﻿function results = opt_layout_fixed_n(spec_in)
% Fixed-count multistage layout search for pyramid TEC layout.
%
% Default run:
%   results = opt_layout_fixed_n();
%
% Custom fixed-n run:
%   spec = struct();
%   spec.stage_count = 5;
%   spec.fixed_n = [322 110 42 14 6];
%   spec.current.I_init = 3.0;
%   spec.targets = struct('DeltaT_target',107,'Qc_target_last',0.5,'ThRes',300);
%   results = opt_layout_fixed_n(spec);
%
% New API:
%   spec.stage_count
%   spec.fixed_n
%   spec.targets.Qc_target_last
%   spec.geometry.L_max_mm
%   spec.geometry.coverage_min
%   spec.geometry.pyramid_gap_min_mm
%   spec.geometry.min_edge_gap_mm
%   spec.soft_prior.fallback_stage_trend
if nargin < 1
    [spec, spec_source] = resolve_runtime_spec_ms();
else
    [spec, spec_source] = resolve_runtime_spec_ms(spec_in);
end

% User-editable entry defaults for no-argument fixed-n layout runs.
entry_cfg = struct();
entry_cfg.stage_count = 5;
entry_cfg.fixed_n = [480, 246, 114, 64, 32];
entry_cfg.current = struct('I_init', 3.8);
entry_cfg.targets = struct('DeltaT_target', 107, 'Qc_target_last', 0.8, 'ThRes', 300);
entry_cfg.plate_k_inplane = 100;
entry_cfg.z_path = struct('enable', false);
spec = merge_struct_recursive_ms(entry_cfg, spec);
spec = apply_default_spec_ms(spec);
spec = prepare_run_output_dir_ms(spec);
spec.use_parallel = ensure_pool_ms(spec);
rng(spec.seed, 'twister');

script_folder = fileparts(mfilename('fullpath'));
if ~contains(path, script_folder)
    addpath(script_folder);
end
mkdir_if_needed_ms(spec.output.output_dir);

fprintf('\n========== optimize_layout_multistage ==========\n');
fprintf('Entry: %s, version=%s\n', mfilename, code_version_stamp_ms());
fprintf('Output: %s\n', spec.output.output_dir);
print_runtime_spec_summary_ms(spec, spec_source);

t0 = tic;
G_template = init_G_params_ms(spec);
G_template.ThRes = spec.targets.ThRes;
G_template.I = spec.current.I_init;
G_template.omega = spec.omega;
G_template.max_inner = spec.max_inner;
G_template.max_outer = spec.max_outer;
G_template.tol_theta = spec.tol_theta;
G_template.tol_g_rel = spec.tol_g_rel;
G_template.Qc_target_last = spec.targets.Qc_target_last;

[ok_run, results_run, fail_reason] = run_layout_pipeline_once_ms(spec, spec_source, G_template, t0);
if ~ok_run
    error('%s', fail_reason);
end
results = results_run;
end

function [ok, results, fail_reason] = run_layout_pipeline_once_ms(spec, spec_source, G_template, t0)
ok = false;
results = struct();
fail_reason = '';

try
    G = G_template;
    soft_prior = init_soft_prior_ms(spec);

    fprintf('[FixedN] using fixed particle counts n=%s, I=%.3f A\n', ...
        vec_to_inline_str_ms(spec.fixed_n), spec.current.I_init);
    count_solution = build_fixed_count_solution_ms(spec, G);
    count_solutions = count_solution;
    G.I = count_solution.I_opt;

    fprintf('[Step2] building coarse jobs ...\n');
    jobs = build_coarse_jobs_ms(spec);
    fprintf('[Step2] coarse jobs: %d\n', numel(jobs));

    warmup_info = struct('enabled', false, 'job_count', 0, 'cand_count_raw', 0, ...
        'cand_count', 0, 'valid_count', 0, 'dedup_removed', 0);
    warmup_dedup_info = struct('enabled', false, 'tol_m', NaN, 'count_in', 0, ...
        'count_out', 0, 'removed', 0);
    if spec.warmup.enable && ~isempty(jobs)
        warmup_info.enabled = true;
        fprintf('[Warmup] screening removed; using fallback prior.\n');
    else
        fprintf('[Warmup] disabled or no jobs. Using fallback prior.\n');
    end

    spec.soft_prior = soft_prior;
    [jobs, jobs_prune_info] = prune_jobs_by_soft_prior_ms(jobs, soft_prior, spec, 'Step2');
    fprintf(['[Step2-Prune] jobs_before=%d, jobs_after=%d, pruned=%d, floor_added=%d, ' ...
        'lock_enabled=%d, fallback_reason=%s\n'], ...
        jobs_prune_info.count_in, jobs_prune_info.count_out, jobs_prune_info.pruned, ...
        jobs_prune_info.floor_added, jobs_prune_info.lock_enabled, jobs_prune_info.fallback_reason);

    fprintf('[Step2] building coarse geometry candidates ...\n');
    coarse_cands = instantiate_candidate_layouts_ms(jobs, spec, count_solution);
    [coarse_cands, coarse_dedup_info] = deduplicate_candidates_by_geometry_ms(coarse_cands, spec, 'Step2');
    if isempty(coarse_cands)
        error('No geometry-feasible coarse candidate.');
    end
    [coarse_cands, supplement_info] = supplement_geometry_candidates_ms(coarse_cands, spec, count_solution);
    [coarse_cands, candidate_batch_info] = limit_candidates_to_batch_ms(coarse_cands, spec, numel(jobs));
    fprintf(['[CandidateBatch] Candidate batch requested=%d, jobs_after_prune=%d, ' ...
        'geometry_candidates_before_batch=%d, geometry_candidates_after_batch=%d\n'], ...
        candidate_batch_info.requested, candidate_batch_info.jobs_after_prune, ...
        candidate_batch_info.geometry_before_batch, candidate_batch_info.geometry_after_batch);
    if candidate_batch_info.geometry_after_batch < candidate_batch_info.requested
        fprintf(['[CandidateBatch-Warning] requested %d geometry-feasible candidates, ' ...
            'but only %d were available after geometry build and dedup.\n'], ...
            candidate_batch_info.requested, candidate_batch_info.geometry_after_batch);
        fprintf(['[CandidateBatch-Warning] FEM parallelism is geometry-limited: ' ...
            'increase parameterized jobs or relax geometry constraints before changing pool_workers.\n']);
    end
    fprintf('[Step2] geometry-feasible candidates: %d\n', numel(coarse_cands));

    t_fem = tic;
    current_calib = empty_current_calib_ms(G.I);
    current_calib.enabled = false;
    G_rank = G;
    fprintf('[Step3A-Full] evaluating full FEM candidates directly ...\n');
    coarse_eval = evaluate_candidates_ms(coarse_cands, G, spec, true);
    all_eval = coarse_eval(:);
    all_valid = all_eval([all_eval.success]);
    if isempty(all_valid)
        summarize_failure_messages_ms(all_eval, 'Step3A-Full');
        write_failed_eval_csv_ms(spec, all_eval, fullfile(spec.output.output_dir, 'all_eval_failed.csv'));
        error('No valid full FEM candidate.');
    end
    all_valid = renormalize_rank_scores_ms(all_valid, spec);
    all_valid = sort_candidates_by_deltaTN_ms(all_valid, spec);
    coarse_valid = all_valid;
    topK = min(max(1, round(spec.output.topK)), numel(all_valid));
    top_eval = all_valid(1:topK);
    runtime_fem_sec = toc(t_fem);

    merged_cands = coarse_cands;

    t_post = tic;
    all_valid_final = all_valid;
    top_postK = min(max(1, round(spec.output.top_postK)), numel(top_eval));
    top_eval_post = top_eval(1:top_postK);
    [top_opt_full, top_uni_full, top_compare, baseline_info] = evaluate_top_with_regular_grid_baseline_ms(top_eval_post, merged_cands, G_rank, spec);
    [top_opt_full, top_uni_full, top_compare, top_eval_post] = ...
        resort_top_post_records_ms(top_opt_full, top_uni_full, top_compare, top_eval_post, spec);
    runtime_post_sec = toc(t_post);

    best_cand = find_candidate_by_id_ms(merged_cands, top_opt_full(1).candidate_id);
    best_full = top_opt_full(1);
    if ~best_full.success
        error('Top-1 full FEM result failed: %s', best_full.message);
    end

    eval_counts = struct();
    eval_counts.candidate_batch_requested = candidate_batch_info.requested;
    eval_counts.candidate_batch_jobs_before = candidate_batch_info.jobs_after_prune;
    eval_counts.candidate_batch_jobs_after = candidate_batch_info.jobs_after_prune;
    eval_counts.jobs_after_prune = candidate_batch_info.jobs_after_prune;
    eval_counts.geometry_candidates_before_batch = candidate_batch_info.geometry_before_batch;
    eval_counts.geometry_candidates_after_batch = candidate_batch_info.geometry_after_batch;
    eval_counts.supplement_info = supplement_info;
    eval_counts.coarse_valid_count = numel(coarse_valid);
    eval_counts.final_valid_count = numel(all_valid_final);
    write_outputs_ms(spec, count_solution, count_solutions, all_eval, coarse_valid, all_valid_final, top_eval, top_eval_post, ...
        best_cand, best_full, top_opt_full, top_uni_full, top_compare, baseline_info, current_calib, eval_counts);

    results = struct();
    results.spec = spec;
    results.spec_source = spec_source;
    results.count_solution = count_solution;
    results.count_solutions = count_solutions;
    results.candidate_batch_info = candidate_batch_info;
    results.all_eval = all_eval;
    results.coarse_valid = coarse_valid;
    results.final_valid = all_valid_final;
    results.top_eval = top_eval;
    results.top_eval_post = top_eval_post;
    results.post_eval_count = numel(top_eval_post);
    results.coarse_eval = coarse_eval;
    results.current_calibration = current_calib;
    results.soft_prior = soft_prior;
    results.warmup_info = warmup_info;
    results.warmup_dedup_info = warmup_dedup_info;
    results.jobs_prune_info = jobs_prune_info;
    results.coarse_dedup_info = coarse_dedup_info;
    results.supplement_info = supplement_info;
    results.top_opt_full = top_opt_full;
    results.top_uniform_full = top_uni_full;
    results.top_compare = top_compare;
    results.baseline_info = baseline_info;
    results.best_candidate = best_cand;
    results.best_metrics = best_full;
    results.runtime_sec = toc(t0);
    results.runtime_fem_sec = runtime_fem_sec;
    results.runtime_post_sec = runtime_post_sec;
    results.output_dir = spec.output.output_dir;

    fprintf('Done. Post Top-1 DeltaTN_mean=%.4f K, DeltaTN_actual=%.4f K\n', ...
        best_full.DeltaTN_mean, best_full.DeltaTN_actual);
    fprintf('Runtime FEM=%.2f s, Top%d post=%.2f s, Total=%.2f s\n', ...
        runtime_fem_sec, top_postK, runtime_post_sec, results.runtime_sec);
    fprintf('===============================================\n\n');
    ok = true;
catch ME
    if ~isempty(ME.stack)
        st = ME.stack(1);
        fail_reason = sprintf('%s (at %s:%d)', ME.message, st.name, st.line);
    else
        fail_reason = ME.message;
    end
end
end

%% ====================== spec/runtime ======================

function [spec, spec_source] = resolve_runtime_spec_ms(spec_in)
if nargin < 1 || isempty(spec_in)
    spec_source = 'default';
    spec = struct();
    return;
end
if ~isstruct(spec_in)
    error('Input must be empty or a struct.');
end
spec_source = 'custom_struct';
spec = merge_struct_recursive_ms(struct(), spec_in);
end

function out = merge_struct_recursive_ms(base, override)
out = base;
if isempty(override)
    return;
end
fn = fieldnames(override);
for i = 1:numel(fn)
    k = fn{i};
    v = override.(k);
    if isstruct(v) && isfield(out, k) && isstruct(out.(k))
        out.(k) = merge_struct_recursive_ms(out.(k), v);
    else
        out.(k) = v;
    end
end
end

function spec = apply_default_spec_ms(spec)
[spec, ~] = optimize_layout_multistage0411_shared_params(spec, 'layout');
end

function n = normalize_fixed_n_ms(raw_n, N)
if nargin < 2 || N < 1
    error('Invalid stage_count for fixed_n.');
end
if nargin < 1 || isempty(raw_n)
    error('spec.fixed_n is required. Automatic particle-count search has been removed from search_upgrade.');
end
n = raw_n(:).';
if numel(n) ~= N
    error('spec.fixed_n length must equal stage_count=%d.', N);
end
if any(~isfinite(n)) || any(n <= 0) || any(abs(n - round(n)) > 0)
    error('spec.fixed_n must contain positive integer particle counts.');
end
n = round(n);
if any(mod(n, 2) ~= 0)
    error('spec.fixed_n must contain even particle counts for the current TEC pair model.');
end
end

function rec = build_fixed_count_solution_ms(spec, G)
N = spec.stage_count;
n = normalize_fixed_n_ms(spec.fixed_n, N);
ratios = NaN(1, max(0, N-1));
for k = 1:(N-1)
    ratios(k) = n(k+1) / max(n(k), eps);
end
rec = struct();
rec.n = n;
rec.ratios = ratios;
rec.I_opt = spec.current.I_init;
rec.DeltaT_0D = NaN;
rec.C = NaN(1, N);
rec.count_rank = 1;
try
    x0 = linspace(spec.targets.ThRes - 15, spec.targets.ThRes - 120, N).';
    [ok0, C0] = solve_0d_at_current_ms(n, spec.current.I_init, G, spec, x0);
    if ok0
        rec.C = C0;
        rec.DeltaT_0D = G.ThRes - C0(end);
    end
catch
    rec.C = NaN(1, N);
    rec.DeltaT_0D = NaN;
end
end

function spec = prepare_run_output_dir_ms(spec)
base_dir = spec.output.output_dir;
if isempty(base_dir)
    base_dir = fullfile(pwd, 'optimize_layout_multistage_output');
end
mkdir_if_needed_ms(base_dir);
dt_tag = numeric_tag_ms(spec.targets.DeltaT_target);
qc_tag = numeric_tag_ms(spec.targets.Qc_target_last);
stamp = datestr(now, 'yyyymmdd_HHMMSS');
n_tag = numeric_vec_tag_ms(spec.fixed_n);
i_tag = numeric_tag_ms(spec.current.I_init);
k_tag = numeric_tag_ms(spec.plate_k_inplane);
run_name = sprintf('layout_fixedn_N%d_n%s_I%s_k%s_DeltaT%s_Qc%s_%s', ...
    spec.stage_count, n_tag, i_tag, k_tag, dt_tag, qc_tag, stamp);
run_dir = fullfile(base_dir, run_name);
suffix = 1;
while exist(run_dir, 'dir')
    run_dir = fullfile(base_dir, sprintf('%s_%02d', run_name, suffix));
    suffix = suffix + 1;
end
spec.output.base_output_dir = base_dir;
spec.output.output_dir = run_dir;
end

function s = numeric_tag_ms(v)
if ~isfinite(v)
    s = 'NaN';
    return;
end
s = char(string(v));
s = strrep(s, '-', 'm');
s = strrep(s, '.', 'p');
s = strrep(s, '+', '');
if isempty(s)
    s = '0';
end
end

function s = code_version_stamp_ms()
s = 'fixedn64_geometry_supplement_64_c4c2shape_20260519';
end

function s = numeric_vec_tag_ms(v)
if isempty(v)
    s = 'empty';
    return;
end
v = v(:).';
parts = cell(1, numel(v));
for i = 1:numel(v)
    parts{i} = numeric_tag_ms(v(i));
end
s = strjoin(parts, '-');
end

function use_parallel = ensure_pool_ms(spec)
use_parallel = isfield(spec, 'use_parallel') && logical(spec.use_parallel);
if ~use_parallel
    return;
end
nw_target = 96;
if isfield(spec, 'parallel') && isstruct(spec.parallel) && isfield(spec.parallel, 'pool_workers')
    nw_target = max(1, round(spec.parallel.pool_workers));
end
local_cap = NaN;
try
    c = parcluster('local');
    if isprop(c, 'NumWorkers')
        local_cap = max(1, round(c.NumWorkers));
    end
catch
    local_cap = NaN;
end
if isfinite(local_cap) && nw_target > local_cap
    fprintf('[Parallel] requested workers=%d exceeds local profile NumWorkers=%d, clamp to %d\n', ...
        nw_target, local_cap, local_cap);
    nw_target = local_cap;
end
pool = gcp('nocreate');
try
    if isempty(pool)
        parpool('local', nw_target);
    elseif pool.NumWorkers ~= nw_target
        fprintf('[Parallel] reset pool workers: %d -> %d\n', pool.NumWorkers, nw_target);
        delete(pool);
        parpool('local', nw_target);
    end
catch ME
    fprintf('[Parallel-Warning] parpool failed; fallback to serial evaluation: %s\n', ME.message);
    use_parallel = false;
end
end

function print_runtime_spec_summary_ms(spec, spec_source)
fprintf('Config source: %s\n', spec_source);
fprintf('stage_count=%d\n', spec.stage_count);
fprintf('Fixed n: %s\n', vec_to_inline_str_ms(spec.fixed_n));
fprintf('Candidate batch: size=%d, top_bottom_n=%d, convergence_gap_K=%.6f\n', ...
    spec.candidate_batch.size, spec.candidate_batch.top_bottom_n, spec.candidate_batch.convergence_gap_K);
fprintf('Output counts: topK=%d, top_postK=%d, top_plotK=%d\n', ...
    spec.output.topK, spec.output.top_postK, spec.output.top_plotK);
fprintf('Targets: DeltaT_target=%.3f K, Qc_target_last=%.3f W, ThRes=%.3f K\n', ...
    spec.targets.DeltaT_target, spec.targets.Qc_target_last, spec.targets.ThRes);
fprintf('L_max_mm: %s\n', vec_to_inline_str_ms(spec.geometry.L_max_mm));
fprintf('coverage_min: %s\n', vec_to_inline_str_ms(spec.geometry.coverage_min));
fprintf('pyramid_gap_min_mm: %s\n', vec_to_inline_str_ms(spec.geometry.pyramid_gap_min_mm));
fprintf('min_edge_gap_mm: %.3f\n', spec.geometry.min_edge_gap_mm);
fprintf('Layout methods: %s\n', strjoin(spec.layout_methods, ','));
fprintf('Shape explore: enable=%d, modes=%s, stage_anis=%d, jitter=%s\n', ...
    logical(spec.shape_explore.enable), strjoin(spec.shape_explore.extra_stage_modes, ','), ...
    logical(spec.shape_explore.stage_anis_enable), vec_to_inline_str_ms(spec.shape_explore.jitter_ratio_list));
fprintf('Warmup: enable=%d, sample_ratio=%.3f, topN=%d\n', ...
    spec.warmup.enable, spec.warmup.sample_ratio, spec.warmup.topN);
fprintf('Current calibration: disabled; direct full FEM mode.\n');
fprintf('Plot: use_global_axis=%d, show_mesh_edges=%d, view_mode=%s, axis_margin=%.3f\n', ...
    spec.output.plot.use_global_axis, spec.output.plot.show_mesh_edges, ...
    spec.output.plot.view_mode, spec.output.plot.axis_margin_ratio);
fprintf('Plot extras: save_overview=%d, save_stage_separate=%d, annotate_substrate_dims=%d, separate_fig_size_px=%s\n', ...
    spec.output.plot.save_overview, spec.output.plot.save_stage_separate, ...
    spec.output.plot.annotate_substrate_dims, vec_to_inline_str_ms(spec.output.plot.separate_fig_size_px));
fprintf('Soft prior fallback trends: %s\n', strjoin(spec.soft_prior.fallback_stage_trend, '/'));
fprintf('Full FEM mesh (stage-wise): nx=%s, ny=%s\n', ...
    vec_to_inline_str_ms(spec.mesh_nx_stage_full), vec_to_inline_str_ms(spec.mesh_ny_stage_full));
fprintf('Plate in-plane thermal conductivity: %.6g W/(m*K)\n', spec.plate_k_inplane);
fprintf('Parallel: use=%d, pool_workers=%d, eval_min_tasks=%d, block=[%d,%d]\n', ...
    spec.use_parallel, spec.parallel.pool_workers, spec.parallel.eval_min_tasks, ...
    spec.parallel.block_min, spec.parallel.block_max);
fprintf('Poisson seeds: %s\n', vec_to_inline_str_ms(spec.poisson.seed_list));
fprintf('Interstage align weights: %s\n', vec_to_inline_str_ms(spec.interstage.align_weights));
fprintf('Job budget: enable=%d\n', logical(spec.job_budget.enable));
if is_z_path_enabled_ms(spec)
    fprintf('Z-path: enable=1, interfaces=%s, t=%s, Rc=%s, sink(k=%.3g,t=%.3g,Rc=%.3g)\n', ...
        vec_to_inline_str_ms(spec.z_path.k_interfaces), ...
        vec_to_inline_str_ms(spec.z_path.t_interface_effs), ...
        vec_to_inline_str_ms(spec.z_path.Rc_interfaces), ...
        spec.z_path.k_sink, spec.z_path.t_sink_eff, spec.z_path.Rc_sink);
else
    fprintf('Z-path: disabled\n');
end
end

function s = vec_to_inline_str_ms(v)
if isempty(v)
    s = '[]';
    return;
end
if iscell(v)
    s = ['[', strjoin(v, ','), ']'];
    return;
end
s = ['[', strtrim(sprintf('%.6g ', v)), ']'];
end

function s = mat_to_inline_str_ms(M)
if isempty(M)
    s = '[]';
    return;
end
rows = cell(size(M,1), 1);
for i = 1:size(M,1)
    rows{i} = sprintf('[%.6g,%.6g]', M(i,1), M(i,2));
end
s = ['{', strjoin(rows, ','), '}'];
end

function s = stage_templates_to_inline_str_ms(templates)
if isempty(templates)
    s = '';
    return;
end
if isstring(templates)
    templates = cellstr(templates);
end
parts = {};
if iscell(templates) && ~isempty(templates) && all(cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)), templates))
    parts{1} = strjoin(cellfun(@char, templates, 'UniformOutput', false), '/');
else
    for i = 1:numel(templates)
        t = templates{i};
        if isstring(t)
            t = cellstr(t);
        end
        if iscell(t)
            parts{end+1} = strjoin(cellfun(@char, t, 'UniformOutput', false), '/'); %#ok<AGROW>
        end
    end
end
s = strjoin(parts, ' | ');
end

function tf = is_z_path_enabled_ms(spec)
tf = isstruct(spec) && isfield(spec, 'z_path') && isstruct(spec.z_path) && ...
    isfield(spec.z_path, 'enable') && logical(spec.z_path.enable);
end

function methods = normalize_layout_methods_ms(raw_methods)
if isempty(raw_methods)
    methods = {'subset_symmetric', 'fixed_c4_grid', 'hex_c6', 'ring_stage3', 'interstage_aligned', 'poisson_disk'};
    return;
end
if ischar(raw_methods) || (isstring(raw_methods) && isscalar(raw_methods))
    raw_methods = {char(raw_methods)};
elseif isstring(raw_methods)
    raw_methods = cellstr(raw_methods);
end
valid = {'subset_symmetric', 'fixed_c4_grid', 'hex_c6', 'gamma_stage2', ...
    'ring_stage3', 'poisson_disk', 'interstage_aligned', 'shape_explore'};
methods = {};
for i = 1:numel(raw_methods)
    k = lower(strtrim(char(raw_methods{i})));
    if ~any(strcmp(k, valid))
        error('Unknown layout method: %s', k);
    end
    if ~any(strcmp(methods, k))
        methods{end+1} = k; %#ok<AGROW>
    end
end
if isempty(methods)
    methods = {'subset_symmetric', 'fixed_c4_grid', 'hex_c6', 'ring_stage3', 'interstage_aligned', 'poisson_disk'};
end
end

function methods = normalize_layout_methods_ms_legacy_unused(raw_methods)
if isempty(raw_methods)
    methods = {'subset_symmetric', 'hex_c6', 'gamma_stage2', ...
        'ring_stage3', 'poisson_disk', 'interstage_aligned'};
    return;
end
if ischar(raw_methods) || (isstring(raw_methods) && isscalar(raw_methods))
    raw_methods = {char(raw_methods)};
elseif isstring(raw_methods)
    raw_methods = cellstr(raw_methods);
end
valid = {'subset_symmetric', 'hex_c6', 'gamma_stage2', ...
    'ring_stage3', 'poisson_disk', 'interstage_aligned'};
methods = {};
for i = 1:numel(raw_methods)
    k = lower(strtrim(char(raw_methods{i})));
    if ~any(strcmp(k, valid))
        error('Unknown layout method: %s', k);
    end
    if ~any(strcmp(methods, k))
        methods{end+1} = k; %#ok<AGROW>
    end
end
if isempty(methods)
    methods = {'subset_symmetric', 'hex_c6'};
end
end

function modes = normalize_mode_list_ms(raw)
if isempty(raw)
    modes = {'center_dense', 'edge_dense'};
    return;
end
if ischar(raw) || (isstring(raw) && isscalar(raw))
    raw = {char(raw)};
elseif isstring(raw)
    raw = cellstr(raw);
end
valid = {'center_dense', 'edge_dense'};
modes = {};
for i = 1:numel(raw)
    k = lower(strtrim(char(raw{i})));
    if any(strcmp(k, valid)) && ~any(strcmp(modes, k))
        modes{end+1} = k; %#ok<AGROW>
    end
end
if isempty(modes)
    modes = {'center_dense', 'edge_dense'};
end
end

function cfg = normalize_shape_explore_cfg_ms(raw, fallback)
cfg = fallback;
if nargin < 2 || ~isstruct(fallback)
    fallback = struct();
end
if nargin >= 1 && isstruct(raw)
    cfg = merge_struct_recursive_ms(fallback, raw);
end
if ~isfield(cfg, 'enable'), cfg.enable = true; end
if ~isfield(cfg, 'extra_stage_modes'), cfg.extra_stage_modes = {'ring_dense','near_center_dense','corner_dense','band_dense_x','band_dense_y','multi_center'}; end
if ~isfield(cfg, 'stage_anis_enable'), cfg.stage_anis_enable = true; end
if ~isfield(cfg, 'anis_ratio_stage_list'), cfg.anis_ratio_stage_list = [0.65, 1.0, 1.55]; end
if ~isfield(cfg, 'ring_radius_ratio_list'), cfg.ring_radius_ratio_list = [0.30, 0.45, 0.60]; end
if ~isfield(cfg, 'ring_width_ratio_list'), cfg.ring_width_ratio_list = [0.12, 0.20]; end
if ~isfield(cfg, 'band_width_ratio_list'), cfg.band_width_ratio_list = [0.20, 0.35]; end
if ~isfield(cfg, 'corner_bias_list'), cfg.corner_bias_list = [0.35, 0.55]; end
if ~isfield(cfg, 'jitter_ratio_list'), cfg.jitter_ratio_list = [0, 0.08, 0.16]; end
if ~isfield(cfg, 'jitter_seed_list'), cfg.jitter_seed_list = [101, 203, 307, 409]; end
if ~isfield(cfg, 'stage_mode_templates'), cfg.stage_mode_templates = {}; end
cfg.enable = logical(cfg.enable);
cfg.stage_anis_enable = logical(cfg.stage_anis_enable);
cfg.extra_stage_modes = normalize_shape_mode_list_ms(cfg.extra_stage_modes);
cfg.anis_ratio_stage_list = normalize_positive_values_ms(cfg.anis_ratio_stage_list, [0.65, 1.0, 1.55]);
cfg.ring_radius_ratio_list = normalize_bounded_values_ms(cfg.ring_radius_ratio_list, [0.30, 0.45, 0.60], 0.05, 0.95);
cfg.ring_width_ratio_list = normalize_bounded_values_ms(cfg.ring_width_ratio_list, [0.12, 0.20], 0.03, 0.80);
cfg.band_width_ratio_list = normalize_bounded_values_ms(cfg.band_width_ratio_list, [0.20, 0.35], 0.03, 0.90);
cfg.corner_bias_list = normalize_bounded_values_ms(cfg.corner_bias_list, [0.35, 0.55], 0.01, 2.0);
cfg.jitter_ratio_list = normalize_bounded_values_ms(cfg.jitter_ratio_list, [0, 0.08, 0.16], 0, 0.45);
cfg.jitter_seed_list = normalize_seed_list_ms(cfg.jitter_seed_list, [101, 203, 307, 409]);
cfg.stage_mode_templates = normalize_stage_mode_templates_ms(cfg.stage_mode_templates, 5);
end

function modes = normalize_shape_mode_list_ms(raw)
if isempty(raw)
    modes = {'ring_dense','near_center_dense','corner_dense','band_dense_x','band_dense_y','multi_center'};
    return;
end
if ischar(raw) || (isstring(raw) && isscalar(raw))
    raw = {char(raw)};
elseif isstring(raw)
    raw = cellstr(raw);
end
valid = {'ring_dense','near_center_dense','center_quad_dense','corner_dense','band_dense_x','band_dense_y','multi_center'};
modes = {};
for i = 1:numel(raw)
    k = lower(strtrim(char(raw{i})));
    if any(strcmp(k, valid)) && ~any(strcmp(modes, k))
        modes{end+1} = k; %#ok<AGROW>
    end
end
if isempty(modes)
    modes = {'ring_dense','near_center_dense','corner_dense','band_dense_x','band_dense_y','multi_center'};
end
end

function templates = normalize_stage_mode_templates_ms(raw, N)
templates = {};
if nargin < 2 || N < 1 || isempty(raw)
    return;
end
if isstring(raw)
    raw = cellstr(raw);
end
valid = {'center_dense','edge_dense','neutral','ring_dense','near_center_dense', ...
    'center_quad_dense','corner_dense','band_dense_x','band_dense_y','multi_center'};
if iscell(raw) && ~isempty(raw) && all(cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)), raw))
    raw = {raw};
end
if ~iscell(raw)
    return;
end
for i = 1:numel(raw)
    t = raw{i};
    if isstring(t)
        t = cellstr(t);
    end
    if ~(iscell(t) && ~isempty(t))
        continue;
    end
    one = cell(1, N);
    ok = true;
    for k = 1:N
        idx = min(k, numel(t));
        key = lower(strtrim(char(string(t{idx}))));
        if ~any(strcmp(key, valid))
            ok = false;
            break;
        end
        one{k} = key;
    end
    if ok
        templates{end+1} = one; %#ok<AGROW>
    end
end
end

function v = normalize_bounded_values_ms(v_in, fallback, lo, hi)
if nargin < 2 || isempty(fallback)
    fallback = lo;
end
if isempty(v_in) || ~isnumeric(v_in)
    v = fallback;
else
    v = double(v_in(:).');
    v = v(isfinite(v));
    if isempty(v)
        v = fallback;
    end
end
v = unique(min(max(v, lo), hi), 'stable');
end

function modes = normalize_symmetry_mode_list_ms(raw)
modes = {'c4'};
end

function modes = normalize_edge_pattern_mode_list_ms(raw)
if isempty(raw)
    modes = {'free', 'edge_spaced'};
    return;
end
if ischar(raw) || (isstring(raw) && isscalar(raw))
    raw = {char(raw)};
elseif isstring(raw)
    raw = cellstr(raw);
end
valid = {'free', 'edge_spaced', 'edge_clean'};
modes = {};
for i = 1:numel(raw)
    k = lower(strtrim(char(raw{i})));
    if any(strcmp(k, valid)) && ~any(strcmp(modes, k))
        modes{end+1} = k; %#ok<AGROW>
    end
end
if isempty(modes)
    modes = {'free', 'edge_spaced'};
end
end

function trends = normalize_trend_cell_ms(raw, N)
if nargin < 2 || N < 1
    N = 1;
end
valid = {'center_heavy', 'edge_heavy', 'neutral'};
if isempty(raw)
    raw = repmat({'neutral'}, 1, N);
end
if ischar(raw) || (isstring(raw) && isscalar(raw))
    raw = {char(raw)};
elseif isstring(raw)
    raw = cellstr(raw);
end
trends = cell(1, N);
for i = 1:N
    idx = min(i, numel(raw));
    v = lower(strtrim(char(raw{idx})));
    if ~any(strcmp(v, valid))
        v = 'neutral';
    end
    trends{i} = v;
end
end

function v = normalize_positive_values_ms(raw, fallback)
if nargin < 2 || isempty(fallback)
    fallback = 1.0;
end
if isempty(raw)
    raw = fallback;
end
v = unique(raw(:).', 'stable');
v = v(isfinite(v) & v > 0);
if isempty(v)
    v = unique(fallback(:).', 'stable');
end
end

function v = normalize_gamma_list_ms(raw, fallback)
if nargin < 2 || isempty(fallback)
    fallback = [-1, -0.5, 0, 0.5, 1];
end
if isempty(raw)
    raw = fallback;
end
v = unique(raw(:).', 'sorted');
v = v(isfinite(v) & v >= -1 & v <= 1);
if isempty(v)
    v = fallback;
end
end

function v = normalize_seed_list_ms(raw, fallback)
if nargin < 2 || isempty(fallback)
    fallback = [42, 137, 271, 618, 1001];
end
if isempty(raw)
    raw = fallback;
end
v = unique(round(double(raw(:).')), 'stable');
v = v(isfinite(v) & v >= 0);
if isempty(v)
    v = unique(round(double(fallback(:).')), 'stable');
end
end

function v = normalize_align_weights_ms(raw, fallback)
if nargin < 2 || isempty(fallback)
    fallback = [0.2, 0.5, 0.8];
end
if isempty(raw)
    raw = fallback;
end
v = unique(double(raw(:).'), 'stable');
v = v(isfinite(v));
if isempty(v)
    v = fallback;
end
v = min(max(v, 0), 1);
if isempty(v)
    v = [0.5];
end
end

function caps = normalize_job_budget_caps_ms(raw, fallback)
caps = fallback;
if nargin < 2 || ~isstruct(fallback)
    return;
end
if nargin < 1 || ~isstruct(raw)
    return;
end
fn = fieldnames(fallback);
for i = 1:numel(fn)
    k = fn{i};
    if ~isfield(raw, k)
        continue;
    end
    v = raw.(k);
    if isempty(v) || ~isfinite(v)
        continue;
    end
    caps.(k) = max(1, round(v));
end
end

function v = fit_len_vec_ms(v_in, N, fallback)
if nargin < 3
    fallback = 1;
end
if isempty(v_in)
    v_in = fallback;
end
v_in = v_in(:).';
if numel(v_in) >= N
    v = v_in(1:N);
else
    v = [v_in, repmat(v_in(end), 1, N - numel(v_in))];
end
end

function v = scalar_with_default_ms(v_in, fallback)
if nargin < 2 || ~isfinite(fallback)
    fallback = 0;
end
if isempty(v_in)
    v = fallback;
    return;
end
v = v_in(1);
if ~isfinite(v)
    v = fallback;
end
end

function v = first_numeric_ms(v_in, fallback)
if nargin < 2
    fallback = NaN;
end
v = fallback;
if isnumeric(v_in) && ~isempty(v_in) && isfinite(v_in(1))
    v = double(v_in(1));
end
end

function G = init_G_params_ms(spec)
[~, G] = optimize_layout_multistage0411_shared_params(spec, 'layout');
end

function [ok, C] = solve_0d_at_current_ms(n, I, G, spec, x0)
ok = false;
C = NaN(1, numel(n));
N = numel(n);
fun = @(x) residual_0d_chain_ms(x, n, I, G, spec);

% Try fsolve first (if Optimization Toolbox exists).
has_fsolve = (exist('fsolve', 'file') == 2) || (exist('fsolve', 'builtin') == 5);
if has_fsolve
    try
        opts = optimoptions('fsolve', 'Display', 'off', ...
            'MaxFunctionEvaluations', 2000, 'MaxIterations', 600, ...
            'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10);
        [sol, fval, exitflag] = fsolve(fun, x0, opts);
        if all(isfinite(sol)) && all(isfinite(fval))
            C_try = sol(:).';
            rnorm = norm(fval);
            if numel(C_try) == N
                % Strict gate (physics-first).
                if exitflag > 0 && rnorm < 1e-5 && validate_0d_solution_ms(C_try, n, I, G, spec)
                    C = C_try;
                    ok = true;
                    return;
                end
                % Relaxed gate (legacy-compatible): prioritize residual consistency.
                if rnorm < 1e-3 && validate_0d_solution_relaxed_ms(C_try, n, I, G, spec)
                    C = C_try;
                    ok = true;
                    return;
                end
            end
        end
    catch
        % Continue with toolbox-free fallback below.
    end
end

% Toolbox-free fallback: small multi-start Gauss-Newton.
seed_mat = build_0d_seed_matrix_ms(x0, spec.targets.ThRes, N);
best_r = inf;
best_C = NaN(1, N);
for si = 1:size(seed_mat, 1)
    [ok_gn, C_gn, r_gn] = solve_0d_gauss_newton_ms(seed_mat(si,:).', fun);
    if ~all(isfinite(C_gn)) || ~isfinite(r_gn)
        continue;
    end
    C_try = C_gn(:).';
    if numel(C_try) ~= N
        continue;
    end
    if r_gn < best_r
        best_r = r_gn;
        best_C = C_try;
    end
    if ~ok_gn
        continue;
    end
    if r_gn < 1e-5 && validate_0d_solution_ms(C_try, n, I, G, spec)
        C = C_try;
        ok = true;
        return;
    end
    if r_gn < 1e-3 && validate_0d_solution_relaxed_ms(C_try, n, I, G, spec)
        C = C_try;
        ok = true;
        return;
    end
end

if isfinite(best_r) && best_r < 2e-3 && all(isfinite(best_C)) && ...
        validate_0d_solution_relaxed_ms(best_C, n, I, G, spec)
    C = best_C;
    ok = true;
end
end

function seeds = build_0d_seed_matrix_ms(x0, ThRes, N)
base = x0(:).';
if numel(base) ~= N || any(~isfinite(base))
    base = linspace(ThRes - 20, ThRes - 120, N);
end
seeds = [ ...
    base; ...
    linspace(ThRes - 10, ThRes - 150, N); ...
    linspace(ThRes - 30, ThRes - 220, N); ...
    max(base - 20, 1.0); ...
    min(base + 20, ThRes + 60) ...
    ];
[~, ia] = unique(round(seeds * 1e6), 'rows', 'stable');
seeds = seeds(sort(ia), :);
end

function [ok, x_best, r_best] = solve_0d_gauss_newton_ms(x0, fun)
ok = false;
x = x0(:);
N = numel(x);
x_best = x;
r_best = inf;
if N < 1 || any(~isfinite(x))
    return;
end
r = fun(x);
if any(~isfinite(r))
    return;
end
rnorm = norm(r);
r_best = rnorm;
x_best = x;
if rnorm < 1e-3
    ok = true;
    return;
end

lambda = 1e-4;
for it = 1:80
    h = 1e-3;
    J = zeros(N, N);
    valid_jac = true;
    for j = 1:N
        xj = x;
        xj(j) = xj(j) + h;
        rj = fun(xj);
        if any(~isfinite(rj))
            valid_jac = false;
            break;
        end
        J(:,j) = (rj - r) / h;
    end
    if ~valid_jac || any(~isfinite(J(:)))
        break;
    end

    A = J.' * J + lambda * eye(N);
    b = J.' * r;
    delta = -A \ b;
    if any(~isfinite(delta))
        break;
    end

    accepted = false;
    step = 1.0;
    for ls = 1:10
        x_try = x + step * delta;
        r_try = fun(x_try);
        if all(isfinite(r_try))
            rn_try = norm(r_try);
            if rn_try < rnorm
                x = x_try;
                r = r_try;
                rnorm = rn_try;
                if rn_try < r_best
                    r_best = rn_try;
                    x_best = x_try;
                end
                lambda = max(1e-10, lambda / 3);
                accepted = true;
                break;
            end
        end
        step = step * 0.5;
    end

    if ~accepted
        lambda = min(1e6, lambda * 10);
        if lambda >= 1e6
            break;
        end
    end

    if rnorm < 1e-3
        ok = true;
        x_best = x;
        r_best = rnorm;
        return;
    end
    if norm(step * delta) < 1e-7
        break;
    end
end

ok = isfinite(r_best) && r_best < 1e-3;
end

function ok = validate_0d_solution_ms(C, n, I, G, spec)
ok = false;
N = numel(C);
if N < 1 || numel(n) ~= N
    return;
end
if any(~isfinite(C)) || any(C <= 0) || any(C >= spec.targets.ThRes + 1e-6)
    return;
end
if N > 1
    for k = 1:N-1
        if C(k) <= C(k+1) + 1e-9
            return;
        end
    end
end
for k = 1:N
    Tc = C(k);
    if k == 1
        Th = spec.targets.ThRes;
    else
        Th = C(k-1);
    end
    if Tc >= Th - 1e-9
        return;
    end
    Qc_k = te_Qc_onecouple_ms(Tc, Th, I, G, k);
    if ~isfinite(Qc_k) || Qc_k < 0
        return;
    end
end
r = residual_0d_chain_ms(C(:), n, I, G, spec);
if ~all(isfinite(r)) || norm(r) > 1e-4
    return;
end
ok = true;
end

function ok = validate_0d_solution_relaxed_ms(C, n, I, G, spec)
ok = false;
N = numel(C);
if N < 1 || numel(n) ~= N
    return;
end
if any(~isfinite(C))
    return;
end

% Loose but bounded temperature range for fallback acceptance.
T_hi = spec.targets.ThRes + 60.0;
T_lo = max(1.0, spec.targets.ThRes - 450.0);
if any(C < T_lo) || any(C > T_hi)
    return;
end

r = residual_0d_chain_ms(C(:), n, I, G, spec);
if ~all(isfinite(r)) || norm(r) > 1e-3
    return;
end
ok = true;
end

function r = residual_0d_chain_ms(x, n, I, G, spec)
N = numel(n);
if numel(x) ~= N || ~all(isfinite(x))
    r = ones(N,1) * 1e12;
    return;
end
[ok_pair, npair] = particle_counts_to_pair_counts_ms(n);
if ~ok_pair
    r = ones(N,1) * 1e12;
    return;
end
Qc = zeros(N,1);
Qh = zeros(N,1);
for k = 1:N
    Tc = x(k);
    if k == 1
        Th = spec.targets.ThRes;
    else
        Th = x(k-1);
    end
    Qc(k) = te_Qc_onecouple_ms(Tc, Th, I, G, k);
    Qh(k) = te_Qh_onecouple_ms(Th, Tc, I, G, k);
end
if ~all(isfinite([Qc; Qh]))
    r = ones(N,1) * 1e12;
    return;
end

r = zeros(N,1);
r(1) = npair(end) * Qc(end) - spec.targets.Qc_target_last;
for k = 2:N
    r(k) = npair(k) * Qh(k) - npair(k-1) * Qc(k-1);
end
end

function [ok, npair] = particle_counts_to_pair_counts_ms(n_vec)
n_vec = n_vec(:).';
npair = NaN(size(n_vec));
ok = ~isempty(n_vec) && all(isfinite(n_vec)) && all(n_vec >= 2) && ...
    all(abs(n_vec - round(n_vec)) < 1e-9) && all(mod(round(n_vec), 2) == 0);
if ok
    npair = round(n_vec) / 2;
end
end

function pairs = make_particle_pair_indices_ms(n_particles)
if ~(isfinite(n_particles) && n_particles >= 2 && mod(round(n_particles), 2) == 0)
    pairs = zeros(0, 2);
    return;
end
idx = 1:round(n_particles);
pairs = reshape(idx, 2, []).';
end

function vals_pair = pair_average_values_ms(vals_particle, pair_idx)
vals_pair = NaN(size(pair_idx, 1), 1);
for p = 1:size(pair_idx, 1)
    vals_pair(p) = mean(vals_particle(pair_idx(p,:)));
end
end

function F = add_heat_to_footprint_elems_ms(F, fpElems, elem, tri, fp_idx, qj)
if fp_idx < 1 || fp_idx > numel(fpElems)
    return;
end
elist = fpElems{fp_idx};
for ee = elist
    nodes = tri(ee,:);
    F(nodes) = F(nodes) + qj * elem.A(ee) / 3;
end
end

%% ====================== Step2: warmup prior / jobs / candidates ======================

function prior = disable_soft_prior_ms(spec)
if nargin < 1 || ~isstruct(spec) || ~isfield(spec, 'stage_count')
    N = 5;
else
    N = spec.stage_count;
end
prior = struct('enable', false, 'stage_trend_prefer', {normalize_trend_cell_ms({}, N)}, ...
    'fallback_stage_trend', {normalize_trend_cell_ms({}, N)}, 'is_fallback', true, ...
    'stage_vote_gap', zeros(1, N), 'stage_weight_gap', zeros(1, N), ...
    'warmup_topn', 0, 'min_keep_per_method', 0);
prior = prior(1);
end

function prior = init_soft_prior_ms(spec)
N = spec.stage_count;
prior = struct('enable', true, ...
    'stage_trend_prefer', {normalize_trend_cell_ms(spec.soft_prior.fallback_stage_trend, N)}, ...
    'fallback_stage_trend', {normalize_trend_cell_ms(spec.soft_prior.fallback_stage_trend, N)}, ...
    'is_fallback', true, ...
    'stage_vote_gap', zeros(1, N), ...
    'stage_weight_gap', zeros(1, N), ...
    'warmup_topn', 0, ...
    'min_keep_per_method', spec.soft_prior.min_keep_per_method);
prior = prior(1);
end

function prior = infer_soft_prior_from_eval_ms(warm_valid_eval, spec)
N = spec.stage_count;
prior = init_soft_prior_ms(spec);
if isempty(warm_valid_eval)
    return;
end
sorted = sort_candidates_by_deltaTN_ms(warm_valid_eval, spec);
topN = min(spec.warmup.topN, numel(sorted));
if topN <= 0
    return;
end
top_arr = sorted(1:topN);
prior.warmup_topn = topN;

fallback_flags = false(1, N);
for k = 1:N
    [trend_out, vote_gap, weight_gap, is_fallback] = infer_stage_trend_prefer_ms( ...
        top_arr, k, prior.fallback_stage_trend{k}, ...
        spec.soft_prior.vote_margin_min, spec.soft_prior.weight_margin_min);
    prior.stage_trend_prefer{k} = trend_out;
    prior.stage_vote_gap(k) = vote_gap;
    prior.stage_weight_gap(k) = weight_gap;
    fallback_flags(k) = is_fallback;
end
prior.is_fallback = any(fallback_flags);
end

function [trend_out, vote_gap, weight_gap, is_fallback] = infer_stage_trend_prefer_ms(arr, stage_idx, fallback_trend, vote_margin, weight_margin)
trend_out = fallback_trend;
vote_gap = 0;
weight_gap = 0;
is_fallback = true;
if isempty(arr)
    return;
end
modes = {'center_heavy', 'edge_heavy', 'neutral'};
vote = zeros(1, numel(modes));
wvote = zeros(1, numel(modes));
for i = 1:numel(arr)
    md = safe_stage_trend_from_eval_ms(arr(i), stage_idx);
    j = find(strcmp(md, modes), 1, 'first');
    if isempty(j)
        continue;
    end
    vote(j) = vote(j) + 1;
    wvote(j) = wvote(j) + 1 / i;
end
if all(vote == 0)
    return;
end
[vote_max, imax_vote] = max(vote);
vote_other = max(vote(setdiff(1:numel(modes), imax_vote)));
vote_gap = vote_max - vote_other;
[wmax, imax_w] = max(wvote);
wother = max(wvote(setdiff(1:numel(modes), imax_w)));
weight_gap = wmax - wother;
if imax_vote ~= imax_w
    return;
end
if vote_gap > vote_margin && weight_gap >= weight_margin
    trend_out = modes{imax_vote};
    is_fallback = false;
end
end

function trend = safe_stage_trend_from_eval_ms(ev, stage_idx)
trend = 'neutral';
if isfield(ev, 'stage_trends') && numel(ev.stage_trends) >= stage_idx
    trend = lower(strtrim(char(ev.stage_trends{stage_idx})));
    if isempty(trend)
        trend = 'neutral';
    end
end
end

function [jobs_out, info] = prune_jobs_by_soft_prior_ms(jobs, soft_prior, spec, tag)
if nargin < 4 || isempty(tag)
    tag = 'Step2';
end
jobs_out = jobs;
info = struct('count_in', numel(jobs), 'count_out', numel(jobs), 'pruned', 0, ...
    'floor_added', 0, 'min_keep_per_method', 0, 'lock_enabled', false, ...
    'fallback_reason', 'none', 'prefer_trends', {{}} );

if isempty(jobs)
    info.fallback_reason = 'no_jobs';
    return;
end
if ~isstruct(soft_prior) || ~isfield(soft_prior, 'is_fallback') || logical(soft_prior.is_fallback)
    info.fallback_reason = 'prior_fallback';
    return;
end
N = spec.stage_count;
prefer = normalize_trend_cell_ms(soft_prior.stage_trend_prefer, N);
info.prefer_trends = prefer;

keep_mask = false(numel(jobs), 1);
for i = 1:numel(jobs)
    jt = normalize_trend_cell_ms(jobs(i).stage_trends, N);
    keep_mask(i) = true;
    for k = 1:N
        if ~strcmpi(prefer{k}, 'neutral') && ~strcmpi(jt{k}, prefer{k})
            keep_mask(i) = false;
            break;
        end
    end
end

min_keep_per_method = 0;
if isfield(soft_prior, 'min_keep_per_method')
    min_keep_per_method = max(0, round(soft_prior.min_keep_per_method));
end
info.min_keep_per_method = min_keep_per_method;
if min_keep_per_method > 0
    method_list = string({jobs.layout_method});
    u_methods = unique(method_list, 'stable');
    for iu = 1:numel(u_methods)
        idx_m = find(method_list == u_methods(iu));
        need_n = min(min_keep_per_method, numel(idx_m));
        cur_n = sum(keep_mask(idx_m));
        add_n = need_n - cur_n;
        if add_n <= 0
            continue;
        end
        cand_idx = idx_m(~keep_mask(idx_m));
        if isempty(cand_idx)
            continue;
        end
        s = -inf(numel(cand_idx), 1);
        for j = 1:numel(cand_idx)
            jt = normalize_trend_cell_ms(jobs(cand_idx(j)).stage_trends, N);
            match_n = 0;
            for k = 1:N
                if strcmpi(jt{k}, prefer{k})
                    match_n = match_n + 1;
                end
            end
            s(j) = match_n;
        end
        [~, ord] = sort(s, 'descend');
        take = min(add_n, numel(ord));
        add_idx = cand_idx(ord(1:take));
        keep_mask(add_idx) = true;
        info.floor_added = info.floor_added + take;
    end
end

if ~any(keep_mask)
    info.fallback_reason = 'empty_after_lock';
    return;
end
jobs_out = jobs(keep_mask);
info.count_out = numel(jobs_out);
info.pruned = info.count_in - info.count_out;
info.lock_enabled = true;
if info.pruned <= 0
    info.fallback_reason = 'none_pruned';
else
    info.fallback_reason = 'none';
end

fprintf('[%s-Prune] prefer=%s, min_keep_per_method=%d, floor_added=%d\n', ...
    tag, strjoin(prefer, '/'), min_keep_per_method, info.floor_added);
end

function jobs = build_coarse_jobs_ms(spec)
jobs = empty_job_struct_ms();
id = 0;
N = spec.stage_count;
mode_templates = build_mode_templates_ms(N, spec.soft_prior.fallback_stage_trend);
min_spacing_ratio = 1.6;

for lm = 1:numel(spec.layout_methods)
    method_name = spec.layout_methods{lm};

    if strcmpi(method_name, 'poisson_disk')
        seed_list = spec.poisson.seed_list;
        for it = 1:numel(mode_templates)
            base_modes = mode_templates{it};
            for iseed = 1:numel(seed_list)
                id = id + 1;
                jobs(end+1) = pack_job_ms(id, method_name, base_modes, 'c4', 'free', ... %#ok<AGROW>
                    spec.coarse_s_dense_list(1), NaN, NaN, 1.0, NaN, NaN, seed_list(iseed), NaN);
            end
        end
        continue;
    end

    if strcmpi(method_name, 'interstage_aligned')
        align_weights = spec.interstage.align_weights;
        base_modes = repmat({'center_dense'}, 1, N);
        for iw = 1:numel(align_weights)
            id = id + 1;
            jobs(end+1) = pack_job_ms(id, method_name, base_modes, 'c4', 'free', ... %#ok<AGROW>
                spec.coarse_s_dense_list(1), NaN, NaN, 1.0, NaN, NaN, NaN, align_weights(iw));
        end
        continue;
    end

    if strcmpi(method_name, 'fixed_c4_grid')
        for it = 1:numel(mode_templates)
            base_modes = mode_templates{it};
            for idd = 1:numel(spec.coarse_s_dense_list)
                for ia = 1:numel(spec.coarse_anis_ratio_list)
                    id = id + 1;
                    jobs(end+1) = pack_job_ms(id, method_name, base_modes, 'c4', 'free', ... %#ok<AGROW>
                        spec.coarse_s_dense_list(idd), spec.coarse_s_sparse_list(1), ...
                        spec.coarse_expo_list(1), spec.coarse_anis_ratio_list(ia), NaN, NaN);
                end
            end
        end
        continue;
    end

    if strcmpi(method_name, 'hex_c6') && spec.hex.enable
        for it = 1:numel(mode_templates)
            base_modes = mode_templates{it};
            for ish = 1:numel(spec.hex.s_list)
                for ia = 1:numel(spec.coarse_anis_ratio_list)
                    id = id + 1;
                    jobs(end+1) = pack_job_ms(id, method_name, base_modes, 'c4', 'free', ... %#ok<AGROW>
                        spec.hex.s_list(ish), spec.coarse_s_sparse_list(1), ...
                        spec.coarse_expo_list(1), spec.coarse_anis_ratio_list(ia), NaN, NaN);
                end
            end
        end
        continue;
    end

    if strcmpi(method_name, 'gamma_stage2')
        for it = 1:numel(mode_templates)
            base_modes = mode_templates{it};
            for ig = 1:numel(spec.gamma_list)
                for idd = 1:numel(spec.coarse_s_dense_list)
                    for ia = 1:numel(spec.coarse_anis_ratio_list)
                        id = id + 1;
                        jobs(end+1) = pack_job_ms(id, method_name, base_modes, 'c4', 'free', ... %#ok<AGROW>
                            spec.coarse_s_dense_list(idd), spec.coarse_s_sparse_list(1), ...
                            spec.coarse_expo_list(1), spec.coarse_anis_ratio_list(ia), ...
                            spec.gamma_list(ig), NaN);
                    end
                end
            end
        end
        continue;
    end

    if strcmpi(method_name, 'ring_stage3') && spec.ring.enable
        for it = 1:numel(mode_templates)
            base_modes = mode_templates{it};
            for idd = 1:numel(spec.coarse_s_dense_list)
                for ia = 1:numel(spec.coarse_anis_ratio_list)
                    id = id + 1;
                    jobs(end+1) = pack_job_ms(id, method_name, base_modes, 'c4', 'free', ... %#ok<AGROW>
                        spec.coarse_s_dense_list(idd), spec.coarse_s_sparse_list(1), ...
                        spec.coarse_expo_list(1), spec.coarse_anis_ratio_list(ia), NaN, NaN);
                end
            end
        end
        continue;
    end

    for it = 1:numel(mode_templates)
        stage_modes = mode_templates{it};
        for isym = 1:numel(spec.symmetry_mode_list)
            for iedge = 1:numel(spec.edge_pattern_mode_list)
                for idd = 1:numel(spec.coarse_s_dense_list)
                    for iss = 1:numel(spec.coarse_s_sparse_list)
                        if spec.coarse_s_sparse_list(iss) <= spec.coarse_s_dense_list(idd)
                            continue;
                        end
                        if (spec.coarse_s_sparse_list(iss) / spec.coarse_s_dense_list(idd)) < min_spacing_ratio
                            continue;
                        end
                        for ie = 1:numel(spec.coarse_expo_list)
                            for ia = 1:numel(spec.coarse_anis_ratio_list)
                                id = id + 1;
                                jobs(end+1) = pack_job_ms(id, method_name, stage_modes, ... %#ok<AGROW>
                                    spec.symmetry_mode_list{isym}, spec.edge_pattern_mode_list{iedge}, ...
                                    spec.coarse_s_dense_list(idd), spec.coarse_s_sparse_list(iss), ...
                                    spec.coarse_expo_list(ie), spec.coarse_anis_ratio_list(ia), NaN, NaN);
                            end
                        end
                    end
                end
            end
        end
    end
end

if spec.method_mix.enable
    method_combos = build_stage_method_combinations_ms(spec.layout_methods, N);
    s_dense_ref = spec.coarse_s_dense_list(max(1, ceil(numel(spec.coarse_s_dense_list) / 2)));
    s_sparse_ref = spec.coarse_s_sparse_list(end);
    expo_ref = spec.coarse_expo_list(max(1, ceil(numel(spec.coarse_expo_list) / 2)));
    anis_ref = 1.0;
    gamma_ref = spec.gamma_list(max(1, ceil(numel(spec.gamma_list) / 2)));
    for ic = 1:numel(method_combos)
        for it = 1:numel(mode_templates)
            id = id + 1;
            jobs(end+1) = pack_job_ms(id, 'mixed_stage_methods', mode_templates{it}, 'c4', 'free', ... %#ok<AGROW>
                s_dense_ref, s_sparse_ref, expo_ref, anis_ref, gamma_ref, NaN, NaN, NaN, method_combos{ic});
        end
    end
end

if spec.shape_explore.enable
    [shape_jobs, id] = build_shape_explore_jobs_ms(spec, mode_templates, id);
    if ~isempty(shape_jobs)
        jobs = append_job_structs_ms(jobs, shape_jobs);
    end
end

jobs = apply_job_budget_ms(jobs, spec);
end

function [jobs, id] = build_shape_explore_jobs_ms(spec, mode_templates, id)
jobs = empty_job_struct_ms();
N = spec.stage_count;
cfg = spec.shape_explore;
if ~cfg.enable || isempty(cfg.extra_stage_modes)
    return;
end
s_ref = spec.coarse_s_dense_list(max(1, ceil(numel(spec.coarse_s_dense_list) / 2)));
s_sparse_ref = spec.coarse_s_sparse_list(max(1, ceil(numel(spec.coarse_s_sparse_list) / 2)));
expo_ref = spec.coarse_expo_list(max(1, ceil(numel(spec.coarse_expo_list) / 2)));
base_anis_list = unique([1.0, spec.coarse_anis_ratio_list(:).'], 'stable');
if cfg.stage_anis_enable
    stage_anis_pool = cfg.anis_ratio_stage_list;
else
    stage_anis_pool = 1.0;
end
base_templates = normalize_stage_mode_templates_for_N_ms(shape_cfg_field_ms(cfg, 'stage_mode_templates', {}), N);
if isempty(base_templates)
    base_templates = mode_templates;
end
if isempty(base_templates)
    base_templates = {repmat({'center_dense'}, 1, N)};
end

for im = 1:numel(cfg.extra_stage_modes)
    shape_mode = cfg.extra_stage_modes{im};
    radius_list = NaN;
    ring_width_list = NaN;
    band_width_list = NaN;
    corner_bias_list = NaN;
    if any(strcmp(shape_mode, {'ring_dense','multi_center','near_center_dense','center_quad_dense'}))
        radius_list = cfg.ring_radius_ratio_list;
        ring_width_list = cfg.ring_width_ratio_list;
    elseif any(strcmp(shape_mode, {'band_dense_x','band_dense_y'}))
        band_width_list = cfg.band_width_ratio_list;
    elseif strcmp(shape_mode, 'corner_dense')
        corner_bias_list = cfg.corner_bias_list;
    end
    for it = 1:numel(base_templates)
        stage_modes = base_templates{it};
        if ~template_contains_shape_mode_ms(stage_modes, shape_mode)
            for k = max(1, ceil(N/2)):N
                stage_modes{k} = shape_mode;
            end
        end
        for ia = 1:min(2, numel(base_anis_list))
            for isa = 1:numel(stage_anis_pool)
                stage_anis = make_stage_anis_vector_ms(stage_anis_pool(isa), N);
                for ir = 1:numel(radius_list)
                    for iw = 1:numel(ring_width_list)
                        for ib = 1:numel(band_width_list)
                            for ic = 1:numel(corner_bias_list)
                                jitter_take = min(2, numel(cfg.jitter_ratio_list));
                                for ij = 1:jitter_take
                                    seed_idx = mod(id + ij - 1, numel(cfg.jitter_seed_list)) + 1;
                                    shape_cfg = struct( ...
                                        'shape_mode', shape_mode, ...
                                        'stage_anis', stage_anis, ...
                                        'ring_radius_ratio', radius_list(ir), ...
                                        'ring_width_ratio', ring_width_list(iw), ...
                                        'band_width_ratio', band_width_list(ib), ...
                                        'corner_bias', corner_bias_list(ic), ...
                                        'jitter_ratio', cfg.jitter_ratio_list(ij), ...
                                        'jitter_seed', cfg.jitter_seed_list(seed_idx));
                                    id = id + 1;
                                    jobs(end+1) = pack_job_ms(id, 'shape_explore', stage_modes, 'c4', 'free', ... %#ok<AGROW>
                                        s_ref, s_sparse_ref, expo_ref, base_anis_list(ia), NaN, NaN, NaN, NaN, ...
                                        repmat({'shape_explore'}, 1, N), shape_cfg);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
end

function jobs_out = append_job_structs_ms(jobs_base, jobs_add)
if isempty(jobs_base)
    jobs_out = jobs_add(:).';
    return;
end
if isempty(jobs_add)
    jobs_out = jobs_base(:).';
    return;
end
jobs_out = [jobs_base(:).', jobs_add(:).'];
end

function v = shape_cfg_field_ms(cfg, field_name, default_v)
v = default_v;
if isstruct(cfg) && isfield(cfg, field_name)
    v = cfg.(field_name);
end
end

function templates = normalize_stage_mode_templates_for_N_ms(raw, N)
templates = {};
if nargin < 2 || N < 1 || isempty(raw)
    return;
end
if isstring(raw)
    raw = cellstr(raw);
end
if iscell(raw) && ~isempty(raw) && all(cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)), raw))
    raw = {raw};
end
if ~iscell(raw)
    return;
end
valid = {'center_dense','edge_dense','neutral','ring_dense','near_center_dense', ...
    'center_quad_dense','corner_dense','band_dense_x','band_dense_y','multi_center'};
for i = 1:numel(raw)
    t = raw{i};
    if isstring(t)
        t = cellstr(t);
    end
    if ~(iscell(t) && ~isempty(t))
        continue;
    end
    one = cell(1, N);
    ok = true;
    for k = 1:N
        idx = min(k, numel(t));
        key = lower(strtrim(char(string(t{idx}))));
        if ~any(strcmp(key, valid))
            ok = false;
            break;
        end
        one{k} = key;
    end
    if ok
        templates{end+1} = one; %#ok<AGROW>
    end
end
end

function tf = template_contains_shape_mode_ms(stage_modes, shape_mode)
tf = false;
if isempty(stage_modes)
    return;
end
shape_mode = lower(strtrim(char(string(shape_mode))));
for i = 1:numel(stage_modes)
    if strcmp(lower(strtrim(char(string(stage_modes{i})))), shape_mode)
        tf = true;
        return;
    end
end
end

function stage_anis = make_stage_anis_vector_ms(anis_val, N)
anis_val = max(0.1, scalar_with_default_ms(anis_val, 1.0));
stage_anis = ones(1, N);
if N >= 1
    stage_anis(1:2:end) = anis_val;
end
if N >= 2
    stage_anis(2:2:end) = 1 / max(anis_val, eps);
end
end

function combos = build_stage_method_combinations_ms(method_pool, N)
combos = {};
if nargin < 2 || N < 1
    return;
end
if isempty(method_pool)
    return;
end
pool = normalize_layout_methods_ms(method_pool);
M = numel(pool);
total = M ^ N;
combos = cell(total, 1);
for idx = 0:(total - 1)
    code = idx;
    c = cell(1, N);
    for k = 1:N
        pick = mod(code, M) + 1;
        c{k} = pool{pick};
        code = floor(code / M);
    end
    combos{idx + 1} = c;
end
end
function jobs_out = apply_job_budget_ms(jobs, spec)
jobs_out = jobs;
if isempty(jobs)
    return;
end
if ~isfield(spec, 'job_budget') || ~isstruct(spec.job_budget) || ~logical(spec.job_budget.enable)
    return;
end
caps = spec.job_budget.method_caps;
method_list = string({jobs.layout_method});
u_methods = unique(method_list, 'stable');
keep_mask = false(numel(jobs), 1);
for i = 1:numel(u_methods)
    idx = find(method_list == u_methods(i));
    cap_i = get_job_method_cap_ms(caps, char(u_methods(i)), numel(idx));
    if cap_i >= numel(idx)
        keep_mask(idx) = true;
        continue;
    end
    pick_local = stratified_pick_jobs_ms(jobs(idx), cap_i);
    keep_mask(idx(pick_local)) = true;
end
jobs_out = jobs(keep_mask);
fprintf('[JobBudget] enabled=%d, jobs_before=%d, jobs_after=%d, pruned=%d\n', ...
    logical(spec.job_budget.enable), numel(jobs), numel(jobs_out), numel(jobs) - numel(jobs_out));
fprintf('[JobBudget] kept_by_method: %s\n', summarize_job_method_counts_ms(jobs_out));
end

function [jobs_out, info] = limit_jobs_to_candidate_batch_ms(jobs, spec)
jobs_out = jobs;
target_n = max(1, round(spec.candidate_batch.size));
info = struct('requested', target_n, 'count_in', numel(jobs), 'count_out', numel(jobs), 'limited', false);
if isempty(jobs)
    info.count_out = 0;
    return;
end
if numel(jobs) <= target_n
    return;
end
pick_mask = stratified_pick_jobs_ms(jobs, target_n);
if ~any(pick_mask)
    pick_mask = false(numel(jobs), 1);
    pick_mask(1:min(target_n, numel(jobs))) = true;
end
jobs_out = jobs(pick_mask);
info.count_out = numel(jobs_out);
info.limited = true;
end

function [cands_out, info] = limit_candidates_to_batch_ms(cands, spec, jobs_after_prune)
if nargin < 3 || ~isfinite(jobs_after_prune)
    jobs_after_prune = NaN;
end
cands_out = cands;
target_n = max(1, round(spec.candidate_batch.size));
info = struct('requested', target_n, ...
    'jobs_after_prune', jobs_after_prune, ...
    'geometry_before_batch', numel(cands), ...
    'geometry_after_batch', numel(cands), ...
    'count_in', numel(cands), ...
    'count_out', numel(cands), ...
    'limited', false);
if isempty(cands)
    info.geometry_after_batch = 0;
    info.count_out = 0;
    return;
end
if numel(cands) <= target_n
    return;
end
pick_mask = stratified_pick_candidates_ms(cands, target_n);
if ~any(pick_mask)
    pick_mask = false(numel(cands), 1);
    pick_mask(1:min(target_n, numel(cands))) = true;
end
cands_out = cands(pick_mask);
info.geometry_after_batch = numel(cands_out);
info.count_out = numel(cands_out);
info.limited = true;
end

function [cands_out, info] = supplement_geometry_candidates_ms(cands, spec, count_solution)
target_n = max(1, round(spec.candidate_batch.size));
quota = supplement_quota_ms(target_n);
max_attempts = 4 * target_n;
info = struct('enabled', false, 'target', target_n, 'max_attempts', max_attempts, ...
    'initial_count', numel(cands), 'final_count', numel(cands), ...
    'attempted_jobs', 0, 'accepted_total', 0, ...
    'quota_c4_main', quota.c4_main, 'quota_c2', quota.c2_explore, 'quota_shape', quota.shape_edge_anis, ...
    'initial_c4_main', 0, 'initial_c2', 0, 'initial_shape', 0, ...
    'final_c4_main', 0, 'final_c2', 0, 'final_shape', 0, ...
    'accepted_c4_main', 0, 'accepted_c2', 0, 'accepted_shape', 0, ...
    'stop_reason', 'not_needed');
cands_out = cands;
if isempty(cands_out) || numel(cands_out) >= target_n
    info = update_supplement_counts_ms(info, cands_out, 'initial');
    info = update_supplement_counts_ms(info, cands_out, 'final');
    return;
end

info.enabled = true;
info.stop_reason = 'attempt_budget_exhausted';
info = update_supplement_counts_ms(info, cands_out, 'initial');
dcfg = get_candidate_dedup_config_ms(spec);
keys = containers.Map('KeyType', 'char', 'ValueType', 'logical');
for i = 1:numel(cands_out)
    keys(char(make_candidate_geometry_key_ms(cands_out(i), dcfg.tol_m))) = true;
end

group_order = {'c4_main', 'c2_explore', 'shape_edge_anis'};
for ig = 1:numel(group_order)
    group_name = group_order{ig};
    if numel(cands_out) >= target_n || info.attempted_jobs >= max_attempts
        break;
    end
    missing_group = supplement_group_missing_ms(cands_out, quota, group_name);
    if missing_group <= 0
        continue;
    end
    remaining_attempts = min(max_attempts - info.attempted_jobs, ...
        supplement_group_attempt_cap_ms(group_name, missing_group, target_n));
    jobs = build_supplement_jobs_ms(spec, group_name, remaining_attempts);
    if isempty(jobs)
        continue;
    end
    info.attempted_jobs = info.attempted_jobs + numel(jobs);
    cand_try = instantiate_candidate_layouts_ms(jobs, spec, count_solution);
    if isempty(cand_try)
        continue;
    end
    for j = 1:numel(cand_try)
        if numel(cands_out) >= target_n
            break;
        end
        if supplement_group_missing_ms(cands_out, quota, group_name) <= 0
            break;
        end
        k = char(make_candidate_geometry_key_ms(cand_try(j), dcfg.tol_m));
        if isKey(keys, k)
            continue;
        end
        keys(k) = true;
        cands_out(end+1) = cand_try(j); %#ok<AGROW>
        info.accepted_total = info.accepted_total + 1;
        switch group_name
            case 'c4_main'
                info.accepted_c4_main = info.accepted_c4_main + 1;
            case 'c2_explore'
                info.accepted_c2 = info.accepted_c2 + 1;
            case 'shape_edge_anis'
                info.accepted_shape = info.accepted_shape + 1;
        end
    end
end

info.final_count = numel(cands_out);
info = update_supplement_counts_ms(info, cands_out, 'final');
if info.final_count >= target_n
    info.stop_reason = 'target_reached';
elseif info.attempted_jobs >= max_attempts
    info.stop_reason = 'attempt_budget_exhausted';
else
    info.stop_reason = 'candidate_pool_exhausted';
end
fprintf(['[CandidateSupplement] enabled=%d, initial=%d, final=%d, target=%d, ' ...
    'attempted_jobs=%d/%d, accepted=%d, groups(c4=%d,c2=%d,shape=%d), stop=%s\n'], ...
    info.enabled, info.initial_count, info.final_count, info.target, ...
    info.attempted_jobs, info.max_attempts, info.accepted_total, ...
    info.final_c4_main, info.final_c2, info.final_shape, info.stop_reason);
end

function cap = supplement_group_attempt_cap_ms(group_name, missing_group, target_n)
missing_group = max(0, round(missing_group));
switch group_name
    case 'c4_main'
        cap = max(32, 4 * missing_group);
    case 'c2_explore'
        cap = max(24, 4 * missing_group);
    otherwise
        cap = max(32, 4 * missing_group);
end
cap = min(cap, 4 * max(1, round(target_n)));
end

function info = empty_supplement_info_ms()
info = struct('enabled', false, 'target', NaN, 'max_attempts', NaN, ...
    'initial_count', NaN, 'final_count', NaN, ...
    'attempted_jobs', 0, 'accepted_total', 0, ...
    'quota_c4_main', NaN, 'quota_c2', NaN, 'quota_shape', NaN, ...
    'initial_c4_main', NaN, 'initial_c2', NaN, 'initial_shape', NaN, ...
    'final_c4_main', NaN, 'final_c2', NaN, 'final_shape', NaN, ...
    'accepted_c4_main', 0, 'accepted_c2', 0, 'accepted_shape', 0, ...
    'stop_reason', 'unknown');
end

function quota = supplement_quota_ms(target_n)
quota = struct('c4_main', 32, 'c2_explore', 15, 'shape_edge_anis', 17);
if target_n == 64
    return;
end
scale = target_n / 64;
quota.c4_main = max(1, round(32 * scale));
quota.c2_explore = max(0, round(15 * scale));
quota.shape_edge_anis = max(0, target_n - quota.c4_main - quota.c2_explore);
end

function info = update_supplement_counts_ms(info, cands, phase)
cnt = count_candidate_groups_ms(cands);
if strcmpi(phase, 'initial')
    info.initial_c4_main = cnt.c4_main;
    info.initial_c2 = cnt.c2_explore;
    info.initial_shape = cnt.shape_edge_anis;
else
    info.final_c4_main = cnt.c4_main;
    info.final_c2 = cnt.c2_explore;
    info.final_shape = cnt.shape_edge_anis;
    info.final_count = numel(cands);
end
end

function missing = supplement_group_missing_ms(cands, quota, group_name)
cnt = count_candidate_groups_ms(cands);
switch group_name
    case 'c4_main'
        missing = quota.c4_main - cnt.c4_main;
    case 'c2_explore'
        missing = quota.c2_explore - cnt.c2_explore;
    otherwise
        missing = quota.shape_edge_anis - cnt.shape_edge_anis;
end
missing = max(0, missing);
end

function cnt = count_candidate_groups_ms(cands)
cnt = struct('c4_main', 0, 'c2_explore', 0, 'shape_edge_anis', 0, 'other', 0);
for i = 1:numel(cands)
    g = classify_candidate_search_group_ms(cands(i));
    switch g
        case 'c4_main'
            cnt.c4_main = cnt.c4_main + 1;
        case 'c2_explore'
            cnt.c2_explore = cnt.c2_explore + 1;
        case 'shape_edge_anis'
            cnt.shape_edge_anis = cnt.shape_edge_anis + 1;
        otherwise
            cnt.other = cnt.other + 1;
    end
end
end

function group_name = classify_candidate_search_group_ms(cand)
group_name = 'other';
if isstruct(cand) && isfield(cand, 'search_group') && ~isempty(cand.search_group)
    raw = lower(strtrim(char(string(cand.search_group))));
    if any(strcmp(raw, {'c4_main','c2_explore','shape_edge_anis'}))
        group_name = raw;
        return;
    end
end
sym = '';
if isstruct(cand) && isfield(cand, 'symmetry_mode')
    sym = lower(strtrim(char(string(cand.symmetry_mode))));
end
if is_c2_symmetry_mode_ms(sym)
    group_name = 'c2_explore';
    return;
end
layout_method = '';
if isstruct(cand) && isfield(cand, 'layout_method')
    layout_method = lower(strtrim(char(string(cand.layout_method))));
end
shape_mode = '';
if isstruct(cand) && isfield(cand, 'shape_mode')
    shape_mode = lower(strtrim(char(string(cand.shape_mode))));
end
edge_mode = '';
if isstruct(cand) && isfield(cand, 'edge_pattern_mode')
    edge_mode = lower(strtrim(char(string(cand.edge_pattern_mode))));
end
anis = NaN;
if isstruct(cand) && isfield(cand, 'anis_ratio')
    anis = cand.anis_ratio;
end
if strcmp(layout_method, 'shape_explore') || ~isempty(shape_mode) || ...
        (~isempty(edge_mode) && ~strcmp(edge_mode, 'free')) || ...
        (isfinite(anis) && abs(anis - 1.0) > 1e-9)
    group_name = 'shape_edge_anis';
elseif strcmp(sym, 'c4') || isempty(sym)
    group_name = 'c4_main';
end
end

function tf = is_c2_symmetry_mode_ms(symmetry_mode)
s = lower(strtrim(char(string(symmetry_mode))));
tf = any(strcmp(s, {'c2_lr', 'c2_ud'}));
end

function jobs = build_supplement_jobs_ms(spec, group_name, max_jobs)
jobs = empty_job_struct_ms();
if max_jobs <= 0
    return;
end
switch group_name
    case 'c4_main'
        jobs = build_c4_supplement_jobs_ms(spec, max_jobs);
    case 'c2_explore'
        jobs = build_c2_supplement_jobs_ms(spec, max_jobs);
    case 'shape_edge_anis'
        jobs = build_shape_edge_anis_supplement_jobs_ms(spec, max_jobs);
end
end

function jobs = build_c4_supplement_jobs_ms(spec, max_jobs)
jobs = empty_job_struct_ms();
N = spec.stage_count;
templates = {
    repmat({'center_dense'}, 1, N), ...
    supplement_template_ms(N, {'center_dense','center_dense','center_dense','edge_dense','edge_dense'}), ...
    supplement_template_ms(N, {'center_dense','center_dense','edge_dense','edge_dense','edge_dense'})};
s_dense_list = unique([spec.coarse_s_dense_list(:).', 2.1e-3, 2.3e-3, 2.5e-3], 'stable');
s_sparse_list = unique([spec.coarse_s_sparse_list(:).', 3.6e-3, 4.2e-3, 4.8e-3], 'stable');
expo_list = unique([spec.coarse_expo_list(:).', 2.5, 5.0, 7.0], 'stable');
id = supplement_candidate_id_base_ms('c4_main');
for it = 1:numel(templates)
    for idd = 1:numel(s_dense_list)
        for iss = 1:numel(s_sparse_list)
            if s_sparse_list(iss) <= s_dense_list(idd) || (s_sparse_list(iss) / s_dense_list(idd)) < 1.45
                continue;
            end
            for ie = 1:numel(expo_list)
                if numel(jobs) >= max_jobs, return; end
                id = id + 1;
                job = pack_job_ms(id, 'subset_symmetric', templates{it}, 'c4', 'free', ...
                    s_dense_list(idd), s_sparse_list(iss), expo_list(ie), 1.0, NaN, NaN);
                job.search_group = 'c4_main';
                jobs(end+1) = job; %#ok<AGROW>
            end
        end
    end
end
end

function jobs = build_c2_supplement_jobs_ms(spec, max_jobs)
jobs = empty_job_struct_ms();
N = spec.stage_count;
templates = {
    repmat({'center_dense'}, 1, N), ...
    supplement_template_ms(N, {'center_dense','center_dense','center_dense','edge_dense','edge_dense'}), ...
    supplement_template_ms(N, {'center_dense','center_dense','edge_dense','center_dense','edge_dense'})};
sym_list = {'c2_lr', 'c2_ud'};
anis_list = unique([0.70, 0.85, 1.0, 1.20, 1.45, spec.coarse_anis_ratio_list(:).'], 'stable');
s_dense_list = spec.coarse_s_dense_list;
s_sparse_list = spec.coarse_s_sparse_list;
expo_list = unique([2.0, 3.2, 5.0, 8.0], 'stable');
id = supplement_candidate_id_base_ms('c2_explore');
for isym = 1:numel(sym_list)
    for it = 1:numel(templates)
        for ia = 1:numel(anis_list)
            for idd = 1:numel(s_dense_list)
                for iss = 1:numel(s_sparse_list)
                    if s_sparse_list(iss) <= s_dense_list(idd) || (s_sparse_list(iss) / s_dense_list(idd)) < 1.45
                        continue;
                    end
                    for ie = 1:numel(expo_list)
                        if numel(jobs) >= max_jobs, return; end
                        id = id + 1;
                        job = pack_job_ms(id, 'subset_symmetric', templates{it}, sym_list{isym}, 'free', ...
                            s_dense_list(idd), s_sparse_list(iss), expo_list(ie), anis_list(ia), NaN, NaN);
                        job.search_group = 'c2_explore';
                        jobs(end+1) = job; %#ok<AGROW>
                    end
                end
            end
        end
    end
end
end

function jobs = build_shape_edge_anis_supplement_jobs_ms(spec, max_jobs)
jobs = empty_job_struct_ms();
N = spec.stage_count;
id = supplement_candidate_id_base_ms('shape_edge_anis');
shape_modes = {'ring_dense','near_center_dense','band_dense_x','band_dense_y','corner_dense','multi_center'};
stage_templates = build_stage13_shape_templates_ms(N, shape_modes);
stage_anis_pool = unique([0.70, 0.85, 1.00, 1.20, 1.45], 'stable');
s_ref = spec.coarse_s_dense_list(max(1, ceil(numel(spec.coarse_s_dense_list) / 2)));
s_sparse_ref = spec.coarse_s_sparse_list(max(1, ceil(numel(spec.coarse_s_sparse_list) / 2)));
expo_ref = spec.coarse_expo_list(max(1, ceil(numel(spec.coarse_expo_list) / 2)));
for it = 1:numel(stage_templates)
    shape_mode = first_shape_mode_in_template_ms(stage_templates{it});
    for ia = 1:numel(stage_anis_pool)
        if numel(jobs) >= max_jobs, return; end
        shape_cfg = struct('shape_mode', shape_mode, ...
            'stage_anis', make_stage_anis_vector_ms(stage_anis_pool(ia), N), ...
            'ring_radius_ratio', 0.35, 'ring_width_ratio', 0.16, ...
            'band_width_ratio', 0.28, 'corner_bias', 0.45, ...
            'jitter_ratio', 0.04, 'jitter_seed', 700 + id);
        id = id + 1;
        job = pack_job_ms(id, 'shape_explore', stage_templates{it}, 'c4', 'free', ...
            s_ref, s_sparse_ref, expo_ref, 1.0, NaN, NaN, NaN, NaN, ...
            repmat({'shape_explore'}, 1, N), shape_cfg);
        job.search_group = 'shape_edge_anis';
        jobs(end+1) = job; %#ok<AGROW>
    end
end

edge_templates = {
    supplement_template_ms(N, {'center_dense','center_dense','center_dense','edge_dense','edge_dense'}), ...
    supplement_template_ms(N, {'center_dense','edge_dense','center_dense','center_dense','center_dense'})};
edge_modes = {'edge_spaced'};
anis_list = unique([0.65, 0.75, 1.25, 1.55, spec.coarse_anis_ratio_list(:).'], 'stable');
for it = 1:numel(edge_templates)
    for ia = 1:numel(anis_list)
        for iedge = 1:numel(edge_modes)
            if numel(jobs) >= max_jobs, return; end
            id = id + 1;
            job = pack_job_ms(id, 'subset_symmetric', edge_templates{it}, 'c4', edge_modes{iedge}, ...
                s_ref, s_sparse_ref, expo_ref, anis_list(ia), NaN, NaN);
            job.search_group = 'shape_edge_anis';
            jobs(end+1) = job; %#ok<AGROW>
        end
    end
end
end

function templates = build_stage13_shape_templates_ms(N, shape_modes)
templates = {};
for im = 1:numel(shape_modes)
    for stage_idx = 1:min(3, N)
        t = repmat({'center_dense'}, 1, N);
        t{stage_idx} = shape_modes{im};
        if stage_idx < min(3, N)
            t{stage_idx + 1} = shape_modes{im};
        end
        templates{end+1} = t; %#ok<AGROW>
    end
end
end

function mode = first_shape_mode_in_template_ms(t)
mode = 'ring_dense';
for i = 1:numel(t)
    k = lower(strtrim(char(string(t{i}))));
    if any(strcmp(k, {'ring_dense','near_center_dense','center_quad_dense','corner_dense','band_dense_x','band_dense_y','multi_center'}))
        mode = k;
        return;
    end
end
end

function t = supplement_template_ms(N, raw)
t = repmat({'center_dense'}, 1, N);
for i = 1:N
    idx = min(i, numel(raw));
    t{i} = raw{idx};
end
end

function id = supplement_candidate_id_base_ms(group_name)
switch group_name
    case 'c4_main'
        id = 100000;
    case 'c2_explore'
        id = 200000;
    otherwise
        id = 300000;
end
end

function pick_mask = stratified_pick_candidates_ms(cands, cap)
n = numel(cands);
pick_mask = false(n, 1);
if n <= 0 || cap <= 0
    return;
end
if cap >= n
    pick_mask(:) = true;
    return;
end
keys = strings(n, 1);
for i = 1:n
    keys(i) = string(candidate_batch_stratum_key_ms(cands(i)));
end
u = unique(keys, 'stable');
rep_idx = zeros(numel(u), 1);
for i = 1:numel(u)
    idx = find(keys == u(i));
    if isempty(idx)
        continue;
    end
    rep_idx(i) = idx(randi(numel(idx), 1));
end
rep_idx = rep_idx(rep_idx > 0);
if isempty(rep_idx)
    rp = randperm(n, cap);
    pick_mask(rp) = true;
    return;
end
if numel(rep_idx) >= cap
    rp = randperm(numel(rep_idx), cap);
    pick_mask(rep_idx(rp)) = true;
    return;
end
pick_mask(rep_idx) = true;
picked = sum(pick_mask);
idx_rem = find(~pick_mask);
add_n = min(cap - picked, numel(idx_rem));
if add_n > 0
    rp = randperm(numel(idx_rem), add_n);
    pick_mask(idx_rem(rp)) = true;
end
end

function key = candidate_batch_stratum_key_ms(cand)
method = char(string(cand.layout_method));
shape_mode = '';
if isfield(cand, 'shape_mode') && ~isempty(cand.shape_mode)
    shape_mode = char(string(cand.shape_mode));
end
if isempty(shape_mode)
    shape_mode = 'none';
end
stage_method_sig = '';
if isfield(cand, 'stage_methods') && ~isempty(cand.stage_methods)
    stage_method_sig = strjoin(cellfun(@char, cand.stage_methods, 'UniformOutput', false), '/');
end
if isempty(stage_method_sig)
    stage_method_sig = method;
end
anis_bin = round(100 * scalar_with_default_ms(cand.anis_ratio, NaN));
expo_bin = round(100 * scalar_with_default_ms(cand.expo, NaN));
sd_bin = round(1e6 * scalar_with_default_ms(cand.s_dense, NaN));
ss_bin = round(1e6 * scalar_with_default_ms(cand.s_sparse, NaN));
key = sprintf('%s|%s|%s|a%d|e%d|sd%d|ss%d', ...
    method, shape_mode, stage_method_sig, anis_bin, expo_bin, sd_bin, ss_bin);
end

function s = summarize_job_method_counts_ms(jobs)
if isempty(jobs)
    s = 'none';
    return;
end
methods = string({jobs.layout_method});
u = unique(methods, 'stable');
parts = cell(1, numel(u));
for i = 1:numel(u)
    parts{i} = sprintf('%s=%d', char(u(i)), sum(methods == u(i)));
end
s = strjoin(parts, ',');
end

function cap = get_job_method_cap_ms(caps, method_name, fallback_n)
cap = fallback_n;
if nargin < 3 || ~isfinite(fallback_n) || fallback_n < 1
    fallback_n = inf;
end
if ~isstruct(caps)
    return;
end
k = lower(strtrim(char(method_name)));
if isfield(caps, k)
    v = caps.(k);
    if ~isempty(v) && isfinite(v) && v > 0
        cap = round(v);
    end
end
cap = min(cap, fallback_n);
end

function pick_mask = stratified_pick_jobs_ms(job_group, cap)
n = numel(job_group);
pick_mask = false(n, 1);
if n <= 0 || cap <= 0
    return;
end
if cap >= n
    pick_mask(:) = true;
    return;
end
keys = strings(n,1);
for i = 1:n
    keys(i) = string(job_budget_stratum_key_ms(job_group(i)));
end
u = unique(keys, 'stable');
rep_idx = zeros(numel(u), 1);
for i = 1:numel(u)
    idx = find(keys == u(i));
    if isempty(idx)
        continue;
    end
    rep_idx(i) = idx(randi(numel(idx), 1));
end
rep_idx = rep_idx(rep_idx > 0);
if isempty(rep_idx)
    rp = randperm(n, cap);
    pick_mask(rp) = true;
    return;
end

if numel(rep_idx) >= cap
    rp = randperm(numel(rep_idx), cap);
    pick_mask(rep_idx(rp)) = true;
    return;
end

pick_mask(rep_idx) = true;
picked = sum(pick_mask);
idx_rem = find(~pick_mask);
add_n = min(cap - picked, numel(idx_rem));
if add_n > 0
    rp = randperm(numel(idx_rem), add_n);
    pick_mask(idx_rem(rp)) = true;
end
end

function key = job_budget_stratum_key_ms(job)
stage_tag = strjoin(normalize_text_cell_ms(job.stage_modes, numel(job.stage_modes)), '_');
method_tag = '';
if isfield(job, 'stage_methods') && ~isempty(job.stage_methods)
    method_tag = strjoin(normalize_text_cell_ms(job.stage_methods, numel(job.stage_modes)), '_');
end
key = sprintf('%s|%s|%s|%s|m%s|a%s|g%s|p%s|w%s|sh%s|sa%s|rr%s|jr%s', ...
    lower(strtrim(char(job.layout_method))), ...
    lower(strtrim(char(job.symmetry_mode))), ...
    lower(strtrim(char(job.edge_pattern_mode))), ...
    stage_tag, ...
    method_tag, ...
    numeric_tag_ms(job.method_anchor_stage), ...
    numeric_tag_ms(job.gamma), ...
    numeric_tag_ms(job.poisson_seed), ...
    numeric_tag_ms(job.alignment_weight), ...
    lower(strtrim(char(string(job.shape_mode)))), ...
    numeric_tag_ms(first_numeric_ms(job.stage_anis, NaN)), ...
    numeric_tag_ms(job.ring_radius_ratio), ...
    numeric_tag_ms(job.jitter_ratio));
end

function templates = build_mode_templates_ms(N, fallback_trends)
templates = {};
templates{end+1} = repmat({'center_dense'}, 1, N); %#ok<AGROW>
templates{end+1} = repmat({'edge_dense'}, 1, N); %#ok<AGROW>

md = repmat({'center_dense'}, 1, N);
for k = max(1, ceil(N/2)):N
    md{k} = 'edge_dense';
end
templates{end+1} = md; %#ok<AGROW>

md = repmat({'edge_dense'}, 1, N);
for k = 1:max(1, floor(N/2))
    md{k} = 'center_dense';
end
templates{end+1} = md; %#ok<AGROW>

fb = normalize_trend_cell_ms(fallback_trends, N);
md = repmat({'center_dense'}, 1, N);
for k = 1:N
    if strcmpi(fb{k}, 'edge_heavy')
        md{k} = 'edge_dense';
    else
        md{k} = 'center_dense';
    end
end
templates{end+1} = md; %#ok<AGROW>
end

function picked_jobs = sample_warmup_jobs_ms(jobs, spec)
picked_jobs = empty_job_struct_ms();
N = numel(jobs);
if N == 0
    return;
end
ratio = spec.warmup.sample_ratio;
target_total = min(N, max(round(ratio * N), spec.warmup.min_sample_jobs));
keys = cell(N,1);
for i = 1:N
    keys{i} = warmup_job_key_ms(jobs(i));
end
uk = unique(keys, 'stable');
selected = false(N,1);
for i = 1:numel(uk)
    idx = find(strcmp(keys, uk{i}));
    take = min(numel(idx), max(spec.warmup.min_per_stratum, round(ratio * numel(idx))));
    rp = randperm(numel(idx), take);
    selected(idx(rp)) = true;
end
picked_count = sum(selected);
if picked_count < target_total
    remain = find(~selected);
    add_n = min(numel(remain), target_total - picked_count);
    if add_n > 0
        rp = randperm(numel(remain), add_n);
        selected(remain(rp)) = true;
    end
end
picked_jobs = jobs(selected);
end

function k = warmup_job_key_ms(job)
stage_mode_tag = strjoin(job.stage_modes, '_');
stage_method_tag = '';
if isfield(job, 'stage_methods') && ~isempty(job.stage_methods)
    stage_method_tag = strjoin(normalize_text_cell_ms(job.stage_methods, numel(job.stage_modes)), '_');
end
k = sprintf('%s|%s|%s|%s|%s|%s', job.layout_method, stage_mode_tag, stage_method_tag, ...
    job.symmetry_mode, job.edge_pattern_mode, warmup_expo_bucket_ms(job.expo));
end

function b = warmup_expo_bucket_ms(expo)
if ~isfinite(expo)
    b = 'nan';
elseif expo <= 2.0
    b = 'lo';
elseif expo <= 6.0
    b = 'mid';
else
    b = 'hi';
end
end

function s = empty_job_struct_ms()
s = struct('candidate_id', {}, 'layout_method', {}, ...
    'search_group', {}, ...
    'stage_modes', {}, 'stage_methods', {}, 'stage_trends', {}, ...
    'symmetry_mode', {}, 'edge_pattern_mode', {}, ...
    's_dense', {}, 's_sparse', {}, 'expo', {}, 'anis_ratio', {}, ...
    'gamma', {}, 'method_anchor_stage', {}, ...
    'poisson_seed', {}, 'alignment_weight', {}, ...
    'shape_mode', {}, 'stage_anis', {}, 'ring_radius_ratio', {}, ...
    'ring_width_ratio', {}, 'band_width_ratio', {}, 'corner_bias', {}, ...
    'jitter_ratio', {}, 'jitter_seed', {});
end

function job = pack_job_ms(id, layout_method, stage_modes, symmetry_mode, edge_pattern_mode, ...
    s_dense, s_sparse, expo, anis, gamma, method_anchor_stage, poisson_seed, alignment_weight, stage_methods, shape_cfg)
if nargin < 10 || isempty(gamma)
    gamma = NaN;
end
if nargin < 11 || isempty(method_anchor_stage)
    method_anchor_stage = NaN;
end
if nargin < 12 || isempty(poisson_seed)
    poisson_seed = NaN;
end
if nargin < 13 || isempty(alignment_weight)
    alignment_weight = NaN;
end
if nargin < 14 || isempty(stage_methods)
    stage_methods = repmat({char(layout_method)}, 1, numel(stage_modes));
end
if nargin < 15 || ~isstruct(shape_cfg)
    shape_cfg = struct();
end
if isstruct(stage_methods)
    shape_cfg = stage_methods;
    stage_methods = repmat({char(layout_method)}, 1, numel(stage_modes));
end
stage_methods = normalize_text_cell_ms(stage_methods, numel(stage_modes));
stage_trends = cell(size(stage_modes));
for i = 1:numel(stage_modes)
    stage_trends{i} = mode_to_trend_ms(stage_modes{i}, gamma);
end
job = struct('candidate_id', id, 'layout_method', char(layout_method), ...
    'search_group', 'regular', ...
    'stage_modes', {stage_modes}, 'stage_methods', {stage_methods}, 'stage_trends', {stage_trends}, ...
    'symmetry_mode', char(symmetry_mode), 'edge_pattern_mode', char(edge_pattern_mode), ...
    's_dense', s_dense, 's_sparse', s_sparse, 'expo', expo, 'anis_ratio', anis, ...
    'gamma', gamma, 'method_anchor_stage', method_anchor_stage, ...
    'poisson_seed', poisson_seed, 'alignment_weight', alignment_weight, ...
    'shape_mode', shape_cfg_value_ms(shape_cfg, 'shape_mode', ''), ...
    'stage_anis', shape_cfg_value_ms(shape_cfg, 'stage_anis', NaN), ...
    'ring_radius_ratio', shape_cfg_value_ms(shape_cfg, 'ring_radius_ratio', NaN), ...
    'ring_width_ratio', shape_cfg_value_ms(shape_cfg, 'ring_width_ratio', NaN), ...
    'band_width_ratio', shape_cfg_value_ms(shape_cfg, 'band_width_ratio', NaN), ...
    'corner_bias', shape_cfg_value_ms(shape_cfg, 'corner_bias', NaN), ...
    'jitter_ratio', shape_cfg_value_ms(shape_cfg, 'jitter_ratio', NaN), ...
    'jitter_seed', shape_cfg_value_ms(shape_cfg, 'jitter_seed', NaN));
end

function v = shape_cfg_value_ms(cfg, fname, default_v)
v = default_v;
if isstruct(cfg) && isfield(cfg, fname)
    v = cfg.(fname);
end
end

function trend = mode_to_trend_ms(mode_name, gamma)
if nargin < 2
    gamma = NaN;
end
mode_key = lower(strtrim(char(string(mode_name))));
switch mode_key
    case 'center_dense'
        trend = 'center_heavy';
    case 'edge_dense'
        trend = 'edge_heavy';
    case {'ring_dense','near_center_dense','center_quad_dense','corner_dense','band_dense_x','band_dense_y','multi_center','neutral'}
        trend = 'neutral';
    case 'gamma'
        if isfinite(gamma) && gamma < 0
            trend = 'edge_heavy';
        elseif isfinite(gamma) && gamma > 0
            trend = 'center_heavy';
        else
            trend = 'neutral';
        end
    otherwise
        trend = 'neutral';
end
end

function cands = instantiate_candidate_layouts_ms(jobs, spec, count_solution)
N = spec.stage_count;
if isempty(jobs)
    cands = empty_candidate_struct_ms(N);
    return;
end
cells = cell(numel(jobs), 1);
use_par = spec.use_parallel && numel(jobs) > 1;
if use_par
    parfor i = 1:numel(jobs)
        [ok, rec] = instantiate_candidate_job_ms(jobs(i), spec, count_solution);
        if ok
            cells{i} = rec;
        end
    end
else
    for i = 1:numel(jobs)
        [ok, rec] = instantiate_candidate_job_ms(jobs(i), spec, count_solution);
        if ok
            cells{i} = rec;
        end
    end
end
keep = ~cellfun(@isempty, cells);
if ~any(keep)
    cands = empty_candidate_struct_ms(N);
    return;
end
idx = find(keep);
tmpl = candidate_record_template_ms(N);
cands = repmat(tmpl, numel(idx), 1);
for i = 1:numel(idx)
    cands(i) = cells{idx(i)};
end
end

function [ok, rec] = instantiate_candidate_job_ms(job, spec, count_solution)
N = spec.stage_count;
ok = false;
rec = candidate_record_template_ms(N);
hc = resolve_hard_constraints_ms(spec);
if isfield(job, 'symmetry_mode') && is_c2_symmetry_mode_ms(job.symmetry_mode)
    hc.force_c4_only = false;
end

stages = cell(1, N);
for k = 1:N
    method_k = stage_method_for_job_ms(job, k, N);
    mode_k = job.stage_modes{k};
    anis_k = job.anis_ratio;
    if isfield(job, 'stage_anis') && isnumeric(job.stage_anis) && numel(job.stage_anis) >= k && isfinite(job.stage_anis(k))
        anis_k = job.stage_anis(k);
    end
    shape_params_k = job_shape_params_for_stage_ms(job, k, N);
    prev_rects = zeros(0,4);
    prev2_rects = zeros(0,4);
    if k > 1 && isstruct(stages{k-1}) && isfield(stages{k-1}, 'rects')
        prev_rects = stages{k-1}.rects;
    end
    if k > 2 && isstruct(stages{k-2}) && isfield(stages{k-2}, 'rects')
        prev2_rects = stages{k-2}.rects;
    end
    [ok_k, stage_k] = build_stage_layout_ms( ...
        count_solution.n(k), mode_k, job.s_dense, job.s_sparse, ...
        job.expo, anis_k, spec.geometry.L_max(k), ...
        spec.geometry.coverage_min(k), method_k, job.symmetry_mode, ...
        job.edge_pattern_mode, spec, job.gamma, job.poisson_seed, ...
        job.alignment_weight, prev_rects, prev2_rects, k, shape_params_k);
    if ~ok_k
        return;
    end
    stages{k} = stage_k;
end

Lx = zeros(1, N);
Ly = zeros(1, N);
Lx(N) = stages{N}.Lx;
Ly(N) = stages{N}.Ly;
for k = N-1:-1:1
    Lx(k) = max(stages{k}.Lx, Lx(k+1) + spec.geometry.pyramid_gap_min(k));
    Ly(k) = max(stages{k}.Ly, Ly(k+1) + spec.geometry.pyramid_gap_min(k));
end

cov = zeros(1, N);
for k = 1:N
    cov(k) = count_solution.n(k) * spec.fp_w * spec.fp_h / max(Lx(k) * Ly(k), eps);
end
[ok_hc, ~] = validate_candidate_hard_constraints_ms(job.symmetry_mode, stages, Lx, Ly, cov, spec, hc);
if ~ok_hc
    return;
end

rec.candidate_id = job.candidate_id;
rec.layout_method = job.layout_method;
if isfield(job, 'search_group')
    rec.search_group = job.search_group;
end
rec.stage_modes = job.stage_modes;
rec.stage_methods = job.stage_methods;
rec.stage_trends = job.stage_trends;
rec.symmetry_mode = job.symmetry_mode;
rec.edge_pattern_mode = job.edge_pattern_mode;
rec.s_dense = job.s_dense;
rec.s_sparse = job.s_sparse;
rec.expo = job.expo;
rec.anis_ratio = job.anis_ratio;
rec.gamma = job.gamma;
rec.method_anchor_stage = job.method_anchor_stage;
rec.shape_mode = job.shape_mode;
rec.stage_anis = fit_len_vec_ms(job.stage_anis, N, NaN(1, N));
rec.ring_radius_ratio = job.ring_radius_ratio;
rec.ring_width_ratio = job.ring_width_ratio;
rec.band_width_ratio = job.band_width_ratio;
rec.corner_bias = job.corner_bias;
rec.jitter_ratio = job.jitter_ratio;
rec.jitter_seed = job.jitter_seed;
rec.I_opt = count_solution.I_opt;
if isfield(count_solution, 'count_rank')
    rec.count_rank = count_solution.count_rank;
else
    rec.count_rank = NaN;
end
rec.n = count_solution.n;
rec.ratios = count_solution.ratios;
rec.Lx = Lx;
rec.Ly = Ly;
rec.cov = cov;
rec.stage_rects = cell(1, N);
for k = 1:N
    rec.stage_rects{k} = stages{k}.rects;
end
rec.lmax_relaxed = any(cellfun(@(s) logical(s.lmax_relaxed), stages));
rec.lmax_relax_ratio = max([1.0, cellfun(@(s) s.lmax_relax_ratio, stages)]);
rec.spacing_ratio = rec.s_sparse / max(rec.s_dense, eps);
rec.contrast_score = calc_contrast_score_ms(rec.s_dense, rec.s_sparse, rec.expo);
ok = true;
end

function method_k = stage_method_for_job_ms(job, stage_idx, N)
method_k = lower(strtrim(char(job.layout_method)));
if isfield(job, 'stage_methods') && numel(job.stage_methods) >= stage_idx
    mk = lower(strtrim(char(job.stage_methods{stage_idx})));
    if ~isempty(mk)
        method_k = mk;
    end
end
if strcmp(method_k, 'mixed_stage_methods')
    method_k = 'subset_symmetric';
end
end

function shape_params = job_shape_params_for_stage_ms(job, stage_idx, N)
if nargin < 3 || ~isfinite(N)
    N = numel(job.stage_modes);
end
shape_params = struct('shape_mode', '', 'stage_anis', NaN(1, N), ...
    'ring_radius_ratio', NaN, 'ring_width_ratio', NaN, ...
    'band_width_ratio', NaN, 'corner_bias', NaN, ...
    'jitter_ratio', NaN, 'jitter_seed', NaN, 'stage_idx', stage_idx);
fields = {'shape_mode','stage_anis','ring_radius_ratio','ring_width_ratio', ...
    'band_width_ratio','corner_bias','jitter_ratio','jitter_seed'};
for i = 1:numel(fields)
    f = fields{i};
    if isfield(job, f)
        shape_params.(f) = job.(f);
    end
end
if isfield(job, 'layout_method') && strcmpi(job.layout_method, 'shape_explore') && ...
        isfield(job, 'stage_modes') && numel(job.stage_modes) >= stage_idx
    shape_params.shape_mode = lower(strtrim(char(string(job.stage_modes{stage_idx}))));
end
end

function cfg = normalize_hard_constraints_ms(raw, fallback)
cfg = fallback;
if nargin < 2 || ~isstruct(fallback)
    cfg = struct('force_c4_only', true, 'unified_lmax_across_stages', true, ...
        'allow_lmax_relax', false, ...
        'enforce_monotonic_geometry', true, 'enforce_coverage_min', true, ...
        'enforce_spacing_min_edge_gap', true);
end
if nargin < 1 || ~isstruct(raw)
    return;
end
fn = fieldnames(cfg);
for i = 1:numel(fn)
    k = fn{i};
    if isfield(raw, k)
        cfg.(k) = any(logical(raw.(k)));
    end
end
end

function s = empty_candidate_struct_ms(N)
if nargin < 1
    N = 5;
end
s = struct('candidate_id', {}, 'layout_method', {}, ...
    'search_group', {}, ...
    'stage_modes', {}, 'stage_methods', {}, 'stage_trends', {}, ...
    'symmetry_mode', {}, 'edge_pattern_mode', {}, ...
    's_dense', {}, 's_sparse', {}, 'expo', {}, 'anis_ratio', {}, 'gamma', {}, ...
    'method_anchor_stage', {}, ...
    'spacing_ratio', {}, 'contrast_score', {}, 'mode_prior_score', {}, ...
    'lmax_relaxed', {}, 'lmax_relax_ratio', {}, ...
    'I_opt', {}, 'count_rank', {}, ...
    'n', {}, 'ratios', {}, ...
    'Lx', {}, 'Ly', {}, 'cov', {}, ...
    'shape_mode', {}, 'stage_anis', {}, 'ring_radius_ratio', {}, ...
    'ring_width_ratio', {}, 'band_width_ratio', {}, 'corner_bias', {}, ...
    'jitter_ratio', {}, 'jitter_seed', {}, ...
    'stage_rects', {});
if N < 2 %#ok<NASGU>
end
end

function rec = candidate_record_template_ms(N)
if nargin < 1
    N = 5;
end
rec = struct('candidate_id', -1, 'layout_method', '', ...
    'search_group', 'regular', ...
    'stage_modes', {repmat({''}, 1, N)}, ...
    'stage_methods', {repmat({''}, 1, N)}, ...
    'stage_trends', {repmat({'neutral'}, 1, N)}, ...
    'symmetry_mode', '', 'edge_pattern_mode', '', ...
    's_dense', NaN, 's_sparse', NaN, 'expo', NaN, 'anis_ratio', NaN, 'gamma', NaN, ...
    'method_anchor_stage', NaN, ...
    'spacing_ratio', NaN, 'contrast_score', NaN, 'mode_prior_score', NaN, ...
    'lmax_relaxed', false, 'lmax_relax_ratio', 1.0, ...
    'I_opt', NaN, 'count_rank', NaN, ...
    'n', NaN(1, N), 'ratios', NaN(1, N-1), ...
    'Lx', NaN(1, N), 'Ly', NaN(1, N), 'cov', NaN(1, N), ...
    'shape_mode', '', 'stage_anis', NaN(1, N), 'ring_radius_ratio', NaN, ...
    'ring_width_ratio', NaN, 'band_width_ratio', NaN, 'corner_bias', NaN, ...
    'jitter_ratio', NaN, 'jitter_seed', NaN, ...
    'stage_rects', {cell(1, N)});
end

function [cands_out, info] = deduplicate_candidates_by_geometry_ms(cands, spec, tag)
if nargin < 3 || isempty(tag)
    tag = 'Dedup';
end
dcfg = get_candidate_dedup_config_ms(spec);
info = struct('enabled', dcfg.enable, 'tol_m', dcfg.tol_m, ...
    'count_in', numel(cands), 'count_out', numel(cands), 'removed', 0, ...
    'removed_ratio', 0);
cands_out = cands;
if isempty(cands)
    return;
end
if ~dcfg.enable
    fprintf('[%s-Dedup] disabled, candidates=%d\n', tag, numel(cands));
    return;
end
N = numel(cands);
keys = strings(N,1);
for i = 1:N
    keys(i) = string(make_candidate_geometry_key_ms(cands(i), dcfg.tol_m));
end
[~, keep_idx] = unique(keys, 'stable');
if numel(keep_idx) < N
    keep_mask = false(N,1);
    keep_mask(keep_idx) = true;
    cands_out = cands(keep_mask);
end
info.count_out = numel(cands_out);
info.removed = info.count_in - info.count_out;
info.removed_ratio = info.removed / max(info.count_in, 1);
fprintf('[%s-Dedup] tol=%.3g m, in=%d, out=%d, removed=%d (ratio=%.3f)\n', ...
    tag, dcfg.tol_m, info.count_in, info.count_out, info.removed, info.removed_ratio);
end

function dcfg = get_candidate_dedup_config_ms(spec)
dcfg = struct('enable', true, 'tol_m', 1e-9);
if nargin < 1 || ~isstruct(spec) || ~isfield(spec, 'candidate_dedup')
    return;
end
cfg = spec.candidate_dedup;
if isfield(cfg, 'enable')
    dcfg.enable = logical(cfg.enable);
end
if isfield(cfg, 'tol_m')
    dcfg.tol_m = max(cfg.tol_m, 1e-12);
end
end

function key = make_candidate_geometry_key_ms(cand, tol_m)
if nargin < 2 || ~isfinite(tol_m) || tol_m <= 0
    tol_m = 1e-9;
end
if ~is_candidate_geometry_finite_ms(cand)
    key = 'invalid';
    return;
end
N = numel(cand.n);
Lq = int64(round([cand.Lx(:); cand.Ly(:)] / tol_m));
h = zeros(1, N);
for k = 1:N
    h(k) = hash_rects_geom_ms(cand.stage_rects{k}, tol_m);
end
hhex = cell(1, N);
for k = 1:N
    hhex{k} = dec2hex(h(k), 16);
end
key = sprintf('%s|%s|%s', ...
    strtrim(sprintf('%d ', round(cand.n))), ...
    strtrim(sprintf('%d ', Lq)), ...
    strjoin(hhex, '_'));
end

function tf = is_candidate_geometry_finite_ms(cand)
tf = all(isfinite(cand.Lx)) && all(isfinite(cand.Ly));
if ~tf
    return;
end
for k = 1:numel(cand.stage_rects)
    r = cand.stage_rects{k};
    if size(r,2) ~= 4 || ~all(isfinite(r(:)))
        tf = false;
        return;
    end
end
end

function h = hash_rects_geom_ms(rects, tol_m)
offset = uint64(1469598103934665603);
prime = uint64(1099511628211);
h = offset;
if isempty(rects)
    return;
end
rq = int64(round(rects / tol_m));
rq = sortrows(rq, [1,2,3,4]);
u = typecast(rq(:), 'uint64');
for i = 1:numel(u)
    h = bitxor(h, u(i));
    h = h * prime;
end
h = bitxor(h, uint64(size(rq,1)));
end

function [ok, stage] = build_stage_layout_ms(n_target, mode, s_dense, s_sparse, expo, anis, ...
    Lmax, coverage_min_stage, layout_method, symmetry_mode, edge_pattern_mode, ...
    spec, gamma_val, poisson_seed, alignment_weight, prev_rects, prev2_rects, stage_idx, shape_params)
if nargin < 19 || ~isstruct(shape_params)
    shape_params = struct();
end
[ok, stage, fail_code] = build_stage_layout_core_ms(n_target, mode, s_dense, s_sparse, expo, anis, ...
    Lmax, coverage_min_stage, layout_method, symmetry_mode, edge_pattern_mode, ...
    spec, gamma_val, poisson_seed, alignment_weight, prev_rects, prev2_rects, stage_idx, shape_params);
if ok
    stage.lmax_relaxed = false;
    stage.lmax_relax_ratio = 1.0;
    return;
end
if ~(strcmp(fail_code, 'capacity') || strcmp(fail_code, 'edge_pattern'))
    return;
end
hc = resolve_hard_constraints_ms(spec);
if ~isfield(hc, 'allow_lmax_relax') || ~logical(hc.allow_lmax_relax)
    return;
end
Lmax_soft = 1.05 * Lmax;
[ok_soft, stage_soft] = build_stage_layout_core_ms(n_target, mode, s_dense, s_sparse, expo, anis, ...
    Lmax_soft, coverage_min_stage, layout_method, symmetry_mode, edge_pattern_mode, ...
    spec, gamma_val, poisson_seed, alignment_weight, prev_rects, prev2_rects, stage_idx, shape_params);
if ~ok_soft
    ok = false;
    stage = struct();
    return;
end
stage = stage_soft;
stage.lmax_relaxed = true;
stage.lmax_relax_ratio = Lmax_soft / max(Lmax, eps);
ok = true;
end

function [ok, stage, fail_code] = build_stage_layout_core_ms(n_target, mode, s_dense, s_sparse, expo, anis, ...
    Lmax_use, coverage_min_stage, layout_method, symmetry_mode, edge_pattern_mode, ...
    spec, gamma_val, poisson_seed, alignment_weight, prev_rects, prev2_rects, stage_idx, shape_params)
ok = false;
stage = struct();
fail_code = 'unknown';
if nargin < 19 || ~isstruct(shape_params)
    shape_params = struct();
end
if n_target < 1
    fail_code = 'capacity';
    return;
end

if strcmpi(layout_method, 'shape_explore')
    [ok_u, pts] = make_shape_explore_points_ms(n_target, Lmax_use, spec.fp_w, spec.fp_h, ...
        spec.edge_margin_fp, s_dense, anis, spec.geometry.min_edge_gap, mode, shape_params);
    if ~ok_u
        fail_code = 'capacity';
        return;
    end
elseif strcmpi(layout_method, 'hex_c6')
    [ok_u, pts] = make_hex_c6_points_ms(n_target, Lmax_use, spec.fp_w, spec.fp_h, ...
        spec.edge_margin_fp, mode, s_dense, anis, spec.geometry.min_edge_gap);
    if ~ok_u
        fail_code = 'capacity';
        return;
    end
elseif strcmpi(layout_method, 'gamma_stage2')
    [ok_u, pts] = make_radial_gamma_points_ms(n_target, Lmax_use, spec.fp_w, spec.fp_h, ...
        spec.edge_margin_fp, gamma_val, s_dense, anis, spec.geometry.min_edge_gap);
    if ~ok_u
        fail_code = 'capacity';
        return;
    end
elseif strcmpi(layout_method, 'ring_stage3')
    [ok_u, pts] = make_concentric_ring_points_ms(n_target, Lmax_use, spec.fp_w, spec.fp_h, ...
        spec.edge_margin_fp, mode, spec.geometry.min_edge_gap);
    if ~ok_u
        fail_code = 'capacity';
        return;
    end
elseif strcmpi(layout_method, 'poisson_disk')
    seed_use = poisson_seed;
    if ~isfinite(seed_use)
        seed_use = 0;
    end
    seed_use = round(seed_use + 1000 * stage_idx + n_target);
    [ok_u, pts] = make_poisson_disk_points_ms(n_target, Lmax_use, spec.fp_w, spec.fp_h, ...
        spec.edge_margin_fp, spec.geometry.min_edge_gap, seed_use);
    if ~ok_u
        fail_code = 'capacity';
        return;
    end
elseif strcmpi(layout_method, 'interstage_aligned')
    [ok_u, pts] = make_interstage_aligned_points_ms(n_target, Lmax_use, spec.fp_w, spec.fp_h, ...
        spec.edge_margin_fp, spec.geometry.min_edge_gap, prev_rects, prev2_rects, alignment_weight);
    if ~ok_u
        fail_code = 'capacity';
        return;
    end
elseif strcmpi(layout_method, 'fixed_c4_grid')
    [ok_u, pts] = make_regular_fullgrid_points_ms(n_target, Lmax_use, spec.fp_w, spec.fp_h, ...
        spec.edge_margin_fp, s_dense, anis, spec.geometry.min_edge_gap);
    if ~ok_u
        fail_code = 'capacity';
        return;
    end
else
    fp = max(spec.fp_w, spec.fp_h);
    min_step_x = max([spec.min_spacing_manufacture, spec.fp_w + spec.geometry.min_edge_gap, spec.fp_w * 1.01]);
    min_step_y = max([spec.min_spacing_manufacture, spec.fp_h + spec.geometry.min_edge_gap, spec.fp_h * 1.01]);
    sxd = max(s_dense * anis, min_step_x);
    sxs = max(s_sparse * anis, sxd + fp * 0.01);
    syd = max(s_dense / anis, min_step_y);
    sys = max(s_sparse / anis, syd + fp * 0.01);

    if count_points_xy_ms(Lmax_use, mode, sxd, sxs, syd, sys, expo, spec.fp_w, spec.fp_h, spec.edge_margin_fp) < n_target
        fail_code = 'capacity';
        return;
    end

    Llo = max(2.2e-3, fp * 1.2);
    Lhi = Lmax_use;
    for it = 1:24
        Lmid = 0.5 * (Llo + Lhi);
        if count_points_xy_ms(Lmid, mode, sxd, sxs, syd, sys, expo, spec.fp_w, spec.fp_h, spec.edge_margin_fp) >= n_target
            Lhi = Lmid;
        else
            Llo = Lmid;
        end
    end

    ok_sel = false;
    pts = [];
    for Lt = linspace(Lhi, Lmax_use, 8)
        rects = make_gradient_footprints_xy_ms(sxd, sxs, syd, sys, expo, Lt, spec.fp_w, spec.fp_h, mode, spec.edge_margin_fp);
        p = rect_centers_ms(rects);
        [ok_sel, pts] = select_exact_count_symmetric_ms(p, n_target, mode, symmetry_mode, edge_pattern_mode);
        if ok_sel
            break;
        end
    end
    if ~ok_sel
        [ok_sel, pts] = force_feasible_layout_ms(n_target, mode, symmetry_mode, edge_pattern_mode, ...
            Lmax_use, spec.fp_w, spec.fp_h, spec.edge_margin_fp, spec.geometry.min_edge_gap);
    end
    if ~ok_sel
        fail_code = 'edge_pattern';
        return;
    end
end

rects = centers_to_rects_ms(pts, spec.fp_w, spec.fp_h);
[spacing_ok, ~] = validate_stage_rect_spacing_fast_ms(rects, spec.fp_w, spec.fp_h, spec.geometry.min_edge_gap);
if ~spacing_ok
    fail_code = 'spacing';
    return;
end
[Lx, Ly] = infer_plate_size_xy_ms(pts, spec.fp_w, spec.fp_h, spec.edge_margin_fp);
if max(Lx, Ly) > Lmax_use
    fail_code = 'capacity';
    return;
end
cov = n_target * spec.fp_w * spec.fp_h / max(Lx * Ly, eps);
if cov < coverage_min_stage
    fail_code = 'coverage';
    return;
end

stage.rects = rects;
stage.Lx = Lx;
stage.Ly = Ly;
stage.coverage = cov;
fail_code = 'ok';
ok = true;
end

function n = count_points_xy_ms(span, mode, sxd, sxs, syd, sys, expo, fp_w, fp_h, edge_margin_fp)
rects = make_gradient_footprints_xy_ms(sxd, sxs, syd, sys, expo, span, fp_w, fp_h, mode, edge_margin_fp);
n = size(rects, 1);
end

function rects = make_gradient_footprints_xy_ms(sxd, sxs, syd, sys, expo, span, fp_w, fp_h, mode, edge_margin_fp)
fp = max(fp_w, fp_h);
half = span/2 - fp/2 - edge_margin_fp * fp;
if half <= 0
    rects = zeros(0,4);
    return;
end
switch lower(mode)
    case 'center_dense'
        sfx = @(u) sxd + (sxs-sxd) * min(abs(u)./half,1).^expo;
        sfy = @(u) syd + (sys-syd) * min(abs(u)./half,1).^expo;
    case 'edge_dense'
        sfx = @(u) sxs - (sxs-sxd) * min(abs(u)./half,1).^expo;
        sfy = @(u) sys - (sys-syd) * min(abs(u)./half,1).^expo;
    otherwise
        error('Unknown mode: %s', mode);
end
px = generate_half_axis_ms(sfx, half);
py = generate_half_axis_ms(sfy, half);
px = px(:).';
py = py(:).';
if numel(px) > 1
    fx = [-px(end:-1:2), px];
else
    fx = px;
end
if numel(py) > 1
    fy = [-py(end:-1:2), py];
else
    fy = py;
end
[X, Y] = meshgrid(fx, fy);
rects = [X(:)-fp_w/2, X(:)+fp_w/2, Y(:)-fp_h/2, Y(:)+fp_h/2];
end

function pos = generate_half_axis_ms(sfn, half)
pos = 0;
x = 0;
while true
    s = sfn(x);
    x_next = x + s;
    if x_next > half + s*1e-9
        break;
    end
    pos(end+1) = x_next; %#ok<AGROW>
    x = x_next;
end
pos = pos(pos <= half + sfn(half)*1e-9);
end

function [ok, pts] = make_shape_explore_points_ms(n_target, Lmax, fp_w, fp_h, edge_margin_fp, ...
    s_dense, anis_ratio, min_edge_gap, mode, shape_params)
ok = false;
pts = zeros(0,2);
if n_target < 1 || Lmax <= 0
    return;
end
fp = max(fp_w, fp_h);
margin = fp/2 + edge_margin_fp * fp;
halfx = Lmax/2 - margin;
halfy = Lmax/2 - margin;
if halfx <= 0 || halfy <= 0
    return;
end
min_dx = max(fp_w + min_edge_gap, fp_w * 1.01);
min_dy = max(fp_h + min_edge_gap, fp_h * 1.01);
sx = max(s_dense * max(anis_ratio, 0.1), min_dx);
sy = max(s_dense / max(anis_ratio, 0.1), min_dy);
nx = max(1, floor(2 * halfx / sx) + 1);
ny = max(1, floor(2 * halfy / sy) + 1);
if nx * ny < n_target
    sx = min_dx;
    sy = min_dy;
    nx = max(1, floor(2 * halfx / sx) + 1);
    ny = max(1, floor(2 * halfy / sy) + 1);
end
if nx * ny < n_target
    return;
end
xvals = ((0:nx-1) - (nx-1)/2) * sx;
yvals = ((0:ny-1) - (ny-1)/2) * sy;
xvals = xvals(abs(xvals) <= halfx + 1e-12);
yvals = yvals(abs(yvals) <= halfy + 1e-12);
[X, Y] = meshgrid(xvals, yvals);
pool = [X(:), Y(:)];
if size(pool,1) < n_target
    return;
end
shape_mode = lower(strtrim(char(string(mode))));
if isfield(shape_params, 'shape_mode') && ~isempty(shape_params.shape_mode)
    shape_mode = lower(strtrim(char(string(shape_params.shape_mode))));
end
scores = score_shape_pool_ms(pool, shape_mode, halfx, halfy, shape_params);
[ok_sel, pts] = select_exact_count_symmetric_weighted_ms(pool, n_target, scores, 'c4');
if ~ok_sel
    [ok_sel, pts] = select_exact_count_symmetric_ms(pool, n_target, 'center_dense', 'c4', 'free');
end
if ~ok_sel
    return;
end
pts = apply_shape_jitter_ms(pts, shape_params, fp_w, fp_h, min_edge_gap, halfx, halfy);
rects = centers_to_rects_ms(pts, fp_w, fp_h);
[spacing_ok, ~] = validate_stage_rect_spacing_fast_ms(rects, fp_w, fp_h, min_edge_gap);
if ~spacing_ok
    return;
end
ok = true;
end

function scores = score_shape_pool_ms(pool, shape_mode, halfx, halfy, shape_params)
x = pool(:,1);
y = pool(:,2);
rx = abs(x) / max(halfx, eps);
ry = abs(y) / max(halfy, eps);
r = sqrt((x / max(halfx, eps)).^2 + (y / max(halfy, eps)).^2) / sqrt(2);
scores = ones(size(x));
switch lower(shape_mode)
    case 'center_dense'
        scores = 1 ./ (1 + r);
    case 'edge_dense'
        scores = r + 1e-3;
    case 'neutral'
        scores = ones(size(x));
    case 'ring_dense'
        rr = bounded_shape_param_ms(shape_params, 'ring_radius_ratio', 0.45, 0.05, 0.95);
        rw = bounded_shape_param_ms(shape_params, 'ring_width_ratio', 0.18, 0.03, 0.80);
        scores = exp(-((r - rr) / max(rw, eps)).^2);
    case {'near_center_dense','center_quad_dense'}
        rr = bounded_shape_param_ms(shape_params, 'ring_radius_ratio', 0.18, 0.03, 0.45);
        rw = bounded_shape_param_ms(shape_params, 'ring_width_ratio', 0.10, 0.03, 0.35);
        center_penalty = exp(-(r / max(0.035, rw * 0.35)).^2);
        scores = exp(-((r - rr) / max(rw, eps)).^2) .* (1 - center_penalty);
    case 'corner_dense'
        cb = bounded_shape_param_ms(shape_params, 'corner_bias', 0.45, 0.01, 2.0);
        scores = (rx .* ry + 1e-3) .^ cb;
    case 'band_dense_x'
        bw = bounded_shape_param_ms(shape_params, 'band_width_ratio', 0.30, 0.03, 0.90);
        scores = exp(-(ry / max(bw, eps)).^2);
    case 'band_dense_y'
        bw = bounded_shape_param_ms(shape_params, 'band_width_ratio', 0.30, 0.03, 0.90);
        scores = exp(-(rx / max(bw, eps)).^2);
    case 'multi_center'
        rr = bounded_shape_param_ms(shape_params, 'ring_radius_ratio', 0.45, 0.05, 0.95);
        rw = bounded_shape_param_ms(shape_params, 'ring_width_ratio', 0.18, 0.03, 0.80);
        centers = [0 0; rr 0; -rr 0; 0 rr; 0 -rr];
        scores = zeros(size(x));
        xn = x / max(halfx, eps);
        yn = y / max(halfy, eps);
        for i = 1:size(centers,1)
            d2 = (xn - centers(i,1)).^2 + (yn - centers(i,2)).^2;
            scores = max(scores, exp(-d2 / max(rw^2, eps)));
        end
    otherwise
        scores = 1 ./ (1 + r);
end
scores = scores(:);
if ~all(isfinite(scores)) || all(scores <= 0)
    scores = ones(size(x));
end
end

function v = bounded_shape_param_ms(shape_params, fname, default_v, lo, hi)
v = default_v;
if isstruct(shape_params) && isfield(shape_params, fname)
    raw = shape_params.(fname);
    if isnumeric(raw) && ~isempty(raw) && isfinite(raw(1))
        v = raw(1);
    end
end
v = min(max(v, lo), hi);
end

function [ok, pts] = select_exact_count_symmetric_weighted_ms(points, n_target, point_scores, symmetry_mode)
ok = false;
pts = zeros(0,2);
if isempty(points) || n_target < 1 || size(points,1) < n_target
    return;
end
if nargin < 4 || isempty(symmetry_mode)
    symmetry_mode = 'c4';
end
[orbits, ~, ok_orbit] = group_points_by_mirror_orbit_ms(points, symmetry_mode);
if ~ok_orbit || isempty(orbits)
    return;
end
sizes = cellfun(@numel, orbits).';
vals = zeros(numel(orbits), 1);
for i = 1:numel(orbits)
    vals(i) = sum(point_scores(orbits{i}));
end
[ok_dp, pick] = subset_exact_dp_ms(sizes, vals, n_target);
if ~ok_dp
    return;
end
[ok, pts] = points_from_orbit_pick_ms(points, orbits, pick, n_target);
end

function pts_out = apply_shape_jitter_ms(pts, shape_params, fp_w, fp_h, min_edge_gap, halfx, halfy)
pts_out = pts;
jr = bounded_shape_param_ms(shape_params, 'jitter_ratio', 0, 0, 0.45);
if jr <= 0 || isempty(pts)
    return;
end
seed_val = bounded_shape_param_ms(shape_params, 'jitter_seed', 1, -1e9, 1e9);
stage_idx = 1;
if isstruct(shape_params) && isfield(shape_params, 'stage_idx') && isfinite(shape_params.stage_idx)
    stage_idx = shape_params.stage_idx;
end
stream = RandStream('mt19937ar', 'Seed', max(0, round(seed_val + 10007 * stage_idx)));
[orbits, ~, ok_orbit] = group_points_by_mirror_orbit_ms(pts, 'c4');
if ~ok_orbit || isempty(orbits)
    return;
end
for attempt = 1:5
    pts_try = pts;
    for oi = 1:numel(orbits)
        scale = 1 + (rand(stream, 1, 1) - 0.5) * 2 * jr;
        pts_try(orbits{oi},:) = pts(orbits{oi},:) * scale;
    end
    pts_try(:,1) = min(max(pts_try(:,1), -halfx), halfx);
    pts_try(:,2) = min(max(pts_try(:,2), -halfy), halfy);
    rects = centers_to_rects_ms(pts_try, fp_w, fp_h);
    [spacing_ok, ~] = validate_stage_rect_spacing_fast_ms(rects, fp_w, fp_h, min_edge_gap);
    if spacing_ok
        pts_out = pts_try;
        return;
    end
end
end

function [ok, pts] = select_exact_count_symmetric_ms(points, n_target, mode, symmetry_mode, edge_pattern_mode)
ok = false;
pts = zeros(0,2);
if isempty(points) || n_target < 1 || size(points,1) < n_target
    return;
end
[orbits, radii, ok_orbit] = group_points_by_mirror_orbit_ms(points, symmetry_mode);
if ~ok_orbit || isempty(orbits)
    return;
end
sizes = cellfun(@numel, orbits).';
if strcmpi(mode, 'edge_dense')
    vals = radii .* sizes;
else
    vals = (1 ./ (1 + radii)) .* sizes;
end
if nargin < 5 || isempty(edge_pattern_mode)
    edge_pattern_mode = 'free';
end
edge_pattern_mode = lower(strtrim(char(edge_pattern_mode)));
if ~any(strcmp(edge_pattern_mode, {'free', 'edge_spaced', 'edge_clean'}))
    edge_pattern_mode = 'free';
end

[ok_dp0, pick0] = subset_exact_dp_ms(sizes, vals, n_target);
if ~ok_dp0
    return;
end
[ok_pts0, pts0, sel0] = points_from_orbit_pick_ms(points, orbits, pick0, n_target);
if ~ok_pts0
    return;
end
if strcmp(edge_pattern_mode, 'free')
    ok = true;
    pts = pts0;
    return;
end

[side_meta, boundary_mask] = build_edge_pattern_side_meta_ms(points);
if check_edge_pattern_ms(side_meta, sel0, edge_pattern_mode)
    ok = true;
    pts = pts0;
    return;
end

orbit_boundary_touch = zeros(numel(orbits), 1);
for i = 1:numel(orbits)
    orbit_boundary_touch(i) = sum(boundary_mask(orbits{i}));
end
lambda_list = [0.05, 0.1, 0.2, 0.4];
for lam = lambda_list
    vals_adj = vals - lam * orbit_boundary_touch;
    [ok_dp, pick] = subset_exact_dp_ms(sizes, vals_adj, n_target);
    if ~ok_dp
        continue;
    end
    [ok_pts, pts_try, sel_try] = points_from_orbit_pick_ms(points, orbits, pick, n_target);
    if ~ok_pts
        continue;
    end
    if check_edge_pattern_ms(side_meta, sel_try, edge_pattern_mode)
        ok = true;
        pts = pts_try;
        return;
    end
end
end

function [ok, pts, sel] = points_from_orbit_pick_ms(points, orbits, pick, n_target)
ok = false;
pts = zeros(0,2);
sel = zeros(0,1);
for i = 1:numel(pick)
    if pick(i)
        sel = [sel; orbits{i}(:)]; %#ok<AGROW>
    end
end
if numel(sel) ~= n_target
    return;
end
sel = unique(sel, 'stable');
if numel(sel) ~= n_target
    return;
end
pts = sortrows(points(sel,:), [1 2]);
ok = true;
end

function [side_meta, boundary_mask] = build_edge_pattern_side_meta_ms(points)
tol = 1e-9;
qx = int64(round(points(:,1) / tol));
qy = int64(round(points(:,2) / tol));
xmin = min(qx); xmax = max(qx);
ymin = min(qy); ymax = max(qy);
left_ids = find(qx == xmin);
right_ids = find(qx == xmax);
bottom_ids = find(qy == ymin);
top_ids = find(qy == ymax);
[~, il] = sort(qy(left_ids), 'ascend');
[~, ir] = sort(qy(right_ids), 'ascend');
[~, ib] = sort(qx(bottom_ids), 'ascend');
[~, it] = sort(qx(top_ids), 'ascend');
side_meta = struct();
side_meta.left = left_ids(il);
side_meta.right = right_ids(ir);
side_meta.bottom = bottom_ids(ib);
side_meta.top = top_ids(it);
boundary_mask = false(size(points,1),1);
boundary_mask(side_meta.left) = true;
boundary_mask(side_meta.right) = true;
boundary_mask(side_meta.bottom) = true;
boundary_mask(side_meta.top) = true;
end

function ok = check_edge_pattern_ms(side_meta, sel_ids, edge_pattern_mode)
ok = true;
edge_pattern_mode = lower(strtrim(char(edge_pattern_mode)));
side_names = {'left', 'right', 'bottom', 'top'};
for i = 1:numel(side_names)
    ids = side_meta.(side_names{i});
    Ns = numel(ids);
    if Ns <= 0
        continue;
    end
    occ = ismember(ids, sel_ids);
    n_occ = sum(occ);
    if n_occ == 0
        continue;
    end
    if strcmp(edge_pattern_mode, 'edge_clean')
        min_side_fill_ratio = min(0.35, max(0.10, 2 / Ns));
        min_segment_len = max(3, ceil(0.12 * Ns));
    else
        min_side_fill_ratio = min(0.20, max(0.04, 1 / Ns));
        min_segment_len = 2;
    end
    fill_ratio = n_occ / Ns;
    if fill_ratio < min_side_fill_ratio
        ok = false;
        return;
    end
    if ~all_true_runs_ge_ms(occ, min_segment_len)
        ok = false;
        return;
    end
end
end

function ok = all_true_runs_ge_ms(mask, min_len)
ok = true;
if isempty(mask)
    return;
end
cur = 0;
for i = 1:numel(mask)
    if mask(i)
        cur = cur + 1;
    else
        if cur > 0 && cur < min_len
            ok = false;
            return;
        end
        cur = 0;
    end
end
if cur > 0 && cur < min_len
    ok = false;
end
end

function [orbits, radii, ok] = group_points_by_mirror_orbit_ms(points, symmetry_mode)
orbits = {};
radii = zeros(0,1);
ok = true;
if isempty(points)
    return;
end
tol = 1e-9;
qx = int64(round(points(:,1) / tol));
qy = int64(round(points(:,2) / tol));
map = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:numel(qx)
    map(sprintf('%d_%d', qx(i), qy(i))) = i;
end
used = false(size(points,1),1);
for i = 1:size(points,1)
    if used(i)
        continue;
    end
    x = qx(i); y = qy(i);
    switch lower(symmetry_mode)
        case 'c2_lr'
            q = [x y; -x y];
        case 'c2_ud'
            q = [x y; x -y];
        otherwise
            q = [x y; -x y; x -y; -x -y];
    end
    idx = zeros(0,1);
    for k = 1:size(q,1)
        key = sprintf('%d_%d', q(k,1), q(k,2));
        if ~isKey(map, key)
            ok = false;
            return;
        end
        idx(end+1,1) = map(key); %#ok<AGROW>
    end
    idx = unique(idx, 'stable');
    if isempty(idx)
        ok = false;
        return;
    end
    used(idx) = true;
    orbits{end+1,1} = idx; %#ok<AGROW>
    p = points(idx(1),:);
    radii(end+1,1) = hypot(p(1), p(2)); %#ok<AGROW>
end
end

function [ok, pick] = subset_exact_dp_ms(sizes, vals, target)
m = numel(sizes);
pick = false(1, m);
ok = false;
dp = -inf(m+1, target+1);
choose = false(m+1, target+1);
dp(1,1) = 0;
for i = 1:m
    s = sizes(i);
    v = vals(i);
    for c = 0:target
        cur = dp(i, c+1);
        if ~isfinite(cur)
            continue;
        end
        if cur > dp(i+1, c+1)
            dp(i+1, c+1) = cur;
        end
        if c + s <= target
            cand = cur + v;
            if cand > dp(i+1, c+s+1)
                dp(i+1, c+s+1) = cand;
                choose(i+1, c+s+1) = true;
            end
        end
    end
end
if ~isfinite(dp(end, target+1))
    return;
end
c = target;
for i = m:-1:1
    if choose(i+1, c+1)
        pick(i) = true;
        c = c - sizes(i);
    end
end
ok = (c == 0);
end

function [ok, pts] = make_regular_fullgrid_points_ms(n_target, Lmax, fp_w, fp_h, edge_margin_fp, s_dense, anis_ratio, min_edge_gap)
    ok = false; pts = zeros(0,2);
    if n_target < 1, return; end

    % 鎸夌敤鎴蜂紶鍏ョ殑 s_dense 涓?anis_ratio 鐢熸垚瑙勫垯缃戞牸锛屼繚鐣欏悓鍙橀噺浣撶郴鍩虹嚎
    min_dx = max(fp_w + min_edge_gap, fp_w * 1.01);
    min_dy = max(fp_h + min_edge_gap, fp_h * 1.01);
    sx = max(s_dense * anis_ratio, min_dx);
    sy = max(s_dense / anis_ratio, min_dy);

    ny = round(sqrt(n_target * (sx / sy)));
    ny = max(1, ny);
    nx = ceil(n_target / ny);

    vx = ((0:nx-1) - (nx-1)/2) * sx;
    vy = ((0:ny-1) - (ny-1)/2) * sy;
    [X, Y] = meshgrid(vx, vy);
    pool = [X(:), Y(:)];

    [ok_sel, pts] = select_exact_count_symmetric_ms(pool, n_target, 'center_dense', 'c4', 'free');
    if ok_sel, ok = true; end
end

function [ok, pts_best, nx_best, ny_best] = make_uniform_points_adaptive_fixed_box_ms(n_target, Lx_fixed, Ly_fixed, fp_w, fp_h, edge_margin_fp, min_step_scale, gamma_val)
if nargin < 8
    gamma_val = NaN;
end
ok = false;
pts_best = zeros(0,2);
nx_best = 0;
ny_best = 0;
if n_target < 1 || Lx_fixed <= 0 || Ly_fixed <= 0
    return;
end
fp = max(fp_w, fp_h);
step_min = fp * max(1.10, min_step_scale);
halfx = Lx_fixed/2 - fp_w/2 - edge_margin_fp * fp;
halfy = Ly_fixed/2 - fp_h/2 - edge_margin_fp * fp;
if halfx <= 0 || halfy <= 0
    return;
end
nx_cap = max(1, floor((2*halfx) / step_min) + 1);
ny_cap = max(1, floor((2*halfy) / step_min) + 1);

best_score = Inf;
for nx = 1:min(nx_cap, n_target)
    ny = ceil(n_target / nx);
    if ny > ny_cap || ny < 1
        continue;
    end
    if nx > 1
        step_x = (2*halfx)/(nx-1);
        if step_x < step_min
            continue;
        end
    else
        step_x = 0;
    end
    if ny > 1
        step_y = (2*halfy)/(ny-1);
        if step_y < step_min
            continue;
        end
    else
        step_y = 0;
    end
    xvals = ((0:nx-1) - (nx-1)/2) * step_x;
    yvals = ((0:ny-1) - (ny-1)/2) * step_y;
    [ok_sel, pts0] = select_single_side_rows_ms(xvals, yvals, n_target);
    if ~ok_sel
        continue;
    end
    [ok_shift, pts_shift, score_shift] = optimize_uniform_shift_ms(pts0, Lx_fixed, Ly_fixed, fp_w, fp_h);
    if ~ok_shift
        continue;
    end
    ratio_pen = abs(log(max(nx/max(ny,1), 1e-12) / max(Lx_fixed/max(Ly_fixed,1e-12), 1e-12)));
    gamma_pen = 0;
    if isfinite(gamma_val)
        gamma_pen = 0.01 * abs(gamma_val);
    end
    score = score_shift + 0.02 * ratio_pen + gamma_pen;
    if score < best_score
        best_score = score;
        pts_best = pts_shift;
        nx_best = nx;
        ny_best = ny;
        ok = true;
    end
end
end

function [ok, pts] = select_single_side_rows_ms(xvals, yvals, n_target)
ok = false;
pts = zeros(0,2);
xvals = xvals(:).';
yvals = yvals(:).';
nx = numel(xvals);
ny = numel(yvals);
if nx < 1 || ny < 1 || n_target > nx * ny
    return;
end
ny_full = floor(n_target / nx);
rem = n_target - ny_full * nx;
sel = zeros(n_target, 2);
ptr = 1;
if ny_full > 0
    [Xf, Yf] = meshgrid(xvals, yvals(1:ny_full));
    block = [Xf(:), Yf(:)];
    sel(ptr:ptr+size(block,1)-1,:) = block;
    ptr = ptr + size(block,1);
end
if rem > 0
    y_idx = min(ny_full + 1, ny);
    order = center_out_order_ms(nx);
    idx_part = sort(order(1:rem));
    x_part = xvals(idx_part(:));
    y_part = repmat(yvals(y_idx), numel(x_part), 1);
    block = [x_part(:), y_part];
    sel(ptr:ptr+rem-1,:) = block;
    ptr = ptr + rem;
end
if ptr - 1 ~= n_target
    return;
end
pts = sortrows(sel, [1 2]);
if size(unique(pts, 'rows'), 1) ~= n_target
    return;
end
ok = true;
end

function order = center_out_order_ms(n)
if n <= 0
    order = zeros(0,1);
    return;
end
if mod(n,2) == 1
    c = (n + 1) / 2;
    order = c;
    k = 1;
    while numel(order) < n
        if c - k >= 1
            order(end+1) = c - k; %#ok<AGROW>
        end
        if c + k <= n
            order(end+1) = c + k; %#ok<AGROW>
        end
        k = k + 1;
    end
else
    c1 = n/2;
    c2 = c1 + 1;
    order = [c1 c2];
    k = 1;
    while numel(order) < n
        if c1 - k >= 1
            order(end+1) = c1 - k; %#ok<AGROW>
        end
        if c2 + k <= n
            order(end+1) = c2 + k; %#ok<AGROW>
        end
        k = k + 1;
    end
end
order = order(:);
end

function [ok, pts_best, score_best] = optimize_uniform_shift_ms(pts0, Lx, Ly, fp_w, fp_h)
ok = false;
pts_best = zeros(0,2);
score_best = Inf;
if isempty(pts0)
    return;
end
xL = min(pts0(:,1) - fp_w/2);
xR = max(pts0(:,1) + fp_w/2);
yB = min(pts0(:,2) - fp_h/2);
yT = max(pts0(:,2) + fp_h/2);
dx_min = -Lx/2 - xL;
dx_max = Lx/2 - xR;
dy_min = -Ly/2 - yB;
dy_max = Ly/2 - yT;
if dx_min > dx_max || dy_min > dy_max
    return;
end
dx_grid = linspace(dx_min, dx_max, 9);
dy_grid = linspace(dy_min, dy_max, 9);
for ix = 1:numel(dx_grid)
    for iy = 1:numel(dy_grid)
        pts = pts0 + [dx_grid(ix), dy_grid(iy)];
        mleft = min((pts(:,1) - fp_w/2) + Lx/2);
        mright = min(Lx/2 - (pts(:,1) + fp_w/2));
        mbot = min((pts(:,2) - fp_h/2) + Ly/2);
        mtop = min(Ly/2 - (pts(:,2) + fp_h/2));
        if any([mleft, mright, mbot, mtop] < -1e-12)
            continue;
        end
        sym_err = abs(mleft - mright) + abs(mbot - mtop);
        score = sym_err;
        if score < score_best
            score_best = score;
            pts_best = pts;
            ok = true;
        end
    end
end
end

function [ok, pts] = force_feasible_layout_ms(n_target, mode, symmetry_mode, edge_pattern_mode, Lmax, fp_w, fp_h, edge_margin_fp, min_edge_gap)
ok = false;
pts = zeros(0,2);
fp = max(fp_w, fp_h);
base_step = max([fp_w + min_edge_gap, fp_h + min_edge_gap, fp * 1.10]);
step_list = base_step * [1.30, 1.15, 1.00];
for s = step_list
    pool = make_uniform_pool_points_ms(Lmax, s, fp_w, fp_h, edge_margin_fp, min_edge_gap);
    if size(pool,1) < n_target
        continue;
    end
    [ok_sel, pts_sel] = select_exact_count_symmetric_ms(pool, n_target, mode, symmetry_mode, edge_pattern_mode);
    if ok_sel
        ok = true;
        pts = pts_sel;
        return;
    end
end
end

function hc = resolve_hard_constraints_ms(spec)
hc = struct('force_c4_only', true, 'enforce_monotonic_geometry', true, ...
    'enforce_coverage_min', true, 'enforce_spacing_min_edge_gap', true);
if nargin < 1 || ~isstruct(spec) || ~isfield(spec, 'hard_constraints') || ~isstruct(spec.hard_constraints)
    return;
end
src = spec.hard_constraints;
fn = fieldnames(hc);
for i = 1:numel(fn)
    k = fn{i};
    if isfield(src, k)
        hc.(k) = any(logical(src.(k)));
    end
end
end

function [ok, msg] = validate_candidate_hard_constraints_ms(symmetry_mode, stages, Lx, Ly, cov, spec, hc)
ok = true;
msg = '';
N = spec.stage_count;
tol = 1e-12;
if hc.force_c4_only && ~strcmpi(symmetry_mode, 'c4')
    ok = false;
    msg = 'symmetry_not_c4';
    return;
end
for k = 1:N
    relax_ratio = 1.0;
    if numel(stages) >= k && isstruct(stages{k}) && isfield(stages{k}, 'lmax_relax_ratio') && isfinite(stages{k}.lmax_relax_ratio)
        relax_ratio = max(1.0, stages{k}.lmax_relax_ratio);
    end
    if hc.force_c4_only
        if ~is_stage_geometry_c4_ms(stages{k}.rects)
            ok = false;
            msg = sprintf('c4_geom_stage%d', k);
            return;
        end
    end
    if max(Lx(k), Ly(k)) > spec.geometry.L_max(k) * relax_ratio + tol
        ok = false;
        msg = sprintf('size_stage%d', k);
        return;
    end
    if hc.enforce_coverage_min && cov(k) < spec.geometry.coverage_min(k) - tol
        ok = false;
        msg = sprintf('coverage_stage%d', k);
        return;
    end
    if hc.enforce_spacing_min_edge_gap
        [ok_spacing, ~] = validate_stage_rect_spacing_fast_ms(stages{k}.rects, spec.fp_w, spec.fp_h, spec.geometry.min_edge_gap);
        if ~ok_spacing
            ok = false;
            msg = sprintf('spacing_stage%d', k);
            return;
        end
    end
end
if hc.enforce_monotonic_geometry
    for k = 1:(N-1)
        if Lx(k) + tol < (Lx(k+1) + spec.geometry.pyramid_gap_min(k))
            ok = false;
            msg = sprintf('mono_lx_stage%d', k);
            return;
        end
        if Ly(k) + tol < (Ly(k+1) + spec.geometry.pyramid_gap_min(k))
            ok = false;
            msg = sprintf('mono_ly_stage%d', k);
            return;
        end
    end
end
end

function tf = is_stage_geometry_c4_ms(rects)
tf = true;
if isempty(rects)
    return;
end
if size(rects,2) ~= 4 || ~all(isfinite(rects(:)))
    tf = false;
    return;
end
pts = rect_centers_ms(rects);
[~, ~, ok_orbit] = group_points_by_mirror_orbit_ms(pts, 'c4');
tf = logical(ok_orbit);
end

function pts = make_uniform_pool_points_ms(span, step, fp_w, fp_h, edge_margin_fp, min_edge_gap)
fp = max(fp_w, fp_h);
half = span/2 - fp/2 - edge_margin_fp * fp;
if half <= 0
    pts = zeros(0,2);
    return;
end
step = max(step, max([fp_w + min_edge_gap, fp_h + min_edge_gap, fp * 1.10]));
p = 0:step:half;
if isempty(p)
    p = 0;
end
p = p(:).';
if numel(p) > 1
    full = [-p(end:-1:2), p];
else
    full = p;
end
[X,Y] = meshgrid(full, full);
pts = [X(:), Y(:)];
end

function centers = rect_centers_ms(rects)
centers = [0.5*(rects(:,1)+rects(:,2)), 0.5*(rects(:,3)+rects(:,4))];
end

function rects = centers_to_rects_ms(points, w, h)
rects = [points(:,1)-w/2, points(:,1)+w/2, points(:,2)-h/2, points(:,2)+h/2];
end

function [Lx, Ly] = infer_plate_size_xy_ms(points, w, h, edge_margin_fp)
if isempty(points)
    Lx = NaN;
    Ly = NaN;
    return;
end
fp = max(w, h);
margin = edge_margin_fp * fp;
Lx = 2 * (max(abs(points(:,1))) + w/2 + margin);
Ly = 2 * (max(abs(points(:,2))) + h/2 + margin);
end

function [ok, diag] = validate_stage_rect_spacing_fast_ms(rects, fp_w, fp_h, min_edge_gap)
diag = struct('overlap_pair_count', 0, 'violating_pair_count', 0);
ok = true;
if isempty(rects) || size(rects, 1) < 2
    return;
end
centers = rect_centers_ms(rects);
x = centers(:,1);
y = centers(:,2);
[x, ord] = sort(x);
y = y(ord);
dx_req = fp_w + min_edge_gap;
dy_req = fp_h + min_edge_gap;
tol = 1e-12;
n = numel(x);
for i = 1:n-1
    j = i + 1;
    while j <= n && (x(j) - x(i)) < dx_req - tol
        dy = abs(y(j) - y(i));
        if dy < dy_req - tol
            diag.violating_pair_count = diag.violating_pair_count + 1;
            if (x(j) - x(i)) < fp_w - tol && dy < fp_h - tol
                diag.overlap_pair_count = diag.overlap_pair_count + 1;
            end
        end
        j = j + 1;
    end
end
ok = (diag.violating_pair_count == 0) && (diag.overlap_pair_count == 0);
end

function diag = calc_stage_spacing_metrics_ms(rects, fp_w, fp_h, min_edge_gap)
diag = struct('overlap_pair_count', 0, 'violating_pair_count', 0, 'min_edge_gap', inf);
if isempty(rects) || size(rects, 1) < 2
    return;
end
centers = rect_centers_ms(rects);
x = centers(:,1);
y = centers(:,2);
n = numel(x);
dx_req = fp_w + min_edge_gap;
dy_req = fp_h + min_edge_gap;
tol = 1e-12;
for i = 1:n-1
    for j = i+1:n
        dx = abs(x(i) - x(j));
        dy = abs(y(i) - y(j));
        sx = max(0, dx - fp_w);
        sy = max(0, dy - fp_h);
        d_edge = hypot(sx, sy);
        if d_edge < diag.min_edge_gap
            diag.min_edge_gap = d_edge;
        end
        if dx < dx_req - tol && dy < dy_req - tol
            diag.violating_pair_count = diag.violating_pair_count + 1;
        end
        if dx < fp_w - tol && dy < fp_h - tol
            diag.overlap_pair_count = diag.overlap_pair_count + 1;
        end
    end
end
end

function [ok, pts] = make_hex_c6_points_ms(n_target, Lmax, fp_w, fp_h, edge_margin_fp, mode, s_base, anis_ratio, min_edge_gap)
ok = false;
pts = zeros(0,2);
if n_target < 1
    return;
end
fp = max(fp_w, fp_h);
margin = fp/2 + edge_margin_fp * fp;
half = Lmax/2 - margin;
if half <= fp * 0.5
    return;
end
sx = max(s_base * anis_ratio, max(fp_w + min_edge_gap, fp_w * 1.01));
sy = max(s_base / anis_ratio, max(fp_h + min_edge_gap, fp_h * 1.01));
a2y = sy * sqrt(3)/2;
i_max = ceil(half / min(sx, a2y)) + 2;
pts_all = zeros(0,2);
for i = -i_max:i_max
    for j = -i_max:i_max
        px = i*sx + j*sx/2;
        py = j*a2y;
        if abs(px) + fp_w/2 <= Lmax/2 - edge_margin_fp*fp && ...
                abs(py) + fp_h/2 <= Lmax/2 - edge_margin_fp*fp
            pts_all(end+1,:) = [px, py]; %#ok<AGROW>
        end
    end
end
if size(pts_all,1) < n_target
    return;
end
[orbits, radii, ok_g] = group_points_by_c6_orbit_ms(pts_all);
if ~ok_g || isempty(orbits)
    return;
end
sizes = cellfun(@numel, orbits).';
if strcmpi(mode, 'center_dense')
    vals = 1./(1 + radii) .* sizes;
else
    vals = radii .* sizes;
end
[ok_dp, pick] = subset_exact_dp_ms(sizes, vals, n_target);
if ~ok_dp
    return;
end
sel = [];
for ii = 1:numel(pick)
    if pick(ii)
        sel = [sel; orbits{ii}(:)]; %#ok<AGROW>
    end
end
sel = unique(sel);
if numel(sel) ~= n_target
    return;
end
pts = pts_all(sel,:);
ok = true;
end

function [orbits, radii, ok] = group_points_by_c6_orbit_ms(points)
ok = true;
orbits = {};
radii = zeros(0,1);
if isempty(points)
    return;
end
tol = 1e-7;
N = size(points,1);
used = false(N,1);
qi = round(points(:,1)/tol);
qj_arr = round(points(:,2)/tol);
map = containers.Map('KeyType','char','ValueType','double');
for i = 1:N
    map(sprintf('%d_%d', qi(i), qj_arr(i))) = i;
end
c60 = 0.5; s60 = sqrt(3)/2;
R60 = [c60, -s60; s60, c60];
for i = 1:N
    if used(i)
        continue;
    end
    orbit_idx = [];
    p = points(i,:).';
    for k = 0:5
        key = sprintf('%d_%d', round(p(1)/tol), round(p(2)/tol));
        if isKey(map, key)
            orbit_idx(end+1) = map(key); %#ok<AGROW>
        end
        p = R60 * p;
    end
    orbit_idx = unique(orbit_idx(:));
    used(orbit_idx) = true;
    orbits{end+1,1} = orbit_idx; %#ok<AGROW>
    radii(end+1,1) = hypot(points(i,1), points(i,2)); %#ok<AGROW>
end
end

function [ok, pts] = make_radial_gamma_points_ms(n_target, Lmax, fp_w, fp_h, edge_margin_fp, gamma, s_base, anis_ratio, min_edge_gap)
ok = false;
pts = zeros(0,2);
if n_target < 1
    return;
end
expo = 1 + 9 * abs(gamma);
if gamma >= 0
    mode = 'center_dense';
else
    mode = 'edge_dense';
end
sxd = max(s_base * anis_ratio, max(fp_w + min_edge_gap, fp_w * 1.01));
syd = max(s_base / anis_ratio, max(fp_h + min_edge_gap, fp_h * 1.01));
sxs = sxd * 2.0;
sys = syd * 2.0;
fp = max(fp_w, fp_h);
Llo = max(2.5e-3, fp*1.2);
Lhi = Lmax;
if count_points_xy_ms(Lhi, mode, sxd, sxs, syd, sys, expo, fp_w, fp_h, edge_margin_fp) < n_target
    return;
end
for it = 1:24
    Lmid = 0.5*(Llo+Lhi);
    if count_points_xy_ms(Lmid, mode, sxd, sxs, syd, sys, expo, fp_w, fp_h, edge_margin_fp) >= n_target
        Lhi = Lmid;
    else
        Llo = Lmid;
    end
end
for Lt = linspace(Lhi, Lmax, 6)
    rects = make_gradient_footprints_xy_ms(sxd, sxs, syd, sys, expo, Lt, fp_w, fp_h, mode, edge_margin_fp);
    p = rect_centers_ms(rects);
    [ok_sel, pts] = select_exact_count_symmetric_ms(p, n_target, mode, 'c4', 'free');
    if ok_sel
        ok = true;
        return;
    end
end
end

function [ok, pts] = make_concentric_ring_points_ms(n_target, Lmax, fp_w, fp_h, edge_margin_fp, mode, min_edge_gap)
ok = false;
pts = zeros(0,2);
if n_target < 1
    return;
end
fp = max(fp_w, fp_h);
margin = fp/2 + edge_margin_fp * fp;
r_max_sq = Lmax/2 - margin;
if r_max_sq <= fp * 0.5
    return;
end
r_step = max(fp * 1.35, max(fp_w + min_edge_gap, fp_h + min_edge_gap));
K_max = floor(r_max_sq / r_step);
if K_max < 1
    return;
end
rings_r = zeros(K_max+1, 1);
rings_cap = zeros(K_max+1, 1);
rings_r(1) = 0;
rings_cap(1) = 1;
for k = 1:K_max
    rings_r(k+1) = k * r_step;
    rings_cap(k+1) = 6 * k;
end
if strcmpi(mode, 'center_dense')
    order = 1:numel(rings_r);
else
    order = numel(rings_r):-1:1;
end
pts_all = zeros(0,2);
for ri = order
    rk = rings_r(ri);
    cap_k = rings_cap(ri);
    if rk == 0
        pts_all(end+1,:) = [0, 0]; %#ok<AGROW>
    else
        angles = (0:cap_k-1) * 2*pi / cap_k;
        for ai = 1:cap_k
            px = rk * cos(angles(ai));
            py = rk * sin(angles(ai));
            if abs(px) + fp_w/2 <= Lmax/2 - edge_margin_fp*fp && ...
                    abs(py) + fp_h/2 <= Lmax/2 - edge_margin_fp*fp
                pts_all(end+1,:) = [px, py]; %#ok<AGROW>
            end
        end
    end
    if size(pts_all,1) >= n_target
        break;
    end
end
if size(pts_all,1) < n_target
    return;
end
pts = pts_all(1:n_target,:);
ok = true;
end

function [ok, pts] = make_poisson_disk_points_ms(n_target, Lmax, fp_w, fp_h, edge_margin_fp, min_edge_gap, seed_val)
ok = false;
pts = zeros(0,2);
if n_target < 1
    return;
end
if nargin < 7 || ~isfinite(seed_val)
    seed_val = 42;
end
stream = RandStream('mt19937ar', 'Seed', round(seed_val));

fp = max(fp_w, fp_h);
r_min = max(fp + min_edge_gap, fp * 1.05);
margin = fp/2 + edge_margin_fp * fp;
half_x = Lmax/2 - margin;
half_y = Lmax/2 - margin;
if half_x <= r_min || half_y <= r_min
    return;
end

cell_size = r_min / sqrt(2);
nx_cell = max(1, ceil((2*half_x) / cell_size));
ny_cell = max(1, ceil((2*half_y) / cell_size));
grid = zeros(ny_cell, nx_cell);

k_tries = 30;
n_alloc = max(16, n_target * 4);
pts_placed = zeros(n_alloc, 2);
active = zeros(n_alloc, 1);
n_active = 0;
n_placed = 0;

p0x = (rand(stream) * 2 - 1) * half_x * 0.5;
p0y = (rand(stream) * 2 - 1) * half_y * 0.5;
n_placed = n_placed + 1;
pts_placed(n_placed,:) = [p0x, p0y];
[gx0, gy0] = to_grid_index_poisson_ms(p0x, p0y, half_x, half_y, cell_size, nx_cell, ny_cell);
grid(gy0, gx0) = n_placed;
n_active = n_active + 1;
active(n_active) = n_placed;

while n_active > 0 && n_placed < max(n_target * 3, n_target + 8)
    ai = max(1, floor(rand(stream) * n_active) + 1);
    src_idx = active(ai);
    src = pts_placed(src_idx,:);
    found = false;
    for t = 1:k_tries
        ang = rand(stream) * 2 * pi;
        rad = r_min * (1 + rand(stream));
        px = src(1) + rad * cos(ang);
        py = src(2) + rad * sin(ang);
        if abs(px) + fp_w/2 > half_x || abs(py) + fp_h/2 > half_y
            continue;
        end
        if ~is_far_enough_poisson_ms(px, py, pts_placed, n_placed, grid, ...
                half_x, half_y, cell_size, nx_cell, ny_cell, r_min)
            continue;
        end
        n_placed = n_placed + 1;
        if n_placed > size(pts_placed,1)
            pts_placed = [pts_placed; zeros(n_target, 2)]; %#ok<AGROW>
        end
        pts_placed(n_placed,:) = [px, py];
        [gx, gy] = to_grid_index_poisson_ms(px, py, half_x, half_y, cell_size, nx_cell, ny_cell);
        grid(gy, gx) = n_placed;
        n_active = n_active + 1;
        if n_active > numel(active)
            active = [active; zeros(n_target, 1)]; %#ok<AGROW>
        end
        active(n_active) = n_placed;
        found = true;
    end
    if ~found
        active(ai) = active(n_active);
        n_active = n_active - 1;
    end
end

if n_placed < n_target
    return;
end
pts_placed = pts_placed(1:n_placed,:);
r2 = pts_placed(:,1).^2 + pts_placed(:,2).^2;
[~, ord] = sort(r2, 'ascend');
pts = pts_placed(ord(1:n_target), :);
ok = true;
end

function [gx, gy] = to_grid_index_poisson_ms(px, py, half_x, half_y, cell_size, nx_cell, ny_cell)
gx = floor((px + half_x) / cell_size) + 1;
gy = floor((py + half_y) / cell_size) + 1;
gx = min(nx_cell, max(1, gx));
gy = min(ny_cell, max(1, gy));
end

function tf = is_far_enough_poisson_ms(px, py, pts_placed, n_placed, grid, ...
    half_x, half_y, cell_size, nx_cell, ny_cell, r_min)
tf = true;
[gx0, gy0] = to_grid_index_poisson_ms(px, py, half_x, half_y, cell_size, nx_cell, ny_cell);
for dgx = -2:2
    for dgy = -2:2
        gi = gx0 + dgx;
        gj = gy0 + dgy;
        if gi < 1 || gi > nx_cell || gj < 1 || gj > ny_cell
            continue;
        end
        idx = grid(gj, gi);
        if idx < 1 || idx > n_placed
            continue;
        end
        dx = px - pts_placed(idx,1);
        dy = py - pts_placed(idx,2);
        if hypot(dx, dy) < r_min - 1e-12
            tf = false;
            return;
        end
    end
end
end

function [ok, pts] = make_interstage_aligned_points_ms(n_target, Lmax, fp_w, fp_h, ...
    edge_margin_fp, min_edge_gap, prev_stage_rects, prev2_stage_rects, alignment_weight)
ok = false;
pts = zeros(0,2);
if n_target < 1
    return;
end
if nargin < 9 || ~isfinite(alignment_weight)
    alignment_weight = 0.6;
end
alignment_weight = min(max(alignment_weight, 0), 1);

fp = max(fp_w, fp_h);
margin = fp/2 + edge_margin_fp * fp;
half_x = Lmax/2 - margin;
half_y = Lmax/2 - margin;
if half_x <= 0 || half_y <= 0
    return;
end

step = max(fp + min_edge_gap, fp * 1.10);
xs = -half_x:step:half_x;
ys = -half_y:step:half_y;
[X, Y] = meshgrid(xs, ys);
pool = [X(:), Y(:)];
if size(pool,1) < n_target
    return;
end

r = hypot(pool(:,1), pool(:,2));
center_score = 1 - r / max(max(r), 1e-12);

prev_cold = zeros(0,2);
if ~isempty(prev_stage_rects) && size(prev_stage_rects,2) == 4
    prev_cold = rect_centers_ms(prev_stage_rects);
end
prev_hot = zeros(0,2);
if ~isempty(prev2_stage_rects) && size(prev2_stage_rects,2) == 4
    prev_hot = rect_centers_ms(prev2_stage_rects);
elseif ~isempty(prev_cold)
    prev_hot = prev_cold;
end

if isempty(prev_cold)
    score = center_score;
else
    cold_dist = nearest_distance_to_points_ms(pool, prev_cold);
    cold_score = 1 - cold_dist / max(max(cold_dist), 1e-12);
    if isempty(prev_hot)
        hot_score = cold_score;
    else
        hot_dist = nearest_distance_to_points_ms(pool, prev_hot);
        hot_score = 1 - hot_dist / max(max(hot_dist), 1e-12);
    end
    mix_align = 0.6 * cold_score + 0.4 * hot_score;
    score = alignment_weight * mix_align + (1 - alignment_weight) * center_score;
end

[~, ord] = sort(score, 'descend');
dx_req = fp_w + min_edge_gap;
dy_req = fp_h + min_edge_gap;
pts_sel = zeros(n_target, 2);
placed = 0;
for i = 1:numel(ord)
    if placed >= n_target
        break;
    end
    px = pool(ord(i),1);
    py = pool(ord(i),2);
    if placed > 0
        dx = abs(pts_sel(1:placed,1) - px);
        dy = abs(pts_sel(1:placed,2) - py);
        if any(dx < dx_req - 1e-12 & dy < dy_req - 1e-12)
            continue;
        end
    end
    placed = placed + 1;
    pts_sel(placed,:) = [px, py];
end
if placed < n_target
    return;
end
pts = pts_sel;
ok = true;
end

function dmin = nearest_distance_to_points_ms(pool_pts, ref_pts)
if isempty(pool_pts)
    dmin = zeros(0,1);
    return;
end
if isempty(ref_pts)
    dmin = inf(size(pool_pts,1), 1);
    return;
end
dmin = inf(size(pool_pts,1), 1);
for i = 1:size(pool_pts,1)
    dx = ref_pts(:,1) - pool_pts(i,1);
    dy = ref_pts(:,2) - pool_pts(i,2);
    d2 = dx.^2 + dy.^2;
    dmin(i) = sqrt(min(d2));
end
end

%% ====================== Step3: evaluation ======================

function calib = empty_current_calib_ms(I0)
if nargin < 1 || ~isfinite(I0)
    I0 = NaN;
end
calib = struct('enabled', false, 'has_solution', false, ...
    'I0', I0, 'I_cal', I0, 'I_shift', 0, ...
    'I_grid', zeros(0,1), 'top_ids', zeros(0,1), 'peak_I', zeros(0,1), ...
    'peak_score', zeros(0,1), 'seed_source', 'direct_full_fem', 'seed_ids', zeros(0,1));
end

function evals = evaluate_candidates_ms(cands, G, spec, keep_fields)
if nargin < 4
    keep_fields = false;
end
N = numel(cands);
evals = repmat(empty_eval_struct_ms(spec.stage_count), N, 1);
if N == 0
    return;
end
    use_par = spec.use_parallel && N > 2;
    if use_par
        pool = gcp('nocreate');
        if isempty(pool)
            ensure_pool_ms(spec);
            pool = gcp('nocreate');
        end
        if isempty(pool)
            use_par = false;
        else
            nw = max(pool.NumWorkers, 1);
            active_limit = min([64, N, nw]);
            fprintf('[Eval] evaluating %d candidates with %d workers, dynamic_inflight=%d\n', N, nw, active_limit);
            futures = parallel.FevalFuture.empty(0, 1);
            future_candidate_idx = zeros(0, 1);
            future_start_time = zeros(0, 1);
            next_idx = 1;
            completed = 0;
            while next_idx <= N && numel(futures) < active_limit
                cand_id = get_candidate_id_for_log_ms(cands(next_idx));
                fprintf('[Eval] candidate %03d/%03d submitted: id=%g\n', next_idx, N, cand_id);
                futures(end+1, 1) = parfeval(pool, @evaluate_single_candidate_ms, 1, cands(next_idx), G, spec, keep_fields); %#ok<AGROW>
                future_candidate_idx(end+1, 1) = next_idx; %#ok<AGROW>
                future_start_time(end+1, 1) = now; %#ok<AGROW>
                next_idx = next_idx + 1;
            end
            while ~isempty(futures)
                [done_slot, rec_done] = fetchNext(futures);
                cand_idx = future_candidate_idx(done_slot);
                elapsed_s = (now - future_start_time(done_slot)) * 86400;
                evals(cand_idx) = rec_done;
                completed = completed + 1;
                cand_id = get_candidate_id_for_log_ms(cands(cand_idx));
                if rec_done.success
                    fprintf('[Eval] candidate %03d/%03d done: id=%g, time=%.2fs, completed=%d/%d\n', ...
                        cand_idx, N, cand_id, elapsed_s, completed, N);
                else
                    fprintf('[Eval] candidate %03d/%03d failed: id=%g, time=%.2fs, completed=%d/%d, message=%s\n', ...
                        cand_idx, N, cand_id, elapsed_s, completed, N, char(string(rec_done.message)));
                end
                futures(done_slot) = [];
                future_candidate_idx(done_slot) = [];
                future_start_time(done_slot) = [];
                if next_idx <= N
                    cand_id = get_candidate_id_for_log_ms(cands(next_idx));
                    fprintf('[Eval] candidate %03d/%03d submitted: id=%g\n', next_idx, N, cand_id);
                    futures(end+1, 1) = parfeval(pool, @evaluate_single_candidate_ms, 1, cands(next_idx), G, spec, keep_fields); %#ok<AGROW>
                    future_candidate_idx(end+1, 1) = next_idx; %#ok<AGROW>
                    future_start_time(end+1, 1) = now; %#ok<AGROW>
                    next_idx = next_idx + 1;
                end
            end
        end
    end
    if ~use_par
        for i = 1:N
            evals(i) = evaluate_single_candidate_ms(cands(i), G, spec, keep_fields);
        end
    end
end

function cand_id = get_candidate_id_for_log_ms(cand)
cand_id = NaN;
if isstruct(cand) && isfield(cand, 'candidate_id') && ~isempty(cand.candidate_id)
    cand_id = cand.candidate_id;
end
end

function rec = evaluate_single_candidate_ms(cand, G, spec, keep_fields)
if nargin < 4
    keep_fields = false;
end
N = spec.stage_count;
rec = empty_eval_struct_ms(N);
rec.candidate_id = cand.candidate_id;
rec.layout_method = cand.layout_method;
if isfield(cand, 'search_group')
    rec.search_group = cand.search_group;
end
rec.stage_modes = cand.stage_modes;
if isfield(cand, 'stage_methods')
    rec.stage_methods = cand.stage_methods;
end
rec.stage_trends = cand.stage_trends;
rec.symmetry_mode = cand.symmetry_mode;
rec.edge_pattern_mode = cand.edge_pattern_mode;
rec.s_dense = cand.s_dense;
rec.s_sparse = cand.s_sparse;
rec.expo = cand.expo;
rec.anis_ratio = cand.anis_ratio;
rec.gamma = cand.gamma;
rec.method_anchor_stage = cand.method_anchor_stage;
if isfield(cand, 'shape_mode'), rec.shape_mode = cand.shape_mode; end
if isfield(cand, 'stage_anis'), rec.stage_anis = cand.stage_anis; end
if isfield(cand, 'ring_radius_ratio'), rec.ring_radius_ratio = cand.ring_radius_ratio; end
if isfield(cand, 'ring_width_ratio'), rec.ring_width_ratio = cand.ring_width_ratio; end
if isfield(cand, 'band_width_ratio'), rec.band_width_ratio = cand.band_width_ratio; end
if isfield(cand, 'corner_bias'), rec.corner_bias = cand.corner_bias; end
if isfield(cand, 'jitter_ratio'), rec.jitter_ratio = cand.jitter_ratio; end
if isfield(cand, 'jitter_seed'), rec.jitter_seed = cand.jitter_seed; end
rec.spacing_ratio = cand.spacing_ratio;
rec.contrast_score = cand.contrast_score;
rec.mode_prior_score = calc_mode_prior_score_ms(cand.stage_trends, get_soft_prior_ms(spec), spec.stage_count);
rec.lmax_relaxed = cand.lmax_relaxed;
rec.lmax_relax_ratio = cand.lmax_relax_ratio;
if isfield(cand, 'count_rank')
    rec.count_rank = cand.count_rank;
end
rec.n = cand.n;
[ok_pair_count, npair] = particle_counts_to_pair_counts_ms(cand.n);
if ok_pair_count
    rec.npair = npair;
else
    rec.message = 'invalid particle counts: every stage n must be even for couple pairing';
    return;
end
rec.ratios = cand.ratios;
rec.Lx = cand.Lx;
rec.Ly = cand.Ly;
rec.cov = cand.cov;
rec.fp_fill_count = 0;
rec.fp_fill_by_stage = zeros(1, N);
rec.stage_rects = cand.stage_rects;
G_eval = G;
if isfield(cand, 'I_opt') && isfinite(cand.I_opt) && cand.I_opt > 0
    G_eval.I = cand.I_opt;
end
rec.I_opt = G_eval.I;

try
    for k = 1:N
        [ok_spacing, spacing_diag] = validate_stage_rect_spacing_fast_ms(cand.stage_rects{k}, ...
            spec.fp_w, spec.fp_h, spec.geometry.min_edge_gap);
        if ~ok_spacing
            rec.message = sprintf(['stage %d spacing violation: overlap_pairs=%d, ' ...
                'violating_pairs=%d, req_gap=%.3f mm'], ...
                k, spacing_diag.overlap_pair_count, spacing_diag.violating_pair_count, ...
                spec.geometry.min_edge_gap_mm);
            return;
        end
    end

    x0 = linspace(spec.targets.ThRes - 20, spec.targets.ThRes - 120, N).';
    [ok0d, C] = solve_0d_at_current_ms(cand.n, G_eval.I, G_eval, spec, x0);
    if ~ok0d
        rec.message = '0D solve failed';
        return;
    end
    C_init = C(:).';

    G_eval.z_path_runtime = make_zpath_runtime_ms(spec, G.fp_A, N);
    plates = cell(1, N);
    fill_count_total = 0;
    fill_count_by_stage = zeros(1, N);
    for k = 1:N
        [nx_use, ny_use] = resolve_stage_mesh_for_eval_ms(spec, k);
        plates{k} = build_plate_struct_ms(struct( ...
            'Lx', cand.Lx(k), 'Ly', cand.Ly(k), 't', spec.plate_t, ...
            'nx', nx_use, 'ny', ny_use, 'k_inplane', spec.plate_k_inplane));
    end

    for k = 1:N
        if k < N
            [fpIdCool, fpIdHot] = map_elements_to_footprints_ms( ...
                plates{k}.elem.cx, plates{k}.elem.cy, cand.stage_rects{k}, cand.stage_rects{k+1});
        else
            [fpIdCool, ~] = map_elements_to_footprints_ms( ...
                plates{k}.elem.cx, plates{k}.elem.cy, cand.stage_rects{k}, zeros(0,4));
            fpIdHot = zeros(size(fpIdCool));
        end
        plates{k}.fpElemsCool = ids_to_cell_ms(fpIdCool, cand.n(k));
        [plates{k}.fpElemsCool, fill_cool] = fill_empty_fp_elems_nearest_ms( ...
            plates{k}.fpElemsCool, cand.stage_rects{k}, plates{k}.elem.cx, plates{k}.elem.cy);
        fill_count_total = fill_count_total + fill_cool;
        fill_count_by_stage(k) = fill_count_by_stage(k) + fill_cool;
        if any(cellfun(@isempty, plates{k}.fpElemsCool))
            rec.message = sprintf('empty cool-side footprint mapping at stage %d', k);
            return;
        end
        if k < N
            plates{k}.fpElemsHot = ids_to_cell_ms(fpIdHot, cand.n(k+1));
            [plates{k}.fpElemsHot, fill_hot] = fill_empty_fp_elems_nearest_ms( ...
                plates{k}.fpElemsHot, cand.stage_rects{k+1}, plates{k}.elem.cx, plates{k}.elem.cy);
            fill_count_total = fill_count_total + fill_hot;
            fill_count_by_stage(k) = fill_count_by_stage(k) + fill_hot;
            if any(cellfun(@isempty, plates{k}.fpElemsHot))
                rec.message = sprintf('empty hot-side footprint mapping at stage %d', k);
                return;
            end
        else
            plates{k}.fpElemsHot = cell(0,1);
        end
    end

    best = newton_solve_Ns_ms(G_eval, plates, cand.n, C_init, N);
    used_continuation = false;
    if ~best.success
        [best_cont, cont_msg] = newton_solve_with_continuation_ms(G_eval, plates, cand.n, C_init, N);
        if best_cont.success
            best = best_cont;
            used_continuation = true;
        else
            rec.message = sprintf('newton not converged (direct+continuation failed: %s)', cont_msg);
            return;
        end
    end
    if ~best.success
        rec.message = 'newton not converged';
        return;
    end

    fields = cell(1, N);
    rec.stage_spread = zeros(1, N);
    for k = 1:N
        fields{k} = best.theta{k} + best.C(k);
        rec.stage_spread(k) = max(fields{k}) - min(fields{k});
    end
    rec.C = best.C(:).';
    Tmin_last = min(fields{N});
    Tmax_last = max(fields{N});
    rec.DeltaTN_actual = spec.targets.ThRes - Tmin_last;
    rec.TN_min = Tmin_last;
    rec.TN_mean = sum(plates{N}.m .* fields{N}) * plates{N}.inv_sum_m;
    rec.DeltaTN_mean = spec.targets.ThRes - rec.TN_mean;
    rec.TN_maxmin = Tmax_last - Tmin_last;
    rec.fp_fill_count = fill_count_total;
    rec.fp_fill_by_stage = fill_count_by_stage;
    rec.newton_relaxed = logical(best.relaxed);
    rec.newton_rel_max = best.rel_max;
    rec.newton_iters = best.iters;
    rec.newton_convergence_mode = best.convergence_mode;
    if used_continuation
        rec.newton_convergence_mode = ['continuation_', char(string(best.convergence_mode))];
    end
    rec.success = true;
    rec.message = 'ok';
    rec.rank_score = calc_rank_score_ms(rec, spec);

    if keep_fields
        rec.plates = plates;
        rec.Tfields = fields;
    end
catch ME
    rec.success = false;
    rec.message = ME.message;
end
end

function [nx_use, ny_use] = resolve_stage_mesh_for_eval_ms(spec, stage_idx)
nx_use = max(5, round(spec.mesh_nx));
ny_use = max(5, round(spec.mesh_ny));
if isfield(spec, 'mesh_nx_stage_full') && numel(spec.mesh_nx_stage_full) >= stage_idx
    nx_use = max(5, round(spec.mesh_nx_stage_full(stage_idx)));
end
if isfield(spec, 'mesh_ny_stage_full') && numel(spec.mesh_ny_stage_full) >= stage_idx
    ny_use = max(5, round(spec.mesh_ny_stage_full(stage_idx)));
end
end


function [best_out, status_msg] = newton_solve_with_continuation_ms(G, plates, n_vec, C_init, N)
best_out = struct('C', C_init(:).', 'theta', {cell(1, N)}, 'g', NaN(N,1), ...
    'cache', struct([]), 'success', false, 'iters', 0);
status_msg = 'schedule_not_started';
I_target = G.I;
if ~(isfinite(I_target) && I_target > 0)
    status_msg = 'invalid_target_current';
    return;
end
I_sched = unique([max(0.6, 0.55 * I_target), max(0.8, 0.75 * I_target), max(1.0, 0.90 * I_target), I_target]);
I_sched = I_sched(isfinite(I_sched) & I_sched > 0);
if isempty(I_sched)
    status_msg = 'empty_schedule';
    return;
end

C_seed = C_init(:).';
for i = 1:numel(I_sched)
    Gs = G;
    Gs.I = I_sched(i);
    cand_i = newton_solve_Ns_ms(Gs, plates, n_vec, C_seed, N);
    if ~cand_i.success
        status_msg = sprintf('schedule_fail_at_I=%.3fA', I_sched(i));
        return;
    end
    C_seed = cand_i.C(:).';
    best_out = cand_i;
end
status_msg = 'schedule_success';
end

function [plates, fields] = build_surrogate_temperature_fields_ms(cand, C, spec)
N = spec.stage_count;
plates = cell(1, N);
fields = cell(1, N);
for k = 1:N
    p = build_plate_struct_ms(struct('Lx', cand.Lx(k), 'Ly', cand.Ly(k), ...
        'nx', spec.mesh_nx, 'ny', spec.mesh_ny));
    x = p.mesh.x;
    y = p.mesh.y;
    r = sqrt((x / max(cand.Lx(k), eps)).^2 + (y / max(cand.Ly(k), eps)).^2);
    trend = cand.stage_trends{k};
    switch trend
        case 'center_heavy'
            phi = -0.8 * exp(-8 * r.^2);
        case 'edge_heavy'
            phi = -0.8 * (1 - exp(-8 * r.^2));
        otherwise
            phi = -0.4 * exp(-6 * r.^2);
    end
    scale = 0.8 + 0.05 * k;
    fields{k} = C(k) + scale * phi;
    plates{k} = p;
end
end

function best = newton_solve_Ns_ms(G, plates, n_vec, C_init_vec, N)
seeds = build_newton_initial_seeds_ms(C_init_vec, N, G);
best = struct('C', seeds(1,:), 'theta', {cell(1, N)}, 'g', NaN(N,1), ...
    'cache', struct([]), 'success', false, 'iters', 0);
best_gnorm = inf;
for i = 1:size(seeds, 1)
    cand = newton_solve_single_start_ms(G, plates, n_vec, seeds(i,:), N);
    gnorm = inf;
    if isfield(cand, 'g') && ~isempty(cand.g)
        gvec = cand.g(:);
        gvec = gvec(isfinite(gvec));
        if ~isempty(gvec)
            gnorm = norm(gvec);
        end
    end
    if cand.success
        if ~best.success || gnorm < best_gnorm
            best = cand;
            best_gnorm = gnorm;
        end
    elseif ~best.success && gnorm < best_gnorm
        best = cand;
        best_gnorm = gnorm;
    end
end
end

function seeds = build_newton_initial_seeds_ms(C_init_vec, N, G)
base = C_init_vec(:).';
if numel(base) ~= N || any(~isfinite(base))
    base = linspace(G.ThRes - 20, G.ThRes - 120, N);
end
fallback = linspace(G.ThRes - 15, G.ThRes - 130, N);
grad = linspace(0, -12, N);
seeds = [base; fallback; base + grad; base - grad];
mask = ~isfinite(seeds);
if any(mask(:))
    fill = repmat(base, size(seeds,1), 1);
    seeds(mask) = fill(mask);
end
[~, ia] = unique(round(seeds * 1e6), 'rows', 'stable');
seeds = seeds(sort(ia), :);
end

function best = newton_solve_single_start_ms(G, plates, n_vec, C0, N)
C = C0(:).';
theta_guess = cell(1, N);
for k = 1:N
    theta_guess{k} = zeros(plates{k}.mesh.Nn, 1);
end
best = struct('C', C, 'theta', {theta_guess}, 'g', NaN(N,1), ...
    'cache', struct([]), 'success', false, 'iters', 0, ...
    'relaxed', false, 'convergence_mode', 'none', 'rel_max', inf);

eps_C = 0.1;
lambda = 1e-4;
newton_rc_min = 1e-12;
relaxed_rel_factor = 3.0;
relaxed_step_limit = 0.35;
step_norm_cap = max(0.25, 0.18 * sqrt(max(1, N)));
for it = 1:G.max_outer
    [g, theta, cache] = eval_gvec_Ns_ms(C, theta_guess, G, plates, n_vec, N);
    if any(~isfinite(g))
        return;
    end
    best = struct('C', C, 'theta', {theta}, 'g', g, 'cache', cache, 'success', false, 'iters', it, ...
        'relaxed', false, 'convergence_mode', 'none', 'rel_max', inf);

    J = zeros(N, N);
    for j = 1:N
        Cj = C;
        Cj(j) = Cj(j) + eps_C;
        gp = eval_g_only_Ns_ms(Cj, theta, G, plates, n_vec, N);
        if any(~isfinite(gp))
            return;
        end
        J(:,j) = (gp - g) / eps_C;
    end
    if any(~isfinite(J(:)))
        return;
    end

    delta_newton = NaN(N,1);
    [J_use, g_use] = scale_newton_system_ms(J, g);
    rc = rcond(J_use);
    if isfinite(rc) && rc > newton_rc_min
        delta_newton = -J_use \ g_use;
    end
    [delta_lm, lambda, ok_lm] = solve_regularized_lm_step_ms(J_use, g_use, lambda);
    if ~ok_lm
        return;
    end
    if all(isfinite(delta_newton))
        delta = delta_newton;
    else
        delta = delta_lm;
    end
    if any(~isfinite(delta))
        return;
    end
    dnorm = norm(delta);
    if isfinite(dnorm) && dnorm > step_norm_cap
        delta = delta * (step_norm_cap / dnorm);
    end

    accepted = false;
    step = 1.0;
    g_norm0 = norm(g);
    for ls = 1:10
        C_try = C + step * delta(:).';
        if any(~isfinite(C_try))
            step = step * 0.5;
            continue;
        end
        [g_try, theta_try, cache_try] = eval_gvec_Ns_ms(C_try, theta, G, plates, n_vec, N);
        if all(isfinite(g_try))
            g_try_norm = norm(g_try);
            if g_try_norm <= g_norm0 * (1 - 1e-4 * step) || g_try_norm < g_norm0 - 1e-12
                accepted = true;
                break;
            end
        end
        step = step * 0.5;
    end
    if ~accepted
        lambda = min(1e10, lambda * 10);
        if lambda >= 1e10
            return;
        end
        continue;
    end

    C = C_try;
    theta_guess = theta_try;
    g = g_try;
    cache = cache_try;
    rel_max = calc_rel_max_newton_ms(g, cache, N);
    best = struct('C', C, 'theta', {theta_guess}, 'g', g, ...
        'cache', cache, 'success', false, 'iters', it, ...
        'relaxed', false, 'convergence_mode', 'none', 'rel_max', rel_max);
    if rel_max < G.tol_g_rel && norm(step * delta) < 0.1
        [g_cv, theta_cv, cache_cv] = eval_gvec_Ns_ms(C, theta_guess, G, plates, n_vec, N);
        if all(isfinite(g_cv))
            best.C = C;
            best.theta = theta_cv;
            best.g = g_cv;
            best.cache = cache_cv;
            best.rel_max = calc_rel_max_newton_ms(g_cv, cache_cv, N);
        end
        best.success = true;
        best.relaxed = false;
        best.convergence_mode = 'strict';
        return;
    end
    if rel_max < relaxed_rel_factor * G.tol_g_rel && norm(step * delta) < relaxed_step_limit
        [g_cv, theta_cv, cache_cv] = eval_gvec_Ns_ms(C, theta_guess, G, plates, n_vec, N);
        if all(isfinite(g_cv))
            best.C = C;
            best.theta = theta_cv;
            best.g = g_cv;
            best.cache = cache_cv;
            best.rel_max = calc_rel_max_newton_ms(g_cv, cache_cv, N);
        end
        best.success = true;
        best.relaxed = true;
        best.convergence_mode = 'relaxed';
        return;
    end
    lambda = max(1e-12, lambda / 3);
end
end

function [delta_lm, lambda_out, ok] = solve_regularized_lm_step_ms(J, g, lambda_in)
N = size(J, 2);
delta_lm = NaN(N,1);
lambda_out = lambda_in;
ok = false;
if ~(isfinite(lambda_in) && lambda_in > 0)
    lambda_out = 1e-4;
end
w1 = warning('off', 'MATLAB:nearlySingularMatrix');
w2 = warning('off', 'MATLAB:singularMatrix');
cleanup = onCleanup(@() restore_warnings_ms(w1, w2));
for t = 1:8
    A_lm = J.' * J + lambda_out * eye(N);
    b_lm = J.' * g;
    if any(~isfinite(A_lm(:))) || any(~isfinite(b_lm))
        lambda_out = min(1e12, lambda_out * 10);
        continue;
    end
    d_try = -(A_lm \ b_lm);
    if all(isfinite(d_try))
        delta_lm = d_try;
        ok = true;
        return;
    end
    lambda_out = min(1e12, lambda_out * 10);
end
end

function [J_use, g_use] = scale_newton_system_ms(J, g)
J_use = J;
g_use = g;
if isempty(J) || isempty(g)
    return;
end
rscale = max(1, max(abs(J), [], 2));
rscale = max(rscale, abs(g));
rscale = max(rscale, 1e-10);
W = spdiags(1 ./ rscale, 0, numel(rscale), numel(rscale));
J_use = W * J;
g_use = W * g;
end

function restore_warnings_ms(w1, w2)
if isstruct(w1) && isfield(w1, 'state')
    warning(w1.state, 'MATLAB:nearlySingularMatrix');
end
if isstruct(w2) && isfield(w2, 'state')
    warning(w2.state, 'MATLAB:singularMatrix');
end
end

function rel_max = calc_rel_max_newton_ms(g, cache, N)
rel_max = inf;
if ~isfield(cache, 'sumQc') || ~isfield(cache, 'sumQh') || ...
        numel(cache.sumQc) ~= N || numel(cache.sumQh) ~= N
    return;
end
refs = zeros(N,1);
refs(1) = max(abs(cache.sumQc(N)) + 1e-12, 1e-12);
for ridx = 2:N
    j = N - ridx + 2;
    refs(ridx) = max(abs(cache.sumQh(j)) + abs(cache.sumQc(j-1)), 1e-12);
end
rel = abs(g) ./ refs;
if all(isfinite(rel))
    rel_max = max(rel);
end
end

function [g, theta, cache] = eval_gvec_Ns_ms(C_vec, theta_cells_0, G, plates, n_vec, N)
[theta, cache] = solve_theta_given_C_ms(C_vec, theta_cells_0, G, plates, n_vec, N);
g = NaN(N,1);
if ~isfield(cache, 'sumQc') || ~isfield(cache, 'sumQh') || ...
        numel(cache.sumQc) ~= N || numel(cache.sumQh) ~= N
    return;
end
g(1) = cache.sumQc(N) - G.Qc_target_last;
for ridx = 2:N
    j = N - ridx + 2;
    g(ridx) = cache.sumQh(j) - cache.sumQc(j-1);
end
end

function g = eval_g_only_Ns_ms(C_vec, theta_cells_0, G, plates, n_vec, N)
[g, ~, ~] = eval_gvec_Ns_ms(C_vec, theta_cells_0, G, plates, n_vec, N);
end

function [theta, cache] = solve_theta_given_C_ms(C_vec, theta_cells_0, G, plates, n_vec, N)
theta = theta_cells_0;
for k = 1:N
    if k > numel(theta) || isempty(theta{k}) || numel(theta{k}) ~= plates{k}.mesh.Nn
        theta{k} = zeros(plates{k}.mesh.Nn, 1);
    end
end
zcfg = get_zpath_runtime_from_G_ms(G, N);
w_prev = zcfg.step_q_prev_weight;
if ~zcfg.enable
    w_prev = 0;
end
[ok_pair_count, npair] = particle_counts_to_pair_counts_ms(n_vec);
if ~ok_pair_count
    cache = make_nonfinite_cache_ms(zcfg, w_prev, 'invalid_even_particle_counts', N);
    return;
end
pair_idx = cell(1, N);
for k = 1:N
    pair_idx{k} = make_particle_pair_indices_ms(n_vec(k));
end

Qc_prev = cell(1, N);
Qh_prev = cell(1, N);
for k = 1:N
    Qc_prev{k} = zeros(npair(k), 1);
    Qh_prev{k} = zeros(npair(k), 1);
end

for it = 1:G.max_inner
    Tfields = cell(1, N);
    for k = 1:N
        Tfields{k} = theta{k} + C_vec(k);
    end

    Tc_plate = cell(1, N);
    Th_plate = cell(1, N);
    for k = 1:N
        Tc_plate{k} = footprint_avg_ms(Tfields{k}, plates{k}.fpElemsCool, plates{k}.elem, plates{k}.mesh.tri);
        if any(~isfinite(Tc_plate{k}))
            cache = make_nonfinite_cache_ms(zcfg, w_prev, 'nonfinite_plate_temperature', N);
            return;
        end
        if k >= 2
            Th_plate{k} = footprint_avg_ms(Tfields{k-1}, plates{k-1}.fpElemsHot, plates{k-1}.elem, plates{k-1}.mesh.tri);
            if any(~isfinite(Th_plate{k}))
                cache = make_nonfinite_cache_ms(zcfg, w_prev, 'nonfinite_plate_temperature', N);
                return;
            end
        else
            Th_plate{k} = G.ThRes + zeros(size(Tc_plate{k}));
        end
    end

    Tc_pair_plate = cell(1, N);
    Th_pair_plate = cell(1, N);
    Tc_eff = cell(1, N);
    Th_eff = cell(1, N);
    for k = 1:N
        Tc_pair_plate{k} = pair_average_values_ms(Tc_plate{k}, pair_idx{k});
        Th_pair_plate{k} = pair_average_values_ms(Th_plate{k}, pair_idx{k});
        if zcfg.enable
            if k == 1
                Tc_eff{k} = Tc_pair_plate{k} - Qc_prev{k} * zcfg.Rz_interfaces(1);
                Th_eff{k} = G.ThRes + Qh_prev{k} * zcfg.Rz_sink;
            elseif k == N
                Tc_eff{k} = Tc_pair_plate{k} - Qc_prev{k} * zcfg.Rz_interfaces(N-1);
                Th_eff{k} = Th_pair_plate{k} + Qh_prev{k} * zcfg.Rz_interfaces(N-1);
            else
                Tc_eff{k} = Tc_pair_plate{k} - Qc_prev{k} * zcfg.Rz_interfaces(k);
                Th_eff{k} = Th_pair_plate{k} + Qh_prev{k} * zcfg.Rz_interfaces(k-1);
            end
        else
            Tc_eff{k} = Tc_pair_plate{k};
            if k == 1
                Th_eff{k} = G.ThRes + zeros(size(Tc_pair_plate{k}));
            else
                Th_eff{k} = Th_pair_plate{k};
            end
        end
        if any(~isfinite(Tc_eff{k})) || any(~isfinite(Th_eff{k}))
            cache = make_nonfinite_cache_ms(zcfg, w_prev, 'nonfinite_leg_or_sink_temperature', N);
            return;
        end
    end

    Qc = cell(1, N);
    Qh = cell(1, N);
    for k = 1:N
        Qc{k} = zeros(npair(k), 1);
        Qh{k} = zeros(npair(k), 1);
        for j = 1:npair(k)
            Qc{k}(j) = te_Qc_onecouple_ms(Tc_eff{k}(j), Th_eff{k}(j), G.I, G, k);
            Qh{k}(j) = te_Qh_onecouple_ms(Th_eff{k}(j), Tc_eff{k}(j), G.I, G, k);
        end
        if any(~isfinite(Qc{k})) || any(~isfinite(Qh{k}))
            cache = make_nonfinite_cache_ms(zcfg, w_prev, 'nonfinite_q_from_te_model', N);
            return;
        end
    end

    Qc_eff = cell(1, N);
    Qh_eff = cell(1, N);
    for k = 1:N
        Qc_eff{k} = (1 - w_prev) * Qc{k} + w_prev * Qc_prev{k};
        Qh_eff{k} = (1 - w_prev) * Qh{k} + w_prev * Qh_prev{k};
        Qc_prev{k} = Qc_eff{k};
        Qh_prev{k} = Qh_eff{k};
        if any(~isfinite(Qc_eff{k})) || any(~isfinite(Qh_eff{k}))
            cache = make_nonfinite_cache_ms(zcfg, w_prev, 'nonfinite_q_after_relax', N);
            return;
        end
    end

    F = cell(1, N);
    for k = 1:N
        F{k} = zeros(plates{k}.mesh.Nn, 1);

        for p = 1:numel(Qc_eff{k})
            qj = -0.5 * Qc_eff{k}(p) / G.fp_A;
            F{k} = add_heat_to_footprint_elems_ms(F{k}, plates{k}.fpElemsCool, ...
                plates{k}.elem, plates{k}.mesh.tri, pair_idx{k}(p,1), qj);
            F{k} = add_heat_to_footprint_elems_ms(F{k}, plates{k}.fpElemsCool, ...
                plates{k}.elem, plates{k}.mesh.tri, pair_idx{k}(p,2), qj);
        end

        if k < N
            for p = 1:numel(Qh_eff{k+1})
                qj = +0.5 * Qh_eff{k+1}(p) / G.fp_A;
                F{k} = add_heat_to_footprint_elems_ms(F{k}, plates{k}.fpElemsHot, ...
                    plates{k}.elem, plates{k}.mesh.tri, pair_idx{k+1}(p,1), qj);
                F{k} = add_heat_to_footprint_elems_ms(F{k}, plates{k}.fpElemsHot, ...
                    plates{k}.elem, plates{k}.mesh.tri, pair_idx{k+1}(p,2), qj);
            end
        end

        if any(~isfinite(F{k}))
            cache = make_nonfinite_cache_ms(zcfg, w_prev, 'nonfinite_rhs', N);
            return;
        end
    end

    dmax = 0;
    for k = 1:N
        sol = plates{k}.Aaug \ [F{k}; 0];
        if any(~isfinite(sol))
            cache = make_nonfinite_cache_ms(zcfg, w_prev, 'nonfinite_linear_solution', N);
            return;
        end
        theta_new = sol(1:end-1);
        theta_new = theta_new - (plates{k}.m' * theta_new) * plates{k}.inv_sum_m;
        theta_next = theta{k} + G.omega * (theta_new - theta{k});
        dk = max(abs(theta_next - theta{k}));
        if ~isfinite(dk)
            cache = make_nonfinite_cache_ms(zcfg, w_prev, 'nonfinite_theta_update', N);
            return;
        end
        dmax = max(dmax, dk);
        theta{k} = theta_next;
    end

    cache = struct();
    cache.sumQc = zeros(N,1);
    cache.sumQh = zeros(N,1);
    cache.step_q_prev_weight = w_prev;
    cache.Rz_interfaces = zcfg.Rz_interfaces(:).';
    cache.Rz_sink = zcfg.Rz_sink;
    for k = 1:N
        cache.sumQc(k) = sum(Qc_eff{k});
        cache.sumQh(k) = sum(Qh_eff{k});
    end

    if dmax < G.tol_theta
        return;
    end
end
end

function cache = make_nonfinite_cache_ms(zcfg, w_prev, reason, N)
if nargin < 4 || isempty(N)
    N = 5;
end
cache = struct();
cache.abort_reason = reason;
cache.sumQc = NaN(N,1);
cache.sumQh = NaN(N,1);
cache.step_q_prev_weight = w_prev;
cache.Rz_interfaces = zcfg.Rz_interfaces(:).';
cache.Rz_sink = zcfg.Rz_sink;
end

function rec = empty_eval_struct_ms(N)
if nargin < 1
    N = 5;
end
rec = struct('success', false, 'message', '', ...
    'candidate_id', -1, 'layout_method', '', ...
    'search_group', 'regular', ...
    'stage_modes', {repmat({''}, 1, N)}, ...
    'stage_methods', {repmat({''}, 1, N)}, ...
    'stage_trends', {repmat({'neutral'}, 1, N)}, ...
    'symmetry_mode', '', 'edge_pattern_mode', '', ...
    's_dense', NaN, 's_sparse', NaN, 'expo', NaN, 'anis_ratio', NaN, 'gamma', NaN, 'method_anchor_stage', NaN, ...
    'shape_mode', '', 'stage_anis', NaN(1,N), 'ring_radius_ratio', NaN, ...
    'ring_width_ratio', NaN, 'band_width_ratio', NaN, 'corner_bias', NaN, ...
    'jitter_ratio', NaN, 'jitter_seed', NaN, ...
    'spacing_ratio', NaN, 'contrast_score', NaN, 'mode_prior_score', NaN, ...
    'lmax_relaxed', false, 'lmax_relax_ratio', 1.0, ...
    'count_rank', NaN, ...
    'n', NaN(1,N), 'npair', NaN(1,N), 'ratios', NaN(1,max(0,N-1)), 'I_opt', NaN, ...
    'Lx', NaN(1,N), 'Ly', NaN(1,N), 'cov', NaN(1,N), ...
    'fp_fill_count', 0, 'fp_fill_by_stage', zeros(1,N), ...
    'newton_relaxed', false, 'newton_rel_max', NaN, 'newton_iters', NaN, 'newton_convergence_mode', 'none', ...
    'C', NaN(1,N), ...
    'DeltaTN_actual', NaN, 'TN_min', NaN, 'TN_mean', NaN, 'DeltaTN_mean', NaN, 'TN_maxmin', NaN, ...
    'stage_spread', NaN(1,N), 'rank_score', NaN, ...
    'plates', {cell(1,N)}, 'Tfields', {cell(1,N)}, ...
    'stage_rects', {cell(1,N)});
end

function s = calc_contrast_score_ms(s_dense, s_sparse, expo)
if ~(isfinite(s_dense) && s_dense > 0 && isfinite(s_sparse) && isfinite(expo))
    s = -inf;
    return;
end
s = (s_sparse / s_dense)^2 * expo;
end

function prior = get_soft_prior_ms(spec)
prior = disable_soft_prior_ms(spec);
prior = prior(1);
if ~isstruct(spec) || ~isfield(spec, 'soft_prior') || ~isstruct(spec.soft_prior)
    return;
end
sp = spec.soft_prior;
sp = sp(1);
if isfield(sp, 'enable')
    v_enable = sp.enable;
    if isempty(v_enable)
        % keep default
    elseif islogical(v_enable)
        prior.enable = any(v_enable(:));
    elseif isnumeric(v_enable)
        prior.enable = any(v_enable(:) ~= 0);
    elseif ischar(v_enable)
        s_enable = lower(strtrim(v_enable));
        prior.enable = any(strcmp(s_enable, {'1', 'true', 'yes', 'on'}));
    elseif isa(v_enable, 'string') && isscalar(v_enable)
        s_enable = lower(strtrim(char(v_enable)));
        prior.enable = any(strcmp(s_enable, {'1', 'true', 'yes', 'on'}));
    end
end
if isfield(sp, 'stage_trend_prefer')
    prior.stage_trend_prefer = normalize_trend_cell_ms(sp.stage_trend_prefer, spec.stage_count);
end
end

function s = calc_mode_prior_score_ms(stage_trends, prior, N)
if nargin < 3
    N = numel(stage_trends);
end
if nargin < 2 || isempty(prior) || ~isfield(prior, 'enable') || ~prior.enable
    s = NaN;
    return;
end
s = 0;
tr = normalize_trend_cell_ms(stage_trends, N);
pr = normalize_trend_cell_ms(prior.stage_trend_prefer, N);
for k = 1:N
    if strcmpi(tr{k}, pr{k})
        s = s + 1;
    end
end
end

function sorted = sort_candidates_by_deltaTN_ms(arr, spec)
if nargin < 2
    spec = struct();
end
if isempty(arr)
    sorted = arr;
    return;
end
score = [arr.DeltaTN_mean].';
score(~isfinite(score)) = -inf;
actual = [arr.DeltaTN_actual].';
actual(~isfinite(actual)) = -inf;
rank = [arr.rank_score].';
rank(~isfinite(rank)) = -inf;
spread = [arr.TN_maxmin].';
spread_score = -spread;
spread_score(~isfinite(spread)) = -inf;
[~, idx] = sortrows([actual, score, rank, spread_score, [arr.candidate_id].'], [-1, -2, -3, -4, 5]);
sorted = arr(idx);
end

function cand = find_candidate_by_id_ms(cands, cid)
idx = find([cands.candidate_id] == cid, 1, 'first');
if isempty(idx)
    error('Candidate id %d not found.', cid);
end
cand = cands(idx);
end

function s = calc_rank_score_ms(ev, spec)
if ~ev.success || ~isfinite(ev.DeltaTN_actual) || ~isfinite(ev.DeltaTN_mean)
    s = -inf;
    return;
end
w1 = 0.6; w2 = 0.3; w3 = 0.1;
if isstruct(spec) && isfield(spec, 'rank')
    if isfield(spec.rank, 'w_actual'), w1 = spec.rank.w_actual; end
    if isfield(spec.rank, 'w_mean'), w2 = spec.rank.w_mean; end
    if isfield(spec.rank, 'w_spread'), w3 = spec.rank.w_spread; end
end
dt_ref = max(abs(spec.targets.DeltaT_target), 1.0);
spread_ref = max(0.20 * dt_ref, 1.0);

s_actual = ev.DeltaTN_actual / dt_ref;
s_mean = ev.DeltaTN_mean / dt_ref;
spread = ev.TN_maxmin;
if ~isfinite(spread) || spread < 0
    spread = 0;
end
s_spread = spread / spread_ref;
s = w1 * s_actual + w2 * s_mean - w3 * s_spread;
end

function arr = renormalize_rank_scores_ms(arr, spec)
if isempty(arr) || ~isfield(arr, 'success')
    return;
end
ok_mask = logical([arr.success]);
n_ok = sum(ok_mask);
if n_ok < 2
    return;
end

ok_idx = find(ok_mask);
dt_a = arrayfun(@(e) e.DeltaTN_actual, arr(ok_idx));
dt_m = arrayfun(@(e) e.DeltaTN_mean, arr(ok_idx));
sp = arrayfun(@(e) ifelse_finite_ms(e.TN_maxmin, 0), arr(ok_idx));

da_min = min(dt_a); da_max = max(dt_a);
dm_min = min(dt_m); dm_max = max(dt_m);
sp_min = min(sp); sp_max = max(sp);

w1 = spec.rank.w_actual;
w2 = spec.rank.w_mean;
w3 = spec.rank.w_spread;

for ii = 1:n_ok
    idx = ok_idx(ii);
    ev = arr(idx);
    s_a = norm01_ms(ev.DeltaTN_actual, da_min, da_max);
    s_m = norm01_ms(ev.DeltaTN_mean, dm_min, dm_max);
    sp_v = ifelse_finite_ms(ev.TN_maxmin, 0);
    s_sp = 1 - norm01_ms(sp_v, sp_min, sp_max);
    arr(idx).rank_score = w1 * s_a + w2 * s_m + w3 * s_sp;
end
end

function v = norm01_ms(x, lo, hi)
range = hi - lo;
if range < 1e-12
    v = 0.5;
    return;
end
v = max(0.0, min(1.0, (x - lo) / range));
end

function v = ifelse_finite_ms(x, default_v)
if isfinite(x) && x >= 0
    v = x;
else
    v = default_v;
end
end

%% ====================== Top post + regular-grid baseline ======================

function [opt_full, uni_full, compare_rows, baseline_info] = evaluate_top_with_regular_grid_baseline_ms(top_eval, all_cands, G, spec)
N = numel(top_eval);
opt_full = repmat(empty_eval_struct_ms(spec.stage_count), N, 1);
uni_full = repmat(empty_eval_struct_ms(spec.stage_count), N, 1);
compare_rows = repmat(empty_compare_row_ms(), N, 1);
baseline_info = empty_baseline_info_ms();
if N == 0
    return;
end

for i = 1:N
    opt_full(i) = top_eval(i);
end

best_cand = find_candidate_by_id_ms(all_cands, top_eval(1).candidate_id);
baseline_info.target_total_n = sum(best_cand.n);
[ok_uni, uni_cand, msg] = build_regular_grid_baseline_candidate_ms(best_cand, spec);
baseline_info.control_candidate_count = double(ok_uni);
u_ref = empty_eval_struct_ms(spec.stage_count);

if ok_uni
    eval_try = evaluate_single_candidate_ms(uni_cand, G, spec, true);
    baseline_info = append_baseline_attempt_ms(baseline_info, 'regular_grid_baseline', NaN, ...
        uni_cand, msg, true, eval_try.success, eval_try.message);
    if eval_try.success
        u_ref = eval_try;
        baseline_info.selected_source = 'regular_grid_baseline';
        baseline_info.selected_scale = NaN;
        baseline_info.selected_n = eval_try.n;
        baseline_info.message = 'regular grid baseline FEM success';
    else
        u_ref = eval_try;
        u_ref.message = sprintf('regular grid baseline FEM failed: %s', u_ref.message);
        baseline_info.selected_source = 'none';
        baseline_info.selected_scale = NaN;
        baseline_info.selected_n = u_ref.n;
        baseline_info.message = u_ref.message;
    end
else
    baseline_info = append_baseline_attempt_ms(baseline_info, 'regular_grid_baseline', NaN, ...
        uni_cand, msg, false, false, 'baseline_build_failed');
    u_ref = build_failed_eval_stub_from_candidate_ms(uni_cand, spec.stage_count, ...
        sprintf('regular grid baseline build failed: %s', msg), 'regular_grid_baseline');
    baseline_info.selected_source = 'none';
    baseline_info.selected_scale = NaN;
    baseline_info.selected_n = u_ref.n;
    baseline_info.message = u_ref.message;
end
baseline_info.success = logical(u_ref.success);

for i = 1:N
    uni_full(i) = u_ref;
    compare_rows(i) = build_compare_row_ms(i, opt_full(i), uni_full(i));
    compare_rows(i).baseline_build_success = ok_uni;
    compare_rows(i).baseline_fem_success = logical(u_ref.success);
    compare_rows(i).baseline_failure_reason = char(string(baseline_info.message));
end
end

function info = empty_baseline_info_ms()
info = struct('success', false, 'selected_source', 'none', 'selected_scale', NaN, ...
    'target_total_n', NaN, 'selected_n', zeros(1,0), 'message', '', ...
    'control_candidate_count', 0, ...
    'attempts', repmat(empty_baseline_attempt_row_ms(), 0, 1));
end

function row = empty_baseline_attempt_row_ms()
row = struct('attempt_id', NaN, 'source', '', 'scale', NaN, ...
    'candidate_id', NaN, 'n', '', 'build_message', '', ...
    'eval_ran', false, 'eval_success', false, 'eval_message', '');
end

function info = append_baseline_attempt_ms(info, source, scale, cand, build_msg, eval_ran, eval_success, eval_msg)
row = empty_baseline_attempt_row_ms();
row.attempt_id = numel(info.attempts) + 1;
row.source = char(string(source));
row.scale = scale;
if isstruct(cand) && isfield(cand, 'candidate_id')
    row.candidate_id = cand.candidate_id;
end
if isstruct(cand) && isfield(cand, 'n') && ~isempty(cand.n)
    row.n = vec_to_inline_str_ms(cand.n);
end
row.build_message = char(string(build_msg));
row.eval_ran = logical(eval_ran);
row.eval_success = logical(eval_success);
row.eval_message = char(string(eval_msg));
info.attempts = [info.attempts; row];
end

function rec = build_failed_eval_stub_from_candidate_ms(cand, N, msg, layout_method)
rec = empty_eval_struct_ms(N);
rec.success = false;
rec.message = char(string(msg));
if nargin >= 4 && ~isempty(layout_method)
    rec.layout_method = char(string(layout_method));
end
if ~isstruct(cand)
    return;
end
if isfield(cand, 'candidate_id'), rec.candidate_id = cand.candidate_id; end
if isfield(cand, 'layout_method') && isempty(layout_method), rec.layout_method = cand.layout_method; end
if isfield(cand, 'stage_modes'), rec.stage_modes = cand.stage_modes; end
if isfield(cand, 'stage_methods'), rec.stage_methods = cand.stage_methods; end
if isfield(cand, 'stage_trends'), rec.stage_trends = cand.stage_trends; end
if isfield(cand, 'symmetry_mode'), rec.symmetry_mode = cand.symmetry_mode; end
if isfield(cand, 'edge_pattern_mode'), rec.edge_pattern_mode = cand.edge_pattern_mode; end
if isfield(cand, 's_dense'), rec.s_dense = cand.s_dense; end
if isfield(cand, 's_sparse'), rec.s_sparse = cand.s_sparse; end
if isfield(cand, 'expo'), rec.expo = cand.expo; end
if isfield(cand, 'anis_ratio'), rec.anis_ratio = cand.anis_ratio; end
if isfield(cand, 'gamma'), rec.gamma = cand.gamma; end
if isfield(cand, 'method_anchor_stage'), rec.method_anchor_stage = cand.method_anchor_stage; end
if isfield(cand, 'shape_mode'), rec.shape_mode = cand.shape_mode; end
if isfield(cand, 'stage_anis'), rec.stage_anis = cand.stage_anis; end
if isfield(cand, 'ring_radius_ratio'), rec.ring_radius_ratio = cand.ring_radius_ratio; end
if isfield(cand, 'ring_width_ratio'), rec.ring_width_ratio = cand.ring_width_ratio; end
if isfield(cand, 'band_width_ratio'), rec.band_width_ratio = cand.band_width_ratio; end
if isfield(cand, 'corner_bias'), rec.corner_bias = cand.corner_bias; end
if isfield(cand, 'jitter_ratio'), rec.jitter_ratio = cand.jitter_ratio; end
if isfield(cand, 'jitter_seed'), rec.jitter_seed = cand.jitter_seed; end
if isfield(cand, 'spacing_ratio'), rec.spacing_ratio = cand.spacing_ratio; end
if isfield(cand, 'contrast_score'), rec.contrast_score = cand.contrast_score; end
if isfield(cand, 'lmax_relaxed'), rec.lmax_relaxed = cand.lmax_relaxed; end
if isfield(cand, 'lmax_relax_ratio'), rec.lmax_relax_ratio = cand.lmax_relax_ratio; end
if isfield(cand, 'n'), rec.n = cand.n; end
if isfield(cand, 'ratios'), rec.ratios = cand.ratios; end
if isfield(cand, 'Lx'), rec.Lx = cand.Lx; end
if isfield(cand, 'Ly'), rec.Ly = cand.Ly; end
if isfield(cand, 'cov'), rec.cov = cand.cov; end
if isfield(cand, 'stage_rects'), rec.stage_rects = cand.stage_rects; end
end

function [ok, cand_base, msg] = build_regular_grid_baseline_candidate_ms(cand_ref, spec)
ok = false;
msg = '';
cand_base = baseline_failed_candidate_shell_ms(cand_ref.n, spec.stage_count, 'regular_grid_baseline');
N = spec.stage_count;
if ~isfield(cand_ref, 'n') || numel(cand_ref.n) ~= N || any(~isfinite(cand_ref.n)) || any(cand_ref.n < 1)
    msg = 'invalid reference n for regular grid baseline';
    return;
end
[ok, cand_base, msg] = build_regular_grid_baseline_from_counts_ms(cand_ref, spec, round(cand_ref.n(:).'));
if ok
    msg = sprintf('regular grid baseline geometry ok (n=%s)', vec_to_inline_str_ms(cand_base.n));
end
end

function cand = baseline_failed_candidate_shell_ms(n_vec, N, layout_method)
cand = candidate_record_template_ms(N);
cand.candidate_id = -1;
cand.layout_method = char(string(layout_method));
cand.stage_modes = repmat({'baseline_failed'}, 1, N);
cand.stage_methods = repmat({char(string(layout_method))}, 1, N);
cand.stage_trends = repmat({'neutral'}, 1, N);
cand.symmetry_mode = 'none';
cand.edge_pattern_mode = 'two_sides_to_center';
if numel(n_vec) == N
    cand.n = n_vec(:).';
    cand.ratios = cand.n(2:end) ./ max(cand.n(1:end-1), eps);
end
end

function [ok, cand_out, msg] = build_regular_grid_baseline_from_counts_ms(cand_ref, spec, n_vec)
ok = false;
msg = '';
cand_out = cand_ref;
N = spec.stage_count;
if numel(n_vec) ~= N || any(~isfinite(n_vec)) || any(n_vec < 1)
    msg = 'invalid n vector';
    return;
end

rects_u = cell(1, N);
Lx_u = NaN(1, N);
Ly_u = NaN(1, N);
for k = 1:N
    [ok_k, rects_k, Lx_k, Ly_k] = make_regular_grid_baseline_rects_ms(n_vec(k), spec);
    if ~ok_k
        msg = sprintf('regular grid stage%d layout failed', k);
        return;
    end
    rects_u{k} = rects_k;
    Lx_u(k) = Lx_k;
    Ly_u(k) = Ly_k;
end

Lx = zeros(1, N);
Ly = zeros(1, N);
Lx(N) = Lx_u(N);
Ly(N) = Ly_u(N);
for k = N-1:-1:1
    Lx(k) = max(Lx_u(k), Lx(k+1) + spec.geometry.pyramid_gap_min(k));
    Ly(k) = max(Ly_u(k), Ly(k+1) + spec.geometry.pyramid_gap_min(k));
end
cov = zeros(1, N);
for k = 1:N
    cov(k) = n_vec(k) * spec.fp_w * spec.fp_h / max(Lx(k) * Ly(k), eps);
end
stages_hc = cell(1, N);
for k = 1:N
    stages_hc{k} = struct('rects', rects_u{k}, 'lmax_relax_ratio', 1.0);
end
hc = resolve_hard_constraints_ms(spec);
hc.force_c4_only = false;
[ok_hc, hc_msg] = validate_candidate_hard_constraints_ms('none', stages_hc, Lx, Ly, cov, spec, hc);
if ~ok_hc
    msg = sprintf('regular grid hard-constraint failed: %s', hc_msg);
    return;
end

cand_out.candidate_id = -1;
cand_out.layout_method = 'k25_standard_fixed50_baseline';
cand_out.stage_modes = repmat({'standard_fixed50'}, 1, N);
cand_out.stage_methods = repmat({'k25_standard_fixed50_baseline'}, 1, N);
for k = 1:N
    cand_out.stage_trends{k} = 'neutral';
end
cand_out.symmetry_mode = 'none';
cand_out.edge_pattern_mode = 'two_sides_to_center';
cand_out.n = n_vec(:).';
cand_out.ratios = cand_out.n(2:end) ./ max(cand_out.n(1:end-1), eps);
cand_out.I_opt = NaN;
cand_out.count_rank = NaN;
cand_out.Lx = Lx;
cand_out.Ly = Ly;
cand_out.cov = cov;
cand_out.stage_rects = cell(1, N);
for k = 1:N
    cand_out.stage_rects{k} = rects_u{k};
end
cand_out.method_anchor_stage = NaN;
cand_out.gamma = NaN;
cand_out.lmax_relaxed = false;
cand_out.lmax_relax_ratio = 1.0;
cand_out.s_dense = spec.fp_w * 1.5;
cand_out.s_sparse = cand_out.s_dense;
cand_out.expo = 1.0;
cand_out.anis_ratio = 1.0;
cand_out.spacing_ratio = 1.0;
cand_out.contrast_score = calc_contrast_score_ms(cand_out.s_dense, cand_out.s_sparse, cand_out.expo);
cand_out.mode_prior_score = calc_mode_prior_score_ms(cand_out.stage_trends, get_soft_prior_ms(spec), N);
ok = true;
end

function [ok, rects, Lx_box, Ly_box] = make_regular_grid_baseline_rects_ms(n_target, spec)
ok = false;
rects = zeros(0,4);
Lx_box = NaN;
Ly_box = NaN;
if n_target < 1
    return;
end
gap_x = 0.5 * spec.fp_w;
gap_y = 0.5 * spec.fp_h;
margin_x = 0.5 * spec.fp_w;
margin_y = 0.5 * spec.fp_h;
pitch_x = spec.fp_w + gap_x;
pitch_y = spec.fp_h + gap_y;

% Baseline follows optimize_layout_multistage0411_shared_n regular layout:
% no C4 orbit selection, just a deterministic 50%-gap grid filled from two sides.
[nx, ny] = choose_regular_grid_baseline_shape_ms(n_target, spec.fp_w, spec.fp_h);
pts = fill_regular_grid_two_sides_to_center_ms(nx, ny, n_target, pitch_x, pitch_y);
if size(pts, 1) ~= n_target
    return;
end
rects = centers_to_rects_ms(pts, spec.fp_w, spec.fp_h);

xmin = min(rects(:,1));
xmax = max(rects(:,2));
ymin = min(rects(:,3));
ymax = max(rects(:,4));
Lx_box = (xmax - xmin) + 2 * margin_x;
Ly_box = (ymax - ymin) + 2 * margin_y;
if ~(isfinite(Lx_box) && isfinite(Ly_box) && Lx_box > 0 && Ly_box > 0)
    return;
end
ok = true;
end

function [nx_best, ny_best] = choose_regular_grid_baseline_shape_ms(n_target, fp_w, fp_h)
nx_best = 1;
ny_best = n_target;
best_score = inf;
for nx = 1:n_target
    ny = ceil(n_target / nx);
    Lx = nx * fp_w + max(nx - 1, 0) * 0.5 * fp_w + fp_w;
    Ly = ny * fp_h + max(ny - 1, 0) * 0.5 * fp_h + fp_h;
    aspect_pen = abs(log(max(Lx / max(Ly, eps), eps)));
    waste_pen = (nx * ny - n_target) / max(n_target, 1);
    score = aspect_pen + 0.08 * waste_pen;
    if score < best_score
        best_score = score;
        nx_best = nx;
        ny_best = ny;
    end
end
end

function pts = fill_regular_grid_two_sides_to_center_ms(nx, ny, n_target, pitch_x, pitch_y)
xvals = ((0:nx-1) - (nx-1)/2) * pitch_x;
yvals = ((0:ny-1) - (ny-1)/2) * pitch_y;
row_order = two_sides_to_center_order_ms(ny);
col_order = two_sides_to_center_order_ms(nx);
pts = zeros(n_target, 2);
p = 0;
rem_first = mod(n_target, nx);
row_ptr = 1;
if rem_first > 0
    y = yvals(row_order(row_ptr));
    cols = sort(col_order(1:rem_first));
    for ic = 1:numel(cols)
        p = p + 1;
        pts(p,:) = [xvals(cols(ic)), y];
    end
    row_ptr = row_ptr + 1;
end
while p < n_target && row_ptr <= numel(row_order)
    y = yvals(row_order(row_ptr));
    take = min(nx, n_target - p);
    cols = sort(col_order(1:take));
    for ic = 1:take
        p = p + 1;
        pts(p,:) = [xvals(cols(ic)), y];
    end
    row_ptr = row_ptr + 1;
end
if p < n_target
    pts = pts(1:p,:);
else
    pts = sortrows(pts, [2 1]);
end
end

function order = two_sides_to_center_order_ms(n)
left = 1;
right = n;
order = zeros(1, n);
p = 0;
while left <= right
    p = p + 1;
    order(p) = left;
    left = left + 1;
    if right >= left
        p = p + 1;
        order(p) = right;
        right = right - 1;
    end
end
order = order(1:p);
end

function row = empty_compare_row_ms()
row = struct('rank', NaN, 'candidate_id', NaN, ...
    'opt_success', false, 'uni_success', false, ...
    'baseline_candidate_id', -1, ...
    'baseline_build_success', false, 'baseline_fem_success', false, 'baseline_failure_reason', '', ...
    'opt_DeltaTN_mean', NaN, 'uni_DeltaTN_mean', NaN, ...
    'opt_DeltaTN_actual', NaN, 'uni_DeltaTN_actual', NaN, ...
    'd_DeltaTN_mean', NaN, 'd_DeltaTN_actual', NaN);
end

function row = build_compare_row_ms(rank_idx, opt_rec, uni_rec)
row = empty_compare_row_ms();
row.rank = rank_idx;
row.candidate_id = opt_rec.candidate_id;
row.opt_success = opt_rec.success;
row.uni_success = uni_rec.success;
if opt_rec.success
    row.opt_DeltaTN_mean = opt_rec.DeltaTN_mean;
    row.opt_DeltaTN_actual = opt_rec.DeltaTN_actual;
end
if uni_rec.success
    row.baseline_candidate_id = uni_rec.candidate_id;
    row.uni_DeltaTN_mean = uni_rec.DeltaTN_mean;
    row.uni_DeltaTN_actual = uni_rec.DeltaTN_actual;
end
if opt_rec.success && uni_rec.success
    row.d_DeltaTN_mean = opt_rec.DeltaTN_mean - uni_rec.DeltaTN_mean;
    row.d_DeltaTN_actual = opt_rec.DeltaTN_actual - uni_rec.DeltaTN_actual;
end
end

function [opt_sorted, uni_sorted, compare_sorted, top_eval_sorted] = ...
    resort_top_post_records_ms(top_opt_full, top_uni_full, top_compare, top_eval_post, spec)
if nargin < 5
    spec = struct();
end
opt_sorted = top_opt_full;
uni_sorted = top_uni_full;
compare_sorted = top_compare;
top_eval_sorted = top_eval_post;
if isempty(top_opt_full)
    return;
end
opt_sorted = sort_candidates_by_deltaTN_ms(top_opt_full, spec);
cid_sorted = [opt_sorted.candidate_id].';
N = numel(cid_sorted);
uni_sorted = repmat(empty_eval_struct_ms(spec.stage_count), N, 1);
top_eval_sorted = repmat(empty_eval_struct_ms(spec.stage_count), N, 1);
for i = 1:N
    j = find([top_opt_full.candidate_id].' == cid_sorted(i), 1, 'first');
    if isempty(j)
        continue;
    end
    if numel(top_uni_full) >= j
        uni_sorted(i) = top_uni_full(j);
    end
    if numel(top_eval_post) >= j
        top_eval_sorted(i) = top_eval_post(j);
    end
end
compare_sorted = repmat(empty_compare_row_ms(), N, 1);
for i = 1:N
    compare_sorted(i) = build_compare_row_ms(i, opt_sorted(i), uni_sorted(i));
end
end

%% ====================== Output ======================

function write_outputs_ms(spec, count_solution, count_solutions, all_eval, coarse_valid, final_valid, top_eval, top_eval_post, ...
    best_cand, best_full, top_opt_full, top_uni_full, top_compare, baseline_info, current_calib, eval_counts)
out_dir = spec.output.output_dir;
mkdir_if_needed_ms(out_dir);

    save(fullfile(out_dir, 'results_data.mat'), 'spec', 'count_solution', 'count_solutions', 'all_eval', 'coarse_valid', 'final_valid', ...
    'top_eval', 'top_eval_post', 'best_cand', 'best_full', 'top_opt_full', 'top_uni_full', 'top_compare', ...
    'baseline_info', 'current_calib', 'eval_counts');

if nargin < 16 || ~isstruct(eval_counts)
    eval_counts = struct('coarse_valid_count', numel(coarse_valid), ...
        'final_valid_count', numel(final_valid));
end
if ~isfield(eval_counts, 'candidate_batch_jobs_before'), eval_counts.candidate_batch_jobs_before = NaN; end
if ~isfield(eval_counts, 'candidate_batch_jobs_after'), eval_counts.candidate_batch_jobs_after = NaN; end
if ~isfield(eval_counts, 'jobs_after_prune'), eval_counts.jobs_after_prune = eval_counts.candidate_batch_jobs_after; end
if ~isfield(eval_counts, 'geometry_candidates_before_batch'), eval_counts.geometry_candidates_before_batch = NaN; end
if ~isfield(eval_counts, 'geometry_candidates_after_batch'), eval_counts.geometry_candidates_after_batch = NaN; end
if ~isfield(eval_counts, 'coarse_valid_count'), eval_counts.coarse_valid_count = numel(coarse_valid); end
if ~isfield(eval_counts, 'final_valid_count'), eval_counts.final_valid_count = numel(final_valid); end
if ~isfield(eval_counts, 'supplement_info'), eval_counts.supplement_info = empty_supplement_info_ms(); end

coarse_valid = coarse_valid([coarse_valid.success]);
final_valid = final_valid([final_valid.success]);
assert_rank_scores_finite_ms(coarse_valid, 'coarse_valid');
assert_rank_scores_finite_ms(final_valid, 'final_valid');
assert_rank_scores_finite_ms(top_eval, 'top_eval');
assert_rank_scores_finite_ms(top_eval_post, 'top_eval_post');

if ~isempty(coarse_valid)
    tbl_all = struct_to_table_rows_ms(strip_large_fields_ms(coarse_valid, spec.stage_count));
    writetable(tbl_all, fullfile(out_dir, 'all_eval_valid_coarse.csv'));
end
if ~isempty(final_valid)
    tbl_final = struct_to_table_rows_ms(strip_large_fields_ms(final_valid, spec.stage_count));
    writetable(tbl_final, fullfile(out_dir, 'all_eval_valid_final.csv'));
end
write_failed_eval_csv_ms(spec, all_eval, fullfile(out_dir, 'all_eval_failed.csv'));
if ~isempty(top_compare)
    writetable(struct_to_table_rows_ms(top_compare), fullfile(out_dir, 'baseline_compare_summary.csv'));
end
write_eval_failure_summary_csv_ms(all_eval, fullfile(out_dir, 'eval_failure_summary.csv'));
write_candidate_batch_composition_csv_ms(all_eval, fullfile(out_dir, 'candidate_batch_composition.csv'));
write_candidate_repro_manifest_ms(spec, count_solution, top_eval_post, top_opt_full, top_uni_full, ...
    fullfile(out_dir, 'candidate_repro_manifest.csv'));
top_bottom_info = write_top_bottom_candidate_reports_ms(spec, count_solution, final_valid, all_eval, ...
    fullfile(out_dir, 'top10_candidates.csv'), fullfile(out_dir, 'bottom10_candidates.csv'));

stage_metrics_rows = build_stage_metrics_summary_rows_ms(top_opt_full, spec);
if ~isempty(stage_metrics_rows)
    writetable(struct_to_table_rows_ms(stage_metrics_rows), fullfile(out_dir, 'stage_metrics_summary.csv'));
end
baseline_stage_metrics_rows = repmat(empty_stage_metrics_row_ms(), 0, 1);
if ~isempty(top_uni_full)
    uni_ref = top_uni_full(1);
    if isstruct(uni_ref) && has_valid_candidate_geometry_ms(uni_ref, spec.stage_count)
        baseline_stage_metrics_rows = build_stage_metrics_summary_rows_ms(uni_ref, spec);
        if ~isempty(baseline_stage_metrics_rows)
            writetable(struct_to_table_rows_ms(baseline_stage_metrics_rows), fullfile(out_dir, 'baseline_stage_metrics_summary.csv'));
        end
    end
end
write_reasonability_report_ms(fullfile(out_dir, 'reasonability_report.txt'), spec, best_full, stage_metrics_rows);

plot_axis_diag_rows = repmat(empty_plot_axis_diag_row_ms(), 0, 1);
baseline_plot_axis_diag_rows = repmat(empty_plot_axis_diag_row_ms(), 0, 1);
plot_cfg = resolve_plot_cfg_ms(spec);

% Generate detailed plots for the ranked top candidates without expanding post/baseline work.
top_plot_recs = top_eval(:);
top_plot_n = min(spec.output.top_plotK, numel(top_plot_recs));
if top_plot_n > 0
    top_plot_recs = top_plot_recs(1:top_plot_n);
    % Compute unified color axis for each stage across all top candidates
    N = spec.stage_count;
    caxis_by_stage = NaN(N, 2);
    for k = 1:N
        all_temps = [];
        for i = 1:numel(top_plot_recs)
            if ~isempty(top_plot_recs(i).Tfields) && numel(top_plot_recs(i).Tfields) >= k && ...
                    ~isempty(top_plot_recs(i).Tfields{k})
                all_temps = [all_temps; top_plot_recs(i).Tfields{k}(:)]; %#ok<AGROW>
            end
        end
        if ~isempty(all_temps)
            caxis_by_stage(k, :) = [min(all_temps), max(all_temps)];
        end
    end

    plot_axis_diag_cells = cell(top_plot_n, 1);
    parfor i = 1:top_plot_n
        top_dir = fullfile(out_dir, sprintf('top%d', i));
        mkdir_if_needed_ms(top_dir);
        layout_info = empty_plot_axis_info_ms(N);
        temp_info = empty_plot_axis_info_ms(N);
        diag_tag = '';

        if plot_cfg.save_overview
            % Save overview layout/temperature plots
            layout_info = save_layout_plot_ms(top_plot_recs(i), fullfile(top_dir, 'layout.png'), spec);
            temp_info = save_temperature_plot_ms(top_plot_recs(i), fullfile(top_dir, 'temperature.png'), caxis_by_stage, spec);
            diag_tag = 'overview';
        end
        if plot_cfg.save_stage_separate
            % Save per-stage separate layout/temperature plots
            layout_sep_info = save_layout_stage_plots_ms(top_plot_recs(i), top_dir, spec);
            temp_sep_info = save_temperature_stage_plots_ms(top_plot_recs(i), top_dir, caxis_by_stage, spec);
            if ~plot_cfg.save_overview
                layout_info = layout_sep_info;
                temp_info = temp_sep_info;
                diag_tag = 'separate';
            end
        end
        if any(isfinite(layout_info.target_axes_mm(:))) || any(isfinite(temp_info.target_axes_mm(:)))
            plot_axis_diag_cells{i} = build_plot_axis_diag_rows_ms(i, layout_info, temp_info, diag_tag);
        end
    end
    for i = 1:top_plot_n
        if ~isempty(plot_axis_diag_cells{i})
            plot_axis_diag_rows = [plot_axis_diag_rows; plot_axis_diag_cells{i}]; %#ok<AGROW>
        end
    end
    if ~isempty(plot_axis_diag_rows)
        writetable(struct_to_table_rows_ms(plot_axis_diag_rows), fullfile(out_dir, 'plot_axis_diagnostics.csv'));
    end

    if ~isempty(top_uni_full)
        baseline_rec = top_uni_full(1);
        baseline_dir = fullfile(out_dir, 'baseline');
        mkdir_if_needed_ms(baseline_dir);
        layout_b = empty_plot_axis_info_ms(N);
        temp_b = empty_plot_axis_info_ms(N);
        baseline_diag_tag = '';
        has_geom = has_valid_candidate_geometry_ms(baseline_rec, spec.stage_count);
        if has_geom
            write_stage_layout_rects_csv_ms(baseline_rec, fullfile(baseline_dir, 'layout_stage_rects.csv'), spec.stage_count);
            if plot_cfg.save_overview
                layout_b = save_layout_plot_ms(baseline_rec, fullfile(baseline_dir, 'layout.png'), spec);
                baseline_diag_tag = 'baseline_overview';
            end
            if plot_cfg.save_stage_separate
                layout_bs = save_layout_stage_plots_ms(baseline_rec, baseline_dir, spec);
                if ~plot_cfg.save_overview
                    layout_b = layout_bs;
                    baseline_diag_tag = 'baseline_separate';
                end
            end
        else
            fid_g = fopen(fullfile(baseline_dir, 'geometry_unavailable.txt'), 'w');
            if fid_g > 0
                fprintf(fid_g, 'Baseline geometry unavailable.\n');
                fprintf(fid_g, 'reason=%s\n', baseline_rec.message);
                fclose(fid_g);
            end
        end
        if baseline_rec.success
            if plot_cfg.save_overview
                temp_b = save_temperature_plot_ms(baseline_rec, fullfile(baseline_dir, 'temperature.png'), caxis_by_stage, spec);
            end
            if plot_cfg.save_stage_separate
                temp_bs = save_temperature_stage_plots_ms(baseline_rec, baseline_dir, caxis_by_stage, spec);
                if ~plot_cfg.save_overview
                    temp_b = temp_bs;
                end
            end
            if isempty(baseline_diag_tag)
                baseline_diag_tag = 'baseline_temp_only';
            end
        else
            fid_b = fopen(fullfile(baseline_dir, 'temperature_unavailable.txt'), 'w');
            if fid_b > 0
                fprintf(fid_b, 'Baseline FEM failed, temperature fields unavailable.\n');
                fprintf(fid_b, 'reason=%s\n', baseline_rec.message);
                fclose(fid_b);
            end
            if isempty(baseline_diag_tag)
                baseline_diag_tag = 'baseline_failed';
            end
        end
        if any(isfinite(layout_b.target_axes_mm(:))) || any(isfinite(temp_b.target_axes_mm(:)))
            baseline_plot_axis_diag_rows = build_plot_axis_diag_rows_ms(1, layout_b, temp_b, baseline_diag_tag);
            writetable(struct_to_table_rows_ms(baseline_plot_axis_diag_rows), fullfile(out_dir, 'baseline_plot_axis_diagnostics.csv'));
        end
    end
end

write_baseline_status_txt_ms(fullfile(out_dir, 'baseline_status.txt'), baseline_info, top_uni_full, spec.stage_count);
write_llm_analysis_packet_ms(fullfile(out_dir, 'llm_analysis_packet.md'), spec, count_solution, ...
    top_bottom_info, eval_counts, baseline_info);
write_llm_analysis_packet_ms(fullfile(out_dir, 'all_candidates_analysis_packet.md'), spec, count_solution, ...
    top_bottom_info, eval_counts, baseline_info);

fid = fopen(fullfile(out_dir, 'summary.txt'), 'w');
if fid > 0
    fprintf(fid, 'optimize_layout_multistage summary\n');
    fprintf(fid, 'entry=%s\n', mfilename);
    fprintf(fid, 'code_version_stamp=%s\n', code_version_stamp_ms());
    fprintf(fid, 'stage_count=%d\n', spec.stage_count);
    fprintf(fid, 'fixed_n_enable=1\n');
    fprintf(fid, 'plate_k_inplane=%.6g\n', spec.plate_k_inplane);
    fprintf(fid, 'ceramic_plate_k_inplane_W_mK=%.12g\n', spec.plate_k_inplane);
    fprintf(fid, 'particle_width_mm=%.12g\n', spec.fp_w * 1e3);
    fprintf(fid, 'particle_length_mm=%.12g\n', spec.fp_h * 1e3);
    fprintf(fid, 'particle_height_stage_mm=%s\n', vec_to_inline_str_ms(spec.fp_t_stage * 1e3));
    fprintf(fid, 'fp_w_mm=%.12g\n', spec.fp_w * 1e3);
    fprintf(fid, 'fp_h_mm=%.12g\n', spec.fp_h * 1e3);
    fprintf(fid, 'fp_t_stage_mm=%s\n', vec_to_inline_str_ms(spec.fp_t_stage * 1e3));
    fprintf(fid, 'particle_area_mm2=%.12g\n', spec.fp_w * spec.fp_h * 1e6);
    fprintf(fid, 'particle_volume_stage_mm3=%s\n', vec_to_inline_str_ms(spec.fp_w * spec.fp_h * spec.fp_t_stage * 1e9));
    fprintf(fid, 'current_I_init_A=%.12g\n', spec.current.I_init);
    fprintf(fid, 'candidate_batch_size=%d\n', spec.candidate_batch.size);
    fprintf(fid, 'topK=%d\n', spec.output.topK);
    fprintf(fid, 'top_postK=%d\n', spec.output.top_postK);
    fprintf(fid, 'top_plotK=%d\n', spec.output.top_plotK);
    fprintf(fid, 'top_bottom_n=%d\n', spec.candidate_batch.top_bottom_n);
    fprintf(fid, 'jobs_after_prune=%d\n', eval_counts.jobs_after_prune);
    fprintf(fid, 'geometry_candidates_before_batch=%d\n', eval_counts.geometry_candidates_before_batch);
    fprintf(fid, 'geometry_candidates_after_batch=%d\n', eval_counts.geometry_candidates_after_batch);
    sinfo = eval_counts.supplement_info;
    fprintf(fid, 'candidate_supplement_enabled=%d\n', logical(sinfo.enabled));
    fprintf(fid, 'candidate_supplement_target=%d\n', sinfo.target);
    fprintf(fid, 'candidate_supplement_attempted_jobs=%d\n', sinfo.attempted_jobs);
    fprintf(fid, 'candidate_supplement_max_attempts=%d\n', sinfo.max_attempts);
    fprintf(fid, 'candidate_supplement_initial_count=%d\n', sinfo.initial_count);
    fprintf(fid, 'candidate_supplement_final_count=%d\n', sinfo.final_count);
    fprintf(fid, 'candidate_supplement_stop_reason=%s\n', char(string(sinfo.stop_reason)));
    fprintf(fid, 'candidate_supplement_quota_c4_c2_shape=[%d %d %d]\n', ...
        sinfo.quota_c4_main, sinfo.quota_c2, sinfo.quota_shape);
    fprintf(fid, 'candidate_supplement_final_c4_c2_shape=[%d %d %d]\n', ...
        sinfo.final_c4_main, sinfo.final_c2, sinfo.final_shape);
    fprintf(fid, 'shape_explore_enable=%d\n', logical(spec.shape_explore.enable));
    fprintf(fid, 'shape_explore_modes=%s\n', strjoin(spec.shape_explore.extra_stage_modes, ','));
    fprintf(fid, 'shape_explore_stage_anis=%s\n', vec_to_inline_str_ms(spec.shape_explore.anis_ratio_stage_list));
    fprintf(fid, 'shape_explore_jitter=%s\n', vec_to_inline_str_ms(spec.shape_explore.jitter_ratio_list));
    if isfield(spec.shape_explore, 'stage_mode_templates') && ~isempty(spec.shape_explore.stage_mode_templates)
        fprintf(fid, 'shape_explore_stage_mode_templates=%s\n', stage_templates_to_inline_str_ms(spec.shape_explore.stage_mode_templates));
    end
    fprintf(fid, 'n=%s\n', vec_to_inline_str_ms(count_solution.n));
    fprintf(fid, 'particle_count_n=%s\n', vec_to_inline_str_ms(count_solution.n));
    fprintf(fid, 'particle_count_total=%d\n', sum(count_solution.n));
    fprintf(fid, 'ratios=%s\n', vec_to_inline_str_ms(count_solution.ratios));
    fprintf(fid, 'I_opt=%.6f A\n', count_solution.I_opt);
    fprintf(fid, 'current_I_eval_A=%.12g\n', count_solution.I_opt);
    fprintf(fid, 'fixed_n_0d_DeltaT=%.6f K\n', count_solution.DeltaT_0D);
    fprintf(fid, 'DeltaT_target=%.3f K, Qc_target_last=%.3f W\n', spec.targets.DeltaT_target, spec.targets.Qc_target_last);
    fprintf(fid, 'L_max_mm=%s\n', vec_to_inline_str_ms(spec.geometry.L_max_mm));
    fprintf(fid, 'coverage_min=%s\n', vec_to_inline_str_ms(spec.geometry.coverage_min));
    fprintf(fid, 'pyramid_gap_min_mm=%s\n', vec_to_inline_str_ms(spec.geometry.pyramid_gap_min_mm));
    fprintf(fid, 'min_edge_gap_mm=%.6f\n', spec.geometry.min_edge_gap_mm);
    if is_z_path_enabled_ms(spec)
        fprintf(fid, 'z_path.enable=1\n');
        fprintf(fid, 'z_path.k_interfaces=%s\n', vec_to_inline_str_ms(spec.z_path.k_interfaces));
        fprintf(fid, 'z_path.t_interface_effs=%s\n', vec_to_inline_str_ms(spec.z_path.t_interface_effs));
        fprintf(fid, 'z_path.Rc_interfaces=%s\n', vec_to_inline_str_ms(spec.z_path.Rc_interfaces));
        fprintf(fid, 'z_path.k_sink=%.6g, t_sink_eff=%.6g, Rc_sink=%.6g\n', ...
            spec.z_path.k_sink, spec.z_path.t_sink_eff, spec.z_path.Rc_sink);
    else
        fprintf(fid, 'z_path=disabled\n');
    end
    fprintf(fid, 'soft_prior(stage trends)=%s, fallback=%d, warmup_topn=%d\n', ...
        strjoin(spec.soft_prior.stage_trend_prefer, '/'), logical(spec.soft_prior.is_fallback), spec.soft_prior.warmup_topn);
    fprintf(fid, 'plot.use_global_axis=%d, show_mesh_edges=%d, view_mode=%s, axis_margin_ratio=%.3f\n', ...
        spec.output.plot.use_global_axis, spec.output.plot.show_mesh_edges, ...
        spec.output.plot.view_mode, spec.output.plot.axis_margin_ratio);
    fprintf(fid, 'plot.save_overview=%d, save_stage_separate=%d, annotate_substrate_dims=%d, separate_fig_size_px=%s\n', ...
        spec.output.plot.save_overview, spec.output.plot.save_stage_separate, ...
        spec.output.plot.annotate_substrate_dims, vec_to_inline_str_ms(spec.output.plot.separate_fig_size_px));
    fprintf(fid, 'coarse_valid_count=%d\n', eval_counts.coarse_valid_count);
    fprintf(fid, 'final_valid_count=%d\n', eval_counts.final_valid_count);
    if ~isempty(plot_axis_diag_rows)
        fprintf(fid, 'plot_axis_diagnostics_csv=%s\n', 'plot_axis_diagnostics.csv');
    end
    if ~isempty(stage_metrics_rows)
        fprintf(fid, 'stage_metrics_summary_csv=%s\n', 'stage_metrics_summary.csv');
    end
    if ~isempty(baseline_stage_metrics_rows)
        fprintf(fid, 'baseline_stage_metrics_summary_csv=%s\n', 'baseline_stage_metrics_summary.csv');
    end
    if ~isempty(baseline_plot_axis_diag_rows)
        fprintf(fid, 'baseline_plot_axis_diagnostics_csv=%s\n', 'baseline_plot_axis_diagnostics.csv');
    end
    if ~isempty(top_compare)
        fprintf(fid, 'baseline_compare_summary_csv=%s\n', 'baseline_compare_summary.csv');
    end
    fprintf(fid, 'eval_failure_summary_csv=%s\n', 'eval_failure_summary.csv');
    fprintf(fid, 'candidate_batch_composition_csv=%s\n', 'candidate_batch_composition.csv');
    fprintf(fid, 'candidate_repro_manifest_csv=%s\n', 'candidate_repro_manifest.csv');
    fprintf(fid, 'top10_candidates_csv=%s\n', 'top10_candidates.csv');
    fprintf(fid, 'bottom10_candidates_csv=%s\n', 'bottom10_candidates.csv');
    fprintf(fid, 'llm_analysis_packet_md=%s\n', 'llm_analysis_packet.md');
    fprintf(fid, 'all_candidates_analysis_packet_md=%s\n', 'all_candidates_analysis_packet.md');
    fprintf(fid, 'baseline_status_txt=%s\n', 'baseline_status.txt');
    if isstruct(baseline_info)
        fprintf(fid, 'baseline_success=%d, baseline_source=%s, baseline_scale=%.6f\n', ...
            logical(baseline_info.success), char(string(baseline_info.selected_source)), baseline_info.selected_scale);
    end
    fprintf(fid, 'reasonability_report_txt=%s\n', 'reasonability_report.txt');
    if isstruct(current_calib) && isfield(current_calib, 'has_solution')
        fprintf(fid, 'Current calibration: enabled=%d, solved=%d, I0=%.6f A, Ical=%.6f A, shift=%.6f A\n', ...
            logical(current_calib.enabled), logical(current_calib.has_solution), ...
            current_calib.I0, current_calib.I_cal, current_calib.I_shift);
    end
    if best_full.success
        fprintf(fid, 'Top1 DeltaTN_mean=%.6f K\n', best_full.DeltaTN_mean);
        fprintf(fid, 'Top1 DeltaTN_actual=%.6f K\n', best_full.DeltaTN_actual);
    end
    fclose(fid);
end
end

function info = write_top_bottom_candidate_reports_ms(spec, count_solution, final_valid, all_eval, top_csv, bottom_csv)
topN = max(1, round(spec.candidate_batch.top_bottom_n));
info = struct('top_count', 0, 'bottom_count', 0, 'valid_count', 0, 'failed_count', 0, ...
    'top_best_delta', NaN, 'top_worst_delta', NaN, 'bottom_best_delta', NaN, ...
    'bottom_worst_delta', NaN, 'top_bottom_gap_K', NaN, 'converged_like', false);
if nargin < 5
    return;
end
if isempty(final_valid)
    valid = repmat(empty_eval_struct_ms(spec.stage_count), 0, 1);
else
    valid = final_valid([final_valid.success]);
end
info.valid_count = numel(valid);
if ~isempty(all_eval) && isfield(all_eval, 'success')
    info.failed_count = sum(~[all_eval.success]);
end
if isempty(valid)
    return;
end
valid = sort_candidates_by_deltaTN_ms(valid, spec);
top_take = min(topN, numel(valid));
top_recs = valid(1:top_take);
bottom_recs = valid(max(1, numel(valid) - top_take + 1):numel(valid));
bottom_recs = flipud(bottom_recs(:));
top_rows = build_candidate_rank_rows_ms(top_recs, 'top', count_solution, spec);
bottom_rows = build_candidate_rank_rows_ms(bottom_recs, 'bottom', count_solution, spec);
try
    if ~isempty(top_rows)
        writetable(struct_to_table_rows_ms(top_rows), top_csv);
    end
catch ME
    warning('Failed to write %s: %s', top_csv, ME.message);
end
try
    if ~isempty(bottom_rows)
        writetable(struct_to_table_rows_ms(bottom_rows), bottom_csv);
    end
catch ME
    warning('Failed to write %s: %s', bottom_csv, ME.message);
end
info.top_count = numel(top_recs);
info.bottom_count = numel(bottom_recs);
info.top_best_delta = top_recs(1).DeltaTN_actual;
info.top_worst_delta = top_recs(end).DeltaTN_actual;
info.bottom_best_delta = bottom_recs(1).DeltaTN_actual;
info.bottom_worst_delta = bottom_recs(end).DeltaTN_actual;
info.top_bottom_gap_K = info.top_best_delta - info.bottom_worst_delta;
info.converged_like = isfinite(info.top_bottom_gap_K) && info.top_bottom_gap_K <= spec.candidate_batch.convergence_gap_K;
end

function rows = build_candidate_rank_rows_ms(recs, source, count_solution, spec)
rows = repmat(empty_repro_manifest_row_ms(), 0, 1);
for i = 1:numel(recs)
    rows = [rows; build_repro_manifest_row_ms(i, source, spec.seed, count_solution, recs(i), recs(i), spec)]; %#ok<AGROW>
end
end

function write_llm_analysis_packet_ms(out_md, spec, count_solution, top_bottom_info, eval_counts, baseline_info)
fid = fopen(out_md, 'w');
if fid < 0
    return;
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# Fixed-n Layout Search Analysis Packet\n\n');
fprintf(fid, '## Run Goal\n');
fprintf(fid, '- Find the layout with maximum `DeltaTN_actual` under fixed particle counts.\n');
fprintf(fid, '- Fixed n: `%s`\n', vec_to_inline_str_ms(count_solution.n));
fprintf(fid, '- Current used for batch evaluation starts from `spec.current.I_init = %.6f A`.\n', spec.current.I_init);
fprintf(fid, '- Candidate batch target: `%d` complete multistage layouts.\n\n', spec.candidate_batch.size);
fprintf(fid, '## Parameter Space\n');
fprintf(fid, '- layout_methods: `%s`\n', strjoin(spec.layout_methods, ','));
fprintf(fid, '- mode_list: `%s`\n', strjoin(spec.mode_list, ','));
fprintf(fid, '- s_dense_mm: `%s`\n', vec_to_inline_str_ms(spec.coarse_s_dense_list * 1e3));
fprintf(fid, '- s_sparse_mm: `%s`\n', vec_to_inline_str_ms(spec.coarse_s_sparse_list * 1e3));
fprintf(fid, '- expo: `%s`\n', vec_to_inline_str_ms(spec.coarse_expo_list));
fprintf(fid, '- anis_ratio: `%s`\n', vec_to_inline_str_ms(spec.coarse_anis_ratio_list));
fprintf(fid, '- gamma: `%s`\n', vec_to_inline_str_ms(spec.gamma_list));
fprintf(fid, '- interstage_align_weights: `%s`\n\n', vec_to_inline_str_ms(spec.interstage.align_weights));
fprintf(fid, '- shape_explore_enable: `%d`\n', logical(spec.shape_explore.enable));
fprintf(fid, '- shape_modes: `%s`\n', strjoin(spec.shape_explore.extra_stage_modes, ','));
fprintf(fid, '- shape_stage_anis: `%s`\n', vec_to_inline_str_ms(spec.shape_explore.anis_ratio_stage_list));
fprintf(fid, '- ring_radius_ratio: `%s`\n', vec_to_inline_str_ms(spec.shape_explore.ring_radius_ratio_list));
fprintf(fid, '- ring_width_ratio: `%s`\n', vec_to_inline_str_ms(spec.shape_explore.ring_width_ratio_list));
fprintf(fid, '- band_width_ratio: `%s`\n', vec_to_inline_str_ms(spec.shape_explore.band_width_ratio_list));
fprintf(fid, '- corner_bias: `%s`\n', vec_to_inline_str_ms(spec.shape_explore.corner_bias_list));
fprintf(fid, '- jitter_ratio: `%s`\n\n', vec_to_inline_str_ms(spec.shape_explore.jitter_ratio_list));
if isfield(spec.shape_explore, 'stage_mode_templates') && ~isempty(spec.shape_explore.stage_mode_templates)
    fprintf(fid, '- stage_mode_templates: `%s`\n\n', stage_templates_to_inline_str_ms(spec.shape_explore.stage_mode_templates));
end
fprintf(fid, '## Result Files\n');
fprintf(fid, '- Primary full-candidate analysis table: `all_eval_valid_final.csv`\n');
fprintf(fid, '- Coarse-stage valid candidates: `all_eval_valid_coarse.csv`\n');
fprintf(fid, '- Top candidates for extreme-sample reference only: `top10_candidates.csv`\n');
fprintf(fid, '- Bottom candidates for extreme-sample reference only: `bottom10_candidates.csv`\n');
fprintf(fid, '- Failure summary: `eval_failure_summary.csv`\n');
fprintf(fid, '- Candidate group composition: `candidate_batch_composition.csv`\n');
fprintf(fid, '- Reproduction manifest: `candidate_repro_manifest.csv`\n\n');
fprintf(fid, '## Batch Summary\n');
fprintf(fid, '- jobs_after_prune: `%d`\n', eval_counts.jobs_after_prune);
fprintf(fid, '- geometry_candidates_before_batch: `%d`\n', eval_counts.geometry_candidates_before_batch);
fprintf(fid, '- geometry_candidates_after_batch: `%d`\n', eval_counts.geometry_candidates_after_batch);
if isfield(eval_counts, 'supplement_info') && isstruct(eval_counts.supplement_info)
    sinfo = eval_counts.supplement_info;
    fprintf(fid, '- candidate_supplement_enabled: `%d`\n', logical(sinfo.enabled));
    fprintf(fid, '- candidate_supplement_attempted_jobs: `%d/%d`\n', sinfo.attempted_jobs, sinfo.max_attempts);
    fprintf(fid, '- candidate_supplement_final_c4_c2_shape: `[%d %d %d]`\n', ...
        sinfo.final_c4_main, sinfo.final_c2, sinfo.final_shape);
    fprintf(fid, '- candidate_supplement_stop_reason: `%s`\n', char(string(sinfo.stop_reason)));
end
fprintf(fid, '- valid_count: `%d`\n', top_bottom_info.valid_count);
fprintf(fid, '- failed_count: `%d`\n', top_bottom_info.failed_count);
fprintf(fid, '- top_best_DeltaTN_actual_K: `%.6f`\n', top_bottom_info.top_best_delta);
fprintf(fid, '- top_worst_DeltaTN_actual_K: `%.6f`\n', top_bottom_info.top_worst_delta);
fprintf(fid, '- bottom_worst_DeltaTN_actual_K: `%.6f`\n', top_bottom_info.bottom_worst_delta);
fprintf(fid, '- top_to_bottom_gap_K: `%.6f`\n', top_bottom_info.top_bottom_gap_K);
fprintf(fid, '- convergence_gap_threshold_K: `%.6f`\n', spec.candidate_batch.convergence_gap_K);
fprintf(fid, '- approximate_parameter_insensitive: `%d`\n\n', logical(top_bottom_info.converged_like));
if isstruct(baseline_info)
    fprintf(fid, '## Baseline\n');
    fprintf(fid, '- baseline_success: `%d`\n', logical(baseline_info.success));
    fprintf(fid, '- baseline_source: `%s`\n\n', char(string(baseline_info.selected_source)));
end
fprintf(fid, '## Prompt For Next-Round Analysis\n');
fprintf(fid, 'Use `all_eval_valid_final.csv` as the primary dataset and analyze every valid FEM candidate, not only top/bottom rows. Use `top10_candidates.csv` and `bottom10_candidates.csv` only as extreme-sample references. Identify which layout parameters consistently improve `DeltaTN_actual`, which shape parameters create useful geometry diversity, which parameters correlate with failed or weak candidates, and propose the next 64-candidate parameter batch.\n\n');
fprintf(fid, '## Fixed Operating Constraints For Next Layout Round\n');
fprintf(fid, '- Keep `spec.fixed_n = %s` unless particle-count changes are listed separately as post-analysis hypotheses.\n', vec_to_inline_str_ms(count_solution.n));
fprintf(fid, '- Keep `spec.current.I_init = %.6f A` and do not mix current optimization into layout conclusions.\n', spec.current.I_init);
fprintf(fid, '- Keep `spec.current_calib.enable = false` for the next layout search.\n');
fprintf(fid, '- `all_eval_valid_final.csv` is the main AI tuning input; top5 plots remain visual diagnostics only.\n\n');
fprintf(fid, '## Next-Round Parameter Template\n');
fprintf(fid, '```matlab\n');
fprintf(fid, 'spec.fixed_n = %s;\n', vec_to_inline_str_ms(count_solution.n));
fprintf(fid, 'spec.current.I_init = %.6f;\n', spec.current.I_init);
fprintf(fid, 'spec.current_calib.enable = false;\n');
fprintf(fid, 'spec.candidate_batch.size = 64;\n');
fprintf(fid, 'spec.layout_methods = {...};\n');
fprintf(fid, 'spec.coarse_s_dense_list = [...] * 1e-3;\n');
fprintf(fid, 'spec.coarse_s_sparse_list = [...] * 1e-3;\n');
fprintf(fid, 'spec.coarse_expo_list = [...];\n');
fprintf(fid, 'spec.coarse_anis_ratio_list = [...];\n');
fprintf(fid, 'spec.gamma_list = [...];\n');
fprintf(fid, 'spec.interstage.align_weights = [...];\n');
fprintf(fid, 'spec.shape_explore.enable = true;\n');
fprintf(fid, 'spec.shape_explore.extra_stage_modes = {...};\n');
fprintf(fid, 'spec.shape_explore.anis_ratio_stage_list = [...];\n');
fprintf(fid, 'spec.shape_explore.ring_radius_ratio_list = [...];\n');
fprintf(fid, 'spec.shape_explore.band_width_ratio_list = [...];\n');
fprintf(fid, 'spec.shape_explore.jitter_ratio_list = [...];\n');
fprintf(fid, '```\n');
end

function arr = strip_large_fields_ms(arr, N)
if nargin < 2
    N = 5;
end
if isempty(arr)
    return;
end
for i = 1:numel(arr)
    arr(i).plates = cell(1, N);
    arr(i).Tfields = cell(1, N);
end
end

function assert_rank_scores_finite_ms(arr, label)
if nargin < 2 || isempty(label)
    label = 'eval_set';
end
if isempty(arr)
    return;
end
for i = 1:numel(arr)
    if ~arr(i).success
        continue;
    end
    if isfinite(arr(i).DeltaTN_actual) && isfinite(arr(i).DeltaTN_mean) && ~isfinite(arr(i).rank_score)
        error('Invalid rank_score detected in %s (candidate_id=%d).', label, arr(i).candidate_id);
    end
end
end

function row = empty_stage_metrics_row_ms()
row = struct('rank', NaN, 'candidate_id', NaN, 'stage', NaN, ...
    'n', NaN, 'Lx_mm', NaN, 'Ly_mm', NaN, 'L_max_mm', NaN, 'size_ok', false, ...
    'coverage', NaN, 'coverage_min', NaN, 'coverage_ok', false, ...
    'gap_req_mm', NaN, 'gap_lx_to_next_mm', NaN, 'gap_ly_to_next_mm', NaN, 'gap_ok', false, ...
    'overlap_pair_count', NaN, 'min_edge_gap_mm', NaN, 'min_edge_gap_ok', false, ...
    'stage_spread_K', NaN, 'spread_ratio_to_prev', NaN, 'spread_monotonic_ok', false, ...
    'DeltaTN_actual_K', NaN, 'DeltaT_target_K', NaN, 'DeltaT_abs_error_K', NaN, ...
    'TN_min_K', NaN, 'TN_mean_K', NaN, 'TN_maxmin_K', NaN);
end

function rows = build_stage_metrics_summary_rows_ms(top_opt_full, spec)
rows = repmat(empty_stage_metrics_row_ms(), 0, 1);
if isempty(top_opt_full)
    return;
end
N = spec.stage_count;
tol = 1e-9;
for i = 1:numel(top_opt_full)
    rec = top_opt_full(i);
    tmpl = empty_stage_metrics_row_ms();
    rank_rows = repmat(tmpl, N, 1);
    for k = 1:N
        rank_rows(k).rank = i;
        rank_rows(k).candidate_id = rec.candidate_id;
        rank_rows(k).stage = k;
        rank_rows(k).n = safe_stage_value_ms(rec, 'n', k, NaN);
        rank_rows(k).Lx_mm = safe_stage_value_ms(rec, 'Lx', k, NaN) * 1e3;
        rank_rows(k).Ly_mm = safe_stage_value_ms(rec, 'Ly', k, NaN) * 1e3;
        rank_rows(k).L_max_mm = safe_stage_value_ms(spec.geometry, 'L_max_mm', k, NaN);
        rank_rows(k).size_ok = isfinite(rank_rows(k).Lx_mm) && isfinite(rank_rows(k).Ly_mm) && ...
            isfinite(rank_rows(k).L_max_mm) && rank_rows(k).Lx_mm <= rank_rows(k).L_max_mm + tol && ...
            rank_rows(k).Ly_mm <= rank_rows(k).L_max_mm + tol;
        rank_rows(k).coverage = safe_stage_value_ms(rec, 'cov', k, NaN);
        rank_rows(k).coverage_min = safe_stage_value_ms(spec.geometry, 'coverage_min', k, NaN);
        rank_rows(k).coverage_ok = isfinite(rank_rows(k).coverage) && isfinite(rank_rows(k).coverage_min) && ...
            rank_rows(k).coverage >= rank_rows(k).coverage_min - tol;
        spacing_diag = calc_stage_spacing_metrics_ms(rec.stage_rects{k}, spec.fp_w, spec.fp_h, spec.geometry.min_edge_gap);
        rank_rows(k).overlap_pair_count = spacing_diag.overlap_pair_count;
        if isfinite(spacing_diag.min_edge_gap)
            rank_rows(k).min_edge_gap_mm = spacing_diag.min_edge_gap * 1e3;
        else
            rank_rows(k).min_edge_gap_mm = NaN;
        end
        rank_rows(k).min_edge_gap_ok = spacing_diag.violating_pair_count == 0;
        rank_rows(k).stage_spread_K = safe_stage_value_ms(rec, 'stage_spread', k, NaN);
        if k > 1
            prev_spread = rank_rows(k-1).stage_spread_K;
            if isfinite(prev_spread) && abs(prev_spread) > eps && isfinite(rank_rows(k).stage_spread_K)
                rank_rows(k).spread_ratio_to_prev = rank_rows(k).stage_spread_K / prev_spread;
            end
            rank_rows(k).spread_monotonic_ok = isfinite(prev_spread) && isfinite(rank_rows(k).stage_spread_K) && ...
                rank_rows(k).stage_spread_K <= prev_spread + tol;
        else
            rank_rows(k).spread_monotonic_ok = true;
        end
        if k < N
            rank_rows(k).gap_req_mm = safe_stage_value_ms(spec.geometry, 'pyramid_gap_min_mm', k, NaN);
            rank_rows(k).gap_lx_to_next_mm = rank_rows(k).Lx_mm - safe_stage_value_ms(rec, 'Lx', k+1, NaN) * 1e3;
            rank_rows(k).gap_ly_to_next_mm = rank_rows(k).Ly_mm - safe_stage_value_ms(rec, 'Ly', k+1, NaN) * 1e3;
            rank_rows(k).gap_ok = isfinite(rank_rows(k).gap_req_mm) && isfinite(rank_rows(k).gap_lx_to_next_mm) && ...
                isfinite(rank_rows(k).gap_ly_to_next_mm) && rank_rows(k).gap_lx_to_next_mm >= rank_rows(k).gap_req_mm - tol && ...
                rank_rows(k).gap_ly_to_next_mm >= rank_rows(k).gap_req_mm - tol;
        else
            rank_rows(k).gap_ok = true;
        end
        rank_rows(k).DeltaTN_actual_K = rec.DeltaTN_actual;
        rank_rows(k).DeltaT_target_K = spec.targets.DeltaT_target;
        rank_rows(k).DeltaT_abs_error_K = abs(rec.DeltaTN_actual - spec.targets.DeltaT_target);
        rank_rows(k).TN_min_K = rec.TN_min;
        rank_rows(k).TN_mean_K = rec.TN_mean;
        rank_rows(k).TN_maxmin_K = rec.TN_maxmin;
    end
    rows = [rows; rank_rows]; %#ok<AGROW>
end
end

function write_reasonability_report_ms(out_txt, spec, best_full, stage_rows)
fid = fopen(out_txt, 'w');
if fid < 0
    return;
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Five-stage result reasonability report\n');
fprintf(fid, 'DeltaT_target=%.6f K\n', spec.targets.DeltaT_target);
if isempty(stage_rows) || ~best_full.success
    fprintf(fid, 'Status=insufficient_data\n');
    return;
end
rows_top1 = stage_rows([stage_rows.rank] == 1);
if isempty(rows_top1)
    fprintf(fid, 'Status=missing_rank1_rows\n');
    return;
end
target = spec.targets.DeltaT_target;
abs_err = abs(best_full.DeltaTN_actual - target);
size_ok = all([rows_top1.size_ok]);
coverage_ok = all([rows_top1.coverage_ok]);
if numel(rows_top1) > 1
    gap_ok = all([rows_top1(1:end-1).gap_ok]);
else
    gap_ok = all([rows_top1.gap_ok]);
end
spread_ok = all([rows_top1.spread_monotonic_ok]);
spacing_ok = all([rows_top1.min_edge_gap_ok]);
overlap_total = sum([rows_top1.overlap_pair_count]);
fprintf(fid, 'Top1 candidate_id=%d\n', best_full.candidate_id);
fprintf(fid, 'Top1 DeltaTN_actual=%.6f K\n', best_full.DeltaTN_actual);
fprintf(fid, 'Top1 DeltaTN_mean=%.6f K\n', best_full.DeltaTN_mean);
fprintf(fid, 'Top1 |DeltaTN_actual-DeltaT_target|=%.6f K\n', abs_err);
if best_full.DeltaTN_actual < 0
    fprintf(fid, 'ALERT: DeltaTN_actual is negative (%.6f K), far from target direction.\n', best_full.DeltaTN_actual);
end
fprintf(fid, 'Constraint checks: size_ok=%d, coverage_ok=%d, gap_ok=%d, spread_monotonic_ok=%d, min_edge_gap_ok=%d, overlap_total=%d\n', ...
    logical(size_ok), logical(coverage_ok), logical(gap_ok), logical(spread_ok), logical(spacing_ok), overlap_total);
if abs_err <= 5 && size_ok && coverage_ok && gap_ok && spread_ok && spacing_ok && overlap_total == 0
    verdict = 'reasonable';
else
    verdict = 'needs_review';
end
fprintf(fid, 'Verdict=%s\n', verdict);
fprintf(fid, 'Top1 stage metrics (mm/K):\n');
for k = 1:numel(rows_top1)
    fprintf(fid, ['  Stage %d: n=%d, Lx=%.3f, Ly=%.3f, cov=%.6f, spread=%.6f, ' ...
        'gap_lx=%.3f, gap_ly=%.3f, overlap_pairs=%d, min_edge_gap=%.3f mm\n'], ...
        rows_top1(k).stage, rows_top1(k).n, rows_top1(k).Lx_mm, rows_top1(k).Ly_mm, ...
        rows_top1(k).coverage, rows_top1(k).stage_spread_K, ...
        rows_top1(k).gap_lx_to_next_mm, rows_top1(k).gap_ly_to_next_mm, ...
        rows_top1(k).overlap_pair_count, rows_top1(k).min_edge_gap_mm);
end
end

function v = safe_stage_value_ms(s, fname, idx, default_v)
if nargin < 4
    default_v = NaN;
end
v = default_v;
if ~isstruct(s) || ~isfield(s, fname)
    return;
end
x = s.(fname);
if isempty(x) || numel(x) < idx
    return;
end
v = x(idx);
end

function row = empty_plot_axis_diag_row_ms()
row = struct('rank', NaN, 'stage', NaN, 'plot_kind', '', ...
    'substrate_lx_mm', NaN, 'substrate_ly_mm', NaN, ...
    'layout_target_xmin_mm', NaN, 'layout_target_xmax_mm', NaN, ...
    'layout_target_ymin_mm', NaN, 'layout_target_ymax_mm', NaN, ...
    'layout_actual_xmin_mm', NaN, 'layout_actual_xmax_mm', NaN, ...
    'layout_actual_ymin_mm', NaN, 'layout_actual_ymax_mm', NaN, ...
    'temp_target_xmin_mm', NaN, 'temp_target_xmax_mm', NaN, ...
    'temp_target_ymin_mm', NaN, 'temp_target_ymax_mm', NaN, ...
    'temp_actual_xmin_mm', NaN, 'temp_actual_xmax_mm', NaN, ...
    'temp_actual_ymin_mm', NaN, 'temp_actual_ymax_mm', NaN);
end

function rows = build_plot_axis_diag_rows_ms(rank_idx, layout_info, temp_info, plot_kind)
if nargin < 4
    plot_kind = '';
end
N = size(layout_info.target_axes_mm, 1);
tmpl = empty_plot_axis_diag_row_ms();
rows = repmat(tmpl, N, 1);
for k = 1:N
    rows(k).rank = rank_idx;
    rows(k).stage = k;
    rows(k).plot_kind = plot_kind;
    if size(layout_info.substrate_mm, 2) >= 2
        rows(k).substrate_lx_mm = layout_info.substrate_mm(k,1);
        rows(k).substrate_ly_mm = layout_info.substrate_mm(k,2);
    end
    if size(layout_info.target_axes_mm, 2) == 4
        rows(k).layout_target_xmin_mm = layout_info.target_axes_mm(k,1);
        rows(k).layout_target_xmax_mm = layout_info.target_axes_mm(k,2);
        rows(k).layout_target_ymin_mm = layout_info.target_axes_mm(k,3);
        rows(k).layout_target_ymax_mm = layout_info.target_axes_mm(k,4);
    end
    if size(layout_info.actual_axes_mm, 2) == 4
        rows(k).layout_actual_xmin_mm = layout_info.actual_axes_mm(k,1);
        rows(k).layout_actual_xmax_mm = layout_info.actual_axes_mm(k,2);
        rows(k).layout_actual_ymin_mm = layout_info.actual_axes_mm(k,3);
        rows(k).layout_actual_ymax_mm = layout_info.actual_axes_mm(k,4);
    end
    if size(temp_info.target_axes_mm, 2) == 4
        rows(k).temp_target_xmin_mm = temp_info.target_axes_mm(k,1);
        rows(k).temp_target_xmax_mm = temp_info.target_axes_mm(k,2);
        rows(k).temp_target_ymin_mm = temp_info.target_axes_mm(k,3);
        rows(k).temp_target_ymax_mm = temp_info.target_axes_mm(k,4);
    end
    if size(temp_info.actual_axes_mm, 2) == 4
        rows(k).temp_actual_xmin_mm = temp_info.actual_axes_mm(k,1);
        rows(k).temp_actual_xmax_mm = temp_info.actual_axes_mm(k,2);
        rows(k).temp_actual_ymin_mm = temp_info.actual_axes_mm(k,3);
        rows(k).temp_actual_ymax_mm = temp_info.actual_axes_mm(k,4);
    end
end
end

function tbl = struct_to_table_rows_ms(rows)
if isempty(rows)
    tbl = table();
    return;
end
if isscalar(rows)
    tbl = struct2table(rows, 'AsArray', true);
else
    tbl = struct2table(rows);
end
end

function summarize_failure_messages_ms(eval_arr, tag)
if nargin < 2 || isempty(tag)
    tag = 'Eval';
end
if isempty(eval_arr) || ~isfield(eval_arr, 'success')
    fprintf('[%s] no evaluation records for failure summary.\n', tag);
    return;
end
failed = eval_arr(~[eval_arr.success]);
if isempty(failed)
    fprintf('[%s] no failed candidate.\n', tag);
    return;
end
msgs = string({failed.message});
msgs(strlength(msgs) == 0) = "empty_message";
[u, ~, ic] = unique(msgs, 'stable');
cnt = accumarray(ic, 1);
[cnt_sorted, ord] = sort(cnt, 'descend');
u_sorted = u(ord);
fprintf('[%s] failure histogram (top %d/%d):\n', tag, min(10, numel(u_sorted)), numel(u_sorted));
for i = 1:min(10, numel(u_sorted))
    fprintf('  - %s : %d\n', char(u_sorted(i)), cnt_sorted(i));
end
end

function write_failed_eval_csv_ms(spec, eval_arr, out_csv)
if nargin < 3 || isempty(out_csv) || isempty(eval_arr)
    return;
end
if ~isfield(eval_arr, 'success')
    return;
end
failed = eval_arr(~[eval_arr.success]);
if isempty(failed)
    return;
end
rows = build_failed_rows_ms(failed, spec.stage_count);
if isempty(rows)
    return;
end
try
    tbl = struct_to_table_rows_ms(rows);
    writetable(tbl, out_csv);
catch ME
    warning('Failed to write %s: %s', out_csv, ME.message);
end
end

function rows = build_failed_rows_ms(failed, N)
rows = struct('candidate_id', {}, 'message', {}, 'layout_method', {}, ...
    'symmetry_mode', {}, 'edge_pattern_mode', {}, ...
    'stage_modes', {}, 'stage_trends', {}, ...
    'n', {}, 'Lx', {}, 'Ly', {}, ...
    'I_opt', {}, 's_dense', {}, 's_sparse', {}, 'expo', {}, ...
    'fp_fill_count', {}, ...
    'anis_ratio', {}, 'rank_score', {}, 'DeltaTN_mean', {}, 'DeltaTN_actual', {});
if isempty(failed)
    return;
end
rows = repmat(rows, numel(failed), 1);
for i = 1:numel(failed)
    rows(i).candidate_id = failed(i).candidate_id;
    rows(i).message = failed(i).message;
    rows(i).layout_method = failed(i).layout_method;
    rows(i).symmetry_mode = failed(i).symmetry_mode;
    rows(i).edge_pattern_mode = failed(i).edge_pattern_mode;
    rows(i).stage_modes = strjoin(normalize_text_cell_ms(failed(i).stage_modes, N), '/');
    rows(i).stage_trends = strjoin(normalize_text_cell_ms(failed(i).stage_trends, N), '/');
    rows(i).n = vec_to_inline_str_ms(failed(i).n);
    rows(i).Lx = vec_to_inline_str_ms(failed(i).Lx);
    rows(i).Ly = vec_to_inline_str_ms(failed(i).Ly);
    rows(i).I_opt = failed(i).I_opt;
    rows(i).s_dense = failed(i).s_dense;
    rows(i).s_sparse = failed(i).s_sparse;
    rows(i).expo = failed(i).expo;
    rows(i).fp_fill_count = failed(i).fp_fill_count;
    rows(i).anis_ratio = failed(i).anis_ratio;
    rows(i).rank_score = failed(i).rank_score;
    rows(i).DeltaTN_mean = failed(i).DeltaTN_mean;
    rows(i).DeltaTN_actual = failed(i).DeltaTN_actual;
end
end

function write_eval_failure_summary_csv_ms(eval_arr, out_csv)
rows = struct('category', {}, 'count', {}, 'fraction_of_failed', {});
if nargin < 2 || isempty(out_csv) || isempty(eval_arr) || ~isfield(eval_arr, 'success')
    return;
end
failed = eval_arr(~[eval_arr.success]);
if isempty(failed)
    rows(1).category = 'none';
    rows(1).count = 0;
    rows(1).fraction_of_failed = 0;
else
    cats = cell(numel(failed), 1);
    for i = 1:numel(failed)
        cats{i} = classify_eval_failure_ms(failed(i).message);
    end
    [u, ~, ic] = unique(cats, 'stable');
    cnt = accumarray(ic, 1);
    rows = repmat(struct('category', '', 'count', 0, 'fraction_of_failed', 0), numel(u), 1);
    for i = 1:numel(u)
        rows(i).category = u{i};
        rows(i).count = cnt(i);
        rows(i).fraction_of_failed = cnt(i) / max(1, numel(failed));
    end
end
try
    writetable(struct_to_table_rows_ms(rows), out_csv);
catch ME
    warning('Failed to write %s: %s', out_csv, ME.message);
end
end

function category = classify_eval_failure_ms(msg)
s = lower(strtrim(char(string(msg))));
if isempty(s)
    category = 'empty_message';
elseif contains(s, 'spacing')
    category = 'spacing';
elseif contains(s, 'hard-constraint') || contains(s, 'coverage') || contains(s, 'size_stage') || ...
        contains(s, 'mono_') || contains(s, 'c4')
    category = 'hard_constraint';
elseif contains(s, '0d')
    category = '0d';
elseif contains(s, 'newton')
    category = 'newton';
elseif contains(s, 'footprint')
    category = 'footprint_mapping';
elseif contains(s, 'baseline')
    category = 'baseline';
else
    category = 'other';
end
end

function write_candidate_batch_composition_csv_ms(eval_arr, out_csv)
if nargin < 2 || isempty(out_csv) || isempty(eval_arr)
    return;
end
groups = {'c4_main','c2_explore','shape_edge_anis','other'};
rows = repmat(struct('search_group', '', 'count', 0, 'success_count', 0, ...
    'failed_count', 0, 'mean_DeltaTN_actual', NaN, 'best_DeltaTN_actual', NaN), numel(groups), 1);
for ig = 1:numel(groups)
    g = groups{ig};
    mask = false(numel(eval_arr), 1);
    for i = 1:numel(eval_arr)
        mask(i) = strcmp(classify_candidate_search_group_ms(eval_arr(i)), g);
    end
    arr = eval_arr(mask);
    rows(ig).search_group = g;
    rows(ig).count = numel(arr);
    if isempty(arr)
        continue;
    end
    if isfield(arr, 'success')
        ok = [arr.success];
    else
        ok = false(size(arr));
    end
    rows(ig).success_count = sum(ok);
    rows(ig).failed_count = numel(arr) - rows(ig).success_count;
    vals = [arr(ok).DeltaTN_actual];
    vals = vals(isfinite(vals));
    if ~isempty(vals)
        rows(ig).mean_DeltaTN_actual = mean(vals);
        rows(ig).best_DeltaTN_actual = max(vals);
    end
end
try
    writetable(struct_to_table_rows_ms(rows), out_csv);
catch ME
    warning('Failed to write %s: %s', out_csv, ME.message);
end
end

function write_candidate_repro_manifest_ms(spec, count_solution, top_eval_post, top_opt_full, top_uni_full, out_csv)
if nargin < 6 || isempty(out_csv)
    return;
end
rows = repmat(empty_repro_manifest_row_ms(), 0, 1);
for i = 1:numel(top_opt_full)
    src = 'post_top';
    seed_i = spec.seed;
    pre_rec = empty_eval_struct_ms(spec.stage_count);
    if numel(top_eval_post) >= i
        pre_rec = top_eval_post(i);
    end
    rows = [rows; build_repro_manifest_row_ms(i, src, seed_i, count_solution, pre_rec, top_opt_full(i), spec)]; %#ok<AGROW>
end
if ~isempty(top_uni_full)
    rows = [rows; build_repro_manifest_row_ms(1, 'baseline', spec.seed, count_solution, ...
        empty_eval_struct_ms(spec.stage_count), top_uni_full(1), spec)]; %#ok<AGROW>
end
if isempty(rows)
    return;
end
try
    writetable(struct_to_table_rows_ms(rows), out_csv);
catch ME
    warning('Failed to write %s: %s', out_csv, ME.message);
end
end

function row = empty_repro_manifest_row_ms()
row = struct('rank', NaN, 'source', '', 'seed', NaN, 'candidate_id', NaN, ...
    'geometry_hash', '', 'layout_method', '', 'search_group', '', 'stage_modes', '', 'stage_methods', '', ...
    'n', '', 'count_rank', NaN, 'count_solution_n', '', 'I_A', NaN, 'success', false, ...
    'DeltaTN_actual', NaN, 'DeltaTN_mean', NaN, 'message', '', ...
    'Lx_mm', '', 'Ly_mm', '','coverage', '', ...
    's_dense_mm', NaN, 's_sparse_mm', NaN, 'expo', NaN, 'anis_ratio', NaN, ...
    'shape_mode', '', 'stage_anis', '', 'ring_radius_ratio', NaN, ...
    'ring_width_ratio', NaN, 'band_width_ratio', NaN, 'corner_bias', NaN, ...
    'jitter_ratio', NaN, 'jitter_seed', NaN);
end

function row = build_repro_manifest_row_ms(rank_idx, source, seed_i, count_solution, pre_rec, rec, spec)
row = empty_repro_manifest_row_ms();
row.rank = rank_idx;
row.source = char(string(source));
row.seed = seed_i;
row.candidate_id = rec.candidate_id;
row.geometry_hash = geometry_hash_from_eval_ms(rec, spec.stage_count);
row.layout_method = rec.layout_method;
if isfield(rec, 'search_group'), row.search_group = char(string(rec.search_group)); end
row.stage_modes = strjoin(normalize_text_cell_ms(rec.stage_modes, spec.stage_count), '/');
row.stage_methods = strjoin(normalize_text_cell_ms(rec.stage_methods, spec.stage_count), '/');
row.n = vec_to_inline_str_ms(rec.n);
if isfield(rec, 'count_rank')
    row.count_rank = rec.count_rank;
end
row.count_solution_n = vec_to_inline_str_ms(count_solution.n);
row.I_A = rec.I_opt;
row.success = logical(rec.success);
row.DeltaTN_actual = rec.DeltaTN_actual;
row.DeltaTN_mean = rec.DeltaTN_mean;
row.message = rec.message;
row.Lx_mm = vec_to_inline_str_ms(rec.Lx * 1e3);
row.Ly_mm = vec_to_inline_str_ms(rec.Ly * 1e3);
row.coverage = vec_to_inline_str_ms(rec.cov);
row.s_dense_mm = rec.s_dense * 1e3;
row.s_sparse_mm = rec.s_sparse * 1e3;
row.expo = rec.expo;
row.anis_ratio = rec.anis_ratio;
if isfield(rec, 'shape_mode'), row.shape_mode = char(string(rec.shape_mode)); end
if isfield(rec, 'stage_anis'), row.stage_anis = vec_to_inline_str_ms(rec.stage_anis); end
if isfield(rec, 'ring_radius_ratio'), row.ring_radius_ratio = rec.ring_radius_ratio; end
if isfield(rec, 'ring_width_ratio'), row.ring_width_ratio = rec.ring_width_ratio; end
if isfield(rec, 'band_width_ratio'), row.band_width_ratio = rec.band_width_ratio; end
if isfield(rec, 'corner_bias'), row.corner_bias = rec.corner_bias; end
if isfield(rec, 'jitter_ratio'), row.jitter_ratio = rec.jitter_ratio; end
if isfield(rec, 'jitter_seed'), row.jitter_seed = rec.jitter_seed; end
if isstruct(pre_rec) && isfield(pre_rec, 'candidate_id') && pre_rec.candidate_id == rec.candidate_id && isfinite(pre_rec.I_opt)
    row.I_A = rec.I_opt;
end
end

function h = geometry_hash_from_eval_ms(rec, N)
parts = cell(1, N + 4);
parts{1} = vec_to_inline_str_ms(rec.n);
parts{2} = vec_to_inline_str_ms(rec.Lx);
parts{3} = vec_to_inline_str_ms(rec.Ly);
parts{4} = rec.layout_method;
for k = 1:N
    if isfield(rec, 'stage_rects') && numel(rec.stage_rects) >= k && ~isempty(rec.stage_rects{k})
        parts{4+k} = mat_to_inline_str_ms(round(rec.stage_rects{k} * 1e9) / 1e9);
    else
        parts{4+k} = '';
    end
end
s = strjoin(parts, '|');
h = stable_text_hash_ms(s);
end

function h = stable_text_hash_ms(s)
txt = uint8(char(string(s)));
v = uint32(2166136261);
prime = uint32(16777619);
for i = 1:numel(txt)
    v = bitxor(v, uint32(txt(i)));
    v = uint32(mod(uint64(v) * uint64(prime), uint64(2)^32));
end
h = lower(dec2hex(v, 8));
end

function c = normalize_text_cell_ms(raw, N)
if nargin < 2
    N = numel(raw);
end
c = repmat({''}, 1, N);
if isempty(raw)
    return;
end
for i = 1:N
    idx = min(i, numel(raw));
    if idx >= 1 && ~isempty(raw{idx})
        c{i} = char(string(raw{idx}));
    else
        c{i} = '';
    end
end
end

function tf = has_valid_candidate_geometry_ms(rec, N)
tf = isstruct(rec) && isfield(rec, 'stage_rects') && numel(rec.stage_rects) >= N;
if ~tf
    return;
end
for k = 1:N
    rk = rec.stage_rects{k};
    if isempty(rk)
        continue;
    end
    if size(rk,2) ~= 4 || ~all(isfinite(rk(:)))
        tf = false;
        return;
    end
end
end

function write_stage_layout_rects_csv_ms(rec, out_csv, N)
if nargin < 3 || ~isfinite(N) || N < 1
    N = 5;
end
if ~has_valid_candidate_geometry_ms(rec, N)
    return;
end
rows = struct('stage', {}, 'particle_idx', {}, ...
    'x1_mm', {}, 'x2_mm', {}, 'y1_mm', {}, 'y2_mm', {}, ...
    'cx_mm', {}, 'cy_mm', {});
for k = 1:N
    rects = rec.stage_rects{k};
    if isempty(rects)
        continue;
    end
    for i = 1:size(rects,1)
        r = struct();
        r.stage = k;
        r.particle_idx = i;
        r.x1_mm = rects(i,1) * 1e3;
        r.x2_mm = rects(i,2) * 1e3;
        r.y1_mm = rects(i,3) * 1e3;
        r.y2_mm = rects(i,4) * 1e3;
        r.cx_mm = 0.5 * (r.x1_mm + r.x2_mm);
        r.cy_mm = 0.5 * (r.y1_mm + r.y2_mm);
        rows = [rows; r]; %#ok<AGROW>
    end
end
if isempty(rows)
    return;
end
writetable(struct_to_table_rows_ms(rows), out_csv);
end

function write_baseline_status_txt_ms(out_txt, baseline_info, top_uni_full, N)
fid = fopen(out_txt, 'w');
if fid < 0
    return;
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'baseline status\n');
if nargin >= 2 && isstruct(baseline_info)
    fprintf(fid, 'success=%d\n', logical(baseline_info.success));
    fprintf(fid, 'selected_source=%s\n', char(string(baseline_info.selected_source)));
    fprintf(fid, 'selected_scale=%.6f\n', baseline_info.selected_scale);
    fprintf(fid, 'target_total_n=%.0f\n', baseline_info.target_total_n);
    if isfield(baseline_info, 'selected_n') && ~isempty(baseline_info.selected_n)
        fprintf(fid, 'selected_n=%s\n', vec_to_inline_str_ms(baseline_info.selected_n));
    end
    if isfield(baseline_info, 'message')
        fprintf(fid, 'message=%s\n', char(string(baseline_info.message)));
    end
end
if nargin >= 3 && ~isempty(top_uni_full)
    u = top_uni_full(1);
    fprintf(fid, 'eval_success=%d\n', logical(u.success));
    fprintf(fid, 'eval_layout_method=%s\n', u.layout_method);
    fprintf(fid, 'eval_message=%s\n', u.message);
    if nargin >= 4 && has_valid_candidate_geometry_ms(u, N)
        fprintf(fid, 'eval_n=%s\n', vec_to_inline_str_ms(u.n));
    end
end
if nargin >= 2 && isfield(baseline_info, 'attempts') && ~isempty(baseline_info.attempts)
    fprintf(fid, 'attempt_count=%d\n', numel(baseline_info.attempts));
    for i = 1:numel(baseline_info.attempts)
        a = baseline_info.attempts(i);
        fprintf(fid, 'attempt%d: source=%s, scale=%.6f, candidate_id=%d, n=%s, build=\"%s\", eval_ran=%d, eval_success=%d, eval_msg=\"%s\"\n', ...
            i, a.source, a.scale, a.candidate_id, a.n, a.build_message, ...
            logical(a.eval_ran), logical(a.eval_success), a.eval_message);
    end
end
end

%% ====================== FEM/mesh helpers ======================

function plate = build_plate_struct_ms(P)
plate.Lx = P.Lx;
plate.Ly = P.Ly;
if isfield(P, 't') && isfinite(P.t) && P.t > 0
    plate.t = P.t;
else
    plate.t = 1e-3;
end
plate.nx = P.nx;
plate.ny = P.ny;
[x, y, tri] = build_structured_tri_mesh_ms(P.Lx, P.Ly, P.nx, P.ny);
if isfield(P, 'k_inplane') && isfinite(P.k_inplane) && P.k_inplane > 0
    k_inplane = P.k_inplane;
else
    k_inplane = 170;
end
elem = precompute_tri_geom_ms(x, y, tri);
m = lumped_area_weights_ms(tri, elem.A, numel(x));
plate.mesh = struct('x', x, 'y', y, 'tri', tri, 'Nn', numel(x), 'Ne', size(tri,1));
plate.elem = elem;
plate.m = m;
plate.A_total = sum(m);
plate.inv_sum_m = 1 / max(plate.A_total, eps);
plate.K = assemble_stiffness_const_ms(k_inplane, plate);
reg_param = max(1e-12, 1e-10 * max(diag(plate.K)));
plate.Aaug = [plate.K, plate.m; plate.m', -reg_param];
end

function [x, y, tri] = build_structured_tri_mesh_ms(Lx, Ly, nx, ny)
xv = linspace(-Lx/2, Lx/2, nx);
yv = linspace(-Ly/2, Ly/2, ny);
[X, Y] = meshgrid(xv, yv);
x = X(:);
y = Y(:);
node = @(i,j) (j-1)*nx + i;
tri = zeros(2*(nx-1)*(ny-1),3);
t = 0;
for j = 1:ny-1
    for i = 1:nx-1
        n1 = node(i,j);
        n2 = node(i+1,j);
        n3 = node(i,j+1);
        n4 = node(i+1,j+1);
        t = t + 1; tri(t,:) = [n1 n2 n4];
        t = t + 1; tri(t,:) = [n1 n4 n3];
    end
end
end

function elem = precompute_tri_geom_ms(x, y, tri)
Ne = size(tri,1);
A = zeros(Ne,1);
b = zeros(Ne,3);
c = zeros(Ne,3);
cx = zeros(Ne,1);
cy = zeros(Ne,1);
for e = 1:Ne
    id = tri(e,:);
    x1 = x(id(1)); y1 = y(id(1));
    x2 = x(id(2)); y2 = y(id(2));
    x3 = x(id(3)); y3 = y(id(3));
    A(e) = 0.5 * abs(det([x2-x1, x3-x1; y2-y1, y3-y1]));
    b(e,:) = [y2-y3, y3-y1, y1-y2];
    c(e,:) = [x3-x2, x1-x3, x2-x1];
    cx(e) = (x1 + x2 + x3) / 3;
    cy(e) = (y1 + y2 + y3) / 3;
end
elem = struct('A', A, 'b', b, 'c', c, 'cx', cx, 'cy', cy);
end

function m = lumped_area_weights_ms(tri, Aelem, Nn)
m = zeros(Nn,1);
Ne = size(tri,1);
for e = 1:Ne
    m(tri(e,:)) = m(tri(e,:)) + Aelem(e)/3;
end
end

function K = assemble_stiffness_const_ms(k_const, P)
Ne = P.mesh.Ne;
tri = P.mesh.tri;
Iind = zeros(9*Ne,1);
Jind = zeros(9*Ne,1);
Vval = zeros(9*Ne,1);
idx = 0;
for e = 1:Ne
    nodes = tri(e,:);
    Ae = P.elem.A(e);
    b = P.elem.b(e,:);
    c = P.elem.c(e,:);
    Ke = (k_const * P.t) * ((b.'*b + c.'*c) / (4*Ae));
    for a = 1:3
        for bb = 1:3
            idx = idx + 1;
            Iind(idx) = nodes(a);
            Jind(idx) = nodes(bb);
            Vval(idx) = Ke(a,bb);
        end
    end
end
K = sparse(Iind(1:idx), Jind(1:idx), Vval(1:idx), P.mesh.Nn, P.mesh.Nn);
end

function [fpId_cool, fpId_hot] = map_elements_to_footprints_ms(cx, cy, rect_cool, rect_hot)
Ne = numel(cx);
fpId_cool = zeros(Ne,1);
fpId_hot = zeros(Ne,1);
for p = 1:size(rect_cool,1)
    r = rect_cool(p,:);
    inside = (cx>=r(1) & cx<=r(2) & cy>=r(3) & cy<=r(4));
    fpId_cool(inside) = p;
end
for p = 1:size(rect_hot,1)
    r = rect_hot(p,:);
    inside = (cx>=r(1) & cx<=r(2) & cy>=r(3) & cy<=r(4));
    fpId_hot(inside) = p;
end
end

function fpElems = ids_to_cell_ms(fpId, nfp)
fpElems = cell(nfp,1);
Ne = numel(fpId);
for e = 1:Ne
    p = fpId(e);
    if p > 0
        fpElems{p}(end+1) = e; %#ok<AGROW>
    end
end
end

function [fpElems, fill_count] = fill_empty_fp_elems_nearest_ms(fpElems, rects, cx, cy)
fill_count = 0;
if isempty(fpElems) || isempty(rects)
    return;
end
for p = 1:min(numel(fpElems), size(rects,1))
    if ~isempty(fpElems{p})
        continue;
    end
    rcx = 0.5 * (rects(p,1) + rects(p,2));
    rcy = 0.5 * (rects(p,3) + rects(p,4));
    [~, idx] = min((cx-rcx).^2 + (cy-rcy).^2);
    if ~isempty(idx) && isfinite(idx)
        fpElems{p} = idx;
        fill_count = fill_count + 1;
    end
end
end

function Tfp = footprint_avg_ms(Tnode, fpElems, elem, tri)
nf = numel(fpElems);
Tfp = zeros(nf,1);
for p = 1:nf
    elist = fpElems{p};
    if isempty(elist)
        Tfp(p) = NaN;
        continue;
    end
    Ae = elem.A(elist);
    Te = zeros(numel(elist),1);
    for ii = 1:numel(elist)
        e = elist(ii);
        Te(ii) = mean(Tnode(tri(e,:)));
    end
    Tfp(p) = sum(Ae .* Te) / sum(Ae);
end
if any(isnan(Tfp))
    valid = Tfp(~isnan(Tfp));
    if isempty(valid)
        fill_v = 0;
    else
        fill_v = mean(valid);
    end
    Tfp(isnan(Tfp)) = fill_v;
end
end

%% ====================== TE model ======================

function Qc = te_Qc_onecouple_ms(Tc, Th, Icur, G, stage_idx)
if nargin < 5 || isempty(stage_idx)
    stage_idx = 1;
end
[L_leg, A_leg] = get_stage_leg_geom_ms(G, stage_idx);
[alpha_avg, k_avg, rho_avg] = avg_props_ms(Tc, Th, G);
Kcond = k_avg * A_leg / L_leg;
Relec = G.Rc + rho_avg * L_leg / A_leg;
Qc = alpha_avg*Icur*Tc - 0.5*(Icur^2)*Relec - Kcond*(Th-Tc);
end
function Qh = te_Qh_onecouple_ms(Th, Tc, Icur, G, stage_idx)
if nargin < 5 || isempty(stage_idx)
    stage_idx = 1;
end
[L_leg, A_leg] = get_stage_leg_geom_ms(G, stage_idx);
[alpha_avg, k_avg, rho_avg] = avg_props_ms(Tc, Th, G);
Kcond = k_avg * A_leg / L_leg;
Relec = G.Rc + rho_avg * L_leg / A_leg;
Qh = alpha_avg*Icur*Th + 0.5*(Icur^2)*Relec - Kcond*(Th-Tc);
end

function [L_leg, A_leg] = get_stage_leg_geom_ms(G, stage_idx)
if nargin < 2 || isempty(stage_idx) || ~isfinite(stage_idx)
    stage_idx = 1;
end
A_leg = G.A_leg;
if isfield(G, 'L_leg_stage') && ~isempty(G.L_leg_stage)
    idx = min(max(1, round(stage_idx)), numel(G.L_leg_stage));
    L_leg = G.L_leg_stage(idx);
elseif isfield(G, 'fp_t_stage') && ~isempty(G.fp_t_stage)
    idx = min(max(1, round(stage_idx)), numel(G.fp_t_stage));
    L_leg = G.fp_t_stage(idx);
else
    L_leg = G.L_leg;
end
L_leg = max(L_leg, 1e-12);
A_leg = max(A_leg, 1e-18);
end
function [alpha_avg, k_avg, rho_avg] = avg_props_ms(T1, T2, G)
if abs(T2-T1) < 1e-12
    alpha_avg = G.alpha_fun(T1);
    k_avg = G.k_leg_fun(T1);
    rho_avg = G.rho_fun(T1);
    return;
end
Ta = min(T1, T2);
Tb = max(T1, T2);
dT = Tb - Ta;
alpha_avg = integral_poly_ms(G.material.alpha_coeffs, Ta, Tb) / dT;
k_avg = integral_poly_ms(G.material.k_coeffs, Ta, Tb) / dT;
rho_avg = integral_poly_ms(G.material.rho_coeffs, Ta, Tb) / dT;
end

function F = integral_poly_ms(coeffs, Ta, Tb)
a4 = coeffs(1); a3 = coeffs(2); a2 = coeffs(3); a1 = coeffs(4); a0 = coeffs(5);
F_Ta = a4*Ta^5/5 + a3*Ta^4/4 + a2*Ta^3/3 + a1*Ta^2/2 + a0*Ta;
F_Tb = a4*Tb^5/5 + a3*Tb^4/4 + a2*Tb^3/3 + a1*Tb^2/2 + a0*Tb;
F = F_Tb - F_Ta;
end

function zcfg = make_zpath_runtime_ms(spec, fp_A, N)
if nargin < 3 || ~isfinite(N) || N < 2
    N = 2;
end
if nargin < 2 || ~isfinite(fp_A) || fp_A <= 0
    fp_A = 1e-6;
end
nif = max(1, N-1);
zcfg = struct( ...
    'enable', false, ...
    'k_interfaces', 170 * ones(1, nif), ...
    't_interface_effs', 1e-3 * ones(1, nif), ...
    'Rc_interfaces', zeros(1, nif), ...
    'k_sink', 170, ...
    't_sink_eff', 1e-3, ...
    'Rc_sink', 0, ...
    'step_fp_iters', 2, ...
    'step_fp_relax', 0.6, ...
    'step_q_prev_weight', 0.85, ...
    'Rz_interfaces', zeros(1, nif), ...
    'Rz_sink', 0);

if nargin < 1 || ~isstruct(spec)
    return;
end
if isfield(spec, 'plate_t') && isfinite(spec.plate_t) && spec.plate_t >= 0
    zcfg.t_interface_effs = spec.plate_t * ones(1, nif);
    zcfg.t_sink_eff = spec.plate_t;
end
if ~isfield(spec, 'z_path') || ~isstruct(spec.z_path)
    return;
end

zp = spec.z_path;
if isfield(zp, 'enable')
    zcfg.enable = logical(zp.enable);
end
if isfield(zp, 'k_interfaces') && ~isempty(zp.k_interfaces)
    zcfg.k_interfaces = max(fit_len_vec_ms(zp.k_interfaces, nif, zcfg.k_interfaces), 1e-9);
end
if isfield(zp, 't_interface_effs') && ~isempty(zp.t_interface_effs)
    zcfg.t_interface_effs = max(fit_len_vec_ms(zp.t_interface_effs, nif, zcfg.t_interface_effs), 0);
end
if isfield(zp, 'Rc_interfaces') && ~isempty(zp.Rc_interfaces)
    zcfg.Rc_interfaces = max(fit_len_vec_ms(zp.Rc_interfaces, nif, zcfg.Rc_interfaces), 0);
end
if isfield(zp, 'k_sink'), zcfg.k_sink = max(zp.k_sink, 1e-9); end
if isfield(zp, 't_sink_eff'), zcfg.t_sink_eff = max(zp.t_sink_eff, 0); end
if isfield(zp, 'Rc_sink'), zcfg.Rc_sink = max(zp.Rc_sink, 0); end
if isfield(zp, 'step_fp_iters')
    zcfg.step_fp_iters = max(1, round(zp.step_fp_iters));
end
if isfield(zp, 'step_fp_relax')
    zcfg.step_fp_relax = min(max(zp.step_fp_relax, 0), 1);
end
if isfield(zp, 'step_q_prev_weight')
    zcfg.step_q_prev_weight = min(max(zp.step_q_prev_weight, 0), 1);
end

if zcfg.enable
    zcfg.Rz_interfaces = zcfg.t_interface_effs ./ (zcfg.k_interfaces * fp_A) + zcfg.Rc_interfaces;
    zcfg.Rz_sink = zcfg.t_sink_eff / (zcfg.k_sink * fp_A) + zcfg.Rc_sink;
else
    zcfg.Rz_interfaces = zeros(1, nif);
    zcfg.Rz_sink = 0;
end
end

function zcfg = get_zpath_runtime_from_G_ms(G, N)
if nargin < 2 || ~isfinite(N) || N < 2
    N = 2;
end
nif = max(1, N-1);
if isstruct(G) && isfield(G, 'z_path_runtime') && isstruct(G.z_path_runtime)
    zcfg = G.z_path_runtime;
    if ~isfield(zcfg, 'Rz_interfaces') || isempty(zcfg.Rz_interfaces)
        zcfg.Rz_interfaces = zeros(1, nif);
    end
    zcfg.Rz_interfaces = fit_len_vec_ms(zcfg.Rz_interfaces, nif, zeros(1, nif));
    if ~isfield(zcfg, 'Rz_sink') || ~isfinite(zcfg.Rz_sink)
        zcfg.Rz_sink = 0;
    end
    if ~isfield(zcfg, 'enable')
        zcfg.enable = false;
    end
    if ~isfield(zcfg, 'step_q_prev_weight')
        zcfg.step_q_prev_weight = 0.85;
    end
    return;
end
zcfg = struct('enable', false, 'Rz_interfaces', zeros(1, nif), 'Rz_sink', 0, ...
    'step_fp_iters', 2, 'step_fp_relax', 0.6, 'step_q_prev_weight', 0.85);
end

%% ====================== utilities ======================

function blk = choose_parfor_block_size_ms(N, nw, blk_min, blk_max, target_blocks_per_worker)
if nargin < 5
    target_blocks_per_worker = 2;
end
if nargin < 4
    blk_max = 96;
end
if nargin < 3
    blk_min = 16;
end
if nargin < 2 || nw < 1
    nw = 1;
end
target_blocks_per_worker = max(1, round(target_blocks_per_worker));
if N <= 2 * nw
    blk = N;
    return;
end
target = ceil(N / (target_blocks_per_worker * nw));
blk = min(blk_max, max(blk_min, target));
end

function v = make_even_ge_ms(x)
v = ceil(x);
if mod(v,2) ~= 0
    v = v + 1;
end
end

function v = make_even_le_ms(x)
v = floor(x);
if mod(v,2) ~= 0
    v = v - 1;
end
end

%% ====================== Plotting functions ======================

function info = empty_plot_axis_info_ms(N)
info = struct('target_axes_mm', NaN(N, 4), 'actual_axes_mm', NaN(N, 4), 'substrate_mm', NaN(N, 2));
end

function info = save_layout_plot_ms(eval_rec, out_png, spec)
N = spec.stage_count;
plot_cfg = resolve_plot_cfg_ms(spec);
global_axis_mm = make_global_axis_limits_mm_ms(spec);
info = empty_plot_axis_info_ms(N);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 800]);

for k = 1:N
    ax = subplot(N, 1, k);
    [axis_target_k, axis_actual_k, substrate_k] = render_layout_stage_on_axis_ms(ax, eval_rec, k, plot_cfg, global_axis_mm);
    info.target_axes_mm(k,:) = axis_target_k;
    info.actual_axes_mm(k,:) = axis_actual_k;
    info.substrate_mm(k,:) = substrate_k;
end

sgtitle('Particle Layout');
exportgraphics(fig, out_png, 'Resolution', 150);
close(fig);
end

function info = save_layout_stage_plots_ms(eval_rec, out_dir, spec)
N = spec.stage_count;
plot_cfg = resolve_plot_cfg_ms(spec);
global_axis_mm = make_global_axis_limits_mm_ms(spec);
info = empty_plot_axis_info_ms(N);
for k = 1:N
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [100 100 plot_cfg.separate_fig_size_px(1) plot_cfg.separate_fig_size_px(2)]);
    ax = axes(fig);
    [axis_target_k, axis_actual_k, substrate_k] = render_layout_stage_on_axis_ms(ax, eval_rec, k, plot_cfg, global_axis_mm);
    info.target_axes_mm(k,:) = axis_target_k;
    info.actual_axes_mm(k,:) = axis_actual_k;
    info.substrate_mm(k,:) = substrate_k;
    out_png = fullfile(out_dir, sprintf('layout_stage%d.png', k));
    exportgraphics(fig, out_png, 'Resolution', 180);
    close(fig);
end
end

function info = save_temperature_plot_ms(eval_rec, out_png, caxis_by_stage, spec)
N = spec.stage_count;
if nargin < 4
    spec = struct('geometry', struct('L_max', ones(1,N)));
end
plot_cfg = resolve_plot_cfg_ms(spec);
global_axis_mm = make_global_axis_limits_mm_ms(spec);
info = empty_plot_axis_info_ms(N);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1400 1000]);
colormap(fig, 'parula');

for k = 1:N
    ax = subplot(N, 1, k);
    [axis_target_k, axis_actual_k, substrate_k] = render_temperature_stage_on_axis_ms( ...
        ax, eval_rec, k, caxis_by_stage, plot_cfg, global_axis_mm);
    info.target_axes_mm(k,:) = axis_target_k;
    info.actual_axes_mm(k,:) = axis_actual_k;
    info.substrate_mm(k,:) = substrate_k;
end

sgtitle('Temperature Distribution');
exportgraphics(fig, out_png, 'Resolution', 220);
close(fig);
end

function info = save_temperature_stage_plots_ms(eval_rec, out_dir, caxis_by_stage, spec)
N = spec.stage_count;
plot_cfg = resolve_plot_cfg_ms(spec);
global_axis_mm = make_global_axis_limits_mm_ms(spec);
info = empty_plot_axis_info_ms(N);
for k = 1:N
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [100 100 plot_cfg.separate_fig_size_px(1) plot_cfg.separate_fig_size_px(2)]);
    colormap(fig, 'parula');
    ax = axes(fig);
    [axis_target_k, axis_actual_k, substrate_k] = render_temperature_stage_on_axis_ms( ...
        ax, eval_rec, k, caxis_by_stage, plot_cfg, global_axis_mm);
    info.target_axes_mm(k,:) = axis_target_k;
    info.actual_axes_mm(k,:) = axis_actual_k;
    info.substrate_mm(k,:) = substrate_k;
    out_png = fullfile(out_dir, sprintf('temperature_stage%d.png', k));
    exportgraphics(fig, out_png, 'Resolution', 220);
    close(fig);
end
end

function [axis_target_mm, axis_actual_mm, substrate_mm] = render_layout_stage_on_axis_ms(ax, eval_rec, k, plot_cfg, global_axis_mm)
hold(ax, 'on');
if isfield(eval_rec, 'stage_rects') && numel(eval_rec.stage_rects) >= k && ~isempty(eval_rec.stage_rects{k})
    rects = eval_rec.stage_rects{k};
    for r = 1:size(rects, 1)
        rectangle(ax, 'Position', [rects(r,1)*1e3, rects(r,3)*1e3, (rects(r,2)-rects(r,1))*1e3, (rects(r,4)-rects(r,3))*1e3], ...
            'EdgeColor', 'b', 'LineWidth', 1);
    end
end
[axis_target_mm, axis_actual_mm, substrate_mm] = apply_axis_and_substrate_ms(ax, eval_rec, k, plot_cfg, global_axis_mm, []);
stage_n = safe_stage_value_ms(eval_rec, 'n', k, NaN);
if isfinite(stage_n)
    n_tag = sprintf('%d', round(stage_n));
else
    n_tag = 'NA';
end
title(ax, sprintf('Stage %d Layout (n=%s)', k, n_tag));
xlabel(ax, 'x / mm');
ylabel(ax, 'y / mm');
grid(ax, 'on');
hold(ax, 'off');
end

function [axis_target_mm, axis_actual_mm, substrate_mm] = render_temperature_stage_on_axis_ms( ...
    ax, eval_rec, k, caxis_by_stage, plot_cfg, global_axis_mm)
edge_color = 'none';
if plot_cfg.show_mesh_edges
    edge_color = 'k';
end
hold(ax, 'on');
has_plate = isfield(eval_rec, 'plates') && numel(eval_rec.plates) >= k && ~isempty(eval_rec.plates{k});
has_field = isfield(eval_rec, 'Tfields') && numel(eval_rec.Tfields) >= k && ~isempty(eval_rec.Tfields{k});
plate_k = [];
if has_plate
    plate_k = eval_rec.plates{k};
end
if has_plate && has_field
    Tk = eval_rec.Tfields{k};
    Pk = eval_rec.plates{k};
    [ok_render, render_msg] = draw_temperature_patch_ms(ax, Pk, Tk, edge_color, plot_cfg.view_mode);
    if strcmp(plot_cfg.view_mode, '3d')
        view(ax, 3);
    else
        view(ax, 2);
    end
    if ok_render
        colorbar(ax);
        if size(caxis_by_stage,1) >= k && size(caxis_by_stage,2) == 2 && all(isfinite(caxis_by_stage(k,:)))
            caxis(ax, caxis_by_stage(k,:));
        end
        title(ax, sprintf('Stage %d: Tmin=%.2f K, Tmax=%.2f K, \\DeltaTspread=%.2f K', ...
            k, min(Tk), max(Tk), safe_stage_value_ms(eval_rec, 'stage_spread', k, NaN)));
    else
        title(ax, sprintf('Stage %d: temperature render failed (%s)', k, render_msg));
    end
else
    title(ax, sprintf('Stage %d: No temperature data', k));
end
[axis_target_mm, axis_actual_mm, substrate_mm] = apply_axis_and_substrate_ms(ax, eval_rec, k, plot_cfg, global_axis_mm, plate_k);
xlabel(ax, 'x / mm');
ylabel(ax, 'y / mm');
grid(ax, 'on');
hold(ax, 'off');
end

function [ok, msg] = draw_temperature_patch_ms(ax, plate_k, Tk, edge_color, view_mode)
ok = false;
msg = 'invalid_input';
if nargin < 5 || isempty(view_mode)
    view_mode = '2d';
end
if isempty(plate_k) || ~isfield(plate_k, 'mesh') || isempty(plate_k.mesh)
    msg = 'missing_mesh';
    return;
end
if ~isfield(plate_k.mesh, 'tri') || ~isfield(plate_k.mesh, 'x') || ~isfield(plate_k.mesh, 'y')
    msg = 'incomplete_mesh';
    return;
end
tri = plate_k.mesh.tri;
x = plate_k.mesh.x(:) * 1e3;
y = plate_k.mesh.y(:) * 1e3;
if size(tri,2) ~= 3 || isempty(tri)
    msg = 'invalid_triangles';
    return;
end
if numel(Tk) ~= numel(x)
    msg = 'node_value_size_mismatch';
    return;
end
if ~all(isfinite(x)) || ~all(isfinite(y)) || ~all(isfinite(Tk(:)))
    msg = 'nonfinite_values';
    return;
end
if strcmp(view_mode, '3d')
    verts = [x, y, Tk(:)];
else
    verts = [x, y];
end
patch('Faces', tri, 'Vertices', verts, ...
    'FaceVertexCData', Tk(:), ...
    'FaceColor', 'interp', ...
    'EdgeColor', edge_color, ...
    'Parent', ax);
ok = true;
msg = 'ok';
end

function [axis_target_mm, axis_actual_mm, substrate_mm] = apply_axis_and_substrate_ms( ...
    ax, eval_rec, stage_idx, plot_cfg, global_axis_mm, plate_k)
axis_target_mm = resolve_stage_axis_limits_mm_ms(eval_rec, stage_idx, plot_cfg, global_axis_mm, plate_k);
xlim(ax, axis_target_mm(1:2));
ylim(ax, axis_target_mm(3:4));
axis(ax, 'equal');
[Lx_m, Ly_m] = resolve_stage_substrate_size_m_ms(eval_rec, stage_idx, plate_k);
substrate_mm = [Lx_m, Ly_m] * 1e3;
draw_substrate_overlay_ms(ax, substrate_mm, plot_cfg);
xl = xlim(ax);
yl = ylim(ax);
axis_actual_mm = [xl(1), xl(2), yl(1), yl(2)];
end

function draw_substrate_overlay_ms(ax, substrate_mm, plot_cfg)
if numel(substrate_mm) < 2 || ~all(isfinite(substrate_mm)) || any(substrate_mm <= 0)
    return;
end
rectangle(ax, 'Position', [-substrate_mm(1)/2, -substrate_mm(2)/2, substrate_mm(1), substrate_mm(2)], ...
    'EdgeColor', [0.1 0.1 0.1], 'LineStyle', '--', 'LineWidth', 1.0);
if plot_cfg.annotate_substrate_dims
    txt = sprintf('Substrate: Lx=%.2f mm, Ly=%.2f mm', substrate_mm(1), substrate_mm(2));
    text(ax, 0.01, 0.98, txt, 'Units', 'normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'Color', [0.1 0.1 0.1], 'BackgroundColor', 'w', 'Margin', 1, 'FontSize', 9);
end
end

function cfg = resolve_plot_cfg_ms(spec)
cfg = struct('use_global_axis', false, 'show_mesh_edges', false, ...
    'view_mode', '2d', 'axis_margin_ratio', 0.05, ...
    'save_overview', true, 'save_stage_separate', true, ...
    'annotate_substrate_dims', true, 'separate_fig_size_px', [900, 900]);
if nargin < 1 || ~isstruct(spec) || ~isfield(spec, 'output') || ~isstruct(spec.output) || ...
        ~isfield(spec.output, 'plot') || ~isstruct(spec.output.plot)
    return;
end
pc = spec.output.plot;
if isfield(pc, 'use_global_axis')
    cfg.use_global_axis = any(logical(pc.use_global_axis(:)));
end
if isfield(pc, 'show_mesh_edges')
    cfg.show_mesh_edges = any(logical(pc.show_mesh_edges(:)));
end
if isfield(pc, 'view_mode')
    cfg.view_mode = normalize_plot_view_mode_ms(pc.view_mode);
end
if isfield(pc, 'axis_margin_ratio') && isfinite(pc.axis_margin_ratio)
    cfg.axis_margin_ratio = min(max(pc.axis_margin_ratio, 0), 0.5);
end
if isfield(pc, 'save_overview')
    cfg.save_overview = any(logical(pc.save_overview(:)));
end
if isfield(pc, 'save_stage_separate')
    cfg.save_stage_separate = any(logical(pc.save_stage_separate(:)));
end
if isfield(pc, 'annotate_substrate_dims')
    cfg.annotate_substrate_dims = any(logical(pc.annotate_substrate_dims(:)));
end
if isfield(pc, 'separate_fig_size_px')
    cfg.separate_fig_size_px = normalize_plot_size_px_ms(pc.separate_fig_size_px);
end
end

function mode = normalize_plot_view_mode_ms(v)
mode = '2d';
if nargin < 1 || isempty(v)
    return;
end
if isstring(v)
    v = char(v(1));
elseif ~ischar(v)
    return;
end
s = lower(strtrim(v));
if strcmp(s, '3d')
    mode = '3d';
end
end

function fig_size = normalize_plot_size_px_ms(v)
fig_size = [900, 900];
if nargin < 1 || isempty(v)
    return;
end
raw = double(v(:).');
if numel(raw) == 1
    raw = [raw, raw];
end
if numel(raw) < 2
    return;
end
if all(isfinite(raw(1:2)))
    fig_size = round(raw(1:2));
end
fig_size = max(fig_size, [300, 300]);
end

function axis_mm = make_global_axis_limits_mm_ms(spec)
half_mm = 40;
if nargin >= 1 && isstruct(spec) && isfield(spec, 'geometry') && isstruct(spec.geometry) && ...
        isfield(spec.geometry, 'L_max_mm') && ~isempty(spec.geometry.L_max_mm)
    vmax = max(spec.geometry.L_max_mm(:));
    if isfinite(vmax) && vmax > 0
        half_mm = 0.5 * vmax;
    end
end
axis_mm = [-half_mm, half_mm, -half_mm, half_mm];
end

function axis_mm = resolve_stage_axis_limits_mm_ms(eval_rec, stage_idx, plot_cfg, global_axis_mm, plate_k)
if nargin < 5
    plate_k = [];
end
if plot_cfg.use_global_axis
    axis_mm = global_axis_mm;
    return;
end

[Lx_m, Ly_m] = resolve_stage_substrate_size_m_ms(eval_rec, stage_idx, plate_k);

if ~(isfinite(Lx_m) && Lx_m > 0)
    Lx_m = (global_axis_mm(2) - global_axis_mm(1)) * 1e-3;
end
if ~(isfinite(Ly_m) && Ly_m > 0)
    Ly_m = (global_axis_mm(4) - global_axis_mm(3)) * 1e-3;
end

half_x_mm = 0.5 * Lx_m * 1e3 * (1 + plot_cfg.axis_margin_ratio);
half_y_mm = 0.5 * Ly_m * 1e3 * (1 + plot_cfg.axis_margin_ratio);
half_x_mm = max(half_x_mm, 1e-3);
half_y_mm = max(half_y_mm, 1e-3);
axis_mm = [-half_x_mm, half_x_mm, -half_y_mm, half_y_mm];
end

function [Lx_m, Ly_m] = resolve_stage_substrate_size_m_ms(eval_rec, stage_idx, plate_k)
Lx_m = NaN;
Ly_m = NaN;
if ~isempty(plate_k) && isstruct(plate_k)
    if isfield(plate_k, 'Lx') && isfinite(plate_k.Lx) && plate_k.Lx > 0
        Lx_m = plate_k.Lx;
    end
    if isfield(plate_k, 'Ly') && isfinite(plate_k.Ly) && plate_k.Ly > 0
        Ly_m = plate_k.Ly;
    end
end
if ~(isfinite(Lx_m) && Lx_m > 0) && isfield(eval_rec, 'Lx') && numel(eval_rec.Lx) >= stage_idx
    if isfinite(eval_rec.Lx(stage_idx)) && eval_rec.Lx(stage_idx) > 0
        Lx_m = eval_rec.Lx(stage_idx);
    end
end
if ~(isfinite(Ly_m) && Ly_m > 0) && isfield(eval_rec, 'Ly') && numel(eval_rec.Ly) >= stage_idx
    if isfinite(eval_rec.Ly(stage_idx)) && eval_rec.Ly(stage_idx) > 0
        Ly_m = eval_rec.Ly(stage_idx);
    end
end
if ~(isfinite(Lx_m) && Lx_m > 0 && isfinite(Ly_m) && Ly_m > 0) && ...
        isfield(eval_rec, 'stage_rects') && numel(eval_rec.stage_rects) >= stage_idx && ...
        ~isempty(eval_rec.stage_rects{stage_idx})
    rects = eval_rec.stage_rects{stage_idx};
    if size(rects,2) == 4 && all(isfinite(rects(:)))
        Lx_m = max(rects(:,2)) - min(rects(:,1));
        Ly_m = max(rects(:,4)) - min(rects(:,3));
    end
end
end

function mkdir_if_needed_ms(p)
if ~exist(p, 'dir')
    mkdir(p);
end
end
