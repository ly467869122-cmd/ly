﻿function results = opt_stage_leglength_fixed_layout(spec_in)
% Optimize per-stage thermoelectric leg length for one fixed candidate or regular-grid layout.
% This is the per-stage mode: spec.fp_t_stage can become [L1 L2 ... LN].
% The public field spec.height_opt.height_grid_mm is retained and means
% candidate thermoelectric leg lengths / particle heights in mm.
%
% Candidate per-stage leg length optimization:
%   spec = struct();
%   spec.source.results_mat = 'D:\...\results_data.mat';
%   spec.target.layout_mode = 'candidate';
%   spec.target.candidate_id = 23;
%   spec.height_opt.height_grid_mm = [1.2 1.5 1.8 2.1];
%   results_H = opt_stage_leglength_fixed_layout(spec);
%
% Regular-grid per-stage leg length optimization:
%   spec = struct();
%   spec.target.layout_mode = 'regular_grid';
%   spec.stage_count = 5;
%   spec.fixed_n = [322 110 42 14 6];
%   spec.current.I_init = 3.0;
%   spec.height_opt.height_grid_mm = [1.2 1.5 1.8 2.1];
%   results_H = opt_stage_leglength_fixed_layout(spec);
%
% Use optimized current:
%   spec.current.I_fixed = results_I.best_I;
%   results_H = opt_stage_leglength_fixed_layout(spec);
%
% Or load optimized current from result MAT:
%   spec.current.source_result_mat = 'D:\...\current_opt_result.mat';
%   results_H = opt_stage_leglength_fixed_layout(spec);
%
% Required input for candidate mode:
%   spec.source.results_mat              path to results_data.mat from layout search
%   spec.target.candidate_id             candidate_id to optimize
%
% Optional input:
%   spec.target.layout_mode              'candidate' (default) or 'regular_grid'
%   spec.current.I_fixed                 explicit fixed current in A
%   spec.current.source_result_mat       current_opt_result.mat containing scan.best_I
%   spec.height_opt.height_grid_mm       candidate leg lengths in mm for each stage
%   spec.height_opt.rounds               beam coordinate rounds, default 2
%   spec.height_opt.beamK                retained beam count after each stage, default 2
%   spec.height_opt.qc_abs_tol           default 1e-3 W
%   spec.height_opt.qc_rel_tol           default 1e-3
%   spec.height_opt.dedup_enable         skip repeated fp_t_stage trials, default true
%   spec.height_opt.dedup_tol_m          height-combination key tolerance in m, default 1e-9
%   spec.use_parallel                    auto-start/reuse parallel pool, default true
%   spec.parallel.pool_workers           requested parallel workers, default 64
%   spec.output.output_dir               default sibling run directory beside results_data.mat
if nargin < 1 || isempty(spec_in)
    spec_in = struct();
end
spec = apply_height_opt_defaults_ms(spec_in);
mkdir_if_needed_ms(spec.output.output_dir);

fprintf('\n========== opt_stage_leglength_fixed_layout ==========' );
if is_regular_grid_mode_ms(spec)
    fprintf('\nSource MAT: <regular_grid generated from spec>\n');
else
    fprintf('\nSource MAT: %s\n', spec.source.results_mat);
end
fprintf('Layout mode: %s\n', spec.target.layout_mode);
fprintf('Target candidate_id=%d\n', spec.target.candidate_id);
if isempty(spec.height_opt.height_grid_mm)
    fprintf('height_grid_mm=<auto from base_fp_t_stage * [0.70 0.85 1.00 1.15 1.30]>\n');
else
    fprintf('height_grid_mm=%s\n', vec_to_inline_str_ms(spec.height_opt.height_grid_mm));
end
fprintf('rounds=%d, beamK=%d\n', spec.height_opt.rounds, spec.height_opt.beamK);
fprintf('Output: %s\n', spec.output.output_dir);

if is_regular_grid_mode_ms(spec)
    [base_spec, cand, source_name] = build_regular_grid_source_ms(spec);
else
    data = load(spec.source.results_mat);
    [base_spec, cand, source_name] = select_candidate_from_results_ms(data, spec.target.candidate_id);
end
[fixed_I, current_source] = resolve_fixed_current_ms(cand, base_spec, spec);
base_spec = merge_height_opt_runtime_spec_ms(base_spec, spec, fixed_I);

[base_spec, G] = optimize_layout_multistage0411_shared_params(base_spec, 'current');
G.I = fixed_I;
parallel_requested = logical(base_spec.use_parallel);
base_spec = ensure_height_opt_parallel_pool_ms(base_spec);
spec.use_parallel = base_spec.use_parallel;
spec.parallel = base_spec.parallel;
fprintf('Parallel: requested=%d, used=%d, workers=%d\n', ...
    parallel_requested, logical(base_spec.use_parallel), ...
    get_parallel_worker_count_ms());
fprintf('Dedup: enabled=%d, tol_m=%.3g\n', ...
    logical(spec.height_opt.dedup_enable), spec.height_opt.dedup_tol_m);

scan = run_height_beam_search_for_candidate_ms(cand, G, base_spec, spec, fixed_I);
scan.current_source = current_source;
write_height_opt_outputs_ms(scan, cand, base_spec, spec, source_name);

results = struct();
results.spec = spec;
results.base_spec = base_spec;
results.source_name = source_name;
results.candidate = cand;
results.scan_rows = scan.rows;
results.scan_records = scan.records;
results.beam_history = scan.beam_history;
results.best = scan.best;
results.best_fp_t_stage = scan.best_fp_t_stage;
results.fixed_I = fixed_I;
results.current_source = current_source;
results.has_feasible_solution = scan.has_feasible_solution;
results.output_dir = spec.output.output_dir;
fprintf('Done. has_feasible_solution=%d, fixed_I=%.6g A, best_DeltaTN=%.6g K, best_fp_t_stage_mm=%s\n', ...
    scan.has_feasible_solution, fixed_I, scan.best_DeltaTN_actual, vec_to_inline_str_ms(scan.best_fp_t_stage * 1e3));
fprintf('============================================================\n\n');
end

function spec = apply_height_opt_defaults_ms(spec)
if ~isstruct(spec), error('spec_in must be a struct.'); end
if ~isfield(spec, 'source') || ~isstruct(spec.source), spec.source = struct(); end
if ~isfield(spec, 'target') || ~isstruct(spec.target), spec.target = struct(); end
if ~isfield(spec, 'height_opt') || ~isstruct(spec.height_opt), spec.height_opt = struct(); end
if ~isfield(spec, 'output') || ~isstruct(spec.output), spec.output = struct(); end
if ~isfield(spec, 'parallel') || ~isstruct(spec.parallel), spec.parallel = struct(); end
if ~isfield(spec, 'current') || ~isstruct(spec.current), spec.current = struct(); end
if ~isfield(spec.target, 'layout_mode') || isempty(spec.target.layout_mode)
    spec.target.layout_mode = 'candidate';
end
spec.target.layout_mode = lower(strtrim(char(string(spec.target.layout_mode))));
if ~any(strcmp(spec.target.layout_mode, {'candidate', 'regular_grid'}))
    error('spec.target.layout_mode must be ''candidate'' or ''regular_grid''.');
end
if ~isfield(spec, 'use_parallel') || isempty(spec.use_parallel)
    spec.use_parallel = true;
end
spec.use_parallel = any(logical(spec.use_parallel(:)));
if is_regular_grid_mode_ms(spec)
    spec.target.candidate_id = -1;
    if ~isfield(spec.source, 'results_mat')
        spec.source.results_mat = '';
    end
else
    if ~isfield(spec.source, 'results_mat') || isempty(spec.source.results_mat)
        error('spec.source.results_mat is required for candidate mode.');
    end
    spec.source.results_mat = char(string(spec.source.results_mat));
    if exist(spec.source.results_mat, 'file') ~= 2
        error('results_data.mat not found: %s', spec.source.results_mat);
    end
    if ~isfield(spec.target, 'candidate_id') || isempty(spec.target.candidate_id) || ~isfinite(spec.target.candidate_id)
        error('spec.target.candidate_id is required for candidate mode.');
    end
    spec.target.candidate_id = round(double(spec.target.candidate_id));
end
if ~isfield(spec.height_opt, 'rounds') || isempty(spec.height_opt.rounds) || ~isfinite(spec.height_opt.rounds)
    spec.height_opt.rounds = 2;
end
if ~isfield(spec.height_opt, 'beamK') || isempty(spec.height_opt.beamK) || ~isfinite(spec.height_opt.beamK)
    spec.height_opt.beamK = 2;
end
if ~isfield(spec.height_opt, 'qc_abs_tol') || isempty(spec.height_opt.qc_abs_tol) || ~isfinite(spec.height_opt.qc_abs_tol)
    spec.height_opt.qc_abs_tol = 1e-3;
end
if ~isfield(spec.height_opt, 'qc_rel_tol') || isempty(spec.height_opt.qc_rel_tol) || ~isfinite(spec.height_opt.qc_rel_tol)
    spec.height_opt.qc_rel_tol = 1e-3;
end
if ~isfield(spec.height_opt, 'dedup_enable') || isempty(spec.height_opt.dedup_enable)
    spec.height_opt.dedup_enable = true;
end
if ~isfield(spec.height_opt, 'dedup_tol_m') || isempty(spec.height_opt.dedup_tol_m) || ~isfinite(spec.height_opt.dedup_tol_m)
    spec.height_opt.dedup_tol_m = 1e-9;
end
spec.height_opt.rounds = max(1, round(double(spec.height_opt.rounds)));
spec.height_opt.beamK = max(1, round(double(spec.height_opt.beamK)));
spec.height_opt.qc_abs_tol = max(0, double(spec.height_opt.qc_abs_tol));
spec.height_opt.qc_rel_tol = max(0, double(spec.height_opt.qc_rel_tol));
spec.height_opt.dedup_enable = any(logical(spec.height_opt.dedup_enable(:)));
spec.height_opt.dedup_tol_m = max(eps, double(spec.height_opt.dedup_tol_m));
if ~isfield(spec.parallel, 'pool_workers') || isempty(spec.parallel.pool_workers) || ~isfinite(spec.parallel.pool_workers)
    spec.parallel.pool_workers = 64;
end
spec.parallel.pool_workers = max(1, round(double(spec.parallel.pool_workers)));
if isfield(spec.height_opt, 'height_grid_mm') && ~isempty(spec.height_opt.height_grid_mm)
    h = unique(double(spec.height_opt.height_grid_mm(:).'));
    h = h(isfinite(h) & h > 0);
    if isempty(h)
        error('spec.height_opt.height_grid_mm must contain at least one positive finite height.');
    end
    spec.height_opt.height_grid_mm = h;
else
    spec.height_opt.height_grid_mm = [];
end
if ~isfield(spec.output, 'output_dir') || isempty(spec.output.output_dir)
    if is_regular_grid_mode_ms(spec)
        base_dir = fullfile(pwd, 'stage_leglength_opt_regular_grid');
        target_tag = 'regular_grid';
    else
        src_dir = fileparts(spec.source.results_mat);
        base_dir = fullfile(src_dir, sprintf('stage_leglength_opt_candidate%d', spec.target.candidate_id));
        target_tag = sprintf('candidate%d', spec.target.candidate_id);
    end
    stamp = datestr(now, 'yyyymmdd_HHMMSS');
    if isempty(spec.height_opt.height_grid_mm)
        h_tag = 'auto';
    else
        h_tag = numeric_vec_tag_ms(spec.height_opt.height_grid_mm);
    end
    run_name = sprintf('stage_leglength_opt_%s_h%s_%s', target_tag, h_tag, stamp);
    spec.output.output_dir = fullfile(base_dir, run_name);
end
spec.output.output_dir = char(string(spec.output.output_dir));
end

function tf = is_regular_grid_mode_ms(spec)
tf = isfield(spec, 'target') && isstruct(spec.target) && ...
    isfield(spec.target, 'layout_mode') && strcmpi(char(string(spec.target.layout_mode)), 'regular_grid');
end

function spec = apply_current_opt_defaults_ms(spec)
if ~isstruct(spec), error('spec_in must be a struct.'); end
if ~isfield(spec, 'source') || ~isstruct(spec.source), spec.source = struct(); end
if ~isfield(spec, 'target') || ~isstruct(spec.target), spec.target = struct(); end
if ~isfield(spec, 'current_opt') || ~isstruct(spec.current_opt), spec.current_opt = struct(); end
if ~isfield(spec, 'output') || ~isstruct(spec.output), spec.output = struct(); end
if ~isfield(spec.source, 'results_mat') || isempty(spec.source.results_mat)
    error('spec.source.results_mat is required.');
end
spec.source.results_mat = char(string(spec.source.results_mat));
if exist(spec.source.results_mat, 'file') ~= 2
    error('results_data.mat not found: %s', spec.source.results_mat);
end
if ~isfield(spec.target, 'candidate_id') || isempty(spec.target.candidate_id) || ~isfinite(spec.target.candidate_id)
    error('spec.target.candidate_id is required.');
end
spec.target.candidate_id = round(double(spec.target.candidate_id));
if ~isfield(spec.current_opt, 'I_list') || isempty(spec.current_opt.I_list)
    spec.current_opt.I_list = 1:0.2:6;
end
spec.current_opt.I_list = unique(double(spec.current_opt.I_list(:).'));
spec.current_opt.I_list = spec.current_opt.I_list(isfinite(spec.current_opt.I_list) & spec.current_opt.I_list > 0);
if isempty(spec.current_opt.I_list)
    error('spec.current_opt.I_list must contain at least one positive finite current.');
end
if ~isfield(spec.current_opt, 'qc_abs_tol') || isempty(spec.current_opt.qc_abs_tol) || ~isfinite(spec.current_opt.qc_abs_tol)
    spec.current_opt.qc_abs_tol = 1e-3;
end
if ~isfield(spec.current_opt, 'qc_rel_tol') || isempty(spec.current_opt.qc_rel_tol) || ~isfinite(spec.current_opt.qc_rel_tol)
    spec.current_opt.qc_rel_tol = 1e-3;
end
spec.current_opt.qc_abs_tol = max(0, double(spec.current_opt.qc_abs_tol));
spec.current_opt.qc_rel_tol = max(0, double(spec.current_opt.qc_rel_tol));
if ~isfield(spec.output, 'output_dir') || isempty(spec.output.output_dir)
    src_dir = fileparts(spec.source.results_mat);
    base_dir = fullfile(src_dir, sprintf('current_opt_candidate%d', spec.target.candidate_id));
    stamp = datestr(now, 'yyyymmdd_HHMMSS');
    i_min_tag = numeric_tag_ms(min(spec.current_opt.I_list));
    i_max_tag = numeric_tag_ms(max(spec.current_opt.I_list));
    i_count = numel(spec.current_opt.I_list);
    run_name = sprintf('current_opt_candidate%d_I%s-%s_n%d_%s', ...
        spec.target.candidate_id, i_min_tag, i_max_tag, i_count, stamp);
    spec.output.output_dir = fullfile(base_dir, run_name);
end
spec.output.output_dir = char(string(spec.output.output_dir));
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

function [base_spec, cand, source_name] = build_regular_grid_source_ms(spec)
[base_spec, ~] = optimize_layout_multistage0411_shared_params(spec, 'layout');
base_spec.target = spec.target;
base_spec.height_opt = spec.height_opt;
base_spec.output = spec.output;
cand = build_regular_grid_candidate_from_spec_ms(base_spec);
source_name = 'regular_grid';
end

function cand = build_regular_grid_candidate_from_spec_ms(spec)
N = spec.stage_count;
n_vec = normalize_even_n_vec_ms(spec.fixed_n, N);
cand = empty_candidate_struct_ms(N);
cand.candidate_id = -1;
cand.layout_method = 'k25_standard_fixed50_baseline';
cand.stage_modes = repmat({'standard_fixed50'}, 1, N);
cand.stage_methods = repmat({'k25_standard_fixed50_baseline'}, 1, N);
cand.stage_trends = repmat({'neutral'}, 1, N);
cand.symmetry_mode = 'none';
cand.edge_pattern_mode = 'two_sides_to_center';
cand.n = n_vec;
cand.ratios = n_vec(2:end) ./ max(n_vec(1:end-1), eps);
cand.I_opt = spec.current.I_init;
cand.count_rank = NaN;

rects_u = cell(1, N);
Lx_u = NaN(1, N);
Ly_u = NaN(1, N);
for k = 1:N
    [rects_u{k}, Lx_u(k), Ly_u(k)] = make_regular_grid_stage_rects_ms(n_vec(k), spec);
end

gap = zeros(1, max(0, N-1));
if isfield(spec, 'geometry') && isstruct(spec.geometry) && isfield(spec.geometry, 'pyramid_gap_min')
    gap_raw = double(spec.geometry.pyramid_gap_min(:).');
    gap(1:min(numel(gap), numel(gap_raw))) = gap_raw(1:min(numel(gap), numel(gap_raw)));
end
cand.Lx(N) = Lx_u(N);
cand.Ly(N) = Ly_u(N);
for k = N-1:-1:1
    cand.Lx(k) = max(Lx_u(k), cand.Lx(k+1) + gap(k));
    cand.Ly(k) = max(Ly_u(k), cand.Ly(k+1) + gap(k));
end
for k = 1:N
    cand.stage_rects{k} = rects_u{k};
    cand.cov(k) = n_vec(k) * spec.fp_w * spec.fp_h / max(cand.Lx(k) * cand.Ly(k), eps);
end
cand.s_dense = spec.fp_w * 1.5;
cand.s_sparse = cand.s_dense;
cand.expo = 1.0;
cand.anis_ratio = 1.0;
cand.gamma = NaN;
cand.method_anchor_stage = NaN;
cand.shape_mode = '';
cand.stage_anis = NaN(1, N);
cand.ring_radius_ratio = NaN;
cand.ring_width_ratio = NaN;
cand.band_width_ratio = NaN;
cand.corner_bias = NaN;
cand.jitter_ratio = NaN;
cand.jitter_seed = NaN;
cand.spacing_ratio = 1.0;
cand.contrast_score = 1.0;
cand.mode_prior_score = 0;
cand.lmax_relaxed = false;
cand.lmax_relax_ratio = 1.0;
end

function cand = empty_candidate_struct_ms(N)
cand = struct('candidate_id', NaN, 'layout_method', '', ...
    'stage_modes', {cell(1,N)}, 'stage_methods', {cell(1,N)}, ...
    'stage_trends', {cell(1,N)}, 'symmetry_mode', '', 'edge_pattern_mode', '', ...
    's_dense', NaN, 's_sparse', NaN, 'expo', NaN, 'anis_ratio', NaN, ...
    'gamma', NaN, 'method_anchor_stage', NaN, 'shape_mode', '', ...
    'stage_anis', NaN(1,N), 'ring_radius_ratio', NaN, 'ring_width_ratio', NaN, ...
    'band_width_ratio', NaN, 'corner_bias', NaN, 'jitter_ratio', NaN, ...
    'jitter_seed', NaN, 'spacing_ratio', NaN, 'contrast_score', NaN, ...
    'mode_prior_score', NaN, 'lmax_relaxed', false, 'lmax_relax_ratio', 1.0, ...
    'count_rank', NaN, 'n', NaN(1,N), 'ratios', NaN(1,max(0,N-1)), ...
    'I_opt', NaN, 'Lx', NaN(1,N), 'Ly', NaN(1,N), 'cov', NaN(1,N), ...
    'stage_rects', {cell(1,N)});
end

function n_vec = normalize_even_n_vec_ms(n_raw, N)
n_vec = round(double(n_raw(:).'));
if numel(n_vec) ~= N || any(~isfinite(n_vec)) || any(n_vec < 2)
    error('spec.fixed_n must contain %d positive particle counts for regular_grid mode.', N);
end
if any(mod(n_vec, 2) ~= 0)
    error('spec.fixed_n must contain even particle counts for the TEC pair model.');
end
end

function [rects, Lx_box, Ly_box] = make_regular_grid_stage_rects_ms(n_target, spec)
gap_x = 0.5 * spec.fp_w;
gap_y = 0.5 * spec.fp_h;
margin_x = 0.5 * spec.fp_w;
margin_y = 0.5 * spec.fp_h;
pitch_x = spec.fp_w + gap_x;
pitch_y = spec.fp_h + gap_y;
[nx, ny] = choose_regular_grid_shape_ms(n_target, spec.fp_w, spec.fp_h);
pts = fill_regular_grid_two_sides_to_center_ms(nx, ny, n_target, pitch_x, pitch_y);
if size(pts, 1) ~= n_target
    error('regular_grid stage layout failed for n=%d.', n_target);
end
rects = centers_to_rects_ms(pts, spec.fp_w, spec.fp_h);
Lx_box = max(rects(:,2)) - min(rects(:,1)) + 2 * margin_x;
Ly_box = max(rects(:,4)) - min(rects(:,3)) + 2 * margin_y;
end

function [nx_best, ny_best] = choose_regular_grid_shape_ms(n_target, fp_w, fp_h)
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
    for ic = 1:nx
        p = p + 1;
        if p > n_target
            break;
        end
        pts(p,:) = [xvals(col_order(ic)), y];
    end
    row_ptr = row_ptr + 1;
end
end

function order = two_sides_to_center_order_ms(n)
lo = 1;
hi = n;
order = zeros(1, n);
p = 0;
while lo <= hi
    p = p + 1;
    order(p) = lo;
    lo = lo + 1;
    if lo <= hi
        p = p + 1;
        order(p) = hi;
        hi = hi - 1;
    end
end
end

function rects = centers_to_rects_ms(pts, fp_w, fp_h)
rects = [pts(:,1) - fp_w/2, pts(:,1) + fp_w/2, pts(:,2) - fp_h/2, pts(:,2) + fp_h/2];
end

function [base_spec, cand, source_name] = select_candidate_from_results_ms(data, candidate_id)
if ~isfield(data, 'spec') || ~isstruct(data.spec)
    error('results_data.mat does not contain a valid spec struct.');
end
base_spec = data.spec;
sources = {'top_opt_full', 'final_valid', 'all_eval', 'best_cand'};
for si = 1:numel(sources)
    name = sources{si};
    if ~isfield(data, name) || isempty(data.(name))
        continue;
    end
    arr = data.(name);
    if isstruct(arr) && isfield(arr, 'candidate_id')
        idx = find([arr.candidate_id] == candidate_id, 1, 'first');
        if ~isempty(idx)
            cand = arr(idx);
            source_name = name;
            if ~has_valid_candidate_geometry_ms(cand, base_spec.stage_count)
                error('candidate_id %d was found in %s, but it does not contain reusable stage_rects geometry.', ...
                    candidate_id, name);
            end
            return;
        end
    end
end
error('candidate_id %d not found in top_opt_full/final_valid/all_eval/best_cand.', candidate_id);
end

function base_spec = merge_current_opt_runtime_spec_ms(base_spec, spec)
if ~isfield(base_spec, 'output') || ~isstruct(base_spec.output), base_spec.output = struct(); end
if ~isfield(base_spec.output, 'plot') || ~isstruct(base_spec.output.plot), base_spec.output.plot = struct(); end
if ~isfield(base_spec, 'current') || ~isstruct(base_spec.current), base_spec.current = struct(); end
base_spec.output.output_dir = spec.output.output_dir;
base_spec.output.plot.save_overview = true;
base_spec.output.plot.save_stage_separate = true;
base_spec.output.plot.use_global_axis = true;
base_spec.output.plot.show_mesh_edges = false;
base_spec.output.plot.view_mode = '2d';
base_spec.current.I_init = spec.current_opt.I_list(1);
if ~isfield(base_spec, 'use_parallel'), base_spec.use_parallel = false; end
base_spec.use_parallel = false;
if ~isfield(base_spec, 'parallel') || ~isstruct(base_spec.parallel), base_spec.parallel = struct(); end
if ~isfield(base_spec.parallel, 'eval_min_tasks'), base_spec.parallel.eval_min_tasks = inf; end
if ~isfield(base_spec, 'hard_constraints') || ~isstruct(base_spec.hard_constraints), base_spec.hard_constraints = struct(); end
if ~isfield(base_spec.hard_constraints, 'enforce_spacing_min_edge_gap'), base_spec.hard_constraints.enforce_spacing_min_edge_gap = true; end
end

function [fixed_I, source] = resolve_fixed_current_ms(cand, base_spec, spec)
fixed_I = NaN;
source = 'none';
if isstruct(spec) && isfield(spec, 'current') && isstruct(spec.current)
    if isfield(spec.current, 'I_fixed') && ~isempty(spec.current.I_fixed) && ...
            isnumeric(spec.current.I_fixed) && isfinite(spec.current.I_fixed(1)) && spec.current.I_fixed(1) > 0
        fixed_I = spec.current.I_fixed(1);
        source = 'spec.current.I_fixed';
        return;
    end
    if isfield(spec.current, 'source_result_mat') && ~isempty(spec.current.source_result_mat)
        src_mat = char(string(spec.current.source_result_mat));
        if exist(src_mat, 'file') ~= 2
            error('spec.current.source_result_mat not found: %s', src_mat);
        end
        data = load(src_mat);
        if isfield(data, 'scan') && isstruct(data.scan) && isfield(data.scan, 'best_I') && ...
                ~isempty(data.scan.best_I) && isfinite(data.scan.best_I(1)) && data.scan.best_I(1) > 0
            fixed_I = data.scan.best_I(1);
            source = 'spec.current.source_result_mat:scan.best_I';
            return;
        end
        error('spec.current.source_result_mat does not contain a positive scan.best_I: %s', src_mat);
    end
end
if isstruct(cand) && isfield(cand, 'I_opt') && ~isempty(cand.I_opt) && ...
        isnumeric(cand.I_opt) && isfinite(cand.I_opt(1)) && cand.I_opt(1) > 0
    fixed_I = cand.I_opt(1);
    source = 'cand.I_opt';
    return;
end
if isstruct(base_spec) && isfield(base_spec, 'current') && isstruct(base_spec.current) && ...
        isfield(base_spec.current, 'I_init') && ~isempty(base_spec.current.I_init) && ...
        isfinite(base_spec.current.I_init(1)) && base_spec.current.I_init(1) > 0
    fixed_I = base_spec.current.I_init(1);
    source = 'base_spec.current.I_init';
    return;
end
error('Unable to resolve a positive fixed current from I_fixed, source_result_mat, candidate.I_opt, or base_spec.current.I_init.');
end

function base_spec = merge_height_opt_runtime_spec_ms(base_spec, spec, fixed_I)
if ~isfield(base_spec, 'output') || ~isstruct(base_spec.output), base_spec.output = struct(); end
if ~isfield(base_spec.output, 'plot') || ~isstruct(base_spec.output.plot), base_spec.output.plot = struct(); end
if ~isfield(base_spec, 'current') || ~isstruct(base_spec.current), base_spec.current = struct(); end
base_spec.output.output_dir = spec.output.output_dir;
base_spec.output.plot.save_overview = true;
base_spec.output.plot.save_stage_separate = true;
base_spec.output.plot.use_global_axis = true;
base_spec.output.plot.show_mesh_edges = false;
base_spec.output.plot.view_mode = '2d';
base_spec.current.I_init = fixed_I;
base_spec = apply_explicit_leglength_override_ms(base_spec, spec);
base_spec.use_parallel = logical(spec.use_parallel);
if ~isfield(base_spec, 'parallel') || ~isstruct(base_spec.parallel), base_spec.parallel = struct(); end
if isfield(spec, 'parallel') && isstruct(spec.parallel)
    base_spec.parallel = merge_struct_recursive_ms(base_spec.parallel, spec.parallel);
end
if ~isfield(base_spec.parallel, 'eval_min_tasks'), base_spec.parallel.eval_min_tasks = inf; end
if ~isfield(base_spec, 'hard_constraints') || ~isstruct(base_spec.hard_constraints), base_spec.hard_constraints = struct(); end
if ~isfield(base_spec.hard_constraints, 'enforce_spacing_min_edge_gap'), base_spec.hard_constraints.enforce_spacing_min_edge_gap = true; end
end

function base_spec = apply_explicit_leglength_override_ms(base_spec, spec)
if ~isstruct(spec)
    return;
end
N = 5;
if isfield(base_spec, 'stage_count') && isfinite(base_spec.stage_count)
    N = max(1, round(base_spec.stage_count));
elseif isfield(spec, 'stage_count') && isfinite(spec.stage_count)
    N = max(1, round(spec.stage_count));
end
if isfield(spec, 'fp_t_stage') && ~isempty(spec.fp_t_stage)
    base_spec.fp_t_stage = fit_len_vec_ms(double(spec.fp_t_stage), N, 1.5e-3);
    base_spec.fp_t = base_spec.fp_t_stage;
elseif isfield(spec, 'fp_t') && ~isempty(spec.fp_t)
    base_spec.fp_t_stage = fit_len_vec_ms(double(spec.fp_t), N, 1.5e-3);
    base_spec.fp_t = base_spec.fp_t_stage;
end
end

function scan = run_height_beam_search_for_candidate_ms(cand, G, base_spec, spec, fixed_I)
N = base_spec.stage_count;
base_fp_t = resolve_base_fp_t_stage_ms(base_spec, N);
height_grid_by_stage_mm = resolve_height_grid_by_stage_mm_ms(spec.height_opt.height_grid_mm, base_fp_t, N);
height_grid_flat_mm = unique([height_grid_by_stage_mm{:}]);
qc_tol = max(spec.height_opt.qc_abs_tol, spec.height_opt.qc_rel_tol * abs(base_spec.targets.Qc_target_last));
beamK = spec.height_opt.beamK;
rounds = spec.height_opt.rounds;

beam = empty_height_beam_ms(N);
beam(1).rank = 1;
beam(1).fp_t_stage = base_fp_t;
beam(1).score_row = empty_height_scan_row_ms(N);
beam(1).source_scan_idx = 0;

rows = repmat(empty_height_scan_row_ms(N), 0, 1);
records = repmat(empty_eval_struct_ms(N), 0, 1);
beam_history = repmat(empty_beam_history_row_ms(N), 0, 1);
scan_idx = 0;
dedup_enabled = logical(spec.height_opt.dedup_enable);
dedup_tol_m = spec.height_opt.dedup_tol_m;
seen_keys = containers.Map('KeyType', 'char', 'ValueType', 'logical');
duplicate_trial_count = 0;
evaluated_trial_count = 0;

for round_idx = 1:rounds
    fprintf('[HeightBeam] round %d/%d start, beam_count=%d\n', round_idx, rounds, numel(beam));
    for stage_idx = 1:N
        trials = expand_height_trials_for_stage_ms(beam, height_grid_by_stage_mm{stage_idx}, round_idx, stage_idx, N);
        if isempty(trials)
            continue;
        end
        raw_trial_count = numel(trials);
        [trials, skipped_duplicate, seen_keys] = filter_height_trials_unique_ms(trials, seen_keys, dedup_enabled, dedup_tol_m);
        duplicate_trial_count = duplicate_trial_count + skipped_duplicate;
        fprintf('[HeightBeam] round %d stage %d: trials=%d, unique=%d, skipped_duplicate=%d\n', ...
            round_idx, stage_idx, raw_trial_count, numel(trials), skipped_duplicate);
        if isempty(trials)
            continue;
        end
        [trial_rows, trial_records] = evaluate_height_trial_batch_ms(trials, cand, base_spec, fixed_I, qc_tol, base_fp_t, scan_idx);
        scan_idx = scan_idx + numel(trial_rows);
        evaluated_trial_count = evaluated_trial_count + numel(trial_rows);
        rows = [rows; trial_rows(:)]; %#ok<AGROW>
        records = [records; trial_records(:)]; %#ok<AGROW>
        [~, ord] = rank_height_scan_rows_ms(trial_rows);
        keep_n = min(beamK, numel(ord));
        next_beam = repmat(empty_height_beam_ms(N), keep_n, 1);
        for k = 1:keep_n
            src_idx = ord(k);
            next_beam(k).rank = k;
            next_beam(k).fp_t_stage = trial_rows(src_idx).fp_t_stage_numeric;
            next_beam(k).score_row = trial_rows(src_idx);
            next_beam(k).source_scan_idx = trial_rows(src_idx).scan_idx;
            beam_history(end+1,1) = build_beam_history_row_ms(round_idx, stage_idx, k, trial_rows(src_idx), base_fp_t); %#ok<AGROW>
        end
        beam = next_beam;
    end
end

[~, all_ord] = rank_height_scan_rows_ms(rows);
scan = struct();
scan.rows = rows;
scan.records = records;
scan.beam_history = beam_history;
scan.qc_tol = qc_tol;
scan.height_grid_by_stage_mm = height_grid_by_stage_mm;
scan.height_grid_mm = height_grid_flat_mm;
scan.base_fp_t_stage = base_fp_t;
scan.fixed_I = fixed_I;
scan.parallel_used = logical(base_spec.use_parallel) && get_parallel_worker_count_ms() > 0;
scan.parallel_workers = get_parallel_worker_count_ms();
scan.dedup_enabled = dedup_enabled;
scan.dedup_tol_m = dedup_tol_m;
scan.duplicate_trial_count = duplicate_trial_count;
scan.evaluated_trial_count = evaluated_trial_count;
scan.has_feasible_solution = ~isempty(rows) && any([rows.target_met]);
scan.best = empty_eval_struct_ms(N);
scan.best_fp_t_stage = NaN(1, N);
scan.best_DeltaTN_actual = NaN;
if ~isempty(all_ord)
    if scan.has_feasible_solution
        feasible_ord = all_ord([rows(all_ord).target_met]);
        best_row_idx = feasible_ord(1);
    else
        best_row_idx = all_ord(1);
    end
    scan.best_fp_t_stage = rows(best_row_idx).fp_t_stage_numeric;
    scan.best_DeltaTN_actual = rows(best_row_idx).DeltaTN_actual;
    scan.best = records(best_row_idx);
end
end

function base_fp_t = resolve_base_fp_t_stage_ms(base_spec, N)
if isfield(base_spec, 'fp_t_stage') && ~isempty(base_spec.fp_t_stage)
    base_fp_t = fit_len_vec_ms(double(base_spec.fp_t_stage), N, 1.5e-3);
elseif isfield(base_spec, 'fp_t') && ~isempty(base_spec.fp_t)
    base_fp_t = fit_len_vec_ms(double(base_spec.fp_t), N, 1.5e-3);
else
    base_fp_t = 1.5e-3 * ones(1, N);
end
base_fp_t = max(1e-9, base_fp_t(:).');
end

function grids = resolve_height_grid_by_stage_mm_ms(height_grid_mm, base_fp_t, N)
grids = cell(1, N);
if ~isempty(height_grid_mm)
    h = unique(double(height_grid_mm(:).'));
    h = h(isfinite(h) & h > 0);
    for k = 1:N
        grids{k} = h;
    end
    return;
end
factors = [0.70, 0.85, 1.00, 1.15, 1.30];
for k = 1:N
    hk = unique(base_fp_t(k) * 1e3 * factors);
    grids{k} = hk(isfinite(hk) & hk > 0);
end
end

function beam = empty_height_beam_ms(N)
beam = struct('rank', NaN, 'fp_t_stage', NaN(1,N), ...
    'score_row', empty_height_scan_row_ms(N), 'source_scan_idx', NaN);
end

function trials = expand_height_trials_for_stage_ms(beam, height_grid_mm, round_idx, stage_idx, N)
trials = struct('round_idx', {}, 'stage_idx', {}, 'parent_beam_rank', {}, ...
    'fp_t_stage', {}, 'changed_stage_height_mm', {});
if isempty(beam) || isempty(height_grid_mm)
    return;
end
for b = 1:numel(beam)
    base_h = beam(b).fp_t_stage;
    for hi = 1:numel(height_grid_mm)
        fp_t = fit_len_vec_ms(base_h, N, base_h);
        fp_t(stage_idx) = height_grid_mm(hi) * 1e-3;
        t = struct();
        t.round_idx = round_idx;
        t.stage_idx = stage_idx;
        t.parent_beam_rank = b;
        t.fp_t_stage = fp_t;
        t.changed_stage_height_mm = height_grid_mm(hi);
        trials(end+1,1) = t; %#ok<AGROW>
    end
end
end

function [rows, records] = evaluate_height_trial_batch_ms(trials, cand, base_spec, fixed_I, qc_tol, base_fp_t, scan_idx0)
N = base_spec.stage_count;
rows = repmat(empty_height_scan_row_ms(N), numel(trials), 1);
records = repmat(empty_eval_struct_ms(N), numel(trials), 1);
use_par = logical(base_spec.use_parallel) && numel(trials) >= 2 && has_active_parallel_pool_ms();
if use_par
    parfor i = 1:numel(trials)
        [records(i), rows(i)] = evaluate_height_trial_one_ms(trials(i), cand, base_spec, fixed_I, qc_tol, base_fp_t, scan_idx0 + i);
    end
else
    for i = 1:numel(trials)
        [records(i), rows(i)] = evaluate_height_trial_one_ms(trials(i), cand, base_spec, fixed_I, qc_tol, base_fp_t, scan_idx0 + i);
    end
end
end

function [rec, row] = evaluate_height_trial_one_ms(trial, cand, base_spec, fixed_I, qc_tol, base_fp_t, scan_idx)
base_i = base_spec;
base_i.fp_t_stage = trial.fp_t_stage(:).';
base_i.fp_t = base_i.fp_t_stage;
base_i.current.I_init = fixed_I;
[base_i, G_i] = optimize_layout_multistage0411_shared_params(base_i, 'current');
G_i.I = fixed_I;
cand_i = cand;
cand_i.I_opt = fixed_I;
rec = evaluate_single_candidate_ms(cand_i, G_i, base_i, true);
row = build_height_scan_row_ms(scan_idx, trial, rec, base_i, qc_tol, base_fp_t, fixed_I);
end

function row = empty_height_scan_row_ms(N)
if nargin < 1 || ~isfinite(N), N = 5; end
row = struct('scan_idx', NaN, 'round_idx', NaN, 'stage_idx', NaN, ...
    'parent_beam_rank', NaN, 'candidate_id', NaN, ...
    'fp_t_stage_mm', '[]', 'fp_t_stage_numeric', NaN(1,N), ...
    'changed_stage_height_mm', NaN, 'I_A', NaN, ...
    'success', false, 'target_met', false, 'DeltaTN_actual', NaN, ...
    'DeltaT_target', NaN, 'Qc_last_total', NaN, 'Qc_target_last', NaN, ...
    'Qc_error', NaN, 'Qc_tol', NaN, 'TN_min', NaN, 'TN_mean', NaN, ...
    'DeltaTN_mean', NaN, 'TN_maxmin', NaN, 'newton_rel_max', NaN, ...
    'newton_iters', NaN, 'sum_height_mm', NaN, 'distance_from_base_mm', NaN, ...
    'message', '');
end

function row = build_height_scan_row_ms(scan_idx, trial, rec, spec, qc_tol, base_fp_t, fixed_I)
N = spec.stage_count;
row = empty_height_scan_row_ms(N);
row.scan_idx = scan_idx;
row.round_idx = trial.round_idx;
row.stage_idx = trial.stage_idx;
row.parent_beam_rank = trial.parent_beam_rank;
row.candidate_id = rec.candidate_id;
row.fp_t_stage_numeric = trial.fp_t_stage(:).';
row.fp_t_stage_mm = vec_to_inline_str_ms(row.fp_t_stage_numeric * 1e3);
row.changed_stage_height_mm = trial.changed_stage_height_mm;
row.I_A = fixed_I;
row.success = logical(rec.success);
row.DeltaTN_actual = rec.DeltaTN_actual;
row.DeltaT_target = spec.targets.DeltaT_target;
if isfield(rec, 'Qc_last_total')
    row.Qc_last_total = rec.Qc_last_total;
elseif isfield(rec, 'sumQc_stage') && numel(rec.sumQc_stage) >= N
    row.Qc_last_total = rec.sumQc_stage(N);
end
row.Qc_target_last = spec.targets.Qc_target_last;
row.Qc_error = row.Qc_last_total - row.Qc_target_last;
row.Qc_tol = qc_tol;
row.TN_min = rec.TN_min;
row.TN_mean = rec.TN_mean;
row.DeltaTN_mean = rec.DeltaTN_mean;
row.TN_maxmin = rec.TN_maxmin;
row.newton_rel_max = rec.newton_rel_max;
row.newton_iters = rec.newton_iters;
row.sum_height_mm = sum(row.fp_t_stage_numeric) * 1e3;
row.distance_from_base_mm = norm((row.fp_t_stage_numeric(:).' - base_fp_t(:).') * 1e3);
row.message = char(string(rec.message));
row.target_met = row.success && isfinite(row.DeltaTN_actual) && row.DeltaTN_actual >= spec.targets.DeltaT_target && ...
    isfinite(row.Qc_last_total) && abs(row.Qc_error) <= qc_tol;
end

function [ranked_rows, ord] = rank_height_scan_rows_ms(rows)
ranked_rows = rows;
ord = [];
if isempty(rows)
    return;
end
target = double([rows.target_met].');
dt = [rows.DeltaTN_actual].';
qe = abs([rows.Qc_error].');
sh = [rows.sum_height_mm].';
db = [rows.distance_from_base_mm].';
dt(~isfinite(dt)) = -inf;
qe(~isfinite(qe)) = inf;
sh(~isfinite(sh)) = inf;
db(~isfinite(db)) = inf;
[~, ord] = sortrows([-target, -dt, qe, sh, db]);
ranked_rows = rows(ord);
end

function tf = has_active_parallel_pool_ms()
tf = false;
if ~((exist('gcp', 'file') == 2) || (exist('gcp', 'builtin') == 5))
    return;
end
try
    pool = gcp('nocreate');
    tf = ~isempty(pool) && pool.NumWorkers > 0;
catch
    tf = false;
end
end

function spec = ensure_height_opt_parallel_pool_ms(spec)
if ~isfield(spec, 'parallel') || ~isstruct(spec.parallel)
    spec.parallel = struct();
end
if ~isfield(spec.parallel, 'pool_workers') || isempty(spec.parallel.pool_workers) || ~isfinite(spec.parallel.pool_workers)
    spec.parallel.pool_workers = 64;
end
if ~isfield(spec, 'use_parallel') || isempty(spec.use_parallel)
    spec.use_parallel = true;
end
use_requested = any(logical(spec.use_parallel(:)));
spec.parallel.pool_used = false;
spec.parallel.pool_workers_actual = 0;
if ~use_requested
    spec.use_parallel = false;
    return;
end
has_parallel = ((exist('parpool', 'file') == 2) || (exist('parpool', 'builtin') == 5)) && ...
    ((exist('gcp', 'file') == 2) || (exist('gcp', 'builtin') == 5));
if ~has_parallel
    warning('Height optimization parallel requested but Parallel Computing Toolbox functions are unavailable; falling back to serial.');
    spec.use_parallel = false;
    return;
end
worker_target = max(1, round(double(spec.parallel.pool_workers)));
try
    cluster = parcluster('local');
    if isprop(cluster, 'NumWorkers')
        worker_target = min(worker_target, max(1, round(cluster.NumWorkers)));
    end
catch
end
try
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('local', worker_target);
    elseif pool.NumWorkers ~= worker_target
        fprintf('[HeightBeam] existing parallel pool workers=%d, requested=%d; using existing pool.\n', ...
            pool.NumWorkers, worker_target);
    end
    if isempty(pool)
        spec.use_parallel = false;
    else
        spec.use_parallel = true;
        spec.parallel.pool_used = true;
        spec.parallel.pool_workers_actual = pool.NumWorkers;
    end
catch ME_pool
    warning('Height optimization parallel pool unavailable: %s. Falling back to serial.', ME_pool.message);
    spec.use_parallel = false;
end
end

function nw = get_parallel_worker_count_ms()
nw = 0;
if ~((exist('gcp', 'file') == 2) || (exist('gcp', 'builtin') == 5))
    return;
end
try
    pool = gcp('nocreate');
    if ~isempty(pool) && isprop(pool, 'NumWorkers')
        nw = pool.NumWorkers;
    end
catch
    nw = 0;
end
end

function [unique_trials, skipped_count, seen_keys] = filter_height_trials_unique_ms(trials, seen_keys, dedup_enabled, dedup_tol_m)
skipped_count = 0;
if isempty(trials) || ~dedup_enabled
    unique_trials = trials;
    return;
end
keep = false(numel(trials), 1);
for i = 1:numel(trials)
    key = height_trial_key_ms(trials(i).fp_t_stage, dedup_tol_m);
    if isKey(seen_keys, key)
        skipped_count = skipped_count + 1;
    else
        seen_keys(key) = true;
        keep(i) = true;
    end
end
unique_trials = trials(keep);
end

function key = height_trial_key_ms(fp_t_stage, tol_m)
if nargin < 2 || ~(isfinite(tol_m) && tol_m > 0)
    tol_m = 1e-9;
end
v = round(fp_t_stage(:).' ./ tol_m) .* tol_m;
key = sprintf('%.9g_', v);
end

function row = empty_beam_history_row_ms(N)
if nargin < 1 || ~isfinite(N), N = 5; end
row = struct('round_idx', NaN, 'stage_idx', NaN, 'beam_rank', NaN, ...
    'fp_t_stage_mm', '[]', 'fp_t_stage_numeric', NaN(1,N), ...
    'I_A', NaN, 'target_met', false, 'DeltaTN_actual', NaN, ...
    'Qc_error', NaN, 'sum_height_mm', NaN, ...
    'distance_from_base_mm', NaN, 'source_scan_idx', NaN);
end

function row = build_beam_history_row_ms(round_idx, stage_idx, beam_rank, src, base_fp_t)
N = numel(src.fp_t_stage_numeric);
row = empty_beam_history_row_ms(N);
row.round_idx = round_idx;
row.stage_idx = stage_idx;
row.beam_rank = beam_rank;
row.fp_t_stage_numeric = src.fp_t_stage_numeric;
row.fp_t_stage_mm = src.fp_t_stage_mm;
row.I_A = src.I_A;
row.target_met = src.target_met;
row.DeltaTN_actual = src.DeltaTN_actual;
row.Qc_error = src.Qc_error;
row.sum_height_mm = src.sum_height_mm;
row.distance_from_base_mm = norm((src.fp_t_stage_numeric(:).' - base_fp_t(:).') * 1e3);
row.source_scan_idx = src.scan_idx;
end

function write_height_opt_outputs_ms(scan, cand, base_spec, spec, source_name)
out_dir = spec.output.output_dir;
mkdir_if_needed_ms(out_dir);
rows_out = rmfield_if_exists_ms(scan.rows, 'fp_t_stage_numeric');
beam_out = rmfield_if_exists_ms(scan.beam_history, 'fp_t_stage_numeric');
writetable(struct_to_table_rows_ms(rows_out), fullfile(out_dir, 'height_scan_results.csv'));
writetable(struct_to_table_rows_ms(beam_out), fullfile(out_dir, 'beam_history.csv'));
save(fullfile(out_dir, 'height_opt_result.mat'), 'scan', 'cand', 'base_spec', 'spec', 'source_name');
write_best_height_summary_ms(fullfile(out_dir, 'best_height_summary.txt'), scan, cand, base_spec, spec, source_name);
if scan.has_feasible_solution
    best_dir = fullfile(out_dir, 'best');
    mkdir_if_needed_ms(best_dir);
    best_spec = base_spec;
    best_spec.fp_t_stage = scan.best_fp_t_stage;
    best_spec.fp_t = scan.best_fp_t_stage;
    write_stage_layout_rects_csv_ms(scan.best, fullfile(best_dir, 'layout_stage_rects.csv'), base_spec.stage_count);
    caxis_by_stage = make_caxis_by_stage_ms(scan.best, base_spec.stage_count);
    save_layout_plot_ms(scan.best, fullfile(best_dir, 'layout.png'), best_spec);
    save_layout_stage_plots_ms(scan.best, best_dir, best_spec);
    save_temperature_plot_ms(scan.best, fullfile(best_dir, 'temperature.png'), caxis_by_stage, best_spec);
    save_temperature_stage_plots_ms(scan.best, best_dir, caxis_by_stage, best_spec);
else
    write_text_file_ms(fullfile(out_dir, 'best_unavailable.txt'), ...
        'No particle-height trial met DeltaT_target and Qc_last_total equality tolerance.');
end
end

function write_best_height_summary_ms(out_txt, scan, cand, base_spec, spec, source_name)
fid = fopen(out_txt, 'w');
if fid < 0, warning('Failed to write %s', out_txt); return; end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'fixed-layout particle-height optimization summary\n');
fprintf(fid, 'source_results_mat=%s\n', spec.source.results_mat);
fprintf(fid, 'source_candidate_array=%s\n', char(string(source_name)));
if isfield(spec, 'target') && isstruct(spec.target) && isfield(spec.target, 'layout_mode')
    fprintf(fid, 'target_layout_mode=%s\n', char(string(spec.target.layout_mode)));
end
fprintf(fid, 'candidate_id=%d\n', cand.candidate_id);
if isfield(cand, 'layout_method'), fprintf(fid, 'layout_method=%s\n', char(string(cand.layout_method))); end
if isfield(cand, 'stage_modes'), fprintf(fid, 'stage_modes=%s\n', vec_to_inline_str_ms(cand.stage_modes)); end
if isfield(cand, 'stage_methods'), fprintf(fid, 'stage_methods=%s\n', vec_to_inline_str_ms(cand.stage_methods)); end
if isfield(cand, 'n'), fprintf(fid, 'particle_count_n=%s\n', vec_to_inline_str_ms(cand.n)); end
if isfield(cand, 'n'), fprintf(fid, 'particle_count_total=%d\n', sum(cand.n)); end
fprintf(fid, 'plate_k_inplane=%.6g\n', base_spec.plate_k_inplane);
fprintf(fid, 'ceramic_plate_k_inplane_W_mK=%.12g\n', base_spec.plate_k_inplane);
fprintf(fid, 'particle_width_mm=%.12g\n', base_spec.fp_w * 1e3);
fprintf(fid, 'particle_length_mm=%.12g\n', base_spec.fp_h * 1e3);
fprintf(fid, 'particle_height_stage_mm=%s\n', vec_to_inline_str_ms(scan.best_fp_t_stage * 1e3));
fprintf(fid, 'base_particle_height_stage_mm=%s\n', vec_to_inline_str_ms(scan.base_fp_t_stage * 1e3));
fprintf(fid, 'fp_w_mm=%.12g\n', base_spec.fp_w * 1e3);
fprintf(fid, 'fp_h_mm=%.12g\n', base_spec.fp_h * 1e3);
fprintf(fid, 'fp_t_stage_mm=%s\n', vec_to_inline_str_ms(scan.best_fp_t_stage * 1e3));
fprintf(fid, 'base_fp_t_stage_mm=%s\n', vec_to_inline_str_ms(scan.base_fp_t_stage * 1e3));
fprintf(fid, 'particle_area_mm2=%.12g\n', base_spec.fp_w * base_spec.fp_h * 1e6);
fprintf(fid, 'particle_volume_stage_mm3=%s\n', vec_to_inline_str_ms(base_spec.fp_w * base_spec.fp_h * scan.best_fp_t_stage * 1e9));
fprintf(fid, 'fixed_I_A=%.12g\n', scan.fixed_I);
fprintf(fid, 'current_I_fixed_A=%.12g\n', scan.fixed_I);
if isfield(scan, 'current_source')
    fprintf(fid, 'current_source=%s\n', char(string(scan.current_source)));
end
fprintf(fid, 'height_grid_mm=%s\n', vec_to_inline_str_ms(scan.height_grid_mm));
fprintf(fid, 'rounds=%d\n', spec.height_opt.rounds);
fprintf(fid, 'beamK=%d\n', spec.height_opt.beamK);
fprintf(fid, 'qc_abs_tol=%.12g\n', spec.height_opt.qc_abs_tol);
fprintf(fid, 'qc_rel_tol=%.12g\n', spec.height_opt.qc_rel_tol);
fprintf(fid, 'qc_effective_tol=%.12g\n', scan.qc_tol);
fprintf(fid, 'parallel_used=%d\n', logical(scan.parallel_used));
fprintf(fid, 'parallel_workers=%d\n', scan.parallel_workers);
fprintf(fid, 'dedup_enabled=%d\n', logical(scan.dedup_enabled));
fprintf(fid, 'dedup_tol_m=%.12g\n', scan.dedup_tol_m);
fprintf(fid, 'evaluated_trial_count=%d\n', scan.evaluated_trial_count);
fprintf(fid, 'duplicate_trial_count=%d\n', scan.duplicate_trial_count);
fprintf(fid, 'has_feasible_solution=%d\n', scan.has_feasible_solution);
fprintf(fid, 'best_fp_t_stage_mm=%s\n', vec_to_inline_str_ms(scan.best_fp_t_stage * 1e3));
if scan.has_feasible_solution
    b = scan.best;
    fprintf(fid, 'best_DeltaTN_actual=%.12g\n', b.DeltaTN_actual);
    fprintf(fid, 'best_Qc_last_total=%.12g\n', b.Qc_last_total);
    fprintf(fid, 'best_Qc_error=%.12g\n', b.Qc_last_total - base_spec.targets.Qc_target_last);
    fprintf(fid, 'best_TN_min=%.12g\n', b.TN_min);
    fprintf(fid, 'best_TN_mean=%.12g\n', b.TN_mean);
    fprintf(fid, 'best_message=%s\n', b.message);
else
    fprintf(fid, 'best_message=no height trial met both DeltaT and Qc equality tolerance\n');
end
end

function s = rmfield_if_exists_ms(s, fname)
if isstruct(s) && isfield(s, fname)
    s = rmfield(s, fname);
end
end

function scan = run_current_scan_for_candidate_ms(cand, G, base_spec, spec)
I_list = spec.current_opt.I_list;
N = numel(I_list);
records = repmat(empty_eval_struct_ms(base_spec.stage_count), N, 1);
rows = repmat(empty_current_scan_row_ms(), N, 1);
qc_tol = max(spec.current_opt.qc_abs_tol, spec.current_opt.qc_rel_tol * abs(base_spec.targets.Qc_target_last));
for i = 1:N
    Ii = I_list(i);
    cand_i = cand;
    cand_i.I_opt = Ii;
    G_i = G;
    G_i.I = Ii;
    fprintf('[CurrentScan] %d/%d: I=%.6g A ...\n', i, N, Ii);
    rec = evaluate_single_candidate_ms(cand_i, G_i, base_spec, true);
    records(i) = rec;
    rows(i) = build_current_scan_row_ms(i, rec, base_spec, qc_tol);
    fprintf('[CurrentScan] I=%.6g success=%d target_met=%d DeltaTN=%.6g Qc_last=%.6g msg=%s\n', ...
        Ii, rows(i).success, rows(i).target_met, rows(i).DeltaTN_actual, rows(i).Qc_last_total, rows(i).message);
end
met = [rows.target_met];
scan = struct();
scan.rows = rows;
scan.records = records;
scan.qc_tol = qc_tol;
scan.has_feasible_solution = any(met);
scan.best = empty_eval_struct_ms(base_spec.stage_count);
scan.best_I = NaN;
scan.best_DeltaTN_actual = NaN;
if scan.has_feasible_solution
    dt = [rows.DeltaTN_actual];
    dt(~met | ~isfinite(dt)) = -inf;
    [~, idx] = max(dt);
    scan.best = records(idx);
    scan.best_I = rows(idx).I_A;
    scan.best_DeltaTN_actual = rows(idx).DeltaTN_actual;
end
end

function row = empty_current_scan_row_ms()
row = struct('scan_idx', NaN, 'candidate_id', NaN, 'I_A', NaN, 'success', false, ...
    'target_met', false, 'DeltaTN_actual', NaN, 'DeltaT_target', NaN, ...
    'Qc_last_total', NaN, 'Qc_target_last', NaN, 'Qc_error', NaN, 'Qc_tol', NaN, ...
    'TN_min', NaN, 'TN_mean', NaN, 'DeltaTN_mean', NaN, 'TN_maxmin', NaN, ...
    'newton_rel_max', NaN, 'newton_iters', NaN, 'message', '');
end

function row = build_current_scan_row_ms(scan_idx, rec, spec, qc_tol)
row = empty_current_scan_row_ms();
row.scan_idx = scan_idx;
row.candidate_id = rec.candidate_id;
row.I_A = rec.I_opt;
row.success = logical(rec.success);
row.DeltaTN_actual = rec.DeltaTN_actual;
row.DeltaT_target = spec.targets.DeltaT_target;
if isfield(rec, 'Qc_last_total')
    row.Qc_last_total = rec.Qc_last_total;
end
row.Qc_target_last = spec.targets.Qc_target_last;
row.Qc_error = row.Qc_last_total - row.Qc_target_last;
row.Qc_tol = qc_tol;
row.TN_min = rec.TN_min;
row.TN_mean = rec.TN_mean;
row.DeltaTN_mean = rec.DeltaTN_mean;
row.TN_maxmin = rec.TN_maxmin;
row.newton_rel_max = rec.newton_rel_max;
row.newton_iters = rec.newton_iters;
row.message = char(string(rec.message));
row.target_met = row.success && isfinite(row.DeltaTN_actual) && row.DeltaTN_actual >= spec.targets.DeltaT_target && ...
    isfinite(row.Qc_last_total) && abs(row.Qc_error) <= qc_tol;
end

function write_current_opt_outputs_ms(scan, cand, base_spec, spec, source_name)
out_dir = spec.output.output_dir;
mkdir_if_needed_ms(out_dir);
writetable(struct_to_table_rows_ms(scan.rows), fullfile(out_dir, 'current_scan_results.csv'));
save(fullfile(out_dir, 'current_opt_result.mat'), 'scan', 'cand', 'base_spec', 'spec', 'source_name');
write_best_current_summary_ms(fullfile(out_dir, 'best_current_summary.txt'), scan, cand, base_spec, spec, source_name);
save_current_vs_deltaTN_plot_ms(scan.rows, fullfile(out_dir, 'current_vs_deltaTN.png'), base_spec);
if scan.has_feasible_solution
    best_dir = fullfile(out_dir, 'best');
    mkdir_if_needed_ms(best_dir);
    write_stage_layout_rects_csv_ms(scan.best, fullfile(best_dir, 'layout_stage_rects.csv'), base_spec.stage_count);
    caxis_by_stage = make_caxis_by_stage_ms(scan.best, base_spec.stage_count);
    save_layout_plot_ms(scan.best, fullfile(best_dir, 'layout.png'), base_spec);
    save_layout_stage_plots_ms(scan.best, best_dir, base_spec);
    save_temperature_plot_ms(scan.best, fullfile(best_dir, 'temperature.png'), caxis_by_stage, base_spec);
    save_temperature_stage_plots_ms(scan.best, best_dir, caxis_by_stage, base_spec);
else
    write_text_file_ms(fullfile(out_dir, 'best_unavailable.txt'), ...
        'No successful current point met DeltaT_target and Qc_last_total equality tolerance.');
end
end

function write_best_current_summary_ms(out_txt, scan, cand, base_spec, spec, source_name)
fid = fopen(out_txt, 'w');
if fid < 0, warning('Failed to write %s', out_txt); return; end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'fixed-layout current optimization summary\n');
fprintf(fid, 'source_results_mat=%s\n', spec.source.results_mat);
fprintf(fid, 'source_candidate_array=%s\n', char(string(source_name)));
fprintf(fid, 'candidate_id=%d\n', cand.candidate_id);
if isfield(cand, 'layout_method'), fprintf(fid, 'layout_method=%s\n', char(string(cand.layout_method))); end
if isfield(cand, 'n'), fprintf(fid, 'particle_count_n=%s\n', vec_to_inline_str_ms(cand.n)); end
if isfield(cand, 'n'), fprintf(fid, 'particle_count_total=%d\n', sum(cand.n)); end
fprintf(fid, 'plate_k_inplane=%.6g\n', base_spec.plate_k_inplane);
fprintf(fid, 'ceramic_plate_k_inplane_W_mK=%.12g\n', base_spec.plate_k_inplane);
fprintf(fid, 'particle_width_mm=%.12g\n', base_spec.fp_w * 1e3);
fprintf(fid, 'particle_length_mm=%.12g\n', base_spec.fp_h * 1e3);
fprintf(fid, 'particle_height_stage_mm=%s\n', vec_to_inline_str_ms(base_spec.fp_t_stage * 1e3));
fprintf(fid, 'fp_w_mm=%.12g\n', base_spec.fp_w * 1e3);
fprintf(fid, 'fp_h_mm=%.12g\n', base_spec.fp_h * 1e3);
fprintf(fid, 'fp_t_stage_mm=%s\n', vec_to_inline_str_ms(base_spec.fp_t_stage * 1e3));
fprintf(fid, 'particle_area_mm2=%.12g\n', base_spec.fp_w * base_spec.fp_h * 1e6);
fprintf(fid, 'particle_volume_stage_mm3=%s\n', vec_to_inline_str_ms(base_spec.fp_w * base_spec.fp_h * base_spec.fp_t_stage * 1e9));
fprintf(fid, 'current_I_init_A=%.12g\n', base_spec.current.I_init);
fprintf(fid, 'I_list=%s\n', vec_to_inline_str_ms(spec.current_opt.I_list));
fprintf(fid, 'DeltaT_target=%.12g\n', base_spec.targets.DeltaT_target);
fprintf(fid, 'Qc_target_last=%.12g\n', base_spec.targets.Qc_target_last);
fprintf(fid, 'qc_abs_tol=%.12g\n', spec.current_opt.qc_abs_tol);
fprintf(fid, 'qc_rel_tol=%.12g\n', spec.current_opt.qc_rel_tol);
fprintf(fid, 'qc_effective_tol=%.12g\n', scan.qc_tol);
fprintf(fid, 'has_feasible_solution=%d\n', scan.has_feasible_solution);
if scan.has_feasible_solution
    b = scan.best;
    fprintf(fid, 'best_I_A=%.12g\n', scan.best_I);
    fprintf(fid, 'current_I_best_A=%.12g\n', scan.best_I);
    fprintf(fid, 'best_DeltaTN_actual=%.12g\n', b.DeltaTN_actual);
    fprintf(fid, 'best_Qc_last_total=%.12g\n', b.Qc_last_total);
    fprintf(fid, 'best_Qc_error=%.12g\n', b.Qc_last_total - base_spec.targets.Qc_target_last);
    fprintf(fid, 'best_TN_min=%.12g\n', b.TN_min);
    fprintf(fid, 'best_TN_mean=%.12g\n', b.TN_mean);
    fprintf(fid, 'best_message=%s\n', b.message);
else
    fprintf(fid, 'best_message=no current point met both DeltaT and Qc equality tolerance\n');
end
end

function save_current_vs_deltaTN_plot_ms(rows, out_png, spec)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 560]);
ax = axes(fig);
hold(ax, 'on');
I = [rows.I_A];
dt = [rows.DeltaTN_actual];
ok = [rows.success] & isfinite(I) & isfinite(dt);
met = [rows.target_met] & ok;
if any(ok)
    plot(ax, I(ok), dt(ok), '-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'success');
end
if any(~ok & isfinite(I))
    plot(ax, I(~ok & isfinite(I)), zeros(1, sum(~ok & isfinite(I))), 'x', 'Color', [0.6 0.6 0.6], 'DisplayName', 'failed');
end
if any(met)
    plot(ax, I(met), dt(met), 'o', 'MarkerSize', 8, 'LineWidth', 1.8, 'DisplayName', 'target met');
    [~, local_idx] = max(dt(met));
    met_idx = find(met);
    best_idx = met_idx(local_idx);
    plot(ax, I(best_idx), dt(best_idx), 'p', 'MarkerSize', 14, 'LineWidth', 2.0, 'DisplayName', 'best');
end
yline(ax, spec.targets.DeltaT_target, '--', 'DeltaT target', 'LineWidth', 1.2);
xlabel(ax, 'Current I / A');
ylabel(ax, 'DeltaTN actual / K');
title(ax, 'Current scan: I vs DeltaTN actual');
grid(ax, 'on');
legend(ax, 'Location', 'best');
hold(ax, 'off');
exportgraphics(fig, out_png, 'Resolution', 180);
close(fig);
end

function caxis_by_stage = make_caxis_by_stage_ms(rec, N)
caxis_by_stage = NaN(N, 2);
for k = 1:N
    if isfield(rec, 'Tfields') && numel(rec.Tfields) >= k && ~isempty(rec.Tfields{k})
        Tk = rec.Tfields{k};
        if any(isfinite(Tk(:)))
            caxis_by_stage(k,:) = [min(Tk(:)), max(Tk(:))];
        end
    end
end
end

function write_text_file_ms(out_txt, txt)
fid = fopen(out_txt, 'w');
if fid < 0, warning('Failed to write %s', out_txt); return; end
fprintf(fid, '%s\n', char(string(txt)));
fclose(fid);
end

function c = normalize_trend_cell_ms(raw, N)
c = repmat({'neutral'}, 1, N);
if isempty(raw), return; end
if ischar(raw) || (isstring(raw) && isscalar(raw))
    raw = repmat({char(string(raw))}, 1, N);
elseif isstring(raw)
    raw = cellstr(raw(:).');
end
if ~iscell(raw), return; end
for i = 1:N
    idx = min(i, numel(raw));
    if idx >= 1 && ~isempty(raw{idx})
        c{i} = char(string(raw{idx}));
    end
end
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

function v = logical_field_ms(s, fname, default_v)
if nargin < 3
    default_v = false;
end
v = default_v;
if ~isstruct(s) || ~isfield(s, fname) || isempty(s.(fname))
    return;
end
try
    v = any(logical(s.(fname)(:)));
catch
    v = default_v;
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
[~, G] = optimize_layout_multistage0411_shared_params(spec, 'current');
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

function centers = rect_centers_ms(rects)
centers = [0.5*(rects(:,1)+rects(:,2)), 0.5*(rects(:,3)+rects(:,4))];
end

function calib = empty_current_calib_ms(I0)
if nargin < 1 || ~isfinite(I0)
    I0 = NaN;
end
calib = struct('enabled', false, 'has_solution', false, ...
    'I0', I0, 'I_cal', I0, 'I_shift', 0, ...
    'I_grid', zeros(0,1), 'top_ids', zeros(0,1), 'peak_I', zeros(0,1), ...
    'peak_score', zeros(0,1), 'seed_source', 'direct_full_fem', 'seed_ids', zeros(0,1));
end

function use_parallel = ensure_pool_ms(~)
use_parallel = false;
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
use_par = spec.use_parallel && N >= spec.parallel.eval_min_tasks;
if use_par
    pool = gcp('nocreate');
    if isempty(pool)
        ensure_pool_ms(spec);
        pool = gcp('nocreate');
    end
    nw = max(pool.NumWorkers, 1);
    blk = choose_parfor_block_size_ms(N, nw, spec.parallel.block_min, spec.parallel.block_max, spec.parallel.target_blocks_per_worker);
    i0 = 1;
    while i0 <= N
        i1 = min(N, i0 + blk - 1);
        chunk = cands(i0:i1);
        chunk_eval = repmat(empty_eval_struct_ms(spec.stage_count), numel(chunk), 1);
        parfor k = 1:numel(chunk)
            chunk_eval(k) = evaluate_single_candidate_ms(chunk(k), G, spec, keep_fields);
        end
        evals(i0:i1) = chunk_eval;
        i0 = i1 + 1;
    end
else
    for i = 1:N
        evals(i) = evaluate_single_candidate_ms(cands(i), G, spec, keep_fields);
    end
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
    if isfield(best, 'cache') && isfield(best.cache, 'sumQc') && numel(best.cache.sumQc) >= N
        rec.Qc_last_total = best.cache.sumQc(N);
    end
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
best_out = newton_solve_Ns_ms(G, plates, n_vec, C_init, N);
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
    'DeltaTN_actual', NaN, 'Qc_last_total', NaN, 'TN_min', NaN, 'TN_mean', NaN, 'DeltaTN_mean', NaN, 'TN_maxmin', NaN, ...
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
if nargin < 1 || ~isstruct(spec) || ~isfield(spec, 'stage_count')
    N = 5;
else
    N = spec.stage_count;
end
prior = struct('enable', false, 'stage_trend_prefer', {normalize_trend_cell_ms({}, N)});
if ~isstruct(spec) || ~isfield(spec, 'soft_prior') || ~isstruct(spec.soft_prior)
    return;
end
sp = spec.soft_prior(1);
if isfield(sp, 'enable') && ~isempty(sp.enable)
    prior.enable = any(logical(sp.enable(:)));
end
if isfield(sp, 'stage_trend_prefer')
    prior.stage_trend_prefer = normalize_trend_cell_ms(sp.stage_trend_prefer, N);
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

%% ====================== Top post + uniform baseline ======================

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
    k_interfaces = zp.k_interfaces;
elseif isfield(zp, 'k_interface')
    k_interfaces = zp.k_interface;
else
    k_interfaces = zcfg.k_interfaces;
end
if isfield(zp, 't_interface_effs') && ~isempty(zp.t_interface_effs)
    t_interfaces = zp.t_interface_effs;
elseif isfield(zp, 't_interface_eff')
    t_interfaces = zp.t_interface_eff;
else
    t_interfaces = zcfg.t_interface_effs;
end
if isfield(zp, 'Rc_interfaces') && ~isempty(zp.Rc_interfaces)
    Rc_interfaces = zp.Rc_interfaces;
elseif isfield(zp, 'Rc_interface')
    Rc_interfaces = zp.Rc_interface;
else
    Rc_interfaces = zcfg.Rc_interfaces;
end
zcfg.k_interfaces = max(fit_len_vec_ms(k_interfaces, nif, zcfg.k_interfaces), 1e-9);
zcfg.t_interface_effs = max(fit_len_vec_ms(t_interfaces, nif, zcfg.t_interface_effs), 0);
zcfg.Rc_interfaces = max(fit_len_vec_ms(Rc_interfaces, nif, zcfg.Rc_interfaces), 0);
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
