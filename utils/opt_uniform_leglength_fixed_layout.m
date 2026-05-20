﻿function results = opt_uniform_leglength_fixed_layout(spec_in)
% Scan uniform thermoelectric leg length for one fixed candidate or regular-grid layout.
% This is the uniform mode: every stage uses the same leg length for each scan point.
% Internally each point calls opt_stage_leglength_fixed_layout with one height value.
%
% Candidate uniform leg length scan:
%   spec = struct();
%   spec.source.results_mat = 'D:\...\results_data.mat';
%   spec.target.layout_mode = 'candidate';
%   spec.target.candidate_id = 23;
%   spec.leg_opt.leg_length_mm = 0.8:0.1:2.4;
%   results_L = opt_uniform_leglength_fixed_layout(spec);
%
% Regular-grid uniform leg length scan:
%   spec = struct();
%   spec.target.layout_mode = 'regular_grid';
%   spec.stage_count = 5;
%   spec.fixed_n = [322 110 42 14 6];
%   spec.current.I_init = 3.0;
%   spec.leg_opt.leg_length_mm = 0.8:0.1:2.4;
%   results_L = opt_uniform_leglength_fixed_layout(spec);
%
% Continue from optimized current:
%   spec.current.I_fixed = results_I.best_I;
%   results_L = opt_uniform_leglength_fixed_layout(spec);
%
% Or load optimized current from result MAT:
%   spec.current.source_result_mat = 'D:\...\current_opt_result.mat';
%   results_L = opt_uniform_leglength_fixed_layout(spec);
if nargin < 1 || isempty(spec_in)
    spec_in = struct();
end
spec = apply_leg_opt_defaults_ms(spec_in);
mkdir_if_needed_ms(spec.output.output_dir);
if is_regular_grid_mode_ms(spec)
    spec = prepare_regular_grid_source_ms(spec);
end

fprintf('\n========== opt_uniform_leglength_fixed_layout ==========\n');
fprintf('Source MAT: %s\n', spec.source.results_mat);
fprintf('Layout mode: %s\n', spec.target.layout_mode);
fprintf('Target candidate_id=%d\n', spec.target.candidate_id);
fprintf('leg_length_mm=%s\n', vec_to_inline_str_ms(spec.leg_opt.leg_length_mm));
fprintf('Output: %s\n', spec.output.output_dir);

L_mm = spec.leg_opt.leg_length_mm(:).';
rows = repmat(empty_leg_scan_row_ms(), numel(L_mm), 1);
height_results = cell(numel(L_mm), 1);
spec = ensure_leg_scan_parallel_pool_ms(spec, numel(L_mm));
use_par = logical(spec.use_parallel) && numel(L_mm) > 1 && has_active_parallel_pool_ms();
fprintf('Parallel: requested=%d, used=%d, workers=%d\n', ...
    logical(spec.parallel.requested), use_par, get_parallel_worker_count_ms());

if use_par
    parfor i = 1:numel(L_mm)
        [height_results{i}, rows(i)] = evaluate_leg_length_point_ms(i, L_mm(i), L_mm, spec, true);
    end
else
    for i = 1:numel(L_mm)
        [height_results{i}, rows(i)] = evaluate_leg_length_point_ms(i, L_mm(i), L_mm, spec, false);
    end
end

[best_row, best_idx] = pick_best_leg_scan_row_ms(rows);
results = struct();
results.spec = spec;
results.rows = rows;
results.height_results = height_results;
results.best_row = best_row;
results.best_idx = best_idx;
results.has_feasible_solution = any([rows.target_met]);
results.output_dir = spec.output.output_dir;

writetable(struct2table(rows), fullfile(spec.output.output_dir, 'leg_length_scan_results.csv'));
save(fullfile(spec.output.output_dir, 'leg_length_scan_result.mat'), ...
    'results', 'rows', 'height_results', 'spec', '-v7.3');
write_leg_scan_summary_ms(fullfile(spec.output.output_dir, 'best_leg_length_summary.txt'), results);
save_leg_length_vs_deltaTN_plot_ms(rows, fullfile(spec.output.output_dir, 'leg_length_vs_DeltaTN.png'), spec);

fprintf('Done. has_feasible_solution=%d, best_leg_length_mm=%.6g, best_DeltaTN=%.6g K\n', ...
    results.has_feasible_solution, best_row.leg_length_mm, best_row.DeltaTN_actual);
fprintf('===========================================================\n\n');
end

function spec = apply_leg_opt_defaults_ms(spec)
if ~isstruct(spec), error('spec_in must be a struct.'); end
if ~isfield(spec, 'source') || ~isstruct(spec.source), spec.source = struct(); end
if ~isfield(spec, 'target') || ~isstruct(spec.target), spec.target = struct(); end
if ~isfield(spec, 'leg_opt') || ~isstruct(spec.leg_opt), spec.leg_opt = struct(); end
if ~isfield(spec, 'output') || ~isstruct(spec.output), spec.output = struct(); end
if ~isfield(spec, 'parallel') || ~isstruct(spec.parallel), spec.parallel = struct(); end
if ~isfield(spec, 'use_parallel') || isempty(spec.use_parallel)
    spec.use_parallel = true;
end
spec.use_parallel = any(logical(spec.use_parallel(:)));
if ~isfield(spec.parallel, 'pool_workers') || isempty(spec.parallel.pool_workers) || ~isfinite(spec.parallel.pool_workers)
    spec.parallel.pool_workers = 64;
end
spec.parallel.pool_workers = max(1, round(double(spec.parallel.pool_workers)));
spec.parallel.requested = spec.use_parallel;
if ~isfield(spec.target, 'layout_mode') || isempty(spec.target.layout_mode)
    spec.target.layout_mode = 'candidate';
end
spec.target.layout_mode = lower(strtrim(char(string(spec.target.layout_mode))));
if ~any(strcmp(spec.target.layout_mode, {'candidate', 'regular_grid'}))
    error('spec.target.layout_mode must be ''candidate'' or ''regular_grid''.');
end
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
if ~isfield(spec.leg_opt, 'leg_length_mm') || isempty(spec.leg_opt.leg_length_mm)
    spec.leg_opt.leg_length_mm = 0.8:0.1:2.4;
end
L = unique(double(spec.leg_opt.leg_length_mm(:).'));
L = L(isfinite(L) & L > 0);
if isempty(L)
    error('spec.leg_opt.leg_length_mm must contain at least one positive finite value.');
end
spec.leg_opt.leg_length_mm = L;
if ~isfield(spec.output, 'output_dir') || isempty(spec.output.output_dir)
    if is_regular_grid_mode_ms(spec)
        base_dir = fullfile(pwd, 'leg_length_scan_regular_grid');
        target_tag = 'regular_grid';
    else
        src_dir = fileparts(spec.source.results_mat);
        base_dir = fullfile(src_dir, sprintf('leg_length_scan_candidate%d', spec.target.candidate_id));
        target_tag = sprintf('candidate%d', spec.target.candidate_id);
    end
    stamp = datestr(now, 'yyyymmdd_HHMMSS');
    run_name = sprintf('leg_length_scan_%s_L%s_%s', target_tag, numeric_vec_tag_ms(L), stamp);
    spec.output.output_dir = fullfile(base_dir, run_name);
end
spec.output.output_dir = char(string(spec.output.output_dir));
end

function tf = is_regular_grid_mode_ms(spec)
tf = isfield(spec, 'target') && isstruct(spec.target) && ...
    isfield(spec.target, 'layout_mode') && strcmpi(char(string(spec.target.layout_mode)), 'regular_grid');
end

function spec = ensure_leg_scan_parallel_pool_ms(spec, task_count)
if nargin < 2 || isempty(task_count) || ~isfinite(task_count)
    task_count = 0;
end
spec.parallel.pool_used = false;
spec.parallel.pool_workers_actual = 0;
if ~logical(spec.use_parallel) || task_count <= 1
    spec.use_parallel = false;
    return;
end
has_parallel = ((exist('parpool', 'file') == 2) || (exist('parpool', 'builtin') == 5)) && ...
    ((exist('gcp', 'file') == 2) || (exist('gcp', 'builtin') == 5));
if ~has_parallel
    warning('Leg-length scan parallel requested but Parallel Computing Toolbox functions are unavailable; falling back to serial.');
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
        fprintf('[LegScan] existing parallel pool workers=%d, requested=%d; using existing pool.\n', ...
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
    warning('Leg-length scan parallel pool unavailable: %s. Falling back to serial.', ME_pool.message);
    spec.use_parallel = false;
end
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

function [ri, row] = evaluate_leg_length_point_ms(i, Li, L_mm, spec, outer_parallel)
run_spec = spec;
run_spec.height_opt = build_height_spec_for_leg_length_ms(spec, Li);
run_spec.output.output_dir = fullfile(spec.output.output_dir, sprintf('leg_%smm', numeric_tag_ms(Li)));
if outer_parallel
    run_spec.use_parallel = false;
end
fprintf('[LegScan] %d/%d: L=%.6g mm ...\n', i, numel(L_mm), Li);
ri = opt_stage_leglength_fixed_layout(run_spec);
row = build_leg_scan_row_ms(i, Li, ri, spec);
fprintf('[LegScan] L=%.6g mm success=%d target_met=%d DeltaTN=%.6g Qc_last=%.6g\n', ...
    Li, row.success, row.target_met, row.DeltaTN_actual, row.Qc_last_total);
end

function height_opt = build_height_spec_for_leg_length_ms(spec, L_mm)
height_opt = struct();
height_opt.height_grid_mm = L_mm;
height_opt.rounds = 1;
height_opt.beamK = 1;
height_opt.dedup_enable = true;
height_opt.dedup_tol_m = 1e-12;
if isfield(spec.leg_opt, 'qc_abs_tol')
    height_opt.qc_abs_tol = spec.leg_opt.qc_abs_tol;
end
if isfield(spec.leg_opt, 'qc_rel_tol')
    height_opt.qc_rel_tol = spec.leg_opt.qc_rel_tol;
end
end

function spec = prepare_regular_grid_source_ms(spec)
[base_spec, ~] = optimize_layout_multistage0411_shared_params(spec, 'layout');
if isfield(spec, 'current') && isstruct(spec.current)
    if ~isfield(base_spec, 'current') || ~isstruct(base_spec.current)
        base_spec.current = struct();
    end
    if isfield(spec.current, 'I_init') && ~isempty(spec.current.I_init)
        base_spec.current.I_init = spec.current.I_init;
    end
    if isfield(spec.current, 'I_fixed') && ~isempty(spec.current.I_fixed)
        base_spec.current.I_fixed = spec.current.I_fixed;
    end
    if isfield(spec.current, 'source_result_mat') && ~isempty(spec.current.source_result_mat)
        base_spec.current.source_result_mat = spec.current.source_result_mat;
    end
end
base_spec.target = spec.target;
base_spec.leg_opt = spec.leg_opt;
base_spec.output = spec.output;
best_cand = build_regular_grid_candidate_from_spec_ms(base_spec);
source_dir = fullfile(spec.output.output_dir, '_regular_grid_source');
mkdir_if_needed_ms(source_dir);
results_mat = fullfile(source_dir, 'results_data.mat');
spec = base_spec;
spec.source.results_mat = results_mat;
best_cand.I_opt = spec.current.I_init;
save(results_mat, 'spec', 'best_cand', '-v7.3');
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

function row = empty_leg_scan_row_ms()
row = struct('scan_idx', NaN, 'leg_length_mm', NaN, 'candidate_id', NaN, ...
    'I_A', NaN, 'current_source', '', 'success', false, 'target_met', false, ...
    'DeltaTN_actual', NaN, 'DeltaT_target', NaN, ...
    'Qc_last_total', NaN, 'Qc_target_last', NaN, 'Qc_error', NaN, ...
    'TN_min', NaN, 'TN_mean', NaN, 'DeltaTN_mean', NaN, ...
    'TN_maxmin', NaN, 'newton_rel_max', NaN, 'newton_iters', NaN, ...
    'source_output_dir', '', 'message', '');
end

function row = build_leg_scan_row_ms(scan_idx, L_mm, ri, spec)
row = empty_leg_scan_row_ms();
row.scan_idx = scan_idx;
row.leg_length_mm = L_mm;
row.source_output_dir = ri.output_dir;
if isfield(ri, 'fixed_I')
    row.I_A = ri.fixed_I;
end
if isfield(ri, 'current_source')
    row.current_source = char(string(ri.current_source));
end
rows = ri.scan_rows;
idx = find_uniform_leg_row_ms(rows, L_mm);
if isempty(idx)
    row.message = 'uniform leg-length row not found in height scan output';
    return;
end
r = rows(idx);
row.candidate_id = r.candidate_id;
if isfinite(r.I_A)
    row.I_A = r.I_A;
end
row.success = logical(r.success);
row.target_met = logical(r.target_met);
row.DeltaTN_actual = r.DeltaTN_actual;
row.DeltaT_target = r.DeltaT_target;
row.Qc_last_total = r.Qc_last_total;
row.Qc_target_last = r.Qc_target_last;
row.Qc_error = r.Qc_error;
row.TN_min = r.TN_min;
row.TN_mean = r.TN_mean;
row.DeltaTN_mean = r.DeltaTN_mean;
row.TN_maxmin = r.TN_maxmin;
row.newton_rel_max = r.newton_rel_max;
row.newton_iters = r.newton_iters;
row.message = char(string(r.message));
if ~row.target_met && isfield(ri, 'best') && isstruct(ri.best) && isfield(ri.best, 'message') && ~isempty(ri.best.message)
    row.message = char(string(ri.best.message));
end
if isfield(spec.leg_opt, 'accept_success_only') && logical(spec.leg_opt.accept_success_only)
    row.target_met = row.success;
end
end

function idx = find_uniform_leg_row_ms(rows, L_mm)
idx = [];
if isempty(rows)
    return;
end
target_m = L_mm * 1e-3;
tol_m = 1e-10;
for i = 1:numel(rows)
    if isfield(rows(i), 'fp_t_stage_numeric')
        v = rows(i).fp_t_stage_numeric;
        if ~isempty(v) && all(isfinite(v)) && max(abs(v(:).' - target_m)) <= tol_m
            idx = i;
            return;
        end
    end
end
idx = numel(rows);
end

function [best_row, best_idx] = pick_best_leg_scan_row_ms(rows)
best_row = empty_leg_scan_row_ms();
best_idx = NaN;
if isempty(rows)
    return;
end
met = [rows.target_met].';
succ = [rows.success].';
dt = [rows.DeltaTN_actual].';
dt(~isfinite(dt)) = -inf;
if any(met)
    score = dt;
    score(~met) = -inf;
elseif any(succ)
    score = dt;
    score(~succ) = -inf;
else
    score = dt;
end
[~, best_idx] = max(score);
best_row = rows(best_idx);
end

function save_leg_length_vs_deltaTN_plot_ms(rows, out_png, spec)
L = [rows.leg_length_mm];
dt = [rows.DeltaTN_actual];
ok = [rows.success] & isfinite(dt);
met = [rows.target_met] & isfinite(dt);
[L, ord] = sort(L);
dt = dt(ord);
ok = ok(ord);
met = met(ord);
fig = figure('Visible', 'off', 'Color', 'w');
cleanup = onCleanup(@() close(fig));
ax = axes(fig);
hold(ax, 'on');
grid(ax, 'on');
if any(ok)
    plot(ax, L(ok), dt(ok), '-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'success');
end
if any(~ok)
    plot(ax, L(~ok), zeros(1, sum(~ok)), 'x', 'Color', [0.6 0.6 0.6], 'DisplayName', 'failed');
end
if any(met)
    plot(ax, L(met), dt(met), 'o', 'MarkerSize', 8, 'LineWidth', 1.8, 'DisplayName', 'target met');
end
if isfield(spec, 'targets') && isstruct(spec.targets) && isfield(spec.targets, 'DeltaT_target')
    yline(ax, spec.targets.DeltaT_target, '--', 'DeltaT target', 'LineWidth', 1.2);
end
xlabel(ax, 'Leg length / mm');
ylabel(ax, 'DeltaTN actual / K');
title(ax, 'Uniform leg length scan: leg length vs DeltaTN actual');
legend(ax, 'Location', 'best');
exportgraphics(fig, out_png, 'Resolution', 200);
end

function write_leg_scan_summary_ms(out_txt, results)
fid = fopen(out_txt, 'w');
if fid < 0
    warning('Failed to write %s', out_txt);
    return;
end
cleanup = onCleanup(@() fclose(fid));
b = results.best_row;
fprintf(fid, 'fixed-layout uniform leg-length scan summary\n');
fprintf(fid, 'source_results_mat=%s\n', results.spec.source.results_mat);
if isfield(results.spec, 'target') && isstruct(results.spec.target) && isfield(results.spec.target, 'layout_mode')
    fprintf(fid, 'target_layout_mode=%s\n', char(string(results.spec.target.layout_mode)));
end
fprintf(fid, 'candidate_id=%d\n', results.spec.target.candidate_id);
fprintf(fid, 'leg_length_grid_mm=%s\n', vec_to_inline_str_ms(results.spec.leg_opt.leg_length_mm));
fprintf(fid, 'has_feasible_solution=%d\n', results.has_feasible_solution);
fprintf(fid, 'best_leg_length_mm=%.12g\n', b.leg_length_mm);
fprintf(fid, 'best_DeltaTN_actual=%.12g\n', b.DeltaTN_actual);
fprintf(fid, 'best_I_A=%.12g\n', b.I_A);
fprintf(fid, 'fixed_I_A=%.12g\n', b.I_A);
fprintf(fid, 'current_source=%s\n', char(string(b.current_source)));
fprintf(fid, 'best_Qc_last_total=%.12g\n', b.Qc_last_total);
fprintf(fid, 'best_Qc_error=%.12g\n', b.Qc_error);
fprintf(fid, 'best_success=%d\n', logical(b.success));
fprintf(fid, 'best_target_met=%d\n', logical(b.target_met));
fprintf(fid, 'best_message=%s\n', b.message);
end

function mkdir_if_needed_ms(path_str)
if exist(path_str, 'dir') ~= 7
    mkdir(path_str);
end
end

function s = vec_to_inline_str_ms(v)
if iscell(v)
    parts = cellfun(@(x) char(string(x)), v, 'UniformOutput', false);
    s = ['[', strjoin(parts, ' '), ']'];
    return;
end
v = double(v(:).');
if isempty(v)
    s = '[]';
else
    s = ['[', strjoin(arrayfun(@(x) sprintf('%.9g', x), v, 'UniformOutput', false), ' '), ']'];
end
end

function s = numeric_vec_tag_ms(v)
parts = arrayfun(@numeric_tag_ms, double(v(:).'), 'UniformOutput', false);
s = strjoin(parts, '_');
end

function s = numeric_tag_ms(v)
if ~isfinite(v)
    s = 'nan';
    return;
end
s = regexprep(sprintf('%.6g', v), '[^\dA-Za-z]+', 'p');
s = regexprep(s, '^p+', '');
if isempty(s)
    s = '0';
end
end
