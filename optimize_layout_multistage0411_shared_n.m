function results = optimize_layout_multistage0411_shared_n(spec_in)
% Modified 2026-05-12 20:51:42 +08:00:
%   Allow intermediate-stage temperature reversal under high heat load, complete FEM Qc target checks,
%   and add lightweight DeltaT_at_Qc metric reporting.
% Generic multistage optimizer (default N=5) for pyramid TEC layout.
% New API (non-compatible with legacy three-stage fields):
%   spec.stage_count
%   spec.targets.Qc_target_last
%   spec.geometry.L_max_mm
%   spec.geometry.coverage_min
%   spec.geometry.pyramid_gap_min_mm
%   spec.geometry.min_edge_gap_mm
if nargin < 1
    [spec, spec_source] = resolve_runtime_spec_ms();
else
    [spec, spec_source] = resolve_runtime_spec_ms(spec_in);
end
spec = apply_default_spec_ms(spec);
spec = apply_metric_mode_targets_ms(spec);
spec = prepare_run_output_dir_ms(spec);
spec.k25_debug.enable = true;
rng(spec.seed, 'twister');

script_folder = fileparts(mfilename('fullpath'));
if ~contains(path, script_folder)
    addpath(script_folder);
end
mkdir_if_needed_ms(spec.output.output_dir);
spec = ensure_k25_parallel_pool_ms(spec);

fprintf('\n========== optimize_layout_multistage ==========\n');
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
G_template = apply_low_k_solver_stabilization_ms(G_template, spec);

[ok_run, results_run, fail_reason] = run_k25_single_candidate_debug_ms(spec, spec_source, G_template, t0);
if ~ok_run
    error('%s', fail_reason);
end
results = results_run;
end

function [ok, results, fail_reason] = run_k25_single_candidate_debug_ms(spec, spec_source, G_template, t0)
ok = false;
results = struct();
fail_reason = '';
try
    G = G_template;
    G.I = spec.current.I_init;

    fprintf('[K25-Debug] particle-count mode enabled.\n');
    if has_explicit_k25_n_ms(spec)
        fprintf('[K25-Debug] fixed n mode: fixed I=%.6g A; n-current scan is not used.\n', G.I);
    elseif is_n_current_scan_enabled_ms(spec)
        fprintf('[K25-Debug] 0D n-current scan enabled: I_list=%s; selected n/I is based on 0D target checks.\n', ...
            vec_to_inline_str_ms(get_n_current_scan_I_list_ms(spec)));
    else
        fprintf('[K25-Debug] n-current scan disabled: fixed I=%.6g A.\n', G.I);
    end
    if spec.k25_debug.search_only
        fprintf('[K25-Debug] search_only=1: final FEM/layout temperature solve will be skipped.\n');
    end

    t_eval = tic;
    particle_search_rows = repmat(empty_particle_search_row_ms(), 0, 1);
    top_n_current = repmat(empty_top_n_current_row_ms(spec.stage_count), 0, 1);
    ratio_ref_rec = empty_eval_struct_ms(spec.stage_count);
    ratio_ref_cand = candidate_record_template_ms(spec.stage_count);
    ratio_ref_row = empty_particle_search_row_ms();
    if has_explicit_k25_n_ms(spec)
        n_vec = normalize_even_count_vector_ms(spec.k25_debug.n, spec.stage_count);
        fprintf('[K25-Debug] fixed n=%s, I_init=%.3f A, k=%.3g W/mK\n', ...
            vec_to_inline_str_ms(n_vec), G.I, spec.plate_k_inplane);
        [rec, cand] = evaluate_k25_n_vector_ms(n_vec, G, spec, true);
        particle_search_rows = build_particle_search_row_ms(rec, n_vec, 'fixed', 1, spec, toc(t_eval));
        top_n_current = build_top_n_current_rows_ms(particle_search_rows, spec);
    else
        fprintf('[K25-Debug] particle-count search enabled, I_init=%.3f A, k=%.3g W/mK\n', ...
            G.I, spec.plate_k_inplane);
        [rec, cand, n_vec, particle_search_rows, ratio_ref_rec, ratio_ref_cand, ratio_ref_row, top_n_current] = ...
            run_particle_count_search_ms(G, spec);
    end
    total_n = sum(n_vec);
    scan_rows = build_particle_search_row_ms(rec, n_vec, 'selected', 1, spec, NaN);
    eval_sec = toc(t_eval);
    row = build_k25_single_result_row_ms(rec, spec, total_n, eval_sec);
    if spec.k25_debug.search_only && ~has_explicit_k25_n_ms(spec)
        fprintf('[K25-Debug] 0D search done: success=%d, 0D_target_met=%d, I_opt_0d=%.4f A, DeltaT_0d=%.4f K, Qc_0d=%.6g W, elapsed=%.2f s, msg=%s\n', ...
            row.success, row.target_met, row.I_A, row.DeltaTN_actual, row.Qc_last_total, row.elapsed_sec, row.message);
    elseif ~has_explicit_k25_n_ms(spec)
        selected_0d_target_met = false;
        selected_0d_I = NaN;
        selected_0d_dt = NaN;
        selected_0d_qc = NaN;
        if isstruct(rec)
            if isfield(rec, 'selected_0d_target_met'), selected_0d_target_met = logical(rec.selected_0d_target_met); end
            if isfield(rec, 'selected_0d_I_opt'), selected_0d_I = rec.selected_0d_I_opt; end
            if isfield(rec, 'selected_0d_DeltaTN_actual'), selected_0d_dt = rec.selected_0d_DeltaTN_actual; end
            if isfield(rec, 'selected_0d_Qc_last_total'), selected_0d_qc = rec.selected_0d_Qc_last_total; end
        end
        fprintf('[K25-Debug] 0D selection: 0D_target_met=%d, I_opt_0d=%.4f A, DeltaT_0d=%.4f K, Qc_0d=%.6g W.\n', ...
            selected_0d_target_met, selected_0d_I, selected_0d_dt, selected_0d_qc);
        fprintf('[K25-Debug] FEM diagnostic done: success=%d, FEM_target_met_diagnostic=%d, I_FEM=%.4f A, DeltaT_FEM=%.4f K, elapsed=%.2f s, msg=%s\n', ...
            row.success, row.target_met, row.I_A, row.DeltaTN_actual, row.elapsed_sec, row.message);
    else
        fprintf('[K25-Debug] FEM done: success=%d, target_met=%d, I=%.4f A, DeltaTN_actual=%.4f, elapsed=%.2f s, msg=%s\n', ...
            row.success, row.target_met, row.I_A, row.DeltaTN_actual, row.elapsed_sec, row.message);
    end

    out_dir = spec.output.output_dir;
    mkdir_if_needed_ms(out_dir);
    write_scalar_struct_csv_ms(row, fullfile(out_dir, 'k25_single_layout_result.csv'));
    write_particle_search_results_csv_ms(particle_search_rows, fullfile(out_dir, 'particle_search_results.csv'));
    write_top_n_current_csv_ms(top_n_current, fullfile(out_dir, 'top_n_current.csv'));
    write_particle_search_results_csv_ms(build_particle_search_row_ms(rec, n_vec, 'best', 1, spec, eval_sec), ...
        fullfile(out_dir, 'best_particle_solution.csv'));
    if isstruct(ratio_ref_row) && isfinite(ratio_ref_row.total_n)
        write_particle_search_results_csv_ms(ratio_ref_row, fullfile(out_dir, 'ratio_reference_solution.csv'));
        write_best_vs_ratio_reference_csv_ms(rec, ratio_ref_rec, n_vec, ratio_ref_row.n_vec_numeric, ...
            fullfile(out_dir, 'best_vs_ratio_reference.csv'), spec);
    end
    write_stage_layout_rects_csv_ms(rec, fullfile(out_dir, 'layout_stage_rects.csv'), spec.stage_count);
    try
        if spec.k25_debug.search_only && ~has_explicit_k25_n_ms(spec)
            write_text_file_ms(fullfile(out_dir, 'fem_skipped_search_only.txt'), ...
                sprintf('search_only=1: 0D particle-count search completed; final FEM/layout temperature solve was skipped.\nn=%s\nDeltaTN_0d=%.9g\n', ...
                vec_to_inline_str_ms(n_vec), rec.DeltaTN_actual));
        else
	            save_layout_plot_ms(rec, fullfile(out_dir, 'layout.png'), spec);
	            save_layout_stage_plots_ms(rec, out_dir, spec);
	            if has_temperature_fields_ms(rec, spec.stage_count)
	                caxis_by_stage = make_caxis_by_stage_from_eval_ms(rec, spec.stage_count);
	                save_temperature_plot_ms(rec, fullfile(out_dir, 'temperature.png'), caxis_by_stage, spec);
	                save_temperature_stage_plots_ms(rec, out_dir, caxis_by_stage, spec);
	                if ~rec.success
	                    fail_txt = sprintf(['FEM temperature fields were exported for diagnostics, but the result was rejected by validation.\n' ...
	                        'reason=%s\n'], rec.message);
	                    if isfield(rec, 'newton_failure_detail') && ~isempty(rec.newton_failure_detail)
	                        fail_txt = sprintf('%s%s\n', fail_txt, rec.newton_failure_detail);
	                    end
	                    write_text_file_ms(fullfile(out_dir, 'temperature_diagnostic_note.txt'), fail_txt);
	                end
	            else
	                fail_txt = sprintf('FEM failed or was rejected by validation; accepted temperature fields unavailable.\nreason=%s\n', rec.message);
	                write_text_file_ms(fullfile(out_dir, 'temperature_unavailable.txt'), ...
	                    fail_txt);
	            end
	            if isstruct(ratio_ref_row) && isfinite(ratio_ref_row.total_n)
	                ratio_dir = fullfile(out_dir, 'ratio_reference');
	                mkdir_if_needed_ms(ratio_dir);
	                save_layout_plot_ms(ratio_ref_rec, fullfile(ratio_dir, 'layout.png'), spec);
	                save_layout_stage_plots_ms(ratio_ref_rec, ratio_dir, spec);
	                if has_temperature_fields_ms(ratio_ref_rec, spec.stage_count)
	                    caxis_ratio = make_caxis_by_stage_from_eval_ms(ratio_ref_rec, spec.stage_count);
	                    save_temperature_plot_ms(ratio_ref_rec, fullfile(ratio_dir, 'temperature.png'), caxis_ratio, spec);
	                    save_temperature_stage_plots_ms(ratio_ref_rec, ratio_dir, caxis_ratio, spec);
	                end
	            end
        end
    catch ME_plot
        warning('K25 debug plot export failed: %s', ME_plot.message);
    end
    metric_results = compute_metric_postprocess_ms(rec, cand, ratio_ref_rec, ratio_ref_cand, G, spec);
    write_metric_outputs_ms(metric_results, out_dir);
    write_k25_fast_summary_ms(fullfile(out_dir, 'k25_debug_summary.txt'), spec, spec_source, ...
        cand, rec, row, t0, metric_results);
    write_k25_debug_result_mat_ms(fullfile(out_dir, 'k25_debug_result.mat'), ...
        spec, cand, rec, row, scan_rows, particle_search_rows, ratio_ref_rec, ratio_ref_cand, top_n_current, metric_results);

    results = struct();
    results.spec = spec;
    results.spec_source = spec_source;
    results.count_solution = struct('n', n_vec, 'I_opt', rec.I_opt);
    results.best_n = n_vec;
    results.best_I = rec.I_opt;
    if ~isempty(top_n_current)
        results.best_n = top_n_current(1).n;
        results.best_I = top_n_current(1).I_opt;
        results.count_solution = struct('n', results.best_n, 'I_opt', results.best_I);
    end
    results.top_n_current = top_n_current;
    results.debug_candidate = cand;
    results.k25_result = rec;
    results.k25_result_row = row;
    results.current_scan_rows = scan_rows;
    results.particle_search_rows = particle_search_rows;
    results.ratio_reference_result = ratio_ref_rec;
    results.ratio_reference_candidate = ratio_ref_cand;
    results.metric_results = metric_results;
    results.output_dir = out_dir;
    results.runtime_sec = toc(t0);
    fprintf('[K25-Debug] wrote single-layout outputs to %s\n', out_dir);
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

function n_vec = allocate_k25_standard_counts_ms(total_n, N)
ratio_cold_to_hot = [1, 3, 7, 21, 55];
ratio_cold_to_hot = fit_len_vec_ms(ratio_cold_to_hot, N, ratio_cold_to_hot);
ratio_hot_to_cold = ratio_cold_to_hot(end:-1:1);
ratio_hot_to_cold = max(ratio_hot_to_cold(:).', eps);
raw = total_n * ratio_hot_to_cold / sum(ratio_hot_to_cold);
n_vec = max(2, round(raw));
n_vec = rebalance_counts_to_total_ms(n_vec, raw, total_n, true);
end

function tf = has_explicit_k25_n_ms(spec)
tf = isstruct(spec) && isfield(spec, 'k25_debug') && isstruct(spec.k25_debug) && ...
    isfield(spec.k25_debug, 'n') && ~isempty(spec.k25_debug.n);
end

function n_vec = normalize_even_count_vector_ms(n_in, N)
n_vec = fit_len_vec_ms(double(n_in), N, 2 * ones(1, N));
n_vec = max(2, 2 * round(n_vec / 2));
end

function [rec, cand] = evaluate_k25_n_vector_ms(n_vec, G, spec, keep_fields)
if nargin < 4
    keep_fields = false;
end
spec_try = spec;
spec_try.k25_debug.total_n = sum(n_vec);
[ok_cand, cand, cand_msg] = build_k25_standard_candidate_ms(spec_try, n_vec, G.I);
if ~ok_cand
    rec = empty_eval_struct_ms(spec.stage_count);
    rec.n = n_vec(:).';
    rec.message = sprintf('particle-count candidate build failed: %s', cand_msg);
    return;
end
[ok_geom, geom_msg] = validate_k25_count_geometry_ms(cand, spec_try);
if ~ok_geom
    rec = empty_eval_struct_ms(spec.stage_count);
    rec.n = n_vec(:).';
    rec.Lx = cand.Lx;
    rec.Ly = cand.Ly;
    rec.cov = cand.cov;
    rec.stage_rects = cand.stage_rects;
    rec.message = geom_msg;
    return;
end
rec = evaluate_single_candidate_ms(cand, G, spec_try, keep_fields);
end

function [rec, cand] = evaluate_k25_n_vector_0d_ms(n_vec, G, spec)
if is_n_current_scan_enabled_ms(spec)
    [rec, cand] = evaluate_k25_n_vector_0d_current_scan_ms(n_vec, G, spec);
else
    [rec, cand] = evaluate_k25_n_vector_0d_at_current_ms(n_vec, G, spec, G.I);
end
end

function [best_rec, best_cand] = evaluate_k25_n_vector_0d_current_scan_ms(n_vec, G, spec)
N = spec.stage_count;
I_list = get_n_current_scan_I_list_ms(spec);
qc_tol = get_n_current_scan_qc_tol_ms(spec);
best_rec = empty_eval_struct_ms(N);
best_cand = candidate_record_template_ms(N);
best_found = false;
for ii = 1:numel(I_list)
    [rec_i, cand_i] = evaluate_k25_n_vector_0d_at_current_ms(n_vec, G, spec, I_list(ii));
    rec_i.Qc_tol = qc_tol;
    rec_i.Qc_target_last = spec.targets.Qc_target_last;
    if isfield(rec_i, 'sumQc_stage') && numel(rec_i.sumQc_stage) >= N
        rec_i.Qc_last_total = rec_i.sumQc_stage(N);
    end
    if isfield(rec_i, 'Qc_last_total') && isfinite(rec_i.Qc_last_total)
        rec_i.Qc_error = rec_i.Qc_last_total - spec.targets.Qc_target_last;
    else
        rec_i.Qc_error = NaN;
    end
    rec_i.target_met = logical(rec_i.success) && isfinite(rec_i.DeltaTN_actual) && ...
        rec_i.DeltaTN_actual >= spec.targets.DeltaT_target && ...
        isfinite(rec_i.Qc_error) && abs(rec_i.Qc_error) <= qc_tol;
    if rec_i.target_met && (~best_found || is_better_n_current_scan_rec_ms(rec_i, best_rec))
        best_rec = rec_i;
        best_cand = cand_i;
        best_found = true;
    elseif ~best_found && rec_i.success && (~best_rec.success || is_better_n_current_scan_rec_ms(rec_i, best_rec))
        best_rec = rec_i;
        best_cand = cand_i;
    end
end
if best_found
    best_rec.message = '0D current scan ok';
elseif isstruct(best_rec) && best_rec.success
    best_rec.message = '0D current scan found no current meeting DeltaT and Qc tolerance';
end
end

function tf = is_better_n_current_scan_rec_ms(a, b)
a_dt = safe_numeric_field_ms(a, 'DeltaTN_actual', -inf);
b_dt = safe_numeric_field_ms(b, 'DeltaTN_actual', -inf);
if abs(a_dt - b_dt) > 1e-12
    tf = a_dt > b_dt;
    return;
end
a_qe = abs(safe_numeric_field_ms(a, 'Qc_error', inf));
b_qe = abs(safe_numeric_field_ms(b, 'Qc_error', inf));
if abs(a_qe - b_qe) > 1e-12
    tf = a_qe < b_qe;
    return;
end
tf = safe_numeric_field_ms(a, 'I_opt', inf) < safe_numeric_field_ms(b, 'I_opt', inf);
end

function [rec, cand] = evaluate_k25_n_vector_0d_at_current_ms(n_vec, G, spec, I_use)
G_eval = G;
G_eval.I = I_use;
spec_try = spec;
spec_try.k25_debug.total_n = sum(n_vec);
N = spec.stage_count;
cand = candidate_record_template_ms(N);
cand.candidate_id = -1;
cand.layout_method = '0d_count_only';
cand.stage_modes = repmat({'0d_count_only'}, 1, N);
cand.stage_methods = repmat({'0d_count_only'}, 1, N);
cand.stage_trends = repmat({'neutral'}, 1, N);
cand.symmetry_mode = 'none';
cand.edge_pattern_mode = 'none';
cand.n = n_vec(:).';
cand.ratios = cand.n(2:end) ./ max(cand.n(1:end-1), eps);
cand.I_opt = G_eval.I;
rec = empty_eval_struct_ms(N);
rec.n = n_vec(:).';
rec.I_opt = G_eval.I;
[ok_pair, npair] = particle_counts_to_pair_counts_ms(n_vec);
if ok_pair
    rec.npair = npair;
else
    rec.message = '0D candidate failed: invalid even particle counts';
    return;
end
if any(rec.n(1:end-1) < rec.n(2:end))
    rec.message = '0D candidate failed: particle counts must be nonincreasing';
    return;
end
rec.candidate_id = cand.candidate_id;
rec.layout_method = cand.layout_method;
rec.stage_modes = cand.stage_modes;
rec.stage_methods = cand.stage_methods;
rec.stage_trends = cand.stage_trends;
rec.symmetry_mode = cand.symmetry_mode;
rec.edge_pattern_mode = cand.edge_pattern_mode;
rec.s_dense = NaN;
rec.s_sparse = NaN;
rec.expo = NaN;
rec.anis_ratio = NaN;
rec.method_anchor_stage = NaN;
rec.spacing_ratio = NaN;
rec.contrast_score = NaN;
rec.mode_prior_score = NaN;
rec.ratios = cand.ratios;

x0 = linspace(spec.targets.ThRes - 20, spec.targets.ThRes - 120, spec.stage_count).';
[ok0d, C] = solve_0d_at_current_ms(n_vec, G_eval.I, G_eval, spec_try, x0);
if ~ok0d
    rec.message = '0D solve failed';
    return;
end
C = C(:).';
rec.C = C;
rec.stage_Tc_mean = C;
rec.stage_Th_mean = [spec.targets.ThRes, C(1:end-1)];
rec.stage_min_Th_minus_Tc = rec.stage_Th_mean - rec.stage_Tc_mean;
rec.has_temperature_reversal = any(rec.stage_min_Th_minus_Tc <= 0);
rec.sumQc_stage = NaN(1, N);
rec.sumQh_stage = NaN(1, N);
rec.Qc_pair_min = NaN(1, N);
rec.Qh_pair_min = NaN(1, N);
tol_neg_qc = get_qc_pair_negative_tol_ms(spec_try);
for k = 1:N
    Tc = C(k);
    if k == 1
        Th = spec.targets.ThRes;
    else
        Th = C(k-1);
    end
    Qc_pair = te_Qc_onecouple_ms(Tc, Th, G_eval.I, G_eval, k);
    Qh_pair = te_Qh_onecouple_ms(Th, Tc, G_eval.I, G_eval, k);
    if ~(isfinite(Qc_pair) && isfinite(Qh_pair))
        rec.message = sprintf('0D nonfinite TE heat at stage %d', k);
        return;
    end
    rec.sumQc_stage(k) = npair(k) * Qc_pair;
    rec.sumQh_stage(k) = npair(k) * Qh_pair;
    rec.Qc_pair_min(k) = Qc_pair;
    rec.Qh_pair_min(k) = Qh_pair;
    if Qc_pair < -tol_neg_qc
        rec.message = sprintf('0D stage%d_negative_Qc_pair', k);
        rec.physics_failure_reason = sprintf('stage%d_negative_Qc_pair', k);
        return;
    end
end
rec.DeltaTN_actual = spec.targets.ThRes - C(N);
rec.Qc_last_total = rec.sumQc_stage(N);
rec.Qc_target_last = spec.targets.Qc_target_last;
rec.Qc_error = rec.Qc_last_total - spec.targets.Qc_target_last;
rec.Qc_tol = get_n_current_scan_qc_tol_ms(spec);
rec.target_met = rec.DeltaTN_actual >= spec.targets.DeltaT_target && ...
    isfinite(rec.Qc_error) && abs(rec.Qc_error) <= rec.Qc_tol;
rec.TN_min = C(N);
rec.TN_mean = C(N);
rec.DeltaTN_mean = rec.DeltaTN_actual;
rec.TN_maxmin = 0;
rec.newton_convergence_mode = '0d';
rec.physics_valid = true;
rec.physics_failure_reason = 'ok';
rec.success = true;
rec.message = '0D ok';
rec.rank_score = rec.DeltaTN_actual;
end

function [ok, msg] = validate_k25_count_geometry_ms(cand, spec)
ok = false;
msg = '';
N = spec.stage_count;
if numel(cand.n) ~= N || any(~isfinite(cand.n)) || any(mod(cand.n, 2) ~= 0)
    msg = 'invalid particle counts';
    return;
end
if any(cand.n(1:end-1) < cand.n(2:end))
    msg = 'particle counts must be nonincreasing from hot to cold stages';
    return;
end
if any(cand.Lx > spec.geometry.L_max + 1e-12) || any(cand.Ly > spec.geometry.L_max + 1e-12)
    msg = sprintf('layout exceeds L_max: Lx_mm=%s, Ly_mm=%s, Lmax_mm=%s', ...
        vec_to_inline_str_ms(cand.Lx * 1e3), vec_to_inline_str_ms(cand.Ly * 1e3), ...
        vec_to_inline_str_ms(spec.geometry.L_max_mm));
    return;
end
if any(cand.cov + 1e-12 < spec.geometry.coverage_min)
    msg = sprintf('coverage below minimum: coverage=%s, coverage_min=%s', ...
        vec_to_inline_str_ms(cand.cov), vec_to_inline_str_ms(spec.geometry.coverage_min));
    return;
end
if any(cand.Lx(1:end-1) + 1e-12 < cand.Lx(2:end)) || ...
        any(cand.Ly(1:end-1) + 1e-12 < cand.Ly(2:end))
    msg = 'stage substrate size must be nonincreasing from hot to cold stages';
    return;
end
ok = true;
msg = 'ok';
end

function [best_rec, best_cand, best_n, rows, ratio_rec, ratio_cand, ratio_row, top_n_current] = run_particle_count_search_ms(G, spec)
N = spec.stage_count;
ps = spec.k25_debug.particle_search;
rows = repmat(empty_particle_search_row_ms(), 0, 1);
top_n_current = repmat(empty_top_n_current_row_ms(N), 0, 1);
ratio_rec = empty_eval_struct_ms(N);
ratio_cand = candidate_record_template_ms(N);
ratio_row = empty_particle_search_row_ms();
best_rec = empty_eval_struct_ms(N);
best_cand = candidate_record_template_ms(N);
best_n = NaN(1, N);
cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
eval_idx = 0;

n5_est = estimate_top_stage_particle_count_ms(G, spec);
n5_start = resolve_n5_min_for_search_ms(n5_est, spec);
ps.n5_estimate = n5_est;
ps.n5_min_effective = n5_start;
ps.n5_max_effective = resolve_n5_max_for_search_ms(n5_est, n5_start, spec);
spec.k25_debug.particle_search = ps;
total_step = max(2, 2 * round(ps.total_step / 2));
min_total = normalize_even_total_count_ms(2 * N, N);
min_total = max(min_total, n5_start + 2 * (N - 1));
if isfield(ps, 'total_n_min') && ~isempty(ps.total_n_min) && isfinite(ps.total_n_min)
    min_total = max(min_total, 2 * round(ps.total_n_min / 2));
end
min_total = total_step * ceil(min_total / total_step);
max_total = max(min_total, 2 * round(ps.total_n_max / 2));

fprintf(['[K25-Debug] 0D template search: total_n=[%d:%d:%d], block=%d, chunk=%d, ' ...
    'n5_est=%d, n5_min=%d, n5_max=%d, adjacent_min_fraction=%.4g\n'], ...
    min_total, total_step, max_total, ps.total_block_size, ps.par_chunk_size, ...
    n5_est, n5_start, ps.n5_max_effective, ps.adjacent_min_fraction);

best_0d_rec = empty_eval_struct_ms(N);
best_0d_cand = candidate_record_template_ms(N);
topK_target = get_n_current_scan_topK_ms(spec);
stop_search = false;
search_stop_reason = 'exhausted_total_n_range';
search_stop_total_n = NaN;
total_cursor = min_total;
while total_cursor <= max_total && ~stop_search
    total_block = total_cursor:total_step:min(max_total, total_cursor + total_step * (ps.total_block_size - 1));
    for tb = 1:numel(total_block)
        total_try = total_block(tb);
        [ok_try, rec_try, cand_try, n_try, cache, eval_idx, rows] = evaluate_0d_total_candidate_set_ms( ...
            total_try, n5_start, 'template_0d', cache, eval_idx, G, spec, rows);
        top_n_current = build_top_n_current_rows_ms(rows, spec);
        if ok_try
            if isempty(top_n_current)
                best_0d_rec = rec_try;
                best_0d_cand = cand_try;
                best_n = n_try;
            else
                best_n = top_n_current(1).n;
                [best_0d_rec, best_0d_cand] = evaluate_k25_n_vector_0d_ms(best_n, G, spec);
            end
        end
        if numel(top_n_current) >= topK_target
            stop_search = true;
            search_stop_reason = 'topK_feasible_found';
            search_stop_total_n = total_try;
            fprintf('[K25-Debug] stopping 0D search: found topK=%d feasible n+I rows at total_n=%d.\n', ...
                topK_target, total_try);
            break;
        end
    end
    total_cursor = total_block(end) + total_step;
end

top_n_current = build_top_n_current_rows_ms(rows, spec);
if isempty(top_n_current)
    best_rec.message = sprintf('0D template search found no target-met candidate for total_n=[%d,%d]', ...
        min_total, max_total);
    best_rec.particle_search_stop_reason = search_stop_reason;
    best_rec.particle_search_stop_total_n = search_stop_total_n;
    best_rec.particle_search_topK_target = topK_target;
    best_rec.particle_search_topK_found = 0;
    best_n = NaN(1, N);
    return;
end
best_n = top_n_current(1).n;
[best_0d_rec, best_0d_cand] = evaluate_k25_n_vector_0d_ms(best_n, G, spec);
G_best = G;
if isfield(best_0d_rec, 'I_opt') && isfinite(best_0d_rec.I_opt) && best_0d_rec.I_opt > 0
    G_best.I = best_0d_rec.I_opt;
end

if spec.k25_debug.search_only
    best_rec = best_0d_rec;
    best_cand = best_0d_cand;
    best_rec.message = 'ok: 0D search-only result; FEM skipped';
    if isfinite(sum(best_n))
        n_ratio = allocate_ratio_reference_counts_ms(sum(best_n), N);
        [ratio_rec, ratio_cand] = evaluate_k25_n_vector_0d_ms(n_ratio, G, spec);
        if ratio_rec.success && ratio_rec.target_met
            ratio_rec.message = 'ok: 0D ratio-reference result; FEM skipped';
        end
        ratio_row = build_particle_search_row_ms(ratio_rec, n_ratio, 'ratio_reference_0d', numel(rows) + 1, spec, NaN);
    end
else
    [best_rec, best_cand] = evaluate_k25_n_vector_ms(best_n, G_best, spec, true);
    if ~best_rec.success && has_valid_candidate_geometry_ms(best_0d_rec, N)
        best_rec.n = best_0d_rec.n;
        best_rec.npair = best_0d_rec.npair;
        best_rec.Lx = best_0d_rec.Lx;
        best_rec.Ly = best_0d_rec.Ly;
        best_rec.cov = best_0d_rec.cov;
        best_rec.stage_rects = best_0d_rec.stage_rects;
    end

    if isfinite(sum(best_n))
        n_ratio = allocate_ratio_reference_counts_ms(sum(best_n), N);
        [ratio_rec, ratio_cand] = evaluate_k25_n_vector_ms(n_ratio, G_best, spec, true);
        ratio_row = build_particle_search_row_ms(ratio_rec, n_ratio, 'ratio_reference', numel(rows) + 1, spec, NaN);
    end
end
if isstruct(best_rec)
    best_rec.selection_basis = '0D';
    best_rec.fem_is_acceptance_criterion = false;
    best_rec.selected_0d_success = logical(safe_numeric_field_ms(best_0d_rec, 'success', false));
    best_rec.selected_0d_target_met = logical(get_target_met_ms(best_0d_rec, spec));
    best_rec.selected_0d_I_opt = safe_numeric_field_ms(best_0d_rec, 'I_opt', NaN);
    best_rec.selected_0d_DeltaTN_actual = safe_numeric_field_ms(best_0d_rec, 'DeltaTN_actual', NaN);
    best_rec.selected_0d_Qc_last_total = safe_numeric_field_ms(best_0d_rec, 'Qc_last_total', NaN);
    best_rec.selected_0d_Qc_error = safe_numeric_field_ms(best_0d_rec, 'Qc_error', NaN);
    best_rec.selected_0d_Qc_tol = safe_numeric_field_ms(best_0d_rec, 'Qc_tol', NaN);
    best_rec.particle_search_n5_estimate = n5_est;
    best_rec.particle_search_n5_min = n5_start;
    best_rec.particle_search_n5_max_effective = ps.n5_max_effective;
    if isfield(ps, 'total_n_min')
        best_rec.particle_search_total_n_min = ps.total_n_min;
    end
    best_rec.particle_search_total_step = total_step;
    best_rec.particle_search_use_ratio_prune = logical(ps.use_ratio_prune);
    best_rec.particle_search_adjacent_min_fraction = ps.adjacent_min_fraction;
    best_rec.particle_search_stop_reason = search_stop_reason;
    best_rec.particle_search_stop_total_n = search_stop_total_n;
    best_rec.particle_search_topK_target = topK_target;
    best_rec.particle_search_topK_found = numel(top_n_current);
end
end

function [ok, best_rec, best_cand, best_n, cache, eval_idx, rows] = evaluate_0d_total_candidate_set_ms( ...
    total_n, n5_min, phase, cache, eval_idx, G, spec, rows)
ok = false;
best_rec = empty_eval_struct_ms(spec.stage_count);
best_cand = candidate_record_template_ms(spec.stage_count);
best_n = NaN(1, spec.stage_count);
n_mat = generate_0d_count_candidates_for_total_ms(total_n, n5_min, spec);
if isempty(n_mat)
    fprintf('[K25-Debug] 0D total_n=%d: candidates=0\n', total_n);
    return;
end
fprintf('[K25-Debug] 0D total_n=%d: candidates=%d\n', total_n, size(n_mat, 1));
[rec_list, cache, eval_idx, rows] = eval_particle_candidate_matrix_chunked_ms( ...
    cache, n_mat, phase, eval_idx, G, spec, rows);
for i = 1:numel(rec_list)
    rec_i = rec_list{i};
    if ~is_particle_target_met_ms(rec_i, spec)
        continue;
    end
    n_i = n_mat(i,:);
    if ~ok || is_better_particle_solution_ms(rec_i, n_i, best_rec, best_n, spec)
        [~, cand_i] = evaluate_k25_n_vector_0d_ms(n_i, G, spec);
        best_rec = rec_i;
        best_cand = cand_i;
        best_n = n_i;
        ok = true;
    end
end
end

function n_mat = generate_0d_count_candidates_for_total_ms(total_n, n5_min, spec)
N = spec.stage_count;
ps = spec.k25_debug.particle_search;
total_n = normalize_even_total_count_ms(total_n, N);
n5_min = max(2, 2 * ceil(n5_min / 2));
if total_n < n5_min + 2 * (N - 1)
    n_mat = zeros(0, N);
    return;
end
n5_est = resolve_n5_estimate_for_search_ms(n5_min, spec);
n5_max = min(total_n - 2 * (N - 1), resolve_n5_max_for_search_ms(n5_est, n5_min, spec));
n5_max = max(n5_min, 2 * floor(n5_max / 2));
if n5_max < n5_min
    n_mat = zeros(0, N);
    return;
end

rows = zeros(0, N);
templates = [
    55, 21, 7, 3, 1;
    42, 38, 13, 6, 2;
    32, 24, 12, 6, 2;
    24, 20, 12, 7, 3;
    16, 14, 10, 7, 4;
    10, 9, 7, 5, 3;
    6, 5, 4, 3, 2;
    1, 1, 1, 1, 1];
templates = templates(:, 1:N);
scales = [0.70, 0.85, 1.00, 1.18, 1.40];
for ti = 1:size(templates, 1)
    base = templates(ti,:);
    for si = 1:numel(scales)
        w = base .^ scales(si);
        raw = total_n * w / sum(w);
        n = max(2, 2 * round(raw / 2));
        n = enforce_nonincreasing_even_counts_ms(n, N);
        n(end) = max(n(end), n5_min);
        n = enforce_nonincreasing_even_counts_ms(n, N);
        n = rebalance_counts_to_total_preserve_order_ms(n, total_n);
        if is_valid_0d_template_candidate_ms(n, total_n, n5_min, spec)
            rows(end+1,:) = n; %#ok<AGROW>
        end
    end
end

n5_points = max(3, min(31, max(3, round(ps.candidate_budget / 20))));
n5_vals = unique(2 * round(linspace(n5_min, n5_max, n5_points) / 2));
n5_vals = n5_vals(isfinite(n5_vals) & n5_vals >= n5_min & n5_vals <= n5_max);
weights4 = [
    55, 21, 7, 3;
    42, 38, 13, 6;
    32, 24, 12, 6;
    24, 20, 12, 7;
    16, 14, 10, 7;
    10, 9, 7, 5;
    1, 1, 1, 1];
for i = 1:numel(n5_vals)
    n5 = n5_vals(i);
    rem_total = total_n - n5;
    if rem_total < 2 * (N - 1)
        continue;
    end
    for wi = 1:size(weights4, 1)
        raw4 = rem_total * weights4(wi,1:N-1) / sum(weights4(wi,1:N-1));
        n = [max(2, 2 * round(raw4 / 2)), n5];
        n = enforce_nonincreasing_even_counts_ms(n, N);
        n = rebalance_counts_to_total_preserve_order_ms(n, total_n);
        if is_valid_0d_template_candidate_ms(n, total_n, n5_min, spec)
            rows(end+1,:) = n; %#ok<AGROW>
        end
    end
end

if isempty(rows)
    n_mat = zeros(0, N);
    return;
end
rows = unique(rows, 'rows', 'stable');
score = abs(rows(:,end) - n5_min) + 0.001 * sum(abs(diff(rows, 1, 2)), 2);
[~, ord] = sort(score, 'ascend');
ord = ord(1:min(numel(ord), ps.candidate_budget));
n_mat = rows(ord,:);
end

function tf = is_valid_0d_template_candidate_ms(n, total_n, n5_min, spec)
tf = ~isempty(n) && all(isfinite(n)) && all(n >= 2) && ...
    all(mod(round(n), 2) == 0) && sum(n) == total_n && ...
    n(end) >= n5_min && all(n(1:end-1) >= n(2:end)) && ...
    all(n <= spec.k25_debug.particle_search.stage_max_count) && ...
    n(end) <= resolve_n5_max_for_search_ms(resolve_n5_estimate_for_search_ms(n5_min, spec), n5_min, spec) && ...
    passes_particle_ratio_prune_ms(n, spec);
end

function [rec, cand, row, cache, eval_idx, rows] = eval_particle_candidate_cached_ms( ...
    cache, n_vec, phase, eval_idx, G, spec, rows)
key = vec_to_count_key_ms(n_vec);
if isKey(cache, key)
    pack = cache(key);
    rec = pack.rec;
    cand = pack.cand;
    row = pack.row;
    return;
end
eval_idx = eval_idx + 1;
t_one = tic;
[rec, cand] = evaluate_k25_n_vector_0d_ms(n_vec, G, spec);
row = build_particle_search_row_ms(rec, n_vec, phase, eval_idx, spec, toc(t_one));
cache(key) = struct('rec', rec, 'cand', cand, 'row', row);
rows(end+1,1) = row;
end

function [rec_list, cache, eval_idx, rows] = eval_particle_candidate_matrix_chunked_ms( ...
    cache, n_mat, phase, eval_idx, G, spec, rows)
n_mat = unique(round(n_mat), 'rows', 'stable');
n_total = size(n_mat, 1);
rec_list = cell(n_total, 1);
if n_total < 1
    return;
end
chunk_size = max(1, round(spec.k25_debug.particle_search.par_chunk_size));
for i0 = 1:chunk_size:n_total
    i1 = min(n_total, i0 + chunk_size - 1);
    [rec_chunk, cache, eval_idx, rows] = eval_particle_candidate_batch_cached_ms( ...
        cache, n_mat(i0:i1,:), phase, eval_idx, G, spec, rows);
    rec_list(i0:i1) = rec_chunk;
end
end

function [rec_list, cache, eval_idx, rows] = eval_particle_candidate_batch_cached_ms( ...
    cache, n_mat, phase, eval_idx, G, spec, rows)
n_mat = round(n_mat);
n_batch = size(n_mat, 1);
rec_list = cell(n_batch, 1);
pending_idx = [];
for i = 1:n_batch
    key = vec_to_count_key_ms(n_mat(i,:));
    if isKey(cache, key)
        pack = cache(key);
        rec_list{i} = pack.rec;
    else
        pending_idx(end+1) = i; %#ok<AGROW>
    end
end
if isempty(pending_idx)
    return;
end
pending_n = n_mat(pending_idx, :);
n_pending = size(pending_n, 1);
base_idx = eval_idx;
eval_idx = eval_idx + n_pending;
rec_pending = cell(n_pending, 1);
cand_pending = cell(n_pending, 1);
row_pending = repmat(empty_particle_search_row_ms(), n_pending, 1);
use_par = logical(spec.use_parallel) && n_pending >= 2 && get_parallel_worker_count_ms() > 0;
if use_par
    parfor ii = 1:n_pending
        t_one = tic;
        [rec_i, cand_i] = evaluate_k25_n_vector_0d_ms(pending_n(ii,:), G, spec);
        rec_pending{ii} = rec_i;
        cand_pending{ii} = cand_i;
        row_pending(ii) = build_particle_search_row_ms(rec_i, pending_n(ii,:), phase, base_idx + ii, spec, toc(t_one));
    end
else
    for ii = 1:n_pending
        t_one = tic;
        [rec_i, cand_i] = evaluate_k25_n_vector_0d_ms(pending_n(ii,:), G, spec);
        rec_pending{ii} = rec_i;
        cand_pending{ii} = cand_i;
        row_pending(ii) = build_particle_search_row_ms(rec_i, pending_n(ii,:), phase, base_idx + ii, spec, toc(t_one));
    end
end
for ii = 1:n_pending
    idx = pending_idx(ii);
    key = vec_to_count_key_ms(pending_n(ii,:));
    cache(key) = struct('rec', rec_pending{ii}, 'cand', cand_pending{ii}, 'row', row_pending(ii));
    rec_list{idx} = rec_pending{ii};
    rows(end+1,1) = row_pending(ii); %#ok<AGROW>
end
end

function n_vec = enforce_nonincreasing_even_counts_ms(n_vec, N)
n_vec = fit_len_vec_ms(n_vec, N, 2 * ones(1, N));
n_vec = max(2, 2 * round(n_vec / 2));
for k = N-1:-1:1
    n_vec(k) = max(n_vec(k), n_vec(k+1));
end
end

function n5 = estimate_top_stage_particle_count_ms(G, spec)
N = spec.stage_count;
if spec.targets.Qc_target_last <= 0
    n5 = 2;
    return;
end
Tc_est = spec.targets.ThRes - spec.targets.DeltaT_target;
Th_est = spec.targets.ThRes - spec.targets.DeltaT_target * max(N - 1, 1) / max(N, 1);
if ~(isfinite(Tc_est) && isfinite(Th_est) && Th_est > Tc_est)
    Th_est = spec.targets.ThRes;
end
Qc_pair_est = te_Qc_onecouple_ms(Tc_est, Th_est, G.I, G, N);
if ~(isfinite(Qc_pair_est) && Qc_pair_est > 0)
    Qc_pair_est = max(spec.targets.Qc_target_last, 1e-6);
end
npair5 = max(1, ceil(spec.targets.Qc_target_last / max(Qc_pair_est, eps)));
n5 = 2 * npair5;
n5 = max(n5, 2);
end

function n5_min = resolve_n5_min_for_search_ms(n5_est, spec)
n5_est = max(2, 2 * ceil(n5_est / 2));
factor = 0.7;
if isstruct(spec) && isfield(spec, 'k25_debug') && isstruct(spec.k25_debug) && ...
        isfield(spec.k25_debug, 'particle_search') && isstruct(spec.k25_debug.particle_search) && ...
        isfield(spec.k25_debug.particle_search, 'n5_min_factor')
    factor = scalar_with_default_ms(spec.k25_debug.particle_search.n5_min_factor, factor);
end
factor = min(max(factor, 0.1), 1.0);
n5_min = max(2, 2 * floor(factor * n5_est / 2));
end

function n5_est = resolve_n5_estimate_for_search_ms(n5_min, spec)
n5_est = n5_min;
if isstruct(spec) && isfield(spec, 'k25_debug') && isstruct(spec.k25_debug) && ...
        isfield(spec.k25_debug, 'particle_search') && isstruct(spec.k25_debug.particle_search) && ...
        isfield(spec.k25_debug.particle_search, 'n5_estimate') && ...
        ~isempty(spec.k25_debug.particle_search.n5_estimate) && ...
        isfinite(spec.k25_debug.particle_search.n5_estimate)
    n5_est = spec.k25_debug.particle_search.n5_estimate;
end
n5_est = max(n5_min, 2 * ceil(n5_est / 2));
end

function n5_max = resolve_n5_max_for_search_ms(n5_est, n5_min, spec)
n5_est = max(2, 2 * ceil(n5_est / 2));
n5_min = max(2, 2 * ceil(n5_min / 2));
ps = spec.k25_debug.particle_search;
if isfield(ps, 'n5_max') && ~isempty(ps.n5_max) && isfinite(ps.n5_max)
    n5_max = max(n5_min, 2 * round(ps.n5_max / 2));
    return;
end
factor = 2.5;
margin = 20;
if isfield(ps, 'n5_max_auto_factor')
    factor = max(1, scalar_with_default_ms(ps.n5_max_auto_factor, factor));
end
if isfield(ps, 'n5_max_auto_margin')
    margin = max(0, 2 * round(scalar_with_default_ms(ps.n5_max_auto_margin, margin) / 2));
end
n5_max = max(n5_est + margin, 2 * ceil(factor * n5_est / 2));
n5_max = max(n5_min, 2 * round(n5_max / 2));
end

function tf = passes_particle_ratio_prune_ms(n_vec, spec)
tf = true;
if ~isstruct(spec) || ~isfield(spec, 'k25_debug') || ...
        ~isfield(spec.k25_debug, 'particle_search') || ...
        ~isstruct(spec.k25_debug.particle_search)
    return;
end
ps = spec.k25_debug.particle_search;
if ~isfield(ps, 'use_ratio_prune') || ~logical(ps.use_ratio_prune)
    return;
end
frac = 0.1;
if isfield(ps, 'adjacent_min_fraction')
    frac = scalar_with_default_ms(ps.adjacent_min_fraction, frac);
end
frac = min(max(frac, 0), 1);
if frac <= 0
    return;
end
n_vec = n_vec(:).';
if numel(n_vec) < 2 || any(~isfinite(n_vec)) || any(n_vec <= 0)
    tf = false;
    return;
end
tf = all(n_vec(2:end) ./ max(n_vec(1:end-1), eps) >= frac - 1e-12);
end

function tf = is_particle_target_met_ms(rec, spec)
tf = get_target_met_ms(rec, spec);
end

function tf = get_target_met_ms(rec, spec)
if isstruct(rec) && isfield(rec, 'target_met') && ~isempty(rec.target_met)
    tf = logical(rec.target_met);
    return;
end
tf = isstruct(rec) && isfield(rec, 'success') && rec.success && ...
    isfield(rec, 'DeltaTN_actual') && isfinite(rec.DeltaTN_actual) && ...
    rec.DeltaTN_actual >= spec.targets.DeltaT_target;
end

function tf = is_better_particle_solution_ms(rec, n_vec, best_rec, best_n, spec)
tf = false;
if ~is_particle_target_met_ms(rec, spec)
    return;
end
if nargin < 4 || isempty(best_n) || any(~isfinite(best_n))
    tf = true;
    return;
end
sum_new = sum(n_vec);
sum_old = sum(best_n);
if sum_new < sum_old
    tf = true;
elseif sum_new == sum_old
    margin_new = rec.DeltaTN_actual - spec.targets.DeltaT_target;
    margin_old = -inf;
    if isstruct(best_rec) && isfield(best_rec, 'DeltaTN_actual') && isfinite(best_rec.DeltaTN_actual)
        margin_old = best_rec.DeltaTN_actual - spec.targets.DeltaT_target;
    end
    tf = margin_new > margin_old;
end
end

function n_vec = allocate_ratio_reference_counts_ms(total_n, N)
ratio = [55, 21, 7, 3, 1];
ratio = fit_len_vec_ms(ratio, N, ratio);
raw = total_n * ratio / sum(ratio);
n_vec = max(2, 2 * round(raw / 2));
n_vec = rebalance_counts_to_total_ms(n_vec, raw, total_n, true);
n_vec = enforce_nonincreasing_even_counts_ms(n_vec, N);
n_vec = rebalance_counts_to_total_preserve_order_ms(n_vec, total_n);
end

function n_vec = rebalance_counts_to_total_preserve_order_ms(n_vec, total_n)
n_vec = max(2, 2 * round(n_vec(:).' / 2));
total_n = 2 * round(total_n / 2);
N = numel(n_vec);
guard = 0;
while sum(n_vec) < total_n && guard < 10000
    room = [inf, n_vec(1:end-1) - n_vec(2:end)];
    idx = find(room >= 2, 1, 'last');
    if isempty(idx)
        idx = 1;
    end
    n_vec(idx) = n_vec(idx) + 2;
    for k = N-1:-1:1
        n_vec(k) = max(n_vec(k), n_vec(k+1));
    end
    guard = guard + 1;
end
guard = 0;
while sum(n_vec) > total_n && guard < 10000
    floor_next = [n_vec(2:end), 2];
    room = n_vec - floor_next;
    idx = find(room >= 2, 1, 'first');
    if isempty(idx)
        break;
    end
    n_vec(idx) = n_vec(idx) - 2;
    guard = guard + 1;
end
end

function key = vec_to_count_key_ms(n_vec)
key = strtrim(sprintf('%d_', round(n_vec)));
end

function row = empty_particle_search_row_ms()
row = struct('eval_index', NaN, 'phase', '', 'total_n', NaN, ...
    'n', '[]', 'n_vec_numeric', NaN(1,0), 'npair', '[]', ...
    'success', false, 'physics_valid', false, 'target_met', false, ...
    'target_met_0d', false, 'C', '[]', 'DeltaTN_0d', NaN, ...
    'DeltaTN_actual', NaN, 'DeltaTN_mean', NaN, 'TN_min', NaN, ...
    'I_A', NaN, 'Qc_last_total', NaN, 'Qc_target_last', NaN, ...
    'Qc_error', NaN, 'Qc_tol', NaN, 'message', '', 'physics_failure_reason', '', ...
    'Lx_mm', '[]', 'Ly_mm', '[]', 'coverage', '[]', ...
    'Qc_pair_min', '[]', 'Qh_pair_min', '[]', 'elapsed_sec', NaN);
end

function row = build_particle_search_row_ms(rec, n_vec, phase, eval_index, spec, elapsed_sec)
row = empty_particle_search_row_ms();
row.eval_index = eval_index;
row.phase = char(string(phase));
row.n_vec_numeric = n_vec(:).';
row.n = vec_to_inline_str_ms(n_vec);
row.total_n = sum(n_vec);
[ok_pair, npair] = particle_counts_to_pair_counts_ms(n_vec);
if ok_pair
    row.npair = vec_to_inline_str_ms(npair);
end
row.target_met = is_particle_target_met_ms(rec, spec);
row.target_met_0d = false;
row.elapsed_sec = elapsed_sec;
if ~isstruct(rec)
    return;
end
if isfield(rec, 'success'), row.success = logical(rec.success); end
if isfield(rec, 'physics_valid'), row.physics_valid = logical(rec.physics_valid); end
if isfield(rec, 'target_met'), row.target_met = logical(rec.target_met); end
if isfield(rec, 'C'), row.C = vec_to_inline_str_ms(rec.C); end
if isfield(rec, 'DeltaTN_actual'), row.DeltaTN_actual = rec.DeltaTN_actual; end
if isfield(rec, 'newton_convergence_mode') && strcmpi(char(string(rec.newton_convergence_mode)), '0d')
    row.DeltaTN_0d = row.DeltaTN_actual;
    row.target_met_0d = row.target_met;
end
if isfield(rec, 'DeltaTN_mean'), row.DeltaTN_mean = rec.DeltaTN_mean; end
if isfield(rec, 'TN_min'), row.TN_min = rec.TN_min; end
if isfield(rec, 'I_opt'), row.I_A = rec.I_opt; end
if isfield(rec, 'Qc_last_total'), row.Qc_last_total = rec.Qc_last_total; end
if isfield(rec, 'Qc_target_last'), row.Qc_target_last = rec.Qc_target_last; end
if isfield(rec, 'Qc_error'), row.Qc_error = rec.Qc_error; end
if isfield(rec, 'Qc_tol'), row.Qc_tol = rec.Qc_tol; end
if isfield(rec, 'message'), row.message = char(string(rec.message)); end
if isfield(rec, 'physics_failure_reason'), row.physics_failure_reason = char(string(rec.physics_failure_reason)); end
if isfield(rec, 'Lx'), row.Lx_mm = vec_to_inline_str_ms(rec.Lx * 1e3); end
if isfield(rec, 'Ly'), row.Ly_mm = vec_to_inline_str_ms(rec.Ly * 1e3); end
if isfield(rec, 'cov'), row.coverage = vec_to_inline_str_ms(rec.cov); end
if isfield(rec, 'Qc_pair_min'), row.Qc_pair_min = vec_to_inline_str_ms(rec.Qc_pair_min); end
if isfield(rec, 'Qh_pair_min'), row.Qh_pair_min = vec_to_inline_str_ms(rec.Qh_pair_min); end
end

function write_particle_search_results_csv_ms(rows, out_csv)
if isempty(rows)
    return;
end
rows_out = rmfield_if_exists_ms(rows, 'n_vec_numeric');
write_struct_rows_csv_ms(rows_out, out_csv);
end

function row = empty_top_n_current_row_ms(N)
if nargin < 1 || ~isfinite(N) || N < 1
    N = 5;
end
row = struct('rank', NaN, 'n', NaN(1,N), 'n_text', '[]', ...
    'I_opt', NaN, 'DeltaTN_actual', NaN, 'Qc_last_total', NaN, ...
    'Qc_error', NaN, 'total_n', NaN, 'score', NaN);
end

function top_rows = build_top_n_current_rows_ms(rows, spec)
N = spec.stage_count;
top_rows = repmat(empty_top_n_current_row_ms(N), 0, 1);
if isempty(rows)
    return;
end
valid = [rows.target_met] & isfinite([rows.I_A]) & ...
    isfinite([rows.DeltaTN_actual]) & isfinite([rows.Qc_error]);
if ~any(valid)
    return;
end
idx = find(valid);
dt = [rows(idx).DeltaTN_actual].';
qe = abs([rows(idx).Qc_error].');
tn = [rows(idx).total_n].';
ii = [rows(idx).I_A].';
[~, ord] = sortrows([tn, qe, -dt, ii]);
idx = idx(ord);
topK = min(get_n_current_scan_topK_ms(spec), numel(idx));
top_rows = repmat(empty_top_n_current_row_ms(N), topK, 1);
for r = 1:topK
    src = rows(idx(r));
    top_rows(r).rank = r;
    top_rows(r).n = fit_count_vec_for_top_row_ms(src.n_vec_numeric, N);
    top_rows(r).n_text = vec_to_inline_str_ms(top_rows(r).n);
    top_rows(r).I_opt = src.I_A;
    top_rows(r).DeltaTN_actual = src.DeltaTN_actual;
    top_rows(r).Qc_last_total = src.Qc_last_total;
    top_rows(r).Qc_error = src.Qc_error;
    top_rows(r).total_n = src.total_n;
    top_rows(r).score = -src.total_n;
end
end

function write_top_n_current_csv_ms(top_rows, out_csv)
if isempty(top_rows)
    return;
end
rows_out = top_rows;
for i = 1:numel(rows_out)
    rows_out(i).n = rows_out(i).n_text;
end
write_struct_rows_csv_ms(rows_out, out_csv);
end

function n = fit_count_vec_for_top_row_ms(v, N)
n = NaN(1, N);
if isnumeric(v) && ~isempty(v)
    m = min(N, numel(v));
    n(1:m) = v(1:m);
end
end

function tf = is_n_current_scan_enabled_ms(spec)
tf = isfield(spec, 'n_current_scan') && isstruct(spec.n_current_scan) && ...
    isfield(spec.n_current_scan, 'enable') && any(logical(spec.n_current_scan.enable(:)));
end

function I_list = get_n_current_scan_I_list_ms(spec)
I_list = 1.0:0.2:8.0;
if isfield(spec, 'n_current_scan') && isstruct(spec.n_current_scan) && ...
        isfield(spec.n_current_scan, 'I_list') && ~isempty(spec.n_current_scan.I_list)
    I_list = spec.n_current_scan.I_list;
end
I_list = unique(double(I_list(:).'));
I_list = I_list(isfinite(I_list) & I_list > 0);
if isempty(I_list)
    I_list = 1.0:0.2:8.0;
end
end

function topK = get_n_current_scan_topK_ms(spec)
topK = 3;
if isfield(spec, 'n_current_scan') && isstruct(spec.n_current_scan) && ...
        isfield(spec.n_current_scan, 'topK') && ~isempty(spec.n_current_scan.topK)
    topK = spec.n_current_scan.topK(1);
end
if ~(isfinite(topK) && topK >= 1)
    topK = 3;
end
topK = max(1, round(topK));
end

function qc_tol = get_n_current_scan_qc_tol_ms(spec)
abs_tol = 1e-3;
rel_tol = 1e-3;
if isfield(spec, 'n_current_scan') && isstruct(spec.n_current_scan)
    if isfield(spec.n_current_scan, 'qc_abs_tol') && ~isempty(spec.n_current_scan.qc_abs_tol)
        abs_tol = spec.n_current_scan.qc_abs_tol(1);
    end
    if isfield(spec.n_current_scan, 'qc_rel_tol') && ~isempty(spec.n_current_scan.qc_rel_tol)
        rel_tol = spec.n_current_scan.qc_rel_tol(1);
    end
end
abs_tol = max(0, double(abs_tol));
rel_tol = max(0, double(rel_tol));
qc_tol = max(abs_tol, rel_tol * abs(spec.targets.Qc_target_last));
end

function v = safe_numeric_field_ms(s, fname, fallback)
v = fallback;
if isstruct(s) && isfield(s, fname) && ~isempty(s.(fname)) && (isnumeric(s.(fname)) || islogical(s.(fname)))
    x = s.(fname);
    if isfinite(x(1))
        v = x(1);
    end
end
end

function s = rmfield_if_exists_ms(s, fname)
if isstruct(s) && isfield(s, fname)
    s = rmfield(s, fname);
end
end

function write_best_vs_ratio_reference_csv_ms(best_rec, ratio_rec, n_best, n_ratio, out_csv, spec)
rows = repmat(empty_compare_solution_row_ms(), 2, 1);
rows(1) = build_compare_solution_row_ms('best', best_rec, n_best, spec);
rows(2) = build_compare_solution_row_ms('ratio_reference', ratio_rec, n_ratio, spec);
write_struct_rows_csv_ms(rows, out_csv);
end

function write_struct_rows_csv_ms(rows, out_csv)
if nargin < 2 || isempty(out_csv) || ~isstruct(rows) || isempty(rows)
    return;
end
fn = fieldnames(rows);
vals = cell(numel(rows), numel(fn));
for r = 1:numel(rows)
    for c = 1:numel(fn)
        vals{r,c} = scalar_csv_value_ms(rows(r).(fn{c}));
    end
end
T = cell2table(vals, 'VariableNames', fn);
writetable(T, out_csv);
end

function row = empty_compare_solution_row_ms()
row = struct('label', '', 'n', '[]', 'total_n', NaN, 'success', false, ...
    'target_met', false, 'DeltaTN_actual', NaN, 'DeltaTN_mean', NaN, ...
    'TN_min', NaN, 'TN_mean', NaN, 'I_A', NaN, 'Lx_mm', '[]', ...
    'Ly_mm', '[]', 'coverage', '[]', 'message', '');
end

function row = build_compare_solution_row_ms(label, rec, n_vec, spec)
row = empty_compare_solution_row_ms();
row.label = char(string(label));
row.n = vec_to_inline_str_ms(n_vec);
row.total_n = sum(n_vec);
row.target_met = is_particle_target_met_ms(rec, spec);
if ~isstruct(rec)
    return;
end
if isfield(rec, 'success'), row.success = logical(rec.success); end
if isfield(rec, 'target_met'), row.target_met = logical(rec.target_met); end
if isfield(rec, 'DeltaTN_actual'), row.DeltaTN_actual = rec.DeltaTN_actual; end
if isfield(rec, 'DeltaTN_mean'), row.DeltaTN_mean = rec.DeltaTN_mean; end
if isfield(rec, 'TN_min'), row.TN_min = rec.TN_min; end
if isfield(rec, 'TN_mean'), row.TN_mean = rec.TN_mean; end
if isfield(rec, 'I_opt'), row.I_A = rec.I_opt; end
if isfield(rec, 'Lx'), row.Lx_mm = vec_to_inline_str_ms(rec.Lx * 1e3); end
if isfield(rec, 'Ly'), row.Ly_mm = vec_to_inline_str_ms(rec.Ly * 1e3); end
if isfield(rec, 'cov'), row.coverage = vec_to_inline_str_ms(rec.cov); end
if isfield(rec, 'message'), row.message = char(string(rec.message)); end
end

function total_n = normalize_even_total_count_ms(total_n_in, N)
if nargin < 2 || ~isfinite(N)
    N = 1;
end
min_total = 2 * max(1, round(N));
if nargin < 1 || isempty(total_n_in) || ~isfinite(total_n_in)
    total_n = min_total;
else
    total_n = max(min_total, round(total_n_in));
end
if mod(total_n, 2) ~= 0
    total_n = total_n + 1;
end
end

function [ok, cand, msg] = build_k25_standard_candidate_ms(spec, n_vec, I_use)
ok = false;
msg = '';
N = spec.stage_count;
cand = candidate_record_template_ms(N);
if numel(n_vec) ~= N || sum(n_vec) ~= spec.k25_debug.total_n
    msg = 'invalid k25 standard count vector';
    return;
end
rects = cell(1, N);
Lx_raw = NaN(1, N);
Ly_raw = NaN(1, N);
for k = 1:N
    [ok_k, rects_k, Lx_k, Ly_k] = make_k25_standard_rects_ms(n_vec(k), spec);
    if ~ok_k
        msg = sprintf('stage%d standard layout failed', k);
        return;
    end
    rects{k} = rects_k;
    Lx_raw(k) = Lx_k;
    Ly_raw(k) = Ly_k;
end
Lx = Lx_raw;
Ly = Ly_raw;
for k = N-1:-1:1
    Lx(k) = max(Lx(k), Lx(k+1) + spec.geometry.pyramid_gap_min(k));
    Ly(k) = max(Ly(k), Ly(k+1) + spec.geometry.pyramid_gap_min(k));
end
cov = zeros(1, N);
for k = 1:N
    cov(k) = n_vec(k) * spec.fp_w * spec.fp_h / max(Lx(k) * Ly(k), eps);
end
cand.candidate_id = 1;
cand.layout_method = 'k25_standard_fixed50';
cand.stage_modes = repmat({'standard_fixed50'}, 1, N);
cand.stage_methods = repmat({'k25_standard_fixed50'}, 1, N);
cand.stage_trends = repmat({'neutral'}, 1, N);
cand.symmetry_mode = 'none';
cand.edge_pattern_mode = 'two_sides_to_center';
cand.s_dense = max(spec.fp_w, spec.fp_h) * 1.5;
cand.s_sparse = cand.s_dense;
cand.expo = 1.0;
cand.anis_ratio = 1.0;
cand.method_anchor_stage = NaN;
cand.spacing_ratio = 1.0;
cand.contrast_score = calc_contrast_score_ms(cand.s_dense, cand.s_sparse, cand.expo);
cand.mode_prior_score = NaN;
cand.lmax_relaxed = false;
cand.lmax_relax_ratio = 1.0;
cand.n = n_vec(:).';
cand.ratios = cand.n(2:end) ./ max(cand.n(1:end-1), eps);
cand.I_opt = I_use;
cand.Lx = Lx;
cand.Ly = Ly;
cand.cov = cov;
cand.stage_rects = rects;
ok = true;
end

function [ok, rects, Lx_box, Ly_box] = make_k25_standard_rects_ms(n_target, spec)
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
[nx, ny] = choose_k25_standard_grid_shape_ms(n_target, spec.fp_w, spec.fp_h);
pts = fill_k25_grid_two_sides_to_center_ms(nx, ny, n_target, pitch_x, pitch_y);
if size(pts,1) ~= n_target
    return;
end
rects = centers_to_rects_ms(pts, spec.fp_w, spec.fp_h);
xmin = min(rects(:,1)); xmax = max(rects(:,2));
ymin = min(rects(:,3)); ymax = max(rects(:,4));
Lx_box = (xmax - xmin) + 2 * margin_x;
Ly_box = (ymax - ymin) + 2 * margin_y;
ok = isfinite(Lx_box) && isfinite(Ly_box) && Lx_box > 0 && Ly_box > 0;
end

function [nx_best, ny_best] = choose_k25_standard_grid_shape_ms(n_target, fp_w, fp_h)
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

function pts = fill_k25_grid_two_sides_to_center_ms(nx, ny, n_target, pitch_x, pitch_y)
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

function row = build_k25_single_result_row_ms(rec, spec, total_n, elapsed_sec)
row = struct('total_n', total_n, 'n', vec_to_inline_str_ms(rec.n), ...
    'npair', vec_to_inline_str_ms(rec.npair), ...
    'k_inplane', spec.plate_k_inplane, 'I_A', rec.I_opt, ...
    'candidate_id', rec.candidate_id, 'success', logical(rec.success), ...
    'target_met', logical(get_target_met_ms(rec, spec)), ...
    'message', char(string(rec.message)), ...
    'newton_failure_detail', char(string(rec.newton_failure_detail)), ...
    'physics_valid', logical(rec.physics_valid), ...
    'physics_failure_reason', char(string(rec.physics_failure_reason)), ...
    'stage_Tc_mean', vec_to_inline_str_ms(rec.stage_Tc_mean), ...
    'stage_Th_mean', vec_to_inline_str_ms(rec.stage_Th_mean), ...
    'stage_min_Th_minus_Tc', vec_to_inline_str_ms(rec.stage_min_Th_minus_Tc), ...
    'has_temperature_reversal', logical(rec.has_temperature_reversal), ...
    'sumQc_stage', vec_to_inline_str_ms(rec.sumQc_stage), ...
    'sumQh_stage', vec_to_inline_str_ms(rec.sumQh_stage), ...
    'Qc_pair_min', vec_to_inline_str_ms(rec.Qc_pair_min), ...
    'Qh_pair_min', vec_to_inline_str_ms(rec.Qh_pair_min), ...
    'DeltaTN_actual', rec.DeltaTN_actual, 'DeltaTN_mean', rec.DeltaTN_mean, ...
    'Qc_last_total', rec.Qc_last_total, 'Qc_target_last', rec.Qc_target_last, ...
    'Qc_error', rec.Qc_error, 'Qc_tol', rec.Qc_tol, ...
    'TN_min', rec.TN_min, 'TN_mean', rec.TN_mean, 'TN_maxmin', rec.TN_maxmin, ...
    'newton_iters', rec.newton_iters, 'newton_relaxed', logical(rec.newton_relaxed), ...
    'newton_rel_max', rec.newton_rel_max, ...
    'newton_mode', char(string(rec.newton_convergence_mode)), ...
    'Lx_mm', vec_to_inline_str_ms(rec.Lx * 1e3), ...
    'Ly_mm', vec_to_inline_str_ms(rec.Ly * 1e3), ...
    'coverage', vec_to_inline_str_ms(rec.cov), ...
    'elapsed_sec', elapsed_sec);
end

function caxis_by_stage = make_caxis_by_stage_from_eval_ms(rec, N)
caxis_by_stage = NaN(N, 2);
if ~isfield(rec, 'Tfields') || isempty(rec.Tfields)
    return;
end
for k = 1:N
    if numel(rec.Tfields) >= k && ~isempty(rec.Tfields{k})
        tk = rec.Tfields{k};
        caxis_by_stage(k,:) = [min(tk(:)), max(tk(:))];
    end
end
end

function tf = has_temperature_fields_ms(rec, N)
tf = false;
if nargin < 2 || ~isfinite(N)
    N = 1;
end
if ~isstruct(rec) || ~isfield(rec, 'Tfields') || numel(rec.Tfields) < N
    return;
end
for k = 1:N
    if isempty(rec.Tfields{k}) || ~all(isfinite(rec.Tfields{k}(:)))
        return;
    end
end
tf = true;
end

function write_text_file_ms(path_out, txt)
fid = fopen(path_out, 'w');
if fid < 0
    return;
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', txt);
end

function write_scalar_struct_csv_ms(row, out_csv)
if nargin < 2 || isempty(out_csv) || ~isstruct(row)
    return;
end
row1 = row(1);
fn = fieldnames(row1);
for i = 1:numel(fn)
    v = row1.(fn{i});
    if isnumeric(v) || islogical(v)
        if isempty(v)
            row1.(fn{i}) = '';
        elseif ~isscalar(v)
            row1.(fn{i}) = vec_to_inline_str_ms(v);
        end
    elseif isstring(v)
        if isscalar(v)
            row1.(fn{i}) = char(v);
        else
            row1.(fn{i}) = vec_to_inline_str_ms(cellstr(v(:).'));
        end
    elseif iscell(v)
        row1.(fn{i}) = vec_to_inline_str_ms(v);
    end
end
try
    T = struct2table(row1, 'AsArray', true);
catch
    vals = cell(1, numel(fn));
    for i = 1:numel(fn)
        vals{i} = {scalar_csv_value_ms(row1.(fn{i}))};
    end
    T = table(vals{:}, 'VariableNames', fn);
end
writetable(T, out_csv);
end

function v = scalar_csv_value_ms(v_in)
if isnumeric(v_in) || islogical(v_in)
    if isempty(v_in)
        v = '';
    elseif isscalar(v_in)
        v = sprintf('%.15g', double(v_in));
    else
        v = vec_to_inline_str_ms(v_in);
    end
elseif isstring(v_in)
    if isempty(v_in)
        v = '';
    elseif isscalar(v_in)
        v = char(v_in);
    else
        v = vec_to_inline_str_ms(cellstr(v_in(:).'));
    end
elseif iscell(v_in)
    v = vec_to_inline_str_ms(v_in);
elseif ischar(v_in)
    v = v_in;
else
    v = char(string(v_in));
end
end

function metric_results = compute_metric_postprocess_ms(best_rec, best_cand, ratio_rec, ratio_cand, G, spec)
metric_results = struct();
metric_results.summary_rows = repmat(empty_metric_summary_row_ms(), 0, 1);
metric_results.qmax_rows = repmat(empty_qmax_scan_row_ms(), 0, 1);
if ~isstruct(spec) || ~isfield(spec, 'metrics') || ~isstruct(spec.metrics)
    return;
end
labels = {'best', 'ratio_reference'};
recs = {best_rec, ratio_rec};
for i = 1:numel(labels)
    rec_i = recs{i};
    row = empty_metric_summary_row_ms();
    row.label = labels{i};
    if isstruct(rec_i) && isfield(rec_i, 'n')
        row.n = vec_to_inline_str_ms(rec_i.n);
        row.total_n = sum(rec_i.n);
        row.DeltaT_at_Qc_K = rec_i.DeltaTN_actual;
        row.I_at_DeltaT_at_Qc_A = rec_i.I_opt;
        row.Qc_at_DeltaT_at_Qc_W = safe_numeric_field_ms(rec_i, 'Qc_last_total', NaN);
        row.has_temperature_reversal = logical(safe_numeric_field_ms(rec_i, 'has_temperature_reversal', false));
    else
        row.message = 'result unavailable; metric post-processing skipped';
    end
    metric_results.summary_rows(end+1,1) = row; %#ok<AGROW>
end
end

function write_metric_outputs_ms(metric_results, out_dir)
if nargin < 2 || isempty(out_dir) || ~isstruct(metric_results)
    return;
end
if isfield(metric_results, 'summary_rows') && ~isempty(metric_results.summary_rows)
    write_struct_rows_csv_ms(metric_results.summary_rows, fullfile(out_dir, 'metric_summary.csv'));
end
if isfield(metric_results, 'qmax_rows') && ~isempty(metric_results.qmax_rows)
    write_struct_rows_csv_ms(metric_results.qmax_rows, fullfile(out_dir, 'qmax_scan.csv'));
end
end

function row = empty_metric_summary_row_ms()
row = struct('label', '', 'n', '[]', 'total_n', NaN, ...
    'DeltaT_at_Qc_K', NaN, 'I_at_DeltaT_at_Qc_A', NaN, 'Qc_at_DeltaT_at_Qc_W', NaN, ...
    'DeltaTmax_K', NaN, 'I_at_DeltaTmax_A', NaN, ...
    'Qmax_W', NaN, 'I_at_Qmax_A', NaN, 'DeltaT_at_Qmax_K', NaN, ...
    'Qmax_success', false, 'has_temperature_reversal', false, ...
    'has_temperature_reversal_DeltaTmax', false, 'has_temperature_reversal_Qmax', false, ...
    'message', '');
end

function row = empty_qmax_scan_row_ms()
row = struct('label', '', 'phase', '', 'I_A', NaN, 'Qc_W', NaN, ...
    'success', false, 'feasible', false, 'DeltaTN_actual', NaN, ...
    'Qc_last_total', NaN, 'Qc_error', NaN, 'has_temperature_reversal', false, ...
    'message', '');
end

function write_k25_debug_result_mat_ms(out_mat, spec, cand, rec, row, scan_rows, particle_search_rows, ratio_ref_rec, ratio_ref_cand, top_n_current, metric_results)
if nargin < 7 || isempty(out_mat)
    return;
end
if nargin < 10
    top_n_current = repmat(empty_top_n_current_row_ms(spec.stage_count), 0, 1);
end
if nargin < 11
    metric_results = struct();
end
try
    save(out_mat, 'spec', 'cand', 'rec', 'row', 'scan_rows', 'particle_search_rows', ...
        'ratio_ref_rec', 'ratio_ref_cand', 'top_n_current', 'metric_results', '-v7.3');
catch ME_v73
    try
        save(out_mat, 'spec', 'cand', 'rec', 'row', 'scan_rows', 'particle_search_rows', ...
            'ratio_ref_rec', 'ratio_ref_cand', 'top_n_current', 'metric_results');
    catch ME_save
        warning('K25 debug MAT export failed: %s; fallback failed: %s', ME_v73.message, ME_save.message);
    end
end
end

function write_k25_fast_summary_ms(out_txt, spec, spec_source, cand, rec, row, t0, metric_results)
if nargin < 8
    metric_results = struct();
end
fid = fopen(out_txt, 'w');
if fid < 0
    return;
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'k25 fast fixed-layout FEM debug summary\n');
fprintf(fid, 'spec_source=%s\n', spec_source);
fprintf(fid, 'seed=%d\n', spec.seed);
fprintf(fid, 'stage_count=%d\n', spec.stage_count);
fprintf(fid, 'total_n=%d\n', row.total_n);
fprintf(fid, 'n=%s\n', row.n);
fprintf(fid, 'npair=%s\n', row.npair);
fprintf(fid, 'ratio_cold_to_hot=1:3:7:21:55\n');
fprintf(fid, 'layout_method=%s\n', cand.layout_method);
fprintf(fid, 'spacing_factor=0.5\n');
fprintf(fid, 'edge_margin_factor=0.5\n');
fprintf(fid, 'plate_k_inplane=%.6g\n', spec.plate_k_inplane);
fprintf(fid, 'ceramic_plate_k_inplane_W_mK=%.12g\n', spec.plate_k_inplane);
fprintf(fid, 'particle_count_n=%s\n', row.n);
fprintf(fid, 'particle_count_total=%d\n', row.total_n);
fprintf(fid, 'particle_width_mm=%.12g\n', spec.fp_w * 1e3);
fprintf(fid, 'particle_length_mm=%.12g\n', spec.fp_h * 1e3);
fprintf(fid, 'particle_height_stage_mm=%s\n', vec_to_inline_str_ms(spec.fp_t_stage * 1e3));
fprintf(fid, 'fp_w_mm=%.12g\n', spec.fp_w * 1e3);
fprintf(fid, 'fp_h_mm=%.12g\n', spec.fp_h * 1e3);
fprintf(fid, 'fp_t_stage_mm=%s\n', vec_to_inline_str_ms(spec.fp_t_stage * 1e3));
fprintf(fid, 'particle_area_mm2=%.12g\n', spec.fp_w * spec.fp_h * 1e6);
fprintf(fid, 'particle_volume_stage_mm3=%s\n', vec_to_inline_str_ms(spec.fp_w * spec.fp_h * spec.fp_t_stage * 1e9));
fprintf(fid, 'current_I_init_A=%.12g\n', spec.current.I_init);
fprintf(fid, 'I_selected=%.6f A\n', row.I_A);
fprintf(fid, 'current_I_selected_A=%.12g\n', row.I_A);
if has_explicit_k25_n_ms(spec)
    fprintf(fid, 'particle_mode=fixed_n\n');
else
    fprintf(fid, 'particle_mode=min_total_search\n');
    fprintf(fid, 'particle_search_only=%d\n', logical(spec.k25_debug.search_only));
    fprintf(fid, 'particle_search_algorithm=template_0d_pruned\n');
	    fprintf(fid, 'n_current_scan_select_rule=%s\n', spec.n_current_scan.select_rule);
	    fprintf(fid, 'n_current_scan_topK=%d\n', get_n_current_scan_topK_ms(spec));
        fprintf(fid, 'selection_basis=0D\n');
        fprintf(fid, 'fem_is_acceptance_criterion=0\n');
        if isfield(rec, 'selected_0d_target_met')
            fprintf(fid, 'selected_0d_target_met=%d\n', logical(rec.selected_0d_target_met));
        end
        if isfield(rec, 'selected_0d_I_opt')
            fprintf(fid, 'selected_0d_I_opt=%.9g\n', rec.selected_0d_I_opt);
        end
        if isfield(rec, 'selected_0d_DeltaTN_actual')
            fprintf(fid, 'selected_0d_DeltaTN_actual=%.9g\n', rec.selected_0d_DeltaTN_actual);
        end
        if isfield(rec, 'selected_0d_Qc_last_total')
            fprintf(fid, 'selected_0d_Qc_last_total=%.9g\n', rec.selected_0d_Qc_last_total);
        end
        if isfield(rec, 'selected_0d_Qc_error')
            fprintf(fid, 'selected_0d_Qc_error=%.9g\n', rec.selected_0d_Qc_error);
        end
	    fprintf(fid, 'particle_search_stop_rule=min_total_topK\n');
    if isfield(rec, 'particle_search_stop_reason')
        fprintf(fid, 'particle_search_stop_reason=%s\n', rec.particle_search_stop_reason);
    end
    if isfield(rec, 'particle_search_stop_total_n')
        fprintf(fid, 'particle_search_stop_total_n=%.9g\n', rec.particle_search_stop_total_n);
    end
    if isfield(rec, 'particle_search_topK_found')
        fprintf(fid, 'particle_search_topK_found=%d\n', rec.particle_search_topK_found);
    end
    fprintf(fid, 'particle_search_mode=%s\n', spec.k25_debug.particle_search.mode);
    fprintf(fid, 'particle_search_strict_min_total=%d\n', logical(spec.k25_debug.particle_search.strict_min_total));
    fprintf(fid, 'particle_search_min_total_grid=template_on_total_step\n');
    if isfield(rec, 'particle_search_n5_estimate')
        fprintf(fid, 'particle_search_n5_estimate=%d\n', rec.particle_search_n5_estimate);
    end
    fprintf(fid, 'particle_search_n5_min_factor=%.9g\n', spec.k25_debug.particle_search.n5_min_factor);
    if isfield(rec, 'particle_search_n5_min')
        fprintf(fid, 'particle_search_n5_min=%d\n', rec.particle_search_n5_min);
    end
    if isfield(rec, 'particle_search_n5_max_effective')
        fprintf(fid, 'particle_search_n5_max_effective=%d\n', rec.particle_search_n5_max_effective);
    end
    if isempty(spec.k25_debug.particle_search.n5_max)
        fprintf(fid, 'particle_search_n5_max=user_empty_auto\n');
    else
        fprintf(fid, 'particle_search_n5_max=%d\n', spec.k25_debug.particle_search.n5_max);
    end
    fprintf(fid, 'particle_search_stage_max_count=%d\n', spec.k25_debug.particle_search.stage_max_count);
    if isfield(spec.k25_debug.particle_search, 'total_n_min')
        fprintf(fid, 'particle_search_total_n_min=%d\n', spec.k25_debug.particle_search.total_n_min);
    end
    fprintf(fid, 'particle_search_total_n_max=%d\n', spec.k25_debug.particle_search.total_n_max);
    fprintf(fid, 'particle_search_total_step=%d\n', spec.k25_debug.particle_search.total_step);
    fprintf(fid, 'particle_search_total_block_size=%d\n', spec.k25_debug.particle_search.total_block_size);
    fprintf(fid, 'particle_search_par_chunk_size=%d\n', spec.k25_debug.particle_search.par_chunk_size);
    fprintf(fid, 'particle_search_use_ratio_prune=%d\n', logical(spec.k25_debug.particle_search.use_ratio_prune));
    fprintf(fid, 'particle_search_adjacent_min_fraction=%.9g\n', spec.k25_debug.particle_search.adjacent_min_fraction);
    fprintf(fid, 'particle_search_candidate_budget=%d\n', spec.k25_debug.particle_search.candidate_budget);
    fprintf(fid, 'particle_search_boundary_check_count=%d\n', spec.k25_debug.particle_search.boundary_check_count);
    fprintf(fid, 'particle_search_seed_ratios_hot_to_cold=%s\n', ...
        vec_to_inline_str_ms(spec.k25_debug.particle_search.seed_ratios));
    if isfinite(row.total_n)
        fprintf(fid, 'particle_search_template_failed_total_n_below=%d\n', max(0, row.total_n - spec.k25_debug.particle_search.total_step));
    end
end
if isfield(spec, 'k25_debug') && isfield(spec.k25_debug, 'physics') && ...
        isfield(spec.k25_debug.physics, 'Qc_pair_negative_tol')
    fprintf(fid, 'Qc_pair_negative_tol=%.9g\n', spec.k25_debug.physics.Qc_pair_negative_tol);
end
fprintf(fid, 'mesh_nx_stage_full=%s\n', vec_to_inline_str_ms(spec.mesh_nx_stage_full));
fprintf(fid, 'mesh_ny_stage_full=%s\n', vec_to_inline_str_ms(spec.mesh_ny_stage_full));
fprintf(fid, 'parallel_use=%d\n', logical(spec.use_parallel));
fprintf(fid, 'candidate_Lx_mm=%s\n', row.Lx_mm);
fprintf(fid, 'candidate_Ly_mm=%s\n', row.Ly_mm);
fprintf(fid, 'coverage=%s\n', row.coverage);
fprintf(fid, 'success=%d\n', row.success);
if ~has_explicit_k25_n_ms(spec)
    target_met_acceptance = false;
    if isfield(rec, 'selected_0d_target_met')
        target_met_acceptance = logical(rec.selected_0d_target_met);
    end
    fprintf(fid, 'target_met=%d\n', target_met_acceptance);
    fprintf(fid, 'acceptance_target_met=%d\n', target_met_acceptance);
    fprintf(fid, 'fem_target_met_diagnostic=%d\n', row.target_met);
else
    fprintf(fid, 'target_met=%d\n', row.target_met);
end
fprintf(fid, 'message=%s\n', row.message);
if isfield(rec, 'newton_failure_detail') && ~isempty(rec.newton_failure_detail)
    fprintf(fid, 'newton_failure_detail=%s\n', rec.newton_failure_detail);
end
fprintf(fid, 'physics_valid=%d\n', logical(row.physics_valid));
fprintf(fid, 'physics_failure_reason=%s\n', row.physics_failure_reason);
fprintf(fid, 'stage_Tc_mean=%s\n', row.stage_Tc_mean);
fprintf(fid, 'stage_Th_mean=%s\n', row.stage_Th_mean);
fprintf(fid, 'stage_min_Th_minus_Tc=%s\n', row.stage_min_Th_minus_Tc);
fprintf(fid, 'has_temperature_reversal=%d\n', logical(row.has_temperature_reversal));
fprintf(fid, 'sumQc_stage=%s\n', row.sumQc_stage);
fprintf(fid, 'sumQh_stage=%s\n', row.sumQh_stage);
fprintf(fid, 'Qc_pair_min=%s\n', row.Qc_pair_min);
fprintf(fid, 'Qh_pair_min=%s\n', row.Qh_pair_min);
fprintf(fid, 'DeltaTN_actual=%.9g\n', row.DeltaTN_actual);
fprintf(fid, 'DeltaTN_mean=%.9g\n', row.DeltaTN_mean);
fprintf(fid, 'Qc_last_total=%.9g\n', row.Qc_last_total);
fprintf(fid, 'Qc_target_last=%.9g\n', row.Qc_target_last);
fprintf(fid, 'Qc_error=%.9g\n', row.Qc_error);
fprintf(fid, 'Qc_tol=%.9g\n', row.Qc_tol);
fprintf(fid, 'newton_iters=%.9g\n', row.newton_iters);
fprintf(fid, 'newton_relaxed=%d\n', logical(row.newton_relaxed));
fprintf(fid, 'newton_rel_max=%.9g\n', row.newton_rel_max);
fprintf(fid, 'newton_mode=%s\n', row.newton_mode);
fprintf(fid, 'fem_elapsed_sec=%.6f\n', row.elapsed_sec);
if isstruct(metric_results) && isfield(metric_results, 'summary_rows') && ~isempty(metric_results.summary_rows)
    fprintf(fid, 'metric_summary_csv=metric_summary.csv\n');
    for mi = 1:numel(metric_results.summary_rows)
        mr = metric_results.summary_rows(mi);
        fprintf(fid, 'metric_%d_label=%s, Qmax_W=%.9g, I_at_Qmax_A=%.9g, DeltaTmax_K=%.9g, DeltaT_at_Qc_K=%.9g\n', ...
            mi, char(string(mr.label)), mr.Qmax_W, mr.I_at_Qmax_A, mr.DeltaTmax_K, mr.DeltaT_at_Qc_K);
    end
end
fprintf(fid, 'total_runtime_sec=%.6f\n', toc(t0));
if ~rec.success
    fprintf(fid, 'failure_note=FEM failed or was rejected before accepted temperature output.\n');
end
end

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

function spec = apply_metric_mode_targets_ms(spec)
if ~isstruct(spec) || ~isfield(spec, 'metrics') || ~isstruct(spec.metrics) || ...
        ~isfield(spec.metrics, 'mode')
    return;
end
mode = lower(strtrim(char(string(spec.metrics.mode))));
switch mode
    case 'deltatmax'
        spec.targets.Qc_target_last = 0;
    case 'qmax'
        spec.targets.DeltaT_target = 0;
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

function spec = ensure_k25_parallel_pool_ms(spec)
if ~isstruct(spec) || ~isfield(spec, 'parallel') || ~isstruct(spec.parallel)
    spec.use_parallel = false;
    return;
end

use_requested = logical(spec.use_parallel);
if ~use_requested
    spec.use_parallel = false;
    spec.parallel.pool_used = false;
    spec.parallel.pool_workers_actual = 0;
    return;
end

has_parallel = ((exist('parpool', 'file') == 2) || (exist('parpool', 'builtin') == 5)) && ...
    ((exist('gcp', 'file') == 2) || (exist('gcp', 'builtin') == 5));
if ~has_parallel
    warning('K25 particle search parallel requested but Parallel Computing Toolbox functions are unavailable; falling back to serial.');
    spec.use_parallel = false;
    spec.parallel.pool_used = false;
    spec.parallel.pool_workers_actual = 0;
    return;
end

worker_target = max(1, round(scalar_with_default_ms(spec.parallel.pool_workers, 64)));
try
    cluster = parcluster('local');
    if isprop(cluster, 'NumWorkers')
        worker_target = min(worker_target, max(1, cluster.NumWorkers));
    end
catch
end

try
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('local', worker_target);
    elseif pool.NumWorkers ~= worker_target
        fprintf('[K25-Debug] existing parallel pool workers=%d, requested=%d; using existing pool.\n', ...
            pool.NumWorkers, worker_target);
    end
    if isempty(pool)
        spec.use_parallel = false;
        spec.parallel.pool_used = false;
        spec.parallel.pool_workers_actual = 0;
    else
        spec.use_parallel = true;
        spec.parallel.pool_used = true;
        spec.parallel.pool_workers_actual = pool.NumWorkers;
    end
catch ME_pool
    warning('K25 particle search parallel pool unavailable: %s. Falling back to serial.', ME_pool.message);
    spec.use_parallel = false;
    spec.parallel.pool_used = false;
    spec.parallel.pool_workers_actual = 0;
end
end

function nw = get_parallel_worker_count_ms()
nw = 0;
if ~((exist('gcp', 'file') == 2) || (exist('gcp', 'builtin') == 5))
    return;
end
try
    pool = gcp('nocreate');
    if ~isempty(pool)
        nw = max(0, pool.NumWorkers);
    end
catch
    nw = 0;
end
end

function spec = apply_default_spec_ms(spec)
[spec, ~] = optimize_layout_multistage0411_shared_params(spec, 'n');
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
i_tag = numeric_tag_ms(spec.current.I_init);
k_tag = numeric_tag_ms(spec.plate_k_inplane);
run_name = sprintf('n_search_N%d_I%s_k%s_DeltaT%s_Qc%s_%s', ...
    spec.stage_count, i_tag, k_tag, dt_tag, qc_tag, stamp);
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

function print_runtime_spec_summary_ms(spec, spec_source)
fprintf('Config source: %s\n', spec_source);
fprintf('stage_count=%d\n', spec.stage_count);
fprintf('Targets: DeltaT_target=%.3f K, Qc_target_last=%.3f W, ThRes=%.3f K\n', ...
    spec.targets.DeltaT_target, spec.targets.Qc_target_last, spec.targets.ThRes);
fprintf('L_max_mm: %s\n', vec_to_inline_str_ms(spec.geometry.L_max_mm));
fprintf('coverage_min: %s\n', vec_to_inline_str_ms(spec.geometry.coverage_min));
fprintf('pyramid_gap_min_mm: %s\n', vec_to_inline_str_ms(spec.geometry.pyramid_gap_min_mm));
fprintf('min_edge_gap_mm: %.3f\n', spec.geometry.min_edge_gap_mm);
if has_explicit_k25_n_ms(spec)
    fprintf('Particle mode: fixed n=%s, fixed_I=%.3f A, plate_k=%.3g W/mK\n', ...
        vec_to_inline_str_ms(spec.k25_debug.n), spec.current.I_init, spec.k25_debug.plate_k_inplane);
else
    fprintf('Particle mode: minimum total search, fixed_I=%.3f A, plate_k=%.3g W/mK\n', ...
        spec.current.I_init, spec.k25_debug.plate_k_inplane);
    fprintf('Particle search_only: %d\n', spec.k25_debug.search_only);
    if isempty(spec.k25_debug.particle_search.n5_max)
        n5_max_txt = 'auto';
    else
        n5_max_txt = sprintf('%d', spec.k25_debug.particle_search.n5_max);
    end
    if isfield(spec.k25_debug.particle_search, 'total_n_min')
        total_n_min_txt = sprintf('%d', spec.k25_debug.particle_search.total_n_min);
    else
        total_n_min_txt = 'auto';
    end
    fprintf(['Particle search: mode=%s, strict=%d, n5_min_factor=%.4g, n5_max=%s, stage_max_count=%d, total_n_min=%s, total_n_max=%d, ' ...
        'total_step=%d, block=%d, chunk=%d, ratio_prune=%d, adjacent_min_fraction=%.4g, seed_ratios=%s\n'], ...
        spec.k25_debug.particle_search.mode, spec.k25_debug.particle_search.strict_min_total, ...
        spec.k25_debug.particle_search.n5_min_factor, n5_max_txt, spec.k25_debug.particle_search.stage_max_count, ...
        total_n_min_txt, spec.k25_debug.particle_search.total_n_max, spec.k25_debug.particle_search.total_step, ...
        spec.k25_debug.particle_search.total_block_size, spec.k25_debug.particle_search.par_chunk_size, ...
        spec.k25_debug.particle_search.use_ratio_prune, spec.k25_debug.particle_search.adjacent_min_fraction, ...
        vec_to_inline_str_ms(spec.k25_debug.particle_search.seed_ratios));
end
fprintf('Physics tolerances: Qc_pair_negative_tol=%.6g W\n', ...
    spec.k25_debug.physics.Qc_pair_negative_tol);
fprintf('Low-k solver: enable=%d, threshold=%.3g, omega=%.3g, max_inner=%d, max_outer=%d, fracs=%s\n', ...
    spec.k25_debug.low_k_solver.enable, spec.k25_debug.low_k_solver.k_threshold, ...
    spec.k25_debug.low_k_solver.omega, spec.k25_debug.low_k_solver.max_inner, ...
    spec.k25_debug.low_k_solver.max_outer, vec_to_inline_str_ms(spec.k25_debug.low_k_solver.continuation_fracs));
fprintf('Plot: use_global_axis=%d, show_mesh_edges=%d, view_mode=%s, axis_margin=%.3f\n', ...
    spec.output.plot.use_global_axis, spec.output.plot.show_mesh_edges, ...
    spec.output.plot.view_mode, spec.output.plot.axis_margin_ratio);
fprintf('Plot extras: save_overview=%d, save_stage_separate=%d, annotate_substrate_dims=%d, separate_fig_size_px=%s\n', ...
    spec.output.plot.save_overview, spec.output.plot.save_stage_separate, ...
    spec.output.plot.annotate_substrate_dims, vec_to_inline_str_ms(spec.output.plot.separate_fig_size_px));
fprintf('Full FEM mesh (stage-wise): nx=%s, ny=%s\n', ...
    vec_to_inline_str_ms(spec.mesh_nx_stage_full), vec_to_inline_str_ms(spec.mesh_ny_stage_full));
fprintf('Parallel: use=%d, pool_workers=%d\n', spec.use_parallel, spec.parallel.pool_workers);
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

function fracs = normalize_continuation_fracs_ms(v_in, fallback)
if nargin < 2 || isempty(fallback)
    fallback = [0.75, 0.875, 0.95, 1.0];
end
if isempty(v_in)
    v_in = fallback;
end
fracs = double(v_in(:).');
fracs = fracs(isfinite(fracs) & fracs > 0 & fracs <= 1);
if isempty(fracs)
    fracs = double(fallback(:).');
end
fracs = unique(sort([fracs, 1.0]));
end

function G = apply_low_k_solver_stabilization_ms(G, spec)
if nargin < 2 || ~isstruct(spec) || ~isfield(spec, 'k25_debug') || ...
        ~isstruct(spec.k25_debug) || ~isfield(spec.k25_debug, 'low_k_solver') || ...
        ~isstruct(spec.k25_debug.low_k_solver)
    return;
end
lks = spec.k25_debug.low_k_solver;
if ~isfield(lks, 'enable') || ~logical(lks.enable) || ...
        ~isfield(spec, 'plate_k_inplane') || ~isfinite(spec.plate_k_inplane) || ...
        spec.plate_k_inplane > lks.k_threshold
    return;
end
G.omega = min(G.omega, lks.omega);
G.max_inner = max(G.max_inner, lks.max_inner);
G.max_outer = max(G.max_outer, lks.max_outer);
end

function fracs = get_current_continuation_fracs_ms(spec)
fracs = [0.75, 0.875, 0.95, 1.0];
if nargin < 1 || ~isstruct(spec) || ~isfield(spec, 'k25_debug') || ...
        ~isstruct(spec.k25_debug) || ~isfield(spec.k25_debug, 'low_k_solver') || ...
        ~isstruct(spec.k25_debug.low_k_solver)
    return;
end
lks = spec.k25_debug.low_k_solver;
if isfield(lks, 'enable') && logical(lks.enable) && ...
        isfield(spec, 'plate_k_inplane') && isfinite(spec.plate_k_inplane) && ...
        spec.plate_k_inplane <= lks.k_threshold && ...
        isfield(lks, 'continuation_fracs')
    fracs = normalize_continuation_fracs_ms(lks.continuation_fracs, fracs);
end
end

function G = init_G_params_ms(spec)
[~, G] = optimize_layout_multistage0411_shared_params(spec, 'n');
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
[ok_pair, ~] = particle_counts_to_pair_counts_ms(n);
if ~ok_pair
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
[ok_pair, ~] = particle_counts_to_pair_counts_ms(n);
if ~ok_pair
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
Qc_pair = zeros(N,1);
Qh_pair = zeros(N,1);
for k = 1:N
    Tc = x(k);
    if k == 1
        Th = spec.targets.ThRes;
    else
        Th = x(k-1);
    end
    Qc_pair(k) = te_Qc_onecouple_ms(Tc, Th, I, G, k);
    Qh_pair(k) = te_Qh_onecouple_ms(Th, Tc, I, G, k);
end
if ~all(isfinite([Qc_pair; Qh_pair]))
    r = ones(N,1) * 1e12;
    return;
end

r = zeros(N,1);
for k = 1:N-1
    r(k) = npair(k+1) * Qh_pair(k+1) - npair(k) * Qc_pair(k);
end
r(N) = npair(N) * Qc_pair(N) - spec.targets.Qc_target_last;
end

function rec = candidate_record_template_ms(N)
if nargin < 1
    N = 5;
end
rec = struct('candidate_id', -1, 'layout_method', '', ...
    'stage_modes', {repmat({''}, 1, N)}, ...
    'stage_methods', {repmat({''}, 1, N)}, ...
    'stage_trends', {repmat({'neutral'}, 1, N)}, ...
    'symmetry_mode', '', 'edge_pattern_mode', '', ...
    's_dense', NaN, 's_sparse', NaN, 'expo', NaN, 'anis_ratio', NaN, 'gamma', NaN, ...
    'method_anchor_stage', NaN, ...
    'spacing_ratio', NaN, 'contrast_score', NaN, 'mode_prior_score', NaN, ...
    'lmax_relaxed', false, 'lmax_relax_ratio', 1.0, ...
    'n', NaN(1, N), 'ratios', NaN(1, N-1), ...
    'Lx', NaN(1, N), 'Ly', NaN(1, N), 'cov', NaN(1, N), ...
    'stage_rects', {cell(1, N)});
end

function centers = rect_centers_ms(rects)
centers = [0.5*(rects(:,1)+rects(:,2)), 0.5*(rects(:,3)+rects(:,4))];
end

function rects = centers_to_rects_ms(points, w, h)
rects = [points(:,1)-w/2, points(:,1)+w/2, points(:,2)-h/2, points(:,2)+h/2];
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
rec.spacing_ratio = cand.spacing_ratio;
rec.contrast_score = cand.contrast_score;
rec.mode_prior_score = NaN;
rec.lmax_relaxed = cand.lmax_relaxed;
rec.lmax_relax_ratio = cand.lmax_relax_ratio;
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
rec.I_opt = G.I;

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
    [ok0d, C] = solve_0d_at_current_ms(cand.n, G.I, G, spec, x0);
    if ~ok0d
        rec.message = '0D solve failed';
        return;
    end
    C_init = C(:).';

    G_eval = G;
    G_eval.vertical_runtime = struct('enable', false, ...
        'Rz_interfaces', zeros(1, max(1, N-1)), 'Rz_sink', 0, ...
        'step_fp_iters', 1, 'step_fp_relax', 1.0, 'step_q_prev_weight', 0);
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
        [best_cont, cont_msg, cont_attempts] = newton_solve_with_continuation_ms(G_eval, plates, cand.n, C_init, N, spec);
        if best_cont.success
            best = best_cont;
            used_continuation = true;
        else
            rec = copy_newton_diag_to_eval_ms(rec, best_cont, 'continuation_failed');
            rec.newton_failure_detail = sprintf('%s | continuation_status=%s | %s', ...
                format_newton_attempt_detail_ms('direct', G_eval.I, best), ...
                cont_msg, format_newton_attempts_ms(cont_attempts));
            rec.message = sprintf('newton not converged (direct+continuation failed: %s)', cont_msg);
            return;
        end
    end
    if ~best.success
        rec = copy_newton_diag_to_eval_ms(rec, best, 'direct_failed');
        rec.newton_failure_detail = format_newton_attempt_detail_ms('direct', G_eval.I, best);
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
    if isfield(best, 'cache') && isstruct(best.cache)
        if isfield(best.cache, 'sumQc'), rec.sumQc_stage = best.cache.sumQc(:).'; end
        if isfield(best.cache, 'sumQh'), rec.sumQh_stage = best.cache.sumQh(:).'; end
        if isfield(best.cache, 'Qc_pair_min'), rec.Qc_pair_min = best.cache.Qc_pair_min(:).'; end
        if isfield(best.cache, 'Qh_pair_min'), rec.Qh_pair_min = best.cache.Qh_pair_min(:).'; end
    end
    [physics_ok, physics_diag] = validate_fem_physics_ms(fields, plates, cand.n, G_eval, spec);
    rec.physics_valid = physics_ok;
    rec.physics_failure_reason = physics_diag.reason;
    rec.stage_Tc_mean = physics_diag.stage_Tc_mean;
    rec.stage_Th_mean = physics_diag.stage_Th_mean;
    rec.stage_min_Th_minus_Tc = physics_diag.stage_min_Th_minus_Tc;
    rec.has_temperature_reversal = logical(physics_diag.has_temperature_reversal);
    rec.sumQc_stage = physics_diag.sumQc_stage;
    rec.sumQh_stage = physics_diag.sumQh_stage;
    rec.Qc_pair_min = physics_diag.Qc_pair_min;
    rec.Qh_pair_min = physics_diag.Qh_pair_min;
    if numel(rec.sumQc_stage) >= N
        rec.Qc_last_total = rec.sumQc_stage(N);
    end
    rec.Qc_target_last = spec.targets.Qc_target_last;
    rec.Qc_error = rec.Qc_last_total - spec.targets.Qc_target_last;
    rec.Qc_tol = get_n_current_scan_qc_tol_ms(spec);
    rec.target_met = rec.DeltaTN_actual >= spec.targets.DeltaT_target && ...
        isfinite(rec.Qc_error) && abs(rec.Qc_error) <= rec.Qc_tol;
    if keep_fields
        rec.plates = plates;
        rec.Tfields = fields;
    end
    if ~physics_ok
        rec.success = false;
        rec.message = sprintf('FEM physics validation failed: %s', physics_diag.reason);
        physics_txt = format_physics_diag_ms(physics_diag);
        if isempty(rec.newton_failure_detail)
            rec.newton_failure_detail = physics_txt;
        else
            rec.newton_failure_detail = sprintf('%s | %s', rec.newton_failure_detail, physics_txt);
        end
        rec.rank_score = -inf;
        return;
    end
    rec.success = true;
    rec.message = 'ok';
    rec.rank_score = rec.DeltaTN_actual;

    if keep_fields
        rec.plates = plates;
        rec.Tfields = fields;
    end
catch ME
    rec.success = false;
    rec.message = ME.message;
end
end

function [ok, npair] = particle_counts_to_pair_counts_ms(n_vec)
n_vec = n_vec(:).';
npair = NaN(size(n_vec));
ok = ~isempty(n_vec) && all(isfinite(n_vec)) && all(n_vec >= 2) && all(abs(n_vec - round(n_vec)) < 1e-9) && ...
    all(mod(round(n_vec), 2) == 0);
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

function v_pair = pair_average_values_ms(v_particle, pairs)
if isempty(pairs)
    v_pair = zeros(0, 1);
    return;
end
v_particle = v_particle(:);
v_pair = 0.5 * (v_particle(pairs(:,1)) + v_particle(pairs(:,2)));
end

function F = add_heat_to_footprint_elems_ms(F, fpElems, elem, tri, fp_idx, q_flux)
elist = fpElems{fp_idx};
for ee = elist
    nodes = tri(ee,:);
    F(nodes) = F(nodes) + q_flux * elem.A(ee) / 3;
end
end

function [ok, diag] = validate_fem_physics_ms(fields, plates, n_vec, G, spec)
N = numel(n_vec);
ok = false;
diag = struct('reason', '', 'stage_Tc_mean', NaN(1,N), ...
    'stage_Th_mean', NaN(1,N), 'stage_min_Th_minus_Tc', NaN(1,N), ...
    'has_temperature_reversal', false, ...
    'sumQc_stage', NaN(1,N), 'sumQh_stage', NaN(1,N), ...
    'Qc_pair_min', NaN(1,N), 'Qh_pair_min', NaN(1,N));
if numel(fields) < N || numel(plates) < N
    diag.reason = 'missing_stage_fields';
    return;
end
tol_neg_qc = get_qc_pair_negative_tol_ms(spec);
[ok_pair, npair] = particle_counts_to_pair_counts_ms(n_vec);
if ~ok_pair
    diag.reason = 'invalid_even_particle_counts';
    return;
end
for k = 1:N
    if isempty(fields{k}) || ~all(isfinite(fields{k}(:)))
        diag.reason = sprintf('stage%d_nonfinite_temperature_field', k);
        return;
    end
    if ~isfield(plates{k}, 'fpElemsCool') || isempty(plates{k}.fpElemsCool)
        diag.reason = sprintf('stage%d_missing_cool_footprints', k);
        return;
    end
    Tc_plate = footprint_avg_ms(fields{k}, plates{k}.fpElemsCool, plates{k}.elem, plates{k}.mesh.tri);
    if k == 1
        Th_plate = G.ThRes + zeros(size(Tc_plate));
    else
        if ~isfield(plates{k-1}, 'fpElemsHot') || isempty(plates{k-1}.fpElemsHot)
            diag.reason = sprintf('stage%d_missing_hot_footprints', k);
            return;
        end
        Th_plate = footprint_avg_ms(fields{k-1}, plates{k-1}.fpElemsHot, ...
            plates{k-1}.elem, plates{k-1}.mesh.tri);
    end
    if numel(Tc_plate) ~= n_vec(k) || numel(Th_plate) ~= n_vec(k)
        diag.reason = sprintf('stage%d_footprint_count_mismatch', k);
        return;
    end
    if any(~isfinite(Tc_plate)) || any(~isfinite(Th_plate))
        diag.reason = sprintf('stage%d_nonfinite_footprint_temperature', k);
        return;
    end
    margin = Th_plate(:) - Tc_plate(:);
    diag.stage_Tc_mean(k) = mean(Tc_plate);
    diag.stage_Th_mean(k) = mean(Th_plate);
    diag.stage_min_Th_minus_Tc(k) = min(margin);
    diag.has_temperature_reversal = diag.has_temperature_reversal || any(margin <= 0);
    pairs = make_particle_pair_indices_ms(n_vec(k));
    Tc_pair = pair_average_values_ms(Tc_plate, pairs);
    Th_pair = pair_average_values_ms(Th_plate, pairs);
    Qc_pair = NaN(npair(k), 1);
    Qh_pair = NaN(npair(k), 1);
    for j = 1:npair(k)
        Qc_pair(j) = te_Qc_onecouple_ms(Tc_pair(j), Th_pair(j), G.I, G, k);
        Qh_pair(j) = te_Qh_onecouple_ms(Th_pair(j), Tc_pair(j), G.I, G, k);
        if ~(isfinite(Qc_pair(j)) && isfinite(Qh_pair(j)))
            diag.reason = sprintf('stage%d_nonfinite_te_heat', k);
            return;
        end
    end
    diag.sumQc_stage(k) = sum(Qc_pair);
    diag.sumQh_stage(k) = sum(Qh_pair);
    diag.Qc_pair_min(k) = min(Qc_pair);
    diag.Qh_pair_min(k) = min(Qh_pair);
    if any(Qc_pair < -tol_neg_qc)
        diag.reason = sprintf('stage%d_negative_Qc_pair', k);
        return;
    end
end
diag.reason = 'ok';
ok = true;
end

function tol = get_qc_pair_negative_tol_ms(spec)
tol = 1e-3;
if nargin < 1 || ~isstruct(spec) || ~isfield(spec, 'k25_debug') || ...
        ~isstruct(spec.k25_debug) || ~isfield(spec.k25_debug, 'physics') || ...
        ~isstruct(spec.k25_debug.physics) || ...
        ~isfield(spec.k25_debug.physics, 'Qc_pair_negative_tol')
    return;
end
tol_in = spec.k25_debug.physics.Qc_pair_negative_tol;
if ~isempty(tol_in) && isnumeric(tol_in) && isfinite(tol_in(1))
    tol = max(0, tol_in(1));
end
end

function txt = format_physics_diag_ms(diag)
if nargin < 1 || ~isstruct(diag)
    txt = 'physics_diag_unavailable';
    return;
end
txt = sprintf(['physics_reason=%s; stage_Tc_mean=%s; stage_Th_mean=%s; ' ...
    'stage_min_Th_minus_Tc=%s; has_temperature_reversal=%d; sumQc_stage=%s; sumQh_stage=%s; ' ...
    'Qc_pair_min=%s; Qh_pair_min=%s'], ...
    char(string(diag.reason)), vec_to_inline_str_ms(diag.stage_Tc_mean), ...
    vec_to_inline_str_ms(diag.stage_Th_mean), vec_to_inline_str_ms(diag.stage_min_Th_minus_Tc), ...
    logical(diag.has_temperature_reversal), ...
    vec_to_inline_str_ms(diag.sumQc_stage), vec_to_inline_str_ms(diag.sumQh_stage), ...
    vec_to_inline_str_ms(diag.Qc_pair_min), vec_to_inline_str_ms(diag.Qh_pair_min));
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

function [best_out, status_msg, attempts] = newton_solve_with_continuation_ms(G, plates, n_vec, C_init, N, spec)
status_msg = 'schedule_not_started';
attempts = repmat(empty_newton_attempt_ms(), 0, 1);
best_out = struct('C', C_init(:).', 'theta', {cell(1, N)}, 'g', NaN(N,1), ...
    'cache', struct('abort_reason', status_msg), 'success', false, 'iters', 0, ...
    'relaxed', false, 'convergence_mode', 'continuation_not_started', 'rel_max', inf);
I_target = G.I;
if ~(isfinite(I_target) && I_target > 0)
    status_msg = 'invalid_target_current';
    best_out.cache.abort_reason = status_msg;
    best_out.convergence_mode = status_msg;
    return;
end
if nargin >= 6 && isstruct(spec)
    G = apply_low_k_solver_stabilization_ms(G, spec);
    frac_sched = get_current_continuation_fracs_ms(spec);
else
    frac_sched = [0.75, 0.875, 0.95, 1.0];
end
I_sched = I_target * frac_sched;
I_sched = I_sched(isfinite(I_sched) & I_sched > 0);
if isempty(I_sched)
    status_msg = 'empty_schedule';
    best_out.cache.abort_reason = status_msg;
    best_out.convergence_mode = status_msg;
    return;
end

C_seed = C_init(:).';
for i = 1:numel(I_sched)
    Gs = G;
    Gs.I = I_sched(i);
    if nargin >= 6 && isstruct(spec)
        [ok_seed, C_0d] = solve_0d_at_current_ms(n_vec, Gs.I, Gs, spec, C_seed(:));
        if ok_seed
            C_seed = C_0d(:).';
        end
    end
    cand_i = newton_solve_Ns_ms(Gs, plates, n_vec, C_seed, N);
    attempts(end+1,1) = pack_newton_attempt_ms(sprintf('continuation_%d', i), I_sched(i), cand_i); %#ok<AGROW>
    if ~cand_i.success
        status_msg = sprintf('schedule_fail_at_I=%.3fA', I_sched(i));
        best_out = cand_i;
        return;
    end
    C_seed = cand_i.C(:).';
    best_out = cand_i;
end
status_msg = 'schedule_success';
end

function rec = copy_newton_diag_to_eval_ms(rec, best, mode)
if nargin < 3
    mode = '';
end
if ~isstruct(best)
    return;
end
if isfield(best, 'iters'), rec.newton_iters = best.iters; end
if isfield(best, 'relaxed'), rec.newton_relaxed = logical(best.relaxed); end
if isfield(best, 'rel_max'), rec.newton_rel_max = best.rel_max; end
if isfield(best, 'convergence_mode')
    rec.newton_convergence_mode = char(string(best.convergence_mode));
end
if ~isempty(mode)
    rec.newton_convergence_mode = char(string(mode));
end
end

function attempt = empty_newton_attempt_ms()
attempt = struct('label', '', 'I_A', NaN, 'success', false, 'iters', NaN, ...
    'rel_max', NaN, 'g_norm', NaN, 'mode', '', 'abort_reason', '');
end

function attempt = pack_newton_attempt_ms(label, I_A, best)
attempt = empty_newton_attempt_ms();
attempt.label = char(string(label));
attempt.I_A = I_A;
if ~isstruct(best)
    return;
end
if isfield(best, 'success'), attempt.success = logical(best.success); end
if isfield(best, 'iters'), attempt.iters = best.iters; end
if isfield(best, 'rel_max'), attempt.rel_max = best.rel_max; end
if isfield(best, 'convergence_mode'), attempt.mode = char(string(best.convergence_mode)); end
if isfield(best, 'g') && ~isempty(best.g)
    gv = best.g(:);
    gv = gv(isfinite(gv));
    if ~isempty(gv)
        attempt.g_norm = norm(gv);
    end
end
if isfield(best, 'cache') && isstruct(best.cache) && isfield(best.cache, 'abort_reason')
    attempt.abort_reason = char(string(best.cache.abort_reason));
end
end

function txt = format_newton_attempt_detail_ms(label, I_A, best)
attempt = pack_newton_attempt_ms(label, I_A, best);
txt = sprintf('%s(I=%.6g,success=%d,iters=%.6g,rel_max=%.6g,g_norm=%.6g,mode=%s,abort=%s)', ...
    attempt.label, attempt.I_A, attempt.success, attempt.iters, attempt.rel_max, ...
    attempt.g_norm, attempt.mode, attempt.abort_reason);
end

function txt = format_newton_attempts_ms(attempts)
if isempty(attempts)
    txt = 'continuation_attempts=[]';
    return;
end
parts = cell(numel(attempts), 1);
for i = 1:numel(attempts)
    a = attempts(i);
    parts{i} = sprintf('%s(I=%.6g,success=%d,iters=%.6g,rel_max=%.6g,g_norm=%.6g,mode=%s,abort=%s)', ...
        a.label, a.I_A, a.success, a.iters, a.rel_max, a.g_norm, a.mode, a.abort_reason);
end
txt = ['continuation_attempts=[', strjoin(parts, '; '), ']'];
end

function best = newton_solve_Ns_ms(G, plates, n_vec, C_init_vec, N)
seeds = build_newton_initial_seeds_ms(C_init_vec, N, G);
best = struct('C', seeds(1,:), 'theta', {cell(1, N)}, 'g', NaN(N,1), ...
    'cache', struct([]), 'success', false, 'iters', 0);
best_gnorm = inf;
has_failed_candidate = false;
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
    elseif ~best.success && (~has_failed_candidate || gnorm < best_gnorm)
        best = cand;
        best_gnorm = gnorm;
        has_failed_candidate = true;
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
        best.theta = theta;
        best.g = g;
        best.cache = cache;
        best.iters = it;
        best.rel_max = calc_rel_max_newton_ms(g, cache, N);
        best.convergence_mode = 'nonfinite_g';
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
for k = 1:N-1
    refs(k) = max(abs(cache.sumQh(k+1)) + abs(cache.sumQc(k)), 1e-12);
end
refs(N) = max(abs(cache.sumQc(N)) + 1e-12, 1e-12);
rel = abs(g) ./ refs;
if all(isfinite(rel))
    rel_max = max(rel);
end
end

function [g, theta, cache] = eval_gvec_Ns_ms(C_vec, theta_cells_0, G, plates, n_vec, N)
[theta, cache] = solve_theta_given_C_ms(C_vec, theta_cells_0, G, plates, n_vec, N);
g = NaN(N,1);
if isfield(cache, 'inner_converged') && ~logical(cache.inner_converged)
    return;
end
if ~isfield(cache, 'sumQc') || ~isfield(cache, 'sumQh') || ...
        numel(cache.sumQc) ~= N || numel(cache.sumQh) ~= N
    return;
end
for k = 1:N-1
    g(k) = cache.sumQh(k+1) - cache.sumQc(k);
end
g(N) = cache.sumQc(N) - G.Qc_target_last;
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
zcfg = get_vertical_runtime_from_G_ms(G, N);
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
    cache.Qc_pair_min = NaN(N,1);
    cache.Qh_pair_min = NaN(N,1);
    cache.step_q_prev_weight = w_prev;
    cache.Rz_interfaces = zcfg.Rz_interfaces(:).';
    cache.Rz_sink = zcfg.Rz_sink;
    cache.inner_converged = false;
    cache.inner_iters = it;
    cache.inner_dmax = dmax;
    for k = 1:N
        cache.sumQc(k) = sum(Qc_eff{k});
        cache.sumQh(k) = sum(Qh_eff{k});
        cache.Qc_pair_min(k) = min(Qc_eff{k});
        cache.Qh_pair_min(k) = min(Qh_eff{k});
    end

    if dmax < G.tol_theta
        cache.inner_converged = true;
        return;
    end
end
if ~exist('cache', 'var') || ~isstruct(cache)
    cache = make_nonfinite_cache_ms(zcfg, w_prev, 'inner_theta_not_converged', N);
else
    cache.abort_reason = 'inner_theta_not_converged';
    cache.inner_converged = false;
    cache.inner_iters = G.max_inner;
    if ~isfield(cache, 'inner_dmax')
        cache.inner_dmax = NaN;
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
cache.Qc_pair_min = NaN(N,1);
cache.Qh_pair_min = NaN(N,1);
cache.step_q_prev_weight = w_prev;
cache.Rz_interfaces = zcfg.Rz_interfaces(:).';
cache.Rz_sink = zcfg.Rz_sink;
cache.inner_converged = false;
cache.inner_iters = NaN;
cache.inner_dmax = NaN;
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
    'spacing_ratio', NaN, 'contrast_score', NaN, 'mode_prior_score', NaN, ...
    'lmax_relaxed', false, 'lmax_relax_ratio', 1.0, ...
    'n', NaN(1,N), 'npair', NaN(1,N), 'ratios', NaN(1,max(0,N-1)), 'I_opt', NaN, ...
    'Lx', NaN(1,N), 'Ly', NaN(1,N), 'cov', NaN(1,N), ...
    'fp_fill_count', 0, 'fp_fill_by_stage', zeros(1,N), ...
    'newton_relaxed', false, 'newton_rel_max', NaN, 'newton_iters', NaN, 'newton_convergence_mode', 'none', ...
    'newton_failure_detail', '', ...
    'physics_valid', false, 'physics_failure_reason', '', ...
    'stage_Tc_mean', NaN(1,N), 'stage_Th_mean', NaN(1,N), 'stage_min_Th_minus_Tc', NaN(1,N), ...
    'has_temperature_reversal', false, ...
    'sumQc_stage', NaN(1,N), 'sumQh_stage', NaN(1,N), ...
    'Qc_pair_min', NaN(1,N), 'Qh_pair_min', NaN(1,N), ...
    'C', NaN(1,N), ...
    'DeltaTN_actual', NaN, 'target_met', false, 'Qc_last_total', NaN, 'Qc_target_last', NaN, 'Qc_error', NaN, 'Qc_tol', NaN, ...
    'TN_min', NaN, 'TN_mean', NaN, 'DeltaTN_mean', NaN, 'TN_maxmin', NaN, ...
    'stage_spread', NaN(1,N), 'rank_score', NaN, 'kept_from_fast', false, ...
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

function n_bal = rebalance_counts_to_total_ms(n_vec, n_raw, target_total, force_even)
n_bal = n_vec(:).';
target_total = round(target_total);
if ~(isfinite(target_total) && target_total >= numel(n_bal))
    return;
end
if nargin < 4
    force_even = false;
end
if force_even
    n_bal = max(2, n_bal);
    n_bal = 2 * round(n_bal / 2);
else
    n_bal = max(1, n_bal);
end
if ~isfinite(sum(n_bal))
    return;
end
if nargin < 2 || isempty(n_raw) || numel(n_raw) ~= numel(n_bal)
    n_raw = n_bal;
end
n_raw = n_raw(:).';

if force_even && mod(target_total, 2) ~= 0
    target_total = target_total - 1;
end
step_unit = 1;
if force_even
    step_unit = 2;
end
delta = target_total - sum(n_bal);
if delta == 0
    return;
end
if delta > 0
    [~, ord] = sort(n_raw - n_bal, 'descend');
else
    [~, ord] = sort(n_bal - n_raw, 'descend');
end
ord = ord(:).';
if isempty(ord)
    ord = 1:numel(n_bal);
end
ptr = 1;
while delta ~= 0
    idx = ord(ptr);
    if delta > 0
        n_bal(idx) = n_bal(idx) + step_unit;
        delta = delta - step_unit;
    else
        floor_n = 1;
        if force_even
            floor_n = 2;
        end
        if n_bal(idx) - step_unit >= floor_n
            n_bal(idx) = n_bal(idx) - step_unit;
            delta = delta + step_unit;
        end
    end
    ptr = ptr + 1;
    if ptr > numel(ord)
        ptr = 1;
    end
    if numel(ord) == 1 && delta ~= 0
        break;
    end
end
end

function v = safe_stage_value_ms(s, fname, idx, default_v)
v = default_v;
if ~isstruct(s) || ~isfield(s, fname)
    return;
end
arr = s.(fname);
if numel(arr) < idx
    return;
end
vv = arr(idx);
if isnumeric(vv) && isfinite(vv)
    v = vv;
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
write_struct_rows_csv_ms(rows, out_csv);
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

function zcfg = get_vertical_runtime_from_G_ms(G, N)
if nargin < 2 || ~isfinite(N) || N < 2
    N = 2;
end
nif = max(1, N-1);
if isstruct(G) && isfield(G, 'vertical_runtime') && isstruct(G.vertical_runtime)
    zcfg = G.vertical_runtime;
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

function infos = save_layout_stage_plots_ms(eval_rec, out_dir, spec)
N = spec.stage_count;
plot_cfg = resolve_plot_cfg_ms(spec);
plot_cfg.use_global_axis = false;
global_axis_mm = make_global_axis_limits_mm_ms(spec);
fig_size = normalize_plot_size_px_ms(spec.output.plot.separate_fig_size_px);
infos = repmat(empty_plot_axis_info_ms(1), N, 1);
for k = 1:N
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 fig_size(1) fig_size(2)]);
    ax = axes(fig);
    [axis_target_k, axis_actual_k, substrate_k] = render_layout_stage_on_axis_ms(ax, eval_rec, k, plot_cfg, global_axis_mm);
    infos(k).target_axes_mm = axis_target_k;
    infos(k).actual_axes_mm = axis_actual_k;
    infos(k).substrate_mm = substrate_k;
    exportgraphics(fig, fullfile(out_dir, sprintf('layout_stage%d.png', k)), 'Resolution', 180);
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

function infos = save_temperature_stage_plots_ms(eval_rec, out_dir, caxis_by_stage, spec)
N = spec.stage_count;
plot_cfg = resolve_plot_cfg_ms(spec);
plot_cfg.use_global_axis = false;
global_axis_mm = make_global_axis_limits_mm_ms(spec);
fig_size = normalize_plot_size_px_ms(spec.output.plot.separate_fig_size_px);
infos = repmat(empty_plot_axis_info_ms(1), N, 1);
for k = 1:N
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 fig_size(1) fig_size(2)]);
    colormap(fig, 'parula');
    ax = axes(fig);
    [axis_target_k, axis_actual_k, substrate_k] = render_temperature_stage_on_axis_ms( ...
        ax, eval_rec, k, caxis_by_stage, plot_cfg, global_axis_mm);
    infos(k).target_axes_mm = axis_target_k;
    infos(k).actual_axes_mm = axis_actual_k;
    infos(k).substrate_mm = substrate_k;
    exportgraphics(fig, fullfile(out_dir, sprintf('temperature_stage%d.png', k)), 'Resolution', 220);
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
