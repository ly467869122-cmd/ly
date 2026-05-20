% Modified 2026-05-12 20:51:42 +08:00:
%   Add spec.metrics defaults for DeltaT_at_Qc, DeltaTmax, and Qmax post-processing.
function [spec, G] = optimize_layout_multistage0411_shared_params(spec_in, profile)
% Shared parameter entry for the staged layout/count/current scripts.
% All structural and material defaults live here. Caller overrides are merged
% through spec_in, then normalized into a single spec/G pair.

if nargin < 1 || isempty(spec_in)
    spec_in = struct();
end
if nargin < 2 || isempty(profile)
    profile = 'layout';
end
if ~isstruct(spec_in)
    error('spec_in must be a struct.');
end
profile = lower(strtrim(char(string(profile))));
if ~any(strcmp(profile, {'n', 'layout', 'current'}))
    error('profile must be one of: n, layout, current.');
end

N0 = 5;
if isfield(spec_in, 'stage_count') && ~isempty(spec_in.stage_count) && isfinite(spec_in.stage_count)
    N0 = max(2, round(double(spec_in.stage_count)));
end

defaults = build_shared_defaults_ms(N0, profile);
spec = merge_struct_recursive_ms(defaults, spec_in);
spec.stage_count = max(2, round(double(spec.stage_count)));
N = spec.stage_count;

spec = normalize_shared_spec_ms(spec, defaults, N, profile);
G = build_shared_G_ms(spec);
end

function defaults = build_shared_defaults_ms(N, profile)
base_fixed_n = [480, 246, 114, 64, 32];

defaults = struct();
defaults.seed = 20260407;
defaults.stage_count = N;
defaults.fp_w = 1.5e-3;
defaults.fp_h = 1.5e-3;
defaults.fp_t_stage = fit_len_vec_ms([2.1, 1.3, 1.1, 1.1, 1.1] * 1e-3, N, 1.1e-3);
defaults.fp_t = defaults.fp_t_stage;
defaults.edge_margin_fp = 0.5;
defaults.min_spacing_manufacture = 2.20e-3;
defaults.mesh_nx = 61;
defaults.mesh_ny = 61;
defaults.mesh_nx_stage_full = fit_len_vec_ms([81, 81, 61, 61, 61], N, 61);
defaults.mesh_ny_stage_full = fit_len_vec_ms([81, 81, 61, 61, 61], N, 61);
defaults.plate_t = 1e-3;
defaults.plate_k_inplane = 170;
defaults.fixed_n = fit_len_vec_ms(base_fixed_n, N, base_fixed_n);
defaults.omega = 0.7;
defaults.max_inner = 60;
defaults.max_outer = 30;
defaults.tol_theta = 1e-3;
defaults.tol_g_rel = 1e-3;
defaults.use_parallel = true;
defaults.eval_context = 'full';

defaults.targets = struct('DeltaT_target', 107, 'Qc_target_last', 0.8, 'ThRes', 300);
defaults.geometry = struct( ...
    'L_max_mm', fit_len_vec_ms([80, 80, 55, 55, 55], N, 55), ...
    'coverage_min', fit_len_vec_ms([0.20, 0.20, 0.20, 0.20, 0.20], N, 0.20), ...
    'pyramid_gap_min_mm', fit_len_vec_ms([0, 0, 0, 0], max(1, N-1), 0), ...
    'min_edge_gap_mm', 0.7);

defaults.material = struct( ...
    'alpha_coeffs', [-8.5940e-14, 8.6090e-11, -3.3765e-08, 7.2220e-06, -3.2951e-04], ...
    'k_coeffs', [-2.22485e-09, 2.13815e-06, -0.0007046819, 0.088625, 0.314], ...
    'rho_coeffs', [1.5671e-15, -1.3921e-12, 6.1640e-10, -4.4950e-08, 4.7288e-06], ...
    'Rc', 0.0);

defaults.scan = struct('I_coarse', 1.0:0.5:8.0, 'I_fine_step', 0.1, 'I_fine_half_window', 0.5);
defaults.current = struct('I_init', 3.8);
defaults.current_opt = struct('I_list', 1:0.2:6, 'qc_abs_tol', 1e-3, 'qc_rel_tol', 0);
defaults.n_current_scan = struct('enable', true, 'I_list', 2.0:0.2:5.4, ...
    'topK', 3, 'qc_abs_tol', 1e-3, 'qc_rel_tol', 0, ...
    'select_rule', 'min_total');
defaults.metrics = struct('mode', 'deltaT_at_Qc', 'enable_qmax_post', true, ...
    'qmax', struct('enable', true, 'fem_topI', 3, 'qc_min', 0, ...
    'qc_max_auto_factor', 2.0, 'tol_W', 1e-3, 'max_iter', 20));
defaults.current_calib = struct('enable', true, 'topN', 4, 'left_span', 1.0, ...
    'right_span', 0.3, 'step', 0.1, 'fallback_seed_count', 24, ...
    'fallback_keep_success_min', 1);
defaults.warmup = struct('enable', true, 'sample_ratio', 0.25, 'min_per_stratum', 1, ...
    'topN', 16, 'min_sample_jobs', 1200);
defaults.soft_prior = struct('fallback_stage_trend', {{'center_heavy','center_heavy','neutral','edge_heavy','edge_heavy'}}, ...
    'vote_margin_min', 1, 'weight_margin_min', 0.15, 'min_keep_per_method', 6);
defaults.parallel = struct('eval_min_tasks', 2, 'block_min', 16, 'block_max', 128, ...
    'target_blocks_per_worker', 2, 'pool_workers', 64, ...
    'step1_enable', true, 'step1_min_tasks', 24, ...
    'pool_used', false, 'pool_workers_actual', 0);
defaults.output = struct('topK', 10, 'top_postK', 5, 'top_plotK', 5, ...
    'output_dir', fullfile(pwd, 'optimize_layout_multistage_output'), ...
    'plot', struct('use_global_axis', true, 'show_mesh_edges', false, ...
    'view_mode', '2d', 'axis_margin_ratio', 0.05, ...
    'save_overview', true, 'save_stage_separate', true, ...
    'annotate_substrate_dims', true, 'separate_fig_size_px', [900, 900]));

defaults.mode_list = {'center_dense', 'edge_dense'};
defaults.layout_methods = {'subset_symmetric', 'fixed_c4_grid', 'hex_c6', 'ring_stage3', 'interstage_aligned', 'poisson_disk'};
defaults.symmetry_mode_list = {'c4'};
defaults.edge_pattern_mode_list = {'free', 'edge_spaced'};
defaults.coarse_s_dense_list = [2.20, 2.40, 2.60, 2.80] * 1e-3;
defaults.coarse_s_sparse_list = [2.20, 2.70, 3.00, 3.40, 3.90, 4.40, 4.90] * 1e-3;
defaults.coarse_expo_list = [2.0, 3.2, 4.5, 6.0, 8.0];
defaults.coarse_anis_ratio_list = [0.65, 0.75, 1.00, 1.35, 1.55];
defaults.gamma_list = [-1.0, -0.5, 0.0, 0.5, 1.0];
defaults.hex = struct('enable', true, 's_list', [2.20, 2.6] * 1e-3);
defaults.ring = struct('enable', true);
defaults.candidate_dedup = struct('enable', true, 'tol_m', 1e-9);
defaults.rank = struct('w_actual', 0.6, 'w_mean', 0.3, 'w_spread', 0.1);
defaults.poisson = struct('seed_list', [42, 137, 271, 618, 1001]);
defaults.interstage = struct('align_weights', [0.15, 0.35, 0.55, 0.75, 0.90]);
defaults.hard_constraints = struct('force_c4_only', true, 'unified_lmax_across_stages', false, ...
    'allow_lmax_relax', false, 'enforce_monotonic_geometry', true, ...
    'enforce_coverage_min', true, 'enforce_spacing_min_edge_gap', true);
defaults.uniform_baseline = struct('gap_factor_x', 0.5, 'gap_factor_y', 0.5, ...
    'edge_margin_factor_x', 0.5, 'edge_margin_factor_y', 0.5);
defaults.control_baseline = struct('enable', true, 'ratio_cold_to_hot', [1, 3, 7, 21, 55], ...
    'anchor_last_stage_count', false, 'force_even', true, 'scale_list', 1.00, ...
    'gap_factor_x', 0.5, 'gap_factor_y', 0.5, ...
    'edge_margin_factor_x', 0.5, 'edge_margin_factor_y', 0.5, ...
    'fallback_to_uniform', false);
defaults.method_mix = struct('enable', false);
defaults.candidate_batch = struct('size', 64, 'top_bottom_n', 10, 'convergence_gap_K', 0.5);
defaults.shape_explore = struct('enable', true, ...
    'extra_stage_modes', {{'ring_dense','near_center_dense','corner_dense','band_dense_x','band_dense_y','multi_center'}}, ...
    'stage_anis_enable', true, 'anis_ratio_stage_list', [0.65, 0.85, 1.00, 1.25, 1.55], ...
    'ring_radius_ratio_list', [0.25, 0.35, 0.45, 0.60], 'ring_width_ratio_list', [0.10, 0.16, 0.22], ...
    'band_width_ratio_list', [0.18, 0.28, 0.38], 'corner_bias_list', [0.30, 0.45, 0.60], ...
    'jitter_ratio_list', [0, 0.04, 0.08, 0.12], 'jitter_seed_list', [101, 203, 307, 409, 503, 607], ...
    'stage_mode_templates', {{ ...
        {'center_dense','center_dense','edge_dense','edge_dense','edge_dense'}, ...
        {'center_dense','center_dense','center_dense','edge_dense','edge_dense'}, ...
        {'center_dense','center_dense','ring_dense','ring_dense','edge_dense'}, ...
        {'center_dense','near_center_dense','ring_dense','edge_dense','edge_dense'}, ...
        {'near_center_dense','center_dense','neutral','ring_dense','edge_dense'} ...
    }});
defaults.job_budget = struct('enable', true, 'method_caps', struct( ...
    'subset_symmetric', 360, 'fixed_c4_grid', 120, 'hex_c6', 120, ...
    'gamma_stage2', 160, 'ring_stage3', 120, 'poisson_disk', 160, ...
    'interstage_aligned', 96, 'mixed_stage_methods', 320, 'shape_explore', 720));
defaults.z_path = struct('enable', false);

defaults.k25_debug = struct('enable', true, 'total_n', 608, 'plate_k_inplane', defaults.plate_k_inplane, ...
    'search_only', false, 'n', [], ...
    'particle_search', struct('n5_max', [], 'stage_max_count', 900, ...
    'mode', 'template_0d_pruned', 'strict_min_total', false, ...
    'total_n_min', 700, 'total_n_max', 2000, 'candidate_budget', 1000, ...
    'total_block_size', 4, 'par_chunk_size', 2000, ...
    'total_step', 4, 'n5_min_factor', 0.7, ...
    'n5_max_auto_factor', 2.5, 'n5_max_auto_margin', 20, ...
    'use_ratio_prune', true, 'adjacent_min_fraction', 0.1, ...
    'boundary_check_count', 2, 'stop_after_best_margin', 120, ...
    'seed_ratios', [55, 21, 7, 3, 1]), ...
    'physics', struct('Qc_pair_negative_tol', 1e-3), ...
    'low_k_solver', struct('enable', true, 'k_threshold', 80, ...
    'omega', 0.35, 'max_inner', 160, 'max_outer', 60, ...
    'continuation_fracs', [0.25, 0.50, 0.75, 0.875, 0.95, 1.0]));

if strcmp(profile, 'n')
    defaults.output.top_postK = defaults.output.topK;
    defaults.output.plot.use_global_axis = false;
elseif strcmp(profile, 'layout')
    defaults.current_calib.enable = false;
    defaults.output.top_postK = 5;
end
end

function spec = normalize_shared_spec_ms(spec, defaults, N, profile)
spec.fp_w = max(1e-12, scalar_with_default_ms(spec.fp_w, defaults.fp_w));
spec.fp_h = max(1e-12, scalar_with_default_ms(spec.fp_h, defaults.fp_h));
spec.fp_t_stage = resolve_fp_t_stage_ms(spec, defaults.fp_t_stage, N);
spec.fp_t = spec.fp_t_stage;
spec.edge_margin_fp = max(0, scalar_with_default_ms(spec.edge_margin_fp, defaults.edge_margin_fp));
spec.min_spacing_manufacture = max(0, scalar_with_default_ms(spec.min_spacing_manufacture, defaults.min_spacing_manufacture));
spec.mesh_nx = max(5, round(scalar_with_default_ms(spec.mesh_nx, defaults.mesh_nx)));
spec.mesh_ny = max(5, round(scalar_with_default_ms(spec.mesh_ny, defaults.mesh_ny)));
spec.mesh_nx_stage_full = max(5, round(fit_len_vec_ms(spec.mesh_nx_stage_full, N, defaults.mesh_nx_stage_full)));
spec.mesh_ny_stage_full = max(5, round(fit_len_vec_ms(spec.mesh_ny_stage_full, N, defaults.mesh_ny_stage_full)));
spec.plate_t = max(0, scalar_with_default_ms(spec.plate_t, defaults.plate_t));
spec.plate_k_inplane = max(1e-12, scalar_with_default_ms(spec.plate_k_inplane, defaults.plate_k_inplane));
if isfield(spec, 'k25_debug') && isstruct(spec.k25_debug)
    spec.k25_debug.plate_k_inplane = spec.plate_k_inplane;
end

spec.targets.DeltaT_target = scalar_with_default_ms(spec.targets.DeltaT_target, defaults.targets.DeltaT_target);
spec.targets.Qc_target_last = scalar_with_default_ms(spec.targets.Qc_target_last, defaults.targets.Qc_target_last);
spec.targets.ThRes = scalar_with_default_ms(spec.targets.ThRes, defaults.targets.ThRes);
spec.geometry.L_max_mm = fit_len_vec_ms(spec.geometry.L_max_mm, N, defaults.geometry.L_max_mm);
spec.geometry.coverage_min = fit_len_vec_ms(spec.geometry.coverage_min, N, defaults.geometry.coverage_min);
spec.geometry.pyramid_gap_min_mm = fit_len_vec_ms(spec.geometry.pyramid_gap_min_mm, max(1, N-1), defaults.geometry.pyramid_gap_min_mm);
spec.geometry.min_edge_gap_mm = max(0, scalar_with_default_ms(spec.geometry.min_edge_gap_mm, defaults.geometry.min_edge_gap_mm));
spec.geometry.L_max = spec.geometry.L_max_mm * 1e-3;
spec.geometry.pyramid_gap_min = spec.geometry.pyramid_gap_min_mm * 1e-3;
spec.geometry.min_edge_gap = spec.geometry.min_edge_gap_mm * 1e-3;

spec.fixed_n = normalize_even_vector_ms(fit_len_vec_ms(spec.fixed_n, N, defaults.fixed_n), N);
spec.omega = min(max(scalar_with_default_ms(spec.omega, defaults.omega), 0.01), 1.0);
spec.max_inner = max(1, round(scalar_with_default_ms(spec.max_inner, defaults.max_inner)));
spec.max_outer = max(1, round(scalar_with_default_ms(spec.max_outer, defaults.max_outer)));
spec.tol_theta = max(1e-12, scalar_with_default_ms(spec.tol_theta, defaults.tol_theta));
spec.tol_g_rel = max(1e-12, scalar_with_default_ms(spec.tol_g_rel, defaults.tol_g_rel));
spec.use_parallel = any(logical(spec.use_parallel(:)));

spec.current.I_init = max(1e-9, scalar_with_default_ms(spec.current.I_init, defaults.current.I_init));
if isfield(spec, 'scan') && isstruct(spec.scan) && isfield(spec.scan, 'I_coarse') && ~isempty(spec.scan.I_coarse)
    spec.scan.I_coarse = spec.scan.I_coarse(isfinite(spec.scan.I_coarse) & spec.scan.I_coarse > 0);
    if isempty(spec.scan.I_coarse), spec.scan.I_coarse = defaults.scan.I_coarse; end
    spec.scan.I_fine_step = max(1e-6, scalar_with_default_ms(spec.scan.I_fine_step, defaults.scan.I_fine_step));
    spec.scan.I_fine_half_window = max(0, scalar_with_default_ms(spec.scan.I_fine_half_window, defaults.scan.I_fine_half_window));
    spec.current.I_init = max(min(spec.current.I_init, max(spec.scan.I_coarse)), min(spec.scan.I_coarse));
end

spec.parallel.eval_min_tasks = max(2, round(scalar_with_default_ms(spec.parallel.eval_min_tasks, defaults.parallel.eval_min_tasks)));
spec.parallel.block_min = max(1, round(scalar_with_default_ms(spec.parallel.block_min, defaults.parallel.block_min)));
spec.parallel.block_max = max(spec.parallel.block_min, round(scalar_with_default_ms(spec.parallel.block_max, defaults.parallel.block_max)));
spec.parallel.target_blocks_per_worker = max(1, round(scalar_with_default_ms(spec.parallel.target_blocks_per_worker, defaults.parallel.target_blocks_per_worker)));
spec.parallel.pool_workers = max(1, round(scalar_with_default_ms(spec.parallel.pool_workers, defaults.parallel.pool_workers)));
spec.parallel.step1_enable = logical_field_ms(spec.parallel, 'step1_enable', defaults.parallel.step1_enable);
spec.parallel.step1_min_tasks = max(1, round(scalar_with_default_ms(spec.parallel.step1_min_tasks, defaults.parallel.step1_min_tasks)));

spec.output.topK = max(1, round(scalar_with_default_ms(spec.output.topK, defaults.output.topK)));
spec.output.top_postK = max(1, min(spec.output.topK, round(scalar_with_default_ms(spec.output.top_postK, defaults.output.top_postK))));
spec.output.top_plotK = max(0, min(spec.output.topK, round(scalar_with_default_ms(spec.output.top_plotK, defaults.output.top_plotK))));
spec.output.plot.use_global_axis = logical_field_ms(spec.output.plot, 'use_global_axis', defaults.output.plot.use_global_axis);
spec.output.plot.show_mesh_edges = logical_field_ms(spec.output.plot, 'show_mesh_edges', defaults.output.plot.show_mesh_edges);
spec.output.plot.view_mode = normalize_plot_view_mode_ms(spec.output.plot.view_mode);
spec.output.plot.axis_margin_ratio = min(max(scalar_with_default_ms(spec.output.plot.axis_margin_ratio, defaults.output.plot.axis_margin_ratio), 0), 0.5);
spec.output.plot.save_overview = logical_field_ms(spec.output.plot, 'save_overview', defaults.output.plot.save_overview);
spec.output.plot.save_stage_separate = logical_field_ms(spec.output.plot, 'save_stage_separate', defaults.output.plot.save_stage_separate);
spec.output.plot.annotate_substrate_dims = logical_field_ms(spec.output.plot, 'annotate_substrate_dims', defaults.output.plot.annotate_substrate_dims);
spec.output.plot.separate_fig_size_px = normalize_plot_size_px_ms(spec.output.plot.separate_fig_size_px);

spec.z_path = normalize_z_path_ms(spec.z_path, defaults.z_path, N, spec.plate_t);

if isfield(spec, 'current_opt') && isstruct(spec.current_opt)
    spec.current_opt.I_list = normalize_positive_values_ms(spec.current_opt.I_list, defaults.current_opt.I_list);
    spec.current_opt.qc_abs_tol = max(0, scalar_with_default_ms(spec.current_opt.qc_abs_tol, defaults.current_opt.qc_abs_tol));
    spec.current_opt.qc_rel_tol = max(0, scalar_with_default_ms(spec.current_opt.qc_rel_tol, defaults.current_opt.qc_rel_tol));
end
if isfield(spec, 'n_current_scan') && isstruct(spec.n_current_scan)
    spec.n_current_scan.enable = logical_field_ms(spec.n_current_scan, 'enable', defaults.n_current_scan.enable);
    spec.n_current_scan.I_list = normalize_positive_values_ms(spec.n_current_scan.I_list, defaults.n_current_scan.I_list);
    spec.n_current_scan.topK = max(1, round(scalar_with_default_ms(spec.n_current_scan.topK, defaults.n_current_scan.topK)));
    spec.n_current_scan.qc_abs_tol = max(0, scalar_with_default_ms(spec.n_current_scan.qc_abs_tol, defaults.n_current_scan.qc_abs_tol));
    spec.n_current_scan.qc_rel_tol = max(0, scalar_with_default_ms(spec.n_current_scan.qc_rel_tol, defaults.n_current_scan.qc_rel_tol));
    spec.n_current_scan.select_rule = char(string(spec.n_current_scan.select_rule));
else
    spec.n_current_scan = defaults.n_current_scan;
end

if isfield(spec, 'metrics') && isstruct(spec.metrics)
    spec.metrics = normalize_metrics_ms(spec.metrics, defaults.metrics);
else
    spec.metrics = defaults.metrics;
end

if isfield(spec, 'k25_debug') && isstruct(spec.k25_debug)
    spec.k25_debug.enable = logical_field_ms(spec.k25_debug, 'enable', defaults.k25_debug.enable);
    spec.k25_debug.search_only = logical_field_ms(spec.k25_debug, 'search_only', defaults.k25_debug.search_only);
    spec.k25_debug.total_n = normalize_even_total_count_ms(scalar_with_default_ms(spec.k25_debug.total_n, defaults.k25_debug.total_n), N);
    if isfield(spec.k25_debug, 'n') && ~isempty(spec.k25_debug.n)
        spec.k25_debug.n = normalize_even_vector_ms(fit_len_vec_ms(spec.k25_debug.n, N, defaults.fixed_n), N);
    end
    spec.k25_debug.particle_search = normalize_particle_search_ms(spec.k25_debug.particle_search, defaults.k25_debug.particle_search, N, spec.k25_debug.total_n);
    spec.k25_debug.physics.Qc_pair_negative_tol = max(0, scalar_with_default_ms(spec.k25_debug.physics.Qc_pair_negative_tol, defaults.k25_debug.physics.Qc_pair_negative_tol));
    spec.k25_debug.low_k_solver = normalize_low_k_solver_ms(spec.k25_debug.low_k_solver, defaults.k25_debug.low_k_solver);
end

if ~strcmp(profile, 'n') && isfield(spec, 'current_calib') && isstruct(spec.current_calib)
    spec.current_calib.enable = logical_field_ms(spec.current_calib, 'enable', defaults.current_calib.enable);
    spec.current_calib.topN = max(1, round(scalar_with_default_ms(spec.current_calib.topN, defaults.current_calib.topN)));
    spec.current_calib.left_span = max(0, scalar_with_default_ms(spec.current_calib.left_span, defaults.current_calib.left_span));
    spec.current_calib.right_span = max(0, scalar_with_default_ms(spec.current_calib.right_span, defaults.current_calib.right_span));
    spec.current_calib.step = max(1e-6, scalar_with_default_ms(spec.current_calib.step, defaults.current_calib.step));
    spec.current_calib.fallback_seed_count = max(1, round(scalar_with_default_ms(spec.current_calib.fallback_seed_count, defaults.current_calib.fallback_seed_count)));
    spec.current_calib.fallback_keep_success_min = max(1, round(scalar_with_default_ms(spec.current_calib.fallback_keep_success_min, defaults.current_calib.fallback_keep_success_min)));
end

spec.mode_list = normalize_text_cell_ms(spec.mode_list, defaults.mode_list);
spec.layout_methods = normalize_text_cell_ms(spec.layout_methods, defaults.layout_methods);
spec.symmetry_mode_list = normalize_text_cell_ms(spec.symmetry_mode_list, defaults.symmetry_mode_list);
spec.edge_pattern_mode_list = normalize_text_cell_ms(spec.edge_pattern_mode_list, defaults.edge_pattern_mode_list);
spec.coarse_s_dense_list = normalize_positive_values_ms(spec.coarse_s_dense_list, defaults.coarse_s_dense_list);
spec.coarse_s_sparse_list = normalize_positive_values_ms(spec.coarse_s_sparse_list, defaults.coarse_s_sparse_list);
spec.coarse_expo_list = normalize_positive_values_ms(spec.coarse_expo_list, defaults.coarse_expo_list);
spec.coarse_anis_ratio_list = normalize_positive_values_ms(spec.coarse_anis_ratio_list, defaults.coarse_anis_ratio_list);
spec.gamma_list = normalize_numeric_values_ms(spec.gamma_list, defaults.gamma_list);

spec.material.alpha_coeffs = normalize_coeffs_ms(spec.material.alpha_coeffs, defaults.material.alpha_coeffs);
spec.material.k_coeffs = normalize_coeffs_ms(spec.material.k_coeffs, defaults.material.k_coeffs);
spec.material.rho_coeffs = normalize_coeffs_ms(spec.material.rho_coeffs, defaults.material.rho_coeffs);
spec.material.Rc = max(0, scalar_with_default_ms(spec.material.Rc, defaults.material.Rc));
end

function metrics = normalize_metrics_ms(metrics, defaults)
metrics = merge_struct_recursive_ms(defaults, metrics);
mode = lower(strtrim(char(string(metrics.mode))));
valid_modes = {'deltat_at_qc', 'deltatmax', 'qmax'};
if ~any(strcmp(mode, valid_modes))
    mode = defaults.mode;
end
if strcmp(mode, 'deltat_at_qc')
    mode = 'deltaT_at_Qc';
end
metrics.mode = mode;
metrics.enable_qmax_post = logical_field_ms(metrics, 'enable_qmax_post', defaults.enable_qmax_post);
if ~isfield(metrics, 'qmax') || ~isstruct(metrics.qmax)
    metrics.qmax = defaults.qmax;
else
    metrics.qmax = merge_struct_recursive_ms(defaults.qmax, metrics.qmax);
end
metrics.qmax.enable = logical_field_ms(metrics.qmax, 'enable', defaults.qmax.enable);
metrics.qmax.fem_topI = max(1, round(scalar_with_default_ms(metrics.qmax.fem_topI, defaults.qmax.fem_topI)));
metrics.qmax.qc_min = max(0, scalar_with_default_ms(metrics.qmax.qc_min, defaults.qmax.qc_min));
metrics.qmax.qc_max_auto_factor = max(1.01, scalar_with_default_ms(metrics.qmax.qc_max_auto_factor, defaults.qmax.qc_max_auto_factor));
metrics.qmax.tol_W = max(1e-12, scalar_with_default_ms(metrics.qmax.tol_W, defaults.qmax.tol_W));
metrics.qmax.max_iter = max(1, round(scalar_with_default_ms(metrics.qmax.max_iter, defaults.qmax.max_iter)));
end

function z = normalize_z_path_ms(raw, defaults, N, plate_t)
if nargin < 1 || ~isstruct(raw)
    raw = struct();
end
if nargin < 2 || ~isstruct(defaults)
    defaults = struct('enable', false);
end
if nargin < 3 || ~isfinite(N)
    N = 2;
end
if nargin < 4 || ~isfinite(plate_t)
    plate_t = 1e-3;
end
legacy_fields = {'k_interface', 't_interface_eff', 'Rc_interface'};
for i = 1:numel(legacy_fields)
    if isfield(raw, legacy_fields{i})
        error('Legacy z_path.%s is not supported. Use vector fields k_interfaces, t_interface_effs, and Rc_interfaces.', legacy_fields{i});
    end
end
enabled = logical_field_ms(raw, 'enable', defaults.enable);
if ~enabled
    z = struct('enable', false);
    return;
end
required_fields = {'k_interfaces', 't_interface_effs', 'Rc_interfaces'};
for i = 1:numel(required_fields)
    if ~isfield(raw, required_fields{i}) || isempty(raw.(required_fields{i}))
        error('z_path.%s is required when z_path.enable is true.', required_fields{i});
    end
end
nif = max(1, N - 1);
z = struct();
z.enable = true;
z.k_interfaces = max(fit_len_vec_ms(raw.k_interfaces, nif, 170 * ones(1, nif)), 1e-9);
z.t_interface_effs = max(0, fit_len_vec_ms(raw.t_interface_effs, nif, plate_t * ones(1, nif)));
z.Rc_interfaces = max(0, fit_len_vec_ms(raw.Rc_interfaces, nif, zeros(1, nif)));
z.k_sink = max(1e-9, scalar_with_default_ms(field_or_default_ms(raw, 'k_sink', 170), 170));
z.t_sink_eff = max(0, scalar_with_default_ms(field_or_default_ms(raw, 't_sink_eff', plate_t), plate_t));
z.Rc_sink = max(0, scalar_with_default_ms(field_or_default_ms(raw, 'Rc_sink', 0), 0));
z.step_fp_iters = max(1, round(scalar_with_default_ms(field_or_default_ms(raw, 'step_fp_iters', 2), 2)));
z.step_fp_relax = min(max(scalar_with_default_ms(field_or_default_ms(raw, 'step_fp_relax', 0.6), 0.6), 0), 1);
z.step_q_prev_weight = min(max(scalar_with_default_ms(field_or_default_ms(raw, 'step_q_prev_weight', 0.85), 0.85), 0), 1);
end

function G = build_shared_G_ms(spec)
G.alpha_coeffs = spec.material.alpha_coeffs;
G.k_coeffs = spec.material.k_coeffs;
G.rho_coeffs = spec.material.rho_coeffs;
G.material = spec.material;
G.alpha_fun = @(T) polyval(G.alpha_coeffs, T);
G.k_leg_fun = @(T) polyval(G.k_coeffs, T);
G.rho_fun = @(T) polyval(G.rho_coeffs, T);
G.fp_w = spec.fp_w;
G.fp_h = spec.fp_h;
G.fp_t_stage = spec.fp_t_stage(:).';
G.fp_t = G.fp_t_stage;
G.fp_A = G.fp_w * G.fp_h;
G.k_plate_fun = @(T) spec.plate_k_inplane + 0*T;
G.A_leg = G.fp_A;
G.L_leg_stage = G.fp_t_stage;
G.L_leg = G.L_leg_stage;
G.Rc = spec.material.Rc;
G.I = spec.current.I_init;
G.ThRes = spec.targets.ThRes;
G.omega = spec.omega;
G.max_inner = spec.max_inner;
G.max_outer = spec.max_outer;
G.tol_theta = spec.tol_theta;
G.tol_g_rel = spec.tol_g_rel;
G.Qc_target_last = spec.targets.Qc_target_last;
end

function out = merge_struct_recursive_ms(base, override)
out = base;
if isempty(override), return; end
fn = fieldnames(override);
for i = 1:numel(fn)
    k = fn{i};
    v = override.(k);
    if isstruct(v) && isfield(out, k) && isstruct(out.(k)) && isscalar(v) && isscalar(out.(k))
        out.(k) = merge_struct_recursive_ms(out.(k), v);
    else
        out.(k) = v;
    end
end
end

function v = resolve_fp_t_stage_ms(spec, fallback, N)
if isfield(spec, 'fp_t_stage') && ~isempty(spec.fp_t_stage)
    raw = spec.fp_t_stage;
elseif isfield(spec, 'fp_t') && ~isempty(spec.fp_t)
    raw = spec.fp_t;
else
    raw = fallback;
end
v = fit_len_vec_ms(raw, N, fallback);
v = max(v, 1e-12);
end

function v = fit_len_vec_ms(v_in, N, fallback)
if nargin < 3 || isempty(fallback), fallback = 1; end
if isempty(v_in), v_in = fallback; end
v_in = double(v_in(:).');
fallback = double(fallback(:).');
if isempty(v_in) || any(~isfinite(v_in))
    v_in = fallback;
end
if numel(v_in) >= N
    v = v_in(1:N);
elseif numel(v_in) == 1
    v = repmat(v_in, 1, N);
else
    v = [v_in, repmat(v_in(end), 1, N - numel(v_in))];
end
end

function v = scalar_with_default_ms(v_in, fallback)
v = fallback;
if isnumeric(v_in) && ~isempty(v_in) && isfinite(v_in(1))
    v = double(v_in(1));
end
end

function v = field_or_default_ms(s, fname, fallback)
if isstruct(s) && isfield(s, fname)
    v = s.(fname);
else
    v = fallback;
end
end

function v = normalize_even_vector_ms(v, N)
v = fit_len_vec_ms(v, N, 2);
v = max(2, 2 * round(v / 2));
end

function total_n = normalize_even_total_count_ms(total_n, N)
total_n = max(2 * N, 2 * round(total_n / 2));
end

function tf = logical_field_ms(s, fname, fallback)
tf = fallback;
if isstruct(s) && isfield(s, fname) && ~isempty(s.(fname))
    tf = any(logical(s.(fname)(:)));
end
end

function vals = normalize_positive_values_ms(raw, fallback)
vals = double(raw(:).');
vals = vals(isfinite(vals) & vals > 0);
if isempty(vals), vals = fallback; end
vals = unique(vals, 'stable');
end

function vals = normalize_numeric_values_ms(raw, fallback)
vals = double(raw(:).');
vals = vals(isfinite(vals));
if isempty(vals), vals = fallback; end
vals = unique(vals, 'stable');
end

function coeffs = normalize_coeffs_ms(raw, fallback)
coeffs = double(raw(:).');
if numel(coeffs) ~= 5 || any(~isfinite(coeffs))
    coeffs = fallback;
end
end

function c = normalize_text_cell_ms(raw, fallback)
if isempty(raw), raw = fallback; end
if ischar(raw) || (isstring(raw) && isscalar(raw))
    c = {char(string(raw))};
elseif isstring(raw)
    c = cellstr(raw(:).');
elseif iscell(raw)
    c = raw(:).';
else
    c = fallback;
end
c = c(~cellfun(@isempty, c));
if isempty(c), c = fallback; end
for i = 1:numel(c)
    c{i} = char(string(c{i}));
end
end

function mode = normalize_plot_view_mode_ms(v)
mode = '2d';
if isempty(v), return; end
s = lower(strtrim(char(string(v))));
if strcmp(s, '3d'), mode = '3d'; end
end

function fig_size = normalize_plot_size_px_ms(v)
fig_size = [900, 900];
if isnumeric(v) && ~isempty(v)
    raw = double(v(:).');
    if numel(raw) == 1, raw = [raw, raw]; end
    if numel(raw) >= 2 && all(isfinite(raw(1:2)))
        fig_size = max(round(raw(1:2)), [300, 300]);
    end
end
end

function ps = normalize_particle_search_ms(ps, defaults, N, total_n)
ps = merge_struct_recursive_ms(defaults, ps);
if isfield(ps, 'n5_max') && ~isempty(ps.n5_max) && isfinite(ps.n5_max)
    ps.n5_max = max(2, 2 * round(ps.n5_max / 2));
else
    ps.n5_max = [];
end
ps.stage_max_count = max(2, 2 * round(scalar_with_default_ms(ps.stage_max_count, defaults.stage_max_count) / 2));
if ~isempty(ps.n5_max), ps.stage_max_count = max(ps.n5_max, ps.stage_max_count); end
ps.total_n_min = max(2 * N, 2 * round(scalar_with_default_ms(ps.total_n_min, defaults.total_n_min) / 2));
ps.total_n_max = max(total_n, 2 * round(scalar_with_default_ms(ps.total_n_max, defaults.total_n_max) / 2));
ps.total_n_max = max(ps.total_n_min, ps.total_n_max);
ps.candidate_budget = max(20, round(scalar_with_default_ms(ps.candidate_budget, defaults.candidate_budget)));
ps.total_block_size = max(1, round(scalar_with_default_ms(ps.total_block_size, defaults.total_block_size)));
ps.par_chunk_size = max(1, round(scalar_with_default_ms(ps.par_chunk_size, defaults.par_chunk_size)));
ps.total_step = max(2, 2 * round(scalar_with_default_ms(ps.total_step, defaults.total_step) / 2));
ps.n5_min_factor = min(max(scalar_with_default_ms(ps.n5_min_factor, defaults.n5_min_factor), 0.1), 1.0);
ps.n5_max_auto_factor = max(1, scalar_with_default_ms(ps.n5_max_auto_factor, defaults.n5_max_auto_factor));
ps.n5_max_auto_margin = max(0, 2 * round(scalar_with_default_ms(ps.n5_max_auto_margin, defaults.n5_max_auto_margin) / 2));
ps.use_ratio_prune = logical_field_ms(ps, 'use_ratio_prune', defaults.use_ratio_prune);
ps.adjacent_min_fraction = min(max(scalar_with_default_ms(ps.adjacent_min_fraction, defaults.adjacent_min_fraction), 0), 1);
ps.boundary_check_count = max(0, round(scalar_with_default_ms(ps.boundary_check_count, defaults.boundary_check_count)));
ps.stop_after_best_margin = max(0, round(scalar_with_default_ms(ps.stop_after_best_margin, defaults.stop_after_best_margin)));
ps.seed_ratios = max(fit_len_vec_ms(ps.seed_ratios, N, defaults.seed_ratios), 1);
ps.mode = char(string(ps.mode));
ps.strict_min_total = logical_field_ms(ps, 'strict_min_total', defaults.strict_min_total);
end

function lks = normalize_low_k_solver_ms(lks, defaults)
lks = merge_struct_recursive_ms(defaults, lks);
lks.enable = logical_field_ms(lks, 'enable', defaults.enable);
lks.k_threshold = max(1e-9, scalar_with_default_ms(lks.k_threshold, defaults.k_threshold));
lks.omega = min(max(scalar_with_default_ms(lks.omega, defaults.omega), 0.05), 1.0);
lks.max_inner = max(1, round(scalar_with_default_ms(lks.max_inner, defaults.max_inner)));
lks.max_outer = max(1, round(scalar_with_default_ms(lks.max_outer, defaults.max_outer)));
lks.continuation_fracs = normalize_numeric_values_ms(lks.continuation_fracs, defaults.continuation_fracs);
lks.continuation_fracs = lks.continuation_fracs(lks.continuation_fracs > 0 & lks.continuation_fracs <= 1);
if isempty(lks.continuation_fracs), lks.continuation_fracs = defaults.continuation_fracs; end
if all(abs(lks.continuation_fracs - 1) > 1e-12), lks.continuation_fracs(end+1) = 1.0; end
end
