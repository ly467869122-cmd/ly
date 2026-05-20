function params = default_params()
% DEFAULT_PARAMS
% 多级热电制冷器件统一参数配置文件
% 后续所有主程序都从这里读取参数，避免参数散落在各个脚本中

%% 基本目标参数
params.target.DeltaT = 110;        % 目标制冷温差 K
params.target.Qc_last = 0.8;       % 最高级冷端制冷量 W
params.target.Th_res = 300;        % 热端热沉温度 K

%% 器件级数
params.stage.count = 5;

%% 电流参数
params.current.I_min = 1.0;        % 最小电流 A
params.current.I_max = 6.0;        % 最大电流 A
params.current.I_step = 0.1;       % 电流扫描步长 A
params.current.I_default = 4.0;    % 默认电流 A

%% 热电腿几何参数
params.leg.A_leg = 1e-6;           % 单腿截面积 m^2, 1mm × 1mm
params.leg.L_leg = 1e-3;           % 腿长 m, 1mm
params.leg.Rc = 0;                 % 接触电阻，先取0

%% 基板参数
params.substrate.k_AlN = 170;      % AlN热导率 W/(m·K)
params.substrate.thickness = 1e-3; % 基板厚度 m

%% FEM参数
params.fem.mesh_size_mm = 0.5;     % 网格尺寸 mm
params.fem.tol_theta = 1e-3;       % 内层温度场收敛阈值 K
params.fem.max_iter_inner = 100;   % 内层最大迭代次数
params.fem.max_iter_outer = 80;    % 外层最大迭代次数
params.fem.damping = 0.5;          % 阻尼系数

%% 优化参数
params.opt.topK = 10;              % 保留前10个候选
params.opt.objective = 'DeltaT';   % 当前目标：温差最小
params.opt.use_parallel = true;    % 是否并行计算

%% 输出参数
params.output.save_figures = true;
params.output.save_mat = true;
params.output.save_log = true;

end