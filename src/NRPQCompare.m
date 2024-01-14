% 使用Matpower读取case39数据
mpc = case39();
baseMVA = mpc.baseMVA;

% 使用快速分解法进行潮流计算
options = mpoption('pf.alg', 'FDXB');
tic;
result_fast_decoupled = runpf(mpc, options);
time_fast_decoupled = toc;

% 使用牛顿拉夫逊法进行潮流计算
options = mpoption('pf.alg', 'NR');
tic;
result_newton_raphson = runpf(mpc, options);
time_newton_raphson = toc;

% 输出结果
disp('快速分解法结果：');
disp(result_fast_decoupled);

disp('牛顿拉夫逊法结果：');
disp(result_newton_raphson);

% 输出迭代信息
disp('快速分解法迭代信息：');
disp(result_fast_decoupled.iterations);
disp(['计算时间：', num2str(time_fast_decoupled), '秒']);

disp('牛顿拉夫逊法迭代信息：');
disp(result_newton_raphson.iterations);
disp(['计算时间：', num2str(time_newton_raphson), '秒']);
