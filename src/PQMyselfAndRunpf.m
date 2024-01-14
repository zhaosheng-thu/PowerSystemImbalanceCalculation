clear
mpc = case39();
% mpc = case33bw();
maxIterations = 50;
epsilon = 1e-5;

% 运行迭代算法
tic;  % 记录开始时间
[angleData, nodeVoltage, iteration] = PQCalculatePoweImbalance(mpc, maxIterations, epsilon);
toc;  % 记录结束时间

% 将角度从弧度转换为度
angleData = rad2deg(angleData);

% 输出迭代总次数
fprintf('迭代总次数：%d\n', iteration);

% FDXB算法
options = mpoption('pf.alg', 'FDXB', 'pf.tol', epsilon);
tic;
results_FDXB = runpf(mpc, options);
toc;

disp(['FDXB次数：', num2str(results_FDXB.iterations)]);
disp(['FDXB时间：', num2str(results_FDXB.et)]);

% FDBX算法
options = mpoption('pf.alg', 'FDBX', 'pf.tol', epsilon);
tic;
results_FDBX = runpf(mpc, options);
toc;

disp(['FDBX次数：', num2str(results_FDBX.iterations)]);
disp(['FDBX时间：', num2str(results_FDBX.et)]);
