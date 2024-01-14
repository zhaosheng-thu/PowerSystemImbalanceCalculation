clear;
mpc = case33bw();

% 设置有功负荷为4140kW
mpc.bus(33, 3) = 4140;

% 使用FDXB方法求解
options_fdxb = mpoption('pf.alg', 'FDXB');
tic;
result_fast_decoupled = runpf(mpc, options_fdxb);
time_fdxb = toc;

% 使用NR方法求解
options_nr = mpoption('pf.alg', 'NR');
tic;
result_newton_raphson = runpf(mpc, options_nr);
time_nr = toc;

% 输出结果和时间
fprintf('NR方法结果:\n');
disp(result_newton_raphson);
fprintf('NR方法执行时间: %.4f秒\n', time_nr);

fprintf('FDXB方法结果:\n');
disp(result_fast_decoupled);
fprintf('FDXB方法执行时间: %.4f秒\n', time_fdxb);

% 输出电压幅值
fprintf('NR方法电压幅值:\n');
disp(result_newton_raphson.bus(:, 8));

fprintf('FDXB方法电压幅值:\n');
disp(result_fast_decoupled.bus(:, 8));
% 绘制NR方法的电压幅值图像
figure;
plot(result_newton_raphson.bus(:, 8), 'ro-', 'LineWidth', 2);
title('NR方法电压幅值');
xlabel('节点编号');
ylabel('电压幅值');
grid on;

% 绘制FDXB方法的电压幅值图像
figure;
plot(result_fast_decoupled.bus(:, 8), 'bo-', 'LineWidth', 2);
title('FDXB方法电压幅值');
xlabel('节点编号');
ylabel('电压幅值');
grid on;

maxIterations = 100;
epsilon = 1e-5;
% 拉格朗日乘子法求解
[PInj, QInj, dPInj, dQInj, dPInj_J, dQInj_J, nodeVoltage, angleDelta, iteration] = LarCalculatePoweImbalance(mpc, maxIterations, epsilon);
