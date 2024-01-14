clear;
% 读取数据
mpc = case39();
maxIterations = 50;
epsilon = 1e-4;

% 调用主函数
[PInj, QInj, dPInj, dQInj, dPInj_J, dQInj_J, nodeVoltage, angleDelta, iteration] = NRCalculatePowerImbalance(mpc, maxIterations, epsilon);

% 打印结果
fprintf('迭代总次数：%d\n', iteration);
result = runpf(mpc);
U = result.bus(:, 2);
voltageDifference = nodeVoltage' - U;
angleDelta1 = result.bus(:, 3);
angleDelta = rad2deg(angleDelta);
angleDifference = angleDelta' - angleDelta1;

disp('电压相角差值：');
disp(angleDifference);

disp('电压幅值差值（保留六位小数）：');
disp(num2str(voltageDifference, '%.6f'));