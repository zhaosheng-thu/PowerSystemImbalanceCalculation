% 导入39节点系统数据
define_constants;
system_data = case39;

% 获取系统数据并初始化节点导纳矩阵 Y
bus_data = system_data.bus;
line_data = system_data.branch;
num_buses = length(bus_data);
num_lines = length(line_data);
Y = zeros(num_buses, num_buses);

% 构建节点导纳矩阵 Y
for k = 1:num_lines
    from_bus = line_data(k, F_BUS);
    to_bus = line_data(k, T_BUS);
    
    % 提取支路参数
    resistance = line_data(k, BR_R);
    reactance = line_data(k, BR_X);
    susceptance = line_data(k, BR_B);

    % 线路导纳
    impedance = complex(resistance, reactance);
    admittance = 1 / impedance;

       
    % 变压器
    tap_ratio = line_data(k, TAP);
    phase_shift = deg2rad(line_data(k, SHIFT));

    if susceptance ~= 0
        admittance = admittance + 1i * susceptance / 2;
    end

    % 考虑变压器变比和移相器
    if tap_ratio ~= 0
        admittance = admittance / (tap_ratio * exp(1i * phase_shift));
    end

    % 自导纳
    Y(from_bus, from_bus) = Y(from_bus, from_bus) + admittance;
    Y(to_bus, to_bus) = Y(to_bus, to_bus) + admittance;
    
    % 互导纳
    Y(from_bus, to_bus) = Y(from_bus, to_bus) - admittance;
    Y(to_bus, from_bus) = Y(to_bus, from_bus) - admittance;
end

% makeYbus 函数
[Ybus, ~, ~] = makeYbus(system_data.baseMVA, system_data.bus, system_data.branch);

difference = Y - Ybus;

% 计算误差
absolute_error = norm(difference, 'fro');
rmse = sqrt(mean(difference(:).^2));

disp(['绝对误差为: ' num2str(absolute_error)]);
disp(['均方根误差为: ' num2str(rmse)]);

% 可视化误差矩阵的热力图
figure;
heatmap(abs(difference), 'Colormap', parula, 'ColorbarVisible', 'on');
title('误差矩阵的热力图');
xlabel('节点编号');
ylabel('节点编号');
