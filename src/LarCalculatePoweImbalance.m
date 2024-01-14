function [PInj, QInj, nodeVoltage, angleDelta, iteration] = LarCalculatePoweImbalance(mpc, maxIterations, epsilon)
    baseMVA = mpc.baseMVA;
    busData = mpc.bus;
    genData = mpc.gen;
    branchData = mpc.branch;
    nodeTypes = mpc.bus(:, 2);
    numNodes = length(nodeTypes);

    % 导纳矩阵
    Ybus = makeYbus(mpc);
    G = real(Ybus);
    B = imag(Ybus);

    % 节点电压和相角
    nodeVoltage = busData(:, 8)';
    angleDelta = busData(:, 9)';
    angleDelta = deg2rad(angleDelta);

    PGen = zeros(numNodes, 1);
    QGen = zeros(numNodes, 1);

    for i = 1:length(genData(:, 1))
        PGen(genData(i, 1)) = genData(i, 2);
        QGen(genData(i, 1)) = genData(i, 3);
    end

    % 初始潮流计算
    [PInj, QInj, ~, ~, ~, ~, ~, ~] = calculatePowerInjections(G, B, nodeVoltage, angleDelta, PGen, QGen, baseMVA, nodeTypes);

    % 迭代求解潮流
    iteration = 1;

    while iteration < maxIterations
        % 求解梯度和雅可比矩阵
        [dPInj, dQInj, dPInj_J, dQInj_J, ~, ~, ~, ~] = calculatePowerInjections(G, B, nodeVoltage, angleDelta, PGen, QGen, baseMVA, nodeTypes);

        % 判断是否收敛
        if (max(abs(dQInj(nodeTypes == 1))) <= epsilon) && (max(abs(dPInj(nodeTypes ~= 3))) <= epsilon)
            disp('潮流收敛');
            break
        end

        % 更新节点电压和相角
        [nodeVoltage, angleDelta] = updateVoltageAngles(nodeVoltage, angleDelta, dPInj_J, dQInj_J);
        iteration = iteration + 1;
    end
end

function [PInj, QInj, dPInj, dQInj, dPInj_J, dQInj_J, nodeVoltage, angleDelta] = calculatePowerInjections(G, B, nodeVoltage, angleDelta, PGen, QGen, baseMVA, nodeTypes)
    numNodes = length(nodeTypes);
    indexNot3 = nodeTypes ~= 3;
    indexIs1 = nodeTypes == 1;

    % 计算有功和无功注入功率
    PLine = busData(:, 3);
    QLine = busData(:, 4);
    PInj = (PGen - PLine) / baseMVA;
    QInj = (QGen - QLine) / baseMVA;

    % 计算雅可比矩阵
    [dPInj, dQInj, dPInj_J, dQInj_J] = calculateJacobian(G, B, nodeVoltage, angleDelta, PGen, QGen, baseMVA, nodeTypes);

    % 更新节点电压和相角
    [nodeVoltage, angleDelta] = updateVoltageAngles(nodeVoltage, angleDelta, dPInj_J, dQInj_J);
end

function [dPInj, dQInj, dPInj_J, dQInj_J] = calculateJacobian(G, B, nodeVoltage, angleDelta, PGen, QGen, baseMVA, nodeTypes)
    numNodes = length(nodeTypes);
    indexNot3 = find(nodeTypes ~= 3);
    indexIs1 = find(nodeTypes == 1);

    % 计算梯度
    [dPInj, dQInj] = calculatePowerImbalance(G, B, nodeVoltage, angleDelta, PGen, QGen, baseMVA, nodeTypes);

    % 计算雅可比矩阵
    dPInj_J = dPInj(indexNot3);
    dQInj_J = dQInj(indexIs1);
end

function [dPInj, dQInj] = calculatePowerImbalance(G, B, nodeVoltage, angleDelta, PGen, QGen, baseMVA, nodeTypes)
    numNodes = length(nodeTypes);
    indexNot3 = find(nodeTypes ~= 3);
    indexIs1 = find(nodeTypes == 1);

    % 计算无功不平衡
    QInj = zeros(1, numNodes);
    for i = 1:numNodes
        for j = 1:numNodes
            QInjM(j) = nodeVoltage(i) * nodeVoltage(j) * (G(i, j) * sin(angleDelta(i) - angleDelta(j)) - B(i, j) * cos(angleDelta(i) - angleDelta(j)));
        end
        QInjSum(i) = sum(QInjM);
    end
    dQInj = QGen - QInjSum;

    % 计算有功不平衡
    PInj = zeros(1, numNodes);
    for i = 1:numNodes
        for j = 1:numNodes
            PInjH(j) = nodeVoltage(i) * nodeVoltage(j) * (G(i, j) * cos(angleDelta(i) - angleDelta(j)) + B(i, j) * sin(angleDelta(i) - angleDelta(j)));
        end
        PInjSum(i) = sum(PInjH);
    end
    dPInj = PGen - PInjSum;
end

function [nodeVoltage, angleDelta] = updateVoltageAngles(nodeVoltage, angleDelta, dPInj_J, dQInj_J)
    % 更新节点电压和相角
    PQNode = find(nodeTypes ~= 3);
    PVNode = find(nodeTypes == 1);

    % 雅可比矩阵的计算
    J = [dPInj_J; dQInj_J];

    % 节点电压修正量
    dUdelta = (-inv(J) * [dPInj_J dQInj_J]')';
    dDelta = dUdelta(1:length(PQNode));
    dU = (nodeVoltage(PQNode) .* dUdelta(length(PQNode)+1:end));
    
    % 更新节点电压和相角
    nodeVoltage = nodeVoltage + dU;
    angleDelta(PQNode) = angleDelta(PQNode) + dDelta;
end
