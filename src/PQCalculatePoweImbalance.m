function [angleDelta, nodeVoltage, iteration] = PQCalculatePoweImbalance(mpc, maxIterations, epsilon)
    baseMVA = mpc.baseMVA;
    busData = mpc.bus;
    genData = mpc.gen;
    branchData = mpc.branch;
    nodeTypes = busData(:, 2);
    numNodes = length(nodeTypes);
    indexNot3 = nodeTypes ~= 3;
    indexIs1 = nodeTypes == 1;

    % 求导纳矩阵
    YBus = makeYbus(mpc);
    G = real(YBus);
    B = imag(YBus);

    PLine = busData(:, 3);
    QLine = busData(:, 4);
    PGen = zeros(numNodes, 1);
    QGen = zeros(numNodes, 1);

    nodeVoltage = busData(:, 8)';
    angleDelta = busData(:, 9)';
    angleDelta = (deg2rad(angleDelta));

    for i = 1:length(genData(:, 1))
        PGen(genData(i, 1)) = genData(i, 2);
        QGen(genData(i, 1)) = genData(i, 3);
    end
    PInj = zeros(1, numNodes);
    QInj = zeros(1, numNodes);
    for i = 1:numNodes
        PInj(i) = PGen(i) - PLine(i);
        QInj(i) = QGen(i) - QLine(i);
    end
    PInj = PInj / baseMVA;
    QInj = QInj / baseMVA;

    % 迭代求解潮流
    iteration = 1;
    
    while iteration < maxIterations
        [dPInj, dQInj, ~, ~, nodeVoltage, angleDelta] = calculateJacobian(PInj, QInj, G, B, nodeVoltage, angleDelta, nodeTypes);
        if (max(abs(dQInj(indexIs1))) <= epsilon) && (max(abs(dPInj(indexNot3))) <= epsilon)
            disp('潮流计算收敛');
            break
        end
        iteration = iteration + 1;
    end
    % angleDelta = rad2deg(angleDelta);
end

function [dPInj, dQInj, dPInj_J, dQInj_J, nodeVoltage, angleDelta] = calculateJacobian(PInj, QInj, G, B, nodeVoltage, angleDelta, nodeTypes)
    numNodes = length(nodeTypes); % 总节点数
    indexNot3 = find(nodeTypes ~= 3);
    indexIs1 = find(nodeTypes == 1);
    for i = 1:numNodes
        for j = 1:numNodes
            PInjM(j) = nodeVoltage(i) * nodeVoltage(j) * (G(i, j) * cos(angleDelta(i) - angleDelta(j)) + B(i, j) * sin(angleDelta(i) - angleDelta(j)));
        end
        PInjSum(i) = sum(PInjM);
    end
    dPInj = PInj - PInjSum;

    for i = 1:numNodes
        for j = 1:numNodes
            QInjM(j) = nodeVoltage(i) * nodeVoltage(j) * (G(i, j) * sin(angleDelta(i) - angleDelta(j)) - B(i, j) * cos(angleDelta(i) - angleDelta(j)));
        end
        QInjSum(i) = sum(QInjM);
    end
    dQInj = QInj - QInjSum;
    % 雅可比矩阵
    dPInj_J = dPInj(indexNot3);
    dQInj_J = dQInj(indexIs1);

    % 迭代求解节点电压修正量
    PQNode = indexNot3;
    PVNode = indexIs1;

    for i = 1:numNodes - 1
        for j = 1:numNodes - 1
            B1(i, j) = B(indexNot3(i), indexNot3(j));
        end
    end

    for i = 1:length(indexIs1)
        for j = 1:length(indexIs1)
            B2(i, j) = B(indexIs1(i), indexIs1(j));
        end
    end
    % case 39
    d_delta = -inv(B1) * (dPInj_J ./ [nodeVoltage(1:30) nodeVoltage(32:39)])';
    d_delta = d_delta';
    angleDelta(1:30) = angleDelta(1:30) + d_delta(1:30);
    angleDelta(32:39) = angleDelta(32:39) + d_delta(31:38);
    dU = -inv(B2) * (dQInj_J ./ nodeVoltage(1:29))';
    nodeVoltage(1:29) = nodeVoltage(1:29) + dU';
end
