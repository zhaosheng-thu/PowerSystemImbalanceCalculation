function [PInj, QInj, dPInj, dQInj, dPInj_J, dQInj_J, nodeVoltage, angleDelta, iteration] = NRCalculatePowerImbalance(mpc, maxIterations, epsilon)
    baseMVA = mpc.baseMVA;
    busData = mpc.bus;
    genData = mpc.gen;
    branchData = mpc.branch;
    nodeTypes = mpc.bus(:, 2);
    numNodes = length(nodeTypes); % 总节点数
    indexNot3 = nodeTypes ~= 3;
    indexIs1 = nodeTypes == 1;
    % 导纳矩阵
    Ybus = makeYbus(mpc);
    G = real(Ybus);
    B = imag(Ybus);

    % 节点电压初始值
    PLine = busData(:, 3);
    QLine = busData(:, 4);
    PGen = zeros(numNodes, 1);
    QGen = zeros(numNodes, 1);

    % 节点电压和相角
    nodeVoltage = busData(:, 8)';
    angleDelta = busData(:, 9)';
    angleDelta = deg2rad(angleDelta);

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
        % 求雅可比矩阵
        [dPInj, dQInj, dPInj_J, dQInj_J, nodeVoltage, angleDelta] = calculateJacobian(PInj, QInj, G, B, nodeVoltage, angleDelta, nodeTypes);

        % 判断是否收敛
        if (max(abs(dQInj(indexIs1))) <= epsilon) && (max(abs(dPInj(indexNot3))) <= epsilon)
            disp('潮流收敛');
            break
        end
        iteration = iteration + 1;
    end

    
end

function [dPInj, dQInj, dPInj_J, dQInj_J, nodeVoltage, angleDelta] = calculateJacobian(PInj, QInj, G, B, nodeVoltage, angleDelta, nodeTypes)
    numNodes = length(nodeTypes); % 总节点数
    indexNot3 = find(nodeTypes ~= 3);
    indexIs1 = find(nodeTypes == 1);
    % 无功不平衡量
    for i = 1:numNodes
        for j = 1:numNodes
            QInjM(j) = nodeVoltage(i) * nodeVoltage(j) * (G(i, j) * sin(angleDelta(i) - angleDelta(j)) - B(i, j) * cos(angleDelta(i) - angleDelta(j)));
        end
        QInjSum(i) = sum(QInjM);
    end
    dQInj = QInj - QInjSum;

    % 有功不平衡
    for i = 1:numNodes
        for j = 1:numNodes
            PInjH(j) = nodeVoltage(i) * nodeVoltage(j) * (G(i, j) * cos(angleDelta(i) - angleDelta(j)) + B(i, j) * sin(angleDelta(i) - angleDelta(j)));
        end
        PInjSum(i) = sum(PInjH);
    end
    dPInj = PInj - PInjSum;

    % 雅可比矩阵
    dPInj_J = dPInj(indexNot3);
    dQInj_J = dQInj(indexIs1);

    % 迭代求解节点电压修正量
    PQNode = indexNot3;
    PVNode = indexIs1;
    
    % 雅可比矩阵的计算
    [H, N, K, L] = calculateJacobianBlocks(nodeVoltage, angleDelta, G, B, PQNode, PVNode, numNodes, PInjSum, QInjSum, nodeTypes);
    J = [H N; K L]; 

    % 节点电压修正量 (case39)
    dPQ = [dPInj_J dQInj_J]';
    dUdelta = (-inv(J) * dPQ)';
    dDelta = dUdelta(1:38);
    dDelta = [dDelta(1:30) 0 dDelta(31:38)];
    dU = (nodeVoltage(1:29) .* dUdelta(39:67));
    dU = [dU zeros(1, 10)];
    nodeVoltage = nodeVoltage + dU;
    angleDelta = angleDelta + dDelta;
end

function [H, N, K, L] = calculateJacobianBlocks(nodeVoltage, angleDelta, G, B, PQNode, PVNode, ~, PInjSum, QInjSum, nodeTypes)
    numNodes = length(nodeTypes);
    numIs1 = sum(nodeTypes == 1);
    numIs3 = sum(nodeTypes == 3);
    numNodes = numNodes - numIs3;
    H = zeros(numNodes, numNodes);
    N = zeros(numNodes, numIs1);
    K = zeros(numIs1, numNodes);
    L = zeros(numIs1, numIs1);
    
    % i==j时
    for i = 1:numIs1
        N(i, i) = -nodeVoltage(i).^2 * G(i, i) - PInjSum(i);
    end
    
    for i = 1:numIs1
        K(i, i) = nodeVoltage(i).^2 * G(i, i) - PInjSum(i);
    end
    
    for i = 1:numIs1
        L(i, i) = nodeVoltage(i).^2 * B(i, i) - QInjSum(i);
    end

    for i = 1:numNodes
        H(i, i) = nodeVoltage(PQNode(i)).^2 * B(PQNode(i), PQNode(i)) + QInjSum(PQNode(i));
    end
    
    % i!=j时
    for i = 1:numNodes
        for j = 1:numNodes
            H(i, j) = -nodeVoltage(PQNode(i)) * nodeVoltage(PQNode(j)) * (G(PQNode(i), PQNode(j)) * sin(angleDelta(PQNode(i)) - angleDelta(PQNode(j))) - B(PQNode(i), PQNode(j)) * cos(angleDelta(PQNode(i)) - angleDelta(PQNode(j))));
        end
    end

    for i = 1:numNodes
        for j = 1:numIs1
            N(i, j) = -nodeVoltage(PQNode(i)) * nodeVoltage(PVNode(j)) * (G(PQNode(i), PVNode(j)) * cos(angleDelta(PQNode(i)) - angleDelta(PVNode(j))) + B(PQNode(i), PVNode(j)) * sin(angleDelta(PQNode(i)) - angleDelta(PVNode(j))));
        end
    end
    
    for i = 1:numIs1
        for j = 1:numNodes
            K(i, j) = nodeVoltage(PVNode(i)) * nodeVoltage(PQNode(j)) * (G(PVNode(i), PQNode(j)) * cos(angleDelta(PVNode(i)) - angleDelta(PQNode(j))) + B(PVNode(i), PQNode(j)) * sin(angleDelta(PVNode(i)) - angleDelta(PQNode(j))));
        end
    end
    
    for i = 1:numIs1
        for j = 1:numIs1
            L(i, j) = -nodeVoltage(PVNode(i)) * nodeVoltage(PVNode(j)) * (G(PVNode(i), PVNode(j)) * sin(angleDelta(PVNode(i)) - angleDelta(PVNode(j))) - B(PVNode(i), PVNode(j)) * cos(angleDelta(PVNode(i)) - angleDelta(PVNode(j))));
        end
    end
    
end
