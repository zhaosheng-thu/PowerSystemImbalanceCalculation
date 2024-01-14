## PowerSystemImbalanceCalculation
 $$
PowerSystem (20220543-0) Program Homework 1
$$ 
- 选取了`matpower`的case39和case33的测例，实现了使用`matpower`库函数以及手动编写的潮流计算matlab代码，包括了NR法和PQ分解法
- 实现了对于修改一些节点的 PQVδ 等数据后的网络的收敛迭代
- 若收敛性差，可以使用拉格朗日乘子法做优化（目前有bug）