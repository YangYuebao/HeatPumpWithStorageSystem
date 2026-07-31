"""
通用动态规划求解器核心模块

提供与具体问题无关的动态规划求解算法，包括：
- forwardSolve: 单步正向递推
- dpSolve: 给定初始温度求解最优路径
- ExhaustiveSolver: 穷举搜索所有初始温度
- GoldenRatioSolver: 黄金分割搜索最优初始温度
"""
module DPSolverCore

#export forwardSolve, dpSolve, ExhaustiveSolver, GoldenRatioSolver

"""
正向递推求解单时间层

参数:
- VForward: 当前时间层的最优成本向量 [nT]
- C: 当前时间层的状态转移成本矩阵 [nT, nT]，C[i,j] 表示从状态 i 转移到状态 j 的成本
- nT: 状态数

返回:
- VForwardNew: 下一时刻的最优成本向量 [nT]
- lastTsIndex: 记录每个状态的前驱状态索引 [nT]
"""
function forwardSolve(VForward::Vector{Float64}, C::Matrix{Float64}, nT::Int)
	VForwardNew = fill(Inf, nT)
	lastTsIndex = zeros(Int, nT)

	for j in 1:nT
		for i in 1:nT
			if VForward[i] + C[i, j] < VForwardNew[j]
				VForwardNew[j] = VForward[i] + C[i, j]
				lastTsIndex[j] = i
			end
		end
	end

	return VForwardNew, lastTsIndex
end

"""
正向递推求解单时间层（带边界限制优化版）

参数:
- VForward: 当前时间层的最优成本向量 [nT]
- C: 当前时间层的状态转移成本矩阵 [nT, nT]
- nT: 状态数
- TsDecreaseIndex: 温度下降最多偏移的索引数
- TsIncreaseIndex: 温度上升最多偏移的索引数

返回:
- VForwardNew: 下一时刻的最优成本向量 [nT]
- lastTsIndex: 记录每个状态的前驱状态索引 [nT]

适用场景: 当状态转移存在物理约束时（如温度变化幅度有限），可以限制搜索范围
"""
function forwardSolve(VForward::Vector{Float64}, C::Matrix{Float64}, nT::Int,
	TsDecreaseIndex::Int, TsIncreaseIndex::Int)
	VForwardNew = fill(Inf, nT)
	lastTsIndex = zeros(Int, nT)

	for j in 1:nT
		start_i = max(1, j - TsIncreaseIndex)
		end_i = min(nT, j + TsDecreaseIndex)

		for i in start_i:end_i
			if VForward[i] + C[i, j] < VForwardNew[j]
				VForwardNew[j] = VForward[i] + C[i, j]
				lastTsIndex[j] = i
			end
		end
	end

	return VForwardNew, lastTsIndex
end

"""
给定初始温度索引求解最优路径

参数:
- C: 状态转移成本张量 [nt, nT, nT]，C[i, :, :] 表示第 i 个时间段的转移成本矩阵
- start_idx: 初始温度状态索引（1-based）
- nT: 状态数
- nt: 时间段数

返回:
- total_cost: 总最优成本
- TsIndexList: 各时间层的状态索引列表 [nt+1]，首尾相同（周期边界）

注意:
- 假设周期边界条件：第 nt+1 个时刻的状态必须等于第 1 个时刻的状态
- C 的维度：[时间段数, 起始状态数, 结束状态数]
"""
function dpSolve(C::Array{Float64, 3}, start_idx::Int, nT::Int, nt::Int)
	TsTransitionMatrix = zeros(Int, nT, nt)

	VForward = C[1, start_idx, :]
	TsTransitionMatrix[:, 1] .= start_idx

	for i in 2:nt
		VForward, TsTransitionMatrix[:, i] = forwardSolve(VForward, C[i, :, :], nT)
	end

	total_cost = VForward[start_idx]

	TsIndexList = Vector{Int}(undef, nt + 1)
	TsIndexList[nt+1] = start_idx

	for i in nt:-1:2
		TsIndexList[i] = TsTransitionMatrix[TsIndexList[i+1], i]
	end
	TsIndexList[1] = start_idx

	return total_cost, TsIndexList
end

"""
给定初始温度索引求解最优路径（带边界限制版）

参数:
- C: 状态转移成本张量 [nt, nT, nT]
- start_idx: 初始温度状态索引
- nT: 状态数
- nt: 时间段数
- TsDecreaseIndex: 温度下降最多偏移的索引数
- TsIncreaseIndex: 温度上升最多偏移的索引数

返回:
- total_cost: 总最优成本
- TsIndexList: 各时间层的状态索引列表 [nt+1]
"""
function dpSolve(C::Array{Float64, 3}, start_idx::Int, nT::Int, nt::Int,
	TsDecreaseIndex::Int, TsIncreaseIndex::Int)
	TsTransitionMatrix = zeros(Int, nT, nt)

	VForward = C[1, start_idx, :]
	TsTransitionMatrix[:, 1] .= start_idx

	for i in 2:nt
		VForward, TsTransitionMatrix[:, i] = forwardSolve(VForward, C[i, :, :], nT,
			TsDecreaseIndex, TsIncreaseIndex)
	end

	total_cost = VForward[start_idx]

	TsIndexList = Vector{Int}(undef, nt + 1)
	TsIndexList[nt+1] = start_idx

	for i in nt:-1:2
		TsIndexList[i] = TsTransitionMatrix[TsIndexList[i+1], i]
	end
	TsIndexList[1] = start_idx

	return total_cost, TsIndexList
end

"""
穷举搜索所有初始温度，找到最优解

参数:
- C: 状态转移成本张量 [nt, nT, nT]
- nT: 状态数
- nt: 时间段数

返回:
- min_cost: 最小总成本
- min_TsIndexList: 最优路径的状态索引列表 [nt+1]
- best_start_idx: 最优初始温度索引

适用场景: 状态数较少时（如 nT < 200），穷举搜索可以保证找到全局最优解
"""
function ExhaustiveSolver(C::Array{Float64, 3}, nT::Int, nt::Int)
	valueList = zeros(nT)
	TsIndexListMatrix = zeros(Int, nt + 1, nT)

	min_cost = Inf
	best_start_idx = 0

	for j in 1:nT
		cost, TsIndexList = dpSolve(C, j, nT, nt)
		valueList[j] = cost
		TsIndexListMatrix[:, j] = TsIndexList

		if cost < min_cost
			min_cost = cost
			best_start_idx = j
		end
	end

	min_TsIndexList = TsIndexListMatrix[:, best_start_idx]

	return min_cost, min_TsIndexList, best_start_idx
end

"""
穷举搜索所有初始温度（带边界限制版）

参数:
- C: 状态转移成本张量 [nt, nT, nT]
- nT: 状态数
- nt: 时间段数
- TsDecreaseIndex: 温度下降最多偏移的索引数
- TsIncreaseIndex: 温度上升最多偏移的索引数

返回:
- min_cost: 最小总成本
- min_TsIndexList: 最优路径的状态索引列表 [nt+1]
- best_start_idx: 最优初始温度索引
"""
function ExhaustiveSolver(C::Array{Float64, 3}, nT::Int, nt::Int,
	TsDecreaseIndex::Int, TsIncreaseIndex::Int)
	valueList = zeros(nT)
	TsIndexListMatrix = zeros(Int, nt + 1, nT)

	min_cost = Inf
	best_start_idx = 0

	for j in 1:nT
		cost, TsIndexList = dpSolve(C, j, nT, nt, TsDecreaseIndex, TsIncreaseIndex)
		valueList[j] = cost
		TsIndexListMatrix[:, j] = TsIndexList

		if cost < min_cost
			min_cost = cost
			best_start_idx = j
		end
	end

	min_TsIndexList = TsIndexListMatrix[:, best_start_idx]

	return min_cost, min_TsIndexList, best_start_idx
end

"""
黄金分割搜索最优初始温度

参数:
- C: 状态转移成本张量 [nt, nT, nT]
- nT: 状态数
- nt: 时间段数

返回:
- min_cost: 最小总成本
- min_TsIndexList: 最优路径的状态索引列表 [nt+1]
- best_start_idx: 最优初始温度索引

适用场景: 状态数较多时（如 nT > 200），黄金分割搜索可以显著减少计算量
注意:
- 假设目标函数关于初始温度索引是单峰的（或近似单峰）
- 如果目标函数不是单峰的，可能只能找到局部最优
"""
function GoldenRatioSolver(C::Array{Float64, 3}, nT::Int, nt::Int)
	phi = 0.618

	jList = [1, nT - round(Int, phi * (nT - 1)), 1 + round(Int, phi * (nT - 1)), nT]
	valueList = zeros(4)
	TsIndexListMatrix = zeros(Int, nt + 1, 4)

	for (i, j) in enumerate(jList)
		valueList[i], TsIndexListMatrix[:, i] = dpSolve(C, j, nT, nt)
	end

	count = 0
	while (jList[4] - jList[1] >= 4) && count < 100
		if valueList[2] < valueList[3]
			jList[4] = jList[3]
			valueList[4] = valueList[3]
			TsIndexListMatrix[:, 4] = TsIndexListMatrix[:, 3]

			jList[3] = jList[2]
			valueList[3] = valueList[2]
			TsIndexListMatrix[:, 3] = TsIndexListMatrix[:, 2]

			jList[2] = min(max(floor(Int, jList[4] - phi * (jList[4] - jList[1])), jList[1] + 1), jList[3] - 1)
			valueList[2], TsIndexListMatrix[:, 2] = dpSolve(C, jList[2], nT, nt)
		else
			jList[1] = jList[2]
			valueList[1] = valueList[2]
			TsIndexListMatrix[:, 1] = TsIndexListMatrix[:, 2]

			jList[2] = jList[3]
			valueList[2] = valueList[3]
			TsIndexListMatrix[:, 2] = TsIndexListMatrix[:, 3]

			jList[3] = max(min(ceil(Int, jList[1] + phi * (jList[4] - jList[1])), jList[4] - 1), jList[2] + 1)
			valueList[3], TsIndexListMatrix[:, 3] = dpSolve(C, jList[3], nT, nt)
		end
		count += 1
	end

	min_cost, index = findmin(valueList)
	min_TsIndexList = TsIndexListMatrix[:, index]
	best_start_idx = jList[index]

	return min_cost, min_TsIndexList, best_start_idx
end

"""
黄金分割搜索最优初始温度（带边界限制版）

参数:
- C: 状态转移成本张量 [nt, nT, nT]
- nT: 状态数
- nt: 时间段数
- TsDecreaseIndex: 温度下降最多偏移的索引数
- TsIncreaseIndex: 温度上升最多偏移的索引数

返回:
- min_cost: 最小总成本
- min_TsIndexList: 最优路径的状态索引列表 [nt+1]
- best_start_idx: 最优初始温度索引
"""
function GoldenRatioSolver(C::Array{Float64, 3}, nT::Int, nt::Int,
	TsDecreaseIndex::Int, TsIncreaseIndex::Int)
	phi = 0.618

	jList = [1, nT - round(Int, phi * (nT - 1)), 1 + round(Int, phi * (nT - 1)), nT]
	valueList = zeros(4)
	TsIndexListMatrix = zeros(Int, nt + 1, 4)

	for (i, j) in enumerate(jList)
		valueList[i], TsIndexListMatrix[:, i] = dpSolve(C, j, nT, nt,
			TsDecreaseIndex, TsIncreaseIndex)
	end

	count = 0
	while (jList[4] - jList[1] >= 4) && count < 100
		if valueList[2] < valueList[3]
			jList[4] = jList[3]
			valueList[4] = valueList[3]
			TsIndexListMatrix[:, 4] = TsIndexListMatrix[:, 3]

			jList[3] = jList[2]
			valueList[3] = valueList[2]
			TsIndexListMatrix[:, 3] = TsIndexListMatrix[:, 2]

			jList[2] = min(max(floor(Int, jList[4] - phi * (jList[4] - jList[1])), jList[1] + 1), jList[3] - 1)
			valueList[2], TsIndexListMatrix[:, 2] = dpSolve(C, jList[2], nT, nt,
				TsDecreaseIndex, TsIncreaseIndex)
		else
			jList[1] = jList[2]
			valueList[1] = valueList[2]
			TsIndexListMatrix[:, 1] = TsIndexListMatrix[:, 2]

			jList[2] = jList[3]
			valueList[2] = valueList[3]
			TsIndexListMatrix[:, 2] = TsIndexListMatrix[:, 3]

			jList[3] = max(min(ceil(Int, jList[1] + phi * (jList[4] - jList[1])), jList[4] - 1), jList[2] + 1)
			valueList[3], TsIndexListMatrix[:, 3] = dpSolve(C, jList[3], nT, nt,
				TsDecreaseIndex, TsIncreaseIndex)
		end
		count += 1
	end

	min_cost, index = findmin(valueList)
	min_TsIndexList = TsIndexListMatrix[:, index]
	best_start_idx = jList[index]

	return min_cost, min_TsIndexList, best_start_idx
end

"""
带s6最小间隔约束的图DP求解器

在标准DP基础上增加第二维状态——"距离上次状态6的时间"（s6_wait）。
状态空间: (Ts_index, s6_wait_index)，其中 s6_wait ∈ {0, dt, 2dt, ..., y⁺}

转移规则:
- 不进入s6: (i, w) → (j, w+dt)，成本 = Cn[t,i,j]（w达到y时压缩为y⁺）
- 进入s6: (i, y⁺) → (j, 0)，成本 = Cs[t,i,j]（仅当 w == y⁺ 时可用）

周期性边界条件:
  Ts[1] == Ts[nt+1] 且 wait[1] == wait[nt+1]
  通过枚举所有 (start_idx, w0) 组合，对每组单独跑DP，最后选成本最低的

参数:
- Cn: 不含s6的转移成本张量 [nt, nT, nT]
- Cs: 仅s6的转移成本张量 [nt, nT, nT]
- nT: 温度状态数
- nt: 时段数
- y_s6: 两次s6的最小间隔 (h)
- dt: 时间步长 (h)

返回:
- min_cost: 最小总成本
- min_TsIndexList: 最优温度路径 [nt+1]
- min_WaitIndexList: 最优等待时间路径 [nt+1]
- best_start_idx: 最优起始温度索引
- best_wait_idx: 最优起始等待时间索引
"""
function GraphLayerSolver(
    Cn::Array{Float64, 3}, Cs::Array{Float64, 3},
    nT::Int, nt::Int, y_s6::Float64, dt::Float64
)
    # 计算s6_wait状态数
    # wait状态: 0, dt, 2dt, ..., y_s6-dt, y_s6⁺
    ns6 = ceil(Int, y_s6 / dt) + 1  # +1 是 y⁺ 状态
    max_wait_idx = ns6 - 1           # y⁺ 状态的索引

    # 全局最优解
    global_min_cost = Inf
    global_best_start_idx = 0
    global_best_wait_idx = 0
    global_best_TsIndexList = zeros(Int, nt + 1)
    global_best_WaitIndexList = zeros(Int, nt + 1)

    # 枚举所有 (起始温度, 起始等待状态) 组合
    # 对每组单独跑DP，保证周期性: Ts[1]==Ts[nt+1] 且 wait[1]==wait[nt+1]
    for w0 in 0:max_wait_idx
        for start_idx in 1:nT
            # DP前向递推矩阵
            # V[t, w+1, j] 存储从 (start_idx, w0) 出发到达 (j, w) 的最小成本
            V = fill(Inf, nt + 1, ns6, nT)
            # 回溯矩阵: prev[t, w+1, j] = (prev_w+1, prev_i)
            prev = zeros(Int, nt + 1, ns6, nT, 2)

            # 初始化: 从 (start_idx, w0) 出发，第一步转移
            for j in 1:nT
                # 非s6转移: (start_idx, w0) → (j, next_w)
                if Cn[1, start_idx, j] < Inf
                    next_w = min(w0 + 1, max_wait_idx)
                    V[2, next_w + 1, j] = Cn[1, start_idx, j]
                    prev[2, next_w + 1, j, 1] = w0 + 1
                    prev[2, next_w + 1, j, 2] = start_idx
                end
                # s6转移: (start_idx, w0) → (j, 0)，仅当 w0 == max_wait_idx 时可用
                if w0 == max_wait_idx && Cs[1, start_idx, j] < Inf
                    if Cs[1, start_idx, j] < V[2, 1, j]
                        V[2, 1, j] = Cs[1, start_idx, j]
                        prev[2, 1, j, 1] = w0 + 1
                        prev[2, 1, j, 2] = start_idx
                    end
                end
            end

            # 正向递推 t=2..nt
            for t in 2:nt
                for w in 0:max_wait_idx
                    for i in 1:nT
                        if V[t, w + 1, i] >= Inf
                            continue
                        end
                        for j in 1:nT
                            # 非s6转移: wait递增，但不超过max_wait_idx
                            next_w = min(w + 1, max_wait_idx)
                            next_cost = V[t, w + 1, i] + Cn[t, i, j]
                            if Cn[t, i, j] < Inf && next_cost < V[t + 1, next_w + 1, j]
                                V[t + 1, next_w + 1, j] = next_cost
                                prev[t + 1, next_w + 1, j, 1] = w + 1
                                prev[t + 1, next_w + 1, j, 2] = i
                            end

                            # s6转移: 仅当 w == max_wait_idx 时可用
                            if w == max_wait_idx && Cs[t, i, j] < Inf
                                next_cost_s6 = V[t, w + 1, i] + Cs[t, i, j]
                                if next_cost_s6 < V[t + 1, 1, j]
                                    V[t + 1, 1, j] = next_cost_s6
                                    prev[t + 1, 1, j, 1] = w + 1
                                    prev[t + 1, 1, j, 2] = i
                                end
                            end
                        end
                    end
                end
            end

            # 周期边界检查: 必须回到 (start_idx, w0)
            cost = V[nt + 1, w0 + 1, start_idx]
            if cost < global_min_cost
                global_min_cost = cost
                global_best_start_idx = start_idx
                global_best_wait_idx = w0

                # 回溯路径
                TsIndexList = zeros(Int, nt + 1)
                WaitIndexList = zeros(Int, nt + 1)

                TsIndexList[nt + 1] = start_idx
                WaitIndexList[nt + 1] = w0

                cur_i = start_idx
                cur_w = w0

                for t in nt+1:-1:2
                    prev_w = prev[t, cur_w + 1, cur_i, 1]
                    prev_i = prev[t, cur_w + 1, cur_i, 2]
                    TsIndexList[t - 1] = prev_i
                    WaitIndexList[t - 1] = prev_w - 1
                    cur_i = prev_i
                    cur_w = prev_w - 1
                end

                global_best_TsIndexList = TsIndexList
                global_best_WaitIndexList = WaitIndexList
            end
        end
    end

    if global_best_start_idx == 0
        # 不可行，返回Inf
        return Inf, zeros(Int, nt + 1), zeros(Int, nt + 1), 0, 0
    end

    return global_min_cost, global_best_TsIndexList, global_best_WaitIndexList,
           global_best_start_idx, global_best_wait_idx
end

end # module DPSolverCore
