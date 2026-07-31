#=
动态规划精确解模块

本模块实现了基于连续COP函数的动态规划精确求解。
与 DP_InitialSolution.jl 中的初始解生成不同，本模块：
1. 使用连续COP函数（而非分段COP），COP取首尾平均值
2. 以MILP启发式解为起点，通过局部状态空间搜索改善解
3. 不用于MILP热启动，而是作为独立的精确求解方法

三阶段架构:
  阶段1: DP初始解 (DP_InitialSolution.jl) → 分段COP，粗略解
  阶段2: MILP启发式求解 → 用阶段1的初始解热启动，获得考虑状态切换的解
  阶段3: DP精确解 (本文件) → 连续COP，局部搜索改善

参考:
- VaryLoadVaryArea.jl 中 generateAndSolve 函数的局部搜索逻辑
- DPSolverCore.jl 中 ExhaustiveSolver 函数
- DP_InitialSolution.jl 中 single_step 函数的状态逻辑
=#

# 本文件被 include 到 HeatPumpWithStorageSystem 模块内
# MILPModelParameters, SingleStepVariables, single_step, accurate_s7_model, accurate_s8_model
# 等都定义在同一模块内，可直接访问
# DPSolverCore 是独立子模块，ExhaustiveSolver 通过 DC.ExhaustiveSolver 调用

"""
DP精确解参数结构体

字段说明:
- dt: 时间步长 (h)
- nT_local_init: 初始局部状态数（默认3，即当前温度±1个步长）
- nT_local_max: 最大局部状态数（默认5，精度足够后扩展）
- dT_initial: 初始温度步长 (℃)，局部搜索的初始搜索半径
- dT_min: 最小温度步长 (℃)，收敛判据
- max_iter: 最大迭代次数
- y_s6: 两次状态6（放热→切换直供）的最小间隔 (h)，0表示无约束
  引入此约束可防止DP频繁切换压缩机模式（级间冷却作弊）
"""
struct DP_PRECISE_PARAMS
	dt::Float64
	nT_local_init::Int
	nT_local_max::Int
	dT_initial::Float64
	dT_min::Float64
	max_iter::Int
	y_s6::Float64
end

"""
连续COP变量结构体

与 SingleStepVariables 的区别:
- COP1v, COP2v, COP3v, COPwv 从 Float64 改为 Function
- 删除 COP_lambda（用 COP1_func 代替）
- 删除 initialize（精确解始终为 false）
- 删除 price（电价在外部处理）

COP函数接口: func(Ts::Float64) -> Float64
在固定时间层上，环境温度Tair已固定，COP只与蓄热温度Ts有关。
首尾平均COP在 single_step_continous 内部计算: (func(Ts_start) + func(Ts_end)) / 2
"""
struct SingleStepVariablesContinuous
	heat_load::Float64           # 当前时段热负荷 (kW)
	COP_ca::Float64              # 热泵直接供热模式 COP（首尾平均后定值）
	COP1_func::Function          # COP1(Ts) -> Float64，模式1/切换点COP
	COP2_func::Function          # COP2(Ts) -> Float64，蓄热供热COP
	COP3_func::Function          # COP3(Ts) -> Float64，热泵储热COP
	COPw_func::Function          # COPw(Ts) -> Float64，水蒸气压缩机COP
	dt::Float64                  # 时间步长 (h)
end

"""
DP精确解结果结构体
"""
struct DPSolution
	Ts::Vector{Float64}          # 温度轨线 [nt+1]，首尾一致
	states::Vector{Int}           # 状态序列 [nt]
	P1::Vector{Float64}          # 模式1总功率 [nt]
	P2::Vector{Float64}          # 模式2总功率 [nt]
	P3::Vector{Float64}          # 模式3总功率 [nt]
	Pe_l::Vector{Float64}        # 电锅炉供热功率 [nt]
	Pe_s::Vector{Float64}        # 电锅炉储热功率 [nt]
	cost::Float64                # 总经济成本
	converged::Bool              # 是否收敛（dT_local <= dT_min）
end


"""
使用连续COP函数计算单步状态转移成本（精确解版本）

与 single_step 的区别:
1. COP值通过调用COP函数并取首尾平均得到，而非直接使用分段定值
2. s7/s8中COP_lambda取前半程平均: (COP1_func(Ts_start) + COP1_func(T_lambda)) / 2
3. 不包含 initialize 分支（精确解始终使用完整逻辑）
4. 复用 accurate_s7_model 和 accurate_s8_model（它们接收Float64）

参数:
- Ts_start: 起始蓄热温度 (℃)
- Ts_end: 结束蓄热温度 (℃)
- variables: SingleStepVariablesContinuous 结构体
- params: MILPModelParameters 结构体

返回: 与 single_step 相同的13个值

禁止状态参数:
- forbidden_states: 禁止使用的状态集合，用于分别计算Cn（不含s6）和Cs（仅s6）
  当集合为空时（默认），所有8种状态均可选
"""
function single_step_continous(
	Ts_start::Float64,
	Ts_end::Float64,
	variables::SingleStepVariablesContinuous,
	params::MILPModelParameters;
	forbidden_states::Set{Int} = Set{Int}(),
)
	heat_load = variables.heat_load
	COP_ca = variables.COP_ca
	COP1_func = variables.COP1_func
	COP2_func = variables.COP2_func
	COP3_func = variables.COP3_func
	COPw_func = variables.COPw_func
	dt = variables.dt

	C_heatpump = params.C_heatpump
	C_boiler = params.C_boiler
	Tcmax = params.Tcmax
	dTs = params.dTs
	Tuse = params.Tuse
	cpm = params.cpm

	# 计算首尾平均COP（一般状态使用）
	COP1v = (COP1_func(Ts_start) + COP1_func(Ts_end)) / 2
	COP2v = (COP2_func(Ts_start) + COP2_func(Ts_end)) / 2
	COP3v = (COP3_func(Ts_start) + COP3_func(Ts_end)) / 2
	COPwv = (COPw_func(Ts_start) + COPw_func(Ts_end)) / 2

	best_power = 9999.0
	best_state = 0
	P1s_best = 9999.0
	P2s_best = 9999.0
	P3s_best = 9999.0
	Pesl_best = 9999.0
	Pess_best = 9999.0
	P1e_best = 9999.0
	P2e_best = 9999.0
	P3e_best = 9999.0
	Peel_best = 9999.0
	Pees_best = 9999.0
	lambda_best = 0.0

	current_power = 9999.0
	P1s = 9999.0
	P2s = 9999.0
	P3s = 9999.0
	Pesl = 9999.0
	Pess = 9999.0
	P1e = 9999.0
	P2e = 9999.0
	P3e = 9999.0
	Peel = 9999.0
	Pees = 9999.0
	lambda = 0.0

	# 计算中间变量
	storage_load = cpm * (Ts_end - Ts_start) / dt

	# s1 热泵供热模式
	if storage_load >= 0 && !(1 in forbidden_states)
		P1s = min(heat_load, C_heatpump) / COP_ca
		Pesl = heat_load - P1s * COP_ca
		Pess = storage_load
		current_power = P1s + Pesl + Pess
		if current_power < best_power && C_boiler >= Pess + Pesl >= 0
			best_power = current_power
			best_state = 1
			P1s_best = P1s
			P2s_best = 0.0
			P3s_best = 0.0
			Pesl_best = Pesl
			Pess_best = Pess
			P1e_best = 0.0
			P2e_best = 0.0
			P3e_best = 0.0
			Peel_best = 0.0
			Pees_best = 0.0
			lambda_best = 1.0
		end
	end

	# s2 蓄热供热模式
	if storage_load < 0 && !(2 in forbidden_states)
		P2s = -storage_load / (COP2v - 1)
		Pesl = heat_load - P2s * COP2v
		current_power = P2s + Pesl
		if current_power < best_power && C_boiler >= Pesl >= 0
			best_power = current_power
			best_state = 2
			P1s_best = 0.0
			P2s_best = P2s
			P3s_best = 0.0
			Pesl_best = Pesl
			Pess_best = 0.0
			P1e_best = 0.0
			P2e_best = 0.0
			P3e_best = 0.0
			Peel_best = 0.0
			Pees_best = 0.0
			lambda_best = 1.0
		end
	end

	# s3 热泵储热模式
	if storage_load > 0 && Ts_end <= Tcmax - dTs && !(3 in forbidden_states)
		P3s = min(C_heatpump, storage_load) / COP3v
		Pess = storage_load - P3s * COP3v
		Pesl = heat_load
		lambda = 1.0
		current_power = P3s + Pesl + Pess
		if current_power < best_power && C_boiler >= Pess + Pesl >= 0
			best_power = current_power
			best_state = 3
			P1s_best = 0.0
			P2s_best = 0.0
			P3s_best = P3s
			Pesl_best = Pesl
			Pess_best = Pess
			P1e_best = 0.0
			P2e_best = 0.0
			P3e_best = 0.0
			Peel_best = 0.0
			Pees_best = 0.0
			lambda_best = lambda
		end
	end

	# s4: 计算合并到s6

	# s5 热泵同时供热和储热
	if storage_load > 0 && (heat_load < C_heatpump || storage_load < C_heatpump) && Ts_end <= Tcmax - dTs && !(5 in forbidden_states)
		P3s = storage_load / COP1v
		if P3s * COP1v < C_heatpump
			P1s = min(heat_load, C_heatpump - P3s * COP1v) / COP1v
			Pesl = heat_load - P1s * COP1v
			lambda = 1.0
			current_power = P1s + P3s + Pesl
			if current_power < best_power && C_boiler >= Pesl >= 0
				best_power = current_power
				best_state = 5
				P1s_best = P1s
				P2s_best = 0.0
				P3s_best = P3s
				Pesl_best = Pesl
				Pess_best = 0.0
				P1e_best = 0.0
				P2e_best = 0.0
				P3e_best = 0.0
				Peel_best = 0.0
				Pees_best = 0.0
				lambda_best = lambda
			end
		end
	end

	# s6 蓄热用完后切换为热泵供热
	if storage_load < 0
		P2s = -storage_load / (COP2v - 1)
		flag = P2s * COP2v <= C_heatpump
		P1e = min(heat_load - P2s * COP2v, C_heatpump) / COP_ca
		lambda = P2s * COP2v / heat_load
		Pesl = lambda * heat_load - P2s * COP2v
		Peel = (1 - lambda) * heat_load - P1e * COP_ca
		current_power = P1e + P2s + Pesl + Peel
		current_state = nothing
		if Ts_end >= Tuse + dTs && !(4 in forbidden_states)
			current_state = 4
		elseif !(6 in forbidden_states)
			current_state = 6
		end


		if (0 < lambda < 1 &&
			current_power < best_power &&
			C_boiler >= Pesl + Peel >= 0 && P1e >= 0 && flag &&
			!isnothing(current_state))
			best_power = current_power
			if current_state == 4
				# 模式4同时供热
				best_state = 4
				P1s_best = P1e
				P2s_best = P2s
				P3s_best = 0.0
				Pesl_best = Pesl + Peel
				Pess_best = 0.0
				P1e_best = 0.0
				P2e_best = 0.0
				P3e_best = 0.0
				Peel_best = 0.0
				Pees_best = 0.0
				lambda_best = 1.0
			elseif current_state == 6
				best_state = 6
				P1s_best = 0.0
				P2s_best = P2s
				P3s_best = 0.0
				Pesl_best = Pesl
				Pess_best = 0.0
				P1e_best = P1e
				P2e_best = 0.0
				P3e_best = 0.0
				Peel_best = Peel
				Pees_best = 0.0
				lambda_best = lambda
			end
		end
	end
	
	# s7 储热完成后热泵向工厂供热
	# 前半段热泵储热，后半段热泵供热
	# COP_lambda取前半程平均: (COP1_func(Ts_start) + COP1_func(T_lambda)) / 2
	if storage_load > 0 && (storage_load < C_heatpump || heat_load < C_heatpump) && !(7 in forbidden_states)
		obj = 9999.0
		P3s_lambda = Pesl_lambda = P1e_lambda = Pees_lambda = Peel_lambda = 9999.0
		Pess_lambda = 0.0
		Ts_lambda = Ts_start
		lambda = 0.0

		# 遍历不同的Ts_lambda，找到最优解
		Ts_lambda_list = range(Ts_start, min(Ts_end, Tcmax - dTs), length = 10)
		for T_lambda_temp in Ts_lambda_list
			# 前半程平均COP: 从Ts_start到T_lambda
			COP_lambda_v = (COP1_func(Ts_start) + COP1_func(T_lambda_temp)) / 2
			obj_temp, P3s_lambda_temp, Pesl_lambda_temp, Pess_lambda_temp, P1e_lambda_temp, Pees_lambda_temp, Peel_lambda_temp, lambda_temp = accurate_s7_model(
				COP_lambda_v,
				COP_ca,
				cpm,
				dt,
				Ts_start,
				T_lambda_temp,
				Ts_end,
				heat_load,
				C_heatpump,
				C_boiler,
			)
			if obj_temp < obj
				obj = obj_temp
				P3s_lambda = P3s_lambda_temp
				Pesl_lambda = Pesl_lambda_temp
				Pess_lambda = Pess_lambda_temp
				P1e_lambda = P1e_lambda_temp
				Pees_lambda = Pees_lambda_temp
				Peel_lambda = Peel_lambda_temp
				Ts_lambda = T_lambda_temp
				lambda = lambda_temp
			end
		end

		current_power = P3s_lambda + Pess_lambda + Pesl_lambda + P1e_lambda + Pees_lambda + Peel_lambda
		if current_power < best_power
			best_power = current_power
			best_state = 7
			P1s_best = 0.0
			P2s_best = 0.0
			P3s_best = P3s_lambda
			Pesl_best = Pesl_lambda
			Pess_best = Pess_lambda
			P1e_best = P1e_lambda
			P2e_best = 0.0
			P3e_best = 0.0
			Peel_best = Peel_lambda
			Pees_best = Pees_lambda
			lambda_best = lambda
		end
	end

	# s8 热泵同时供热和储热，蓄热完成后热泵供热
	# 前半段同时供热储热，后半段仅供热
	# COP_lambda取前半程平均: (COP1_func(Ts_start) + COP1_func(T_lambda)) / 2
	if storage_load > 0 && C_heatpump > heat_load && !(8 in forbidden_states)
		obj = 9999.0
		P3s_lambda = P1s_lambda = Pesl_lambda = P1e_lambda = Pees_lambda = Peel_lambda = 9999.0
		Pess_lambda = 0.0   # s8前半段电锅炉不储热，储热由热泵完成
		Ts_lambda = Ts_start
		lambda = 0.0

		# 遍历不同的Ts_lambda，找到最优解
		Ts_lambda_list = range(Ts_start, min(Ts_end, Tcmax - dTs), length = 10)
		for T_lambda_temp in Ts_lambda_list
			# 前半程平均COP: 从Ts_start到T_lambda
			COP_lambda_v = (COP1_func(Ts_start) + COP1_func(T_lambda_temp)) / 2
			obj_temp, P3s_lambda_temp, P1s_lambda_temp, Pesl_lambda_temp, Pess_lambda_temp, P1e_lambda_temp, Pees_lambda_temp, Peel_lambda_temp, lambda_temp = accurate_s8_model(
				COP_lambda_v,
				COP_ca,
				cpm,
				dt,
				Ts_start,
				T_lambda_temp,
				Ts_end,
				heat_load,
				C_heatpump,
				C_boiler,
			)
			if obj_temp < obj
				obj = obj_temp
				P3s_lambda = P3s_lambda_temp
				P1s_lambda = P1s_lambda_temp
				Pesl_lambda = Pesl_lambda_temp
				Pess_lambda = Pess_lambda_temp
				P1e_lambda = P1e_lambda_temp
				Pees_lambda = Pees_lambda_temp
				Peel_lambda = Peel_lambda_temp
				Ts_lambda = T_lambda_temp
				lambda = lambda_temp
			end
		end

		current_power = P3s_lambda + P1s_lambda + Pess_lambda + Pesl_lambda + P1e_lambda + Pees_lambda + Peel_lambda
		if current_power < best_power
			best_power = current_power
			best_state = 8
			P1s_best = P1s_lambda
			P2s_best = 0.0
			P3s_best = P3s_lambda
			Pesl_best = Pesl_lambda
			Pess_best = Pess_lambda
			P1e_best = P1e_lambda
			P2e_best = 0.0
			P3e_best = 0.0
			Peel_best = Peel_lambda
			Pees_best = Pees_lambda
			lambda_best = lambda
		end
	end

	return best_state, best_power, P1s_best, P2s_best, P3s_best, Pesl_best, Pess_best, P1e_best, P2e_best, P3e_best, Peel_best, Pees_best, lambda_best
end


"""
使用连续COP函数计算温度状态转移成本（精确解版本）

该函数是 single_step_continous 的封装层，负责:
1. 从 designParameters 获取连续COP函数
2. 创建 COP 函数（固定时间层上，Tair已固定，COP只与Ts有关）
3. 计算首尾平均后的 COP_ca
4. 调用 single_step_continous 计算最优功率分配

参数:
- Ts_start: 起始蓄热温度 (℃)
- Ts_end: 结束蓄热温度 (℃)
- designParameters: 设计参数结构体（包含 COPWater, COPOverlapFunction 等连续COP函数）
- params: MILPModelParameters 结构体
- sysVariables: 系统变量（含 load, price）
- Tair_start: 起始时刻环境温度 (℃)
- Tair_end: 结束时刻环境温度 (℃)
- t_index: 当前时段索引
- dt: 时间步长 (h)

返回: 与 calculateTransitionCost 相同的13个值（C为纯功率，不乘电价）

COP函数设计:
- COP1_func(Ts): 模式1/切换点COP，使用平均Tair
- COP2_func(Ts): 蓄热供热COP，不依赖Tair
- COP3_func(Ts): 热泵储热COP，使用平均Tair
- COPw_func(Ts): 水蒸气压缩机COP，不依赖Tair

注意: 在固定时间层上，Tair取首尾平均 (Tair_start + Tair_end) / 2，
COP函数内部用此固定Tair，首尾平均在 single_step_continous 中通过
(func(Ts_start) + func(Ts_end)) / 2 实现。
"""
function calculateTransitionCost_continous(
	Ts_start::Float64, Ts_end::Float64,
	designParameters::DesignOptimizeParameters, params::MILPModelParameters,
	sysVariables::SysVariables, Tair_start::Float64, Tair_end::Float64,
	t_index::Int, dt::Float64;
	forbidden_states::Set{Int} = Set{Int}(),
)
	# 温度范围检查：超出可行范围直接返回大数
	# 蓄热温度必须在 [Ts_min, Ts_max] 范围内
	if Ts_start < params.Tsmin || Ts_start > params.Tsmax ||
	   Ts_end < params.Tsmin || Ts_end > params.Tsmax
		return 9999.0, 0, 9999.0, 9999.0, 9999.0, 9999.0, 9999.0,
		9999.0, 9999.0, 9999.0, 9999.0, 9999.0, 0.0
	end

	# 从designParameters获取常量参数
	dTES = designParameters.dT_EvaporationStandard
	Tuse = designParameters.Tuse
	ThMax = designParameters.ThMax
	TCompressorIn = designParameters.TCompressorIn

	# 平均环境温度（固定时间层上使用）
	Tair_avg = (Tair_start + Tair_end) / 2
	# 蒸发温度 Te = Tair - dT_EvaporationStandard
	Te_avg = Tair_avg - dTES
	Te_start = Tair_start - dTES
	Te_end = Tair_end - dTES

	# COP1函数: 模式1/切换点COP
	# COP1(Te, Ts) = COPOverlapFunction(Te, max(Tuse, Ts + dT_EvaporationStandard))
	# 当 Ts + dTES > ThMax 时返回1.0（不可行）
	COP1_func = (Ts::Float64) -> begin
		Ts_cond = max(Tuse, Ts + dTES)
		if Ts_cond > ThMax
			return 1.0
		end
		return designParameters.COPOverlapFunction(Te_avg, Ts_cond)
	end

	# COP2函数: 蓄热供热COP（不依赖Tair）
	# COP2(Ts) = COPWater(Ts - dT_EvaporationStandard, Tuse)
	COP2_func = (Ts::Float64) -> begin
		return designParameters.COPWater(Ts - dTES, Tuse)
	end

	# COP3函数: 热泵储热COP
	# COP3(Te, Ts) = COPOverlapFunction(Te, Ts + dT_EvaporationStandard)
	# 当 Ts + dTES >= ThMax 时返回1.0（不可行）
	COP3_func = (Ts::Float64) -> begin
		Ts_cond = Ts + dTES
		if Ts_cond >= ThMax
			return 1.0
		end
		return designParameters.COPOverlapFunction(Te_avg, Ts_cond)
	end

	# COPw函数: 水蒸气压缩机COP（不依赖Tair）
	# COPw(Ts) = COPWater(TCompressorIn, Ts)
	# 当 Ts > ThMax 时返回1.0（不可行）
	COPw_func = (Ts::Float64) -> begin
		if Ts > ThMax
			return 1.0
		end
		return designParameters.COPWater(TCompressorIn, Ts)
	end

	# COPca: 首尾平均后定值（不依赖Ts）
	# COPca = COPOverlapFunction(Te, Tuse)，其中 Te = Tair - dTES
	COPca = (
		designParameters.COPOverlapFunction(Te_start, Tuse) +
		designParameters.COPOverlapFunction(Te_end, Tuse)
	) / 2

	heat_load = sysVariables.load[t_index]

	variables = SingleStepVariablesContinuous(
		heat_load, COPca,
		COP1_func, COP2_func, COP3_func, COPw_func,
		dt,
	)

	state, C, P1s, P2s, P3s, Pesl, Pess, P1e, P2e, P3e, Peel, Pees, lambda =
		single_step_continous(Ts_start, Ts_end, variables, params; forbidden_states = forbidden_states)

	return C, state, P1s, P2s, P3s, Pesl, Pess, P1e, P2e, P3e, Peel, Pees, lambda
end


"""
为局部状态空间构建成本张量（精确解版本）

与 buildCostTensor 的区别:
1. 每个时刻的状态空间可以不同（TsMatrix[:, t] ≠ TsMatrix[:, t+1]）
2. 使用连续COP函数（通过 calculateTransitionCost_continous）
3. 不使用工况缓存（温度点每次迭代都变）
4. 存储完整的功率矩阵，方便后续从DP结果提取功率

参数:
- TsMatrix: 温度状态矩阵 [nT, nt+1]，每列是一个时刻的状态空间
- nt: 时段数
- dt: 时间步长
- designParameters: 设计参数（连续COP函数来源）
- params: MILPModelParameters 结构体
- sysVariables: 系统变量
- Tair_list: 各时刻环境温度 [nt+1]

返回:
  当 y_s6 > 0 时: (Cn, Cs, stateMatrix, P1sMatrix, ...)
	Cn: 不含s6的成本张量，Cs: 仅s6的成本张量
	同时返回对应的状态和功率矩阵（Cn和Cs共用同一套矩阵）
  当 y_s6 = 0 时: (C, stateMatrix, P1sMatrix, ...)
	C: 全状态成本张量
"""
function buildCostTensor_Precise(
	TsMatrix::Matrix{Float64},
	nt::Int, dt::Float64,
	designParameters::DesignOptimizeParameters, params::MILPModelParameters,
	sysVariables::SysVariables, Tair_list::Vector{Float64};
	y_s6::Float64 = 0.0,
)
	nT = size(TsMatrix, 1)
	use_s6_constraint = y_s6 > 0.0

	# 分配成本张量和功率矩阵（统一返回类型，始终分配所有矩阵）
	Cn = fill(Inf, nt, nT, nT)
	Cs = fill(Inf, nt, nT, nT)  # 始终分配，未启用时保持 Inf
	stateMatrix = zeros(Int, nt, nT, nT)
	P1sMatrix = zeros(nt, nT, nT)
	P2sMatrix = zeros(nt, nT, nT)
	P3sMatrix = zeros(nt, nT, nT)
	P1eMatrix = zeros(nt, nT, nT)
	P2eMatrix = zeros(nt, nT, nT)
	P3eMatrix = zeros(nt, nT, nT)
	PeslMatrix = zeros(nt, nT, nT)
	PessMatrix = zeros(nt, nT, nT)
	PeelMatrix = zeros(nt, nT, nT)
	PeesMatrix = zeros(nt, nT, nT)
	lambdaMatrix = zeros(nt, nT, nT)
	# s6专用矩阵（状态6没有P1s, P3s, Pess, Pees, P2e, P3e）
	P2s_s6_Matrix = zeros(nt, nT, nT)
	Pesl_s6_Matrix = zeros(nt, nT, nT)
	P1e_s6_Matrix = zeros(nt, nT, nT)
	Peel_s6_Matrix = zeros(nt, nT, nT)
	lambda_s6_Matrix = zeros(nt, nT, nT)

	for t in 1:nt
		Tair_start = Tair_list[t]
		Tair_end = Tair_list[t+1]

		for i in 1:nT
			Ts_start = TsMatrix[i, t]
			for j in 1:nT
				Ts_end = TsMatrix[j, t+1]

				# 成本计算：根据是否启用s6约束决定forbidden_states
				# - s6约束启用: Cn禁止s6（分离到Cs）
				# - s6约束禁用: 允许所有状态
				forbidden = use_s6_constraint ? Set{Int}([6]) : Set{Int}()
				power_cost, state, P1s, P2s, P3s, Pesl, Pess,
				P1e, P2e, P3e, Peel, Pees, lambda = calculateTransitionCost_continous(
					Ts_start, Ts_end, designParameters, params,
					sysVariables, Tair_start, Tair_end, t, dt;
					forbidden_states = forbidden,
				)

				if power_cost <= 9999.0
					Cn[t, i, j] = power_cost
					stateMatrix[t, i, j] = state
					P1sMatrix[t, i, j] = P1s
					P2sMatrix[t, i, j] = P2s
					P3sMatrix[t, i, j] = P3s
					P1eMatrix[t, i, j] = P1e
					P2eMatrix[t, i, j] = P2e
					P3eMatrix[t, i, j] = P3e
					PeslMatrix[t, i, j] = Pesl
					PessMatrix[t, i, j] = Pess
					PeelMatrix[t, i, j] = Peel
					PeesMatrix[t, i, j] = Pees
					lambdaMatrix[t, i, j] = lambda
				end

				# 若启用s6约束，额外计算仅s6的转移成本 Cs
				if use_s6_constraint
					power_cost_s6, _, _, P2s_s6, _, Pesl_s6, _,
					P1e_s6, _, _, Peel_s6, _, lambda_s6 = calculateTransitionCost_continous(
						Ts_start, Ts_end, designParameters, params,
						sysVariables, Tair_start, Tair_end, t, dt;
						forbidden_states = Set{Int}([1, 2, 3, 4, 5, 7, 8]),
					)
					if power_cost_s6 <= 9999.0
						Cs[t, i, j] = power_cost_s6
						P2s_s6_Matrix[t, i, j] = P2s_s6
						Pesl_s6_Matrix[t, i, j] = Pesl_s6
						P1e_s6_Matrix[t, i, j] = P1e_s6
						Peel_s6_Matrix[t, i, j] = Peel_s6
						lambda_s6_Matrix[t, i, j] = lambda_s6
					end
				end
			end
		end
	end

    # 功率平滑
    Cn += sysVariables.epsilon*((P1sMatrix + P1eMatrix).^2+(P2sMatrix + P2eMatrix).^2+(P3sMatrix + P3eMatrix).^2+(PeslMatrix + PessMatrix + PeelMatrix + PeesMatrix).^2)
    Cs += sysVariables.epsilon*((P2s_s6_Matrix + P1e_s6_Matrix).^2+(Pesl_s6_Matrix + Peel_s6_Matrix).^2)

	return Cn, Cs, stateMatrix,
		   P1sMatrix, P2sMatrix, P3sMatrix,
		   P1eMatrix, P2eMatrix, P3eMatrix,
		   PeslMatrix, PessMatrix, PeelMatrix, PeesMatrix,
		   lambdaMatrix, P2s_s6_Matrix, Pesl_s6_Matrix, P1e_s6_Matrix, Peel_s6_Matrix, lambda_s6_Matrix
end


"""
使用局部状态空间搜索求解DP精确解

算法流程（参考 VaryLoadVaryArea.jl generateAndSolve）:
1. 以MILP解为起点，生成局部状态空间
2. 构建成本张量，求解DP（保证首尾温度一致）
3. 根据结果自适应调整:
   - 解在边界: 保持步长，重新生成范围
   - 解在内部: 减小步长
   - 连续多次边界: 重置步长
4. 精度足够后增加状态数（3→5）
5. 迭代直到收敛或达最大次数
6. 求解完成后，沿最优路径重新计算各时段的功率和状态

参数:
- TsList_milp: MILP解的温度轨线 [nt+1]
- dp_precise_params: DP精确解参数
- designParameters: 设计参数（连续COP函数来源）
- params: MILPModelParameters 结构体
- sysVariables: 系统变量
- Tair_list: 各时刻环境温度 [nt+1]

返回: DPSolution 结构体
"""
function solvePreciseDP(
	TsList_milp::Vector{Float64},
	dp_precise_params::DP_PRECISE_PARAMS,
	designParameters::DesignOptimizeParameters,
	params::MILPModelParameters,
	sysVariables::SysVariables,
	Tair_list::Vector{Float64},
)
	nt = length(TsList_milp) - 1
	dt = dp_precise_params.dt

	nT = dp_precise_params.nT_local_init
	half_nT = Int((nT - 1) / 2)
	dT_local = dp_precise_params.dT_initial

	use_s6 = dp_precise_params.y_s6 > 0.0

	@info "DP精确解求解开始"
	@info "时段数: $nt, 初始状态数: $nT, 初始步长: $dT_local, s6约束: $use_s6"

	# 初始化局部状态空间: 以MILP解为中心
	TsMatrix = zeros(nT, nt + 1)
	for j in 1:nt+1
		TsMatrix[:, j] = TsList_milp[j] .+ (-half_nT:half_nT) .* dT_local
	end

	TsList = copy(TsList_milp)
	countAll = 0
	countSingleGap = 0
	is_nt_changed = false

	# 记录最优路径信息（存储完整矩阵，便于直接提取功率）
	best_cost = Inf
	best_TsList = copy(TsList_milp)
	best_TsIndexList = fill(half_nT + 1, nt + 1)
	best_WaitIndexList = use_s6 ? zeros(Int, nt + 1) : Int[]
	# 最优功率矩阵（在最优迭代时保存）
	best_stateMatrix = zeros(Int, nt, nT, nT)
	best_P1sMatrix = zeros(nt, nT, nT)
	best_P2sMatrix = zeros(nt, nT, nT)
	best_P3sMatrix = zeros(nt, nT, nT)
	best_P1eMatrix = zeros(nt, nT, nT)
	best_P2eMatrix = zeros(nt, nT, nT)
	best_P3eMatrix = zeros(nt, nT, nT)
	best_PeslMatrix = zeros(nt, nT, nT)
	best_PessMatrix = zeros(nt, nT, nT)
	best_PeelMatrix = zeros(nt, nT, nT)
	best_PeesMatrix = zeros(nt, nT, nT)
	best_lambdaMatrix = zeros(nt, nT, nT)
	# s6专用最优矩阵
	best_P2s_s6_Matrix = zeros(nt, nT, nT)
	best_Pesl_s6_Matrix = zeros(nt, nT, nT)
	best_P1e_s6_Matrix = zeros(nt, nT, nT)
	best_Peel_s6_Matrix = zeros(nt, nT, nT)
	best_lambda_s6_Matrix = zeros(nt, nT, nT)

	total_tensor_time = 0.0
	total_dp_time = 0.0

	while (dT_local > dp_precise_params.dT_min && countAll < dp_precise_params.max_iter)
		# 1. 构建成本张量（连续COP，存储完整功率矩阵，统一返回类型）
		tensor_time = @elapsed begin
			Cn_power, Cs_power, stateMatrix, P1sMatrix, P2sMatrix, P3sMatrix,
			P1eMatrix, P2eMatrix, P3eMatrix, PeslMatrix, PessMatrix, PeelMatrix, PeesMatrix,
			lambdaMatrix, P2s_s6_Matrix, Pesl_s6_Matrix, P1e_s6_Matrix, Peel_s6_Matrix, lambda_s6_Matrix = buildCostTensor_Precise(
				TsMatrix, nt, dt, designParameters, params, sysVariables, Tair_list;
				y_s6 = dp_precise_params.y_s6,
		)
	end
	total_tensor_time += tensor_time

		# 2. 乘以电价得经济成本
		dp_time = @elapsed begin
			if use_s6
				Cn_econ = fill(Inf, nt, nT, nT)
				Cs_econ = fill(Inf, nt, nT, nT)
				for t in 1:nt
					price_t = sysVariables.price[t]
					for i in 1:nT
						for j in 1:nT
							if Cn_power[t, i, j] < Inf
								Cn_econ[t, i, j] = Cn_power[t, i, j] * price_t * dt
							end
							if Cs_power[t, i, j] < Inf
								Cs_econ[t, i, j] = Cs_power[t, i, j] * price_t * dt
							end
						end
					end
				end

				# 图DP求解（带s6间隔约束）
				cost, TsIndexList, WaitIndexList, best_start_idx, best_wait_idx =
					DC.GraphLayerSolver(Cn_econ, Cs_econ, nT, nt, dp_precise_params.y_s6, dt)
			else
				C_economic = fill(Inf, nt, nT, nT)
				for t in 1:nt
					price_t = sysVariables.price[t]
					for i in 1:nT
						for j in 1:nT
							if Cn_power[t, i, j] < Inf
								C_economic[t, i, j] = Cn_power[t, i, j] * price_t * dt
							end
						end
					end
				end

				# 标准DP求解（无s6约束）
				cost, TsIndexList, best_start_idx = DC.ExhaustiveSolver(C_economic, nT, nt)
			end
		end
		total_dp_time += dp_time

		# 3. 提取温度路径
		TsList = [TsMatrix[TsIndexList[j], j] for j in 1:nt+1]

		@info "迭代 $countAll: cost=$(round(cost, digits=4)), dT_local=$(round(dT_local, digits=4)), nT=$nT, 张量耗时=$(round(tensor_time, digits=3))s, DP耗时=$(round(dp_time, digits=3))s"

		# 记录最优结果（存储完整矩阵）
		if cost < best_cost
			best_cost = cost
			best_TsList = copy(TsList)
			best_TsIndexList = copy(TsIndexList)
			if use_s6
				best_WaitIndexList = copy(WaitIndexList)
			end
			# 保存当前迭代的功率矩阵
			best_stateMatrix = copy(stateMatrix)
			best_P1sMatrix = copy(P1sMatrix)
			best_P2sMatrix = copy(P2sMatrix)
			best_P3sMatrix = copy(P3sMatrix)
			best_P1eMatrix = copy(P1eMatrix)
			best_P2eMatrix = copy(P2eMatrix)
			best_P3eMatrix = copy(P3eMatrix)
			best_PeslMatrix = copy(PeslMatrix)
			best_PessMatrix = copy(PessMatrix)
			best_PeelMatrix = copy(PeelMatrix)
			best_PeesMatrix = copy(PeesMatrix)
			best_lambdaMatrix = copy(lambdaMatrix)
			# 保存s6专用矩阵
			best_P2s_s6_Matrix = copy(P2s_s6_Matrix)
			best_Pesl_s6_Matrix = copy(Pesl_s6_Matrix)
			best_P1e_s6_Matrix = copy(P1e_s6_Matrix)
			best_Peel_s6_Matrix = copy(Peel_s6_Matrix)
			best_lambda_s6_Matrix = copy(lambda_s6_Matrix)
		end

		# 4. 增加状态数（3→5）
		if dT_local <= 1.0 && !is_nt_changed
			nT = dp_precise_params.nT_local_max
			half_nT = Int((nT - 1) / 2)
			is_nt_changed = true
			TsMatrix = zeros(nT, nt + 1)
			for j in 1:nt+1
				TsMatrix[:, j] = TsList[j] .+ (-half_nT:half_nT) .* dT_local
			end
			countAll += 1
			continue
		end

		# 5. 自适应调整
		flag_nextgap = true
		for j in 1:nt+1
			if TsIndexList[j] == 1 || TsIndexList[j] == nT
				TsMatrix[:, j] = TsList[j] .+ (-half_nT:half_nT) .* dT_local
				flag_nextgap = false
			end
		end

		if flag_nextgap
			dT_local = dT_local / 2
			for j in 1:nt+1
				TsMatrix[:, j] = TsList[j] .+ (-half_nT:half_nT) .* dT_local
			end
			countSingleGap = 0
		else
			countSingleGap += 1
			if countSingleGap > 6
				dT_local = min(dT_local * 1.5, dp_precise_params.dT_initial)
				for j in 1:nt+1
					TsMatrix[:, j] = TsList[j] .+ (-half_nT:half_nT) .* dT_local
				end
				countSingleGap = 0
			end
		end
		countAll += 1
	end

	@info "DP精确解求解完成"
	@info "最优成本: $(round(best_cost, digits=4)), 迭代次数: $countAll, 最终步长: $(round(dT_local, digits=4))"
	@info "总耗时: 张量构建=$(round(total_tensor_time, digits=2))s, DP求解=$(round(total_dp_time, digits=2))s"

	# 沿最优路径直接从矩阵提取功率和状态（不再重新计算）
	TsList = best_TsList
	TsIndexList = best_TsIndexList

	states = zeros(Int, nt)
	P1_list = zeros(nt)
	P2_list = zeros(nt)
	P3_list = zeros(nt)
	Pe_l_list = zeros(nt)
	Pe_s_list = zeros(nt)

	extract_time = @elapsed begin
		for t in 1:nt
			i = TsIndexList[t]      # 起始温度索引
			j = TsIndexList[t+1]  # 结束温度索引

			# 判断是否使用了s6（wait被重置为0表示使用了s6）
			s6_used = use_s6 && (best_WaitIndexList[t+1] == 0)

			if s6_used
				# s6状态：从s6专用矩阵提取（状态6没有P1s, P3s, Pess, Pees, P2e, P3e）
				states[t] = 6
				P1_list[t] = best_P1e_s6_Matrix[t, i, j]  # P1s_s6=0，只有P1e
				P2_list[t] = best_P2s_s6_Matrix[t, i, j]  # P2e_s6=0，只有P2s
				P3_list[t] = 0.0                          # P3s_s6=P3e_s6=0
				Pe_l_list[t] = best_Pesl_s6_Matrix[t, i, j] + best_Peel_s6_Matrix[t, i, j]
				Pe_s_list[t] = 0.0                        # Pess_s6=Pees_s6=0
			else
				# 普通状态：从普通矩阵提取
				states[t] = best_stateMatrix[t, i, j]
				P1_list[t] = best_P1sMatrix[t, i, j] + best_P1eMatrix[t, i, j]
				P2_list[t] = best_P2sMatrix[t, i, j] + best_P2eMatrix[t, i, j]
				P3_list[t] = best_P3sMatrix[t, i, j] + best_P3eMatrix[t, i, j]
				Pe_l_list[t] = best_PeslMatrix[t, i, j] + best_PeelMatrix[t, i, j]
				Pe_s_list[t] = best_PessMatrix[t, i, j] + best_PeesMatrix[t, i, j]
			end
		end
	end
	@info "功率提取耗时: $(round(extract_time, digits=3))s"

	# 从 best_cost 中扣除功率正则项（正则项仅用于引导DP搜索，不参与实际成本）
	reg_deduct_time = @elapsed begin
		regularization_sum = 0.0
		epsilon = sysVariables.epsilon
		for t in 1:nt
			i = best_TsIndexList[t]
			j = best_TsIndexList[t+1]
			price_t = sysVariables.price[t]

			if use_s6 && (best_WaitIndexList[t+1] == 0)
				# s6路径：使用s6矩阵的正则项
				P2s_s6 = best_P2s_s6_Matrix[t, i, j]
				P1e_s6 = best_P1e_s6_Matrix[t, i, j]
				Pesl_s6 = best_Pesl_s6_Matrix[t, i, j]
				Peel_s6 = best_Peel_s6_Matrix[t, i, j]
				reg_t = epsilon * ((P2s_s6 + P1e_s6)^2 + (Pesl_s6 + Peel_s6)^2)
			else
				# 普通路径：使用普通矩阵的正则项
				P1s = best_P1sMatrix[t, i, j]
				P1e = best_P1eMatrix[t, i, j]
				P2s = best_P2sMatrix[t, i, j]
				P2e = best_P2eMatrix[t, i, j]
				P3s = best_P3sMatrix[t, i, j]
				P3e = best_P3eMatrix[t, i, j]
				Pesl = best_PeslMatrix[t, i, j]
				Pess = best_PessMatrix[t, i, j]
				Peel = best_PeelMatrix[t, i, j]
				Pees = best_PeesMatrix[t, i, j]
				reg_t = epsilon * ((P1s + P1e)^2 + (P2s + P2e)^2 + (P3s + P3e)^2 + (Pesl + Pess + Peel + Pees)^2)
			end
			regularization_sum += reg_t * price_t * dt
		end
		best_cost -= regularization_sum
	end
	@info "正则项扣除耗时: $(round(reg_deduct_time, digits=3))s, 扣除金额: $(round(regularization_sum, digits=4))"

	converged = dT_local <= dp_precise_params.dT_min

	return DPSolution(
		TsList, states,
		P1_list, P2_list, P3_list,
		Pe_l_list, Pe_s_list,
		best_cost, converged,
	)
end
