
# 使用经济性优化需要在项目主目录下切换到calculations环境
#=
用于测试 getMinimumCost 这类函数的计算结果
=#
using Pkg

#Pkg.activate("calculations")
using Plots
using HeatPumpWithStorageSystem

using DataFrames, CSV
using CoolProp
using COPT, JuMP

#=
测试DP生成初始解
=#

situation = "situation31"
#第一步，指定设计条件变量
#项目设计条件
begin
	heatPumpServiceCoff, heatStorageCapacity, maxheatStorageInputHour = 0.4, 3.0, 0.5
	
	#=
	hourlyTariff = zeros(24)
	hourlyTariff[1:8] .= 99.0094
	hourlyTariff[9:16] .= 0.658
	hourlyTariff[17:24] .= 0.3725

	Tair = fill(85.0,25)

	heatConsumptionPower = ones(24)
	=#
	
	hourly_tariff_ori = ones(48)
	p = 1.7#1.7
	pp = p * 1.2
	v = 0.35
	vv = v * 0.8
	hourly_tariff_ori[1:12] .*= v
	hourly_tariff_ori[23:26] .*= v
	hourly_tariff_ori[29:30] .*= pp
	hourly_tariff_ori[31:39] .*= p
	hourly_tariff_ori[40:43] .*= pp
	hourly_tariff_ori[44] *= p
	hourlyTariff = hourly_tariff_ori * 0.7393
	
	heatConsumptionPower = vcat(
		fill(0.0, 16),
		fill(1.0, 8),
		fill(0.0, 2),
		fill(1.0, 8),
		fill(0.0, 14),
	)
	heatConsumptionPower = fill(1.0, 48)
	Tair = fill(85.0, 49)
	
	dt = 0.5
	if length(hourlyTariff) != length(heatConsumptionPower) || length(hourlyTariff) != length(Tair) - 1
		@warn "hourlyTariff, heatConsumptionPower, and Tair must have the same length"
	end

	# 系数
	#heatPumpServiceCoff = 0.5
	maxCOP = 21.0# 最大COP
	eta_s = 0.7# 绝热效率
	workingStartHour = 0                # 生产开始时间
	workingHours = 24                   # 每日工作小时数
	TWaste = 85.0                     # 废热源温度
	#Tair = 25.0                        # 外部环境温度
	TCompressorIn = 115.0
	maxTcHigh = 180.0
	dT_EvaporationStandard = 5.0
	Tsmin = 120.0
	Tsmax = 220.0
	dTRecycleSupply = 5.0
	dTRecycleBackward = 5.0

	# 计算参数
	dT = 0.01
	#dt = 1/2# 时间步长过小会导致初始温度优化的目标不是一个单峰函数


	smoother = 1e-8
	Tuse = PropsSI("T", "P", 0.45e6, "Q", 0, "water") - 273.15
	Tuse = 150.0
	refrigerant = R1233zdE_Water
	sysStruct = RecycleStruct(1, 0, 0)

	segmentHeatLoad = heatConsumptionPower[1:48]
	segmentTariff = hourlyTariff[1:48]
	segmentTair = Tair[1:48]
	dt_list = ones(length(segmentHeatLoad)) * dt
	
	t_list = vcat(0.0, cumsum(dt_list))

	# 生成合并后的热负荷、电价、环境温度
	n_segments = length(dt_list)

	lambda_norm = 1e-4

	T_e = segmentTair .- dT_EvaporationStandard
end

# 经济性参数条件
begin
	# 设备成本
	waterCompresorCost = 1200.0# 水蒸气压缩机单位供热功率成本 元/kW
	lowHPCost = 1450.0# 低温热泵单位供热功率成本 元/kW
	elecHeaterCost = 1000.0# 电极锅炉单位供热功率成本 元/kW
	storageCost = 450e4 / 1000# 承压水蓄热单位体积成本 元/m³


	# 安装费系数
	heatpumpInstallCoff = 2
	elecHeaterInstallCoff = 1.2
	storageInstallCoff = 1.2

	# 维保费用
	heatpumpAnnualCost = 0.1
	elecHeaterAnnualCost = 0.05
	storageAnnualCost = 0.05

	# 折现率
	discountRate = 0.032
	# 年运行天数
	annualDays = 300
	# 运行年数
	lifeYears = 6
end

# 第二步，生成设计参数输入结构体
begin
	designInput = DesignOptimizeInput(;
		hourlyTariff = hourlyTariff,
		Tair = Tair,
		heatConsumptionPower = heatConsumptionPower,
		maxCOP = maxCOP,
		eta_s = eta_s,
		workingStartHour = workingStartHour,
		workingHours = workingHours, Tuse = Tuse,
		TWaste = TWaste,
		TCompressorIn = TCompressorIn,
		maxTcHigh = maxTcHigh,
		dT_EvaporationStandard = dT_EvaporationStandard,
		Tsmin = Tsmin,
		Tsmax = Tsmax,
		dTRecycleSupply = dTRecycleSupply,
		dTRecycleBackward = dTRecycleBackward,
		dT = dT,
		dt = dt,
		smoother = smoother,
		sysStruct = sysStruct,
		refrigerant = refrigerant,
	)
	# 第三步，生成设计参数常数结构体
	designParameters = generateDesignOptimizeParameters(PressedWaterOneStorageOneCompressor(), designInput)


	COP_ca = map(T -> designParameters.COPOverlapFunction(T, Tuse), T_e)



	# 第四步，生成带优化的目标函数
	optimizeFunction, getCount, getControlResult, getParams = generateOperationFunction(PressedWaterOneStorageOneCompressor(), designParameters, designInput)

	# heatPumpServiceCoff和maxheatStorageInputHour只是用来计算params，heatStorageCapacity是预设的常数参数
	params = getParams(heatPumpServiceCoff, heatStorageCapacity, maxheatStorageInputHour)
	cpm = params.cpm# kWh/K
	n = length(dt_list)

	# COP分档信息：
	T_s_1g_list = vcat(120.0:2.0:185.0,220.0)
	T_s_2g_list = vcat(120.0:2.0:185.0,220.0)
	T_s_3g_list = vcat(120.0:2.0:185.0,220.0)
	T_s_wg_list = vcat(120.0:2.0:185.0,220.0)

	m1, m2, m3, mw, COP1v, COP2v, COP3v, COPwv, T1g, T2g, T3g, Twg = getCOP_piecewise_data(
		designParameters,
		T_s_1g_list,
		T_s_2g_list,
		T_s_3g_list,
		T_s_wg_list,
		t_list,
	)

	M_1 = Tsmax - Tsmin
	M_2 = max(maxTcHigh - dT_EvaporationStandard - Tsmin, Tsmax - maxTcHigh + dT_EvaporationStandard)
	M_3 = max(Tuse - Tsmin + dT_EvaporationStandard, Tsmax - Tuse - dT_EvaporationStandard)
	M_4 = max(
		maximum(T_s_1g_list[2:end] - T_s_1g_list[1:end-1]),
		maximum(T_s_2g_list[2:end] - T_s_2g_list[1:end-1]),
		maximum(T_s_3g_list[2:end] - T_s_3g_list[1:end-1]),
		maximum(T_s_wg_list[2:end] - T_s_wg_list[1:end-1]),
	)# 档位温差
	M_5 = heatPumpServiceCoff * maximum(heatConsumptionPower)
end
# 第五步，计算一些经济性参数
begin
	p = 1 / (1 + discountRate)
	# 折线系数
	pc = (1 - p^lifeYears) / discountRate
	finalWaterCompresorCost = waterCompresorCost * (heatpumpInstallCoff + heatpumpAnnualCost * p * (1 - p^lifeYears) / (1 - p))
	finalLowHPCost = lowHPCost * (1 - 1 / designParameters.COPWater(TCompressorIn + dT_EvaporationStandard, Tuse)) * (heatpumpInstallCoff + heatpumpAnnualCost * p * (1 - p^lifeYears) / (1 - p))

	finalElecHeaterCost = elecHeaterCost * (elecHeaterInstallCoff + elecHeaterAnnualCost * p * (1 - p^lifeYears) / (1 - p))

	# 把蓄热的立方米造价转换成蓄热时长造价
	# 1立方米的蓄热能量除以3600秒，得到1立方够用多久
	hour_per_m3 = 1 * 900 * 4.275 * (Tsmax - Tuse) / 3600 # 1立方可以用hour_per_m3小时
	storageHourCost = storageCost / hour_per_m3# 单位时长成本 元/kWh
	finalStorageCost = storageHourCost * (storageInstallCoff + storageAnnualCost * p * (1 - p^lifeYears) / (1 - p))
end

# 第六步，生成双层优化的目标函数
fp = FinanceParameters(
	Plow_cost = finalLowHPCost,       # 低温热泵价格 (¥/kW)
	Phigh_cost = finalWaterCompresorCost,      # 高温热泵价格 (¥/kW)
	Pelec_cost = finalElecHeaterCost,       # 电锅炉价格 (¥/kW)
	Storage_cost = finalStorageCost,     # 蓄热设备价格 (¥/kWh)
	Life_years = lifeYears,          # 设备寿命 (年)
	annual_days = annualDays,        # 年运行天数
	Discount_rate = discountRate,      # 折现率 
)

# 第七步，构建MILP模型参数并求解

# 计算蓄热系统容量 C_storage (kWh)
# C_storage = cpm * (Tsmax - Tuse) / 3600
C_storage_value = heatStorageCapacity

# 创建MILP模型参数
milp_params = MILPModelParameters(
	# 系统设计参数
	Tuse = Tuse,
	Tsmin = Tsmin,
	Tsmax = Tsmax,
	Tcmax = maxTcHigh,
	dTs = dTRecycleSupply,
	Te = T_e[1],
	cpm = cpm,
	dt = dt,
	COPmax = maxCOP,
	#designParameters.COPOverlapFunction(Te, max(designParameters.Tuse, Ts + designParameters.dT_EvaporationStandard))
	COPdesign = map(i->designParameters.COPOverlapFunction(designParameters.TairFunction(t_list[i]) - dT_EvaporationStandard, Tuse), 1:n),
	lambdaNorm = lambda_norm,  # 注意: 使用lambda_norm变量
	smoother = smoother,

	# 设备容量参数
	C_heatpump = heatPumpServiceCoff * maximum(heatConsumptionPower),
	C_boiler = params.PeMax,#,maximum(heatConsumptionPower)#
	C_storage = C_storage_value,
	optimizeCapacity = false,  # 固定容量优化

	# 设备成本系数
	c_heatpump = finalLowHPCost,
	c_storage = finalStorageCost,
	c_boiler = finalElecHeaterCost,

	# 经济性参数
	D = annualDays,
	r = discountRate,
	N = lifeYears,

	# 外部输入曲线
	hourlyTariff = hourlyTariff,
	heatLoad = segmentHeatLoad,  # 使用合并后的热负荷
	Tair = Tair,
	dt_list = dt_list,

	# 时间分段参数
	n_segments = length(dt_list),
	segmentTariff = segmentTariff,  # 使用合并后的电价

	# COP分档数据
	m1 = m1,
	COP1v = COP1v,
	T1g = T1g, m2 = m2,
	COP2v = COP2v,
	T2g = T2g, m3 = m3,
	COP3v = COP3v,
	T3g = T3g, mw = mw,
	COPwv = COPwv,
	Twg = Twg,

	# 定参数COP
	COPca = COP_ca,

	# 大M参数 (文档1.2.1节)
	M1 = M_1,
	M2 = M_2,
	M3 = M_3,
	M4 = M_4,
	M5 = M_5,

	# 求解器选择
	solver = :COPT,  # 使用 COPT 求解器 (SOS1约束)
)

# 收敛数据记录类型：用于记录求解过程中的下界、上界和最优解
begin
	mutable struct ConvergenceData
		upper_bounds::Vector{Float64}          # 全局最优上界历史（单调递降）
		lower_bounds::Vector{Float64}          # 全局最优下界历史（单调递增）
		upper_iterations::Vector{Int}          # 上界对应的迭代次数
		lower_iterations::Vector{Int}          # 下界对应的迭代次数
		global_upper_bound::Float64            # 当前全局最小上界
		global_lower_bound::Float64            # 当前全局最大下界
		best_obj::Ref{Float64}                 # 当前最优目标值
		best_C_heatpump::Ref{Float64}          # 当前最优热泵容量
		best_C_boiler::Ref{Float64}            # 当前最优电加热容量
		best_Ts::Ref{Vector{Float64}}          # 当前最优温度曲线
	end

	# 初始化收敛数据
	convergence_data = ConvergenceData(
		Float64[],                             # upper_bounds
		Float64[],                             # lower_bounds
		Int[],                                 # upper_iterations
		Int[],                                 # lower_iterations
		Inf,                                   # global_upper_bound
		-Inf,                                  # global_lower_bound
		Ref(Inf),                              # best_obj
		Ref(0.0),                              # best_C_heatpump
		Ref(0.0),                              # best_C_boiler
		Ref(Float64[])                         # best_Ts
	)

	# 迭代计数器
	iteration_counter = Ref(0)

	# 定义回调函数：当找到可行解时触发
	function incumbent_callback(cb_data, cb_context, model, conv_data::ConvergenceData, params::MILPModelParameters)
		# 迭代计数器自增
		iteration_counter[] += 1
		# 只处理找到可行解的情况
		if cb_context != COPT.COPT_CBCONTEXT_MIPSOL
			return
		end
		
		current_iter = iteration_counter[]
		
		# 获取全局最优下界（分支定界算法维护的下界）
		best_bnd_ref = Ref{Cdouble}(0.0)
		COPT.COPT_GetCallbackInfo(cb_data, COPT.COPT_CBINFO_BESTBND, best_bnd_ref)
		best_bnd = best_bnd_ref[]
		
		# 获取全局最优上界（当前最优可行解目标值）
		best_obj_ref = Ref{Cdouble}(0.0)
		COPT.COPT_GetCallbackInfo(cb_data, COPT.COPT_CBINFO_BESTOBJ, best_obj_ref)
		best_obj = best_obj_ref[]
		
		# 检查是否有可行解
		has_incumbent_ref = Ref{Cint}(0)
		COPT.COPT_GetCallbackInfo(cb_data, COPT.COPT_CBINFO_HASINCUMBENT, has_incumbent_ref)
		has_incumbent = has_incumbent_ref[] == 1
		
		if !has_incumbent
			return  # 没有可行解，不处理
		end
		
		# 更新全局上下界并记录历史
		lower_bound_updated = false
		upper_bound_updated = false
		
		if best_bnd > conv_data.global_lower_bound
			conv_data.global_lower_bound = best_bnd
			push!(conv_data.lower_bounds, best_bnd)
			push!(conv_data.lower_iterations, current_iter)
			lower_bound_updated = true
		end
		
		if best_obj < conv_data.global_upper_bound
			conv_data.global_upper_bound = best_obj
			push!(conv_data.upper_bounds, best_obj)
			push!(conv_data.upper_iterations, current_iter)
			upper_bound_updated = true
		end
		
		# 只有在上下界任一更新时才绘图
		if !lower_bound_updated && !upper_bound_updated
			return
		end
		
		# 获取可行解的变量值
		COPT.load_callback_variable_primal(cb_data, cb_context)
		
		# 获取底层的 MOI 模型
		moi_model = JuMP.backend(model)
		
		# 获取变量值
		if params.optimizeCapacity
			C_heatpump_idx = JuMP.index(model[:C_heatpump])
			C_boiler_idx = JuMP.index(model[:C_boiler])
			C_heatpump_val = JuMP.MOI.get(moi_model, JuMP.MOI.CallbackVariablePrimal(cb_data), C_heatpump_idx)
			C_boiler_val = JuMP.MOI.get(moi_model, JuMP.MOI.CallbackVariablePrimal(cb_data), C_boiler_idx)
		else
			C_heatpump_val = params.C_heatpump
			C_boiler_val = params.C_boiler
		end
		Ts_idx = JuMP.index.(model[:Ts])
		Ts_val = JuMP.MOI.get.(moi_model, JuMP.MOI.CallbackVariablePrimal(cb_data), Ts_idx)
		
		# 更新最优解记录
		conv_data.best_obj[] = best_obj
		conv_data.best_C_heatpump[] = C_heatpump_val
		conv_data.best_C_boiler[] = C_boiler_val
		conv_data.best_Ts[] = Ts_val
		
		# 计算蓄热容量和时长
		C_storage_val = heatStorageCapacity
		heat_storage_hours = C_heatpump_val > 0 ? C_storage_val : 0
		
		if upper_bound_updated && length(conv_data.best_obj) > 1
			println("iter:", current_iter, "upper_bound:",round(conv_data.upper_bounds[end-1],digits=2), "->",round(best_obj,digits=2))
			println("C_heatpump: ", C_heatpump_val)
			println("hour: ", heat_storage_hours)
			println("C_boiler: ", C_boiler_val)
		end
		if lower_bound_updated && length(conv_data.best_obj) > 1
			println("iter:", current_iter, "lower_bound:",round(conv_data.lower_bounds[end-1],digits=2), "->",round(best_obj,digits=2))
		end

		# 创建子图布局：2行1列
		layout = @layout [a; b]
		
		# 上图：温度曲线
		p1 = plot(t_list, Ts_val,
			title = "Ts vs t",
			xlabel = "t (h)",
			ylabel = "Ts (℃)",
			legend = false,
			color = :blue,
			linewidth = 2,
			ylims = (Tsmin, Tsmax)
		)
		
		# 下图：收敛曲线
		p2 = plot(
			xlabel = "iterations",
			ylabel = "objective value",
			title = "Convergence Curve",
			legend = :bottomright
		)
		if !isempty(conv_data.lower_bounds)
			plot!(p2, conv_data.lower_iterations, conv_data.lower_bounds,
				label = "lower bound",
				color = :red,
				linewidth = 2
			)
		end
		if !isempty(conv_data.upper_bounds)
			plot!(p2, conv_data.upper_iterations, conv_data.upper_bounds,
				label = "upper bound",
				color = :green,
				linewidth = 2
			)
		end
		
		# 合并两张图并显示
		plt = plot(p1, p2, layout = layout, size = (800, 600))
		display(plt)
	end
end
sysVariables = HeatPumpWithStorageSystem.SysVariables(milp_params.heatLoad, milp_params.segmentTariff,1e-2)

# DP参数: n=温度系数量, q=单位热份数, solver_type=求解器类型, dt=时间步长
dp_params = HeatPumpWithStorageSystem.DP_INITIAL_PARAMS(10, Int(1/dt), :Exhaustive, 0.5)

initial,initial_cost, Ts_list = generateInitialSolution_DP(dp_params, milp_params, sysVariables)

model = generate_model(PressedWaterOneStorageOneCompressor_MILP(), milp_params)
set_attribute(model, "TimeLimit", 10)
set_attribute(model, "MipStartMode", 2)
#set_objective_function(model, 0.0)
#fix(model[:Ts][5], 120.0,force=true)
#initial.Ts[:] = precise_solution.Ts[:]
@time result, model = solve_model(PressedWaterOneStorageOneCompressor_MILP(), model, milp_params;
	initial_solution = initial,
	#callback = (cb_data, cb_context, model) -> incumbent_callback(cb_data, cb_context, model, convergence_data, milp_params)
)
# 测试
#=
	fix(model[:Ts][5], 120.0,force=true)
	fix(model[:Ts][1], 220.0,force=true)
	#unfix(model[:Ts][1])
	set_attribute(model, "TimeLimit", 3)

	optimize!(model)
	value.(model[:Ts])
	isFeasible = primal_status(model) in [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]
=#



# 输出结果
println("\n========== MILP求解结果 ==========")
println("求解状态: ", result.isFeasible ? "可行" : "不可行")

if result.isFeasible
	println("\n目标函数值: ", result.objective)
	println("运行成本 (元/天): ", result.C_operation)
	println("正则项: ", result.C_norm)
	println("初投资成本 (元): ", result.C_initial)
	println("总现值 (元): ", result.C_total)

	println("\n设备容量:")
	println("热泵容量 (kW): ", result.C_heatpump_opt)
	println("电锅炉容量 (kW): ", result.C_boiler_opt)
	println("蓄热容量 (kWh): ", C_storage_value)

	println("\n各时段运行状态: ", result.states)

	println("\n功率分布:")
	println("热泵直接供热功率 P1 (kW): ", result.P1)
	println("蓄热供热功率 P2 (kW): ", result.P2)
	println("热泵储热功率 P3 (kW): ", result.P3)
	println("电加热补热功率 P_el (kW): ", result.P_el)
	println("电加热储热功率 P_es (kW): ", result.P_es)

	println("\n蓄热温度 (℃): ", result.Ts)

	println("\nCOP值:")
	println("COP1: ", result.COP1)
	println("COP2: ", result.COP2)
	println("COP3: ", result.COP3)
	println("COPw: ", result.COPw)
else
	println("模型不可行，请检查参数设置！")
end

# ========== 阶段3: DP精确解 ==========
# 使用连续COP函数，以MILP解为起点进行局部搜索改善
if result.isFeasible
    println("\n========== 阶段3: DP精确解求解 ==========")

    # 从MILP解中提取温度轨线
    Ts_milp = result.Ts

    # 生成各时刻的Tair列表
    # Tair_list[t] 对应 t_list[t] 时刻的环境温度
    # 注意: Tair[t] 对应 t_list[t+1] 时刻（即时段t的末尾）
    # 这里用时段起始和结束的Tair
    Tair_list = Float64[]
    for t in 1:length(t_list)
        # t_list[t] 时刻对应的Tair
        # Tair[1] 对应 t_list[2]=dt 时刻，Tair[2] 对应 t_list[3]=2dt 时刻...
        # t_list[1]=0, 对应第一个Tair需要特殊处理
        if t == 1
            push!(Tair_list, Tair[1])  # 起始时刻用第一个Tair
        else
            push!(Tair_list, Tair[t-1])
        end
    end

    # 构建sysVariables（如果没有定义）
    if !@isdefined(sysVariables)
        sysVariables = HeatPumpWithStorageSystem.HSOneStorageOneCompressorMILP.SysVariables(
            segmentHeatLoad, segmentTariff
        )
    end

    # DP精确解参数
    dp_precise_params = HeatPumpWithStorageSystem.DP_PRECISE_PARAMS(
        dt,    # 时间步长
        3,     # 初始局部状态数
        5,     # 最大局部状态数
        4.0,   # 初始温度步长 (℃)
        0.2,   # 最小温度步长 (℃)
        50,    # 最大迭代次数
        4.0    # y_s6: 两次s6最小间隔 (h)，0表示无约束
    )

    # 求解DP精确解
    precise_solution = HeatPumpWithStorageSystem.solvePreciseDP(
        Ts_milp,
        dp_precise_params,
        designParameters,
        milp_params,
        sysVariables,
        Tair_list
    )

    # 输出结果
    println("\nDP精确解结果:")
    println("收敛状态: ", precise_solution.converged ? "已收敛" : "未收敛")
    println("总成本: ", round(precise_solution.cost, digits=4))
    println("温度轨线: ", round.(precise_solution.Ts, digits=2))
    println("状态序列: ", precise_solution.states)
    println("P1 (kW): ", round.(precise_solution.P1, digits=2))
    println("P2 (kW): ", round.(precise_solution.P2, digits=2))
    println("P3 (kW): ", round.(precise_solution.P3, digits=2))
    println("Pe_l (kW): ", round.(precise_solution.Pe_l, digits=2))
    println("Pe_s (kW): ", round.(precise_solution.Pe_s, digits=2))
end
plot([initial.Ts result.Ts precise_solution.Ts],label=["DP initial" "MILP" "DP accurate"]) |> display
println("\n========== 求解结束 ==========")


