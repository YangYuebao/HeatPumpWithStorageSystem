
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
封装了用于设计优化的函数
发布分支design_optimize
=#

situation = "situation29"
#第一步，指定设计条件变量
#项目设计条件
begin
	heatPumpServiceCoff, heatStorageCapacity, maxheatStorageInputHour = 0.4, 6.0, 3.0
	
	#=
	hourlyTariff = zeros(24)
	hourlyTariff[1:8] .= 99.0094
	hourlyTariff[9:16] .= 0.658
	hourlyTariff[17:24] .= 0.3725

	Tair = fill(85.0,25)

	heatConsumptionPower = ones(24)
	=#
	
	hourly_tariff_ori = ones(48)
	p = 4.7#1.7
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

	# 统计每个时间段长度
	idx = 1
	dt_list = [1]
	segmentHeatLoad = [heatConsumptionPower[1]]
	segmentTariff = [hourlyTariff[1]]
	segmentTair = [Tair[1]]
	for i in 2:length(hourlyTariff)
		global idx
		if hourlyTariff[i] == hourlyTariff[i-1] && heatConsumptionPower[i] == heatConsumptionPower[i-1] && Tair[i] == Tair[i-1]
			dt_list[idx] += 1
		else
			push!(dt_list, 1)
			push!(segmentHeatLoad, heatConsumptionPower[i])
			push!(segmentTariff, hourlyTariff[i])
			push!(segmentTair, Tair[i])
			idx += 1
		end
	end
	dt_list *= dt

	inner_divide = 1

	dt_list = repeat(dt_list, inner = inner_divide)/inner_divide
	segmentHeatLoad = repeat(segmentHeatLoad, inner = inner_divide)
	segmentTariff = repeat(segmentTariff, inner = inner_divide)
	segmentTair = repeat(segmentTair, inner = inner_divide)
	
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
	T_s_1g_list = [120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 220.0]
	T_s_2g_list = [120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 220.0]
	T_s_3g_list = [120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 220.0]
	T_s_wg_list = [120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 220.0]

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
# 方式一：生成初值（全程热泵供热）
initial = generateInitialSolution_HeatPumpOnly(milp_params)
model = generate_model(PressedWaterOneStorageOneCompressor_MILP(), milp_params)
#set_attribute(model, "TimeLimit", 30)
#fix(model[:Ts][5], 120.0,force=true)
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
println("开始求解MILP模型...")
# 调优
#=
	#优化求解参数
	MOI.set(model, MOI.RawOptimizerAttribute("TuneMode"), 0)
	MOI.set(model, MOI.RawOptimizerAttribute("TuneMethod"), 0)
	MOI.set(model, MOI.RawOptimizerAttribute("TuneTimeLimit"), 1200.0)
	MOI.set(model, MOI.RawOptimizerAttribute("TuneOutputLevel"), 2)

	optimizer = JuMP.backend(model)
	prob = optimizer.prob
	COPT.COPT_Tune(prob)
	num_results = Ref{Cint}()
	COPT.COPT_GetIntAttr(prob, "TuneResults", num_results)
	println("调优结果数量: ", num_results[])
	COPT.COPT_LoadTuneParam(prob, 0)

	COPT.COPT_WriteTuneParam(prob, 0, joinpath(pwd(),"calculations",situation, "best_tune.par"))
	println("已保存最佳调优结果到 best_tune.par")
	#COPT.COPT_ReadParam(new_prob, "best_tune.par")
	#println("已从 best_tune.par 加载参数")
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




println("\n========== 求解结束 ==========")

# =============================================================
# heat_k17 诊断函数：检查从参数到辅助变量再到原始变量的完整链路
# =============================================================
function diagnose_heat_k17(model, params, result)
    n = params.n_segments
    m1 = params.m1
    m2 = params.m2
    
    println("\n\n========== heat_k17 详细诊断 ==========")
    println("格式说明：heat_k17[i] = sum_{k=1..7} [COPca[i]*u2[k,i] + sum_{j=1..m1} COP1v[j,i]*v2[k,i,j] + sum_{j=1..m2} COP2v[j,i]*v4[k,i,j]]")
    
    for i in 1:1
        println("\n" * "="^80)
        println("时段 $i / $n")
        println("="^80)
        
        println("\n【参数值】")
        println("COPca[$i] = $(params.COPca[i])")
        println("heatLoad[$i] = $(params.heatLoad[i])")
        println("Ts[$i] = $(result.Ts[i])")
        println("实际状态 = $(result.states[i])")
        
        println("\n【COP分档参数】")
        println("COP1v[:, $i] = $(params.COP1v[:, i])")
        println("COP2v[:, $i] = $(params.COP2v[:, i])")
        println("T1g[$i, :] = $(params.T1g[i, :])")
        println("T2g[$i, :] = $(params.T2g[i, :])")
        
        total_heat_k17 = 0.0
        total_from_u2 = 0.0
        total_from_v2 = 0.0
        total_from_v4 = 0.0
        
        for k in 1:7
            s_val = value(model[:s][k, i])
            P_k_s_val = value(model[:P_k_s][k, i, 1])
            u2_val = value(model[:u2][k, i])
            
            println("\n--- 状态 k=$k ---")
            println("s[$k,$i] = $s_val")
            println("P_k_s[$k,$i,1] (起始功率) = $P_k_s_val")
            println("u2[$k,$i] = s * P_k_s = $u2_val")
            println("  验证: u2 = s * P_k_s = $(s_val * P_k_s_val)")
            
            # COPca贡献
            copca_contrib = params.COPca[i] * u2_val
            total_from_u2 += copca_contrib
            println("  COPca贡献 = COPca[$i] * u2[$k,$i] = $(params.COPca[i]) * $u2_val = $copca_contrib")
            
            # v2贡献 (COP1v)
            v2_contrib = 0.0
            for j in 1:m1
                v2_val = value(model[:v2][k, i, j])
                contrib = params.COP1v[j, i] * v2_val
                v2_contrib += contrib
                if v2_val > 1e-6
                    println("  v2[$k,$i,$j] = $v2_val, COP1v[$j,$i] = $(params.COP1v[j,i])")
                    println("    贡献 = $(params.COP1v[j,i]) * $v2_val = $contrib")
                end
            end
            total_from_v2 += v2_contrib
            println("  COP1v总贡献 = $v2_contrib")
            
            # v4贡献 (COP2v)
            v4_contrib = 0.0
            for j in 1:m2
                v4_val = value(model[:v4][k, i, j])
                contrib = params.COP2v[j, i] * v4_val
                v4_contrib += contrib
                if v4_val > 1e-6
                    println("  v4[$k,$i,$j] = $v4_val, COP2v[$j,$i] = $(params.COP2v[j,i])")
                    println("    贡献 = $(params.COP2v[j,i]) * $v4_val = $contrib")
                end
            end
            total_from_v4 += v4_contrib
            println("  COP2v总贡献 = $v4_contrib")
            
            k_total = copca_contrib + v2_contrib + v4_contrib
            total_heat_k17 += k_total
            
            # 详细分解三项
            u2_contrib = params.COPca[i] * u2_val
            v2_sum = 0.0
            v4_sum = 0.0
            for j in 1:m1
                v2_j = value(model[:v2][k, i, j])
                v2_sum += params.COP1v[j, i] * v2_j
            end
            for j in 1:m2
                v4_j = value(model[:v4][k, i, j])
                v4_sum += params.COP2v[j, i] * v4_j
            end
            
            println("  状态 k=$k 供热分解:")
            println("    第1项: COPca[$i] * u2[$k,$i] = $(params.COPca[i]) * $u2_val = $u2_contrib")
            println("    第2项: sum_{j=1..$m1} COP1v[j,$i] * v2[$k,$i,j] = $v2_sum")
            println("    第3项: sum_{j=1..$m2} COP2v[j,$i] * v4[$k,$i,j] = $v4_sum")
            println("    总供热 = $k_total")
        end
        
        println("\n【时段 $i 汇总】")
        println("heat_k17[$i] (计算值) = $total_heat_k17")
        println("heat_k17[$i] (结果值) = $(result.heat_k17[i])")
        println("heat_k8[$i] = $(result.heat_k8[i])")
        println("P_el[$i] (电加热补热) = $(result.P_el[i])")
        println("总供热 = $(result.heat_total[i])")
        println("热负荷需求 = $(params.heatLoad[i])")
        println("供热余量 = $(result.heat_total[i] - params.heatLoad[i])")
        
        # 验证约束
        if result.heat_total[i] >= params.heatLoad[i] - 1e-6
            println("✓ 热负荷约束满足: heat_total >= heatLoad")
        else
            println("✗ 热负荷约束不满足!")
        end
    end
    
    println("\n" * "="^80)
    println("诊断结束")
    println("="^80)
end

#=
if result.isFeasible
    diagnose_heat_k17(model, milp_params, result)
end
=#

