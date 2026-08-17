

# 使用经济性优化需要在项目主目录下切换到calculations环境
#=
situation34: 基于situation32计算多个省市不同分时电价下的结果
=#
using Pkg

#Pkg.activate("calculations")
using Plots
using HeatPumpWithStorageSystem

using DataFrames, CSV
using CoolProp
using COPT, JuMP
using JSON3
using Dates

include(joinpath(pwd(), "tools", "plottool.jl"))
include(joinpath(pwd(), "tools", "ParameterDocGenerator.jl"))
# 分时电价曲线定义与计算顺序（grid_price.jl 在 src 子目录，未被主模块加载，需显式 include）
include(joinpath(pwd(), "src", "systemModels", "GridPrice", "grid_price.jl"))

situation = "situation34"
file_path0 = joinpath(pwd(), "calculations", situation)
if !isdir(file_path0)
	mkdir(file_path0)
end

# 多参数遍历设置
heatPumpServiceCoff_list = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
heatStorageCapacity_list = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
maxheatStorageInputHour_list = [1.0, 2.0, 3.0, 4.0]
grid_price_list = build_grid_price_order()

#=
heatPumpServiceCoff_list = [0.4]
heatStorageCapacity_list = [3.0]
maxheatStorageInputHour_list = [0.5]
=#

# 继续计算标志：跳过已完成的算例
continue_calculate = true
stage_milp = false

# 算例计数器（全局，供内层循环 global 引用）
case_count = 0
n_cases = length(heatStorageCapacity_list) * length(heatPumpServiceCoff_list) * length(maxheatStorageInputHour_list) * length(grid_price_list)

# ===== 已完成电价对象进度管理（CSV 持久化，按城市+季节组织） =====
progress_csv = joinpath(pwd(), "calculations", situation, "completed_cases.csv")
if continue_calculate && isfile(progress_csv)
	completed_df = CSV.read(progress_csv, DataFrame)     # 续算：加载历史进度
else
	completed_df = DataFrame(                             # 全新计算：清空重建
		city            = String[],
		season          = String[],
		grid_price_name = String[],
		status          = String[],               # success / skipped / failed
		timestamp       = String[],
	)
end

"""
记录一条电价对象进度并写回CSV。
按 (city, season) 去重（每城市每季节1条）；跨季对象为其覆盖的每个季节各写一行。
"""
function record_progress!(df, city, season, gp_name, status)
	exists = any((df.city .== city) .& (df.season .== season))
	exists && return
	push!(df, (city, season, gp_name, status, string(now())))
	CSV.write(progress_csv, df)
end

# 预处理：创建文件夹
for gp in grid_price_list
	file_path_grid = joinpath(file_path0, gp.name)
	if !isdir(file_path_grid)
		mkdir(file_path_grid)
	end
	for hs in heatStorageCapacity_list
		hs_str = string(round(hs, digits = 1))
		file_path1 = joinpath(file_path_grid, "storage_" * hs_str)
		if !isdir(file_path1)
			mkdir(file_path1)
		end
	end
end
for gp in grid_price_list
	# 第一步，指定设计条件变量
	begin
		hourly_tariff_ori = gp.price
		dt = gp.time_gap
		n = gp.n

		baseElectricityPrice = 1.0                     # 计算用基准电价（研究电价影响时按1计算）
		analysis_base_price = 0.6                   # 绘图用基准电价
		hourlyTariff = hourly_tariff_ori * baseElectricityPrice
		# 需求曲线固定给定24段（与situation31一致）；
		# 电价精度为0.5h(n=48，江苏)时按 inner=2 展开为48段，使电价与热负荷时段一一对应
		heatConsumptionPower = repeat(vcat(
				fill(0.0, 1),
				fill(1.0, 5),
			), outer = 4)
		if n == 48
			heatConsumptionPower = repeat(heatConsumptionPower, inner = 2)
		end
		Tair = fill(50.0, length(heatConsumptionPower) + 1)

		# 系数
		y_s6 = 4.0   # 两次s6最小间隔 (h)，0表示无约束
		maxCOP = 21.0
		eta_s = 0.7
		workingStartHour = 0
		workingHours = 24
		TWaste = 85.0
		TCompressorIn = 115.0
		maxTcHigh = 180.0
		dT_EvaporationStandard = 5.0
		Tsmin = 120.0
		Tsmax = 220.0
		dTRecycleSupply = 5.0
		dTRecycleBackward = 5.0

		dT = 0.01
		smoother = 1e-8
		Tuse = 150.0
		refrigerant = R1233zdE_Water
		sysStruct = RecycleStruct(1, 0, 0)

		segmentHeatLoad = heatConsumptionPower[1:end]
		segmentTariff = hourlyTariff[1:end]
		segmentTair = Tair[1:end-1]
		dt_list = ones(length(segmentHeatLoad)) * dt

		t_list = vcat(0.0, cumsum(dt_list))

		t_list = vcat(0.0, cumsum(dt_list))

		n_segments = length(dt_list)

		lambda_norm = 1e-4

		T_e = segmentTair .- dT_EvaporationStandard
	end

	# 经济性参数条件
	begin
		waterCompresorCost = 1200.0
		lowHPCost = 1450.0
		elecHeaterCost = 1000.0
		storageCost = 450e4 / 1000

		heatpumpInstallCoff = 2
		elecHeaterInstallCoff = 1.2
		storageInstallCoff = 1.2

		heatpumpAnnualCost = 0.1
		elecHeaterAnnualCost = 0.05
		storageAnnualCost = 0.05

		discountRate = 0.032
		annualDays = 300
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

		designParameters = generateDesignOptimizeParameters(PressedWaterOneStorageOneCompressor(), designInput)

		COP_ca = map(T -> designParameters.COPOverlapFunction(T, Tuse), T_e)

		optimizeFunction, getCount, getControlResult, getParams = generateOperationFunction(PressedWaterOneStorageOneCompressor(), designParameters, designInput)
	end

	# 第三步，计算部分常数
	begin
		# COP分档信息（与situation31一致，2°C间隔）
		T_s_1g_list = vcat(120.0:2.0:185.0, 220.0)
		T_s_2g_list = vcat(120.0:2.0:185.0, 220.0)
		T_s_3g_list = vcat(120.0:2.0:185.0, 220.0)
		T_s_wg_list = vcat(120.0:2.0:185.0, 220.0)

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
		)
	end

	# 第五步，计算一些经济性参数
	begin
		local p = 1 / (1 + discountRate)   # local: 避免与 grid_price.jl 的全局 p 软作用域歧义
		pc = (1 - p^lifeYears) / discountRate
		finalWaterCompresorCost = waterCompresorCost * (heatpumpInstallCoff + heatpumpAnnualCost * p * (1 - p^lifeYears) / (1 - p))
		finalLowHPCost = lowHPCost * (1 - 1 / designParameters.COPWater(TCompressorIn + dT_EvaporationStandard, Tuse)) * (heatpumpInstallCoff + heatpumpAnnualCost * p * (1 - p^lifeYears) / (1 - p))

		finalElecHeaterCost = elecHeaterCost * (elecHeaterInstallCoff + elecHeaterAnnualCost * p * (1 - p^lifeYears) / (1 - p))

		hour_per_m3 = 1 * 900 * 4.275 * (Tsmax - Tuse) / 3600
		storageHourCost = storageCost / hour_per_m3
		finalStorageCost = storageHourCost * (storageInstallCoff + storageAnnualCost * p * (1 - p^lifeYears) / (1 - p))
	end

	# 第六步，生成双层优化的目标函数
	fp = FinanceParameters(
		Plow_cost = finalLowHPCost,
		Phigh_cost = finalWaterCompresorCost,
		Pelec_cost = finalElecHeaterCost,
		Storage_cost = finalStorageCost,
		Life_years = lifeYears,
		annual_days = annualDays,
		Discount_rate = discountRate,
	)

	generateDesignDoc(situation, designInput, designParameters, fp, joinpath(pwd(), "calculations", situation, gp.name, situation * gp.name * "_设计参数.md"))

	# 电价对象级进度统计（Dict避免嵌套for作用域shadow问题）
	stats = Dict("success" => 0, "skipped" => 0, "failed" => 0)

	for heatStorageCapacity in heatStorageCapacity_list
		for heatPumpServiceCoff in heatPumpServiceCoff_list
			for maxheatStorageInputHour in maxheatStorageInputHour_list
				global case_count, n_cases
				case_count += 1

				hs_str = string(round(heatStorageCapacity, digits = 1))
				result_dir = joinpath(pwd(), "calculations", situation, gp.name, "storage_$(hs_str)", "$(round(heatPumpServiceCoff,digits=1))_$(hs_str)_$(round(maxheatStorageInputHour,digits=1)).json")

				if isfile(result_dir) && continue_calculate
					data = JSON3.read(read(result_dir, String))
					if data.status == "success"
						# 被跳过的已完成算例：计入统计（对象级汇总时写进度）
						stats["skipped"] += 1
						println("跳过 算例 $(case_count)/$(n_cases): $(round(heatPumpServiceCoff,digits=1))_$(hs_str)_$(round(maxheatStorageInputHour,digits=1))")
						continue
					end
				end

				println("\n========== 算例 $(case_count)/$(n_cases): heatPumpServiceCoff=$(round(heatPumpServiceCoff,digits=1)), heatStorageCapacity=$(hs_str), maxheatStorageInputHour=$(round(maxheatStorageInputHour,digits=1)) ==========")

				try
					params = getParams(heatPumpServiceCoff, heatStorageCapacity, maxheatStorageInputHour)
				cpm = params.cpm
				n = length(dt_list)

				M_5 = heatPumpServiceCoff * maximum(heatConsumptionPower)
				C_storage_value = heatStorageCapacity

				# 创建MILP模型参数
				milp_params = MILPModelParameters(
					Tuse = Tuse,
					Tsmin = Tsmin,
					Tsmax = Tsmax,
					Tcmax = maxTcHigh,
					dTs = dTRecycleSupply,
					Te = T_e[1],
					cpm = cpm,
					dt = dt,
					COPmax = maxCOP,
					COPdesign = map(i -> designParameters.COPOverlapFunction(designParameters.TairFunction(t_list[i]) - dT_EvaporationStandard, Tuse), 1:n),
					lambdaNorm = lambda_norm,
					smoother = smoother, C_heatpump = heatPumpServiceCoff * maximum(heatConsumptionPower),
					C_boiler = params.PeMax,
					C_storage = C_storage_value,
					optimizeCapacity = false, c_heatpump = finalLowHPCost,
					c_storage = finalStorageCost,
					c_boiler = finalElecHeaterCost, D = annualDays,
					r = discountRate,
					N = lifeYears, hourlyTariff = hourlyTariff,
					heatLoad = segmentHeatLoad,
					Tair = Tair,
					dt_list = dt_list, n_segments = length(dt_list),
					segmentTariff = segmentTariff, m1 = m1,
					COP1v = COP1v,
					T1g = T1g, m2 = m2,
					COP2v = COP2v,
					T2g = T2g, m3 = m3,
					COP3v = COP3v,
					T3g = T3g, mw = mw,
					COPwv = COPwv,
					Twg = Twg, COPca = COP_ca, M1 = M_1,
					M2 = M_2,
					M3 = M_3,
					M4 = M_4,
					M5 = M_5, solver = :COPT,
				)

				# ========== 阶段1: DP初始解 ==========
				println("\n---------- 阶段1: DP初始解 ----------")
				sysVariables = HeatPumpWithStorageSystem.SysVariables(milp_params.heatLoad, milp_params.segmentTariff, 1e-2)
				dp_params = HeatPumpWithStorageSystem.DP_INITIAL_PARAMS(20, Int(5 / dt), :Exhaustive, dt, y_s6)

				initial = nothing
				initial_cost = 9999.0
				try
					@time initial, initial_cost, Ts_list = generateInitialSolution_DP(dp_params, milp_params, sysVariables;
						cop_mode = :continuous, designParameters = designParameters,
					)
					println("DP初始解温度: ", round.(initial.Ts[1:min(5, end)], digits = 2), "...")
				catch e
					println("DP初始解生成失败: ", e)
				end

				# ========== 阶段2: MILP求解 ==========
				if stage_milp
					println("\n---------- 阶段2: MILP求解 ----------")
					model = generate_model(PressedWaterOneStorageOneCompressor_MILP(), milp_params)
					set_attribute(model, "TimeLimit", 60)
					set_attribute(model, "MipStartMode", 2)

					@time result, model = solve_model(PressedWaterOneStorageOneCompressor_MILP(), model, milp_params;
						initial_solution = initial,
					)

					println("MILP求解状态: ", result.isFeasible ? "可行" : "不可行")
					if result.isFeasible
						println("MILP目标函数值: ", round(result.objective, digits = 4))
						println("MILP运行成本: ", round(result.C_operation, digits = 4))
					end
				end
				# ========== 阶段3: DP精确解 ==========
				sysVariables_precise = HeatPumpWithStorageSystem.SysVariables(milp_params.heatLoad, milp_params.segmentTariff, 1e-2)
				precise_solution = nothing
				println("\n---------- 阶段3: DP精确解 ----------")
				Ts_milp = stage_milp && result.isFeasible ? result.Ts : initial.Ts

				Tair_list = Float64[]
				for t in 1:length(t_list)
					if t == 1
						push!(Tair_list, Tair[1])
					else
						push!(Tair_list, Tair[t-1])
					end
				end

				dp_precise_params = HeatPumpWithStorageSystem.DP_PRECISE_PARAMS(
					dt,    # 时间步长
					5,     # 初始局部状态数
					5,     # 最大局部状态数
					4.0,   # 初始温度步长 (℃)
					0.2,   # 最小温度步长 (℃)
					80,    # 最大迭代次数
					y_s6,    # y_s6: 两次s6最小间隔 (h)，0表示无约束
				)

				try
					@time precise_solution = HeatPumpWithStorageSystem.solvePreciseDP(
						Ts_milp,
						dp_precise_params,
						designParameters,
						milp_params,
						sysVariables_precise,
						Tair_list,
					)
					println("DP精确解收敛: ", precise_solution.converged ? "是" : "否")
					println("DP精确解成本: ", round(precise_solution.cost, digits = 4))
				catch e
					println("DP精确解求解失败: ", e)
				end

				# 经济性参数
				annualOperationCost = precise_solution.cost * annualDays
				pw, capitalCost, operatingPV, annuity_pv_factor = totalPresentWorth(
					PressedWaterOneStorageOneCompressor(),
					fp,
					precise_solution.cost,
					heatPumpServiceCoff,
					heatStorageCapacity,
					maxheatStorageInputHour,
				)

				# 保存结果
				result_dict = Dict(
				"status" => initial_cost < 999.0 ? "success" : "infeasible",
				# 电价曲线信息（price 为相对倍数，不乘基准电价；season 为所属季节轮）
				"gridPrice" => Dict(
					"name" => gp.name,                                # 电价曲线名称
					"city" => gp.city,                                # 城市
					"price" => gp.price,                              # 分时电价向量（相对倍数，长度 n）
					"timeGap" => gp.time_gap,                         # 时间间隔 (h)
					"n" => gp.n,                                      # 时段数
					"season" => gp.season,                            # 所属季节（逗号分隔，如 "summer,winter"）
				),
				"heatPumpServiceCoff" => heatPumpServiceCoff,
					"heatStorageCapacity" => heatStorageCapacity,
					"maxheatStorageInputHour" => maxheatStorageInputHour,
					"initial" =>
						initial !== nothing ?
						Dict(
							"cost" => initial_cost,
							"Ts" => initial.Ts,
							"states" => [findfirst(x -> x > 0.5, initial.s[:, i]) for i in 1:length(initial.P_el)],
							"P1" => [sum(initial.P_k_s[k, i, 1] + initial.P_k_e[k, i, 1] for k ∈ 1:8) for i in 1:length(initial.P_el)],
							"P2" => [sum(initial.P_k_s[k, i, 2] + initial.P_k_e[k, i, 2] for k ∈ 1:8) for i in 1:length(initial.P_el)],
							"P3" => [sum(initial.P_k_s[k, i, 3] + initial.P_k_e[k, i, 3] for k ∈ 1:8) for i in 1:length(initial.P_el)],
							"P_el" => initial.P_el,
							"P_es" => initial.P_es,
							"COP1e" => initial.COP1e,
							"COP2e" => initial.COP2e,
							"COP3e" => initial.COP3e,
							"COPwe" => initial.COPwe,
						) : nothing,
					"milp" =>
						stage_milp && result.isFeasible ?
						Dict(
							"objective" => result.objective,
							"C_operation" => result.C_operation,
							"C_initial" => result.C_initial,
							"C_total" => result.C_total,
							"Ts" => result.Ts,
							"states" => result.states,
							"P1" => result.P1,
							"P2" => result.P2,
							"P3" => result.P3,
							"P_el" => result.P_el,
							"P_es" => result.P_es,
						) : nothing,
					"dp_precise" =>
						precise_solution !== nothing ?
						Dict(
							"converged" => precise_solution.converged,
							"cost" => precise_solution.cost,
							"Ts" => precise_solution.Ts,
							"states" => precise_solution.states,
							"P1" => precise_solution.P1,
							"P2" => precise_solution.P2,
							"P3" => precise_solution.P3,
							"Pe_l" => precise_solution.Pe_l,
							"Pe_s" => precise_solution.Pe_s,
						) : nothing,
					"economicResults" => Dict(
						"basePrice" => baseElectricityPrice,   # 计算用基准电价
						"operatingPV" => operatingPV,          # 运行成本现值（与基准电价成正比）
						"dailyOperationCost" => round(precise_solution.cost, digits = 3),
						"annualOperationCost" => annualOperationCost,
						"capitalCost" => capitalCost,
						"totalPresentWorth" => pw,
					),
				)

				open(result_dir, "w+") do f
					JSON3.write(f, result_dict)
				end

				# 算例计算成功：计入统计（对象级汇总时写进度）
				stats["success"] += 1

				println("\n结果已保存: ", result_dir)
				catch e
					# 算例失败：计入统计并继续下一个算例
					println("算例失败: ", e)
					stats["failed"] += 1
				end
			end
		end
	end

	# 对象算例全部处理完毕：汇总写入进度
	# status规则: 有failed→failed；否则有success→success；否则全跳过→skipped
	obj_status = stats["failed"] > 0 ? "failed" :
				 stats["success"] > 0 ? "success" : "skipped"
	# 跨季对象为其覆盖的每个季节各写一行；首个=实际计算季节，其余季节行标记为skipped（被覆盖）
	seasons = split(gp.season, ",")
	for (i, s) in enumerate(seasons)
		row_status = (i == 1) ? obj_status : "skipped"
		record_progress!(completed_df, gp.city, s, gp.name, row_status)
	end

	# ========== 每条电价曲线计算完成后立即绘图（与 situation32 单条电价曲线一致） ==========
	# 各绘图函数假设 file_path0 下存在 storage_* 文件夹，
	# 因此对每个电价对象传入其子目录 {gp.name}/ 即可（plottool.jl 无需修改）
	println("\n" * "="^60)
	println("绘图: $(gp.name) ...")
	println("="^60)

	gp_dir = joinpath(file_path0, gp.name)

	# 绘图参数（与计算口径一致：按 gp 精度展开）
	hourly_tariff_ori = gp.price
	segmentTariff_ori = hourly_tariff_ori[1:end]               # 峰谷时段判断用（相对倍数）
	t_list = vcat(0.0, cumsum(ones(length(segmentTariff_ori)) .* gp.time_gap))

	# 1. DP精确解运行结果绘图
	println("\n--- DP精确解运行结果绘图 ---")
	dp_plotted, dp_skipped, dp_failed = batch_plot_operation_results(
		gp_dir, t_list, segmentTariff_ori;
		result_key = "dp_precise",
		cost_field = "cost",
		pe_l_field = "Pe_l",
		pe_s_field = "Pe_s",
		output_suffix = "_dp",
		output_dir = "dp_precise",
	)
	println("DP精确解: 成功=$(dp_plotted), 跳过=$(dp_skipped), 失败=$(dp_failed)")

	# 2. 各时段能耗模式构成堆叠柱状图（5层：P1/P2/P3/Pe_l/Pe_s）
	println("\n--- 各时段能耗模式构成堆叠柱状图 ---")
	dp_en_plotted, dp_en_skipped, dp_en_failed = batch_plot_mode_energy_stacked(
		gp_dir, t_list;
		result_key = "dp_precise",
		pe_l_field = "Pe_l",
		pe_s_field = "Pe_s",
		output_suffix = "_dp_energy",
		output_dir = "energy_stacked",
	)
	println("DP精确解能耗堆叠图: 成功=$(dp_en_plotted), 跳过=$(dp_en_skipped), 失败=$(dp_en_failed)")

	# 3. 经济性分析绘图（堆叠图+对比图，输出到子文件夹）
	println("\n--- 经济性分析绘图 ---")
	num_stacked, num_comparison = plot_economic_analysis(
		gp_dir;
		params_key = "",              # 参数在JSON顶层
		plot_modes = [:all],          # 绘制全部模式
		use_subfolders = true,        # 使用子文件夹组织
		base_price = 1.0,             # 目标基准电价
	)
	println("经济性分析: 堆叠图 $(num_stacked) 张，对比图 $(num_comparison) 张")

	# 4. 最小总现值包络分析绘图（对电锅炉容量维度取最小）
	println("\n--- 最小总现值包络分析绘图 ---")
	num_envelope = plot_optimal_envelope(
		gp_dir;
		params_key = "",                   # 参数在JSON顶层
		econ_key = "economicResults",      # 总现值在 economicResults.totalPresentWorth
		total_pw_field = "totalPresentWorth",
		base_price = 1.0,                  # 目标基准电价
	)
	println("包络分析绘图完成！共 $(num_envelope) 张图")

	# 5. 基准电价对最优配置影响分析绘图
	println("\n--- 基准电价对最优配置影响分析绘图 ---")
	num_price = plot_optimal_capacity_vs_price(
		gp_dir;
		params_key = "",                   # 参数在JSON顶层
		econ_key = "economicResults",      # 经济数据在 economicResults
		base_prices = 0.58:0.02:2.0,       # 基准电价范围
	)
	println("基准电价影响分析绘图完成！共 $(num_price) 张图")

	# 6. 最优运行成本绘图（横轴=蓄热容量，每条曲线=固定热泵容量）
	println("\n--- 最优运行成本绘图 ---")
	num_oc = plot_optimal_operating_cost(
		gp_dir;
		params_key = "",                   # 参数在JSON顶层
		econ_key = "economicResults",      # 经济数据在 economicResults
		cost_field = "dailyOperationCost", # 每日运行成本
		base_price = 1.0,                  # 目标基准电价（运行成本=日运行成本×基准电价）
	)
	println("最优运行成本绘图完成！共 $(num_oc) 张图")

end
println("\n========== 所有算例计算完成 ==========")
