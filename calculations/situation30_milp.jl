
# 使用经济性优化需要在项目主目录下切换到calculations环境
#=
situation30用于计算不同的热泵容量、电加热功率、存储容量下的运行成本，不需要回调过程，将每种情况的运行优化结果整理成图和json文件，保存到对应目录。具体参考`calculations\situation27_运行结果展示.jl`
此外还有系统设计的参数要输出md文件，参考situation27的150行左右designParameters = generateDesignOptimizeParameters(PressedWaterOneStorageOneCompressor(),designInput)
=#
using Pkg

#Pkg.activate("calculations")
using Plots
using HeatPumpWithStorageSystem

using DataFrames, CSV
using CoolProp
using COPT, JuMP
using JSON3

include(joinpath(pwd(), "tools", "plottool.jl"))
include(joinpath(pwd(), "tools", "ParameterDocGenerator.jl"))

#=
封装了用于设计优化的函数
发布分支design_optimize
=#

situation = "situation30"
file_path0 = joinpath(pwd(), "calculations", situation)
if !isdir(file_path0)
	mkdir(file_path0)
end
#第一步，指定设计条件变量
#项目设计条件
# 跳过heatPumpServiceCoff+maxheatStorageInputHour < 1的工况
heatPumpServiceCoff_list = 0.4:0.2:1.2						# 5
heatStorageCapacity_list = 2.0:1.0:8						# 7
maxheatStorageInputHour_list = [0.5,1.0,1.5,2.5,3.5,4.5]	# 6

for hs in heatStorageCapacity_list
	# 使用 round 确保文件夹名称的一致性，避免浮点精度问题
	hs_str = string(round(hs, digits = 1))
	file_path1 = joinpath(file_path0, "storage_" * hs_str)
	if !isdir(file_path1)
		mkdir(file_path1)
	end
end


begin
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

	heatConsumptionPower = ones(48)

	Tair = fill(85.0, 49)

	dt = 0.5

	# 系数
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
	#Tuse = PropsSI("T", "P", 0.45e6, "Q", 0, "water") - 273.15
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

	dt_list = repeat(dt_list, inner = inner_divide) / inner_divide
	segmentHeatLoad = repeat(segmentHeatLoad, inner = inner_divide)
	segmentTariff = repeat(segmentTariff, inner = inner_divide)
	segmentTair = repeat(segmentTair, inner = inner_divide)

	segmentTariff_ori = segmentTariff/0.7393

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
end

# 第三步，计算部分常数
begin
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

	generateDesignDoc(situation, designInput, designParameters, fp, joinpath(pwd(), "calculations", situation, situation * "_设计参数.md"))

case_count = 0
n_cases = length(heatStorageCapacity_list) * length(heatPumpServiceCoff_list) * length(maxheatStorageInputHour_list)
for heatStorageCapacity in heatStorageCapacity_list
	for heatPumpServiceCoff in heatPumpServiceCoff_list
		for maxheatStorageInputHour in maxheatStorageInputHour_list
			global case_count, n_cases
			case_count += 1

			result_dir = joinpath(pwd(), "calculations", "situation30", "storage_$(round(heatStorageCapacity,digits=1))", "$(round(heatPumpServiceCoff,digits=1))_$(round(heatStorageCapacity,digits=1))_$(round(maxheatStorageInputHour,digits=1)).json")
            if isfile(result_dir)
                data = JSON3.read(read(result_dir, String))
                if data.status == "success" && data.operationResults.gap <= 0.4
                    println("算例  $(case_count)/$(n_cases)  ", round(heatPumpServiceCoff, digits=1), " ", round(heatStorageCapacity, digits=1), " ", round(maxheatStorageInputHour, digits=1))
                    continue
                end
            end

			params = getParams(heatPumpServiceCoff, heatStorageCapacity, maxheatStorageInputHour)
			cpm = params.cpm# kWh/K
			n = length(dt_list)


			M_5 = heatPumpServiceCoff * maximum(heatConsumptionPower)

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
				COPdesign = map(i -> designParameters.COPOverlapFunction(designParameters.TairFunction(t_list[i]) - dT_EvaporationStandard, Tuse), 1:n),
				lambdaNorm = lambda_norm,  # 注意: 使用lambda_norm变量
				smoother = smoother,

				# 设备容量参数
				C_heatpump = heatPumpServiceCoff * maximum(heatConsumptionPower),
				C_boiler = params.PeMax,
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

			# 方式一：生成初值（全程热泵供热）
			initial = generateInitialSolution_HeatPumpOnly(milp_params)
			model = generate_model(PressedWaterOneStorageOneCompressor_MILP(), milp_params)

			@time result, model = solve_model(
				PressedWaterOneStorageOneCompressor_MILP(),
				model,
				milp_params;
				initial_solution = initial,
				#callback = (cb_data, cb_context, model) -> incumbent_callback(cb_data, cb_context, model, convergence_data, milp_params)
			)

			println("算例  $(case_count)/$(n_cases)  ",round(heatPumpServiceCoff, digits=1)," ",round(heatStorageCapacity, digits=1)," ",round(maxheatStorageInputHour, digits=1))
			# 写入结果
			if result.isFeasible
				annualOperationCost = result.C_operation * annualDays
				pw, capitalCost, annuity_pv_factor = totalPresentWorth(
					PressedWaterOneStorageOneCompressor(),
					fp,
					result.C_operation,
					heatPumpServiceCoff,
					heatStorageCapacity,
					maxheatStorageInputHour
				)
				
				pw = isinf(pw) || isnan(pw) ? nothing : round(pw, digits=3)
				capitalCost = isinf(capitalCost) || isnan(capitalCost) ? nothing : round(capitalCost, digits=3)
				annualOperationCost = isinf(annualOperationCost) || isnan(annualOperationCost) ? nothing : round(annualOperationCost, digits=3)
				
				function safe_val(x)
					if x isa AbstractVector
						return [isinf(v) || isnan(v) ? nothing : v for v in x]
					else
						return isinf(x) || isnan(x) ? nothing : x
					end
				end
				
				filename_base = "$(round(heatPumpServiceCoff, digits=1))_$(round(heatStorageCapacity, digits=1))_$(round(maxheatStorageInputHour, digits=1))"
				json_filepath = joinpath(pwd(), "calculations", situation, "storage_" * string(round(heatStorageCapacity, digits=1)), "$(filename_base).json")
				
				result_dict = Dict(
					"status" => "success",
					"caseParameters" => Dict(
						"heatPumpServiceCoff" => heatPumpServiceCoff,
						"heatStorageCapacity" => heatStorageCapacity,
						"maxHeatStorageInputHour" => maxheatStorageInputHour
					),
					"operationResults" => Dict(
						"objective" => safe_val(result.objective),
						"C_operation" => safe_val(result.C_operation),
						"C_norm" => safe_val(result.C_norm),
						"C_initial" => safe_val(result.C_initial),
						"C_total" => safe_val(result.C_total),
						"C_heatpump_opt" => safe_val(result.C_heatpump_opt),
						"C_boiler_opt" => safe_val(result.C_boiler_opt),
						"C_storage" => C_storage_value,
						"states" => safe_val(result.states),
						"P1" => safe_val(result.P1),
						"P2" => safe_val(result.P2),
						"P3" => safe_val(result.P3),
						"P_el" => safe_val(result.P_el),
						"P_es" => safe_val(result.P_es),
						"Ts" => safe_val(result.Ts),
						"COP1" => safe_val(result.COP1),
						"COP2" => safe_val(result.COP2),
						"COP3" => safe_val(result.COP3),
						"COPw" => safe_val(result.COPw),
						"isFeasible" => result.isFeasible,
						"gap" => safe_val(result.gap),
						"best_bound" => safe_val(result.best_bound)
					),
					"economicResults" => Dict(
						"dailyOperationCost" => safe_val(round(result.C_operation, digits=3)),
						"annualOperationCost" => annualOperationCost,
						"capitalCost" => capitalCost,
						"totalPresentWorth" => pw
					)
				)
				
				open(json_filepath, "w") do f
					write(f, JSON3.write(result_dict))
				end
			else
				filename_base = "$(round(heatPumpServiceCoff, digits=1))_$(round(heatStorageCapacity, digits=1))_$(round(maxheatStorageInputHour, digits=1))"
				json_filepath = joinpath(pwd(), "calculations", situation, "storage_" * string(round(heatStorageCapacity, digits=1)), "$(filename_base).json")
				
				result_dict = Dict(
					"status" => "infeasible",
					"caseParameters" => Dict(
						"heatPumpServiceCoff" => heatPumpServiceCoff,
						"heatStorageCapacity" => heatStorageCapacity,
						"maxHeatStorageInputHour" => maxheatStorageInputHour
					),
					"operationResults" => Dict(
						"isFeasible" => false
					),
					"economicResults" => Dict()
				)
				
				open(json_filepath, "w") do f
					write(f, JSON3.write(result_dict))
				end
				
				println("模型不可行，请检查参数设置！")
				println("结果已保存到: ", json_filepath)
			end
		end
	end
end

# 绘图

function batch_plot_results()
	storage_folders = filter(x -> isdir(joinpath(file_path0, x)), readdir(file_path0))
	
	total_plotted = 0
	total_skipped = 0
	total_failed_plot = 0
	
	for folder in storage_folders
		folder_path = joinpath(file_path0, folder)
		json_files = filter(x -> endswith(x, ".json"), readdir(folder_path))
		
		for json_file in json_files
			json_path = joinpath(folder_path, json_file)
			
			try
				data = JSON3.read(read(json_path, String))
				
				if data.status == "infeasible"
					total_skipped += 1
					continue
				end
				
				result = [
					data.operationResults.objective,
					data.operationResults.Ts,
					data.operationResults.P1,
					data.operationResults.P2,
					data.operationResults.P3,
					data.operationResults.P_el,
					zeros(length(data.operationResults.P1))
				]
				
				png_path = joinpath(folder_path, replace(json_file, ".json" => ".png"))
				
				n_time_points = length(data.operationResults.Ts)
				n_segments = length(data.operationResults.P1)
				
				time_list = t_list
				tariff_list = segmentTariff_ori
				
				plt = operation_result_plot(
					time_list,
					tariff_list,
					result;
					w=0.45
				)
				savefig(plt, png_path)
				
				total_plotted += 1
				
			catch e
				@warn "绘图失败: $(json_file)" exception=e
				total_failed_plot += 1
			end
		end
	end
	
	return total_plotted, total_skipped, total_failed_plot
end

total_plotted, total_skipped, total_failed_plot = batch_plot_results()

println("\n" * "="^60)
println("🎉 批量绘图完成！")
println("="^60)
println("   成功绘制: $(total_plotted)")
println("   跳过(不可行算例): $(total_skipped)")
println("   绘图失败: $(total_failed_plot)")
println("="^60)

function plot_economic_analysis()
	all_data = []
	
	storage_folders = filter(x -> isdir(joinpath(file_path0, x)), readdir(file_path0))
	
	for folder in storage_folders
		folder_path = joinpath(file_path0, folder)
		json_files = filter(x -> endswith(x, ".json"), readdir(folder_path))
		
		for json_file in json_files
			json_path = joinpath(folder_path, json_file)
			
			try
				data = JSON3.read(read(json_path, String))
				
				if data.status == "success" && data.economicResults.totalPresentWorth !== nothing && data.economicResults.capitalCost !== nothing
					push!(all_data, (
						heatPumpCapacity = data.caseParameters.heatPumpServiceCoff,
						storageCapacity = data.caseParameters.heatStorageCapacity,
						maxInputHour = data.caseParameters.maxHeatStorageInputHour,
						capitalCost = data.economicResults.capitalCost,
						totalPresentWorth = data.economicResults.totalPresentWorth
					))
				end
			catch e
				@warn "读取JSON失败: $(json_file)" exception=e
			end
		end
	end
	
	println("\n✅ 共收集 $(length(all_data)) 个成功算例")
	
	grouped_by_config = Dict{Tuple{Float64, Float64}, Vector{Any}}()
	
	for item in all_data
		key = (item.storageCapacity, item.maxInputHour)
		if !haskey(grouped_by_config, key)
			grouped_by_config[key] = []
		end
		push!(grouped_by_config[key], item)
	end
	
	println("📈 共有 $(length(grouped_by_config)) 个不同的配置组")
	
	max_total_PW = maximum([item.totalPresentWorth for item in all_data])
	
	for ((storageCap, maxHour), items) in sort(collect(grouped_by_config), by=x->(x[1][1], x[1][2]))
		sort!(items, by=x->x.heatPumpCapacity)
		
		heatPumpCapacities = [item.heatPumpCapacity for item in items]
		capitalCosts = [item.capitalCost for item in items]
		totalPWs = [item.totalPresentWorth for item in items]
		
		operatingCostPW = totalPWs .- capitalCosts
		
		plt = plot(
			size=(800, 600),
			grid=true,
			xlabel="Heat Pump Capacity",
			ylabel="Cost (CNY)",
			title="Economic Analysis\nStorage = $(round(storageCap, digits=1)), Max Input Hour = $(round(maxHour, digits=1)) h",
			legend=:topleft,
			ylim=(0, max_total_PW * 1.05)
		)
		
		plot!(plt,
			heatPumpCapacities,
			capitalCosts,
			label="Capital Cost",
			lw=2,
			color=:blue,
			fillrange=zeros(length(capitalCosts)),
			fillalpha=0.5,
			fillcolor=:blue
		)
		
		plot!(plt,
			heatPumpCapacities,
			totalPWs,
			label="Operating Cost PW",
			lw=2,
			color=:red,
			fillrange=capitalCosts,
			fillalpha=0.5,
			fillcolor=:red
		)
		
		storage_str = string(round(storageCap, digits=1))
		hour_str = string(round(maxHour, digits=1))
		png_path = joinpath(file_path0, "economic_storage_$(storage_str)_hour_$(hour_str).png")
		savefig(plt, png_path)
		
		println("✅ 已保存: economic_storage_$(storage_str)_hour_$(hour_str).png")
	end

	println("\n📊 生成总现值对比图 (不同Max Input Hour)...")

	grouped_by_storage_only = Dict{Float64, Vector{Any}}()

	for item in all_data
		storageCap = item.storageCapacity
		if !haskey(grouped_by_storage_only, storageCap)
			grouped_by_storage_only[storageCap] = []
		end
		push!(grouped_by_storage_only[storageCap], item)
	end

	for (storageCap, items) in sort(collect(grouped_by_storage_only), by=x->x[1])
		hour_groups = Dict{Float64, Vector{Any}}()
		for item in items
			hour = item.maxInputHour
			if !haskey(hour_groups, hour)
				hour_groups[hour] = []
			end
			push!(hour_groups[hour], item)
		end
		
		plt = plot(
			size=(800, 600),
			grid=true,
			xlabel="Heat Pump Capacity",
			ylabel="Total Present Worth (CNY)",
			title="Total Present Worth Comparison\nStorage Capacity = $(round(storageCap, digits=1))",
			legend=:topleft,
			ylim=(0, max_total_PW * 1.05)
		)
		
		color_palette = [:blue, :red, :green, :orange, :purple, :brown, :pink, :gray]
		color_idx = 1
		
		for (hour, hour_items) in sort(collect(hour_groups), by=x->x[1])
			sort!(hour_items, by=x->x.heatPumpCapacity)
			
			heatPumpCapacities = [item.heatPumpCapacity for item in hour_items]
			totalPWs = [item.totalPresentWorth for item in hour_items]
			
			plot!(plt,
				heatPumpCapacities,
				totalPWs,
				label="$(round(hour, digits=1)) h",
				lw=2,
				marker=:circle,
				markersize=4,
				color=color_palette[(color_idx - 1) % length(color_palette) + 1]
			)
			
			color_idx += 1
		end
		
		storage_str = string(round(storageCap, digits=1))
		png_path = joinpath(file_path0, "total_PW_comparison_storage_$(storage_str).png")
		savefig(plt, png_path)
		
		println("✅ 已保存: total_PW_comparison_storage_$(storage_str).png (包含 $(length(hour_groups)) 条曲线)")
	end

	return length(grouped_by_config)
end

num_charts = plot_economic_analysis()
