
# 使用经济性优化需要在项目主目录下切换到calculations环境
#=
]
activate calculations
=#
begin
using Pkg
#Pkg.activate("calculations")
using Plots
using HeatPumpWithStorageSystem
import BlackBoxOptim as bbo
using DataFrames, CSV
using CoolProp
using Dates
using JSON3  # ← 新增：用于结构化数据保存

include(joinpath(pwd(),"temp","plottool.jl"))	# 引入绘图函数
include(joinpath(pwd(),"temp","ParameterDocGenerator.jl"))# 引入参数文档生成函数

#=
封装了用于设计优化的函数
发布分支design_optimize
=#

situation = "situation27_工况4"
#第一步，指定设计条件变量
#项目设计条件
begin
	# 电价，按照等时间隔输入时段起始时刻的电费。如果需要1小时的间隔描述，就为24个；需要半小时就为48个
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

	Tair = vcat(
		fill(26.0, 7),
		fill(27.0, 7),
		fill(26.0, 11),
	)

	#=
	heatConsumptionPower = vcat(
		fill(0.0, 16),
		fill(1.0, 8),
		fill(0.0, 2),
		fill(1.0, 8),
		fill(0.0, 14),
	)
	=#

	heatConsumptionPower = ones(48)

	# 系数
	maxCOP = 21.0						# 最大COP
	eta_s = 0.7							# 绝热效率
	workingStartHour = 0                # 生产开始时间，已无效，这个参数不用修改
	workingHours = 24                   # 每日工作小时数，，已无效，这个参数不用修改
	TWaste = 85.0                     	# 废热源温度℃，是低温热泵的蒸发温度
	TCompressorIn = 115.0				# 水蒸气压缩机吸气温度℃
	maxTcHigh = 180.0					# 水蒸气压缩机最大排气温度℃
	dT_EvaporationStandard = 5.0		# 传热温差℃
	Tsmin = 120.0						# 蓄热最低温度℃
	Tsmax = 220.0						# 蓄热最低温度℃
	dTRecycleSupply=5.0					# 从蓄热到工厂供热的热回收温差℃
    dTRecycleBackward=5.0				# 从工厂回收的冷凝水到蓄热储热的热回收温差℃

	# 计算参数
	dT = 1.0		# 运行优化计算的温度曲线精度
	dt = 0.5		# 运行优化计算的时间间隔
	smoother = 1e-8		# 正则化参数，不用改

	# 用热温度℃
	# Tuse = 150.0
	Tuse = PropsSI("T","P",0.45e6,"Q",0,"water")-273.15
	
	# 选择低温工质与高温工质。R1233zdE_Water表示低温工质是R1233zdE，高温工质是水。
	# 支持的全部工质有：R134a_Water,NH3_Water,R1233zdE_Water
	refrigerant = R1233zdE_Water

	# 系统结构参数，第一个参数表示允许蓄热与热泵同时供热
	# 第二个参数表示允许蓄热到工厂的热回收
	# 第三个表示允许工厂冷凝水到蓄热的热回收
	sysStruct = RecycleStruct(1,1,1)
end
# 经济性参数条件
begin
	# 设备成本
	waterCompresorCost=1200.0	# 水蒸气压缩机单位供热功率成本 元/kW
	lowHPCost = 1450.0			# 低温热泵单位供热功率成本 元/kW
	elecHeaterCost = 1000.0		# 电极锅炉单位供热功率成本 元/kW
	storageCost = 450e4/1000	# 承压水蓄热单位体积成本 元/m³
	
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
	lifeYears = 15
end

#从这里往下直接执行即可
end
# 第二步，生成设计参数输入结构体
designInput = DesignOptimizeInput(;
	hourlyTariff = hourlyTariff,
	Tair=Tair,
	heatConsumptionPower=heatConsumptionPower,
	maxCOP=maxCOP,
	eta_s=eta_s,
	workingStartHour=workingStartHour,
	workingHours=workingHours,
	
	Tuse=Tuse,
	TWaste=TWaste,
	TCompressorIn=TCompressorIn,
	maxTcHigh=maxTcHigh,
	dT_EvaporationStandard=dT_EvaporationStandard,
	Tsmin=Tsmin,
	Tsmax=Tsmax,
	dTRecycleSupply=dTRecycleSupply,
	dTRecycleBackward=dTRecycleBackward,
	dT=dT,
	dt=dt,
	smoother=smoother,
	sysStruct=sysStruct,
	refrigerant=refrigerant
)
# 第三步，生成设计参数常数结构体
designParameters = generateDesignOptimizeParameters(PressedWaterOneStorageOneCompressor(),designInput)
# 第四步，生成带优化的目标函数
optimizeFunction,getCount,getControlResult,getParams = generateOperationFunction(PressedWaterOneStorageOneCompressor(),designParameters,designInput)

# 第五步，计算一些经济性参数
begin
	p=1/(1+discountRate)
	# 折线系数
	pc=(1-p^lifeYears)/discountRate
	finalWaterCompresorCost=waterCompresorCost*(heatpumpInstallCoff+heatpumpAnnualCost*p*(1-p^lifeYears)/(1-p))
	finalLowHPCost = lowHPCost*(1-1/designParameters.COPWater(TCompressorIn+dT_EvaporationStandard,Tuse))*(heatpumpInstallCoff+heatpumpAnnualCost*p*(1-p^lifeYears)/(1-p))

	finalElecHeaterCost = elecHeaterCost*(elecHeaterInstallCoff+elecHeaterAnnualCost*p*(1-p^lifeYears)/(1-p))

	# 把蓄热的立方米造价转换成蓄热时长造价
	# 1立方米的蓄热能量除以3600秒，得到1立方够用多久
	hour_per_m3=1*900*4.275*(Tsmax-Tuse)/3600
	storageHourCost = storageCost/hour_per_m3
	finalStorageCost = storageHourCost*(storageInstallCoff+storageAnnualCost*p*(1-p^lifeYears)/(1-p))
end

fp=FinanceParameters(
	Plow_cost = finalLowHPCost,       # 低温热泵价格 (¥/kW)
	Phigh_cost = finalWaterCompresorCost,      # 高温热泵价格 (¥/kW)
	Pelec_cost = finalElecHeaterCost,       # 电锅炉价格 (¥/kW)
	Storage_cost = finalStorageCost,     # 蓄热设备价格 (¥/kWh)
	Life_years = lifeYears,          # 设备寿命 (年)
	annual_days = annualDays,        # 年运行天数
	Discount_rate = discountRate,      # 折现率 
)

generateDesignDoc(situation, designInput, designParameters, fp, joinpath(pwd(),"calculations",situation,situation*"_设计参数.md"))

# 现在可以调用函数进行优化了
# 1.7594
heatPumpServiceCoff = 0.2
heatStorageCapacity = 0.0
maxheatStorageInputHour=8.0
@time result = optimizeFunction(
	heatPumpServiceCoff,    # 热泵服务系数
	heatStorageCapacity,    # 蓄热容量
	maxheatStorageInputHour    # 蓄热电加热储满时长
)
plt = operation_result_plot(
	0:dt:24,
	hourly_tariff_ori,
	result;
	w=0.45
)


# 工况列表：
heatPumpServiceCoffList = 0.0:0.2:1.2
heatStorageCapacityList = 0.0:2:10.0
maxheatStorageInputHourList = 1.0:2.0:7.0

case_list = [(i,j,k) for i in heatPumpServiceCoffList, j in heatStorageCapacityList, k in maxheatStorageInputHourList]

# 建文件夹
file_path0 = joinpath(pwd(),"calculations",situation)
if !isdir(file_path0)
	mkdir(file_path0)
end
for hs in heatStorageCapacityList
	# 使用 round 确保文件夹名称的一致性，避免浮点精度问题
	hs_str = string(round(hs, digits=1))
	file_path1 = joinpath(file_path0,"storage_"*hs_str)
	if !isdir(file_path1)
		mkdir(file_path1)
	end
end

num_threads = Threads.nthreads()-1
total_cases = length(case_list)
println("🚀 启动高效任务队列模式（原子索引）")
println("   线程数: $(num_threads)")
println("   总算例数: $(total_cases)")
println("   预计开始时间: $(now())")

# 创建线程安全的计数器
const completed_counter = Threads.Atomic{Int}(0)
const failed_counter = Threads.Atomic{Int}(0)
const task_index = Threads.Atomic{Int}(0)  # ← 原子任务索引，指向下一个待处理的任务

# 定义任务函数
function task()
	while true
		# 获取下一个任务索引
		index = Threads.atomic_add!(task_index, 1) + 1  # ← +1 修正: atomic_add! 返回自增前的值
		if index > total_cases
			break
		end
		
		# 获取当前任务
		case = case_list[index]
		
		# 准备文件路径和名称（无需锁保护）
		hs_str = string(round(case[2], digits=1))
		p = joinpath(file_path0,"storage_"*hs_str)
		filename_base = "$(round(case[1], digits=1))_$(hs_str)_$(round(case[3], digits=1))"
		json_filepath = joinpath(p, "$(filename_base).json")

		start_time_task = time()
		
		try
			# 执行优化
			result = optimizeFunction(case...)
			
			# 计算经济性参数
			annualOperationCost = result[1] * annualDays
			pw,capitalCost,annuity_pv_factor = totalPresentWorth(
				PressedWaterOneStorageOneCompressor(),
				fp,# 经济参数
				result[1],          # 每日运行成本（元）
				case...
			)
			compute_time = time() - start_time_task
			
			# ✅ 构建JSON数据结构（无需锁保护）
			result_dict = Dict(
				"status" => "success",
				"caseParameters" => Dict(
					"heatPumpServiceCoff" => case[1],
					"heatStorageCapacity" => case[2],
					"maxHeatStorageInputHour" => case[3]
				),
				"operationResults" => Dict(
					"totalRealCost" => result[1],
					"TsList" => result[2],
					"P1List" => result[3],
					"P2List" => result[4],
					"P3List" => result[5],
					"PeList" => result[6],
					"realCostList" => result[7]
				),
				"economicResults" => Dict(
					"dailyOperationCost" => round(result[1], digits=3),
					"annualOperationCost" => round(annualOperationCost, digits=3),
					"capitalCost" => round(capitalCost, digits=3),
					"totalPresentWorth" => round(pw, digits=3),
					"annuityPvFactor" => round(annuity_pv_factor, digits=3),
					"computeTimeSeconds" => round(compute_time, digits=2)
				),
				"metadata" => Dict(
					"timestamp" => string(now()),
					"situation" => situation,
					"dt" => dt
				)
			)
			
			# ✅ JSON文件写入（无需锁保护，不同文件操作系统处理并发）
			try
				open(json_filepath, "w") do f
					JSON3.pretty(f, result_dict)  # 格式化输出便于阅读
				end
			catch e_save
				@warn "JSON保存失败: $(json_filepath)" exception=e_save
			end
			
			# 增加成功计数并输出进度
			completed = Threads.atomic_add!(completed_counter, 1) + 1
			println("✅ [$(completed)/$(total_cases)] 完成 | 耗时: $(round(compute_time, digits=1))s | 热泵=$(case[1]), 蓄热=$(case[2]), 时长=$(case[3])")

		catch e
			# 捕获并记录错误
			failed = Threads.atomic_add!(failed_counter, 1) + 1
			
			# ✅ 构建失败标记的JSON结构
			error_dict = Dict(
				"status" => "failed",
				"caseParameters" => Dict(
					"heatPumpServiceCoff" => case[1],
					"heatStorageCapacity" => case[2],
					"maxHeatStorageInputHour" => case[3]
				),
				"errorInfo" => Dict(
					"errorType" => string(typeof(e)),
					"errorMessage" => sprint(showerror, e),
					"stackTrace" => sprint(io -> showerror(io, e, catch_backtrace())),
					"possibleCauses" => (
						if occursin("NaN", sprint(showerror, e)) || occursin("Inf", sprint(showerror, e))
							["参数组合导致数学表达式无意义（如除以零、Inf-Inf）",
							 "蓄热容量过小无法平衡负荷",
							 "热泵服务系数不合理"]
						elseif occursin("UndefVarError", sprint(showerror, e))
							["代码中存在未定义的变量或函数",
							 "可能是模块导入问题"]
						else
							["需要查看堆栈跟踪进一步分析"]
						end
					)
				),
				"metadata" => Dict(
					"timestamp" => string(now()),
					"situation" => situation
				)
			)
			
			# ✅ 保存失败信息的JSON文件
			try
				open(json_filepath, "w") do f
					JSON3.pretty(f, error_dict)
				end
			catch e2
				@warn "错误JSON保存失败: $(json_filepath)" exception=e2
			end
			
			println("❌ [$(completed + failed)/$(total_cases)] 失败 | 热泵=$(case[1]), 蓄热=$(case[2]), 时长=$(case[3])")
			println("   错误: ", sprint(showerror, e))
		end
	end
end

# 启动线程
println("\n📊 启动 $(num_threads) 个工作者线程...\n")
start_time = time()
Threads.@threads for _ in 1:num_threads
	task()
end

# 打印结果
completed = completed_counter[]
failed = failed_counter[]
elapsed_time = time() - start_time

println("\n" * "="^60)
println("🎉 计算完成！")
println("="^60)
println("   成功算例: $(completed)")
println("   失败算例: $(failed)")
println("   结束时间: $(now())")
println("   总耗时: $(round(elapsed_time, digits=2))秒")
println("="^60)

# ============================================================================
# 批量绘图阶段 - 在所有计算完成后执行
# ============================================================================
println("\n" * "="^60)
println("🎨 开始批量绘图...")
println("="^60)

# ✅ 直接使用全局电价曲线,无需重新定义
# hourly_tariff_ori 已在脚本开头定义(第30-42行)

# 批量绘图函数(避免软作用域问题)
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
				# 读取JSON
				data = JSON3.read(read(json_path, String))
				
				# 检查状态
				if data.status == "failed"
					total_skipped += 1
					continue
				end
				
				# 提取运行结果
				result = [
					data.operationResults.totalRealCost,
					data.operationResults.TsList,
					data.operationResults.P1List,
					data.operationResults.P2List,
					data.operationResults.P3List,
					data.operationResults.PeList,
					data.operationResults.realCostList
				]
				
				# 生成图片路径
				png_path = joinpath(folder_path, replace(json_file, ".json" => ".png"))
				
				# 绘图并保存
				plt = operation_result_plot(
					0:dt:24,
					hourly_tariff_ori,  # ← 使用全局电价曲线
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

# 执行批量绘图
total_plotted, total_skipped, total_failed_plot = batch_plot_results()

println("\n" * "="^60)
println("🎉 批量绘图完成！")
println("="^60)
println("   成功绘制: $(total_plotted)")
println("   跳过(失败算例): $(total_skipped)")
println("   绘图失败: $(total_failed_plot)")
println("="^60)

# ============================================================================
# 数据整理阶段 - 按蓄热时长分组绘制堆叠折线图
# ============================================================================
println("\n" * "="^60)
println("📊 开始绘制经济性分析图表...")
println("="^60)

# 数据整理与绘图函数
function plot_economic_analysis()
	# 收集所有成功算例的数据
	all_data = []
	
	storage_folders = filter(x -> isdir(joinpath(file_path0, x)), readdir(file_path0))
	
	for folder in storage_folders
		folder_path = joinpath(file_path0, folder)
		json_files = filter(x -> endswith(x, ".json"), readdir(folder_path))
		
		for json_file in json_files
			json_path = joinpath(folder_path, json_file)
			
			try
				data = JSON3.read(read(json_path, String))
				
				if data.status == "success"
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
	
	println("✅ 共收集 $(length(all_data)) 个成功算例")
	
	# 按(蓄热容量, max input hour)分组
	grouped_by_config = Dict{Tuple{Float64, Float64}, Vector{Any}}()
	
	for item in all_data
		key = (item.storageCapacity, item.maxInputHour)
		if !haskey(grouped_by_config, key)
			grouped_by_config[key] = []
		end
		push!(grouped_by_config[key], item)
	end
	
	println("📈 共有 $(length(grouped_by_config)) 个不同的配置组")
	
	# 为每个配置组合绘制一张图
	
	# 先计算所有数据中的最大总现值(用于统一纵轴)
	max_total_PW = maximum([item.totalPresentWorth for item in all_data])
	
	for ((storageCap, maxHour), items) in sort(collect(grouped_by_config), by=x->(x[1][1], x[1][2]))
		# 按热泵容量排序
		sort!(items, by=x->x.heatPumpCapacity)
		
		# 提取数据
		heatPumpCapacities = [item.heatPumpCapacity for item in items]
		capitalCosts = [item.capitalCost for item in items]
		totalPWs = [item.totalPresentWorth for item in items]
		
		# 计算运行费用现值(用于堆叠)
		operatingCostPW = totalPWs .- capitalCosts
		
		# 创建堆叠面积图
		plt = plot(
			size=(800, 600),
			grid=true,
			xlabel="Heat Pump Capacity",
			ylabel="Cost (CNY)",
			title="Economic Analysis\nStorage = $(round(storageCap, digits=1)), Max Input Hour = $(round(maxHour, digits=1)) h",
			legend=:topleft,
			ylim=(0, max_total_PW * 1.05)  # ← 设置统一的纵轴上限(留5%余量)
		)
		
		# 绘制初投资(底层)
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
		
		# 绘制运行费用现值(堆叠在初投资之上)
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
		
		# 保存图片
		storage_str = string(round(storageCap, digits=1))
		hour_str = string(round(maxHour, digits=1))
		png_path = joinpath(file_path0, "economic_storage_$(storage_str)_hour_$(hour_str).png")
		savefig(plt, png_path)
		
		println("✅ 已保存: economic_storage_$(storage_str)_hour_$(hour_str).png")
	end

	# ============================================================================
	# 新增: 按蓄热容量分组,绘制不同max input hour的总现值对比图
	# ============================================================================
	println("\n📊 生成总现值对比图 (不同Max Input Hour)...")

	# 按蓄热容量重新分组
	grouped_by_storage_only = Dict{Float64, Vector{Any}}()

	for item in all_data
		storageCap = item.storageCapacity
		if !haskey(grouped_by_storage_only, storageCap)
			grouped_by_storage_only[storageCap] = []
		end
		push!(grouped_by_storage_only[storageCap], item)
	end

	# 为每个蓄热容量绘制一张对比图
	for (storageCap, items) in sort(collect(grouped_by_storage_only), by=x->x[1])
		# 按max input hour分组
		hour_groups = Dict{Float64, Vector{Any}}()
		for item in items
			hour = item.maxInputHour
			if !haskey(hour_groups, hour)
				hour_groups[hour] = []
			end
			push!(hour_groups[hour], item)
		end
		
		# 创建对比图
		plt = plot(
			size=(800, 600),
			grid=true,
			xlabel="Heat Pump Capacity",
			ylabel="Total Present Worth (CNY)",
			title="Total Present Worth Comparison\nStorage Capacity = $(round(storageCap, digits=1))",
			legend=:topleft,
			ylim=(0, max_total_PW * 1.05)
		)
		
		# 为每个max input hour绘制一条曲线
		color_palette = [:blue, :red, :green, :orange, :purple, :brown, :pink, :gray]
		color_idx = 1
		
		for (hour, hour_items) in sort(collect(hour_groups), by=x->x[1])
			# 按热泵容量排序
			sort!(hour_items, by=x->x.heatPumpCapacity)
			
			# 提取数据
			heatPumpCapacities = [item.heatPumpCapacity for item in hour_items]
			totalPWs = [item.totalPresentWorth for item in hour_items]
			
			# 绘制曲线
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
		
		# 保存图片
		storage_str = string(round(storageCap, digits=1))
		png_path = joinpath(file_path0, "total_PW_comparison_storage_$(storage_str).png")
		savefig(plt, png_path)
		
		println("✅ 已保存: total_PW_comparison_storage_$(storage_str).png (包含 $(length(hour_groups)) 条曲线)")
	end

	return length(grouped_by_config)
end

# 执行数据整理与绘图
num_charts = plot_economic_analysis()

println("\n" * "="^60)
println("🎉 经济性分析图表绘制完成！")
println("="^60)
println("   生成图表数: $(num_charts)")
println("   保存位置: $(file_path0)")
println("="^60)

# ============================================================================
# 最优算例统计分析
# ============================================================================
println("\n" * "="^60)
println("🏆 最优算例统计分析")
println("="^60)

# 重新收集所有成功算例的数据(用于统计)
all_success_data = []

storage_folders = filter(x -> isdir(joinpath(file_path0, x)), readdir(file_path0))

for folder in storage_folders
	folder_path = joinpath(file_path0, folder)
	json_files = filter(x -> endswith(x, ".json"), readdir(folder_path))
	
	for json_file in json_files
		json_path = joinpath(folder_path, json_file)
		
		try
			data = JSON3.read(read(json_path, String))
			
			if data.status == "success"
				push!(all_success_data, (
					fileName = json_file,
					heatPumpCapacity = data.caseParameters.heatPumpServiceCoff,
					storageCapacity = data.caseParameters.heatStorageCapacity,
					maxInputHour = data.caseParameters.maxHeatStorageInputHour,
					capitalCost = data.economicResults.capitalCost,
					totalPresentWorth = data.economicResults.totalPresentWorth,
					dailyOperationCost = data.economicResults.dailyOperationCost,
					annualOperationCost = data.economicResults.annualOperationCost,
					computeTime = data.economicResults.computeTimeSeconds
				))
			end
		catch e
			@warn "读取JSON失败: $(json_file)" exception=e
		end
	end
end

if length(all_success_data) > 0
	# 找出总现值最小的算例
	min_pw_idx = argmin([item.totalPresentWorth for item in all_success_data])
	best_case = all_success_data[min_pw_idx]
	
	println("\n✅ 找到最优算例 (最小总现值):")
	println("-"^60)
	println("   文件名: $(best_case.fileName)")
	println("   配置参数:")
	println("     • Heat Pump Capacity:    $(best_case.heatPumpCapacity)")
	println("     • Storage Capacity:      $(best_case.storageCapacity)")
	println("     • Max Input Hour:        $(best_case.maxInputHour)")
	println()
	println("   经济性指标:")
	println("     • Capital Cost:          ¥$(round(best_case.capitalCost, digits=2))")
	println("     • Total Present Worth:   ¥$(round(best_case.totalPresentWorth, digits=2)) ← 最小值")
	println("     • Daily Operation Cost:  ¥$(round(best_case.dailyOperationCost, digits=2))")
	println("     • Annual Operation Cost: ¥$(round(best_case.annualOperationCost, digits=2))")
	println()
	println("   计算性能:")
	println("     • Compute Time:          $(round(best_case.computeTime, digits=2)) seconds")
	println("-"^60)
	
	
else
	println("\n⚠️  未找到任何成功算例,无法进行统计分析")
end

println("\n" * "="^60)
println("🎊 全部分析完成！")
println("="^60)
