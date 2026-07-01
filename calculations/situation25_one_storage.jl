
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

include(joinpath(pwd(),"temp","plottool.jl"))	# 引入绘图函数
include(joinpath(pwd(),"temp","ParameterDocGenerator.jl"))# 引入参数文档生成函数

#=
封装了用于设计优化的函数
发布分支design_optimize
=#

situation = "situation25_电锅炉加蓄热"
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

	#=
	heatConsumptionPower = vcat(
		fill(0.0, 16),
		fill(1.0, 8),
		fill(0.0, 2),
		fill(1.0, 8),
		fill(0.0, 14),
	)
	=#

	heatConsumptionPower = vcat(
		fill(0.0, 16),
		fill(1.0, 8),
		fill(0.0, 2),
		fill(1.0, 8),
		fill(0.0, 14),
	)
	heatConsumptionPower = ones(48)
	#heatStorageCapacity = 2.0,      # 蓄热量kWh(相变蓄热)
	#PheaterMax::Real = 1.0,               # 热泵最大功率kW

	timePoint = 0:0.5:24 |> collect        # 
	heatStorageOutEfficiency = 0.95    # 蓄热释放效率
	heatStorageInEfficiency = 0.95      # 蓄热充能效率
	heatStorageVelocity = 1.0
	
	Tsmax = 220.0
	Tsmin = 150.0
end
# 经济性参数条件
begin
	# 设备成本
	elecHeaterCost = 1000.0		# 电极锅炉单位供热功率成本 元/kW
	storageCost = 450e4/1000	# 承压水蓄热单位体积成本 元/m³
	
	# 安装费系数
	elecHeaterInstallCoff = 1.2
	storageInstallCoff = 1.2

	# 维保费用
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
osp = generateSystemCoff(OneStorage();
    hourlyTariff=hourlyTariff,
    timePoint=timePoint,
    load=heatConsumptionPower,
    heatStorageOutEfficiency=heatStorageOutEfficiency,    # 蓄热释放效率
    heatStorageInEfficiency=heatStorageInEfficiency,      # 蓄热充能效率
    heatStorageVelocity=heatStorageVelocity
)

#generateAndSolve(OneStorage(),osp)

# 第三步，生成带优化的目标函数,现在可以调用函数进行优化了
optimizeFunction = generateOperationFunction(OneStorage(),osp)

result = optimizeFunction(0.0,0.0)
plot([result.heaterPower result.heatStorage[1:end-1]])

# 第四步，计算一些经济性参数
begin
	p=1/(1+discountRate)
	# 折线系数
	pc=(1-p^lifeYears)/discountRate
	
	finalElecHeaterCost = elecHeaterCost*(elecHeaterInstallCoff+elecHeaterAnnualCost*p*(1-p^lifeYears)/(1-p))

	# 把蓄热的立方米造价转换成蓄热时长造价
	# 1立方米的蓄热能量除以3600秒，得到1立方够用多久
	hour_per_m3=1*900*4.275*(Tsmax-Tsmin)/3600
	storageHourCost = storageCost/hour_per_m3
	finalStorageCost = storageHourCost*(storageInstallCoff+storageAnnualCost*p*(1-p^lifeYears)/(1-p))
end

# 第六步，生成双层优化的目标函数
fp=OneStorageFinanceParameters(
	Pelec_cost = finalElecHeaterCost,       # 电锅炉价格 (¥/kW)
	Storage_cost = finalStorageCost,     # 蓄热设备价格 (¥/kWh)
	Life_years = lifeYears,          # 设备寿命 (年)
	annual_days = annualDays,        # 年运行天数
	Discount_rate = discountRate,      # 折现率 
)

# 生成环境温度向量（如果没有定义，使用默认值）
if !@isdefined(Tair)
    Tair = fill(25.0, length(hourlyTariff))
end

generateDesignDoc(situation, hourlyTariff, heatConsumptionPower, Tair, Tsmin, Tsmax, fp, joinpath(pwd(),"calculations",situation,situation*"_设计参数.md"))

bb_cost = get_bb_cost(OneStorage(),optimizeFunction,fp)

# 第七步，制定收敛性检查的回调函数
# 创建自定义的实时绘图对象
mutable struct FitnessPlot
    fitness_data::Vector{Float64}
	solution_data::Vector{Vector{Float64}}
    iteration_count::Int
	evaluation_count::Int
	#best_fitness::Float64
	#best_solution::Vector{Float64}#无法从优化结果中获取当前解对应的控制参数
    verbose::Bool
end
function monitor_callback(fp::FitnessPlot, opt_controller)
	
	# 当前评估次数
	current_evals = bbo.num_func_evals(opt_controller)
	# 获取当前最佳适应度值（最小值）
    current_best_fitness = bbo.best_fitness(opt_controller)
    # 获取当前最佳解
    current_best_solution = bbo.best_candidate(opt_controller)

	# 记录数据
    push!(fp.fitness_data, current_best_fitness)
	push!(fp.solution_data, current_best_solution)
	fp.iteration_count += 1
	fp.evaluation_count = current_evals
	
	

	if fp.verbose && fp.iteration_count % 50 == 0
		@info "触发回调,最优解：$(round.(current_best_solution,digits=3)), 最佳适应度值：$(current_best_fitness)"
        println("迭代 $(fp.iteration_count): 最小成本 = $(fp.fitness_data[end])")
		plt = plot(fp.fitness_data,xlabel="iteration times",ylabel="least cost",label=:none,title="optimize value vs iteration times")
		display(plt)
    end
end

fitnessPlotController = FitnessPlot([],[[]],0,0,true)

search_range = [
    (0.0, 10.0),  # heatStorageCapacity  kWh
    (0.0, 3.0)]  # PeMax kW

good_guess = [
	[1.0,1.0]
]

@info "开始进行黑盒优化..."
@info "线程数:$(Threads.nthreads())"
@time res = bboptimize(bb_cost,good_guess;
    SearchRange=search_range,
    MaxSteps=1500,      # 最多迭代步数
    NumDimensions=3,
	Method = :adaptive_de_rand_1_bin_radiuslimited,
	PopulationSize = 12,
    TraceInterval=1.0,
    TraceMode=:compact,
	CallbackFunction = oc -> monitor_callback(fitnessPlotController,oc),
	CallbackInterval = 0.0,
	NThreads=4
)

best_x = bbo.best_candidate(res)
best_f = bbo.best_fitness(res)

@info "优化结束"
@info "最优变量" heatStorageCapacity = best_x[1] PeMax = best_x[2]
@info "最小总现值" best_f

#=
[ Info: 优化结束
┌ Info: 最优变量
│   heatStorageCapacity = 0.9500000936010964
└   PeMax = 1.000000094936982     
┌ Info: 最小总现值
└   best_f = 24575.384166055388
=#
result = optimizeFunction(best_x...)
annualOperationCost = result.objective * fp.annual_days
pw,capitalCost,annuity_pv_factor = totalPresentWorth(
	OneStorage(),
	fp,
	result.objective,   # 每日运行成本（元）
	best_x...
)

println("""
日运行费用：$(round(result.objective,digits=3))
年运行费用：$(round(annualOperationCost,digits=3))
初投资费用：$(round(capitalCost,digits=3))
总现值：$(round(pw,digits=3))
折现系数：$(round(annuity_pv_factor,digits=3))
""")