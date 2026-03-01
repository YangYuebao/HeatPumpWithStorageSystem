
# 使用经济性优化需要在项目主目录下切换到calculations环境
#=
]
activate calculations
=#
using Pkg
#Pkg.activate("calculations")
using Plots
using HeatPumpWithStorageSystem
import BlackBoxOptim as bbo
using DataFrames, CSV
using CoolProp

#=
封装了用于设计优化的函数
发布分支design_optimize
=#

situation = "situation22"
#第一步，指定设计条件变量
#项目设计条件
begin
	# 电价，按照等时间隔输入时段起始时刻的电费。如果需要1小时的间隔描述，就为24个；需要半小时就为48个
	hourly_tariff_ori = ones(48)
	p = 1.7
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

	heatConsumptionPower = vcat(
		fill(0.0, 16),
		fill(1.0, 8),
		fill(0.0, 2),
		fill(1.0, 8),
		fill(0.0, 14),
	)

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
	dT = 0.02		# 运行优化计算的温度曲线精度
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
	lifeYears = 6
end

#从这里往下直接执行即可

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
designParameters = generateDesignOptimizeParameters(designInput)
# 第四步，生成带优化的目标函数
optimizeFunction,getCount,getControlResult,getParams = generateOperationFunction(designParameters,designInput)

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
	finalStorageCost = storageCost*(storageInstallCoff+storageAnnualCost*p*(1-p^lifeYears)/(1-p))
end

# 现在可以调用函数进行优化了

# 1.7594
heatPumpServiceCoff = 1.0
heatStorageCapacity = 1.0
maxheatStorageInputHour=1.0
@time result = optimizeFunction(
	heatPumpServiceCoff,    # 热泵服务系数
	heatStorageCapacity,    # 蓄热容量
	maxheatStorageInputHour    # 蓄热电加热储满时长
)

#=
plt=plot(result[2],xlabel="Hour",ylabel="Temerature ℃",title="Heat Storage Temperature")
savefig(plt,"plots/1.0_8.0_4.0.png")
=#

#=
 Info: 最优变量
│   heatPumpServiceCoff = 0.12675378941319726
│   heatStorageCapacity = 0.5095878178083123
└   maxheatStorageInputHour = 7.96062243329475
┌ Info: 最小总现值
└   best_f = 5649.636713178115
=#

# 第六步，生成双层优化的目标函数
fp=FinanceParameters(
	Plow_cost = finalLowHPCost,       # 低温热泵价格 (¥/kW)
	Phigh_cost = finalWaterCompresorCost,      # 高温热泵价格 (¥/kW)
	Pelec_cost = finalElecHeaterCost,       # 电锅炉价格 (¥/kW)
	Storage_cost = finalStorageCost,     # 蓄热设备价格 (¥/kWh)
	Life_years = lifeYears,          # 设备寿命 (年)
	annual_days = annualDays,        # 年运行天数
	Discount_rate = discountRate,      # 折现率 
)
bb_cost = get_bb_cost(optimizeFunction,fp)

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
	
	@info "触发回调,最优解：$(round.(current_best_solution,digits=3)), 最佳适应度值：$(current_best_fitness)"

	if fp.verbose #&& plot.iteration_count % 10 == 0
        println("迭代 $(fp.iteration_count): 最小成本 = $(fp.fitness_data[end])")
		plt = plot(fp.fitness_data,xlabel="iteration times",ylabel="least cost",label=:none,title="optimize value vs iteration times")
		display(plt)
    end
end

fitnessPlotController = FitnessPlot([],[[]],0,0,true)

search_range = [
    (0.1, 1.5),   # heatPumpServiceCoff
    (0.0, 10.0),  # heatStorageCapacity  kWh
    (0.0, 10.0)]  # maxheatStorageInputHour h

good_guess = [
	[1.0,0.0,10.0]
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

bb_cost(0.5,1.0,10.0)

println("调用次数：$(getCount())")

best_x = bbo.best_candidate(res)
best_f = bbo.best_fitness(res)

@info "优化结束"
@info "最优变量" heatPumpServiceCoff = best_x[1]heatStorageCapacity = best_x[2]maxheatStorageInputHour = best_x[3]
@info "最小总现值" best_f
