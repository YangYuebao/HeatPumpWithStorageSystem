
# 使用经济性优化需要在项目主目录下切换到calculations环境
#=
用于测试 getMinimumCost 这类函数的计算结果
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

situation = "situation26_只热泵的运行结果"
#第一步，指定设计条件变量
#项目设计条件
begin
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

	heatConsumptionPower = ones(48)

	# 系数
	#heatPumpServiceCoff = 0.5
	maxCOP = 21.0						# 最大COP
	eta_s = 0.7							# 绝热效率
	workingStartHour = 0                # 生产开始时间
	workingHours = 24                   # 每日工作小时数
	TWaste = 30.0                     	# 废热源温度
	#Tair = 25.0                        # 外部环境温度
	TCompressorIn = 115.0
	maxTcHigh = 180.0
	dT_EvaporationStandard = 5.0
	Tsmin = 120.0
	Tsmax = 220.0
	dTRecycleSupply=5.0
    dTRecycleBackward=5.0

	# 计算参数
	dT = 0.01
	#dt = 1/2# 时间步长过小会导致初始温度优化的目标不是一个单峰函数

	dt = 0.5
	smoother = 1e-8
	Tuse = PropsSI("T","P",0.45e6,"Q",0,"water")-273.15
	Tuse = 150.0
	refrigerant = R1233zdE_Water
	sysStruct = RecycleStruct(1,0,0)
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


# 第六步，算运行优化结果
result = optimizeFunction(1.0,0.0,1e8)

include(joinpath(pwd(),"temp","plottool.jl"))

plt = operation_result_plot(
    0:dt:24,
    hourlyTariff,
    result;
    w=0.45
)

#savefig(plt, "plots/plt_2.png")


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
annualOperationCost = result[1] * annualDays
pw,capitalCost,annuity_pv_factor = totalPresentWorth(
    PressedWaterOneStorageOneCompressor(),
	fp,# 经济参数
	result[1],          # 每日运行成本（元）
	1.0,    # 热泵服务系数
	0.0,    # 蓄热容量kwh
	1e9, # 蓄热电加热储满时长h
)

include(joinpath(pwd(),"temp","ParameterDocGenerator.jl"))
generateDesignDoc(situation, designInput, designParameters, fp, joinpath(pwd(),"calculations",situation,situation*"_设计参数.md"))

println("""
日运行费用：$(round(result[1],digits=3))
年运行费用：$(round(annualOperationCost,digits=3))
初投资费用：$(round(capitalCost,digits=3))
总现值：$(round(pw,digits=3))
折现系数：$(round(annuity_pv_factor,digits=3))
""")

plt = operation_result_plot(
    0:dt:24,
    hourlyTariff,
    result;
    w=0.45
)

savefig(plt, joinpath(pwd(),"calculations",situation,"plot.png"))