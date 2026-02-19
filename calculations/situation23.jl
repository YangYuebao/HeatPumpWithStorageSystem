
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
	hourlyTariff = zeros(24)
	hourlyTariff[1:8] .= 1.0094
	hourlyTariff[9:16] .= 0.658
	hourlyTariff[17:24] .= 0.3725

	Tair = vcat(
		fill(26.0, 7),
		fill(27.0, 7),
		fill(26.0, 11),
	)

	heatConsumptionPower = ones(24)

	# 系数
	#heatPumpServiceCoff = 0.5
	maxCOP = 21.0						# 最大COP
	eta_s = 0.7							# 绝热效率
	workingStartHour = 0                # 生产开始时间
	workingHours = 24                   # 每日工作小时数
	TWaste = 85.0                     	# 废热源温度
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

	dt = 1.0
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
	lifeYears = 6
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


# 第五步，现在可以调用函数进行优化了

# 1.7594



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


heatLoad=1.0
heatPumpServiceCoff = 1.0
heatStorageCapacity = 0.5
maxheatStorageInputHour=10.0

TsStart=175.0
TsEnd=170.0
dt=1.0

params = getParams(heatPumpServiceCoff,heatStorageCapacity,maxheatStorageInputHour)
sysVariables = SystemVariables(
	heatLoad,
	designParameters.COPLowFunction(TWaste, TCompressorIn + dT_EvaporationStandard),
	26.0,
	TWaste,
)

#=
[ Info: 内部检查结束
(0.24808180185848658, true, 0.24634801692900746, 0.0017337849294791472, 0.0, 0.0)
=#
single_result = getMinimumCost(TsStart,TsEnd,dt,params,sysVariables;show=true)

TsList2 = vcat(
	175.0,
	fill(170.0,20),
	fill(175.0,4)
)
TsList1 = fill(120.0,25)
TsList3 = vcat(
	120.0:100/6:220.0,#7个
	fill(220.0,11),
	220.0:-100/6:120.0
)

line_result = getTemperatureLineCost(TsList3;
	COPLowFunction=designParameters.COPLowFunction,
	hourlyTariffFunction=designParameters.hourlyTariffFunction,   # 电价函数
	heatConsumptionPowerFunction=designParameters.heatConsumptionPowerFunction,  # 用热负载函数
	TairFunction=designParameters.TairFunction,# 环境温度函数

	# 总循环参数
	params=params,
	TWaste=TWaste,# 废热回收蒸发器温度

	# 求解参数
	dt = dt,# 时间步长
)

result = optimizeFunction(
	heatPumpServiceCoff,    # 热泵服务系数
	heatStorageCapacity,    # 蓄热容量
	maxheatStorageInputHour    # 蓄热电加热储满时长
)

getTemperatureLineCost(result[2];
	COPLowFunction=designParameters.COPLowFunction,
	hourlyTariffFunction=designParameters.hourlyTariffFunction,   # 电价函数
	heatConsumptionPowerFunction=designParameters.heatConsumptionPowerFunction,  # 用热负载函数
	TairFunction=designParameters.TairFunction,# 环境温度函数

	# 总循环参数
	params=params,
	TWaste=TWaste,# 废热回收蒸发器温度

	# 求解参数
	dt = dt,# 时间步长
)

plot(result[2])#Ts
plot([result[3],result[4],result[5],result[6]])#P1
plot(result[4])#P2
plot(result[5])#P3
plot(result[6])
vscodedisplay(result[2])