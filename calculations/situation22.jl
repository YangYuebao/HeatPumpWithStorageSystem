using Revise
using Pkg
using HeatPumpWithStorageSystem
using DataFrames, CSV

#=
封装了用于设计优化的函数
=#

situation = "situation22"
#第一步，指定设计条件变量
begin
	hourlyTariff = ones(48) * 0.7393
	p = 1.7
	pp = p * 1.2
	v = 0.35
	vv = v * 0.8

	hourlyTariff[1:12] .*= v
	hourlyTariff[23:26] .*= v
	hourlyTariff[29:30] .*= pp
	hourlyTariff[31:39] .*= p
	hourlyTariff[40:43] .*= pp
	hourlyTariff[44] *= p

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
	heatPumpServiceCoff = 0.5
	maxCOP = 21.0# 最大COP
	eta_s = 0.7# 绝热效率
	workingStartHour = 0                # 生产开始时间
	workingHours = 24                   # 每日工作小时数
	TWaste = 30.0                     # 废热源温度
	#Tair = 25.0                        # 外部环境温度
	TCompressorIn = 115.0
	maxTcHigh = 180.0
	dT_EvaporationStandard = 5.0
	Tsmin = 120.0
	Tsmax = 220.0
	dTRecycleSupply=5.0
    dTRecycleBackward=5.0

	# 计算参数
	dT = 0.05
	#dt = 1/2# 时间步长过小会导致初始温度优化的目标不是一个单峰函数

	dt = 0.5
	smoother = 1e-8
	Tuse = 150.0
	refrigerant = NH3_Water
	sysStruct = RecycleStruct(1,0,0)
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
optimizeFunction = generateOperationFunction(designParameters,designInput)


# 第五步，现在可以调用函数进行优化了
heatPumpServiceCoff = 0.5
heatStorageCapacity = 8.0
maxheatStorageInputHour=4.0
result = optimizeFunction(
	heatPumpServiceCoff,    # 热泵服务系数
	heatStorageCapacity,    # 蓄热容量
	maxheatStorageInputHour    # 蓄热电加热储满时长
)


