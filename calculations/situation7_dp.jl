using Revise
using HeatPumpWithStorageSystem

begin
	hourlyTariff = zeros(24)
	hourlyTariff[1:7] .= 0.3340
	hourlyTariff[8:11] .= 0.7393
	hourlyTariff[12:14] .= 1.2360
	hourlyTariff[15:18] .= 0.7393
	hourlyTariff[19:23] .= 1.2360
	hourlyTariff[24] = 0.3340

	or = NH3_Water
	heatStorageCapacity = 3.0
	Tuse = 180.0

	# 系数
	heatPumpServiceCoff = 1.0
	maxCOP = 21# 最大COP
	eta_s = 0.7# 绝热效率
	workingStartHour = 0                # 生产开始时间
	workingHours = 24                   # 每日工作小时数
	TWaste = 55.0                     # 废热源温度
	Tair = 25.0                        # 外部环境温度
	TCompressorIn = 115.0
	maxTcHigh = 180.0
	dT_EvaporationStandard = 5.0
	heatConsumptionPower = 1.0
	# 计算参数
	dT = 0.1
	dt = 1.0# 时间步长过小会导致初始温度优化的目标不是一个单峰函数
	K=2

	COPWater = getCOP(
		TCompressorIn,# 蒸发温度下限,这里是实际设计中的蒸发冷凝温度界限
		maxTcHigh,# 蒸发温度上限
		TCompressorIn,# 冷凝温度下限
		maxTcHigh,# 冷凝温度上限
		refWater,# 工质
		maxCOP,# 最大COP
		eta_s,# 绝热效率
		dT,# 插值步长
	)

	COPOverlap = getOverlapCOP_fixMidTemperature(
		or,
		TCompressorIn + or.midTDifference / 2;
		maxCOP = maxCOP,# 最大COP
		eta_s = eta_s,# 绝热效率
		dT = dT,# 插值步长
	)	
end

#=
			GoldenRatioMethod	
步长1/1	cost 6.040405104886531
步长1/6 cost 6.211648136391782
步长1/10 cost 6.13
步长1/15 cost 6.359773767716163
=#

hourlyTariffFunction, heatConsumptionPowerFunction, TairFunction,
Tuse, TCompressorIn,
dT_EvaporationStandard,
latentHeat, cp_cw, cp_cs,
TcChangeToElec, TWaste,
cpm_h,
TstorageTankMax, PheatPumpMax, PelecHeatMax = generateSystemCoff(PressedWaterDoubleStorageOneCompressor();
	overlapRefrigerant = or,    # 复叠工质
	maxTcHigh = maxTcHigh,                  # 高温热泵冷凝器温度上限
	TCompressorIn = TCompressorIn,              # 中间温度
	TWaste = TWaste,                      # 废热源温度
	Tair = Tair,                        # 外部环境温度
	Tuse = Tuse,                       # 工厂使用温度
	heatStorageCapacity = heatStorageCapacity,  # 蓄热量kWh(承压水蓄热)
	TstorageTankMax = 220.0,            # 蓄热罐的最高温度
	maxheatStorageInputHour = 4,        # 蓄热充满时长
	dT_EvaporationStandard = dT_EvaporationStandard,       # 全蒸温差
	heatConsumptionPower = heatConsumptionPower,         # 每小时用热功率kW
	workingStartHour = workingStartHour,# 生产开始时间
	workingHours = workingHours,                  # 每日工作小时数
	heatPumpServiceCoff = heatPumpServiceCoff,    # 热泵功率服务系数
	hourlyTariff = hourlyTariff,       # 电价向量
)

#=
@time bestValueList,TsMatrix,minCostTest, minTsListTest, P1ListTest, P2ListTest, P3ListTest, PeListTest = HeatPumpWithStorageSystem.testSolve(PressedWaterDoubleStorageOneCompressor(), MinimizeCost(), ConstloadandArea(), ExhaustiveMethod();
	COPOverlap = COPOverlap,
	COPWater = COPWater,
	hourlyTariffFunction = hourlyTariffFunction,
	heatConsumptionPowerFunction = heatConsumptionPowerFunction,
	TairFunction = TairFunction,
	Tuse = Tuse,
	TCompressorIn = TCompressorIn,
	dT_EvaporationStandard = dT_EvaporationStandard,
	latentHeat = latentHeat,
	cp_cw = cp_cw,
	cp_cs = cp_cs,
	TcChangeToElec = TcChangeToElec,
	TWaste = TWaste,
	cpm_h = cpm_h,
	TstorageTankMax = TstorageTankMax,
	PheatPumpMax = PheatPumpMax,
	PelecHeatMax = PelecHeatMax,
	# 求解参数
	dT = dT,# 状态参数高温蓄热温度离散步长
	dt = dt,# 时间步长
	K=K
)

=#
# 2025年3月17日07:14:52
@time minCostTest, minTsListTest, P1ListTest, P2ListTest, P3ListTest, PeListTest = generateAndSolve(PressedWaterDoubleStorageOneCompressor(), MinimizeCost(), ConstloadandArea(), GoldenRatioMethod();
	COPOverlap = COPOverlap,
	COPWater = COPWater,
	hourlyTariffFunction = hourlyTariffFunction,
	heatConsumptionPowerFunction = heatConsumptionPowerFunction,
	TairFunction = TairFunction,
	Tuse = Tuse,
	TCompressorIn = TCompressorIn,
	dT_EvaporationStandard = dT_EvaporationStandard,
	latentHeat = latentHeat,
	cp_cw = cp_cw,
	cp_cs = cp_cs,
	TcChangeToElec = TcChangeToElec,
	TWaste = TWaste,
	cpm_h = cpm_h,
	TstorageTankMax = TstorageTankMax,
	PheatPumpMax = PheatPumpMax,
	PelecHeatMax = PelecHeatMax,
	# 求解参数
	dT = dT,# 状态参数高温蓄热温度离散步长
	dt = dt,# 时间步长
	K=K
)
#= K=2 dt=1.0
170 3.0 6.68601222715751
170 3.75 6.341932166742927
180 3.0 6.692398996519907
=#
using Plots, DataFrames, CSV

tList = 0:dt:(24-dt)

# 看看不同的初始温度下运行费用的变化

# 看看最优运行费用下蓄热温度的变化
plt = plot(tList .+ 0.5, [P1ListTest P2ListTest P3ListTest PeListTest], label = ["P1" "P2" "P3" "Pe"])
plot!(plt, tList, minTsListTest[1:end-1] / 220, label = "Ts")

#=
dt	K	Tuse	Storage		Cost		Ts0
1/3 12	170		3.0			6.68052		132.18
1.0 12	170		3.0			6.68007		137.83
1/2	6	170		3.0			6.68341		140.5
1.0	1	170		3.0			6.69800		130.3
1/2	2	170		3.0			6.70000		139.35
1/2 8	170		3.0			6.68623		134.68
=#

