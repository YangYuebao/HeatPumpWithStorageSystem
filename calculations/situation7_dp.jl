using Revise
using HeatPumpWithStorageSystem

COPOverlap, COPWater,
hourlyTariffFunction, heatConsumptionPowerFunction, TairFunction,
Tuse, TCompressorIn,
dT_EvaporationStandard,
latentHeat, cp_cw, cp_cs,
TcChangeToElec, TWaste,
cpm_h,
TstorageTankMax, PheatPumpMax, PelecHeatMax = generateSystemCoff(PressedWaterDoubleStorageOneCompressor();
	overlapRefrigerant = NH3_Water,    # 复叠工质
	maxTcHigh = 180.0,                  # 高温热泵冷凝器温度上限
	TCompressorIn = 115.0,              # 中间温度
	TWaste = 60.0,                      # 废热源温度
	Tair = 25.0,                        # 外部环境温度
	Tuse = 150.0,                       # 工厂使用温度
	Trecycle = 115.0,                   # 回收冷凝水温度
	heatStorageCapacity = 6.0,          # 蓄热量kWh(承压水蓄热)
	TstorageTankMax = 220.0,            # 蓄热罐的最高温度
	maxheatStorageInputHour = 4,        # 蓄热充满时长
	dT_EvaporationStandard = 5.0,       # 全蒸温差
	heatConsumptionPower = 1.0,         # 每小时用热功率kW
	workingStartHour = 0,               # 生产开始时间
	workingHours = 24,                  # 每日工作小时数
	heatPumpServiceCoff = 1.2,          # 热泵功率服务系数
	elecHeatServiceCoff = 1.5,    # 电锅炉服务系数
	hourlyTariff = fill(0.7, 24),       # 电价向量
)

C, TsDecreaseIndexList, TsIncreaseIndexList=generateAndSolve(PressedWaterDoubleStorageOneCompressor(), MinimizeCost(), ConstloadandArea();
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
	dT = 0.1,# 状态参数高温蓄热温度离散步长
	dt = 1 / 6,# 时间步长
)

using Plots
plot(TsDecreaseIndexList)
plot(TsIncreaseIndexList)

using DataFrames,CSV
df = DataFrame(C,:auto)
CSV.write("C.csv", df)
plot(diag(C,10)[1:end-1])

