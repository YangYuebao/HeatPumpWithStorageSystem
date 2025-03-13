using Revise
using HeatPumpWithStorageSystem
begin
begin
	hourlyTariff = zeros(24)
	hourlyTariff[1:7] .= 0.3340
	hourlyTariff[8:11] .= 0.7393
	hourlyTariff[12:14] .= 1.2360
	hourlyTariff[15:18] .= 0.7393
	hourlyTariff[19:23] .= 1.2360
	hourlyTariff[24] = 0.3340
	workingStartHour = 0                # 生产开始时间
	workingHours = 24                   # 每日工作小时数
	
	heatStorageCapacity=6.0
	heatPumpServiceCoff=1.2
	elecHeatServiceCoff=1.5
	dT=0.1
	dt=3
end

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
	heatStorageCapacity = heatStorageCapacity,          # 蓄热量kWh(承压水蓄热)
	TstorageTankMax = 220.0,            # 蓄热罐的最高温度
	maxheatStorageInputHour = 4,        # 蓄热充满时长
	dT_EvaporationStandard = 5.0,       # 全蒸温差
	heatConsumptionPower = 1.0,         # 每小时用热功率kW
	workingStartHour = workingStartHour,               # 生产开始时间
	workingHours = workingHours,                  # 每日工作小时数
	heatPumpServiceCoff = heatPumpServiceCoff,          # 热泵功率服务系数
	elecHeatServiceCoff = elecHeatServiceCoff,    # 电锅炉服务系数
	hourlyTariff = hourlyTariff,       # 电价向量
)

@time minCostEx,minTsListEx,P1ListEx,P2ListEx,P3ListEx,PeListEx=generateAndSolve(PressedWaterDoubleStorageOneCompressor(), MinimizeCost(), ConstloadandArea(),ExhaustiveMethod();
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
)

@time minCostGo,minTsListGo,P1ListGo,P2ListGo,P3ListGo,PeListGo=generateAndSolve(PressedWaterDoubleStorageOneCompressor(), MinimizeCost(), ConstloadandArea(),GoldenRatioMethod();
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
)

println(minTsListEx-minTsListGo)
println(minCostEx-minCostGo)
end
using Plots,DataFrames,CSV

TsList = TCompressorIn+dT_EvaporationStandard:dT:TstorageTankMax
tList=0:dt:24

# 看看不同的初始温度下运行费用的变化

# 看看最优运行费用下蓄热温度的变化
plot(tList,minTsList,title="Ts",legend=:none)
plot(tList[1:end-1].+0.5,P1List,title="P1",legend=:none)
plot(tList[1:end-1].+0.5,P2List,title="P2",legend=:none)
plot(tList[1:end-1].+0.5,P3List,title="P3",legend=:none)
plot(tList[1:end-1].+0.5,PeList,title="Pe",legend=:none)

plt=plot(tList[1:end-1].+0.5,[P1List P2List P3List PeList],label=["P1" "P2" "P3" "Pe"])
plot!(plt,tList,minTsList/220,label="Ts")

P1List,P2List,P3List,PeList


gridPrice=[hourlyTariffFunction(t) for t in tList]
