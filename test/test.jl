using Revise
using Plots, CSV, DataFrames
using HeatPumpWithStorageSystem

begin

heatConsumptionPower = 1.0     # 每小时用热功率kW
heatStorageCapacity = 6.0      # 蓄热量kWh(承压水蓄热)
TwastCapacity = 0.8
hourlyTariff = zeros(24)
hourlyTariff[1:6] .= 0.3340
hourlyTariff[7:10] .= 0.7393
hourlyTariff[11:13] .= 1.2360
hourlyTariff[14:17] .= 0.7393
hourlyTariff[18:22] .= 1.2360
hourlyTariff[23:24] .= 0.3340



COPSupplyWaste, COPSupplyAir, COPStorageWaste, COPStorageAir,
COPSupplyWaste_g, COPSupplyAir_g, COPStorageWaste_g, COPStorageAir_g,
COPSupplyWaste_h, COPSupplyAir_h, COPStorageWaste_h, COPStorageAir_h,
heatConsumptionPowerList,
heatStorageCapacityConstraint,
heatpumpPowerConstraint,
hourlyTariffList,
heatStorageVelocity,
T4,
cpqml,
cpqmh,
cpm,
kt,
TwastCapacity,
TstorageTankMax = generateSystemCoff(PressedWaterHighStorage();
    TwastCapacity = TwastCapacity,
    hourlyTariff = hourlyTariff,
    heatConsumptionPower=heatConsumptionPower,
    heatStorageCapacity=heatStorageCapacity
)
end

isFesable,
Pl1List,
Pl2List,
Ph1List,
Ph2List,
PList,
TstorageList,
T1List,
T3List,
T4List,
T5List,
COPl1,COPl2,COPh1,COPh2,
costList,
heatStorageList = generateAndSolve(PressedWaterHighStorage(), MinimizeCost();
	COPSupplyWaste, COPSupplyAir, COPStorageWaste, COPStorageAir,
	COPSupplyWaste_g, COPSupplyAir_g, COPStorageWaste_g, COPStorageAir_g,
	COPSupplyWaste_h, COPSupplyAir_h, COPStorageWaste_h, COPStorageAir_h,
	heatConsumptionPowerList,  # 用热负载向量
	heatStorageCapacityConstraint, # 蓄热量约束（最大值）
	heatpumpPowerConstraint,   # 热泵功率约束（最大值）
	hourlyTariffList,   # 电价向量
	heatStorageVelocity,           # 蓄热速率约束
	T4,
	cpqml,
	cpqmh,
	cpm,
	kt,
	TwastCapacity,
	TstorageTankMax,
)


resultdf=DataFrame(
    "时间"=>0:23,
    "低温热泵废热功率Pl1"=>Pl1List,
    "低温热泵空气源功率Pl2"=>Pl2List,
    "高温热泵废热功率Ph1"=>Ph1List,
    "高温热泵空气源功率Ph2"=>Ph2List,
    "低温热泵废热COP"=>COPl1,
    "低温热泵空气源COP"=>COPl2,
    "高温热泵废热COP"=>COPh1,
    "高温热泵空气源COP"=>COPh2,
    "热泵总功率"=>PList,
    "蓄热量"=>heatStorageList,
    "蓄热罐平均温度T"=>TstorageList,
    "蓄热温度T1"=>T1List,
    "低温热泵出水温度T3"=>T3List,
    "用热温度T4"=>T4List,
    "供热回水温度T5"=>T5List,
    "总成本"=>costList,
)

using CSV,DataFrames