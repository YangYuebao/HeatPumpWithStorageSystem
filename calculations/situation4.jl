using Revise
#using Plots, CSV, DataFrames
using HeatPumpWithStorageSystem

begin
hourlyTariff = zeros(24)
hourlyTariff[1:6] .= 0.3340
hourlyTariff[7:10] .= 0.7393
hourlyTariff[11:13] .= 1.2360
hourlyTariff[14:17] .= 0.7393
hourlyTariff[18:22] .= 1.2360
hourlyTariff[23:24] .= 0.3340

COPSupplyWaste, COPSupplyAir, COPStorageWaste, COPStorageAir,
COPSupplyWaste_g,COPSupplyAir_g,COPStorageWaste_g,COPStorageAir_g,
COPSupplyWaste_h,COPSupplyAir_h,COPStorageWaste_h,COPStorageAir_h,
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
TstorageTankMax = generateSystemCoff(HeatPumpStoragePressedWater())
end

isFesable,
Pl1List,
Pl2List,
Ph1List,
Ph2List,
PList,
TstorageList,
costList,
heatStorageList=generateAndSolve(HeatPumpStoragePressedWater(), MinimizeCost();
    COPSupplyWaste, COPSupplyAir, COPStorageWaste, COPStorageAir,
    COPSupplyWaste_g,COPSupplyAir_g,COPStorageWaste_g,COPStorageAir_g,
    COPSupplyWaste_h,COPSupplyAir_h,COPStorageWaste_h,COPStorageAir_h,
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
    TstorageTankMax
)



using Plots
plt3 = plot(1:25, vcat(TstorageList, TstorageList[1]), title="2h", xlabel="Hour", ylabel="kW", legend=:none)
plt4 = plot(1:25, vcat(heatStorageList, heatStorageList[1]), title="2h", xlabel="Hour", ylabel="kWh", legend=:none)