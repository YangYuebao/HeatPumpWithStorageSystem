using Revise
using Plots, CSV, DataFrames
using HeatPumpWithStorageSystem

begin
Tuse=145.0
TChangeToElec=178.0
heatConsumptionPower = 1.0     # 每小时用热功率kW
heatStorageCapacity = 6.0      # 蓄热量kWh(承压水蓄热)
TwastCapacity = 0.8
heatStorageOutEfficiency = 0.0001
workingStartHour=8
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
TstorageTankMax,
heatStorageOutEfficiency = generateSystemCoff(PressedWaterHighStorage();
    Tuse=Tuse,
    TChangeToElec=TChangeToElec,
    TwastCapacity = TwastCapacity,
    hourlyTariff = hourlyTariff,
    heatConsumptionPower=heatConsumptionPower,
    heatStorageCapacity=heatStorageCapacity,
    heatStorageOutEfficiency=heatStorageOutEfficiency,
    workingStartHour=workingStartHour
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
heatConsumptionPowerList,
costList,
heatStorageList = generateAndSolve(PressedWaterHighStorage(), MinimizeCost();
    COPSupplyWaste=COPSupplyWaste,
	COPSupplyAir=COPSupplyAir,
	COPStorageWaste=COPStorageWaste,
	COPStorageAir=COPStorageAir,
	COPSupplyWaste_g=COPSupplyWaste_g,
	COPSupplyAir_g=COPSupplyAir_g,
	COPStorageWaste_g=COPStorageWaste_g,
	COPStorageAir_g=COPStorageAir_g,
	COPSupplyWaste_h=COPSupplyWaste_h,
	COPSupplyAir_h=COPSupplyAir_h,
	COPStorageWaste_h=COPStorageWaste_h,
	COPStorageAir_h=COPStorageAir_h,
	heatConsumptionPowerList=heatConsumptionPowerList,
	heatStorageCapacityConstraint=heatStorageCapacityConstraint,
	heatpumpPowerConstraint=heatpumpPowerConstraint,
	hourlyTariffList=hourlyTariffList,
	heatStorageVelocity=heatStorageVelocity,
	T4=T4,
	cpqml=cpqml,
	cpqmh=cpqmh,
	cpm=cpm,
	kt=kt,
	TwastCapacity=TwastCapacity,
	TstorageTankMax=TstorageTankMax,
	heatStorageOutEfficiency=1e-3
)

@show isFesable


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
    "用热需求"=>heatConsumptionPowerList,
    "蓄热罐平均温度T"=>TstorageList,
    "蓄热温度T1"=>T1List,
    "蓄热温差"=>T1List-TstorageList,
    "低温热泵出水温度T3"=>T3List,
    "用热温度T4"=>T4List,
    "供热回水温度T5"=>T5List,
    "总成本"=>costList,
)



CSV.write(joinpath(pwd(),"calculations","situation5", "result2.csv"),round.(resultdf,digits=2))

#plt3 = plot(1:25, vcat(TstorageList, TstorageList[1]), title = "2h", xlabel = "Hour", ylabel = "kW", legend = :none)
#plt4 = plot(1:25, vcat(heatStorageList, heatStorageList[1]), title = "2h", xlabel = "Hour", ylabel = "kWh", legend = :none)
