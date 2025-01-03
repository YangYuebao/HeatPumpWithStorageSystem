
using Revise
using CSV, DataFrames
#using Plots
using HeatPumpWithStorageSystem
(COPh, COPh_g, COPh_h,
	COPl, COPl_g, COPl_h,
	T13, T9, qm, latenHeat, cp_cw, cp_cs,
	cpqm_l, k1, dTe_l1, dTe_l2, dTc_l, QhRecycle, Tair, TWaste,
	cpm_l, KTloss_l, k6,
	cpqm_m, k2, dTe_h, k3, dTc_h1, dTc_h2, cpqm_h, minTeh,
	cpm_h, KTloss_h, k4, k5, k7, k8,
	TstorageTankMax, heatStorageVelocity, heatStorageCapacityConstraint, heatpumpPowerConstraint,
	hourlyTariffList, heatConsumptionPowerList) = generateSystemCoff(PressedWaterDoubleStorage();
	refrigerantLow = "R134a",    # 供热循环使用制冷剂
	refrigerantHigh = "water",   # 蓄热使用制冷剂
	maxTcHigh = 180.0,  # 高温热泵冷凝器温度上限
	maxTcLow = 80.0,  # 低温热泵冷凝器温度上限
	eta_s = 0.7,                        # 压缩机绝热效率
	Twastein = 80.0,                    # 废气进入温度
	Twasteout = 40.0,                   # 废气排出温度
	Tair = fill(25.0, 24),            # 外部环境温度
	dTair = 5.0,                        # 外部环境温度-蒸发器温度
	dTlc_he = 10.0,  # 高温热泵蒸发器温度与低温热泵冷凝器温度差
	maxTeh = 180.0,  # 高温热泵蒸发器温度上限
	Tuse = 120.0,                       # 工厂使用温度
	Trecycle = 115.0,    # 回收蒸汽温度
	heatStorageCapacity = 6.0,          # 蓄热量kWh(承压水蓄热)
	TstorageTankMax = 220.0,            # 蓄热罐的最高温度
	maxheatStorageInputHour = 4,        # 蓄热充满时长
	dTstorageInput = 5.0,               # 蓄热温差
	kt = 0.5,                           # 蓄热 \Delta T_2 / \Delta T_1
	dT_l = 5.0,    # 低温热泵蒸发器与冷凝器传热温差,也用作低温系统的传热温差
	TwastCapacity = 0.8,                # 废热容量是工厂用热量的倍数
	heatStorageOutEfficiency = 0.0001,  # 蓄热衰减系数K
	dT_h = 5.0,    # 高温热泵蒸发器与冷凝器传热温差
	heatConsumptionPower = 1.0,         # 每小时用热功率kW
	maxCOP = 21.0,                      # 热泵COP上限
	workingStartHour = 8,                # 生产开始时间
	workingHours = 16,                   # 每日工作小时数
	PheatPumpMax = 1.0,                 # 热泵最大功率kW
	hourlyTariff = fill(0.7, 24),     # 电价向量
	COPInterpolateGap = 0.1,    # COP插值时步长
	Te_hStandard =80.0
)

generateAndSolve(PressedWaterDoubleStorage(), MinimizeCost();
	COPSupplyWaste = COPSupplyWaste,
	COPSupplyAir = COPSupplyAir,
	COPStorageWaste = COPStorageWaste,
	COPStorageAir = COPStorageAir,
	COPSupplyWaste_g = COPSupplyWaste_g,
	COPSupplyAir_g = COPSupplyAir_g,
	COPStorageWaste_g = COPStorageWaste_g,
	COPStorageAir_g = COPStorageAir_g,
	COPSupplyWaste_h = COPSupplyWaste_h,
	COPSupplyAir_h = COPSupplyAir_h,
	COPStorageWaste_h = COPStorageWaste_h,
	COPStorageAir_h = COPStorageAir_h,
	heatConsumptionPowerList = heatConsumptionPowerList,
	heatStorageCapacityConstraint = heatStorageCapacityConstraint,
	heatpumpPowerConstraint = heatpumpPowerConstraint,
	hourlyTariffList = hourlyTariffList,
	heatStorageVelocity = heatStorageVelocity,
	T4 = T4,
	cpqml = cpqml,
	cpqmh = cpqmh,
	cpm = cpm,
	kt = kt,
	TwastCapacity = TwastCapacity,
	TstorageTankMax = TstorageTankMax,
	heatStorageOutEfficiency = 1e-3,
)

@show isFesable


resultdf = DataFrame(
	"时间" => 0:23,
	"低温热泵废热功率Pl1" => Pl1List,
	"低温热泵空气源功率Pl2" => Pl2List,
	"高温热泵废热功率Ph1" => Ph1List,
	"高温热泵空气源功率Ph2" => Ph2List,
	"低温热泵废热COP" => COPl1,
	"低温热泵空气源COP" => COPl2,
	"高温热泵废热COP" => COPh1,
	"高温热泵空气源COP" => COPh2,
	"热泵总功率" => PList,
	"蓄热量" => heatStorageList,
	"用热需求" => heatConsumptionPowerList,
	"蓄热罐平均温度T" => TstorageList,
	"蓄热温度T1" => T1List,
	"蓄热温差" => T1List - TstorageList,
	"低温热泵出水温度T3" => T3List,
	"用热温度T4" => T4List,
	"供热回水温度T5" => T5List,
	"总成本" => costList,
)



CSV.write(joinpath(pwd(), "calculations", "situation5", "result2.csv"), round.(resultdf, digits = 2))

#plt3 = plot(1:25, vcat(TstorageList, TstorageList[1]), title = "2h", xlabel = "Hour", ylabel = "kW", legend = :none)
#plt4 = plot(1:25, vcat(heatStorageList, heatStorageList[1]), title = "2h", xlabel = "Hour", ylabel = "kWh", legend = :none)
