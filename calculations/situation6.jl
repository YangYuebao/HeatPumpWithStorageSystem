
using Revise
using CSV, DataFrames,Dates
#using Plots
using HeatPumpWithStorageSystem
using LinearAlgebra
using JuMP

# 参数设置

hourlyTariff = zeros(24)
hourlyTariff[1:7] .= 0.3340
hourlyTariff[8:11] .= 0.7393
hourlyTariff[12:14] .= 1.2360
hourlyTariff[15:18] .= 0.7393
hourlyTariff[19:23] .= 1.2360
hourlyTariff[24] = 0.3340
workingStartHour = 0                # 生产开始时间
workingHours = 24                   # 每日工作小时数

Tair_init=fill(25.0, 24)
dTair_init=5.0
dT_l_init=5.0
TWaste_init=60.0
maxTcLow_init=92.0
dTlc_he_init=13.0
maxTeh_init=180.0
maxTcHigh_init=180
refrigerantHigh="water"
refrigerantLow="R134a"
maxCOP=21
eta_s=0.7
COPInterpolateGap=0.1

#COPh, COPh_g, COPh_h, COPl, COPl_g, COPl_h,
COPh, COPh_g, COPh_h, COPl, COPl_g, COPl_h=HeatPumpWithStorageSystem.getCOPFunction(Tair_init,dTair_init,dT_l_init,TWaste_init,maxTcLow_init,dTlc_he_init,maxTeh_init,maxTcHigh_init,
	refrigerantHigh,
	refrigerantLow,
	maxCOP,
	eta_s,
	COPInterpolateGap,
)
maxHeatStorageHour = 19

filePath0=joinpath(pwd(), "calculations", "situation6","WorkAllDay")
if !isdir(filePath0)
	mkdir(filePath0)
end

dfSummary = DataFrame("蓄热时长" => 0:maxHeatStorageHour)
@time for Tuse in 120:10:180
	filePath=joinpath(filePath0,"Temperature$(Tuse)")
	if !isdir(filePath)
        mkdir(filePath)
    end
	dailyCostList = fill(99.0,maxHeatStorageHour+1)

	@Threads.threads for heatStorageTime in 1:maxHeatStorageHour
		(T10, T9, qm, latentHeat, cp_cs,
			cpqm_l, k1, dTe_l1, dTe_l2, dTc_l, QhRecycle, Tair, TWaste,
			cpm_l, KTloss_l, dT_EvaporationStandard, minTeh,
			cpm_h, KTloss_h,
			TstorageTankMax, heatStorageVelocity, heatpumpPowerConstraint,
			hourlyTariffList, heatConsumptionPowerList,
			TcChangeToElec, maxTcLow) = generateSystemCoff(PressedWaterDoubleStorage();
			maxTcHigh = maxTcHigh_init,  # 高温热泵冷凝器温度上限
			maxTcLow = maxTcLow_init,  # 低温热泵冷凝器温度上限
			TWaste = TWaste_init,                    # 废热源温度
			Tair = Tair_init,            # 外部环境温度
			dTlc_he = dTlc_he_init,  # 高温热泵蒸发器温度与低温热泵冷凝器温度差
			Tuse = Tuse+0.0,                       # 工厂使用温度
			Trecycle = Tuse-5.0,    # 回收蒸汽温度
			heatStorageCapacity = heatStorageTime,          # 蓄热量kWh(承压水蓄热)
			TstorageTankMax = 220.0,            # 蓄热罐的最高温度
			maxheatStorageInputHour = 4,        # 蓄热充满时长
			dTstorageInput = 5.0,               # 蓄热温差
			kt = 0.5,                           # 蓄热 \Delta T_2 / \Delta T_1
			dT_l = dT_l_init,    # 低温热泵蒸发器与冷凝器传热温差,也用作低温系统的传热温差
			TwastCapacity = 0.8,                # 废热容量是工厂用热量的倍数
			heatStorageOutEfficiency = 0.000,  # 蓄热衰减系数K
			heatConsumptionPower = 1e0,         # 每小时用热功率kW
			workingStartHour = workingStartHour,                # 生产开始时间
			workingHours = workingHours,                   # 每日工作小时数
			PheatPumpMax = 1.0,                 # 热泵最大功率kW
			hourlyTariff = hourlyTariff,     # 电价向量
		)


		#=
		蓄热的温度不会低于冷凝水的回水温度，所以实际上温度的扩展范围变化不大
		=#

		isFeasible,dailyCost, dfResult = generateAndSolve(PressedWaterDoubleStorage(), MinimizeCost();
			COPh = COPh, COPh_g = COPh_g, COPh_h = COPh_h, COPl = COPl, COPl_g = COPl_g, COPl_h = COPl_h,
			T10 = T10, T9 = T9, qm = qm, latentHeat = latentHeat, cp_cs = cp_cs,
			cpqm_l = cpqm_l, k1 = k1, dTe_l1 = dTe_l1, dTe_l2 = dTe_l2, dTc_l = dTc_l, QhRecycle = QhRecycle, Tair = Tair, TWaste = TWaste,
			cpm_l = cpm_l, KTloss_l = KTloss_l, dT_EvaporationStandard = dT_EvaporationStandard, minTeh = minTeh,
			cpm_h = cpm_h, KTloss_h = KTloss_h,
			TstorageTankMax = TstorageTankMax, heatStorageVelocity = heatStorageVelocity, heatpumpPowerConstraint = heatpumpPowerConstraint,
			hourlyTariffList = hourlyTariffList, heatConsumptionPowerList = heatConsumptionPowerList,
			TcChangeToElec = TcChangeToElec, maxTcLow = maxTcLow,
		)
		if isFeasible in [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]
			CSV.write(joinpath(filePath, "StorageHour_$(heatStorageTime).csv"), round.(dfResult, digits = 4))
			dailyCostList[heatStorageTime+1] = dailyCost
		end
		println("T:$Tuse, Storage:$heatStorageTime, Cost:$(dailyCostList[heatStorageTime+1])")
	end
	dfSummary[!,"$(Tuse)"]=dailyCostList
end

CSV.write(joinpath(filePath0,"summary$(Dates.format(now(), "HH_MM")).csv"), round.(dfSummary,digits=4))


