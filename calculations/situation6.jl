
using Revise
using CSV, DataFrames
#using Plots
using HeatPumpWithStorageSystem
using LinearAlgebra

# 参数设置
begin
	hourlyTariff = zeros(24)
	hourlyTariff[1:6] .= 0.3340
	hourlyTariff[7:10] .= 0.7393
	hourlyTariff[11:13] .= 1.2360
	hourlyTariff[14:17] .= 0.7393
	hourlyTariff[18:22] .= 1.2360
	hourlyTariff[23:24] .= 0.3340

	(COPh, COPh_g, COPh_h,
		COPl, COPl_g, COPl_h,
		T13, T9, qm, latenHeat, cp_cw, cp_cs,
		cpqm_l, k1, dTe_l1, dTe_l2, dTc_l, QhRecycle, Tair, TWaste,
		cpm_l, KTloss_l, k6,
		cpqm_m, k2, dTe_h, k3, dTc_h1, dTc_h2, cpqm_h, minTeh,
		cpm_h, KTloss_h, k4, k5, k7, k8,
		TstorageTankMax, heatStorageVelocity, heatStorageCapacityConstraint, heatpumpPowerConstraint,
		hourlyTariffList, heatConsumptionPowerList,
		TcChangeToElec,maxTcLow) = generateSystemCoff(PressedWaterDoubleStorage();
		refrigerantLow = "R134a",    # 供热循环使用制冷剂
		refrigerantHigh = "water",   # 蓄热使用制冷剂
		maxTcHigh = 180.0,  # 高温热泵冷凝器温度上限
		maxTcLow = 92.0,  # 低温热泵冷凝器温度上限
		eta_s = 0.7,                        # 压缩机绝热效率
		Twastein = 80.0,                    # 废气进入温度
		Twasteout = 40.0,                   # 废气排出温度
		Tair = fill(25.0, 24),            	# 外部环境温度
		dTair = 5.0,                        # 外部环境温度-蒸发器温度
		dTlc_he = 13.0,  					# 高温热泵蒸发器温度与低温热泵冷凝器温度差
		maxTeh = 180.0,  					# 高温热泵蒸发器温度上限
		Tuse = 120.0,                       # 工厂使用温度
		Trecycle = 115.0,    # 回收蒸汽温度
		heatStorageCapacity = 6.0,          # 蓄热量kWh(承压水蓄热)
		TstorageTankMax = 220.0,            # 蓄热罐的最高温度
		maxheatStorageInputHour = 4,        # 蓄热充满时长
		dTstorageInput = 5.0,               # 蓄热温差
		kt = 0.5,                           # 蓄热 \Delta T_2 / \Delta T_1
		dT_l = 5.0,    # 低温热泵蒸发器与冷凝器传热温差,也用作低温系统的传热温差
		TwastCapacity = 0.8,                # 废热容量是工厂用热量的倍数
		heatStorageOutEfficiency = 0.000,  # 蓄热衰减系数K
		dT_h = 5.0,    # 高温热泵蒸发器与冷凝器传热温差
		heatConsumptionPower = 1e0,         # 每小时用热功率kW
		maxCOP = 21.0,                      # 热泵COP上限
		workingStartHour = 8,                # 生产开始时间
		workingHours = 16,                   # 每日工作小时数
		PheatPumpMax = 1.0,                 # 热泵最大功率kW
		hourlyTariff = hourlyTariff,     # 电价向量
		COPInterpolateGap = 0.1,    # COP插值时步长
	)
end

#=
蓄热的温度不会低于冷凝水的回水温度，所以实际上温度的扩展范围变化不大
=#

generateAndSolve(PressedWaterDoubleStorage(), MinimizeCost();
	COPh = COPh, COPh_g = COPh_g, COPh_h = COPh_h,
	COPl = COPl, COPl_g = COPl_g, COPl_h = COPl_h,
	T13 = T13, T9 = T9, qm = qm, latenHeat = latenHeat, cp_cw = cp_cw, cp_cs = cp_cs,
	cpqm_l = cpqm_l, k1 = k1, dTe_l1 = dTe_l1, dTe_l2 = dTe_l2, dTc_l = dTc_l, QhRecycle = QhRecycle, Tair = Tair, TWaste = TWaste,
	cpm_l = cpm_l, KTloss_l = KTloss_l, k6 = k6,
	cpqm_m = cpqm_m, k2 = k2, dTe_h = dTe_h, k3 = k3, dTc_h1 = dTc_h1, dTc_h2 = dTc_h2, cpqm_h = cpqm_h, minTeh = minTeh,
	cpm_h = cpm_h, KTloss_h = KTloss_h, k4 = k4, k5 = k5, k7 = k7, k8 = k8,
	TstorageTankMax = TstorageTankMax, heatStorageVelocity = heatStorageVelocity, heatStorageCapacityConstraint = heatStorageCapacityConstraint, heatpumpPowerConstraint = heatpumpPowerConstraint,
	hourlyTariffList = hourlyTariffList, heatConsumptionPowerList = heatConsumptionPowerList,
	TcChangeToElec = TcChangeToElec,
	maxTcLow=maxTcLow,
)

#3.901206566133562

