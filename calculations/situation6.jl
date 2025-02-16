begin
	using Revise
	using CSV, DataFrames, Dates
	#using Plots
	using HeatPumpWithStorageSystem
	using LinearAlgebra
	using JuMP
end

# 参数设置
TuseList=120:10:180
function main(TuseList)
#begin
	hourlyTariff = zeros(24)
	hourlyTariff[1:7] .= 0.3340
	hourlyTariff[8:11] .= 0.7393
	hourlyTariff[12:14] .= 1.2360
	hourlyTariff[15:18] .= 0.7393
	hourlyTariff[19:23] .= 1.2360
	hourlyTariff[24] = 0.3340
	workingStartHour = 0                # 生产开始时间
	workingHours = 24                   # 每日工作小时数

	Tair_init = fill(25.0, 24)
	dTair_init = 5.0
	dT_l_init = 5.0
	TWaste_init = 60.0
	maxTcLow_init = 92.0
	dTlc_he_init = 13.0
	maxTeh_init = 180.0
	maxTcHigh_init = 180
	refrigerantHigh = "water"
	refrigerantLow = "R134a"
	maxCOP = 21
	eta_s = 0.7
	COPInterpolateGap = 0.1

	#COPh, COPh_g, COPh_h, COPl, COPl_g, COPl_h,
	COPh, COPh_g, COPh_h, COPl, COPl_g, COPl_h = getCOPFunction(Tair_init, dTair_init, dT_l_init, TWaste_init, maxTcLow_init, dTlc_he_init, maxTeh_init, maxTcHigh_init,
		refrigerantHigh,
		refrigerantLow,
		maxCOP,
		eta_s,
		COPInterpolateGap,
	)
	maxHeatStorageHour = 19

	filePath0 = joinpath(pwd(), "calculations", "situation6", "WorkAllDay")
	filePath2 = joinpath(filePath0,"Summary_24Hour")
	if !isdir(filePath0)
		mkdir(filePath0)
	end

	dfSummary = DataFrame("蓄热时长" => 0:maxHeatStorageHour)

	m = 24# 24小时
	#@Threads.threads 
	for Tuse in TuseList
		filePath = joinpath(filePath0, "Temperature$(Tuse)")
		
		if !isdir(filePath)
			mkdir(filePath)
		end
		if !isdir(filePath2)
			mkdir(filePath2)
		end
		dailyCostList = fill(99.0, maxHeatStorageHour + 1)

		# 初始化初值信息
		reinitialize = true
		qm1_init = qm2_init = qm3_init = T3_init = T8_init = Tc_l_init = P_l1_init = P_l2_init = P_h1_init = P_h2_init = P_h3_init = Qc_l_init = Qs_l_init = Qs_h_init = COPh1_temp = COPh2_temp = COPh3_temp = zeros(24)
		lastHeatStorageTime = 0.0
		objVal=99.0

		dataEconomy=DataFrame("工作时长"=>[],"蓄热时长"=>[],"承压水体积"=>[],"低温热泵总功率"=>[],"高温热泵总功率"=>[],"每天运行费用"=>[])
		dataTech=DataFrame("工作时长"=>[],"蓄热时长"=>[],"承压水体积"=>[],"低温热泵总功率"=>[],"高温热泵总功率"=>[],"每天运行费用"=>[],"低温热泵1功率"=>[],"低温热泵1最低蒸发温度"=>[],"低温热泵1最高冷凝温度"=>[],"低温热泵2功率"=>[],"低温热泵2最低蒸发温度"=>[],"低温热泵2最高冷凝温度"=>[],"高温热泵1功率"=>[],"高温热泵1最低蒸发温度"=>[],"高温热泵1最高冷凝温度"=>[],"高温热泵2功率"=>[],"高温热泵2最低蒸发温度"=>[],"高温热泵2最高冷凝温度"=>[],"高温热泵3功率"=>[],"高温热泵3最低蒸发温度"=>[],"高温热泵3最高冷凝温度"=>[])
		for heatStorageTime in 0:maxHeatStorageHour
			(T10, T9, qm, latentHeat, cp_cs,
				dTe_l1, dTe_l2, dTc_l, QhRecycle, Tair, TWaste,
				cpm_l, KTloss_l, dT_EvaporationStandard, minTeh,
				cpm_h, KTloss_h,
				TstorageTankMax, heatStorageVelocity, heatpumpPowerConstraint,
				hourlyTariffList, heatConsumptionPowerList,
				maxTcLow,
				storageTankVolume) = generateSystemCoff(PressedWaterDoubleStorage();
				maxTcHigh = maxTcHigh_init,  # 高温热泵冷凝器温度上限
				maxTcLow = maxTcLow_init,  # 低温热泵冷凝器温度上限
				TWaste = TWaste_init,                    # 废热源温度
				Tair = Tair_init,            # 外部环境温度
				dTlc_he = dTlc_he_init,  # 高温热泵蒸发器温度与低温热泵冷凝器温度差
				Tuse = Tuse + 0.0,                       # 工厂使用温度
				Trecycle = Tuse - 5.0,    # 回收蒸汽温度
				heatStorageCapacity = heatStorageTime,          # 蓄热量kWh(承压水蓄热)
				TstorageTankMax = 220.0,            # 蓄热罐的最高温度
				maxheatStorageInputHour = 4,        # 蓄热充满时长
				dT_EvaporationStandard = 3.0,
				dT_l = dT_l_init,    # 低温热泵蒸发器与冷凝器传热温差,也用作低温系统的传热温差
				TwastCapacity = 0.8,                # 废热容量是工厂用热量的倍数
				heatStorageOutEfficiency = 0.001,  # 蓄热衰减系数K
				heatConsumptionPower = 1e0,         # 每小时用热功率kW
				workingStartHour = workingStartHour,                # 生产开始时间
				workingHours = workingHours,                   # 每日工作小时数
				PheatPumpMax = 1.0,                 # 热泵最大功率kW
				hourlyTariff = hourlyTariff,     # 电价向量
			)

			# 以上一次的优化结果为基础，调节蓄热罐容量后生成新的初值，并重新计算热泵功率
			if reinitialize == false
				Tair_min = minimum(Tair)# 环境最低温度
				Te_l1 = TWaste - dTe_l1# 低温热泵废热源蒸发器温度
				Te_l2 = Tair .- dTe_l2# 低温热泵空气源蒸发器温度

				# 计算初值
				maxT3 = maxTcLow - dTc_l
				minT3 = minTeh + dT_EvaporationStandard
				# 更新初值的时候主要修改T8和5个功率
				T8_init = (T8_init .- T10) * lastHeatStorageTime / heatStorageTime .+ T10
				# 计算5个COP

				COPh1_init = map(i -> COPh(T3_init[i] - dT_EvaporationStandard, T10), 1:m)
				COPh2_init = map(i -> COPh(T8_init[i] - dT_EvaporationStandard, T10), 1:m)
				COPh3_init = map(i -> COPh(T3_init[i] - dT_EvaporationStandard, T8_init[i] + dT_EvaporationStandard), 1:m)

				P_h1_init = P_h1_init .* COPh1_temp ./ COPh1_init
				P_h2_init = P_h2_init .* COPh2_temp ./ COPh2_init
				P_h3_init = P_h3_init .* COPh3_temp ./ COPh3_init

			end

			# 从上次的最优解出发
			isFeasible, dailyCost, dfResult = generateAndSolve(PressedWaterDoubleStorage(), MinimizeCost();
				COPh = COPh, COPh_g = COPh_g, COPh_h = COPh_h, COPl = COPl, COPl_g = COPl_g, COPl_h = COPl_h,
				T10 = T10, T9 = T9, qm = qm, latentHeat = latentHeat, cp_cs = cp_cs,
				dTe_l1 = dTe_l1, dTe_l2 = dTe_l2, dTc_l = dTc_l, QhRecycle = QhRecycle, Tair = Tair, TWaste = TWaste,
				cpm_l = cpm_l, KTloss_l = KTloss_l, dT_EvaporationStandard = dT_EvaporationStandard, minTeh = minTeh,
				cpm_h = cpm_h, KTloss_h = KTloss_h,
				TstorageTankMax = TstorageTankMax, heatStorageVelocity = heatStorageVelocity, heatpumpPowerConstraint = heatpumpPowerConstraint,
				hourlyTariffList = hourlyTariffList, heatConsumptionPowerList = heatConsumptionPowerList,
				maxTcLow = maxTcLow,
				reinitialize = reinitialize,# 将初值重置为T3取最小值，不使用蓄热的情况
				qm1_init = qm1_init,
				qm2_init = qm2_init,
				qm3_init = qm3_init, T3_init = T3_init,
				T8_init = T8_init,
				Tc_l_init = Tc_l_init, P_l1_init = P_l1_init,
				P_l2_init = P_l2_init, P_h1_init = P_h1_init,
				P_h2_init = P_h2_init,
				P_h3_init = P_h3_init, Qc_l_init = Qc_l_init, Qs_l_init = Qs_l_init,
				Qs_h_init = Qs_h_init,
				aimObjectValue=objVal
			)

			# 从0蓄热出发
			isFeasible2, dailyCost2, dfResult2 = generateAndSolve(PressedWaterDoubleStorage(), MinimizeCost();
				COPh = COPh, COPh_g = COPh_g, COPh_h = COPh_h, COPl = COPl, COPl_g = COPl_g, COPl_h = COPl_h,
				T10 = T10, T9 = T9, qm = qm, latentHeat = latentHeat, cp_cs = cp_cs,
				dTe_l1 = dTe_l1, dTe_l2 = dTe_l2, dTc_l = dTc_l, QhRecycle = QhRecycle, Tair = Tair, TWaste = TWaste,
				cpm_l = cpm_l, KTloss_l = KTloss_l, dT_EvaporationStandard = dT_EvaporationStandard, minTeh = minTeh,
				cpm_h = cpm_h, KTloss_h = KTloss_h,
				TstorageTankMax = TstorageTankMax, heatStorageVelocity = heatStorageVelocity, heatpumpPowerConstraint = heatpumpPowerConstraint,
				hourlyTariffList = hourlyTariffList, heatConsumptionPowerList = heatConsumptionPowerList,
				maxTcLow = maxTcLow,
				reinitialize = true,# 将初值重置为T3取最小值，不使用蓄热的情况
				qm1_init = qm1_init,
				qm2_init = qm2_init,
				qm3_init = qm3_init, T3_init = T3_init,
				T8_init = T8_init,
				Tc_l_init = Tc_l_init, P_l1_init = P_l1_init,
				P_l2_init = P_l2_init, P_h1_init = P_h1_init,
				P_h2_init = P_h2_init,
				P_h3_init = P_h3_init, Qc_l_init = Qc_l_init, Qs_l_init = Qs_l_init,
				Qs_h_init = Qs_h_init,
				aimObjectValue=99.0
			)

			if isFeasible2 in [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT] && dailyCost2 < dailyCost
				isFeasible, dailyCost, dfResult = isFeasible2, dailyCost2, dfResult2
			end

			if isFeasible in [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]
				CSV.write(joinpath(filePath, "StorageHour_$(heatStorageTime).csv"), round.(dfResult, digits = 6))
				dailyCostList[heatStorageTime+1] = dailyCost

				lastHeatStorageTime = heatStorageTime
				reinitialize = false
				objVal=dailyCost
				qm1_init = dfResult[!, "直供流量qm1"] |> Vector
				qm2_init = dfResult[!, "高蓄取热流量qm2"] |> Vector
				qm3_init = dfResult[!, "高蓄蓄热流量qm3"] |> Vector

				T3_init = dfResult[!, "低温罐温度T3"] |> Vector
				T8_init = dfResult[!, "高温罐温度T8"] |> Vector
				Tc_l_init = dfResult[!, "低温热泵冷凝器温度Tc_l"] |> Vector


				P_l1_init = dfResult[!, "低温热回收热泵功率P_l1"] |> Vector
				P_l2_init = dfResult[!, "低温空气源热泵功率P_l2"] |> Vector

				P_h1_init = dfResult[!, "高温供热热泵功率P_h1"] |> Vector
				P_h2_init = dfResult[!, "高温取热热泵功率P_h2"] |> Vector
				P_h3_init = dfResult[!, "高温蓄热热泵功率P_h3"] |> Vector

				COPh1_temp = dfResult[!, "高温供热热泵COPh1"] |> Vector
				COPh2_temp = dfResult[!, "高温取热热泵COPh2"] |> Vector
				COPh3_temp = dfResult[!, "高温蓄热热泵COPh3"] |> Vector

				Qc_l_init = dfResult[!, "低温热泵输出热功率Qc_l"] |> Vector

				Qs_l_init = dfResult[!, "低温罐蓄热量Qs_l"] |> Vector
				Qs_h_init = dfResult[!, "高温罐蓄热量Qs_h"] |> Vector

				#dataEconomy=DataFrame("工作时长"=>[],"蓄热时长"=>[],"承压水体积"=>[],"低温热泵总功率"=>[],"高温热泵总功率"=>[],"每天运行费用"=>[])
				#storageTankVolume=cpm_h
				push!(dataEconomy, [24,heatStorageTime,storageTankVolume,maximum(P_l1_init+P_l2_init),maximum(P_h1_init)+maximum(P_h2_init)+maximum(P_h3_init),dailyCost])

				#dataTech=DataFrame("工作时长"=>[],"蓄热时长"=>[],"承压水体积"=>[],"低温热泵总功率"=>[],"高温热泵总功率"=>[],"每天运行费用"=>[],----------"低温热泵1功率"=>[],"低温热泵1最低蒸发温度"=>[],"低温热泵1最高冷凝温度"=>[],"低温热泵2功率"=>[],"低温热泵2最低蒸发温度"=>[],T"低温热泵2最高冷凝温度"=>[],"高温热泵1功率"=>[],"高温热泵1最低蒸发温度"=>[],"高温热泵1最高冷凝温度"=>[],"高温热泵2功率"=>[],"高温热泵2最低蒸发温度"=>[],"高温热泵2最高冷凝温度"=>[],"高温热泵3功率"=>[],"高温热泵3最低蒸发温度"=>[],"高温热泵3最高冷凝温度"=>[])
				push!(dataTech, [24,heatStorageTime,storageTankVolume,maximum(P_l1_init+P_l2_init),maximum(P_h1_init)+maximum(P_h2_init)+maximum(P_h3_init),dailyCost,
				maximum(P_l1_init),
				minimum(dfResult[:,"Te_l1"]),
				maximum(dfResult[:,"低温热泵冷凝器温度Tc_l"]),
				maximum(P_l2_init),
				minimum(dfResult[:,"Te_l2"]),
				maximum(dfResult[:,"低温热泵冷凝器温度Tc_l"]),
				maximum(P_h1_init),
				minimum(dfResult[:,"Te_h1_h3"]),
				maximum(dfResult[:,"用热温度T10"]),
				maximum(P_h2_init),
				minimum(dfResult[:,"Te_h2"]),
				maximum(dfResult[:,"用热温度T10"]),
				maximum(P_h3_init),
				minimum(dfResult[:,"Te_h1_h3"]),
				min(180.0,maximum(dfResult[:,"Tc_h3"])),
				])
			else
				break
			end
			println("T:$Tuse, Storage:$heatStorageTime, Cost:$(dailyCostList[heatStorageTime+1])")

		end
		dfSummary[!, "$(Tuse)℃"] = dailyCostList
		CSV.write(joinpath(filePath2, "summary$(Tuse)℃.csv"), round.(dataEconomy, digits = 6))
		CSV.write(joinpath(filePath2, "techsummary$(Tuse)℃.csv"), round.(dataTech, digits = 6))
		
	end
	CSV.write(joinpath(filePath0, "summary$(Dates.format(now(), "HH_MM")).csv"), round.(dfSummary, digits = 6))
end

main(TuseList)

