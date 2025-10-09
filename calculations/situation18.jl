using Revise
using Pkg
using HeatPumpWithStorageSystem
using DataFrames, CSV

#=
全天工作
=#

function main()
	situation = "situation18"
	hourlyTariff = ones(48) * 0.7393
	p = 1.7
	pp = p * 1.2
	v = 0.35
	vv = v * 0.8

	hourlyTariff[1:12] .*= v
	hourlyTariff[23:26] .*= v
	hourlyTariff[29:30] .*= pp
	hourlyTariff[31:39] .*= p
	hourlyTariff[40:43] .*= pp
	hourlyTariff[44] *= p

	Tair = vcat(
		fill(26.0, 7),
		fill(27.0, 7),
		fill(26.0, 11),
	)

	#heatConsumptionPower = fill(1.0,48)

	heatConsumptionPower = vcat(
		fill(0.0, 16),
		fill(1.0, 8),
		fill(0.0, 2),
		fill(1.0, 8),
		fill(0.0, 14),
	)
	#=
	heatConsumptionPower = vcat(
		fill(0.0,8),
		fill(1.0, 4),
		fill(0.0,1),
		fill(1.0,4),
		fill(0.0,7)
	)

	heatConsumptionPower = vcat(
		fill(0.0,4),
		fill(1.0, 2),
		fill(1.0,2),
		fill(0.0,4)
	)
	=#

	# 系数
	heatPumpServiceCoff = 1.0
	maxCOP = 21# 最大COP
	eta_s = 0.7# 绝热效率
	workingStartHour = 0                # 生产开始时间
	workingHours = 24                   # 每日工作小时数
	TWaste = 30.0                     # 废热源温度
	#Tair = 25.0                        # 外部环境温度
	TCompressorIn = 115.0
	maxTcHigh = 180.0
	dT_EvaporationStandard = 5.0
	Tsmin = 120

	# 计算参数
	dT = 1e-3
	#dt = 1/2# 时间步长过小会导致初始温度优化的目标不是一个单峰函数

	dt = 0.5
	smoother = 1e-8

	nt = Int(24 / dt + 1)

	# 可调参数：循环工质、用热温度、蓄热容量、热泵服务系数、电锅炉服务系数、中间级温度、废热温度
	overlapRefrigerantList = [
		NH3_Water,
		#R1233zdE_Water
	]

	# 热容的计算列表
	#heatStorageCapacityList = 1.0:1.0:10.0
	heatStorageCapacityList = 0.0:1.0:0.0
	# 用热温度的计算列表
	#TuseList = 130.0:10.0:180.0
	TuseList = 130.0:10.0:180.0

	# 总算例数
	totalCalculationTime = length(overlapRefrigerantList) * length(heatStorageCapacityList) * length(TuseList)

	COPWater = getCOP(
		TCompressorIn,# 蒸发温度下限,这里是实际设计中的蒸发冷凝温度界限
		maxTcHigh,# 蒸发温度上限
		TCompressorIn,# 冷凝温度下限
		maxTcHigh,# 冷凝温度上限
		refWater,# 工质
		maxCOP,# 最大COP
		eta_s,# 绝热效率
		1.0,#dT 插值步长
	)

	column_names_temperature = vcat("蓄热容量", string.(TuseList) .* "℃")
	column_names_eco = string.([:工作时长, :蓄热时长, :承压水体积, :低温热泵总功率, :高温热泵总功率, :电加热功率, :每天运行费用, :总电度])


	# 先生成路径
	filePath0 = joinpath(pwd(), "calculations", situation)
	if !isdir(filePath0)
		mkdir(filePath0)
	end
	filePath0 = joinpath(filePath0, "WorkAllDay")
	if !isdir(filePath0)
		mkdir(filePath0)
	end
	for or in overlapRefrigerantList
		filePathOr = joinpath(filePath0, or.refrigerant)# 存放一个工质的所有计算数据
		if !isdir(filePathOr)
			mkdir(filePathOr)
		end
		# 存放不同温度的经济性指标：工作时长,蓄热时长,承压水体积,低温热泵总功率,高温热泵总功率,电加热功率,每天运行费用,总电度
		filePathEconomic = joinpath(filePathOr, "economic")
		if !isdir(filePathEconomic)
			mkdir(filePathEconomic)
		end
		for Tuse in TuseList
			# 存放系统的运行参数
			filePathTuse = joinpath(filePathOr, "Tuse_" * string(Tuse))
			if !isdir(filePathTuse)
				mkdir(filePathTuse)
			end
		end
	end



	count = 0
	for or in overlapRefrigerantList
		filePathOr = joinpath(filePath0, or.refrigerant)
		filePathEconomic = joinpath(filePathOr, "economic")
		dfCost = DataFrame("蓄热容量" => heatStorageCapacityList)# 不同用热和蓄热的运行成本

		# 复叠循环COP
		COPOverlap = getOverlapCOP_fixMidTemperature(
			or,
			TCompressorIn + or.midTDifference / 2;
			maxCOP = maxCOP,# 最大COP
			eta_s = eta_s,# 绝热效率
			dT = 1.0,# 插值步长
		)

		# 低温热泵COP
		refLow = or.refrigerantLow
		COPLowFunction = getCOP(
			refLow.minTe,# 蒸发温度下限,这里是实际设计中的蒸发冷凝温度界限
			refLow.maxTe,# 蒸发温度上限
			refLow.minTc,# 冷凝温度下限
			refLow.maxTc,# 冷凝温度上限
			refLow,# 工质
			maxCOP,# 最大COP
			eta_s,# 绝热效率
			1.0,# 插值步长
		)

		#Threads.@threads for Tuse in TuseList
		for Tuse in TuseList
			COPLow = COPLowFunction(TWaste, TCompressorIn + or.midTDifference)
			filePathTuse = joinpath(filePathOr, "Tuse_" * string(Tuse))
			dfEconomic = DataFrame([(col => Float64[]) for col in column_names_eco]...)# 经济性指标

			temp = COPWater(TCompressorIn, Tuse)
			temp2 = COPOverlap(TWaste, Tuse)
			PLowMAX = maximum(heatConsumptionPower) * (temp - 1) / (COPLow + temp - 1) / temp2
			PHighMAX = maximum(heatConsumptionPower) * (COPLow) / (COPLow + temp - 1) / temp2

			#=
			无储能的曲线计算时需要修改电加热逻辑
			=#
			temp = PLowMAX + PHighMAX
			hf = generateGridPriceFunction(hourlyTariff, 24)
			#=
			push!(dfEconomic,[
				workingHours,
				0.0,
				0.0,
				PLowMAX,
				PHighMAX,
				0.0,
				sum(temp.*hf.(0:dt:(24-dt)).*heatConsumptionPower)*dt,
				temp*sum(heatConsumptionPower)*dt
			])
			=#

			for heatStorageCapacity in heatStorageCapacityList
				dfOperation = DataFrame()# 运行参数
				COP1_design = COPOverlap(TWaste, Tuse)
				COPWater_design = COPWater(TCompressorIn, Tuse)
				# 生成计算参数
				hourlyTariffFunction, heatConsumptionPowerFunction, TairFunction,
				Tuse, TCompressorIn,
				dT_EvaporationStandard,
				latentHeat, cp_cw, cp_cs,
				TcChangeToElec, TWaste,
				cpm_h,
				TstorageTankMax, PheatPumpMax, PelecHeatMax, PWaterCompressorMax = generateSystemCoff(PressedWaterOneStorageOneCompressor();
					overlapRefrigerant = or,    # 复叠工质
					COP1_design = COP1_design,
					COPWater_design = COPWater_design,
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
				heatLoadList = heatConsumptionPowerFunction.(0:dt:24)
				TairList = TairFunction.(0:dt:24)
				costGridList = hourlyTariffFunction.(0:dt:24)

				# 开始计算，重复5次取最小值
				minCostGo = 9999.0
				minTsListGo = []
				P1ListGo = []
				P2ListGo = []
				P3ListGo = []
				PeListGo = []
				realCostListGo = []
				for _ ∈ 1:1
					minCostGotemp, minTsListGotemp, P1ListGotemp, P2ListGotemp, P3ListGotemp, PeListGotemp,realCostListtemp = generateAndSolve(PressedWaterOneStorageOneCompressor(), MinimizeCost(), VaryLoadVaryArea(), GoldenRatioMethod();
						COPOverlap = COPOverlap,
						COPWater = COPWater,
						COPLowFunction = COPLowFunction,
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
						PWaterCompressorMax = PWaterCompressorMax,
						Tsmin = Tsmin,
						# 求解参数
						dT = dT,# 状态参数高温蓄热温度离散步长
						dt = dt,# 时间步长
						#lambdaPe=lambdaPe,
						smoother = smoother,
					)
					# 发现更小解，就选择更小解
					if minCostGotemp < minCostGo
						minCostGo = minCostGotemp
						minTsListGo = minTsListGotemp
						P1ListGo = P1ListGotemp
						P2ListGo = P2ListGotemp
						P3ListGo = P3ListGotemp
						PeListGo = PeListGotemp
						realCostListGo = realCostListtemp
					end
				end

				# 验证结果是否满足能量守恒
				caseName = or.refrigerant * "_" * string(Tuse) * "_" * string(heatStorageCapacity) * "h"
				params = SystemParameters(
					ThMax = TcChangeToElec,
					Tuse = Tuse,
					dT = dT_EvaporationStandard,
					TCompressorIn = TCompressorIn,
					cpm = cpm_h,
					COPWater = COPWater,
					PhMax = PheatPumpMax,
					PeMax = PelecHeatMax,
				)

				isEnergyConserved = true
				for i ∈ 1:nt-1
					sysVariables = HeatPumpWithStorageSystem.SystemVariables(
						heatLoadList[i],
						COPLowFunction(TWaste, TCompressorIn + dT_EvaporationStandard),
						TairList[i],
						TWaste,
					)
					cost_test, flag_test, P1_test, P2_test, P3_test, Pe_test = getMinimumCost(minTsListGo[i], minTsListGo[i+1], dt, params, sysVariables)

					if !flag_test
						println(caseName, "第$(i)个时间段有误")
						isEnergyConserved = false
					end
					if cost_test * costGridList[i] - realCostListGo[i] > 1e-5
						println(caseName, "第$(i)个时间段的成本有误")
						isEnergyConserved = false
					end
					if P1_test - P1ListGo[i] > 1e-5
						println(caseName, "第$(i)个时间段P1ListGo[i]和P1_test有误")
						isEnergyConserved = false
					end
					if P2_test - P2ListGo[i] > 1e-5
						println(caseName, "第$(i)个时间段P2ListGo[i]和P2_test有误")
						isEnergyConserved = false
					end
					if P3_test - P3ListGo[i] > 1e-5
						println(caseName, "第$(i)个时间段P3ListGo[i]和P3_test有误")
						isEnergyConserved = false
					end
					if Pe_test - PeListGo[i] > 1e-5
						println(caseName, "第$(i)个时间段PeListGo[i]和Pe_test有误")
						isEnergyConserved = false
					end

				end
				if !isEnergyConserved
					println(caseName, "能量守恒有误！")
				end

				x1List = P1ListGo.>0 .|> Int
				x2List = P2ListGo.>0 .|> Int
				x3List = P3ListGo.>0 .|> Int
				storageTankMass = cpm_h / cp_cw * 3600#kg	蓄热水质量
				storageTankVolume = storageTankMass / 900#m^3 蓄热罐体积
				tList = collect(0:dt:(24-dt).+0.5*dt)

				# 这里nt的数值变了
				nt = length(tList)

				COPLowList = fill(1.0, nt)
				for i ∈ 1:nt
					if (P1ListGo[i] > 0 || P3ListGo[i] > 0)
						COPLowList[i] = COPLow
					end
				end

				#=这是用整数规划求解器做的，目前整数规划求解器有问题
				COPh1List = fill(COPWater(TCompressorIn,Tuse),nt)
				COPh2List = [COPWater(Ts-dT_EvaporationStandard,Tuse) for Ts in minTsListGo[1:nt]]
				COPh3List = [COPWater(TCompressorIn, Ts+dT_EvaporationStandard) for Ts in minTsListGo[1:nt]]

				minCOP1COP3 = min.(COPh1List, COPh3List)
				=#
				COPh1List = zeros(nt)
				COPh2List = zeros(nt)
				COPh3List = zeros(nt)
				COPOverlapList = zeros(nt)
				for i ∈ 1:nt
					sysVariables = HeatPumpWithStorageSystem.SystemVariables(
						heatLoadList[i],
						COPLowFunction(TWaste, TCompressorIn + dT_EvaporationStandard),
						TairList[i],
						TWaste,
					)
					_, COPh1List[i], COPh2List[i], COPh3List[i], COPOverlapList[i] = getCOPbyMode(x1List[i], x2List[i], x3List[i], minTsListGo[i], minTsListGo[i+1], params, sysVariables)
				end

				#println(length(P1ListGo)," ",length(COPWaterList)," ",nt)
				PHighList = zeros(nt)
				PLowList = zeros(nt)
				PPumpList = zeros(nt)
				for i ∈ 1:nt
					# 先算模式1和3的高温热泵功率
					PHighList[i] = P1ListGo[i] * COPOverlapList[i] / COPh1List[i] + P3ListGo[i] * COPOverlapList[i] / COPh3List[i]

					# 低温热泵功率只有1，3开启了
					PLowList[i] = P1ListGo[i] + P3ListGo[i] - PHighList[i]

					# 如果同时开启蓄热供热与热泵直接供热，P2不用往里加
					if x2List[i] == 1 && x1List[i] == 0
						PHighList[i] += P2ListGo[i]
					else
						PPumpList[i] = P2ListGo[i]
					end
				end

				dfOperation[!, :时间] = tList
				dfOperation[!, :电价] = costGridList[1:nt]

				dfOperation[!, :环境温度] = TairFunction.(tList)
				dfOperation[!, :蒸发器温度] = fill(TWaste, nt)
				dfOperation[!, :用热温度] = fill(Tuse, nt)
				dfOperation[!, :用热负载] = hf.(tList)
				dfOperation[!, :蓄热温度] = minTsListGo[1:end-1]

				dfOperation[!, :热泵直接供热功率] = P1ListGo
				dfOperation[!, :蓄热供热功率] = P2ListGo
				dfOperation[!, :热泵储热供热功率] = P3ListGo
				dfOperation[!, :电加热储热供热功率] = PeListGo

				dfOperation[!, :低温热泵COP] = COPLowList
				dfOperation[!, :COPh1] = COPh1List
				dfOperation[!, :COPh2] = COPh2List
				dfOperation[!, :COPh3] = COPh3List
				dfOperation[!, :复叠COP] = COPOverlapList

				dfOperation[!, :低温热泵功率] = PLowList
				dfOperation[!, :水蒸气压缩机功率] = PHighList
				dfOperation[!, :水泵功率] = PPumpList

				PLowList + PHighList + PPumpList - P1ListGo - P3ListGo - P2ListGo

				dfOperation[!, :蓄热储入功率反馈] = (minTsListGo[2:end] - minTsListGo[1:end-1]) * cpm_h / dt

				dfOperation[!, :蓄热储入功率计算] = P3ListGo .* COPOverlapList + PeListGo - (P2ListGo + PPumpList) .* (COPh2List .- 1)

				dfOperation[!, :蓄热储入功率反馈误差] = dfOperation[!, :蓄热储入功率反馈] - dfOperation[!, :蓄热储入功率计算]

				PLowMAX = maximum(PLowList)
				PHighMAX = maximum(PHighList)
				temp = PLowMAX + PHighMAX

				# if temp <1
				# 	PLowMAX /= temp
				# 	PHighMAX /= temp
				# end
				#[:工作时长, :蓄热时长, :承压水体积, :低温热泵总功率, :高温热泵总功率, :电加热功率, :每天运行费用, :总电度]

				push!(dfEconomic, [
					workingHours,
					heatStorageCapacity,
					storageTankVolume,
					PLowMAX,
					PHighMAX,
					maximum(PeListGo),
					minCostGo,
					sum(PLowList + PHighList + PeListGo) * dt,
				])

				CSV.write(joinpath(filePathTuse, string(Tuse) * "_" * string(heatStorageCapacity) * ".csv"), round.(dfOperation, digits = 5))

				count += 1

				println(count, "/", totalCalculationTime, " ", "Tuse=", Tuse, ", heatStorageCapacity=", heatStorageCapacity)
			end

			sort!(dfEconomic, "蓄热时长")
			CSV.write(joinpath(filePathEconomic, "经济环境指标" * string(Tuse) * "℃.csv"), round.(dfEconomic, digits = 5))

			dfCost[!, string(Tuse)*"℃"] = dfEconomic[:, "每天运行费用"]
		end
		select!(dfCost, column_names_temperature)
		CSV.write(joinpath(filePathOr, "summary.csv"), round.(dfCost, digits = 5))
	end
end
#2小时20分钟
# 2025年3月17日3时开始计算
@info "开始计算"
@time main()

