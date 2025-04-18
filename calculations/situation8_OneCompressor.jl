using Revise
using Pkg
using HeatPumpWithStorageSystem
using DataFrames, CSV

#=
全天工作
=#

function main()
	situation = "situation8"
	hourlyTariff = zeros(24)
	hourlyTariff[1:7] .= 0.3340
	hourlyTariff[8:11] .= 0.7393
	hourlyTariff[12:14] .= 1.2360
	hourlyTariff[15:18] .= 0.7393
	hourlyTariff[19:23] .= 1.2360
	hourlyTariff[24] = 0.3340
	# 系数
	heatPumpServiceCoff = 1.0
	maxCOP = 21# 最大COP
	eta_s = 0.7# 绝热效率
	workingStartHour = 0                # 生产开始时间
	workingHours = 24                   # 每日工作小时数
	TWaste = 30.0                     # 废热源温度
	Tair = 25.0                        # 外部环境温度
	TCompressorIn = 115.0
	maxTcHigh = 180.0
	dT_EvaporationStandard = 5.0
	heatConsumptionPower = 1.0

	# 计算参数
	dT = 0.1
	dt = 1# 时间步长过小会导致初始温度优化的目标不是一个单峰函数
	K = 2
	#lambdaPe=0.01

	# 可调参数：循环工质、用热温度、蓄热容量、热泵服务系数、电锅炉服务系数、中间级温度、废热温度
	overlapRefrigerantList = [NH3_Water, R1233zdE_Water]
	heatStorageCapacityList = 1.0:1.0:10.0
	TuseList = 120.0:20.0:180.0
	totalCalculationTime=length(overlapRefrigerantList)*length(heatStorageCapacityList)*length(TuseList)

	COPWater = getCOP(
		TCompressorIn,# 蒸发温度下限,这里是实际设计中的蒸发冷凝温度界限
		maxTcHigh,# 蒸发温度上限
		TCompressorIn,# 冷凝温度下限
		maxTcHigh,# 冷凝温度上限
		refWater,# 工质
		maxCOP,# 最大COP
		eta_s,# 绝热效率
		dT,# 插值步长
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



	count=0
	for or in overlapRefrigerantList
		filePathOr = joinpath(filePath0, or.refrigerant)
		filePathEconomic = joinpath(filePathOr, "economic")
		dfCost = DataFrame("蓄热容量" => vcat(0.0,heatStorageCapacityList))# 不同用热和蓄热的运行成本

		# 复叠循环COP
		COPOverlap = getOverlapCOP_fixMidTemperature(
			or,
			TCompressorIn + or.midTDifference / 2;
			maxCOP = maxCOP,# 最大COP
			eta_s = eta_s,# 绝热效率
			dT = dT,# 插值步长
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
			dT,# 插值步长
		)

		Threads.@threads for Tuse in TuseList
			COPLow = COPLowFunction(TWaste, TCompressorIn+or.midTDifference)
			filePathTuse = joinpath(filePathOr, "Tuse_" * string(Tuse))
			dfEconomic = DataFrame([(col => Float64[]) for col in column_names_eco]...)# 经济性指标

			temp=COPWater(TCompressorIn,Tuse)
			temp2=COPOverlap(TWaste,Tuse)
			PLowMAX=heatConsumptionPower*(temp-1)/(COPLow+temp-1)/temp2
			PHighMAX=heatConsumptionPower*(COPLow)/(COPLow+temp-1)/temp2
			
			temp = PLowMAX+PHighMAX
			push!(dfEconomic,[
				workingHours,
				0.0,
				0.0,
				PLowMAX*temp2,
				PHighMAX*temp2,
				0.0,
				temp*sum(hourlyTariff),
				temp*24
			])

			for heatStorageCapacity in heatStorageCapacityList
				dfOperation = DataFrame()# 运行参数
				# 生成计算参数
				hourlyTariffFunction, heatConsumptionPowerFunction, TairFunction,
				Tuse, TCompressorIn,
				dT_EvaporationStandard,
				latentHeat, cp_cw, cp_cs,
				TcChangeToElec, TWaste,
				cpm_h,
				TstorageTankMax, PheatPumpMax, PelecHeatMax = generateSystemCoff(PressedWaterDoubleStorageOneCompressor();
					overlapRefrigerant = or,    # 复叠工质
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

				minCostGo, minTsListGo, P1ListGo, P2ListGo, P3ListGo, PeListGo = generateAndSolve(PressedWaterDoubleStorageOneCompressor(), MinimizeCost(), ConstloadandArea(), GoldenRatioMethod();
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
					K=K,
					#lambdaPe=lambdaPe,
				)

				storageTankMass = cpm_h / cp_cw * 3600#kg	蓄热水质量
				storageTankVolume = storageTankMass / 900#m^3 蓄热罐体积
				tList = collect(0:dt:(24-dt).+0.5*dt)
				nt = length(tList)

				COPLowList = fill(1.0, nt)
				for i ∈ 1:nt
					if (P1ListGo[i] > 0 || P3ListGo[i] > 0)
						COPLowList[i] = COPLow
					end
				end

				COPOverlaplist = fill(COPOverlap(TWaste, Tuse), nt)
				for i ∈ 1:nt
					if P2ListGo[i] > 0
						COPOverlaplist[i] = COPWater((minTsListGo[i+1]+minTsListGo[i+1])/2 - dT_EvaporationStandard, Tuse)
					elseif P3ListGo[i] > 0
						COPOverlaplist[i] = COPOverlap(TWaste, max((minTsListGo[i+1]+minTsListGo[i+1])/2 + dT_EvaporationStandard, Tuse))
					end
				end

				COPWaterList = fill(COPOverlap(TWaste, Tuse), nt)
				for i ∈ 1:nt
					if P2ListGo[i] > 0
						COPWaterList[i] = COPWater(P2ListGo[i] - dT_EvaporationStandard, Tuse)
					elseif P1ListGo[i] > 0 && P3ListGo[i] == 0
						COPWaterList[i] = COPWater(TCompressorIn, Tuse)
					elseif P1ListGo[i] > 0 && P3ListGo[i] > 0
						COPWaterList[i] = COPWater(TCompressorIn, max((minTsListGo[i+1]+minTsListGo[i+1])/2 + dT_EvaporationStandard, Tuse))
					else
						COPWaterList[i] = 1.0
					end
				end
				#println(length(P1ListGo)," ",length(COPWaterList)," ",nt)
				PLowList = (P1ListGo + P3ListGo) .* (COPWaterList .- 1) ./ (COPLowList + COPWaterList .- 1)
				PHighList = P1ListGo + P2ListGo + P3ListGo - PLowList

				dfOperation[!, :时间] = tList

				dfOperation[!, :环境温度] = fill(Tair, nt)
				dfOperation[!, :蒸发器温度] = fill(TWaste, nt)
				dfOperation[!, :用热温度] = fill(Tuse, nt)
				dfOperation[!, :用热负载] = fill(heatConsumptionPower, nt)
				dfOperation[!, :蓄热温度] = minTsListGo[1:end-1]

				dfOperation[!, :热泵直接供热功率] = P1ListGo
				dfOperation[!, :蓄热供热功率] = P2ListGo
				dfOperation[!, :热泵储热供热功率] = P3ListGo
				dfOperation[!, :电加热储热供热功率] = PeListGo


				dfOperation[!, :低温热泵COP] = COPLowList
				dfOperation[!, :水蒸气压缩机COP] = COPWaterList
				dfOperation[!, :复叠COP] = COPOverlaplist

				dfOperation[!, :低温热泵功率] = PLowList
				dfOperation[!, :水蒸气压缩机功率] = PHighList

				PLowMAX=maximum(PLowList)
				PHighMAX=maximum(PHighList)
				temp = PLowMAX+PHighMAX
				if temp <1
					PLowMAX /= temp
					PHighMAX /= temp
				end
				#[:工作时长, :蓄热时长, :承压水体积, :低温热泵总功率, :高温热泵总功率, :电加热功率, :每天运行费用, :总电度]
				push!(dfEconomic, [
					workingHours,
					heatStorageCapacity,
					storageTankVolume,
					PLowMAX,
					PHighMAX,
					PelecHeatMax,
					minCostGo,
					sum(PLowList + PHighList + PeListGo) * dt,
				])
				CSV.write(joinpath(filePathTuse, string(Tuse) * "_" * string(heatStorageCapacity) * ".csv"), round.(dfOperation,digits=5))
				count+=1
				println(count,"/",totalCalculationTime," ","Tuse=",Tuse,", heatStorageCapacity=",heatStorageCapacity)
			end
			
			sort!(dfEconomic, "蓄热时长")
			CSV.write(joinpath(filePathEconomic, "经济环境指标" * string(Tuse) * "℃.csv"), round.(dfEconomic,digits=5))

			dfCost[!, string(Tuse)*"℃"] = dfEconomic[:, "每天运行费用"]
		end
		select!(dfCost, column_names_temperature)
		CSV.write(joinpath(filePathOr, "summary.csv"), round.(dfCost,digits=5))
	end
end
#2小时20分钟
# 2025年3月17日3时开始计算
@info "开始计算"
@time main()
