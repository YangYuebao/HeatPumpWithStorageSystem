
abstract type FunctionGenerateMethod end
struct BackwardGenerate <: FunctionGenerateMethod end
struct LinearGenerate <: FunctionGenerateMethod end
struct BackwardGenerateCycled <: FunctionGenerateMethod end
struct LinearGenerateCycled <: FunctionGenerateMethod end

function functionInterpolationGenerator(list::Vector, duration::Real, ::BackwardGenerate)
	n = length(list)
	dt = duration / n
	function f(x::Real)
		if x >= duration
			x = x % duration
		elseif x < 0
			x = x % duration + duration
		end
		return list[floor(Int, x / dt)+1]
	end
	return f
end
functionInterpolationGenerator(list::Vector, duration::Real, ::BackwardGenerateCycled) = functionInterpolationGenerator(list, duration, BackwardGenerate())

function functionInterpolationGenerator(list::Vector, duration::Real, ::LinearGenerate)
	n = length(list)
	if n == 1
		list = [list[1], list[1]]
	end
	n = 2
	tList = collect(range(0, stop = duration, length = n))
	sitpCOP = linear_interpolation(tList, list)
	function f(x::Real)
		if x >= duration
			x = x % duration
		elseif x < 0
			x = x % duration + duration
		end
		return sitpCOP(x)
	end
	return f
end

function functionInterpolationGenerator(list::Vector, duration::Real, ::LinearGenerateCycled)
	if abs(list[1] - list[end]) < 1e-6
		return functionInterpolationGenerator(list, duration, LinearGenerate())
	else
		@warn "列表首尾不一致，已经自动添加一位"
		return functionInterpolationGenerator(vcat(list, list[1]), duration, LinearGenerate())
	end
end


"""生成电价函数"""
function generateGridPriceFunction(gridPriceList::Vector, duration::Real)
	functionInterpolationGenerator(gridPriceList, duration, BackwardGenerate())
end

"""生成负载需求函数"""
function generateLoadFunction(loadList::Vector, duration::Real)
	functionInterpolationGenerator(loadList, duration, LinearGenerateCycled())
end

"""生成环境温度函数"""
function generateAreaTemperatureFunction(temperatureList::Vector, duration::Real)
	functionInterpolationGenerator(temperatureList, duration, LinearGenerateCycled())
end

"""
输入一定系统结构和工作参数,返回系统计算需要用到的各种向量
"""
function generateSystemCoff(::PressedWaterDoubleStorageOneCompressor;
	overlapRefrigerant::OverlapRefrigerant = NH3_Water,# 复叠工质
	maxTcHigh::Real = 180.0,  # 高温热泵冷凝器温度上限
	TCompressorIn::Real = 115.0,# 中间温度
	TWaste::Real = 60.0,                  # 废热源温度
	Tair::Real = 25.0,          # 外部环境温度
	Tuse::Real = 120.0,                     # 工厂使用温度
	heatStorageCapacity::Real = 6.0,        # 蓄热量kWh(承压水蓄热)
	TstorageTankMax::Real = 220.0,          # 蓄热罐的最高温度
	maxheatStorageInputHour::Real = 4,      # 蓄热充满时长
	dT_EvaporationStandard::Real = 5.0,# 全蒸温差
	heatConsumptionPower::Real = 1.0,       # 每小时用热功率kW
	workingStartHour::Int = 8,              # 生产开始时间
	workingHours::Int = 16,                 # 每日工作小时数
	heatPumpServiceCoff::Real = 1.2,        # 热泵功率服务系数
	hourlyTariff::Vector = fill(0.7, 24),   # 电价向量
)
	# COP函数外部生成
	# 计算工作时间
	temp = workingHours
	workingHours = workingHours % 24
	if workingHours == 0
		workingHours = 24
	end
	if temp > 24
		@warn "工作时长大于24小时,改为$(workingStartHour)h,终止计算"
		return -1
	end

	# 计算热泵基础功率和实际配置功率
	#P1base = heatConsumptionPower / COPOverlap(TWaste, Tuse)
	PheatPumpMax = heatConsumptionPower * heatPumpServiceCoff
	PelecHeatMax = heatConsumptionPower * heatStorageCapacity / maxheatStorageInputHour

	# 生成需求与环境函数
	hourlyTariffFunction = generateGridPriceFunction(hourlyTariff, 24)
	heatConsumptionPowerFunction = generateLoadFunction([heatConsumptionPower], 24)
	TairFunction = generateAreaTemperatureFunction([Tair], 24)

	#=
	heatConsumptionPowerList = zeros(24)
	heatConsumptionPowerList[1:workingHours] .= heatConsumptionPower
	heatConsumptionPowerList = rotateVectorForward(heatConsumptionPowerList, workingStartHour)
	=#

	latentHeat = 2150.0# 汽化潜热kJ/kg
	cp_cw = 4.275# 循环水定压热容kJ/kg
	cp_cs = 2.281# 循环蒸汽定压热容kJ/kg


	# 生成低温热泵参数：cpqm_l,k1,dTe_l1,dTe_l2,dTc_l,QhRecycle
	cpm_h = heatStorageCapacity * heatConsumptionPower / (TstorageTankMax - Tuse)# kWh/K
	storageTankMass = cpm_h / cp_cw * 3600#kg
	storageTankVolume = storageTankMass / 900#m^3

	#Qe_hStandard=1.0
	dT_EvaporationStandard = dT_EvaporationStandard

	# 设备运行约束:TstorageTankMax,heatStorageVelocity,heatStorageCapacityConstraint,heatpumpPowerConstraint
	heatStorageVelocity = heatStorageCapacity / maxheatStorageInputHour * heatConsumptionPower# kWh/h
	heatStorageCapacityConstraint = heatStorageCapacity * heatConsumptionPower# kWh
	heatpumpPowerConstraint = PheatPumpMax# kW


	# 其它参数:hourlyTariffList,heatConsumptionPowerList
	TcChangeToElec = maxTcHigh

	return hourlyTariffFunction, heatConsumptionPowerFunction, TairFunction,
	Tuse, TCompressorIn,
	dT_EvaporationStandard,
	latentHeat, cp_cw, cp_cs,
	TcChangeToElec, TWaste,
	cpm_h,
	TstorageTankMax, PheatPumpMax, PelecHeatMax
end

abstract type SimplifiedType end
struct NoSimplify <: SimplifiedType end
struct ConstloadandArea <: SimplifiedType end

# 计算最优的0点温度时需要考虑优化方法
abstract type OptimizeMethod end
struct ExhaustiveMethod <: OptimizeMethod end# 穷举法（全面的）
struct GoldenRatioMethod <: OptimizeMethod end# 0.618法
struct MomentumMethod <: OptimizeMethod end# 动量法

"""计算状态转移成本矩阵"""
function getStateTransitionCost(::PressedWaterDoubleStorageOneCompressor;
	COPOverlap::Function,
	COPWater::Function,
	heatLoad::Real,#热负荷
	Tair::Real,# 外部环境温度
	# 总循环参数
	Tuse::Real,# 供热蒸汽温度
	dT_EvaporationStandard::Real,#全蒸温差
	latentHeat::Real,# 汽化潜热
	cp_cw::Real,# 循环水定压热容
	TcChangeToElec::Real,
	TWaste::Real,# 废热回收蒸发器温度

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容

	# 设备运行约束
	PheatPumpMax::Real,# 热泵最大功率
	PelecHeatMax::Real,# 电锅炉最大功率
	# 求解参数
	TsList::Vector,# 状态参数高温蓄热温度列表
	dt::Real = 0.1,# 时间步长
)
	COP1 = COPOverlap(TWaste, Tuse)
	COP2 = COPWater
	COP3 = COPOverlap

	nT = length(TsList)# 温度步数	

	#C[i,j]表示温度从TsList[i]到TsList[j]时的最低功率；在负载、环境不变的情况下，C[i,j]是不变的
	C = fill(99.0, nT, nT)# 状态转移矩阵
	P1Matrix = zeros(nT, nT)# 状态转移功率1参数
	P2Matrix = zeros(nT, nT)# 状态转移功率2参数
	P3Matrix = zeros(nT, nT)# 状态转移功率3参数
	PeMatrix = zeros(nT, nT)# 状态转移功率电加热参数
	TsDecreaseIndexList = zeros(Int, nT)# 记录温度下降最多偏移的index
	TsIncreaseIndexList = zeros(Int, nT)# 记录温度上升最多偏移的index

	# 只用热泵供热时的功率
	P1Only = heatLoad / COP1# 只用热泵供热时的功率
	for i ∈ 1:nT
		P1Matrix[i, i] = P1Only
	end

	"""计算蓄热温度相对较低时的电度情况"""
	function powerCalculate_lowTs(Tsaim, Tsstart, dt)
		P3 = min(PheatPumpMax - P1Only, cpm_h / dt / COP1 * (Tsaim - Tsstart))
		Pe = cpm_h / dt * (Tsaim - Tsstart) - P3 * COP1
		flag = Pe <= PelecHeatMax
		return (P1Only + P3 + Pe) * dt, flag, P1Only, P3, Pe
	end

	"""计算蓄热温度相对较高时的电度情况"""
	function powerCalculate_highTs(Tsaim, Tsstart, dt)
		COP3value = COP3(TWaste, (Tsaim + Tsstart) / 2 + dT_EvaporationStandard)
		Ptotal = 999.0
		P1res = 0.0
		P3res = 0.0
		Peres = 0.0
		flag = false
		# 模式1：电热加热
		Pe = cpm_h / dt * (Tsaim - Tsstart)
		if Pe <= PelecHeatMax
			flag = true
			P1res = P1Only
			P3res = 0.0
			Peres = Pe
			Ptotal = P1Only + Pe
		end

		# 模式2：热泵加电热
		P3 = min(
			PheatPumpMax,
			(PheatPumpMax - heatLoad / COP3value) / (1 - cp_cw / latentHeat * (Tsaim - Tsstart)),
			cpm_h / dt / COP3value * (Tsaim - Tsstart),
		)
		P1 = max(0, heatLoad / COP3value - cp_cw / latentHeat * (Tsaim - Tsstart) * P3)
		Pe = cpm_h / dt * (Tsaim - Tsstart) - P3 * COP3value
		if (P1 + P3 + Pe < Ptotal) && (Pe <= PelecHeatMax) && (Tsaim <= TcChangeToElec - dT_EvaporationStandard)
			flag = true
			P1res = P1
			P3res = P3
			Peres = Pe
			Ptotal = P1 + P3 + Pe
		end
		return Ptotal * dt, flag, P1res, P3res, Peres
	end

	for i ∈ 1:nT
		# 先计算温度下降的功率,存在一边电加热一边开2号模式的情况
		# 首先计算该模式在下一刻的最低温度
		TsNextMin = TsList[i]
		# 简单迭代求蓄热的最低温度
		for _ ∈ 1:3
			TsNextMin = TsList[i] - heatLoad * dt / cpm_h * (1 - 1 / COP2(0.5 * (TsNextMin + TsList[i]) - dT_EvaporationStandard, Tuse))
		end
		# 计算最低温度的下标
		ibias=findfirst(x -> x >= TsNextMin, TsList)
		TsDecreaseIndexList[i]=i-ibias
		TsNextMinGrid = TsList[ibias]
		elecP = (TsNextMinGrid - TsNextMin) * cpm_h

		index = i - TsDecreaseIndexList[i]
		#=
		# 纯蓄热供热的工况下不考虑电加热的能量守恒会导致严重的计算错误,但在时间步长减小时由于温度降低幅度减小，又要依附到网格上，导致电加热的补齐效应一直开启。
		P2 = heatLoad / (COP2((TsList[i] + TsList[index]) / 2 - dT_EvaporationStandard, Tuse))
		C[i, index] = P2 * dt
		P2Matrix[i, index] = P2
		=#

		for j ∈ index:i-1
			P2 = heatLoad / (COP2((TsList[i] + TsList[j]) / 2 - dT_EvaporationStandard, Tuse))
			elecP = (TsList[j] - TsNextMin) * cpm_h
			#补热的功率基本上超不了，因为这个工况下从蓄热取走的热量小于工厂的负载，而补热的功率大于等于工厂的负载
			if elecP < PelecHeatMax * dt
				C[i, j] = P2 * dt + elecP
				P2Matrix[i, j] = P2
				PeMatrix[i, j] = elecP / dt
			end
		end

		# 计算温度不变的电度
		C[i, i] = P1Only * dt

		# 计算温度上升的功率，首先计算该模式的转换温度。升温模式有两种工作状态：1.热泵直供，电加热提高蓄热温度；2.热泵按照两者中高的供热
		# 需要一个计算温度低于Tuse和温度高于Tuse时功率的函数
		j = i + 1
		flag = true
		while flag && j < nT
			if TsList[j] <= Tuse - dT_EvaporationStandard
				C1, flag, P1, P3, Pe = powerCalculate_lowTs(TsList[j], TsList[i], dt)
				if flag
					C[i, j] = C1
					P1Matrix[i, j] = P1
					P3Matrix[i, j] = P3
					PeMatrix[i, j] = Pe
				end
				#C[i, j] = flag ? C1 : 99
			elseif TsList[i] < Tuse - dT_EvaporationStandard < TsList[j]
				dt1 = (Tuse - dT_EvaporationStandard - TsList[i]) / (TsList[j] - TsList[i]) * dt
				dt2 = dt - dt1
				C1, flag1, P11, P31, Pe1 = powerCalculate_lowTs(Tuse - dT_EvaporationStandard, TsList[i], dt1)
				C2, flag2, P12, P32, Pe2 = powerCalculate_highTs(TsList[j], Tuse - dT_EvaporationStandard, dt2)
				flag = flag1 && flag2
				if flag
					C[i, j] = C1 + C2
					P1Matrix[i, j] = (P11 * dt1 + P12 * dt2) / dt
					P3Matrix[i, j] = (P31 * dt1 + P32 * dt2) / dt
					PeMatrix[i, j] = (Pe1 * dt1 + Pe2 * dt2) / dt
				end
			elseif TsList[i] >= Tuse - dT_EvaporationStandard
				C1, flag, P1, P3, Pe = powerCalculate_highTs(TsList[j], TsList[i], dt)
				if flag
					C[i, j] = C1
					P1Matrix[i, j] = P1
					P3Matrix[i, j] = P3
					PeMatrix[i, j] = Pe
				end
			end
			j += 1
		end
		TsIncreaseIndexList[i] = flag ? j - 1 - i : j - 2 - i
	end

	return C, P1Matrix, P2Matrix, P3Matrix, PeMatrix, TsDecreaseIndexList, TsIncreaseIndexList
end

"""正向计算一个时间层上的动态规划"""
function forwardSolve(
	valueList::Vector,# 最优成本
	costMatrix::Matrix, # 状态转移成本矩阵
	TsDecreaseIndex::Int,# 温度下降最多偏移的index
	TsIncreaseIndex::Int,# 温度上升最多偏移的index
)
	nT = length(valueList)
	valueListNext = fill(99.0, nT)
	lastTsIndex = zeros(nT)
	for j ∈ 1:nT
		for i ∈ j-TsIncreaseIndex:j+TsDecreaseIndex
			if i < 1 || i > nT
				continue
			end
			if valueList[i] + costMatrix[i, j] < valueListNext[j]
				valueListNext[j] = valueList[i] + costMatrix[i, j]
				lastTsIndex[j] = i
			end
		end
	end
	return valueListNext, lastTsIndex
end

"""反向计算一个时间层上的动态规划"""
function backwardSolve(
	valueList::Vector,# 最优成本
	costMatrix::Matrix, # 状态转移成本矩阵
	TsDecreaseIndexList::Vector{Int},# 温度下降最多偏移的index
	TsIncreaseIndexList::Vector{Int},# 温度上升最多偏移的index
)
	nT = length(valueList)
	valueListNext = fill(99.0, nT)
	lastTsIndex = zeros(nT)
	for j ∈ 1:nT
		for i ∈ j-TsIncreaseIndexList[j]:j+TsDecreaseIndexList[j]
			if valueList[i] + costMatrix[i, j] < valueListNext[j]
				valueListNext[j] = valueList[i] + costMatrix[i, j]
				lastTsIndex[j] = i
			end
		end
	end
	return valueListNext, lastTsIndex
end

"""用动态规划完整求解一个初始温度下的最优运行策略"""
function dpSolve(
	nT::Int,
	nt::Int,
	C::Matrix{Float64},
	hourlyTariffFunction::Function,
	dt::Real,
	j::Int,
	TsDecreaseIndex::Int,# 温度下降最多偏移的index
	TsIncreaseIndex::Int,# 温度上升最多偏移的index
)
	TsTransitionMatrix = zeros(nT, nt - 1)
	# 第一步
	VForward = hourlyTariffFunction(0) * C[j, :]
	TsTransitionMatrix[:, 1] .= j
	for i in 2:nt-2
		VForward, TsTransitionMatrix[:, i] = forwardSolve(
			VForward,
			C * hourlyTariffFunction(dt * (i - 1)),
			TsDecreaseIndex,
			TsIncreaseIndex,
		)
	end

	valueMin = Inf # 在第j个温度下的最优成本
	tempC = C * hourlyTariffFunction(dt * (nt - 2))
	lastTsIndex = 0
	for i ∈ j-TsIncreaseIndex:j+TsDecreaseIndex
		if i < 1 || i > nT
			continue
		end
		if VForward[i] + tempC[i, j] < valueMin
			valueMin = VForward[i] + tempC[i, j]
			lastTsIndex = i
		end
	end

	# 状态回溯
	#valueMin=VForward[j] # 在第j个温度下的最优成本
	TsIndexList = Vector{Int}(undef, nt)
	TsIndexList[nt] = j
	TsIndexList[nt-1] = lastTsIndex
	for i in nt-2:-1:1
		TsIndexList[i] = TsTransitionMatrix[TsIndexList[i+1], i]
	end

	return valueMin, TsIndexList
end

"""给定系统参数,求解系统成本,返回成本、热泵功率向量、蓄热量向量"""
function generateAndSolve(::PressedWaterDoubleStorageOneCompressor, ::MinimizeCost, ::ConstloadandArea, ::ExhaustiveMethod;
	COPOverlap::Function,
	COPWater::Function,
	hourlyTariffFunction::Function,   # 电价函数
	heatConsumptionPowerFunction::Function,  # 用热负载函数
	TairFunction::Function,# 环境温度函数

	# 总循环参数
	Tuse::Real,# 供热蒸汽温度
	TCompressorIn::Real,# 中间级温度
	dT_EvaporationStandard::Real,#全蒸温差
	latentHeat::Real,# 汽化潜热
	cp_cw::Real,# 循环水定压热容
	cp_cs::Real,# 蒸汽定压热容 cp cycled steam
	TcChangeToElec::Real,
	TWaste::Real,# 废热回收蒸发器温度

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容

	# 设备运行约束
	TstorageTankMax::Real,# 蓄热罐的最高温度
	PheatPumpMax::Real,# 热泵最大功率
	PelecHeatMax::Real,# 电锅炉最大功率
	# 求解参数
	dT::Real = 0.1,# 状态参数高温蓄热温度离散步长
	dt::Real = 1 / 6,# 时间步长
)
	heatLoad = heatConsumptionPowerFunction(0.0)
	Tair = TairFunction(0.0)

	TsList = TCompressorIn+dT_EvaporationStandard:dT:TstorageTankMax
	tList = 0:dt:24
	nT = length(TsList)# 温度步数	
	nt = length(tList)# 时间步数
	TsMatrix = zeros(Int, nt, nT)# 存储状态参数：高温蓄热温度
	bestValueList = fill(99.0, nT)

	# 已经生成了C, TsDecreaseIndexList, TsIncreaseIndexList


	C, P1Matrix, P2Matrix, P3Matrix, PeMatrix, TsDecreaseIndexList, TsIncreaseIndexList = getStateTransitionCost(
		PressedWaterDoubleStorageOneCompressor();
		COPOverlap = COPOverlap,
		COPWater = COPWater,
		heatLoad = heatLoad,#热负荷
		Tair = Tair,# 外部环境温度
		# 总循环参数
		Tuse = Tuse,# 供热蒸汽温度
		dT_EvaporationStandard = dT_EvaporationStandard,#全蒸温差
		latentHeat = latentHeat,# 汽化潜热
		cp_cw = cp_cw,# 循环水定压热容
		TcChangeToElec = TcChangeToElec,
		TWaste = TWaste,# 废热回收蒸发器温度

		# 高温蓄热参数
		cpm_h = cpm_h,# 高温蓄热热容

		# 设备运行约束
		PheatPumpMax = PheatPumpMax,# 热泵最大功率
		PelecHeatMax = PelecHeatMax,# 电锅炉最大功率
		# 求解参数
		TsList = collect(TsList),# 状态参数高温蓄热温度列表
		dt = dt,# 时间步长
	)
	TsDecreaseIndex = maximum(TsDecreaseIndexList)
	TsIncreaseIndex = maximum(TsIncreaseIndexList)
	#最多降温TsDecreaseIndex,最多升温TsIncreaseIndex,所以从目标出发要多TsDecreaseIndex个，少TsIncreaseIndex个
	# 先直接正向计算
	# 正向计算步数
	count = 0
	Threads.@threads for j ∈ 1:nT
		# 生成状态记录矩阵
		valueMin, TsIndexList = dpSolve(
			nT,
			nt,
			C,
			hourlyTariffFunction,
			dt,
			j,
			TsDecreaseIndex,# 温度下降最多偏移的index
			TsIncreaseIndex,# 温度上升最多偏移的index
		)


		# 写入bestValueList与TsMatrix
		bestValueList[j] = valueMin
		TsMatrix[:, j] .= TsIndexList
		count += 1
		if count % 50 == 0
			println("$count/$nT")
		end
	end

	minCost, index = findmin(bestValueList)
	minTsList = [TsList[i] for i in TsMatrix[:, index]]

	P1List = [P1Matrix[TsMatrix[i, index], TsMatrix[i+1, index]] for i in 1:nt-1]
	P2List = [P2Matrix[TsMatrix[i, index], TsMatrix[i+1, index]] for i in 1:nt-1]
	P3List = [P3Matrix[TsMatrix[i, index], TsMatrix[i+1, index]] for i in 1:nt-1]
	PeList = [PeMatrix[TsMatrix[i, index], TsMatrix[i+1, index]] for i in 1:nt-1]

	return minCost, minTsList, P1List, P2List, P3List, PeList
end

function generateAndSolve(::PressedWaterDoubleStorageOneCompressor, ::MinimizeCost, ::ConstloadandArea, ::GoldenRatioMethod;
	COPOverlap::Function,
	COPWater::Function,
	hourlyTariffFunction::Function,   # 电价函数
	heatConsumptionPowerFunction::Function,  # 用热负载函数
	TairFunction::Function,# 环境温度函数

	# 总循环参数
	Tuse::Real,# 供热蒸汽温度
	TCompressorIn::Real,# 中间级温度
	dT_EvaporationStandard::Real,#全蒸温差
	latentHeat::Real,# 汽化潜热
	cp_cw::Real,# 循环水定压热容
	cp_cs::Real,# 蒸汽定压热容 cp cycled steam
	TcChangeToElec::Real,
	TWaste::Real,# 废热回收蒸发器温度

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容

	# 设备运行约束
	TstorageTankMax::Real,# 蓄热罐的最高温度
	PheatPumpMax::Real,# 热泵最大功率
	PelecHeatMax::Real,# 电锅炉最大功率
	# 求解参数
	dT::Real = 0.1,# 状态参数高温蓄热温度离散步长
	dt::Real = 1 / 6,# 时间步长
	K::Real = 12,# 温度列表离散分辨率倍数
)
	heatLoad = heatConsumptionPowerFunction(0.0)
	Tair = TairFunction(0.0)

	TsList = TCompressorIn+dT_EvaporationStandard:dT/K:TstorageTankMax
	tList = 0:dt:24
	nT = length(TsList)# 温度步数	
	nt = length(tList)# 时间步数
	TsMatrix = zeros(Int, nt, nT)# 存储状态参数：高温蓄热温度
	bestValueList = fill(99.0, nT)

	# 已经生成了C, TsDecreaseIndexList, TsIncreaseIndexList


	C, P1Matrix, P2Matrix, P3Matrix, PeMatrix, TsDecreaseIndexList, TsIncreaseIndexList = getStateTransitionCost(
		PressedWaterDoubleStorageOneCompressor();
		COPOverlap = COPOverlap,
		COPWater = COPWater,
		heatLoad = heatLoad,#热负荷
		Tair = Tair,# 外部环境温度
		# 总循环参数
		Tuse = Tuse,# 供热蒸汽温度
		dT_EvaporationStandard = dT_EvaporationStandard,#全蒸温差
		latentHeat = latentHeat,# 汽化潜热
		cp_cw = cp_cw,# 循环水定压热容
		TcChangeToElec = TcChangeToElec,
		TWaste = TWaste,# 废热回收蒸发器温度

		# 高温蓄热参数
		cpm_h = cpm_h,# 高温蓄热热容

		# 设备运行约束
		PheatPumpMax = PheatPumpMax,# 热泵最大功率
		PelecHeatMax = PelecHeatMax,# 电锅炉最大功率
		# 求解参数
		TsList = collect(TsList),# 状态参数高温蓄热温度列表
		dt = dt,# 时间步长
	)
	TsDecreaseIndex = maximum(TsDecreaseIndexList)
	TsIncreaseIndex = maximum(TsIncreaseIndexList)
	#最多降温TsDecreaseIndex,最多升温TsIncreaseIndex,所以从目标出发要多TsDecreaseIndex个，少TsIncreaseIndex个
	# 先直接正向计算
	# 正向计算步数
	count = 0
	# 初始化黄金分割法
	phi = 0.618
	jList = [1, nT - round(Int, phi * (nT - 1)), 1 + round(Int, phi * (nT - 1)), nT]
	valueList = zeros(4)
	TsIndexListMatirx = zeros(Int, nt, 4)
	for (i, j) in enumerate(jList)
		valueList[i], TsIndexListMatirx[:, i] = dpSolve(
			nT,
			nt,
			C,
			hourlyTariffFunction,
			dt,
			j,
			TsDecreaseIndex,# 温度下降最多偏移的index
			TsIncreaseIndex,# 温度上升最多偏移的index
		)
	end
	#count = 0
	while (jList[4] - jList[1] >= 4) && count < 100
		if valueList[2] < valueList[3]# 最值在j1到j3之间
			jList[4] = jList[3]
			valueList[4] = valueList[3]
			TsIndexListMatirx[:, 4] = TsIndexListMatirx[:, 3]
			jList[3] = jList[2]
			valueList[3] = valueList[2]
			TsIndexListMatirx[:, 3] = TsIndexListMatirx[:, 2]
			jList[2] = min(max(floor(Int, jList[4] - phi * (jList[4] - jList[1])), jList[1] + 1), jList[3] - 1)
			valueList[2], TsIndexListMatirx[:, 2] = dpSolve(
				nT,
				nt,
				C,
				hourlyTariffFunction,
				dt,
				jList[2],
				TsDecreaseIndex,# 温度下降最多偏移的index
				TsIncreaseIndex,# 温度上升最多偏移的index
			)
		elseif valueList[2] >= valueList[3]
			jList[1] = jList[2]
			valueList[1] = valueList[2]
			TsIndexListMatirx[:, 1] = TsIndexListMatirx[:, 2]
			jList[2] = jList[3]
			valueList[2] = valueList[3]
			TsIndexListMatirx[:, 2] = TsIndexListMatirx[:, 3]
			jList[3] = max(min(ceil(Int, jList[1] + phi * (jList[4] - jList[1])), jList[4] - 1), jList[2] + 1)
			valueList[3], TsIndexListMatirx[:, 3] = dpSolve(
				nT,
				nt,
				C,
				hourlyTariffFunction,
				dt,
				jList[3],
				TsDecreaseIndex,# 温度下降最多偏移的index
				TsIncreaseIndex,# 温度上升最多偏移的index
			)
		end
		#count += 1
		#println("$count/$nT", "j1=$(jList[1]),j2=$(jList[2]),j3=$(jList[3]),j4=$(jList[4])")
	end

	minCost, index = findmin(valueList)
	minTsList = TsIndexListMatirx[:, index]

	P1List = [P1Matrix[minTsList[i], minTsList[i+1]] for i in 1:nt-1]
	P2List = [P2Matrix[minTsList[i], minTsList[i+1]] for i in 1:nt-1]
	P3List = [P3Matrix[minTsList[i], minTsList[i+1]] for i in 1:nt-1]
	PeList = [PeMatrix[minTsList[i], minTsList[i+1]] for i in 1:nt-1]

	return minCost, TsList[minTsList], P1List, P2List, P3List, PeList
end

function generateAndSolve(::PressedWaterDoubleStorageOneCompressor, ::MinimizeCost, ::NoSimplify, ::GoldenRatioMethod;
	COPOverlap::Function,
	COPWater::Function,
	hourlyTariffFunction::Function,   # 电价函数
	heatConsumptionPowerFunction::Function,  # 用热负载函数
	TairFunction::Function,# 环境温度函数

	# 总循环参数
	Tuse::Real,# 供热蒸汽温度
	TCompressorIn::Real,# 中间级温度
	dT_EvaporationStandard::Real,#全蒸温差
	latentHeat::Real,# 汽化潜热
	cp_cw::Real,# 循环水定压热容
	cp_cs::Real,# 蒸汽定压热容 cp cycled steam
	TcChangeToElec::Real,
	TWaste::Real,# 废热回收蒸发器温度

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容

	# 设备运行约束
	TstorageTankMax::Real,# 蓄热罐的最高温度
	PheatPumpMax::Real,# 热泵最大功率
	PelecHeatMax::Real,# 电锅炉最大功率
	# 求解参数
	dT::Real = 0.1,# 状态参数高温蓄热温度离散步长
	dt::Real = 1 / 6,# 时间步长
)
	heatLoad = heatConsumptionPowerFunction(0.0)
	Tair = TairFunction(0.0)

	TsList = TCompressorIn+dT_EvaporationStandard:dT:TstorageTankMax
	tList = 0:dt:24
	nT = length(TsList)# 温度步数	
	nt = length(tList)# 时间步数
	TsMatrix = zeros(Int, nt, nT)# 存储状态参数：高温蓄热温度
	bestValueList = fill(99.0, nT)

	# 已经生成了C, TsDecreaseIndexList, TsIncreaseIndexList


	C, P1Matrix, P2Matrix, P3Matrix, PeMatrix, TsDecreaseIndexList, TsIncreaseIndexList = getStateTransitionCost(
		PressedWaterDoubleStorageOneCompressor();
		COPOverlap = COPOverlap,
		COPWater = COPWater,
		heatLoad = heatLoad,#热负荷
		Tair = Tair,# 外部环境温度
		# 总循环参数
		Tuse = Tuse,# 供热蒸汽温度
		dT_EvaporationStandard = dT_EvaporationStandard,#全蒸温差
		latentHeat = latentHeat,# 汽化潜热
		cp_cw = cp_cw,# 循环水定压热容
		TcChangeToElec = TcChangeToElec,
		TWaste = TWaste,# 废热回收蒸发器温度

		# 高温蓄热参数
		cpm_h = cpm_h,# 高温蓄热热容

		# 设备运行约束
		PheatPumpMax = PheatPumpMax,# 热泵最大功率
		PelecHeatMax = PelecHeatMax,# 电锅炉最大功率
		# 求解参数
		TsList = collect(TsList),# 状态参数高温蓄热温度列表
		dt = dt,# 时间步长
	)
	TsDecreaseIndex = maximum(TsDecreaseIndexList)
	TsIncreaseIndex = maximum(TsIncreaseIndexList)
	#最多降温TsDecreaseIndex,最多升温TsIncreaseIndex,所以从目标出发要多TsDecreaseIndex个，少TsIncreaseIndex个
	# 先直接正向计算
	# 正向计算步数
	count = 0
	# 初始化黄金分割法
	phi = 0.618
	jList = [1, nT - round(Int, phi * (nT - 1)), 1 + round(Int, phi * (nT - 1)), nT]
	valueList = zeros(4)
	TsIndexListMatirx = zeros(Int, nt, 4)
	for (i, j) in enumerate(jList)
		valueList[i], TsIndexListMatirx[:, i] = dpSolve(
			nT,
			nt,
			C,
			hourlyTariffFunction,
			dt,
			j,
			TsDecreaseIndex,# 温度下降最多偏移的index
			TsIncreaseIndex,# 温度上升最多偏移的index
		)
	end
	count = 0
	while (jList[4] - jList[1] >= 4) && count < 100
		if valueList[2] < valueList[3]# 最值在j1到j3之间
			jList[4] = jList[3]
			valueList[4] = valueList[3]
			TsIndexListMatirx[:, 4] = TsIndexListMatirx[:, 3]
			jList[3] = jList[2]
			valueList[3] = valueList[2]
			TsIndexListMatirx[:, 3] = TsIndexListMatirx[:, 2]
			jList[2] = min(max(floor(Int, jList[4] - phi * (jList[4] - jList[1])), jList[1] + 1), jList[3] - 1)
			valueList[2], TsIndexListMatirx[:, 2] = dpSolve(
				nT,
				nt,
				C,
				hourlyTariffFunction,
				dt,
				jList[2],
				TsDecreaseIndex,# 温度下降最多偏移的index
				TsIncreaseIndex,# 温度上升最多偏移的index
			)
		elseif valueList[2] >= valueList[3]
			jList[1] = jList[2]
			valueList[1] = valueList[2]
			TsIndexListMatirx[:, 1] = TsIndexListMatirx[:, 2]
			jList[2] = jList[3]
			valueList[2] = valueList[3]
			TsIndexListMatirx[:, 2] = TsIndexListMatirx[:, 3]
			jList[3] = max(min(ceil(Int, jList[1] + phi * (jList[4] - jList[1])), jList[4] - 1), jList[2] + 1)
			valueList[3], TsIndexListMatirx[:, 3] = dpSolve(
				nT,
				nt,
				C,
				hourlyTariffFunction,
				dt,
				jList[3],
				TsDecreaseIndex,# 温度下降最多偏移的index
				TsIncreaseIndex,# 温度上升最多偏移的index
			)
		end
		count += 1
		println("$count/$nT", "j1=$(jList[1]),j2=$(jList[2]),j3=$(jList[3]),j4=$(jList[4])")
	end

	minCost, index = findmin(valueList)
	minTsList = TsIndexListMatirx[:, index]

	P1List = [P1Matrix[minTsList[i], minTsList[i+1]] for i in 1:nt-1]
	P2List = [P2Matrix[minTsList[i], minTsList[i+1]] for i in 1:nt-1]
	P3List = [P3Matrix[minTsList[i], minTsList[i+1]] for i in 1:nt-1]
	PeList = [PeMatrix[minTsList[i], minTsList[i+1]] for i in 1:nt-1]

	return minCost, TsList[minTsList], P1List, P2List, P3List, PeList
end

function testSolve(::PressedWaterDoubleStorageOneCompressor, ::MinimizeCost, ::ConstloadandArea, ::ExhaustiveMethod;
	COPOverlap::Function,
	COPWater::Function,
	hourlyTariffFunction::Function,   # 电价函数
	heatConsumptionPowerFunction::Function,  # 用热负载函数
	TairFunction::Function,# 环境温度函数

	# 总循环参数
	Tuse::Real,# 供热蒸汽温度
	TCompressorIn::Real,# 中间级温度
	dT_EvaporationStandard::Real,#全蒸温差
	latentHeat::Real,# 汽化潜热
	cp_cw::Real,# 循环水定压热容
	cp_cs::Real,# 蒸汽定压热容 cp cycled steam
	TcChangeToElec::Real,
	TWaste::Real,# 废热回收蒸发器温度

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容

	# 设备运行约束
	TstorageTankMax::Real,# 蓄热罐的最高温度
	PheatPumpMax::Real,# 热泵最大功率
	PelecHeatMax::Real,# 电锅炉最大功率
	# 求解参数
	dT::Real = 0.1,# 状态参数高温蓄热温度离散步长
	dt::Real = 1 / 6,# 时间步长
	K::Real = 12,# 温度列表离散分辨率倍数
)
	heatLoad = heatConsumptionPowerFunction(0.0)
	Tair = TairFunction(0.0)

	TsList = TCompressorIn+dT_EvaporationStandard:dT/K:TstorageTankMax
	tList = 0:dt:24
	nT = length(TsList)# 温度步数	
	nt = length(tList)# 时间步数
	TsMatrix = zeros(Int, nt, nT)# 存储状态参数：高温蓄热温度
	bestValueList = fill(99.0, nT)

	# 已经生成了C, TsDecreaseIndexList, TsIncreaseIndexList


	C, P1Matrix, P2Matrix, P3Matrix, PeMatrix, TsDecreaseIndexList, TsIncreaseIndexList = getStateTransitionCost(
		PressedWaterDoubleStorageOneCompressor();
		COPOverlap = COPOverlap,
		COPWater = COPWater,
		heatLoad = heatLoad,#热负荷
		Tair = Tair,# 外部环境温度
		# 总循环参数
		Tuse = Tuse,# 供热蒸汽温度
		dT_EvaporationStandard = dT_EvaporationStandard,#全蒸温差
		latentHeat = latentHeat,# 汽化潜热
		cp_cw = cp_cw,# 循环水定压热容
		TcChangeToElec = TcChangeToElec,
		TWaste = TWaste,# 废热回收蒸发器温度

		# 高温蓄热参数
		cpm_h = cpm_h,# 高温蓄热热容

		# 设备运行约束
		PheatPumpMax = PheatPumpMax,# 热泵最大功率
		PelecHeatMax = PelecHeatMax,# 电锅炉最大功率
		# 求解参数
		TsList = collect(TsList),# 状态参数高温蓄热温度列表
		dt = dt,# 时间步长
	)
	TsDecreaseIndex = maximum(TsDecreaseIndexList)
	TsIncreaseIndex = maximum(TsIncreaseIndexList)
	#最多降温TsDecreaseIndex,最多升温TsIncreaseIndex,所以从目标出发要多TsDecreaseIndex个，少TsIncreaseIndex个
	# 先直接正向计算
	# 正向计算步数
	count = 0
	Threads.@threads for j ∈ 1:nT
		# 生成状态记录矩阵
		valueMin, TsIndexList = dpSolve(
			nT,
			nt,
			C,
			hourlyTariffFunction,
			dt,
			j,
			TsDecreaseIndex,# 温度下降最多偏移的index
			TsIncreaseIndex,# 温度上升最多偏移的index
		)


		# 写入bestValueList与TsMatrix
		bestValueList[j] = valueMin
		TsMatrix[:, j] .= TsIndexList
		count += 1
		println("$count/$nT")
	end

	minCost, index = findmin(bestValueList)
	minTsList = [TsList[i] for i in TsMatrix[:, index]]

	P1List = [P1Matrix[TsMatrix[i, index], TsMatrix[i+1, index]] for i in 1:nt-1]
	P2List = [P2Matrix[TsMatrix[i, index], TsMatrix[i+1, index]] for i in 1:nt-1]
	P3List = [P3Matrix[TsMatrix[i, index], TsMatrix[i+1, index]] for i in 1:nt-1]
	PeList = [PeMatrix[TsMatrix[i, index], TsMatrix[i+1, index]] for i in 1:nt-1]

	return bestValueList,TsMatrix,minCost, minTsList, P1List, P2List, P3List, PeList
end
