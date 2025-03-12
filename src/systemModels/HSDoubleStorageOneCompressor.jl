
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
		list = [list[1],list[1]]
	end
	n=2
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
	overlapRefrigerant::OverlapRefrigerant = NH3_Water,	# 复叠工质
	maxTcHigh::Real = 180.0,  # 高温热泵冷凝器温度上限
	TCompressorIn::Real = 115.0,# 中间温度
	TWaste::Real = 60.0,                  # 废热源温度
	Tair::Real = 25.0,          # 外部环境温度
	Tuse::Real = 120.0,                     # 工厂使用温度
	Trecycle::Real = 115.0,    # 回收冷凝水温度
	heatStorageCapacity::Real = 6.0,        # 蓄热量kWh(承压水蓄热)
	TstorageTankMax::Real = 220.0,          # 蓄热罐的最高温度
	maxheatStorageInputHour::Real = 4,      # 蓄热充满时长
	dT_EvaporationStandard::Real = 5.0,		# 全蒸温差
	heatConsumptionPower::Real = 1.0,       # 每小时用热功率kW
	workingStartHour::Int = 8,              # 生产开始时间
	workingHours::Int = 16,                 # 每日工作小时数
	heatPumpServiceCoff::Real = 1.2,        # 热泵功率服务系数
	elecHeatServiceCoff::Real = 1.5,		# 电锅炉服务系数
	hourlyTariff::Vector = fill(0.7, 24),   # 电价向量
)
	# 首先生产COP函数
	COPOverlap=getOverlapCOP_fixMidTemperature(
		overlapRefrigerant,
		TCompressorIn+overlapRefrigerant.midTDifference/2;
		maxCOP=21,# 最大COP
		eta_s=0.7,# 绝热效率
		dT=0.1,# 插值步长
	)
	
	COPWater=getCOP(
		Trecycle,# 蒸发温度下限,这里是实际设计中的蒸发冷凝温度界限
		maxTcHigh,# 蒸发温度上限
		Trecycle,# 冷凝温度下限
		maxTcHigh,# 冷凝温度上限
		refWater,# 工质
		21,# 最大COP
		0.7,# 绝热效率
		0.1,# 插值步长
	)

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
	P1base=heatConsumptionPower/COPOverlap(TWaste, Tuse)
	PheatPumpMax=P1base*heatPumpServiceCoff
	PelecHeatMax=heatConsumptionPower*elecHeatServiceCoff

	# 生成需求与环境函数
	hourlyTariffFunction=generateGridPriceFunction(hourlyTariff,24)
	heatConsumptionPowerFunction=generateLoadFunction([heatConsumptionPower], 24)
	TairFunction=generateAreaTemperatureFunction([Tair], 24)

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

	return COPOverlap,COPWater,
	hourlyTariffFunction, heatConsumptionPowerFunction, TairFunction,
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

"""给定系统参数,求解系统成本,返回成本、热泵功率向量、蓄热量向量"""
function generateAndSolve(::PressedWaterDoubleStorageOneCompressor, ::MinimizeCost,::ConstloadandArea;
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
	COP1 = COPOverlap(TWaste, Tuse)
	COP2 = COPWater
	COP3 = COPOverlap


	heatLoad = heatConsumptionPowerFunction(0.0)
	Tair = TairFunction(0.0)

	TsList = TCompressorIn+dT_EvaporationStandard:dT:TstorageTankMax
	tList = 0:dt:24
	nT = length(TsList)# 温度步数	
	nt = length(tList)# 时间步数
	TsMatrix = zeros(nT, nt)# 存储状态参数：高温蓄热温度

	#C[i,j]表示温度从TsList[i]到TsList[j]时的最低功率；在负载、环境不变的情况下，C[i,j]是不变的
	C = fill(99.0, nT, nT)# 状态转移矩阵
	TsDecreaseIndexList = zeros(Int,nT)# 记录温度下降最多偏移的index
	TsIncreaseIndexList = zeros(Int,nT)# 记录温度上升最多偏移的index

	# 只用热泵供热时的功率
	P1Only = heatLoad / COP1# 只用热泵供热时的功率

	"""计算蓄热温度相对较低时的电度情况"""
	function powerCalculate_lowTs(Tsaim, Tsstart, dt)
		P3 = min(PheatPumpMax - P1Only, cpm_h / dt / COP1 * (Tsaim - Tsstart))
		Pe = cpm_h / dt * (Tsaim - Tsstart) - P3 * COP1
		flag = Pe <= PelecHeatMax
		return (P1Only + P3 + Pe) * dt, flag
	end

	"""计算蓄热温度相对较高时的电度情况"""
	function powerCalculate_highTs(Tsaim, Tsstart, dt)
		COP3value = COP3(TWaste, (Tsaim + Tsstart) / 2 + dT_EvaporationStandard)
		Ptotal = 999.0
		flag = false
		# 模式1：电热加热
		Pe = cpm_h / dt * (Tsaim - Tsstart)
		if Pe <= PelecHeatMax
			flag = true
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
		if (Pe <= PelecHeatMax) && (Tsaim <= TcChangeToElec - dT_EvaporationStandard)
			flag = true
			Ptotal = min(Ptotal, P1 + P3 + Pe)
		end
		return Ptotal * dt, flag
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
		TsNextMinGrid = max(ceil(TsNextMin, digits = 1), TsList[1])
		TsDecreaseIndexList[i] = round(Int, (TsList[i] - TsNextMinGrid) * 10)
		elecP = (TsNextMinGrid - TsNextMin) * cpm_h

		index = i - TsDecreaseIndexList[i]
		for j ∈ index:i
			P2 = heatLoad / (COP2((TsList[i] + TsList[j]) / 2 - dT_EvaporationStandard, Tuse))
			elecP = (TsList[j] - TsNextMin) * cpm_h
			#补热的功率基本上超不了，因为这个工况下从蓄热取走的热量小于工厂的负载，而补热的功率大于等于工厂的负载
			if elecP < PelecHeatMax * dt
				C[i, j] = P2 * dt + elecP
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
				C1, flag = powerCalculate_lowTs(TsList[j], TsList[i], dt)
				if flag C[i, j]=C1 end
				#C[i, j] = flag ? C1 : 99
			elseif TsList[i] < Tuse - dT_EvaporationStandard < TsList[j]
				dt1 = (Tuse - dT_EvaporationStandard - TsList[i]) / (TsList[j] - TsList[i]) * dt
				dt2 = dt - dt1
				C1, flag1 = powerCalculate_lowTs(Tuse - dT_EvaporationStandard, TsList[i], dt1)
				C2, flag2 = powerCalculate_highTs(TsList[j], Tuse - dT_EvaporationStandard, dt2)
				flag = flag1 && flag2
				if flag C[i, j]=C1+C2 end
			elseif TsList[i] >= Tuse - dT_EvaporationStandard
				C1, flag = powerCalculate_highTs(TsList[j], TsList[i], dt)
				if flag C[i, j]=C1 end
			end
			j += 1
		end
		TsIncreaseIndexList[i] = flag ? j - 1 - i : j - 2 - i
	end
	
	# 已经生成了C, TsDecreaseIndexList, TsIncreaseIndexList
	TsDecreaseIndex=maximum(TsDecreaseIndexList)
	TsIncreaseIndex=maximum(TsIncreaseIndexList)
	#最多降温TsDecreaseIndex,最多升温TsIncreaseIndex,所以从目标出发要多TsDecreaseIndex个，少TsIncreaseIndex个
	# 先直接正向计算
	# 正向计算步数
	nForward = nt//2
	
	C, TsDecreaseIndexList, TsIncreaseIndexList
end



