#=
变工况与定工况的区别为：热泵配置可能小于标准配置，原先热泵的配置大小足够，蓄热在放热时因为效率总是高于标准工况所以不用考虑需要电热补热；现在需要考虑电热补热。因此：
1.需要考虑热泵供热不足时用电热补热
2.需要给定压缩机的运行功率约束
=#



"""
单个时间层的状态转移成本函数
"""
function getStateTransitionCost_SingleStep(
	::PressedWaterOneStorageOneCompressor;
	COPOverlap::Function,
	Tair::Real,# 外部环境温度
	# 总循环参数

	# 求解参数
	TsListStart::Vector,# 状态参数高温蓄热温度起始值列表
	TsListEnd::Vector,# 状态参数高温蓄热温度结束值列表
	params::SystemParameters,
	sysVariables::SystemVariables,
	dt::Real = 1.0,# 时间步长
)
	Tsmin = params.Tsmin
	Tsmax = params.Tsmax
	#heatLoadPumpMax = PWaterCompressorMax * COP2_design# 热泵设计的供热功率

	nT = length(TsListStart)# 温度步数

	#C[i,j]表示温度从TsList[i]到TsList[j]时的最低功率；在负载、环境不变的情况下，C[i,j]是不变的
	C = fill(9999.0, nT, nT)# 状态转移矩阵
	P1Matrix = zeros(nT, nT)# 状态转移功率1参数
	P2Matrix = zeros(nT, nT)# 状态转移功率2参数
	P3Matrix = zeros(nT, nT)# 状态转移功率3参数
	PeMatrix = zeros(nT, nT)# 状态转移功率电加热参数


	# 只用热泵供热时的功率
	#P1Only = min(heatLoad, heatLoadPumpMax) / COP1# 只用热泵供热时的功率
	#P1hOnly = P1Only * COP1 / COP2_design#只用热泵供热时的压缩机功率



	for (j, Tsaim) in enumerate(TsListEnd)
		#println(j,"/",length(TsListEnd))
		midIndex = findfirst(x -> x >= Tsaim, TsListStart)
		# 大于midIndex的是温度不增部分，小于midIndex的是温度下降部分
		if isnothing(midIndex)
			midIndex = length(TsListStart) + 1
		end

		# 温度不增部分，出现不可行解时跳过
		for i ∈ midIndex:length(TsListStart)
			Tsstart = TsListStart[i]
			if Tsstart < Tsmin || Tsaim > Tsmax
				C[i, j] = 9999.0
				P1Matrix[i, j] = 9999.0
				P3Matrix[i, j] = 9999.0
				PeMatrix[i, j] = 9999.0
				continue
			end
			C[i, j], flag, P1Matrix[i, j], P2Matrix[i, j], P3Matrix[i, j], PeMatrix[i, j] = getMinimumCost(Tsstart, Tsaim, dt, params, sysVariables)
			if !flag
				C[i, j] = 9999.0
				P1Matrix[i, j] = 9999.0
				P3Matrix[i, j] = 9999.0
				PeMatrix[i, j] = 9999.0
				break
			end
		end
		for i ∈ midIndex-1:-1:1
			Tsstart = TsListStart[i]
			if Tsstart < Tsmin || Tsaim > Tsmax
				C[i, j] = 9999.0
				P1Matrix[i, j] = 9999.0
				P3Matrix[i, j] = 9999.0
				PeMatrix[i, j] = 9999.0
				continue
			end
			C[i, j], flag, P1Matrix[i, j], P2Matrix[i, j], P3Matrix[i, j], PeMatrix[i, j] = getMinimumCost(Tsstart, Tsaim, dt, params, sysVariables)
			if !flag
				C[i, j] = 9999.0
				P1Matrix[i, j] = 9999.0
				P3Matrix[i, j] = 9999.0
				PeMatrix[i, j] = 9999.0
				break
			end
		end
	end

	return C, P1Matrix, P2Matrix, P3Matrix, PeMatrix
end

function getStateTransitionCost(::PressedWaterOneStorageOneCompressor, ::VaryLoadVaryArea;
	COPOverlap::Function,
	COPLowFunction::Function,
	heatLoad::Vector,#热负荷
	Tair::Vector,# 外部环境温度
	costGrid::Vector,# 电网电价
	# 总循环参数
	TWaste::Real,# 废热回收蒸发器温度
	params::SystemParameters = params,

	# 求解参数
	TsListStart::Matrix,# 状态参数高温蓄热温度起始值列表
	TsListEnd::Matrix,# 状态参数高温蓄热温度结束值列表
	dt::Real = 1.0,# 时间步长
	smoother::Real = 1e-8,# 状态转移平滑系数
)
	dT_EvaporationStandard = params.dT
	TCompressorIn = params.TCompressorIn

	nt = 24.0 / dt |> Int#环节数
	nT = size(TsListStart, 1)

	C_smoothed = zeros(nt, nT, nT)
	P1Matrix = zeros(nt, nT, nT)
	P2Matrix = zeros(nt, nT, nT)
	P3Matrix = zeros(nt, nT, nT)
	PeMatrix = zeros(nt, nT, nT)

	for i ∈ 1:nt
		sysVariables = SystemVariables(
			heatLoad[i],
			COPLowFunction(TWaste, TCompressorIn + dT_EvaporationStandard),
			Tair[i],
			TWaste,
		)
		# 计算每个时段内的状态转移矩阵
		C_singlestep, P1Matrix[i, :, :], P2Matrix[i, :, :], P3Matrix[i, :, :], PeMatrix[i, :, :] = getStateTransitionCost_SingleStep(
			PressedWaterOneStorageOneCompressor();
			COPOverlap = COPOverlap,
			Tair = Tair[i],
			TsListStart = TsListStart[:, i],
			TsListEnd = TsListEnd[:, i],
			params = params,
			sysVariables = sysVariables,
			dt = dt,
		)
		C_smoothed[i, :, :] = C_singlestep * costGrid[i] + smoother * (P1Matrix[i, :, :] .^ 2 + P2Matrix[i, :, :] .^ 2 + P3Matrix[i, :, :] .^ 2 + PeMatrix[i, :, :] .^ 2)
	end

	return C_smoothed, P1Matrix, P2Matrix, P3Matrix, PeMatrix
end

function generateAndSolve(::PressedWaterOneStorageOneCompressor, ::MinimizeCost, ::VaryLoadVaryArea, ::GoldenRatioMethod;
	COPOverlap::Function,
	COPLowFunction::Function,
	hourlyTariffFunction::Function,   # 电价函数
	heatConsumptionPowerFunction::Function,  # 用热负载函数
	TairFunction::Function,# 环境温度函数

	# 总循环参数
	params::SystemParameters,
	TWaste::Real,# 废热回收蒸发器温度

	# 求解参数
	dT::Real = 0.01,# 状态参数高温蓄热温度离散步长
	dt::Real = 1 / 6,# 时间步长
	smoother::Real = 0.01,
)
	Tuse = params.Tuse
	dT_EvaporationStandard = params.dT
	TCompressorIn = params.TCompressorIn
	cpm_h = params.cpm
	Tsmin = params.Tsmin
	Tsmax = params.Tsmax
	TstorageTankMax = params.Tsmax
	TcChangeToElec = params.ThMax

	#Tuse::Real,# 供热蒸汽温度
	if cpm_h == 0
		#return sum(realCostList), TsList, P1List, P2List, P3List, PeList,realCostList
		nt = Int(24 / dt + 1)
		tList = 0:dt:24
		TsList = fill(Tsmin, nt)
		P1List = zeros(nt - 1)
		PeList = zeros(nt - 1)
		realCostList = zeros(nt - 1)

		heatLoadList = heatConsumptionPowerFunction.(tList)
		TairList = TairFunction.(tList)
		costGridList = hourlyTariffFunction.(tList)

		for i ∈ 1:nt-1
			sysVariables = SystemVariables(
				heatLoadList[i],
				COPLowFunction(TWaste, TCompressorIn + dT_EvaporationStandard),
				TairList[i],
				TWaste,
			)
			cost_test, _, P1List[i], _, _, PeList[i] = getMinimumCost(TsList[i], TsList[i+1], dt, params, sysVariables)

			realCostList[i] = cost_test * costGridList[i]
		end
		return sum(realCostList), TsList, P1List, zeros(nt - 1), zeros(nt - 1), PeList, realCostList
	end

	# 先试算，温差除以时间要小于一个数，默认是10℃/1h=10
	k_dT_to_dt = 10
	#=
	# 生成初始解
	begin
		dt_test = 2 #2小时试算
		tList = 0:dt_test:24
		nt=length(tList)

		dT_origin = dt_test * 5 #试算温度步长

		TsMatrix = repeat(TCompressorIn+dT_EvaporationStandard:dT_origin:TstorageTankMax, 1, length(tList))

		nT = size(TsMatrix, 1)# 温度步数


		heatLoadList = heatConsumptionPowerFunction.(tList)
		TairList = TairFunction.(tList)
		costGridList = hourlyTariffFunction.(tList)

		# 生成初始解
		## 状态转移矩阵的计算
		C, P1Matrix, P2Matrix, P3Matrix, PeMatrix = getStateTransitionCost(
			PressedWaterOneStorageOneCompressor(),
			VaryLoadVaryArea();
			COPOverlap = COPOverlap,
			COPLowFunction = COPLowFunction,
			heatLoad = heatLoadList,# 热负荷
			Tair = TairList,# 外部环境温度
			costGrid = costGridList,# 电网电价
			# 总循环参数
			TWaste = TWaste,# 废热回收蒸发器温度
			params = params,
			# 求解参数
			TsListStart = TsMatrix[:, 1:end-1], # 状态参数高温蓄热温度列表
			TsListEnd = TsMatrix[:, 2:end], # 状态参数高温蓄热温度列表
			dt = dt_test, # 时间步长
			smoother = smoother,
		)

		## 动态规划求解
		cost, TsIndex = GoldenRatioSolver(nT, nt, C)
		TsList_test = map(i -> TsMatrix[TsIndex[i], i], 1:nt)

		# 验证测试首尾是否一致
		if TsList_test[1] != TsList_test[end]
			println("""
			试算首尾Ts不一致
			TsList_test[1] = $(TsList_test[1])
			TsList_test[end] = $(TsList_test[end])
			dT_origin = $dT_origin
			""")
		else
			#println("试算首尾Ts一致")
		end

		# 把TsList从粗时间网格上扩充到细时间网格上。1小时一定是dt的整倍数
		kt = dt_test / dt |> Int
		TsList = zeros((nt - 1) * kt + 1)

		for i ∈ 1:nt-1
			for j ∈ 1:kt
				TsList[(i-1)*kt+j] = TsList_test[i] + (j - 1) / kt * (TsList_test[i+1] - TsList_test[i])
			end
		end
		TsList[end] = TsList_test[end]
	end
	=#
	
	dT_origin = dt * k_dT_to_dt#正式计算的初始温度步长

	# 开始改良解
	tList = 0:dt:24# 正式计算的时间步

	nTList = [3,5]
	changedT = 64*dT
	is_nt_changed = false
	nT = nTList[1]# 温度步数
	half_nT = Int((nT - 1) / 2)
	nt = length(tList)
	#TsList = fill(Tsmax, nt)
	TsMatrix = zeros(nT, nt)

	heatLoadList = heatConsumptionPowerFunction.(tList)
	TairList = TairFunction.(tList)
	costGridList = hourlyTariffFunction.(tList)
	
	countAll = 0
	countSingleGap = 0
	maxcount = 500
	#=
	C = zeros(nt, nT, nT)
	P1Matrix = zeros(nt, nT, nT)
	P2Matrix = zeros(nt, nT, nT)
	P3Matrix = zeros(nt, nT, nT)
	PeMatrix = zeros(nt, nT, nT)
	TsIndex = zeros(nt)
	=#

	# 初始化初值，这一步应该调用一个模型做预测
	TsList = fill(TcChangeToElec+5.0, nt)
	for j ∈ 1:nt
		TsMatrix[:, j] = TsList[j]-half_nT*dT_origin:dT_origin:TsList[j]+half_nT*dT_origin
	end
	
	

	C, P1Matrix, P2Matrix, P3Matrix, PeMatrix = getStateTransitionCost(
		PressedWaterOneStorageOneCompressor(),
		VaryLoadVaryArea();
		COPOverlap = COPOverlap,
		COPLowFunction = COPLowFunction,
		heatLoad = heatLoadList,# 热负荷
		Tair = TairList,# 外部环境温度
		costGrid = costGridList,# 电网电价
		# 总循环参数
		TWaste = TWaste,# 废热回收蒸发器温度
		params = params,
		# 求解参数
		TsListStart = TsMatrix[:, 1:end-1],# 状态参数高温蓄热温度列表
		TsListEnd = TsMatrix[:, 2:end],# 状态参数高温蓄热温度列表
		dt = dt,# 时间步长
		smoother = smoother,
	)
	## 动态规划求解
	cost, TsIndex = ExhaustiveSolver(nT, nt, C)
	TsList = map(i -> TsMatrix[TsIndex[i], i], 1:nt)
	# 初值有问题则采用最原始的初值
	if cost > 1000
		println("初值有问题")
		TsList = fill(TcChangeToElec+5.0, nt)
		for j ∈ 1:nt
			TsMatrix[:, j] = TsList[j]-half_nT*dT_origin:dT_origin:TsList[j]+half_nT*dT_origin
		end
		
		

		C, P1Matrix, P2Matrix, P3Matrix, PeMatrix = getStateTransitionCost(
			PressedWaterOneStorageOneCompressor(),
			VaryLoadVaryArea();
			COPOverlap = COPOverlap,
			COPLowFunction = COPLowFunction,
			heatLoad = heatLoadList,# 热负荷
			Tair = TairList,# 外部环境温度
			costGrid = costGridList,# 电网电价
			# 总循环参数
			TWaste = TWaste,# 废热回收蒸发器温度
			params = params,
			# 求解参数
			TsListStart = TsMatrix[:, 1:end-1],# 状态参数高温蓄热温度列表
			TsListEnd = TsMatrix[:, 2:end],# 状态参数高温蓄热温度列表
			dt = dt,# 时间步长
			smoother = smoother,
		)
		## 动态规划求解
		cost, TsIndex = ExhaustiveSolver(nT, nt, C)
		TsList = map(i -> TsMatrix[TsIndex[i], i], 1:nt)
	end

	while dT_origin > dT && countAll < maxcount
		flag_nextgap = true
		# 精度足够后尝试改变温度步数
		if dT_origin <= changedT && !is_nt_changed
			nT = nTList[2]
			half_nT = Int((nT - 1) / 2)
			TsMatrix = zeros(nT, nt)
			TsIndex = fill(half_nT+1,nt)
			C = zeros(nt, nT, nT)
			P1Matrix = zeros(nt, nT, nT)
			P2Matrix = zeros(nt, nT, nT)
			P3Matrix = zeros(nt, nT, nT)
			PeMatrix = zeros(nt, nT, nT)
			is_nt_changed = true
		end
		
		#=
		for j ∈ 1:nt #范围调整
			if TsIndex[j] == 1 || TsIndex[j] == nT # 温度范围要更小
				TsMatrix[:, j] = TsList[j]-half_nT*dT_origin:dT_origin:(TsList[j]+half_nT*dT_origin+1e-8)
				flag_nextgap = false
			end
		end
		=#
		TsIndexChange = TsIndex.-half_nT# 温度索引的偏移量
		for j ∈ 1:nt #范围调整
			if TsIndex[j] != half_nT # 温度范围要更小
				TsMatrix[:, j] = TsList[j]-half_nT*dT_origin:dT_origin:(TsList[j]+half_nT*dT_origin+1e-8)
				flag_nextgap = false
			end
		end
		# 根据温度索引的偏移量更新状态转移成本矩阵
		for j ∈ 1:nt #范围调整
			# 先复制不需要修改的部分
			# 起始温度上，Ts的索引偏移+1，状态转移上移1行;
			# 结束温度上，Ts的索引偏移+1，状态转移左移1列
			# 上移m行，左移n列，保留[i,end-m+1:end,end-n+1:end]
			C[i,:,:]
		end
		# 如果没有范围调整，那么减小间隔
		if flag_nextgap
			# =dT_origin/2*(1+0.5*(rand()-0.5))
			dT_origin = dT_origin / 2 * (1 + 0.5 * (rand() - 0.5))
			for j ∈ 1:nt
				TsMatrix[:, j] = TsList[j]-half_nT*dT_origin:dT_origin:(TsList[j]+half_nT*dT_origin+1e-8)
			end
			countSingleGap = 0
		else
			countSingleGap += 1
			if countSingleGap > 6
				dT_origin = min(dT_origin * 2 * (1 + 0.5 * (rand() - 0.5)), dt * k_dT_to_dt)
				for j ∈ 1:nt
					TsMatrix[:, j] = TsList[j]-half_nT*dT_origin:dT_origin:(TsList[j]+half_nT*dT_origin+1e-8)
				end
				countSingleGap = 0
			end
		end
		countAll += 1
		#println("countAll:$countAll", " dT_origin:$(round(dT_origin,digits=4))", " cost:$(round(cost,digits=4))")
	end


	P1List = map(i -> P1Matrix[i, TsIndex[i], TsIndex[i+1]], 1:nt-1)
	P2List = map(i -> P2Matrix[i, TsIndex[i], TsIndex[i+1]], 1:nt-1)
	P3List = map(i -> P3Matrix[i, TsIndex[i], TsIndex[i+1]], 1:nt-1)
	PeList = map(i -> PeMatrix[i, TsIndex[i], TsIndex[i+1]], 1:nt-1)
	realCostList = map(i -> C[i, TsIndex[i], TsIndex[i+1]], 1:nt-1)
	realCostList .-= smoother * (sum(abs2.(P1List)) + sum(abs2.(P2List)) + sum(abs2.(P3List)) + sum(abs2.(PeList)))

	return sum(realCostList), TsList, P1List, P2List, P3List, PeList, realCostList
end

"""
返回给定初始温度下的最优解
"""
function dpSolve(::VaryLoadVaryArea;
	C::Array{Float64, 3},
	j::Int,# The index of start temperature
)
	nT = size(C, 2)
	nt = size(C, 1)
	TsTransitionMatrix = zeros(nT, nt)
	# 第一步
	VForward = C[1, j, :]
	TsTransitionMatrix[:, 1] .= j#TsTransitionMatrix[:,i]表示第i+1个时层上的前驱节点

	for i in 2:nt
		VForward, TsTransitionMatrix[:, i] = forwardSolve(
			VaryLoadVaryArea(), VForward, C[i, :, :], nT,
		)
	end

	valueMin = VForward[j]
	# 状态回溯
	#valueMin=VForward[j] # 在第j个温度下的最优成本
	TsIndexList = Vector{Int}(undef, nt + 1)
	TsIndexList[nt+1] = j
	#println("begin nt=",nt," j=",j)
	for i in nt:-1:2
		#println(i," ",j," ",TsIndexList[i+1])
		TsIndexList[i] = TsTransitionMatrix[TsIndexList[i+1], i]
	end
	TsIndexList[1] = j

	return valueMin, TsIndexList
end

function forwardSolve(::VaryLoadVaryArea, VForward::Vector, C::Matrix, nT::Int)
	VForwardNew = fill(99999.0, nT)
	lastTsIndex = zeros(nT)
	for j in 1:nT
		for i in 1:nT
			if VForward[i] + C[i, j] < VForwardNew[j]
				VForwardNew[j] = VForward[i] + C[i, j]
				lastTsIndex[j] = i
			end
		end
	end
	return VForwardNew, lastTsIndex
end

function GoldenRatioSolver(
	nT::Int,
	nt::Int,
	C::Array{Float64, 3},
)
	# 初始化黄金分割法
	phi = 0.618
	jList = [1, nT - round(Int, phi * (nT - 1)), 1 + round(Int, phi * (nT - 1)), nT]
	valueList = zeros(4)
	TsIndexListMatirx = zeros(Int, nt, 4)
	for (i, j) in enumerate(jList)
		temp1, temp2 = dpSolve(
			VaryLoadVaryArea();
			C = C,
			j = j,
		)

		valueList[i], TsIndexListMatirx[:, i] = temp1, temp2

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
				VaryLoadVaryArea();
				C = C,
				j = jList[2],
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
				VaryLoadVaryArea();
				C = C,
				j = jList[3],
			)
		end
		count += 1
	end

	minCost, index = findmin(valueList)
	minTsList = TsIndexListMatirx[:, index]

	return minCost, minTsList
end

function ExhaustiveSolver(
	nT::Int,
	nt::Int,
	C::Array{Float64, 3},
)
	# 初始化黄金分割法
	valueList = zeros(nT)
	jList = 1:nT
	TsIndexListMatirx = zeros(Int, nt, nT)
	for j in jList
		temp1, temp2 = dpSolve(
			VaryLoadVaryArea();
			C = C,
			j = j,
		)
		valueList[j], TsIndexListMatirx[:, j] = temp1, temp2
	end

	minCost, index = findmin(valueList)
	minTsList = TsIndexListMatirx[:, index]

	return minCost, minTsList
end

"""
将矩阵上移m行，左移n列
"""
function moveCost(C_i::Array,m::Int,n::Int)
	
end

