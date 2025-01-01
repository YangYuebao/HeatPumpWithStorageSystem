

"""
基于物性的COP计算
"""
function COPTe_Tc(Te, Tc, eta_s, refrigerant)
	h1 = CoolProp.PropsSI("H", "T", Te + 273.15, "Q", 1, refrigerant)
	s1 = CoolProp.PropsSI("S", "T", Te + 273.15, "Q", 1, refrigerant)
	p2 = CoolProp.PropsSI("P", "T", Tc + 273.15, "Q", 1, refrigerant)
	h2 = CoolProp.PropsSI("H", "S", s1, "P", p2, refrigerant)
	wt = (h2 - h1) / eta_s  # 理论压缩功等于绝热压缩功除以绝热效率
	h3 = CoolProp.PropsSI("H", "T", Tc + 273.15, "Q", 0, refrigerant)
	return (h1 - h3) / wt + 1 # 制热循环效率
end


"""
生成COP的函数、COP的梯度和海森矩阵
"""
function getCOP_g_h(
	minTe::Real,# 蒸发温度下限
	maxTe::Real,# 蒸发温度上限
	minTc::Real,# 冷凝温度下限
	maxTc::Real,# 冷凝温度上限
	refrigerant::String,# 工质
	maxCOP::Real,# 最大COP
	eta_s::Real,# 绝热效率
	dT::Real,# 插值步长
)
	TcChangeToElec = maxTc# 冷凝温度转换到电热的阈值
	# 生成基于物性的COP计算函数
	# 然后进行样条插值,并计算梯度和海森矩阵
	minTe = floor(minTe, digits = 1)
	maxTe = ceil(maxTe, digits = 1)
	minTc = floor(minTc, digits = 1)
	maxTc = ceil(maxTc, digits = 1)
	dT = floor(minTe, digits = 1)
	TeList = minTe:dT:maxTe
	TcList = minTc:dT:maxTc

	iStart = Int((minTe - 10) * 10 + 1)
	iEnd = Int((maxTe - 10) * 10 + 1)
	jStart = Int((minTc - 20) * 10 + 1)
	jEnd = Int((maxTc - 20) * 10 + 1)
	step = Int(dT * 10)

	#Threads.@threads for (j,Tc) in enumerate(TcList)
	if refrigerant in ["R134a"]
		dfTemp = CSV.read(joinpath(pwd(), "src", "refrigerantPropertys", "R134a_10_80_20_100_0.1.csv"), DataFrame)
	elseif refrigerant in ["water", "Water"]
		dfTemp = CSV.read(joinpath(pwd(), "src", "refrigerantPropertys", "water_70_190_70_190_0.1.csv"), DataFrame)
	end
	COPMatrix = Matrix(dfTemp[iStart:step:iEnd, jStart:step:jEnd])

	itpCOP = interpolate(COPMatrix, BSpline(Cubic(Line(OnGrid()))))
	sitpCOP = scale(itpCOP, TeList, TcList)

	function COPfunction(Te, Tc)
		COP = sitpCOP(Te, Tc)
		if Tc > TcChangeToElec
			return 1.0
		end
		if COP > maxCOP || COP <= 0
			return maxCOP
		end
		return COP
	end

	function COPfunction_g(Te, Tc)
		if (Tc > TcChangeToElec) || (COP > maxCOP) || (COP <= 0)
			return zeros(2)
		end
		return Interpolations.gradient(sitpCOPSupplyWaste, Te, Tc)
	end

	function COPfunction_h(Te, Tc)
		if (Tc > TcChangeToElec) || (COP > maxCOP) || (COP <= 0)
			return zeros(2, 2)
		end
		return Interpolations.hessian(sitpCOPSupplyWaste, Te, Tc)
	end

	return COPfunction, COPfunction_g, COPfunction_h
end

"""
计算低温闪蒸需要的焓
"""
function getLatenHeat()

end

"""
设计问题：计算管道的额定流量：
1kW的用热需求,额定的换热温差(5℃),用热温度Tuse为回水温度。假设蒸汽为饱和蒸汽,回水为饱和水,此时流量是多少？可能比实际偏大一些
"""
function qmperkW(Tuse::Float64, Trecycle::Float64)
	hSupply = CoolProp.PropsSI("H", "T", Tuse + 273.15, "Q", 1, "water")
	hRecycle = CoolProp.PropsSI("H", "T", Trecycle + 273.15, "Q", 0, "water")
	return 1 / (hSupply - hRecycle) * 1e3
end


"""
输入一定系统结构和工作参数,返回系统计算需要用到的各种向量
"""
function generateSystemCoff(::PressedWaterDoubleStorage;
	refrigerantLow::String = "R134a",    # 供热循环使用制冷剂
	refrigerantHigh::String = "water",   # 蓄热使用制冷剂
	maxTcHigh::Real = 180.0,  # 高温热泵冷凝器温度上限
	maxTcLow::Real = 80.0,  # 低温热泵冷凝器温度上限
	eta_s::Real = 0.7,                        # 压缩机绝热效率
	Twastein::Real = 80.0,                    # 废气进入温度
	Twasteout::Real = 40.0,                   # 废气排出温度
	dTwaste::Real = 5.0,                      # 废气温度-蒸发器温度
	Tair::Vector = fill(25.0, 24),             # 外部环境温度
	dTair::Real = 5.0,                        # 外部环境温度-蒸发器温度
	Tuse::Real = 120.0,                       # 工厂使用温度
	Trecycle::Real = 115.0,  # 回收蒸汽温度
	heatStorageCapacity::Real = 6.0,          # 蓄热量kWh(承压水蓄热)
	TstorageTankMax::Real = 220.0,            # 蓄热罐的最高温度
	maxheatStorageInputHour::Real = 4,        # 蓄热充满时长
	dTstorageInput::Real = 5.0,               # 蓄热温差
	kt::Real = 0.5,                           # 蓄热 \Delta T_2 / \Delta T_1
	dT_l::Real = 5.0,  # 低温热泵蒸发器与冷凝器传热温差
	TwastCapacity::Real = 0.8,                # 废热容量是工厂用热量的倍数
	heatStorageOutEfficiency::Real = 0.0001,  # 蓄热衰减系数K
	dT_h::Real = 5.0,  # 高温热泵蒸发器与冷凝器传热温差
	dTuse::Real = 5.0,                        # 冷凝器温度-工厂使用温度
	dTuseStandard::Real = 5.0,                # 工厂用热标准温差
	dTstorageStandard::Real = 5.0,            # 工厂蓄热标准温差
	dTstorageOutputStandard::Real = 5.0,      # 工厂蓄热释放标准温差
	heatConsumptionPower::Real = 1.0,         # 每小时用热功率kW
	maxCOP::Real = 21.0,                      # 热泵COP上限
	workingStartHour::Int = 8,                # 生产开始时间
	workingHours::Int = 16,                   # 每日工作小时数
	PheatPumpMax::Real = 1.0,                 # 热泵最大功率kW
	hourlyTariff::Vector = fill(0.7, 24),     # 电价向量
	COPInterpolateGap = 0.1,    # COP插值时步长
)
	temp = workingHours
	workingHours = workingHours % 24
	if workingHours == 0
		workingHours = 24
	end
	if temp > 24
		@warn "工作时长大于24小时,改为$(workingStartHour)h,终止计算"
		return -1
	end

	# 计算一天的热泵蒸发冷凝温度
	TeRecycle = (Twastein + Twasteout) / 2 - dTwaste    # 余热回收蒸发器温度
	TeAirSource = Tair .- dTair                      # 空气源蒸发器温度
	minTe = minimum(TeAirSource)
	maxTe = max(TeRecycle, TeAirSource...)
	minTc = minimum(Tair)

	COPh, COPh_g, COPh_h = getCOP_g_h(minTe, maxTe, minTc, maxTcHigh, refrigerantHigh, maxCOP, eta_s, COPInterpolateGap)
	COPl, COPl_g, COPl_h = getCOP_g_h(minTe, maxTe, Tair, maxTcLow, refrigerantLow, maxCOP, eta_s, COPInterpolateGap)

	# 生成总循环参数:T13::Real,T9::Real,qm::Vector
	heatConsumptionPowerList = zeros(24)
	heatConsumptionPowerList[1:workingHours] .= heatConsumptionPower
	heatConsumptionPowerList = rotateVectorForward(heatConsumptionPowerList, workingStartHour)
	T13 = Tuse
	T9 = Trecycle
	qmStandard = qmperkW(Tuse, Trecycle)
	qm = qmStandard * heatConsumptionPowerList
	latenHeat = 2150.0# 汽化潜热kJ/kg
	cp_cw = 4.275# 循环水定压热容kJ/kg
	cp_cs = 2.281# 循环蒸汽定压热容kJ/kg


	# 生成低温热泵参数：cpqm_l,k1,dTe_l1,dTe_l2,dTc_l,Qhrecycle
	cpm_l = cpm_h = heatStorageCapacity * heatConsumptionPower / (TstorageTankMax - Tuse)
	cpqm_h = heatStorageCapacity * heatConsumptionPower / (maxheatStorageInputHour * dTstorageInput)
	cpqm_m = heatStorageCapacity * heatConsumptionPower / (COPh(Te_hStandard, Tc_hStandard)) / (maxheatStorageInputHour * dTstorageInput)
	cpqm_l = cpqm_m
	k1 = k2 = k3 = k4 = k5 = k6 = k7 = kt
	dTc_l = dTe_l1 = dTe_l2 = dT_l
	Qhrecycle = fill(TwastCapacity * heatConsumptionPower, 24)

	# 生成低温蓄热参数：cpm_l,Tair,cp_cw,KTloss_l
	KTloss_l = KTloss_h = heatStorageOutEfficiency

	# 生成高温热泵参数:cpqm_m,k2,dTe_h,k3,dTc_h1,dTc_h2,cpqm_h,cp_cs
	dTe_h = dTc_h1 = dTc_h2 = dT_h

	# 高温蓄热参数:cpm_h,KTloss_h,k4,k5,k6,k7

	# 设备运行约束:TstorageTankMax,heatStorageVelocity,heatStorageCapacityConstraint,heatpumpPowerConstraint
	heatStorageVelocity = heatStorageCapacity / maxheatStorageInputHour
	heatStorageCapacityConstraint = heatStorageCapacity
	heatpumpPowerConstraint = PheatPumpMax


	# 其它参数:hourlyTariffList,heatConsumptionPowerList
	hourlyTariffList = hourlyTariff


	return COPh, COPh_g, COPh_h,
	COPl, COPl_g, COPl_h,
	T13, T9, qm, latenHeat, cp_cw, cp_cs,
	cpqm_l, k1, dTe_l1, dTe_l2, dTc_l, Qhrecycle,
	cpm_l, Tair, KTloss_l,
	cpqm_m, k2, dTe_h, k3, dTc_h1, dTc_h2, cpqm_h,
	cpm_h, KTloss_h, k4, k5, k6, k7
	TstorageTankMax, heatStorageVelocity, heatStorageCapacityConstraint, heatpumpPowerConstraint,
	hourlyTariffList, heatConsumptionPowerList
end

"""给定系统参数,求解系统成本,返回成本、热泵功率向量、蓄热量向量"""
function generateAndSolve(::PressedWaterDoubleStorage, ::MinimizeCost;
	COPh::Function,
	COPh_g::Function,
	COPh_h::Function,
	COPl::Function,
	COPl_g::Function,
	COPl_h::Function,

	# 总循环参数
	T13::Real,# 供热蒸汽温度
	T9::Real,# 蒸汽冷却循环水回水温度
	qm::Vector,# 总循环水质量流量
	latenHeat::Real,# 汽化潜热
	cp_cw::Real,# 循环水液态cp*qm cycled water
	cp_cs::Real,# 蒸汽定压热容 cp cycled steam

	# 低温热泵参数
	cpqm_l::Real,# 低温热泵换热工质循环cp*qm
	k1::Real,# 低温热泵热汇换热温差比ΔT_2/ΔT_1
	dTe_l1::Real,# 低温热泵蒸发器与废热源的换热温差
	dTe_l2::Real,# 低温热泵蒸发器与空气源的换热温差
	dTc_l::Real,# 低温热泵冷凝器与低温回路的换热温差
	Qhrecycle::Vector,# 废热源最大回收功率

	# 低温蓄热参数
	cpm_l::Real,# 低温蓄热罐热容
	Tair::Real,# 环境温度
	KTloss_l::Real,# 低温蓄热罐热损失系数,是个很小的数

	# 高温热泵参数
	cpqm_m::Real,# 高温热泵蒸发器侧循环工质cp*qm
	k2::Real,# 高温热泵热源换热温差比
	dTe_h::Real,# 高温热泵蒸发器与蓄热罐热源的换热温差
	k3::Real,# 高温热泵热汇换热温差比
	dTc_h1::Real,# 高温热泵冷凝器与高温蓄热罐侧循环工质的换热温差
	dTc_h2::Real,# 高温热泵冷凝器与用热支流qm4+qm1的换热温差
	cpqm_h::Real,# 高温热泵换热工质循环cp*qm

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容
	KTloss_h::Real,# 高温蓄热罐热损失系数,是个很小的数
	k4::Real,# 直接从高温蓄热取热的支路的温差比
	k5::Real,# 从高温罐取热后还需要加热的支路的温差比
	k6::Real,# 直接从高温蓄热取热的支路的温差比
	k7::Real,# 从高温罐取热后还需要加热的支路的温差比

	# 设备运行约束
	TstorageTankMax::Real,# 蓄热罐的最高温度
	heatStorageVelocity::Real,           # 蓄热速率约束
	heatStorageCapacityConstraint::Float64, # 蓄热量约束(最大值)
	heatpumpPowerConstraint::Float64,   # 热泵功率约束(最大值)

	# 其它参数
	hourlyTariffList::Vector{Float64},   # 电价向量
	heatConsumptionPowerList::Vector{Float64},  # 用热负载向量
)
	m = 24
	model = Model(Ipopt.Optimizer)
	set_silent(model)

	# 定义COP函数
	@operator(model, opCOPl, 2, COPl, COPl_g, COPl_h)
	@operator(model, opCOPh, 2, COPh, COPh_g, COPh_h)

	#=
		@variable(model, 0 <= heatPumpPowerl1[1:m] <= heatpumpPowerConstraint)  # 供热回收
		@constraint(model, [i = 1:m], cpqmh * 2 * kt / (1 + kt) * (T1[i] - T[i]) <= heatpumpPowerConstraint)
	=#

	@variable(model, qm1[i = 1:m] >= 0)
	@variable(model, qm2[i = 1:m] >= 0)
	@variable(model, qm3[i = 1:m] >= 0)
	@variable(model, qm4[i = 1:m] >= 0)
	@variable(model, qm5[i = 1:m] >= 0)
	@variable(model, qm6[i = 1:m] >= 0)
	@variable(model, Q[i = 1:m] >= 0)

	@variable(model, heatPumpPowerl1[i = 1:m] >= 0)
	@variable(model, T10[i = 1:m] >= Tair)
	@variable(model, T11[i = 1:m] >= Tair)
	@variable(model, T12[i = 1:m] >= Tair)
	@variable(model, T14[i = 1:m] >= Tair)

	#@constraint(model,[i=1:m],)
	@constraint(model, cons1[i = 1:m],qm[i]==qm1[i]+qm2[i]+qm3[i]+qm4[i]+qm5[i]+qm6[i])
	


	# 目标函数
	@objective(model, Min, sum(hourlyTariffList[i] * (heatPumpPowerl1[i] + heatPumpPowerl2[i] + heatPumpPowerh1[i] + heatPumpPowerh2[i]) for i ∈ 1:m))

	optimize!(model)

	isFesable = primal_status(model) # ==FEASIBLE_POINT ? true : false

	Pl1List = value.(heatPumpPowerl1)
	Pl2List = value.(heatPumpPowerl2)
	Ph1List = value.(heatPumpPowerh1)
	Ph2List = value.(heatPumpPowerh2)

	PList = Pl1List .+ Pl2List .+ Ph1List .+ Ph2List
	TstorageList = value.(T)
	T1List = value.(T1)
	T3List = value.(T3)
	T4List = fill(T4, m)
	T5List = value.(T5)
	COPl1 = COPSupplyWaste.(T3List)
	COPl2 = COPSupplyAir.(T3List)
	COPh1 = COPStorageWaste.(T1List)
	COPh2 = COPStorageAir.(T1List)
	heatConsumptionPowerList = heatConsumptionPowerList
	costList = PList .* hourlyTariffList
	heatStorageList = cpm * (value.(T) .- T4)
	#=
	蓄热的增加和减少有待修改
	=#

	return isFesable, Pl1List, Pl2List, Ph1List, Ph2List, PList, TstorageList, T1List, T3List, T4List, T5List, COPl1, COPl2, COPh1, COPh2, heatConsumptionPowerList, costList, heatStorageList
end
