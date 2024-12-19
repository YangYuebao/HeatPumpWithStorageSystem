

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
	minTe::Real,	# 蒸发温度下限
	maxTe::Real,	# 蒸发温度上限
	minTc::Real,	# 冷凝温度下限
	maxTc::Real,	# 冷凝温度上限
	refrigerant::String,	# 工质
	maxCOP::Real,			# 最大COP
	TcChangeToElec::Real,	# 冷凝温度转换到电热的阈值
	eta_s::Real,			# 绝热效率
	dT::Real				# 插值步长
)
	# 生成基于物性的COP计算函数
	# 然后进行样条插值，并计算梯度和海森矩阵
	TeList=minTe:dT:maxTe
	TcList=minTc:dT:maxTc
	
	COPMatrix = zeros(length(TeList), length(TcList))

	Threads.@threads for (j,Tc) in enumerate(TcList)
		for (i,Te) in enumerate(TeList)
			COPMatrix[i,j] = COPTe_Tc(Te, Tc, eta_s, refrigerant)
		end
	end
	
	itpCOP = interpolate(COPMatrix, BSpline(Cubic(Line(OnGrid()))))
	sitpCOP = scale(itpCOP, TeList,TcList)

	function COPfunction(Te,Tc)
		COP=sitpCOP(Te,Tc)
		if Tc>TcChangeToElec
			return 1.0
		end
		if COP > maxCOP || COP <=0
			return maxCOP
		end
		return COP
	end

	function COPfunction_g(Te,Tc)
		if (Tc>TcChangeToElec) || (COP > maxCOP) || (COP <=0)
			return zeros(2)
		end
		return Interpolations.gradient(sitpCOPSupplyWaste, Te,Tc)
	end

	function COPfunction_h(Te,Tc)
		if (Tc>TcChangeToElec) || (COP > maxCOP) || (COP <=0)
			return zeros(2,2)
		end
		return Interpolations.hessian(sitpCOPSupplyWaste, Te,Tc)
	end

	return COPfunction, COPfunction_g, COPfunction_h
end


"""
输入一定系统结构和工作参数,返回系统计算需要用到的各种向量
"""
function generateSystemCoff(::PressedWaterDoubleStorage;
	refrigerantLow::String = "R134a",  		  # 供热循环使用制冷剂
	refrigerantHigh::String = "water", 		  # 蓄热使用制冷剂
	maxTcHigh::Real = 180.0,				  # 高温热泵冷凝器温度上限
	maxTcLow::Real = 80.0,					  # 低温热泵冷凝器温度上限
	eta_s::Real = 0.7,                        # 压缩机绝热效率
	Twastein::Real = 80.0,                    # 废气进入温度
	Twasteout::Real = 40.0,                   # 废气排出温度
	dTwaste::Real = 5.0,                      # 废气温度-蒸发器温度
	TwastCapacity::Real = 0.8,                # 废热容量是工厂用热量的倍数
	Tair::Vector = fill(25.0,24),             # 外部环境温度
	dTair::Real = 5.0,                        # 外部环境温度-蒸发器温度
	Tuse::Real = 120.0,                       # 工厂使用温度
	dTuse::Real = 5.0,                        # 冷凝器温度-工厂使用温度
	dTstorageInput::Real = 5.0,               # 冷凝器温度-蓄热温度
	dTuseStandard::Real = 5.0,                # 工厂用热标准温差
	dTstorageStandard::Real = 5.0,            # 工厂蓄热标准温差
	dTstorageOutputStandard::Real = 5.0,      # 工厂蓄热释放标准温差
	TstorageTankMax::Real = 220.0,            # 蓄热罐的最高温度
	heatConsumptionPower::Real = 1.0,         # 每小时用热功率kW
	heatStorageCapacity::Real = 6.0,          # 蓄热量kWh(承压水蓄热)
	heatStorageOutEfficiency::Real = 0.0001,  # 蓄热衰减系数K
	maxheatStorageInputHour::Real = 4,        # 蓄热充满时长
    maxCOP::Real = 21.0,                      # 热泵COP上限
	workingStartHour::Int = 8,                # 生产开始时间
    workingHours::Int = 16,                   # 每日工作小时数
	PheatPumpMax::Real = 1.0,                 # 热泵最大功率kW
	hourlyTariff::Vector = fill(0.7, 24),     # 电价向量
	COPInterpolateGap = 0.1,  				  # COP插值时步长
    kt::Real = 0.5,                           # 蓄热 \Delta T_2 / \Delta T_1
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
	minTe=minimum(TeAirSource)
	maxTe=maximum(TeAirSource)

	COPh,COPh_g,COPh_h=getCOP_g_h(minTe,maxTe,Tair,maxTcHigh,refrigerantHigh,maxCOP,TcChangeToElec,eta_s,COPInterpolateGap)

	COPl,COPl_g,COPl_h=getCOP_g_h(minTe,maxTe,Tair,maxTcLow,refrigerantLow,maxCOP,TcChangeToElec,eta_s,COPInterpolateGap)

	# 生成用热负载向量
	heatConsumptionPowerList = zeros(24)
	heatConsumptionPowerList[1:workingHours] .= heatConsumptionPower
	heatConsumptionPowerList = rotateVectorForward(heatConsumptionPowerList, workingStartHour)

	# 生成蓄热量约束
	heatStorageCapacityConstraint = heatStorageCapacity

	# 生成热泵功率约束
	heatpumpPowerConstraint = PheatPumpMax

	# 生成电价向量
	hourlyTariffList = hourlyTariff

	# 生成充热速率
	heatStorageVelocity = heatStorageCapacity / maxheatStorageInputHour

	# 生成用热温度
	T4 = Tuse

	# 计算c_p*q_ml
	cpqml = heatConsumptionPower / dTuseStandard

	# 计算c_p*q_mh
	cpqmh = heatStorageVelocity / dTstorageStandard

	# 计算cpm
	cpm = heatStorageCapacity / (TstorageTankMax - Tuse + dTstorageOutputStandard)

	# 生成kt
	kt = kt

	# 生成废热容量倍数
	TwastCapacity = TwastCapacity

	# 生成蓄热罐最高温度
	TstorageTankMax = TstorageTankMax

	# 蓄热自损失
	heatStorageOutEfficiency=heatStorageOutEfficiency

	return COPSupplyWaste, COPSupplyAir, COPStorageWaste, COPStorageAir,
	COPSupplyWaste_g, COPSupplyAir_g, COPStorageWaste_g, COPStorageAir_g,
	COPSupplyWaste_h, COPSupplyAir_h, COPStorageWaste_h, COPStorageAir_h,
	heatConsumptionPowerList,
	heatStorageCapacityConstraint,
	heatpumpPowerConstraint,
	hourlyTariffList,
	heatStorageVelocity,
	T4,
	cpqml,
	cpqmh,
	cpm,
	kt,
	TwastCapacity,
	TstorageTankMax,
	heatStorageOutEfficiency
end

"""给定系统参数，求解系统成本，返回成本、热泵功率向量、蓄热量向量"""
function generateAndSolve(::PressedWaterDoubleStorage, ::MinimizeCost;
	COPh::Function,
	COPh_g::Function,
	COPh_h::Function,
	COPl::Function,
	COPl_g::Function,
	COPl_h::Function,
	
	# 总循环参数
	T13::Vector,	# 供热蒸汽温度
	T9::Vector,		# 蒸汽冷却循环水回水温度
	qm::Real,		# 总循环水质量流量
	
	# 低温热泵参数
	cpqm_l::Real,	# 低温热泵换热工质循环cp*qm
	k1::Real,		# 低温热泵热汇换热温差比ΔT_2/ΔT_1
	dTe_l1::Real,	# 低温热泵蒸发器与废热源的换热温差
	dTe_l2::Real,	# 低温热泵蒸发器与空气源的换热温差
	dTc_l::Real,	# 低温热泵冷凝器与低温回路的换热温差
	Qhrecycle::Vector,# 废热源最大回收功率

	# 低温蓄热参数
	cpm_l::Real,	# 低温蓄热罐热容
	Tair::Real,		# 环境温度
	cpqm_cw::Real,	# 循环水液态cp*qm cycled water
	KTloss_l::Real,	# 低温蓄热罐热损失系数，是个很小的数

	# 高温热泵参数
	cpqm_m::Real,	# 高温热泵蒸发器侧循环工质cp*qm
	k2::Real,		# 高温热泵热源换热温差比
	dTe_h::Real,	# 高温热泵蒸发器与蓄热罐热源的换热温差
	k3::Real,		# 高温热泵热汇换热温差比
	dTc_h1::Real,	# 高温热泵冷凝器与高温蓄热罐侧循环工质的换热温差
	dTc_h2::Real,	# 高温热泵冷凝器与用热支流qm4+qm1的换热温差
	cpqm_h::Real,	# 高温热泵换热工质循环cp*qm
	cp_cs::Real,	# 蒸汽定压热容 cp cycled steam

	# 高温蓄热参数
	cpm_h::Real,	# 高温蓄热热容
	KTloss_h::Real,	# 高温蓄热罐热损失系数，是个很小的数
	k4::Real,		# 直接从高温蓄热取热的支路的温差比
	k5::Real,		# 从高温罐取热后还需要加热的支路的温差比

	# 设备运行约束
	TstorageTankMax::Real,						# 蓄热罐的最高温度
	heatStorageVelocity::Real,           		# 蓄热速率约束
	heatStorageCapacityConstraint::Float64, 	# 蓄热量约束（最大值）
	heatpumpPowerConstraint::Float64,   		# 热泵功率约束（最大值）
	
	# 其它参数
	hourlyTariffList::Vector{Float64},   		# 电价向量
	heatConsumptionPowerList::Vector{Float64},  # 用热负载向量
)
	model = Model(Ipopt.Optimizer)
	set_silent(model)

	@operator(model, opCOPl, 2, COPl, COPl_g, COPl_h)
	@operator(model, opCOPh, 2, COPh, COPh_g, COPh_h)
	# 建立热泵功率变量与储热量变量，直接约束了变量的范围
	m = 24
	#=
		@variable(model, 0 <= heatPumpPowerl1[1:m] <= heatpumpPowerConstraint)  # 供热回收
		@constraint(model, [i = 1:m], cpqmh * 2 * kt / (1 + kt) * (T1[i] - T[i]) <= heatpumpPowerConstraint)
	=#
	@variable(model,T10[i=1:m]>=Tair)
	@variable(model,T11[i=1:m]>=Tair)
	@variable(model,T12[i=1:m]>=Tair)
	@variable(model,T14[i=1:m]>=Tair)

	@constraint(model,[i=1:m],qqq)

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

	return isFesable,Pl1List, Pl2List, Ph1List, Ph2List, PList,TstorageList,T1List,T3List,T4List,T5List,COPl1,COPl2,COPh1,COPh2,heatConsumptionPowerList, costList, heatStorageList
end
