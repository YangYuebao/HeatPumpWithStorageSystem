

"""
输入一定系统结构和工作参数,返回系统计算需要用到的各种向量
"""
function generateSystemCoff(::PressedWaterDoubleStorageSimplified;
	maxTcHigh::Real = 180.0,  				# 高温热泵冷凝器温度上限
	maxTcLow::Real = 90.0,  				# 低温热泵冷凝器温度上限
	TWaste::Real = 60.0,                  	# 废热源温度
	Tair::Vector = fill(25.0, 24),          # 外部环境温度
	dTlc_he::Real = 10.0,  # 高温热泵蒸发器温度与低温热泵冷凝器温度差
	Tuse::Real = 120.0,                     # 工厂使用温度
	Trecycle::Real = 115.0,    # 回收蒸汽温度
	heatStorageCapacity::Real = 6.0,        # 蓄热量kWh(承压水蓄热)
	TstorageTankMax::Real = 220.0,          # 蓄热罐的最高温度
	maxheatStorageInputHour::Real = 4,      # 蓄热充满时长
	dT_EvaporationStandard::Real = 3.0,# 全蒸温差
	dT_l::Real = 5.0,    # 低温热泵蒸发器与冷凝器传热温差,也用作低温系统的传热温差
	TwastCapacity::Real = 0.8,              # 废热容量是工厂用热量的倍数
	heatStorageOutEfficiency::Real = 0.0001,# 蓄热衰减系数K
	heatConsumptionPower::Real = 1.0,       # 每小时用热功率kW
	workingStartHour::Int = 8,              # 生产开始时间
	workingHours::Int = 16,                 # 每日工作小时数
	PheatPumpMax::Real = 1.0,               # 热泵最大功率kW
	hourlyTariff::Vector = fill(0.7, 24),   # 电价向量
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

	minTeh = maxTcLow - dTlc_he

	# 生成总循环参数:T10::Real,T9::Real,qm::Vector
	heatConsumptionPowerList = zeros(24)
	heatConsumptionPowerList[1:workingHours] .= heatConsumptionPower
	heatConsumptionPowerList = rotateVectorForward(heatConsumptionPowerList, workingStartHour)
	T10 = Tuse
	T9 = Trecycle
	qmStandard = qmperkW(Tuse, Trecycle)
	qm = qmStandard * heatConsumptionPowerList#kg/s
	qm .+= 1e-8
	latentHeat = 2150.0# 汽化潜热kJ/kg
	cp_cw = 4.275# 循环水定压热容kJ/kg
	cp_cs = 2.281# 循环蒸汽定压热容kJ/kg


	# 生成低温热泵参数：cpqm_l,k1,dTe_l1,dTe_l2,dTc_l,QhRecycle
	cpm_l = cpm_h = heatStorageCapacity * heatConsumptionPower / (TstorageTankMax - Tuse)# kWh/K
	storageTankMass=cpm_h/cp_cw*3600#kg
	storageTankVolume=storageTankMass/900#m^3
	#cpqm_h = heatStorageCapacity * heatConsumptionPower / (maxheatStorageInputHour * dTstorageInput)# kWh/(K*h)
	#cpqm_m = cpqm_h * 3.0#*((COPh(Te_hStandard, Tuse + dT_h))-1.0)/ (COPh(Te_hStandard, Tuse + dT_h)) 
	#cpqm_l = cpqm_m# kWh/(K*h)
	#k1 = kt
	dTc_l = dTe_l1 = dTe_l2 = dT_l
	QhRecycle = TwastCapacity * heatConsumptionPowerList# kWh

	# 生成低温蓄热参数：cpm_l,Tair,cp_cw,KTloss_l
	KTloss_l = KTloss_h = heatStorageOutEfficiency
	#Qe_hStandard=1.0
	dT_EvaporationStandard = dT_EvaporationStandard

	# 高温蓄热参数:cpm_h,KTloss_h

	# 设备运行约束:TstorageTankMax,heatStorageVelocity,heatStorageCapacityConstraint,heatpumpPowerConstraint
	heatStorageVelocity = heatStorageCapacity / maxheatStorageInputHour * heatConsumptionPower# kWh/h
	heatStorageCapacityConstraint = heatStorageCapacity * heatConsumptionPower# kWh
	heatpumpPowerConstraint = PheatPumpMax# kW


	# 其它参数:hourlyTariffList,heatConsumptionPowerList
	hourlyTariffList = hourlyTariff
	TcChangeToElec = maxTcHigh

	return (T10, T9, qm, latentHeat, cp_cs,
		dTe_l1, dTe_l2, dTc_l, QhRecycle, Tair, TWaste,
		cpm_l, KTloss_l, dT_EvaporationStandard, minTeh,
		cpm_h, KTloss_h,
		TstorageTankMax, heatStorageVelocity, heatpumpPowerConstraint,
		hourlyTariffList, heatConsumptionPowerList,
		maxTcLow,
		storageTankVolume)
end



"""给定系统参数,求解系统成本,返回成本、热泵功率向量、蓄热量向量"""
function generateAndSolve(::PressedWaterDoubleStorageSimplified, ::MinimizeCost;
	COPh::Function,
	COPh_g::Function,
	COPh_h::Function,

	# 总循环参数
	T10::Real,# 供热蒸汽温度
	T9::Real,# 蒸汽冷却循环水回水温度
	qm::Vector,# 总循环水质量流量
	latentHeat::Real,# 汽化潜热
	cp_cs::Real,# 蒸汽定压热容 cp cycled steam

	# 低温热泵参数
	dTe_l1::Real,# 低温热泵蒸发器与废热源的换热温差
	dTe_l2::Real,# 低温热泵蒸发器与空气源的换热温差
	dTc_l::Real,# 低温热泵冷凝器与低温回路的换热温差
	QhRecycle::Vector,# 废热源最大回收功率
	Tair::Vector,# 环境温度
	TWaste::Real,# 废热回收蒸发器温度

	# 低温蓄热参数
	cpm_l::Real,# 低温蓄热罐热容
	KTloss_l::Real,# 低温蓄热罐热损失系数,是个很小的数
	dT_EvaporationStandard::Real,# 全蒸换热温差
	minTeh::Real,# 高温热泵蒸发器温度下限

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容
	KTloss_h::Real,# 高温蓄热罐热损失系数,是个很小的数

	# 设备运行约束
	TstorageTankMax::Real,# 蓄热罐的最高温度
	heatStorageVelocity::Real,           # 蓄热速率约束
	heatpumpPowerConstraint::Float64,   # 热泵功率约束(最大值)

	# 其它参数
	hourlyTariffList::Vector{Float64},   # 电价向量
	heatConsumptionPowerList::Vector{Float64},  # 用热负载向量

	# 初值
	maxTcLow::Real,# 低温热泵冷凝器温度上限

	reinitialize::Bool,	# 将初值重置为T3取最小值，不使用蓄热的情况
	qm1_init::Vector{Float64},
	qm2_init::Vector{Float64},
	qm3_init::Vector{Float64},

	T3_init::Vector{Float64},
	T8_init::Vector{Float64},
	Tc_l_init::Vector{Float64},

	P_l1_init::Vector{Float64},
	P_l2_init::Vector{Float64},

	P_h1_init::Vector{Float64},
	P_h2_init::Vector{Float64},
	P_h3_init::Vector{Float64},
	
	Qc_l_init::Vector{Float64},
	
	Qs_l_init::Vector{Float64},
	Qs_h_init::Vector{Float64},
	aimObjectValue::Real
)
	m = 24
	model = Model(Ipopt.Optimizer)
	heatpumpPowerConstraint *= 10
	set_silent(model)
	set_attribute(model, "warm_start_init_point", "yes")
	set_attribute(model, "alpha_for_y", "min")
	set_attribute(model, "max_iter", 64000)
	set_attribute(model, "acceptable_iter", 160)
	set_attribute(model, "acceptable_tol", 5e-3)
	set_attribute(model, "tol", 5e-4)

	# 定义COP函数
	@operator(model, opCOPl, 2, COPl, COPl_g, COPl_h)
	@operator(model, opCOPh, 2, COPh, COPh_g, COPh_h)
	
	# 定义变量与初值计算
	# 常数变量
	Tair_min = minimum(Tair)# 环境最低温度
	Te_l1 = TWaste - dTe_l1	# 低温热泵废热源蒸发器温度
	Te_l2 = Tair .- dTe_l2	# 低温热泵空气源蒸发器温度

	# 计算初值
	maxT3=maxTcLow-dTc_l
	minT3=minTeh+dT_EvaporationStandard
	#T3_init=rand(m)*(maxT3-minT3).+minT3
	
	if reinitialize==true
		qm1_init = copy(qm)
		qm2_init = zeros(m)
		qm3_init = zeros(m)

		T3_init=fill(minT3,m)
		T8_init=fill(T10+dT_EvaporationStandard,m)
		Tc_l_init=T3_init.+dTc_l

		COPl2_init=map(i->COPl(Te_l2[i],Tc_l_init[i]),1:m)
		COPh1_init=COPh(minT3-dT_EvaporationStandard,T10)

		P_l1_init = zeros(m)
		P_l2_init = @. qm*(cp_cs*(T3_init-dT_EvaporationStandard-T9)+latentHeat)./(COPl2_init)

		P_h1_init = @. qm*(cp_cs*(T3_init-dT_EvaporationStandard-T9)+latentHeat)/(COPh1_init-1)
		P_h2_init = zeros(m)
		P_h3_init = zeros(m)
		
		Qc_l_init = P_l2_init.*COPl2_init
		
		Qs_l_init = cpm_l*(T3_init.-Tair_min)
		Qs_h_init = cpm_l*(T8_init.-Tair_min)
	end

	#=
	dfInitValue=DataFrame(
		maxT3=fill(maxT3,m),
		minT3=fill(minT3,m),
		maxTcLow=fill(maxTcLow,m),
		minTcLow=fill(minTeh,m),
		T3=T3_init,
		Tc_l=Tc_l_init,
		COPh1=COPh1_init,
		P_h1=P_h1_init,
		COPl2=COPl2_init,
		P_l2=P_l2_init
	)
	CSV.write(joinpath(pwd(), "calculations", "situation6", "initValue.csv"), round.(dfInitValue,digits=4))
	=#

	# 定义变量
	@variable(model, qm1[i = 1:m] >= 0, start = qm1_init[i])# 质量流量:低温蓄热->高温热泵->供热
	@variable(model, qm2[i = 1:m] >= 0, start = qm2_init[i])# 质量流量:高温蓄热->高温热泵->供热
	@variable(model, qm3[i = 1:m] >= 0, start = qm3_init[i])# 质量流量:高温蓄热

	@variable(model, Qc_l[i = 1:m] >= 0,start = Qc_l_init[i])# 低温热泵输出热功率
	
	#@variable(model, TstorageTankMax >= T3[i = 1:m] >= Tair_min, start = maxTcLow)# 低温蓄热温度
	@variable(model, maxT3 >= T3[i = 1:m] >= minT3, start = T3_init[i])# 低温蓄热温度

	@variable(model, TstorageTankMax >= T8[i = 1:m] >= maxTcLow, start = T8_init[i])# 高温蓄热温度
	# T9是常数，工厂回水温度
	# T10是常数，工厂用蒸汽温度
	@variable(model, maxTcLow >= Tc_l[i = 1:m] >= minT3, start = Tc_l_init[i]) # 低温热泵冷凝器温度
	# Te_l1 = TWaste .- dTe_l1 					# 低温热泵废热回收蒸发器温度，确定值直接计算于下方
	# Te_l2 = Tair .- dTe_l2					# 低温热泵空气源蒸发器温度
	@variable(model, heatpumpPowerConstraint >= P_l1[i = 1:m] >= 0,start=P_l1_init[i])# 低温热泵废热源功率
	@variable(model, heatpumpPowerConstraint >= P_l2[i = 1:m] >= 0,start=P_l2_init[i])# 低温热泵空气源功率
	@variable(model, heatpumpPowerConstraint >= P_h1[i = 1:m] >= 0,start=P_h1_init[i])# 高温供热热泵功率
	@variable(model, heatpumpPowerConstraint >= P_h2[i = 1:m] >= 0,start=P_h2_init[i])# 高温取热热泵功率
	@variable(model, heatpumpPowerConstraint >= P_h3[i = 1:m] >= 0,start=P_h3_init[i])# 高温储热热泵功率

	@variable(model, Qs_l[i = 1:m] >= 0,start=Qs_l_init[i])# 低温蓄热量
	@variable(model, Qs_h[i = 1:m] >= 0,start=Qs_h_init[i])# 高温蓄热量
	

	# 总循环回水约束
	# 1. 循环水质量守恒
	@constraint(model, cons01[i = 1:m], qm[i] == qm1[i] + qm2[i])

	# 低温热泵约束
	# 3. 低温热泵冷凝器温度
	@constraint(model, cons03[i = 1:m], Tc_l[i] ==  T3[i] + dTc_l)
	# 4.5. 低温热泵蒸发器温度 低温热泵废热回收蒸发器温度 低温热泵空气源蒸发器温度
	 
	# 6. 低温热泵输出热量功率与COP的关系
	@constraint(model, cons06[i = 1:m], Qc_l[i] == opCOPl(Te_l1, Tc_l[i]) * P_l1[i] + opCOPl(Te_l2[i], Tc_l[i]) * P_l2[i])
	# 7. 废热源热量约束
	@constraint(model, cons07[i = 1:m], (opCOPl(Te_l1, Tc_l[i]) - 1) * P_l1[i] <= QhRecycle[i])

	# 低温蓄热约束
	# 9. 低温蓄热量约束
	@constraint(model, cons09[i = 1:m], Qs_l[i] == cpm_l * (T3[i] - Tair_min))
	# 10.蓄热罐能量守恒
	@constraint(model, cons10_1[i = 1:m-1], cpm_l * (T3[i+1] - T3[i]) == Qc_l[i] - (cp_cs * (T3[i] - dT_EvaporationStandard - T9) + latentHeat) * (qm1[i] + qm3[i]) - KTloss_l * cpm_l * (T3[i+1] - Tair[i+1]))
	@constraint(model, cons10_2, cpm_l * (T3[1] - T3[m]) == Qc_l[m] - (cp_cs * (T3[m] - dT_EvaporationStandard - T9) + latentHeat) * (qm1[m] + qm3[m]) - KTloss_l * cpm_l * (T3[1] - Tair[1]))

	# 11. 温度关系约束——蓄热罐温度>=环境温度
	@constraint(model, cons11[i = 1:m], T3[i] >= Tair[i])

	# 高温部分
	# 13. 供热从低温罐取热量
	@constraint(model, cons13[i = 1:m], (cp_cs * (T3[i] - dT_EvaporationStandard - T9) + latentHeat) * qm1[i] == P_h1[i] * (opCOPh(T3[i] - dT_EvaporationStandard, T10) - 1))

	# 16. 高温蓄热取热约束
	@constraint(model, cons16[i = 1:m], qm2[i] * (cp_cs * (T10 - T9) + latentHeat) == P_h2[i] * (opCOPh(T8[i] - dT_EvaporationStandard, T10) - 1))
	# 17. 高温蓄热储热从低温罐取热
	@constraint(model, cons17[i = 1:m], qm3[i] * (cp_cs * (T3[i] - dT_EvaporationStandard - T9) + latentHeat) == P_h3[i] * (opCOPh(T3[i] - dT_EvaporationStandard, T8[i] + dT_EvaporationStandard) - 1))
	
	# 18.19. 蓄热速率约束
	#@constraint(model, cons18[i = 1:m], Qc_l[i]<=heatStorageVelocity)
	@constraint(model, cons19[i = 1:m], P_h3[i] * opCOPh(T3[i] - dT_EvaporationStandard, T8[i] + dT_EvaporationStandard) <= heatStorageVelocity)

	# 20. 低温罐温度约束
	#@constraint(model, cons20[i = 1:m], T3[i]>=minTeh)

	# 21. 罐体温度约束
	@constraint(model, cons21[i = 1:m], T8[i]>=T3[i])

	# 22. 费用约束
	@constraint(model,cons22,sum(hourlyTariffList[i] * (P_l1[i] + P_l2[i] + P_h1[i] + P_h2[i] + P_h3[i]) for i ∈ 1:m)<=aimObjectValue)

	# 28.高温蓄热量
	@constraint(model, cons28[i = 1:m], Qs_h[i] == cpm_l * (T8[i] - Tair_min))
	# 29.蓄热罐能量守恒
	@constraint(
		model,
		cons29_1[i = 1:m-1],
		cpm_h * (T8[i+1] - T8[i]) == P_h3[i] * opCOPh(T3[i] - dT_EvaporationStandard, T8[i] + dT_EvaporationStandard) - qm2[i] * (cp_cs * (T10 - T9) + latentHeat) - KTloss_h * cpm_h * (T8[i+1] - Tair[i+1])
	)
	@constraint(
		model,
		cons29_2,
		cpm_h * (T8[1] - T8[m]) == P_h3[m] * opCOPh(T3[m] - dT_EvaporationStandard, T8[m] + dT_EvaporationStandard) - qm2[m] * (cp_cs * (T10 - T9) + latentHeat) - KTloss_h * cpm_h * (T8[1] - Tair[1])
	)

	# errReportStart=primal_feasibility_report(model) do v
	# 	return start_value(v)
	# end

	# 目标函数
	@objective(model, Min, sum(hourlyTariffList[i] * (P_l1[i] + P_l2[i] + P_h1[i] + P_h2[i] + P_h3[i]) for i ∈ 1:m))

	optimize!(model)

	isFeasible = primal_status(model) # ==FEASIBLE_POINT ? true : false
	#println("isFeasible: $isFeasible")
	#errReportSolved=primal_feasibility_report(model)


	T4 = value.(T3) .- dT_EvaporationStandard
	T5 = value.(T8) .- dT_EvaporationStandard
	T6 = value.(T8) .+ dT_EvaporationStandard

	Pl1List = value.(P_l1)
	Pl2List = value.(P_l2)
	Ph1List = value.(P_h1)
	Ph2List = value.(P_h2)
	Ph3List = value.(P_h3)
	COPl1 = map(i -> COPl(Te_l1, value(Tc_l[i])), 1:m)
	COPl2 = map(i -> COPl(Te_l2[i], value(Tc_l[i])), 1:m)
	COPh1 = map(i -> COPh(T4[i], value(T10)), 1:m)
	COPh2 = map(i -> COPh(T5[i], value(T10)), 1:m)
	COPh3 = map(i -> COPh(T4[i], T6[i]), 1:m)

	dfResult = DataFrame()
	if isFeasible in [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]
		dfResult = DataFrame(
			"时间" => 0:23,
			"电价" => hourlyTariffList,
			"用热需求" => heatConsumptionPowerList,
			"质量流量_总kg/h" => qm * 3600,
			"直供流量qm1" => value.(qm1) * 3600,
			"高蓄取热流量qm2" => value.(qm2) * 3600,
			"高蓄蓄热流量qm3" => value.(qm3) * 3600,
			"低温罐温度T3" => value.(T3),
			"高温罐温度T8" => value.(T8),
			"低温热泵冷凝器温度Tc_l" => value.(Tc_l),
			"冷凝水温度T9" => fill(T9, 24),
			"用热温度T10" => fill(T10, 24),
			"低温热泵蒸发器温度_热回收TWast" => fill(TWaste, 24),
			"低温热泵蒸发器温度_空气源Tair" => Tair,
			"低温热回收热泵COPl1" => COPl1,
			"低温热回收热泵功率P_l1" => Pl1List,
			"热回收功率上限QhRecycle" => QhRecycle,
			"低温空气源热泵COPl2" => COPl2,
			"低温空气源热泵功率P_l2" => Pl2List,
			"高温供热热泵COPh1" => COPh1,
			"高温供热热泵功率P_h1" => Ph1List,
			"高温取热热泵COPh2" => COPh2,
			"高温取热热泵功率P_h2" => Ph2List,
			"高温蓄热热泵COPh3" => COPh3,
			"高温蓄热热泵功率P_h3" => Ph3List,
			"低温热泵输出热功率Qc_l" => value.(Qc_l),
			"低温罐蓄热量Qs_l" => value.(Qs_l),
			"高温罐蓄热量Qs_h" => value.(Qs_h),
			"热回收热泵制热量" => COPl1 .* Pl1List,
			"空气源热泵制热量" => COPl2 .* Pl2List,
			"低温循环水取热量" => value.(qm1 + qm3) .* (cp_cs * (T4 .- T9) .+ latentHeat),
		)
		#CSV.write(joinpath(pwd(), "calculations", "situation6", "result.csv"), round.(dfResult, digits = 4))
		#@info "The result is saved in $(joinpath(pwd(),"calculations","situation6","result.csv"))"
	end

	#println("cost per day:",objective_value(model))
	return isFeasible, objective_value(model), dfResult#,errReportStart,errReportSolved
end
