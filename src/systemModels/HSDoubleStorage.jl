

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
	dT = floor(dT, digits = 1)
	TeList = minTe:dT:maxTe
	TcList = minTc:dT:maxTc

	step = round(Int,dT * 10)

	#Threads.@threads for (j,Tc) in enumerate(TcList)
	if refrigerant in ["R134a"]
		dfTemp = CSV.read(joinpath(pwd(), "src", "refrigerantPropertys", "R134a_10_80_20_100_0.1.csv"), DataFrame)
		iStart = round(Int,(minTe - 10) * 10 + 1)
		iEnd = round(Int,(maxTe - 10) * 10 + 1)
		jStart = round(Int,(minTc - 20) * 10 + 1)
		jEnd = round(Int,(maxTc - 20) * 10 + 1)
	elseif refrigerant in ["water", "Water"]
		dfTemp = CSV.read(joinpath(pwd(), "src", "refrigerantPropertys", "water_70_190_70_190_0.1.csv"), DataFrame)
		iStart = round(Int,(minTe - 70) * 10 + 1)
		iEnd = round(Int,(maxTe - 70) * 10 + 1)
		jStart = round(Int,(minTc - 70) * 10 + 1)
		jEnd = round(Int,(maxTc - 70) * 10 + 1)
		#@info any(isnan.(dfTemp))
	end
	println(minTe," ", maxTe," ", minTc," ", maxTc)
	println(size(dfTemp))
	println(iStart," ", iEnd," ", jStart," ", jEnd," ", step)
	COPMatrix = Matrix(dfTemp[iStart:step:iEnd, jStart:step:jEnd])*eta_s .+ 1.0

	m,n=size(COPMatrix)
	count=0
	for j=1:n
		for i=1:m
			if COPMatrix[i,j]>=21.0
				COPMatrix[i,j]=21.0
			end
			if COPMatrix[i,j]<=-21.0
				COPMatrix[i,j]=-21.0
			end
			count+=1
		end
	end

	itpCOP = interpolate(COPMatrix, BSpline(Cubic(Line(OnGrid()))))
	sitpCOP = scale(itpCOP, TeList, TcList)

	function COPfunction(x::T ...)::T where {T<:Real}
		COP = sitpCOP(x[1], x[2])
		if x[2] > TcChangeToElec
			return 1.0
		end
		if COP > maxCOP || COP <= 0
			return maxCOP
		end
		return COP
	end

	function COPfunction_g(g::AbstractVector{T}, x::T ...)::Nothing where {T<:Real}
		COP = sitpCOP(x[1], x[2])
		if (x[2] > TcChangeToElec) || (COP > maxCOP) || (COP <= 0)
			g[1],g[2] = zeros(2)
			return nothing
		end
		g[1],g[2] = Interpolations.gradient(sitpCOP, x[1], x[2]) |> Vector
		return nothing
	end

	function COPfunction_h(H::AbstractMatrix{T}, x::T ...)::Nothing where {T<:Real}
		COP = sitpCOP(x[1], x[2])
		if (x[2] > TcChangeToElec) || (COP > maxCOP) || (COP <= 0)
			H[1,1],H[2,1],_,H[2,2]=zeros(4)
			return nothing
		end
		H[1,1],H[2,1],_,H[2,2]=Interpolations.hessian(sitpCOP, x[1], x[2]) |> Matrix

		H=LowerTriangular(H)

		return nothing
	end

	return COPfunction, COPfunction_g, COPfunction_h
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
	refrigerantLow::String = "R134a",    		# 供热循环使用制冷剂
	refrigerantHigh::String = "water",   		# 蓄热使用制冷剂
	maxTcHigh::Real = 180.0,  					# 高温热泵冷凝器温度上限
	maxTcLow::Real = 90.0,  					# 低温热泵冷凝器温度上限
	eta_s::Real = 0.7,                        	# 压缩机绝热效率
	Twastein::Real = 80.0,                    	# 废气进入温度
	Twasteout::Real = 40.0,                   	# 废气排出温度
	Tair::Vector = fill(25.0, 24),            	# 外部环境温度
	dTair::Real = 5.0,                        	# 外部环境温度-蒸发器温度
	dTlc_he::Real = 10.0,					  	# 高温热泵蒸发器温度与低温热泵冷凝器温度差
	maxTeh::Real = 180.0,					  	# 高温热泵蒸发器温度上限
	Tuse::Real = 120.0,                       	# 工厂使用温度
	Trecycle::Real = 115.0,  				  	# 回收蒸汽温度
	heatStorageCapacity::Real = 6.0,          	# 蓄热量kWh(承压水蓄热)
	TstorageTankMax::Real = 220.0,            	# 蓄热罐的最高温度
	maxheatStorageInputHour::Real = 4,        	# 蓄热充满时长
	dTstorageInput::Real = 5.0,               	# 蓄热温差
	kt::Real = 0.5,                           	# 蓄热 \Delta T_2 / \Delta T_1
	dT_l::Real = 5.0,  						  	# 低温热泵蒸发器与冷凝器传热温差,也用作低温系统的传热温差
	TwastCapacity::Real = 0.8,                	# 废热容量是工厂用热量的倍数
	heatStorageOutEfficiency::Real = 0.0001,  	# 蓄热衰减系数K
	dT_h::Real = 5.0,  						  	# 高温热泵蒸发器与冷凝器传热温差
	heatConsumptionPower::Real = 1.0,         	# 每小时用热功率kW
	maxCOP::Real = 21.0,                      	# 热泵COP上限
	workingStartHour::Int = 8,                	# 生产开始时间
	workingHours::Int = 16,                   	# 每日工作小时数
	PheatPumpMax::Real = 1.0,                 	# 热泵最大功率kW
	hourlyTariff::Vector = fill(0.7, 24),     	# 电价向量
	COPInterpolateGap = 0.1,    				# COP插值时步长
	Te_hStandard::Real =80.0					# 高温热泵蒸发器温度标准值
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
	TWaste = (Twastein + Twasteout) / 2 		# 余热回收温度
	TeAirSource = Tair .- dTair                	# 空气源蒸发器温度
	minTel = minimum(TeAirSource)
	maxTel = max(TWaste-dT_l, TeAirSource...)
	minTcl = minimum(Tair)

	minTeh = maxTcLow-dTlc_he
	minTch = minTeh + 0.1

	
	#println("minTeh:", minTeh, " maxTeh:", maxTeh)
	#println("minTch:", minTch, " maxTch:", maxTcHigh)
	COPh, COPh_g, COPh_h = getCOP_g_h(minTeh, maxTeh, minTch, maxTcHigh, refrigerantHigh, maxCOP, eta_s, COPInterpolateGap)
	#println("minTel:", minTel, " maxTel:", maxTel)
	#println("minTcl:", minTcl, " maxTcl:", maxTcLow)
	COPl, COPl_g, COPl_h = getCOP_g_h(minTel, maxTel, minTcl, maxTcLow, refrigerantLow, maxCOP, eta_s, COPInterpolateGap)

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


	# 生成低温热泵参数：cpqm_l,k1,dTe_l1,dTe_l2,dTc_l,QhRecycle
	cpm_l = cpm_h = heatStorageCapacity * heatConsumptionPower / (TstorageTankMax - Tuse)
	cpqm_h = heatStorageCapacity * heatConsumptionPower / (maxheatStorageInputHour * dTstorageInput)
	cpqm_m = heatStorageCapacity * heatConsumptionPower / (COPh(Te_hStandard, Tuse+dT_h)) / (maxheatStorageInputHour * dTstorageInput)
	cpqm_l = cpqm_m
	k1 = k2 = k3 = k4 = k5 = k6 = k7 = k8 = kt
	dTc_l = dTe_l1 = dTe_l2 = dT_l
	QhRecycle = fill(TwastCapacity * heatConsumptionPower, 24)

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


	return (COPh, COPh_g, COPh_h,
	COPl, COPl_g, COPl_h,
	T13, T9, qm, latenHeat, cp_cw, cp_cs,
	cpqm_l, k1, dTe_l1, dTe_l2, dTc_l, QhRecycle, Tair,TWaste,
	cpm_l, KTloss_l, k6,
	cpqm_m, k2, dTe_h, k3, dTc_h1, dTc_h2, cpqm_h,minTeh,
	cpm_h, KTloss_h, k4, k5, k7,k8,
	TstorageTankMax, heatStorageVelocity, heatStorageCapacityConstraint, heatpumpPowerConstraint,
	hourlyTariffList, heatConsumptionPowerList)
end

"""给定系统参数,求解系统成本,返回成本、热泵功率向量、蓄热量向量"""
function generateAndSolve(::PressedWaterDoubleStorage, 	::MinimizeCost;
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
	QhRecycle::Vector,# 废热源最大回收功率
	Tair::Vector,# 环境温度
	TWaste::Real,# 废热回收蒸发器温度

	# 低温蓄热参数
	cpm_l::Real,# 低温蓄热罐热容
	KTloss_l::Real,# 低温蓄热罐热损失系数,是个很小的数
	k6::Real,# 从低温蓄热取热的支路的温差比

	# 高温热泵参数
	cpqm_m::Real,# 高温热泵蒸发器侧循环工质cp*qm
	k2::Real,# 高温热泵热源换热温差比
	dTe_h::Real,# 高温热泵蒸发器与蓄热罐热源的换热温差
	k3::Real,# 高温热泵热汇换热温差比
	dTc_h1::Real,# 高温热泵冷凝器与高温蓄热罐侧循环工质的换热温差
	dTc_h2::Real,# 高温热泵冷凝器与用热支流qm4+qm1的换热温差
	cpqm_h::Real,# 高温热泵换热工质循环cp*qm
	minTeh::Real,# 高温热泵蒸发器温度下限

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容
	KTloss_h::Real,# 高温蓄热罐热损失系数,是个很小的数
	k4::Real,# 直接从高温蓄热取热的支路的温差比
	k5::Real,# 从高温罐取热后还需要加热的支路的温差比
	k7::Real,# 直接从高温蓄热取热的支路的温差比
	k8::Real,# 从高温罐取热后还需要加热的支路的温差比

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
		@variable(model, 0 <= P_l1[1:m] <= heatpumpPowerConstraint)  # 供热回收
		@constraint(model, [i = 1:m], cpqmh * 2 * kt / (1 + kt) * (T1[i] - T[i]) <= heatpumpPowerConstraint)
	=#

	Tair_min= minimum(Tair)					# 环境最低温度
	@variable(model, qm1[i = 1:m] >= 0)	# 质量流量:低温蓄热->高温热泵->供热
	@variable(model, qm2[i = 1:m] >= 0)	# 质量流量:低温蓄热->高温蓄热->供热
	@variable(model, qm3[i = 1:m] >= 0)	# 质量流量:低温蓄热->供热
	@variable(model, qm4[i = 1:m] >= 0)	# 质量流量:低温蓄热->高温蓄热->高温热泵->供热
	@variable(model, qm5[i = 1:m] >= 0)	# 质量流量:高温蓄热->供热
	@variable(model, qm6[i = 1:m] >= 0)	# 质量流量:高温蓄热->高温热泵->供热
	@variable(model, Qc_l[i = 1:m] >= 0)	# 低温热泵输出热功率
	@variable(model, T1[i = 1:m] >= Tair_min)	# 低温热泵出水温度=低温蓄热加热进水温度
	@variable(model, T3[i = 1:m] >= Tair_min)	# 低温蓄热温度
	@variable(model, T5[i = 1:m] >= minTeh)	# 高温热泵蒸发器从低温蓄热取热的温度
	@variable(model, T6[i = 1:m] >= Tair_min)	# 高温热泵冷凝器向高温蓄热储热的温度
	@variable(model, T8[i = 1:m] >= minTeh+0.1)	# 高温蓄热温度
	# T9是常数，工厂回水温度
	@variable(model, T10[i = 1:m] >= T9)	# 低温蓄热供热温度
	@variable(model, T11[i = 1:m] >= T9)	# 高温蓄热供热温度
	@variable(model, T12[i = 1:m] >= minTeh)	# 再热蒸汽排出高温热泵冷凝器温度
	# T13是常数，工厂用蒸汽温度
	@variable(model, T14[i = 1:m] >= minTeh)	# 低温再热蒸汽进入高温热泵冷凝器温度
	@variable(model, T15[i = 1:m] >= minTeh)	# 总再热蒸汽进入高温热泵冷凝器温度
	@variable(model, T16[i = 1:m] >= T9)		# 回水再热蒸汽排出高温蓄热罐温度
	@variable(model, T17[i = 1:m] >= minTeh)	# 回水供热蒸汽排出高温蓄热罐温度
	@variable(model, T18[i = 1:m] >= minTeh)	# 低温供热蒸汽排出高温蓄热罐温度
	@variable(model, Tc_l[i = 1:m] >= Tair_min)	# 低温热泵冷凝器温度
	# Te_l1 = TWaste .- dTe_l1 	# 低温热泵废热回收蒸发器温度，确定值直接计算于下方
	# Te_l2 = Tair .- dTe_l2	# 低温热泵空气源蒸发器温度
	@variable(model, P_l1[i = 1:m] >= 0)	# 低温热泵废热源功率
	@variable(model, P_l2[i = 1:m] >= 0)	# 低温热泵空气源功率
	@variable(model, P_h1[i = 1:m] >= 0)	# 高温蓄热热泵功率
	@variable(model, P_h2[i = 1:m] >= 0)	# 高温供热热泵功率
	@variable(model, Qs_l[i = 1:m] >= 0)	# 低温蓄热量
	@variable(model, Qs_h[i = 1:m] >= 0)	# 高温蓄热量
	@variable(model, 0 <= lambda1[i = 1:m] <=1)	# 低温蓄热罐直接换热的流量比
	@variable(model, 0 <= lambda2[i = 1:m] <=1)	# 高温蓄热罐再热支路的流量比
	@variable(model, 0 <= lambda3[i = 1:m] <=1)	# 高温蓄热罐直供支路的流量比
	@variable(model, 0 <= lambda4[i = 1:m] <=1)	# 高温蓄热罐回水直供支路的流量比
	@variable(model, 0 <= lambda5[i = 1:m] <=1)	# 高温蓄热罐回水再热支路的流量比
	@variable(model, Te_h[i = 1:m] >= minTeh)	# 高温热泵蒸发器温度
	@variable(model, Tc_h1[i = 1:m] >= minTeh+0.1)	# 高温热泵蒸发器温度——蓄热
	@variable(model, Tc_h2[i = 1:m] >= minTeh+0.1)	# 高温热泵蒸发器温度——供热
	@variable(model, Qe_h[i = 1:m] >= 0)	# 高温热泵蒸发器取热功率
	@variable(model, Qc_h1[i = 1:m] >= 0)	# 高温热泵冷凝器蓄热功率
	@variable(model, Qc_h2[i = 1:m] >= 0)	# 高温热泵冷凝器供热功率

	# 总循环回水约束
	# 1. 循环水质量守恒
	@constraint(model, cons01[i = 1:m],qm[i]==qm1[i]+qm2[i]+qm3[i]+qm4[i]+qm5[i]+qm6[i])
	
	# 低温热泵约束
	# 2. 低温热泵输出热量功率=低温蓄热加热功率
	@constraint(model, cons02[i = 1:m],Qc_l[i]==cpqm_l * 2 * k1 / (1 + k1) * (T1[i] - T3[i]))
	# 3. 低温热泵冷凝器温度
	@constraint(model, cons03[i = 1:m],Tc_l[i]==1 / (1 + k1) * T1[i] + k1 / (1+k1) * T3[i] + dTc_l)
	# 4.5. 低温热泵蒸发器温度
	Te_l1 = TWaste - dTe_l1	# 低温热泵废热回收蒸发器温度
	Te_l2 = Tair .- dTe_l2		# 低温热泵空气源蒸发器温度
	# 6. 低温热泵输出热量功率与COP的关系
	@constraint(model, cons06[i = 1:m],Qc_l[i]==opCOPl(Te_l1,Tc_l[i]) * P_l1[i]+opCOPl(Te_l2[i],Tc_l[i]) * P_l2[i])
	# 7. 废热源热量约束
	@constraint(model, cons07[i = 1:m],(opCOPl(Te_l1,Tc_l[i])-1) * P_l1[i] <= QhRecycle[i])
	# 8. 温度关系约束
	@constraint(model, cons08[i = 1:m],T1[i]>=T3[i])

	# 低温蓄热约束
	# 9. 低温蓄热量约束
	@constraint(model, cons09[i = 1:m],Qs_l[i]==cpm_l*(T3[i]-Tair_min))
	# 10.蓄热罐能量守恒
	@constraint(model, cons10_1[i = 1:m-1],cpm_l*(T3[i+1]-T3[i])==(cp_cs*(T10[i]-T9)+latenHeat)*(qm1[i]+qm2[i]+qm3[i]+qm4[i])+Qc_l[i]-Qe_h[i]-KTloss_l*cpm_l*(T3[i+1]-Tair[i+1]))
	@constraint(model, cons10_2,cpm_l*(T3[1]-T3[m])==(cp_cs*(T10[m]-T9)+latenHeat)*(qm1[m]+qm2[m]+qm3[m]+qm4[m])+Qc_l[m]-Qe_h[m]-KTloss_l*cpm_l*(T3[1]-Tair[1]))

	# 11. 温度关系约束——蓄热罐温度>=环境温度
	@constraint(model, cons11[i = 1:m],T3[i]>=Tair[i])
	# 12. 温度关系约束——传热温度
	@constraint(model, cons12[i = 1:m],T10[i]==lambda1[i]*(T9*(1-k6)/(1+k6)+T3[i]*2*k6/(1+k6))+(1-lambda1[i])*T9)

	# 高温热泵约束
	# 13. 高温热泵取热功率=低温蓄热释放功率
	@constraint(model, cons13[i = 1:m],Qe_h[i]==cpqm_m*(T3[i]-T5[i])*2*k2/(1+k2))
	# 14.15. 高温热泵冷凝器温度
	@constraint(model, cons14[i = 1:m],Tc_h1[i]==T6[i]/(1+k3)+T8[i]*k3/(1+k3) + dTc_h1)
	@constraint(model, cons15[i = 1:m],Tc_h2[i]==(T15[i] + T12[i])/2 + dTc_h2)
	# 16. 高温热泵蒸发器温度
	@constraint(model, cons16[i = 1:m],Te_h[i]==T5[i]/(1+k2)+T3[i]*k2/(1+k2) - dTe_h)
	# 17. 混流温度
	@constraint(model, cons17[i = 1:m],(qm1[i]+qm4[i]+qm6[i])*T15[i]==qm1[i]*T10[i]+qm4[i]*T14[i]+qm6[i]*T16[i])
	# 18. 高温热泵输出蓄热功率
	@constraint(model, cons18[i = 1:m],Qc_h1[i]==cpqm_h*(T6[i]-T8[i])*2*k3/(1+k3))
	# 19. 高温热泵输出供热功率
	@constraint(model, cons19[i = 1:m],Qc_h2[i]==cp_cs*(qm1[i]+qm4[i]+qm6[i])*(T12[i]-T15[i]))
	# 20.21. 高温热泵功率与COP
	@constraint(model, cons20[i = 1:m],Qc_h1[i]==opCOPh(Te_h[i],Tc_h1[i])*P_h1[i])
	@constraint(model, cons21[i = 1:m],Qc_h2[i]==opCOPh(Te_h[i],Tc_h2[i])*P_h2[i])
	# 22.~27. 温度约束
	@constraint(model, cons22[i = 1:m],T6[i]>=T8[i])
	@constraint(model, cons23[i = 1:m],T14[i]==lambda2[i]*(T10[i]*(1-k5)/(1+k5)+T8[i]*2*k5/(1+k5))+(1-lambda2[i])*T10[i])
	#@constraint(model, cons24[i = 1:m],T14[i]>=T15[i])
	@constraint(model, cons25[i = 1:m],T16[i]>=T9)
	@constraint(model, cons26[i = 1:m],T15[i]>=T10[i])
	@constraint(model, cons27[i = 1:m],T12[i]>=T15[i])

	# 高温蓄热约束
	# 28.高温蓄热量
	@constraint(model, cons28[i = 1:m],Qs_l[i]==cpm_l*(T8[i]-Tair_min))
	# 29.蓄热罐能量守恒
	@constraint(model, cons29_1[i = 1:m-1],cpm_h*(T8[i+1]-T8[i])==cp_cs*(qm4[i]*(T10[i]-T14[i])+qm2[i]*(T10[i]-T18[i]))-qm5[i]*(cp_cs*(T17[i]-T9)+latenHeat)-qm6[i]*(cp_cs*(T16[i]-T9)+latenHeat)-KTloss_h*cpm_h*(T8[i+1]-Tair[i+1])+Qc_h1[i])
	@constraint(model, cons29_2,cpm_h*(T8[1]-T8[m])==cp_cs*(qm4[m]*(T10[m]-T14[m])+qm2[m]*(T10[m]-T18[m]))-qm5[m]*(cp_cs*(T17[m]-T9)+latenHeat)-qm6[m]*(cp_cs*(T16[m]-T9)+latenHeat)-KTloss_h*cpm_h*(T8[1]-Tair[1])+Qc_h1[m])
	# 30.~36. 温度约束
	@constraint(model, cons30[i = 1:m],T18[i]>=T10[i])
	@constraint(model, cons31[i = 1:m],T14[i]>=T10[i])
	@constraint(model, cons32[i = 1:m],T17[i]>=T9)
	@constraint(model, cons33[i = 1:m],T16[i]>=T9)
	# T10[i]==lambda1[i]*(T9*(1-k6)/(1+k6)+T3[i]*2*k6/(1+k6))+(1-lambda1[i])*T9
	@constraint(model, cons34[i = 1:m],T16[i]==lambda5[i]*(T9*(1-k8)/(1+k8)+T8[i]*2*k8/(1+k8))+(1-lambda5[i])*T9)
	@constraint(model, cons35[i = 1:m],T17[i]==lambda4[i]*(T9*(1-k7)/(1+k7)+T8[i]*2*k7/(1+k7))+(1-lambda4[i])*T9)
	@constraint(model, cons36[i = 1:m],T18[i]==lambda3[i]*(T10[i]*(1-k4)/(1+k4)+T8[i]*2*k4/(1+k4))+(1-lambda3[i])*T10[i])
	# T14在23
	
	# 总循环供蒸汽约束
	# 37. 温度约束
	@constraint(model, cons37[i = 1:m],qm[i]*T13==T10[i]*qm3[i]+T18[i]*qm2[i]+T17[i]*qm5[i]+(qm1[i]+qm4[i]+qm6[i])*T12[i])
	# 目标函数
	@objective(model, Min, sum(hourlyTariffList[i] * (P_l1[i] + P_l2[i] + P_h1[i] + P_h2[i]) for i ∈ 1:m))

	optimize!(model)

	isFesable = primal_status(model) # ==FEASIBLE_POINT ? true : false

	Pl1List = value.(P_l1)
	Pl2List = value.(P_l2)
	Ph1List = value.(P_h1)
	Ph2List = value.(P_h2)

	#=
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
	=#
	return Pl1List,Pl2List
end
