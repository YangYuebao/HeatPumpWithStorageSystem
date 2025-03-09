

"""
输入一定系统结构和工作参数,返回系统计算需要用到的各种向量
"""
function generateSystemCoff(::PressedWaterDoubleStorageOneCompressor;
	maxTcHigh::Real = 180.0,  				# 高温热泵冷凝器温度上限
	maxTcLow::Real = 90.0,  				# 低温热泵冷凝器温度上限
	TWaste::Real = 60.0,                  	# 废热源温度
	Tair::Vector = fill(25.0, 24),          # 外部环境温度
	dTlc_he::Real = 10.0,  					# 高温热泵蒸发器温度与低温热泵冷凝器温度差
	Tuse::Real = 120.0,                     # 工厂使用温度
	Trecycle::Real = 115.0,    				# 回收蒸汽温度
	heatStorageCapacity::Real = 6.0,        # 蓄热量kWh(承压水蓄热)
	TstorageTankMax::Real = 220.0,          # 蓄热罐的最高温度
	maxheatStorageInputHour::Real = 4,      # 蓄热充满时长
	dT_EvaporationStandard::Real = 3.0,		# 全蒸温差
	dT_l::Real = 5.0,    					# 低温热泵蒸发器与冷凝器传热温差,也用作低温系统的传热温差
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

	return (T10, T9, latentHeat, cp_cw, cp_cs,
		dTe_l1, dTe_l2, dTc_l, TwastCapacity, Tair, TWaste,
		cpm_l,  dT_EvaporationStandard, minTeh,
		cpm_h, 
		TstorageTankMax, heatStorageVelocity, heatpumpPowerConstraint,
		hourlyTariffList, heatConsumptionPowerList,
		maxTcLow,
		storageTankVolume)
end



"""给定系统参数,求解系统成本,返回成本、热泵功率向量、蓄热量向量"""
function generateAndSolve(::PressedWaterDoubleStorageOneCompressor, ::MinimizeCost;
	COP1::Function,
	COP2::Function,
	COP3::Function,

	# 总循环参数
	T10::Real,# 供热蒸汽温度
	T9::Real,# 蒸汽冷却循环水回水温度
	latentHeat::Real,# 汽化潜热
	cp_cw::Real,# 循环水定压热容
	cp_cs::Real,# 蒸汽定压热容 cp cycled steam
	

	# 低温热泵参数
	TwastCapacity::Vector,# 废热源最大回收功率
	Tair::Function,# 环境温度
	TWaste::Real,# 废热回收蒸发器温度

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容

	# 设备运行约束
	TstorageTankMax::Real,				# 蓄热罐的最高温度
	heatpumpPowerConstraint::Float64,   # 热泵功率约束(最大值)

	# 其它参数
	hourlyTariffFunction::Function,   # 电价向量
	heatConsumptionPowerFunction::Function,  # 用热负载向量

)
	#hourlyTariffList::Vector{Float64},   # 电价向量
	#heatConsumptionPowerList::Vector{Float64},  # 用热负载向量
	#TairList
end
