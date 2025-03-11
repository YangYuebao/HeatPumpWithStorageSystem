
abstract type FunctionGenerateMethod end
struct BackwardGenerate <: FunctionGenerateMethod end
struct LinearGenerate <: FunctionGenerateMethod end
struct BackwardGenerateCycled <: FunctionGenerateMethod end
struct LinearGenerateCycled <: FunctionGenerateMethod end

function functionInterpolationGenerator(list::Vector,duration::Real,::BackwardGenerate)
	n=length(list)
	dt=duration/n
	function f(x::Real)
		if x >= duration
			x = x % duration
		elseif x < 0
			x = x % duration + duration
		end
		return list[floor(Int, x/dt)+1]
	end
	return f
end
functionInterpolationGenerator(list::Vector,duration::Real,::BackwardGenerateCycled)=Interpolation(list,duration,BackwardGenerate())

function functionInterpolationGenerator(list::Vector,duration::Real,::LinearGenerate)
	n=length(list)
	tList=collect(range(0,stop=duration,length=n))
	sitpCOP = linear_interpolation(tList,list)
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

function functionInterpolationGenerator(list::Vector,duration::Real,::LinearGenerateCycled)
	if abs(list[1]-list[end])<1e-6
		return Interpolation(list,duration,LinearGenerate())
	else
		@warn "列表首尾不一致，已经自动添加一位"
		return Interpolation(vcat(list,list[1]),duration,LinearGenerate())
	end
end


"""生成电价函数"""
function generateGridPriceFunction(gridPriceList::Vector,duration::Real)
	functionInterpolationGenerator(gridPriceList,duration,BackwardGenerate())
end

"""生成负载需求函数"""
function generateLoadFunction(loadList::Vector,duration::Real)
	functionInterpolationGenerator(loadList,duration,LinearGenerateCycled())
end

"""生成环境温度函数"""
function generateAreaTemperatureFunction(temperatureList::Vector,duration::Real)
	functionInterpolationGenerator(temperatureList,duration,LinearGenerateCycled())
end

"""
输入一定系统结构和工作参数,返回系统计算需要用到的各种向量
"""
function generateSystemCoff(::PressedWaterDoubleStorageOneCompressor;
	maxTcHigh::Real = 180.0,  				# 高温热泵冷凝器温度上限
	TCompressorIn::Real								# 中间温度
	TWaste::Real = 60.0,                  	# 废热源温度
	Tair::Vector = fill(25.0, 24),          # 外部环境温度
	dTlc_he::Real = 5.0,  					# 高温热泵蒸发器温度与低温热泵冷凝器温度差
	Tuse::Real = 120.0,                     # 工厂使用温度
	Trecycle::Real = 115.0,    				# 回收蒸汽温度
	heatStorageCapacity::Real = 6.0,        # 蓄热量kWh(承压水蓄热)
	TstorageTankMax::Real = 220.0,          # 蓄热罐的最高温度
	maxheatStorageInputHour::Real = 4,      # 蓄热充满时长
	dT_EvaporationStandard::Real = 5.0,		# 全蒸温差
	dT_l::Real = 5.0,    					# 低温热泵蒸发器与冷凝器传热温差,也用作低温系统的传热温差
	TwastCapacity::Real = 0.8,              # 废热容量是工厂用热量的倍数
	heatConsumptionPower::Real = 1.0,       # 每小时用热功率kW
	workingStartHour::Int = 8,              # 生产开始时间
	workingHours::Int = 16,                 # 每日工作小时数
	PheatPumpMax::Real = 1.2,               # 热泵功率服务系数
	elecHeat::Real=1.5,						# 电锅炉服务系数
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
	#QhRecycle = TwastCapacity * heatConsumptionPowerList# kWh

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

abstract type SimplifiedType end
struct NoSimplify <: SimplifiedType end
struct ConstloadandArea <: SimplifiedType end

"""给定系统参数,求解系统成本,返回成本、热泵功率向量、蓄热量向量"""
function generateAndSolve(::PressedWaterDoubleStorageOneCompressor, ::MinimizeCost;
	COPOverlap::Function,
	COPWater::Function,
	hourlyTariffFunction::Function,   # 电价函数
	heatConsumptionPowerFunction::Function,  # 用热负载函数
	TairFunction::Function,# 环境温度函数

	# 总循环参数
	Tuse::Real,# 供热蒸汽温度
	Tcondense::Real,# 蒸汽循环冷凝水回水温度
	TCompressorIn::Real,# 中间级温度
	dT_EvaporationStandard::Real,#全蒸温差
	latentHeat::Real,# 汽化潜热
	cp_cw::Real,# 循环水定压热容
	cp_cs::Real,# 蒸汽定压热容 cp cycled steam
	
	TwastCapacity::Real,# 废热源最大回收比例
	TWaste::Real,# 废热回收蒸发器温度

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容

	# 设备运行约束
	TstorageTankMax::Real,				# 蓄热罐的最高温度
	PheatPumpMax::Real,					# 热泵最大功率
	PelecHeatMax::Real,					# 电锅炉最大功率
	::ConstloadandArea;

	# 求解参数
	dT::Real=0.1,	# 状态参数高温蓄热温度离散步长
	dt::Real=1/6,	# 时间步长
)
	COP1=COP3=COPOverlap
	COP2=COPWater
	heatLoad=heatConsumptionPowerFunction(0.0)
	Tair=TairFunction(0.0)

	TsList = TCompressorIn+dT_EvaporationStandard:dT:TstorageTankMax
	tList=0:dt:24
	nT=length(TsList)	# 温度步数	
	nt=length(tList)		# 时间步数
	TsMatrix=zeros(nT,nt)	# 存储状态参数：高温蓄热温度
	#=
	"""计算当前温度下蓄热罐最高能升高多少度"""
	function TsIncreaseMax(TsNow,heatLoad)
		Timax=
	end

	"""计算当前温度下蓄热罐最低能降低多少度"""
	function TsDecreaseMax(TsNow,heatLoad)

	end
	=#

	#C[i,j]表示温度从TsList[i]到TsList[j]时的最低功率；在负载、环境不变的情况下，C[i,j]是不变的
	C=fill(99.0,nT,nT)# 状态转移矩阵
	TsDecreaseIndexList=zeros(nT)# 记录温度下降最多偏移的index
	TsIncreaseIndexList=zeros(nT)# 记录温度上升最多偏移的index
	for i=1:nT
		# 先计算温度下降的功率,存在一边电加热一边开2号模式的情况
		# 首先计算该模式在下一刻的最低温度
		TsNextMin=TsList[i]
		# 简单迭代求蓄热的最低温度
		for i=1:3
			TsNextMin=TsList[i]-heatLoad*dt/cpm_h*(1-1/COP2(0.5*(TsNextMin+TsList[i])-dT_EvaporationStandard,Tuse))
		end
		# 计算最低温度的下标
		TsNextMinGrid=max(ceil(TsNextMin,digits=1),TsList[1])
		TsDecreaseIndexList[i]=round(Int,(TsList[i]-TsNextMinGrid)*10)
		elecP=(TsNextMinGrid-TsNextMin)*cpm_h
		
		index = i-TsDecreaseIndexList[i]
		for j = index : i
			P2=heatLoad/(COP2((TsList[i]+TsList[j])/2-dT_EvaporationStandard,Tuse))
			elecP=(TsList[j]-TsNextMin)*cpm_h
			#补热的功率基本上超不了，因为这个工况下从蓄热取走的热量小于工厂的负载，而补热的功率大于等于工厂的负载
			if elecP<PelecHeatMax 
				C[i,j]=P2+elecP
			end
		end

		# 计算温度不变的功率
		P1Only=heatLoad/COP1(TWaste,Tuse)
		C[i,i]=P1Only

		# 计算温度上升的功率，首先计算该模式在下一刻的最高温度
		P3Max=PheatPumpMax-P1Only
		TsNextMax = TsList[i]
		# 当前温度小于用热温度
		if TsList[i]+dT_EvaporationStandard<Tuse
			TsNextMax=(cpm_h/dt*TsList[i]+PelecHeatMax+P3Max*COP3(TWaste,Tuse)*(1+cp_cw/latentHeat*(Tuse-(TsList[i])/2)))/(cpm_h/dt+0.5*cp_cw/latentHeat*P3Max*COP3(TWaste,Tuse))
			#TsNextMax=TsList[i]+dt/cpm_h*(PelecHeatMax+P3Max*COP3(TWaste,Tuse)*(1+cp_cw/latentHeat*(Tuse-(TsList[i]+TsNextMax)/2)))
		else # 当前温度大于等于用热温度
			P1=max(0,heatLoad/COP1(TWaste,(TsList[i]+TsNextMax)*0.5)-cp_cw/latentHeat*(0.5*(TsList[i]+TsNextMax)-Tuse)*P3)
		end

		



		
	end

	"""计算功率的状态转移矩阵"""
	function generateStatusMoveMatrix()

	end


end



# 计算状态转移函数
function 