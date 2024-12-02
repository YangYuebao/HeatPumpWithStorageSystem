
"""输入一定系统结构和工作参数,返回系统计算需要用到的各种向量"""
function generateSystemCoff(::HeatPumpStoragePressedWater;
	refrigerantHeatSupply::String = "water",  # 供热循环使用制冷剂
	refrigerantHeatStorage::String = "water", # 蓄热使用制冷剂
	eta_s::Real = 0.7,                        # 压缩机绝热效率
	Twastein::Real = 80.0,                    # 废气进入温度
	Twasteout::Real = 40.0,                   # 废气排出温度
	dTwaste::Real = 5.0,                      # 废气温度-蒸发器温度
	TwastCapacity::Real = 2,                  # 废热容量是工厂用热量的倍数
	Tair::Real = 25.0,                        # 外部环境温度
	dTair::Real = 5.0,                        # 外部环境温度-蒸发器温度
	Tuse::Real = 120.0,                       # 工厂使用温度
	dTuse::Real = 5.0,                        # 冷凝器温度-工厂使用温度
	dTstorageInput::Real = 5.0,               # 冷凝器温度-蓄热温度
    dTuseStandard::Real = 5.0,                # 工厂用热标准温差
    dTstorageStandard::Real = 5.0,            # 工厂蓄热标准温差
    dTstorageOutputStandard::Real = 5.0,      # 工厂蓄热释放标准温差
    TstorageTankMax::Real = 200.0,            # 蓄热罐的最高温度
	heatConsumptionPower::Real = 10000.0,     # 每小时用热功率kW
	heatStorageCapacity::Real = 20000.0,      # 蓄热量kWh(承压水蓄热)
	heatStorageOutEfficiency::Real = 0.00,    # 蓄热衰减系数K
	maxheatStorageInputHour::Real = 4,        # 蓄热充满时长
	kt::Real = 0.5,                           # 蓄热 \Delta T_2 / \Delta T_1
	TChangeToElec::Real = 140,                # 热泵蓄热温度上限
	workingHours::Int = 16,                   # 每日工作小时数
	workingStartHour::Int = 0,                # 生产开始时间
	PheatPumpMax::Real = 10000.0,             # 热泵最大功率kW
	hourlyTariff::Vector = fill(0.7, 24),     # 电价向量
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
	TeAirSource = Tair - dTair                      # 空气源蒸发器温度

	# 分别生成供热废热源、供热空气源、蓄热废热源、蓄热空气源4个不同循环的COP函数，一阶导数和二阶导数
	function COPTe_Tc(Te, Tc, refrigerant)
		h1 = CoolProp.PropsSI("H", "T", Te + 273.15, "Q", 1, refrigerant)
		s1 = CoolProp.PropsSI("S", "T", Te + 273.15, "Q", 1, refrigerant)
		p2 = CoolProp.PropsSI("P", "T", Tc + 273.15, "Q", 1, refrigerant)
		h2 = CoolProp.PropsSI("H", "S", s1, "P", p2, refrigerant)
		wt = (h2 - h1) / eta_s  # 理论压缩功等于绝热压缩功除以绝热效率
		h3 = CoolProp.PropsSI("H", "T", Tc + 273.15, "Q", 0, refrigerant)
		return (h1 - h3) / wt + 1 # 制热循环效率
	end

    # 供热废热源、供热空气源、蓄热废热源、蓄热空气源4个不同循环的COP函数
	#=
    COPSupplyWaste(T3) = T3+dTstorageInput < TChangeToElec ? COPTe_Tc(TeRecycle, T3 + dTuse, refrigerantHeatSupply) : 1.0
	COPSupplyAir(T3) = T3+dTstorageInput < TChangeToElec ? COPTe_Tc(TeAirSource, T3 + dTuse, refrigerantHeatSupply) : 1.0
	COPStorageWaste(T1) = T1+dTstorageInput < TChangeToElec ? COPTe_Tc(TeRecycle, T1 + dTstorageInput, refrigerantHeatStorage) : 1.0
	COPStorageAir(T1) = T1+dTstorageInput < TChangeToElec ? COPTe_Tc(TeAirSource, T1 + dTstorageInput, refrigerantHeatStorage) : 1.0
    =#

    COPSupplyWaste(T3) =COPTe_Tc(TeRecycle, T3 + dTuse, refrigerantHeatSupply)
	COPSupplyAir(T3) =COPTe_Tc(TeAirSource, T3 + dTuse, refrigerantHeatSupply)
	COPStorageWaste(T1) =COPTe_Tc(TeRecycle, T1 + dTstorageInput, refrigerantHeatStorage)
	COPStorageAir(T1) =COPTe_Tc(TeAirSource, T1 + dTstorageInput, refrigerantHeatStorage)



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
    T4=Tuse

    # 计算c_p*q_ml
    cpqml=heatConsumptionPower/dTuseStandard

    # 计算c_p*q_mh
    cpqmh=heatStorageVelocity/dTstorageStandard

    # 计算cpm
    cpm=heatStorageCapacity/(TstorageTankMax-Tuse+dTstorageOutputStandard)

    # 生成kt
    kt=kt

    # 生成废热容量倍数
    TwastCapacity=TwastCapacity

    # 生成蓄热罐最高温度
    TstorageTankMax=TstorageTankMax

	return 	COPSupplyWaste,COPSupplyAir,COPStorageWaste,COPStorageAir,
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
    TstorageTankMax
end

"""给定系统参数，求解系统成本，返回成本、热泵功率向量、蓄热量向量"""
function generateAndSolve(::HeatPumpStoragePressedWater, ::MinimizeCost;
	COPSupplyWaste::Function,
    COPSupplyAir::Function,
    COPStorageWaste::Function,
    COPStorageAir::Function,
	heatConsumptionPowerList::Vector{Float64},  # 用热负载向量
	heatStorageCapacityConstraint::Float64, # 蓄热量约束（最大值）
	heatpumpPowerConstraint::Float64,   # 热泵功率约束（最大值）
	hourlyTariffList::Vector{Float64},   # 电价向量
	heatStorageVelocity::Real,           # 蓄热速率约束
    T4::Real,
    cpqml::Real,
    cpqmh::Real,
    cpm::Real,
    kt::Real,
    TwastCapacity::Real,
    TstorageTankMax::Real
)
	model = Model(Ipopt.Optimizer)
	set_silent(model)

	# 建立热泵功率变量与储热量变量，直接约束了变量的范围
	m = 24
	@variable(model, 0 <= heatPumpPowerl1[1:m] <= heatpumpPowerConstraint)  # 供热回收
    @variable(model, 0 <= heatPumpPowerl2[1:m] <= heatpumpPowerConstraint)  # 供热空气
    @variable(model, 0 <= heatPumpPowerh1[1:m] <= heatpumpPowerConstraint)  # 蓄热回收
    @variable(model, 0 <= heatPumpPowerh2[1:m] <= heatpumpPowerConstraint)  # 蓄热空气
	@variable(model, 0 <= T[1:m] <= TstorageTankMax)
	@variable(model, 0 <= T1[1:m] <= TstorageTankMax)
    @variable(model, 0 <= T3[1:m] <= T4)
	@variable(model, 0 <= T5[1:m] <= T4)

    # T5
    @constraint(model, [i = 1:m], heatConsumptionPowerList[i] == cpqml*(T4-T5[i]))
    @constraint(model, [i = 1:m], T[i] <= T1[i])
    @constraint(model, [i = 1:m], T5[i] <= T3[i])
    
    # 热平衡
    @constraint(model, [i = 1:m-1], cpm*(T[i+1]-T[i])==cpqmh*2*kt/(1+kt)*(T1[i]-T[i])-cpqml*(T4-T3[i]))
    @constraint(model, cpm*(T[1]-T[m])==cpqmh*2*kt/(1+kt)*(T1[m]-T[m])-cpqml*(T4-T3[m]))

    @constraint(model, [i = 1:m], COPSupplyWaste(T3[i])*heatPumpPowerl1[i]+COPSupplyAir(T3[i])*heatPumpPowerl2[i]==cpqml*(T3[i]-T5[i]))
    @constraint(model, [i = 1:m], COPStorageWaste(T1[i])*heatPumpPowerh1[i]+COPStorageAir(T1[i])*heatPumpPowerh2[i]==cpqmh*2*kt/(1+kt)*(T3[i]-T5[i]))
    @constraint(model, [i = 1:m], COPSupplyWaste(T3[i])*heatPumpPowerl1[i]+COPStorageWaste(T1[i])*heatPumpPowerh1[i] <= heatConsumptionPowerList*TwastCapacity)
	
    # 蓄热速率约束
    @constraint(model, [i = 1:m], cpqmh*2*kt/(1+kt)*(T1[i]-T[i]) <= heatpumpPowerConstraint)

    # 目标函数
	@objective(model, Min, sum(hourlyTariffList[i] * (heatPumpPowerl1[i]+heatPumpPowerl2[i]+heatPumpPowerh1[i]+heatPumpPowerh2[i]) for i ∈ 1:m))
	optimize!(model)



	Pl1List = value.(heatPumpPowerl1)
    Pl2List = value.(heatPumpPowerl2)
    Ph1List = value.(heatPumpPowerh1)
    Ph2List = value.(heatPumpPowerh2)
    PList = Pl1List .+ Pl2List .+ Ph1List .+ Ph2List
	costList = PList .* hourlyTariffList
	heatStorageList = cpm*(value.(T).-T4)
	#=
	蓄热的增加和减少有待修改
	=#

	return Pl1List,Pl2List,Ph1List,Ph2List,PList,costList,heatStorageList
end


