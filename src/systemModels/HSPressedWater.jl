
"""
向量循环向后
"""
function rotateVector(vec::Vector,k::Int)
    k=k%length(vec)
    return vcat(vec[end-k+1:end],vec[1:end-k])        
end

"""输入一定系统结构和工作参数,返回系统计算需要用到的各种向量"""
function generateSystemCoff(::HeatPumpStoragePressedWater;
    refrigerant::String="water",            # 制冷剂
    eta_s::Real=0.7,                        # 压缩机绝热效率
    Twastein::Real=80.0,                    # 废气进入温度
    Twasteout::Real=40.0,                   # 废气排出温度
    dTwaste::Real=5.0,                      # 废气温度-蒸发器温度
    Tair::Real=25.0,                        # 外部环境温度
    dTair::Real=5.0,                        # 外部环境温度-蒸发器温度
    Tuse::Real=180.0,                       # 工厂使用温度
    dTuse::Real=5.0,                        # 冷凝器温度-工厂使用温度
    dTstorageInput::Real=5.0,               # 冷凝器温度-蓄热温度
    heatConsumptionPower::Real=10000.0,     # 每小时用热功率kW
    heatStorageCapacity::Real=20000.0,      # 蓄热量kWh(相变蓄热)
    workingHours::Int=16,                   # 每日工作小时数
    workingStartHour::Int=0,                # 生产开始时间
    PheatPumpMax::Real=10000.0,             # 热泵最大功率kW
    hourlyTariff::Vector=fill(0.7,24),      # 电价向量
    heatStorageOutEfficiency::Real=0.98,    # 蓄热释放效率
    heatStorageInEfficiency::Real=1.0,      # 蓄热充能效率
    maxheatStorageInputHour::Real = 4,      # 蓄热充满时长
    )
    workingStartHour = workingStartHour%24
    # 计算一天的热泵蒸发冷凝温度
    Te=zeros(24)
    Tc=zeros(24)
    Te[1:workingHours] .= (Twastein+Twasteout)/2-dTwaste   # 蒸发器温度是废热平均温度减换热温差
    Tc[1:workingHours] .= Tuse+dTuse                    # 冷凝器温度是工厂使用温度加上用热温差
    Te[workingHours+1:24] .= Tair-dTair                 # 剩余时间蒸发器用外环境温度减换热温差
    Tc[workingHours+1:24] .= Tuse+dTstorageInput+dTuse  # 剩余时间冷凝器用工厂使用温度加上蓄热温差和用热温差

    Tc = rotateVector(Tc,workingStartHour)
    Te = rotateVector(Te,workingStartHour)

    # 计算一天的热泵COP
    if workingHours >24
        @error "工作时长大于24小时"
        return -1
    end
    COP=zeros(24)
    for i=1:24
        # CoolProp.PropsSI("H","T",Te+273.15,"Q",1,rf) - CoolProp.PropsSI("H","T",Tc+273.15,"Q",0,rf) 
        # 用热物性计算给定蒸发冷凝温度下单质工质的循环效率
        h1 = CoolProp.PropsSI("H","T",Te[i]+273.15,"Q",1,refrigerant)
        s1 = CoolProp.PropsSI("S","T",Te[i]+273.15,"Q",1,refrigerant)
        p2 = CoolProp.PropsSI("P","T",Tc[i]+273.15,"Q",1,refrigerant)
        h2 = CoolProp.PropsSI("H","S",s1,"P",p2,refrigerant)
        wt = (h2-h1)/eta_s  # 理论压缩功等于绝热压缩功除以绝热效率
        h3 = CoolProp.PropsSI("H","T",Tc[i]+273.15,"Q",0,refrigerant)
        COP[i] = (h1-h3)/wt + 1 # 制热循环效率
    end

    # 生成用热负载向量
    heatConsumptionPowerList = zeros(24)
    heatConsumptionPowerList[1:workingHours] .= heatConsumptionPower

    heatConsumptionPowerList=rotateVector(heatConsumptionPowerList,workingStartHour)

    # 生成蓄热量约束
    heatStorageCapacityConstraint = heatStorageCapacity

    # 生成泵功率约束
    heatpumpPowerConstraint = PheatPumpMax

    # 生成充放热效率
    heatStorageOutEfficiency=heatStorageOutEfficiency
    heatStorageInEfficiency = heatStorageInEfficiency

    # 生成电价向量
    hourlyTariffList = hourlyTariff

    # 生成充热速率
    heatStorageVelocity = heatStorageCapacity/maxheatStorageInputHour

    return COP,heatConsumptionPowerList,heatStorageCapacityConstraint,heatpumpPowerConstraint,hourlyTariffList,heatStorageOutEfficiency,heatStorageInEfficiency,heatStorageVelocity
end

"""给定系统参数，求解系统成本，返回成本、热泵功率向量、蓄热量向量"""
function generateAndSolve(::HeatPumpStoragePressedWater,::MinimizeCost;
    simulationDays::Int=7,  # 为使程序计算进入一般的工况，需要多仿真一段时间
    COP::Vector{Float64},   # 热泵COP向量
    heatConsumptionPowerList::Vector{Float64},  # 用热负载向量
    heatStorageCapacityConstraint::Float64, # 蓄热量约束（最大值）
    heatpumpPowerConstraint::Float64,   # 热泵功率约束（最大值）
    hourlyTariffList::Vector{Float64},   # 电价向量
    heatStorageOutEfficiency::Real,     # 蓄热释放效率
    heatStorageInEfficiency::Real,      # 蓄热输入效率
    heatStorageVelocity::Real           # 蓄热速率约束
    )
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    # 把一天的数据扩展到simulationDays
    COP=repeat(COP,simulationDays)
    heatConsumptionPowerList=repeat(heatConsumptionPowerList,simulationDays)
    hourlyTariffList=repeat(hourlyTariffList,simulationDays)

    # 建立热泵功率变量与储热量变量，直接约束了变量的范围
    m = 24*simulationDays
    @variable(model, 0<=heatPumpPower[1:m]<=heatpumpPowerConstraint)
    @variable(model, 0<=heatStorage[1:m]<=heatStorageCapacityConstraint)
    @variable(model, 0<=heatStorageIn[1:m-1]<=heatStorageVelocity)
    @variable(model, 0<=heatStorageOut[1:m-1])

    # 蓄热速率约束
    @constraint(model,[i=1:m-1], heatStorage[i+1]==heatStorage[i]+heatStorageIn[i]-heatStorageOut[i])

    # 热平衡约束
    @constraint(model,[i=1:m-1], heatPumpPower[i]*COP[i] + heatStorageOut[i] * heatStorageOutEfficiency == heatConsumptionPowerList[i] + heatStorageIn[i]/heatStorageInEfficiency)

    # 目标函数
    @objective(model, Min, sum(hourlyTariffList[i]*heatPumpPower[i] for i=1:m))
    optimize!(model)

    k=floor(Int,simulationDays/2)
    
    PList = value.(heatPumpPower)
    costList = PList.*hourlyTariffList
    heatStorageList = value.(heatStorage)
    #=
    蓄热的增加和减少有待修改
    =#
    heatStorageOutput = value.(heatStorageOut)*heatStorageOutEfficiency
    heatStorageInput = value.(heatStorageIn)
    heatpumpOutput = PList[k*24+1:(k+1)*24].*COP[1:24]
    
    return  costList[k*24+1:(k+1)*24],PList[k*24+1:(k+1)*24],heatStorageList[k*24+1:(k+1)*24],heatStorageOutput,heatpumpOutput,COP[k*24+1:(k+1)*24]
end


