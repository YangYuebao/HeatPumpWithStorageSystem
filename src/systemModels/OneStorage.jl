
mutable struct OneStorageParameters
    heatStorageCapacity::Real
    PheaterMax::Real

    hourlyTariffFunction::Function
    timePoint::Vector
    loadFunction::Function
    heatStorageOutEfficiency::Real    # 蓄热释放效率
    heatStorageInEfficiency::Real      # 蓄热充能效率
    heatStorageVelocity::Real
end

struct OneStorageResult
    objective::Real
    heaterPower::Vector
    heatStorage::Vector
    heatStorageIn::Vector
    heatStorageOut::Vector
end

OneStorageParameters(;
    heatStorageCapacity::Real,
    PheaterMax::Real,
    hourlyTariffFunction::Function,
    timePoint::Vector,
    loadFunction::Function,
    heatStorageOutEfficiency::Real=0.95,    # 蓄热释放效率
    heatStorageInEfficiency::Real=0.99,      # 蓄热充能效率
    heatStorageVelocity::Real=1.0
)=OneStorageParameters(
    heatStorageCapacity,
    PheaterMax,
    hourlyTariffFunction,
    timePoint,
    loadFunction,
    heatStorageOutEfficiency,    # 蓄热释放效率
    heatStorageInEfficiency,      # 蓄热充能效率
    heatStorageVelocity
)

function Base.getproperty(r::OneStorageParameters, s::Symbol)
	if s === :hourlyTariff
		return r.hourlyTariffFunction.(r.timePoint)
	elseif s === :load
		return r.loadFunction.(r.timePoint)
    elseif s === :timePeriod
        return r.timePoint[2:end]-r.timePoint[1:end-1]
	else
		return getfield(r,s)
	end
end

function generateSystemCoff(::OneStorage;
        heatStorageCapacity::Real = 2.0,      # 蓄热量kWh(相变蓄热)
        PheaterMax::Real = 1.0,               # 热泵最大功率kW
        hourlyTariff::Vector = fill(0.7,24),      # 电价向量
        timePoint::Vector = 0:24,        # 
        load::Vector = fill(1.0,24),
        heatStorageOutEfficiency::Real = 0.95,    # 蓄热释放效率
        heatStorageInEfficiency::Real = 0.99,      # 蓄热充能效率
        heatStorageVelocity::Real = 1.0
    )

    return OneStorageParameters(
        heatStorageCapacity,
        PheaterMax,
        generateGridPriceFunction(hourlyTariff, 24),
        timePoint,
        generateLoadFunction(load, 24),
        heatStorageOutEfficiency,    # 蓄热释放效率
        heatStorageInEfficiency,      # 蓄热充能效率
        heatStorageVelocity
    )
end

"""输入一定系统结构和工作参数,返回系统计算需要用到的各种向量"""
function generateAndSolve(::OneStorage;
    heatStorageCapacity::Real=2.0,      # 蓄热量kWh(相变蓄热)
    PheaterMax::Real=1.0,               # 热泵最大功率kW
    hourlyTariff::Vector=fill(0.7,24),      # 电价向量
    timePeriod::Vector=fill(1.0,24),        # 
    load::Vector=fill(1.0,24),
    heatStorageOutEfficiency::Real=0.95,    # 蓄热释放效率
    heatStorageInEfficiency::Real=0.99,      # 蓄热充能效率
    heatStorageVelocity::Real=1.0
)
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    m=0
    if length(timePeriod) == length(hourlyTariff) == length(load)
        m=length(timePeriod)
    else
        @error("length not match")
        println("""
        length(timePeriod) = $(length(timePeriod))
        length(hourlyTariff) = $(length(hourlyTariff))
        length(load) = $(length(load))
        """)
    end
    
    @variable(model, 0<=heaterPower[1:m]<=PheaterMax)
    @variable(model, 0<=heatStorage[1:m+1]<=heatStorageCapacity)
    @variable(model, 0<=heatStorageIn[1:m]<=heatStorageVelocity)
    @variable(model, 0<=heatStorageOut[1:m])

    # 蓄热速率约束
    @constraint(model,[i=1:m], heatStorage[i+1]==heatStorage[i]+heatStorageIn[i]-heatStorageOut[i])
    @constraint(model,heatStorage[m+1]==heatStorage[1])

    # 热平衡约束
    @constraint(model,[i=1:m], heaterPower[i] + (heatStorageOut[i] * heatStorageOutEfficiency - heatStorageIn[i]/heatStorageInEfficiency)/timePeriod[i] == load[i])

    # 目标函数
    @objective(model, Min, sum(hourlyTariff[i]*heaterPower[i]*timePeriod[i] for i=1:m))
    optimize!(model)

    isFeasible = primal_status(model)
    
    flag = isFeasible in [FEASIBLE_POINT,NEARLY_FEASIBLE_POINT]

    if !flag
        return OneStorageResult(
            9999.0,
            fill(9999.0,m),
            fill(9999.0,m+1),
            fill(9999.0,m),
            fill(9999.0,m)
        )
    end

    return OneStorageResult(objective_value(model),value.(heaterPower),value.(heatStorage),value.(heatStorageIn),value.(heatStorageOut))
end

function generateAndSolve(::OneStorage,osp::OneStorageParameters)
    return generateAndSolve(OneStorage();
        heatStorageCapacity=osp.heatStorageCapacity,      # 蓄热量kWh(相变蓄热)
        PheaterMax=osp.PheaterMax,               # 热泵最大功率kW
        hourlyTariff=osp.hourlyTariff[1:end-1],      # 电价向量
        timePeriod=osp.timePeriod,        # 
        load=osp.load[1:end-1],
        heatStorageOutEfficiency=osp.heatStorageOutEfficiency,    # 蓄热释放效率
        heatStorageInEfficiency= osp.heatStorageInEfficiency,      # 蓄热充能效率
        heatStorageVelocity=osp.heatStorageVelocity,
    )
end
