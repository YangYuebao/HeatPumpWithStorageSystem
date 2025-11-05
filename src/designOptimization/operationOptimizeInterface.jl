#=
考虑一个设计优化问题
用热温度、电价曲线均确定
=#

"""设计优化参数结构体"""
abstract type DesignOptimizeInterface end

"""
    设计优化常量参数，用来描述工厂用热的设计条件
"""
struct DesignOptimizeInput <: DesignOptimizeInterface
    hourlyTariff :: Vector{Float64}
    Tair :: Vector{Float64}
    heatConsumptionPower :: Vector{Float64}
    maxCOP::Float64
    eta_s::Float64
    workingStartHour::Int
	workingHours::Int
    Tuse::Float64
	TWaste::Float64
	TCompressorIn::Float64
	maxTcHigh::Float64
	dT_EvaporationStandard::Float64
	Tsmin::Float64
	Tsmax::Float64
	dTRecycleSupply::Float64
    dTRecycleBackward::Float64

	# 计算参数
	dT::Float64
	dt::Float64
	smoother::Float64
    sysStruct::RecycleStruct
    refrigerant::OverlapRefrigerant
    DesignOptimizeInput(;
        hourlyTariff :: Vector{Float64},
        Tair :: Vector{Float64},
        heatConsumptionPower :: Vector{Float64},
        maxCOP::Float64,
        eta_s::Float64,
        workingStartHour::Int,
        workingHours::Int,
        Tuse::Float64,
        TWaste::Float64,
        TCompressorIn::Float64,
        maxTcHigh::Float64,
        dT_EvaporationStandard::Float64,
        Tsmin::Float64,
        Tsmax::Float64,
        dTRecycleSupply::Float64,
        dTRecycleBackward::Float64,
        dT::Float64,
        dt::Float64,
        smoother::Float64,
        sysStruct::RecycleStruct,
        refrigerant::OverlapRefrigerant
    )=new(
        hourlyTariff,
        Tair,
        heatConsumptionPower,
        maxCOP,
        eta_s,
        workingStartHour,
        workingHours,
        Tuse,
        TWaste,
        TCompressorIn,
        maxTcHigh,
        dT_EvaporationStandard,
        Tsmin,
        Tsmax,
        dTRecycleSupply,
        dTRecycleBackward,
        dT,
        dt,
        smoother,
        sysStruct,
        refrigerant
    )
end

struct DesignOptimizeVariables <: DesignOptimizeInterface
    heatPumpServiceCoff::Float64    # 热泵服务系数
    heatStorageCapacity::Float64    # 蓄热容量
    maxheatStorageInputHour::Float64    # 蓄热电加热储满时长
    DesignOptimizeVariables(;
        heatPumpServiceCoff::Float64,
        heatStorageCapacity::Float64,
        maxheatStorageInputHour::Float64
    )=new(
        heatPumpServiceCoff,
        heatStorageCapacity,
        maxheatStorageInputHour
    )
    DesignOptimizeVariables(x...) = new(x...)
end

struct DesignOptimizeParameters <: DesignOptimizeInterface
    COPOverlapFunction::Function
    COPLowFunction::Function
    hourlyTariffFunction::Function
    heatConsumptionPowerFunction::Function
    TairFunction::Function
    TWaste::Float64
    dT::Float64
    dt::Float64
    smoother::Float64

    # param中的变量
    ThMax::Float64
    Tuse::Float64
    dT_EvaporationStandard::Float64
    TCompressorIn::Float64
    #cpm::Float64   # 由输入确定
    COPWater::Function
    #PhMax::Float64 # 由输入确定
    #PeMax::Float64 # 由输入确定
    cp_cw::Float64
    latentHeat::Float64
    Tsmin::Float64
    Tsmax::Float64
    dTRecycleSupply::Float64
    dTRecycleBackward::Float64
    sysStruct::RecycleStruct
    DesignOptimizeParameters(;
        COPOverlapFunction::Function,
        COPLowFunction::Function,
        hourlyTariffFunction::Function,
        heatConsumptionPowerFunction::Function,
        TairFunction::Function,
        TWaste::Float64,
        dT::Float64,
        dt::Float64,
        smoother::Float64,
        ThMax::Float64,
        Tuse::Float64,
        dT_EvaporationStandard::Float64,
        TCompressorIn::Float64,
        COPWater::Function,
        cp_cw::Float64,
        latentHeat::Float64,
        Tsmin::Float64,
        Tsmax::Float64,
        dTRecycleSupply::Float64,
        dTRecycleBackward::Float64,
        sysStruct::RecycleStruct
    )=new(
        COPOverlapFunction,
        COPLowFunction,
        hourlyTariffFunction,
        heatConsumptionPowerFunction,
        TairFunction,
        TWaste,
        dT,
        dt,
        smoother,
        ThMax,
        Tuse,
        dT_EvaporationStandard,
        TCompressorIn,
        COPWater,
        cp_cw,
        latentHeat,
        Tsmin,
        Tsmax,
        dTRecycleSupply,
        dTRecycleBackward,
        sysStruct
    )
    DesignOptimizeParameters(x...) = new(x...)
end

"""
从DesignOptimizeInput生成DesignOptimizeParameters
"""
function generateDesignOptimizeParameters(input::DesignOptimizeInput)
    maxCOP = input.maxCOP
    eta_s = input.eta_s
    or = input.refrigerant
    TCompressorIn=input.TCompressorIn
    COPOverlapFunction = getOverlapCOP_fixMidTemperature(
		or,
		TCompressorIn + or.midTDifference / 2;
		maxCOP = maxCOP,# 最大COP
		eta_s = eta_s,# 绝热效率
		dT = 1.0,# 插值步长
	)
    COPLowFunction = getCOP(
		or.refrigerantLow.minTe,# 蒸发温度下限,这里是实际设计中的蒸发冷凝温度界限
		or.refrigerantLow.maxTe,# 蒸发温度上限
		or.refrigerantLow.minTc,# 冷凝温度下限
		or.refrigerantLow.maxTc,# 冷凝温度上限
		or.refrigerantLow,# 工质
		maxCOP,# 最大COP
		eta_s,# 绝热效率
		1.0,# 插值步长
	)
    hourlyTariffFunction = generateGridPriceFunction(input.hourlyTariff, 24)
	heatConsumptionPowerFunction = generateLoadFunction(input.heatConsumptionPower, 24)
	TairFunction = generateAreaTemperatureFunction(input.Tair, 24)
    TWaste = input.TWaste
    dT=input.dT
    dt=input.dt
    smoother=input.smoother
    ThMax=input.maxTcHigh
    Tuse=input.Tuse
    dT_EvaporationStandard=input.dT_EvaporationStandard
    
    COPWater = getCOP(
		TCompressorIn,# 蒸发温度下限,这里是实际设计中的蒸发冷凝温度界限
		ThMax,# 蒸发温度上限
		TCompressorIn,# 冷凝温度下限
		ThMax,# 冷凝温度上限
		refWater,# 工质
		maxCOP,# 最大COP
		eta_s,# 绝热效率
		1.0,#dT 插值步长
	)
    
	cp_cw = 4.275# 循环水定压热容kJ/kg
    latentHeat = 2150.0# 汽化潜热kJ/kg
    Tsmin = input.Tsmin
    Tsmax = input.Tsmax
    dTRecycleSupply = input.dTRecycleSupply
    dTRecycleBackward = input.dTRecycleBackward
    sysStruct = input.sysStruct
    return DesignOptimizeParameters(
        COPOverlapFunction,
        COPLowFunction,
        hourlyTariffFunction,
        heatConsumptionPowerFunction,
        TairFunction,
        TWaste,
        dT,
        dt,
        smoother,
        ThMax,
        Tuse,
        dT_EvaporationStandard,
        TCompressorIn,
        COPWater,
        cp_cw,
        latentHeat,
        Tsmin,
        Tsmax,
        dTRecycleSupply,
        dTRecycleBackward,
        sysStruct
    )
end

"""
生成设计优化的目标函数。
算法流程：
1. 输入设计变量
2. 结合设计常量，重新整合成优化问题
"""
function generateOperationFunction(designParameters::DesignOptimizeParameters,designInput::DesignOptimizeInput)
    function operationFunction(
        heatPumpServiceCoff::Float64,    # 热泵服务系数
        heatStorageCapacity::Float64,    # 蓄热容量
        maxheatStorageInputHour::Float64    # 蓄热电加热储满时长
    )
        COPOverlapFunction = designParameters.COPOverlapFunction
        COPLowFunction = designParameters.COPLowFunction
        hourlyTariffFunction = designParameters.hourlyTariffFunction
        heatConsumptionPowerFunction = designParameters.heatConsumptionPowerFunction
        TairFunction = designParameters.TairFunction
        TWaste = designParameters.TWaste
        dT=designParameters.dT
        dt=designParameters.dt
        smoother=designParameters.smoother
        Tuse = designParameters.Tuse
        TCompressorIn = designParameters.TCompressorIn
        maxheatPower = maximum(designInput.heatConsumptionPower)
        COPWater_design = designParameters.COPWater(TCompressorIn,Tuse)
        TstorageTankMax = designParameters.Tsmax
        PhMax = maxheatPower/COPWater_design * heatPumpServiceCoff
        PeMax = maxheatPower * (
            1- heatPumpServiceCoff + 
            heatStorageCapacity / maxheatStorageInputHour
        )
        cpm = heatStorageCapacity * maxheatPower / (TstorageTankMax - Tuse)

        params = SystemParameters(
            ThMax = designParameters.ThMax,
            Tuse = designParameters.Tuse,
            dT = designParameters.dT_EvaporationStandard,
            TCompressorIn = designParameters.TCompressorIn,
            cpm = cpm,#
            COPWater = designParameters.COPWater,
            PhMax = PhMax,#
            PeMax = PeMax,#
            cp_cw = designParameters.cp_cw,
            latentHeat = designParameters.latentHeat,
            Tsmin = designParameters.Tsmin,
            Tsmax = designParameters.Tsmax,
            dTRecycleSupply = designParameters.dTRecycleSupply,
            dTRecycleBackward = designParameters.dTRecycleBackward,
            sysStruct = designParameters.sysStruct
        )
        #=
        println("""
        ThMax: $(params.ThMax)
        Tuse: $(params.Tuse)
        dT: $(params.dT)
        TCompressorIn: $(params.TCompressorIn)
        cpm: $(params.cpm)
        PhMax: $(params.PhMax)
        PeMax: $(params.PeMax)
        Tsmin: $(params.Tsmin)
        Tsmax: $(params.Tsmax)
        dTRecycleSupply: $(params.dTRecycleSupply)
        dTRecycleBackward: $(params.dTRecycleBackward)
        sysStruct: $(params.sysStruct)
        """)
        =#

        return generateAndSolve(PressedWaterOneStorageOneCompressor(), MinimizeCost(), VaryLoadVaryArea(), GoldenRatioMethod();
            COPOverlap = COPOverlapFunction,
            COPLowFunction = COPLowFunction,
            hourlyTariffFunction = hourlyTariffFunction,
            heatConsumptionPowerFunction = heatConsumptionPowerFunction,
            TairFunction = TairFunction,

            params = params,

            TWaste = TWaste,
            # 求解参数
            dT = dT,# 状态参数高温蓄热温度离散步长
            dt = dt,# 时间步长
            #lambdaPe=lambdaPe,
            smoother = smoother,
        )
    end
    function operationFunction(designVariables::DesignOptimizeVariables)
        return operationFunction(
            designVariables.heatPumpServiceCoff,
            designVariables.heatStorageCapacity,
            designVariables.maxheatStorageInputHour
        )
    end

    return operationFunction
end

struct OperationOptimizeResult <: DesignOptimizeInterface
    TsList::Vector{Float64}
    tList::Vector{Float64}
    TairList::Vector{Float64}
    
    refrigerant::OverlapRefrigerant
end
