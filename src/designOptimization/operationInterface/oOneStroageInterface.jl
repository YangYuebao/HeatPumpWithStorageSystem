
"""
从DesignOptimizeInput生成DesignOptimizeParameters
"""
function generateDesignOptimizeParameters(::OneStorage,input::DesignOptimizeInput)
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
生成电加热蓄热设计优化的目标函数。
算法流程：
1. 输入设计变量
2. 结合设计常量，重新整合成优化问题
"""
function generateOperationFunction(::OneStorage,osp::OneStorageParameters)
    function operationFunction(
        heatStorageCapacity::Float64,    # 蓄热容量
        PeMax::Float64    # 蓄热电加热功率
    )
        osp.heatStorageCapacity = heatStorageCapacity
        osp.PheaterMax = PeMax

        return generateAndSolve(OneStorage(),osp)
    end

    return operationFunction
end
