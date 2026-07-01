"""
    PressedWaterOneStorageOneCompressor_MILP 操作优化接口
    生成运行优化的目标函数，供双层优化调用
"""

"""
    生成操作优化函数
    输入：设计变量 -> 输出：每日运行成本
"""
function generateOperationFunction(
    ::PressedWaterOneStorageOneCompressor_MILP,
    designParameters::DesignOptimizeParameters,
    designInput::DesignOptimizeInput
)
    local call_count = 0
    local controlResult = []

    function operationFunction(
        heatPumpServiceCoff::Float64,    # 热泵服务系数
        heatStorageCapacity::Float64,    # 蓄热容量 kWh
        maxheatStorageInputHour::Float64 # 蓄热电加热储满时长 h
    )
        call_count += 1
        println("MILP调用$(call_count):$(round.([heatPumpServiceCoff, heatStorageCapacity, maxheatStorageInputHour], digits=3))")

        Tuse = designParameters.Tuse
        TCompressorIn = designParameters.TCompressorIn
        maxheatPower = maximum(designInput.heatConsumptionPower)
        COPWater_design = designParameters.COPWater(TCompressorIn, Tuse)

        # 计算设备容量
        PhMax = maxheatPower / COPWater_design * heatPumpServiceCoff
        PeMax = max(maxheatPower * (
            1 - heatPumpServiceCoff +
            heatStorageCapacity / maxheatStorageInputHour
        ), 0.0)
        cpm = heatStorageCapacity * maxheatPower / (designParameters.Tsmax - Tuse)

        # 生成时段数据
        n = length(designInput.hourlyTariff)
        dt = designParameters.dt

        # 生成外部输入曲线
        hourlyTariff = [designParameters.hourlyTariffFunction(t) for t in 0:dt:(n*dt-dt)]
        heatLoad = [designParameters.heatConsumptionPowerFunction(t) for t in 0:dt:(n*dt-dt)]
        Tair = [designParameters.TairFunction(t) for t in 0:dt:(n*dt-dt)]

        # COP分档预计算
        Ts_range = [designParameters.Tsmin, designParameters.Tsmax]
        cop_data = computeCOPdiscretization(
            designParameters.COPOverlapFunction,
            designParameters.COPLowFunction,
            designParameters.COPWater,
            designParameters.COPOverlapFunction,  # COP3 also uses overlap COP
            Ts_range,
            Tair,
            Tuse,
            TCompressorIn,
            n,
            10  # 默认分10档
        )

        # 构建MILP参数
        milp_params = MILPModelParameters(
            Tuse,
            designParameters.Tsmin,
            designParameters.Tsmax,
            designParameters.ThMax,
            designParameters.dT_EvaporationStandard,
            TCompressorIn,
            cpm,
            dt,
            COPWater_design,
            COPWater_design,  # COPdesign
            designParameters.smoother,
            designParameters.smoother,
            hourlyTariff,
            heatLoad,
            Tair,
            cop_data.m1, cop_data.COP1v, cop_data.T1g,
            cop_data.m2, cop_data.COP2v, cop_data.T2g,
            cop_data.m3, cop_data.COP3v, cop_data.T3g,
            cop_data.mw, cop_data.COPwv, cop_data.Twg,
            cop_data.COP_ca
        )

        result = generateAndSolve(PressedWaterOneStorageOneCompressor_MILP(), milp_params)

        controlResult = result
        return result
    end

    function operationFunction(designVariables::DesignOptimizeVariables)
        return operationFunction(
            designVariables.heatPumpServiceCoff,
            designVariables.heatStorageCapacity,
            designVariables.maxheatStorageInputHour
        )
    end

    function getParams(
        heatPumpServiceCoff::Float64,
        heatStorageCapacity::Float64,
        maxheatStorageInputHour::Float64
    )
        Tuse = designParameters.Tuse
        TCompressorIn = designParameters.TCompressorIn
        maxheatPower = maximum(designInput.heatConsumptionPower)
        COPWater_design = designParameters.COPWater(TCompressorIn, Tuse)

        PhMax = maxheatPower / COPWater_design * heatPumpServiceCoff
        PeMax = max(maxheatPower * (
            1 - heatPumpServiceCoff +
            heatStorageCapacity / maxheatStorageInputHour
        ), 0.0)
        cpm = heatStorageCapacity * maxheatPower / (designParameters.Tsmax - Tuse)

        return SystemParameters(
            ThMax = designParameters.ThMax,
            Tuse = designParameters.Tuse,
            dT = designParameters.dT_EvaporationStandard,
            TCompressorIn = TCompressorIn,
            cpm = cpm,
            COPWater = designParameters.COPWater,
            PhMax = PhMax,
            PeMax = PeMax,
            cp_cw = designParameters.cp_cw,
            latentHeat = designParameters.latentHeat,
            Tsmin = designParameters.Tsmin,
            Tsmax = designParameters.Tsmax,
            dTRecycleSupply = designParameters.dTRecycleSupply,
            dTRecycleBackward = designParameters.dTRecycleBackward,
            sysStruct = designParameters.sysStruct
        )
    end

    return operationFunction, () -> (call_count), () -> (controlResult), getParams
end