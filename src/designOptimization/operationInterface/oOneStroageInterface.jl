

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
