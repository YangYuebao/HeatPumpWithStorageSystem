module HeatPumpWithStorageSystem



using JuMP, HiGHS, Ipopt, COPT # 优化接口JuMP，优化器HiGHS
using LinearAlgebra
using Interpolations
using CoolProp# 物性库
using CSV,DataFrames
using BlackBoxOptim

abstract type EnergySystem end# 计算系统形式

"""相变蓄热"""
struct HeatPumpStoragePhaseChange <: EnergySystem end# 热泵+蓄热系统
"""承压水蓄热，蓄热温度高"""
struct PressedWaterHighStorage <: EnergySystem end# 热泵+蓄热系统
"""承压水蓄热，中高温双蓄热"""
struct PressedWaterDoubleStorage <: EnergySystem end# 双蓄热系统
struct PressedWaterDoubleStorageSimplified <: EnergySystem end# 双蓄热系统，简化了管路
struct PressedWaterOneStorageOneCompressor <: EnergySystem end# 蓄热系统，只有一个压缩机
struct PressedWaterOneStorageOneCompressor_MILP <: EnergySystem end# 蓄热系统，只有一个压缩机
struct OneStorage <: EnergySystem end# 电锅炉加蓄热


export HeatPumpStoragePhaseChange,
	PressedWaterHighStorage,
	PressedWaterDoubleStorage,
	PressedWaterDoubleStorageSimplified,
	PressedWaterDoubleStorageOneCompressor,
	PressedWaterOneStorageOneCompressor,
	#PressedWaterOneStorageOneCompressor_MILP,
	OneStorage,
	SystemStructure,RecycleStruct

abstract type SystemObjectiveType end
struct MinimizeCost <: SystemObjectiveType end# 最小化成本
struct MinimizeEnergyCost <: SystemObjectiveType end# 最小化能耗

export MinimizeCost, MinimizeEnergyCost

#=
dirname = joinpath(pwd(), "src", "refrigerantPropertys")
for file in readdir(dirname)
	include(joinpath(dirname, file))
end
=#
dirname = joinpath(pwd(), "src", "systemModels")
for file in readdir(dirname)
	srcpath = joinpath(dirname, file)
	if isfile(srcpath)
		include(srcpath)
	end
end

include(joinpath(pwd(), "src", "DPSolver", "DPSolverCore.jl"))
using .DPSolverCore
const DC = DPSolverCore

export Refrigerant,OverlapRefrigerant
export readCOPFile,generateCOPFile,getCOP,getCOP_g_h
export generateCOP,getOverlapCOP_fixMidTemperature

# 预设的制冷剂
export refR134a,refWater,refNH3,refR1233zdE
# 预设的复叠循环
export R134a_Water,NH3_Water,R1233zdE_Water

#双蓄系统生成COP函数
export getCOPFunction

# 运行优化计算程序
export generateSystemCoff, generateAndSolve,getStateTransitionCost,getStateTransitionCost_SingleStep,getTemperatureLineCost

# 根据向量生成函数，用于生成价格函数，负载函数，区域温度函数
export generateGridPriceFunction,generateLoadFunction,generateAreaTemperatureFunction

# 导出动态规划求解单压缩机系统的类型
export NoSimplify,ConstloadandArea,VaryLoadVaryArea

# 求解单压缩机系统的优化方法
export ExhaustiveMethod,GoldenRatioMethod,MomentumMethod
export ExhaustiveSolver,GoldenRatioSolver


# 定义宏用于展开结构体字段
export @unpackParameters
# 混合整数线性规划导出系统参数结构体
export SystemParameters,SystemVariables
export getMinimumCost,getMinimumCost_MILP	#导出
export getCOPbyMode
#export MILPModelParameters, MILPModelResult
#export computeCOPdiscretization

"""双层优化问题的经济性参数"""
abstract type AbstractFinanceParameters end

"""设计优化参数结构体"""
abstract type DesignOptimizeInterface end

dirname = joinpath(pwd(), "src", "designOptimization","designInterface")
for file in readdir(dirname)
	srcpath = joinpath(dirname, file)
    println("design include:",file)
	if isfile(srcpath)
		include(srcpath)
	end
end

dirname = joinpath(pwd(), "src", "designOptimization","operationInterface")

for file in readdir(dirname)
	srcpath = joinpath(dirname, file)
    println("operation include:",file)
	if isfile(srcpath)
		include(srcpath)
	end
end

# DP精确解模块依赖 DesignOptimizeParameters（定义在 operationInterface 中），
# 因此必须在 operationInterface 之后加载
include(joinpath(pwd(), "src","systemModels","HSOneStorageOneCompressorMILP","DP_PreciseSolution.jl"))

# 导出用来与设计优化联合的类型和函数
export DesignOptimizeInput,DesignOptimizeVariables,DesignOptimizeParameters
export generateDesignOptimizeParameters,generateOperationFunction

export FinanceParameters,OneStorageFinanceParameters
#,MILPFinanceParameters
export totalPresentWorth, get_bb_cost,bboptimize

end # module HeatPumpWithStorageSystem
