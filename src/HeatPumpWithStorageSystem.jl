module HeatPumpWithStorageSystem

using JuMP, HiGHS, Ipopt # 优化接口JuMP，优化器HiGHS
using LinearAlgebra
using Interpolations
using CoolProp# 物性库
using CSV,DataFrames

abstract type EnergySystem end# 计算系统形式

"""相变蓄热"""
struct HeatPumpStoragePhaseChange <: EnergySystem end# 热泵+蓄热系统
"""承压水蓄热，蓄热温度高"""
struct PressedWaterHighStorage <: EnergySystem end# 热泵+蓄热系统
"""承压水蓄热，中高温双蓄热"""
struct PressedWaterDoubleStorage <: EnergySystem end# 双蓄热系统
struct PressedWaterDoubleStorageSimplified <: EnergySystem end# 双蓄热系统，简化了管路
struct PressedWaterOneStorageOneCompressor <: EnergySystem end# 蓄热系统，只有一个压缩机

abstract type OOTest <: EnergySystem end
struct OO000 <: OOTest end
struct OO001 <: OOTest end
struct OO010 <: OOTest end
struct OO011 <: OOTest end
struct OO100 <: OOTest end
struct OO101 <: OOTest end
struct OO110 <: OOTest end
struct OO111 <: OOTest end



export HeatPumpStoragePhaseChange,
	PressedWaterHighStorage,
	PressedWaterDoubleStorage,
	PressedWaterDoubleStorageSimplified,
	PressedWaterDoubleStorageOneCompressor,
	PressedWaterOneStorageOneCompressor,
	OOTest,
  OO000, OO001, OO010, OO011, OO100, OO101

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

export Refrigerant,OverlapRefrigerant
export readCOPFile,generateCOPFile,getCOP,getCOP_g_h
export generateCOP,getOverlapCOP_fixMidTemperature

# 预设的制冷剂
export refR134a,refWater,refNH3,refR1233zdE
# 预设的复叠循环
export R134a_Water,NH3_Water,R1233zdE_Water

#双蓄系统生成COP函数
export getCOPFunction

# 古老的计算程序
export generateSystemCoff, generateAndSolve,getStateTransitionCost,getStateTransitionCost_SingleStep

# 根据向量生成函数，用于生成价格函数，负载函数，区域温度函数
export generateGridPriceFunction,generateLoadFunction,generateAreaTemperatureFunction

# 导出动态规划求解单压缩机系统的类型
export NoSimplify,ConstloadandArea,VaryLoadVaryArea

# 求解单压缩机系统的优化方法
export ExhaustiveMethod,GoldenRatioMethod,MomentumMethod
export GoldenRatioSolver


# 定义宏用于展开结构体字段
export @unpackParameters
# 混合整数线性规划导出系统参数结构体
export SystemParameters
export getMinimumCost,getMinimumCost_MILP	#导出
end # module HeatPumpWithStorageSystem
