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
struct PressedWaterDoubleStorageOneCompressor <: EnergySystem end# 蓄热系统，只有一个压缩机

export HeatPumpStoragePhaseChange,
	PressedWaterHighStorage,
	PressedWaterDoubleStorage,
	PressedWaterDoubleStorageSimplified,
	PressedWaterDoubleStorageOneCompressor

#双蓄系统生成COP函数
export getCOPFunction

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
	include(joinpath(dirname, file))
end

export generateSystemCoff, generateAndSolve

end # module HeatPumpWithStorageSystem
