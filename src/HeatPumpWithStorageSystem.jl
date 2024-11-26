module HeatPumpWithStorageSystem

using JuMP, HiGHS# 优化接口JuMP，优化器HiGHS
using CoolProp# 物性库

abstract type EnergySystem end# 计算系统形式
struct HeatPumpStoragePhaseChange <: EnergySystem end# 热泵+蓄热系统
struct HeatPumpStoragePressedWater <: EnergySystem end# 热泵+蓄热系统

export HeatPumpStoragePhaseChange,
	HeatPumpStoragePressedWater

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
