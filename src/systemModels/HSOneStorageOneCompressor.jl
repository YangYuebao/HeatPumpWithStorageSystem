
abstract type FunctionGenerateMethod end
struct BackwardGenerate <: FunctionGenerateMethod end
struct LinearGenerate <: FunctionGenerateMethod end
struct BackwardGenerateCycled <: FunctionGenerateMethod end
struct LinearGenerateCycled <: FunctionGenerateMethod end

# 系统类型
abstract type SimplifiedType end
struct VaryLoadVaryArea <: SimplifiedType end	# 无简化
struct ConstloadandArea <: SimplifiedType end

# 计算最优的0点温度时需要考虑优化方法
abstract type OptimizeMethod end
struct ExhaustiveMethod <: OptimizeMethod end# 穷举法
struct GoldenRatioMethod <: OptimizeMethod end# 0.618法
struct MomentumMethod <: OptimizeMethod end# 动量法

function functionInterpolationGenerator(list::Vector, duration::Real, ::BackwardGenerate)
	n = length(list)
	dt = duration / n
	function f(x::Real)
		if x >= duration
			x = x % duration
		elseif x < 0
			x = x % duration + duration
		end
		return list[floor(Int, x / dt)+1]
	end
	return f
end
functionInterpolationGenerator(list::Vector, duration::Real, ::BackwardGenerateCycled) = functionInterpolationGenerator(list, duration, BackwardGenerate())

function functionInterpolationGenerator(list::Vector, duration::Real, ::LinearGenerate)
	n = length(list)
	if n == 1
		list = [list[1], list[1]]
		n=2
	end
	#n = 2
	tList = collect(range(0, stop = duration, length = n))
	sitpCOP = linear_interpolation(tList, list)
	function f(x::Real)
		if x >= duration
			x = x % duration
		elseif x < 0
			x = x % duration + duration
		end
		return sitpCOP(x)
	end
	return f
end

function functionInterpolationGenerator(list::Vector, duration::Real, ::LinearGenerateCycled)
	if abs(list[1] - list[end]) < 1e-6
		return functionInterpolationGenerator(list, duration, LinearGenerate())
	else
		@warn "列表首尾不一致，已经自动添加一位"
		return functionInterpolationGenerator(vcat(list, list[1]), duration, LinearGenerate())
	end
end


"""生成电价函数"""
function generateGridPriceFunction(gridPriceList::Vector, duration::Real)
	functionInterpolationGenerator(gridPriceList, duration, BackwardGenerate())
end

"""生成负载需求函数"""
function generateLoadFunction(loadList::Vector, duration::Real)
	functionInterpolationGenerator(loadList, duration, LinearGenerateCycled())
end

"""生成环境温度函数"""
function generateAreaTemperatureFunction(temperatureList::Vector, duration::Real)
	functionInterpolationGenerator(temperatureList, duration, LinearGenerateCycled())
end

"""
输入一定系统结构和工作参数,返回系统计算需要用到的各种向量
"""
function generateSystemCoff(::PressedWaterOneStorageOneCompressor;
	overlapRefrigerant::OverlapRefrigerant = NH3_Water,# 复叠工质
	COP1_design::Real = 2.1301295025490354,
	COPWater_design::Real = 3.0,
	maxTcHigh::Real = 180.0,  # 高温热泵冷凝器温度上限
	TCompressorIn::Real = 115.0,# 中间温度
	TWaste::Real = 60.0,                  # 废热源温度
	Tair::Vector = fill(25.0,24),          # 外部环境温度
	Tuse::Real = 120.0,                     # 工厂使用温度
	heatStorageCapacity::Real = 6.0,        # 蓄热量kWh(承压水蓄热)
	TstorageTankMax::Real = 220.0,          # 蓄热罐的最高温度
	maxheatStorageInputHour::Real = 4,      # 蓄热充满时长
	dT_EvaporationStandard::Real = 5.0,# 全蒸温差
	heatConsumptionPower::Vector = fill(1.0,24),       # 每小时用热功率kW
	workingStartHour::Int = 8,              # 生产开始时间
	workingHours::Int = 16,                 # 每日工作小时数
	heatPumpServiceCoff::Real = 1.2,        # 热泵功率服务系数
	hourlyTariff::Vector = fill(0.7, 24),   # 电价向量
)
	# COP函数外部生成
	# 计算工作时间
	temp = workingHours
	workingHours = workingHours % 24
	if workingHours == 0
		workingHours = 24
	end
	if temp > 24
		@warn "工作时长大于24小时,改为$(workingStartHour)h,终止计算"
		return -1
	end

	# 计算热泵基础功率和实际配置功率
	maxheatPower=maximum(heatConsumptionPower)
	P1base =  maxheatPower / COP1_design
	PheatPumpMax = P1base * heatPumpServiceCoff
	PelecHeatMax = maxheatPower * heatStorageCapacity / maxheatStorageInputHour
	PWaterCompressorMax = maxheatPower/COPWater_design

	# 生成需求与环境函数
	hourlyTariffFunction = generateGridPriceFunction(hourlyTariff, 24)
	heatConsumptionPowerFunction = generateLoadFunction(heatConsumptionPower, 24)
	TairFunction = generateAreaTemperatureFunction(Tair, 24)

	#=
	heatConsumptionPowerList = zeros(24)
	heatConsumptionPowerList[1:workingHours] .= heatConsumptionPower
	heatConsumptionPowerList = rotateVectorForward(heatConsumptionPowerList, workingStartHour)
	=#

	latentHeat = 2150.0# 汽化潜热kJ/kg
	cp_cw = 4.275# 循环水定压热容kJ/kg
	cp_cs = 2.281# 循环蒸汽定压热容kJ/kg


	# 生成低温热泵参数：cpqm_l,k1,dTe_l1,dTe_l2,dTc_l,QhRecycle
	cpm_h = heatStorageCapacity * maxheatPower / (TstorageTankMax - Tuse)# kWh/K
	storageTankMass = cpm_h / cp_cw * 3600#kg
	storageTankVolume = storageTankMass / 900#m^3

	#Qe_hStandard=1.0
	dT_EvaporationStandard = dT_EvaporationStandard

	# 设备运行约束:TstorageTankMax,heatStorageVelocity,heatStorageCapacityConstraint,heatpumpPowerConstraint
	heatStorageVelocity = heatStorageCapacity / maxheatStorageInputHour * heatConsumptionPower# kWh/h
	heatStorageCapacityConstraint = heatStorageCapacity * heatConsumptionPower# kWh
	heatpumpPowerConstraint = PheatPumpMax# kW


	# 其它参数:hourlyTariffList,heatConsumptionPowerList
	TcChangeToElec = maxTcHigh

	return hourlyTariffFunction, heatConsumptionPowerFunction, TairFunction,
	Tuse, TCompressorIn,
	dT_EvaporationStandard,
	latentHeat, cp_cw, cp_cs,
	TcChangeToElec, TWaste,
	cpm_h,
	TstorageTankMax, PheatPumpMax, PelecHeatMax,
	PWaterCompressorMax
end

# 导入定需求定环境的代码
include(joinpath(pwd(),"src","systemModels","HSOneStorageOneCompressor","ConstantLoadConstantArea.jl"))

# 导入变需求变环境的代码
include(joinpath(pwd(),"src","systemModels","HSOneStorageOneCompressor","VaryLoadVaryArea.jl"))



