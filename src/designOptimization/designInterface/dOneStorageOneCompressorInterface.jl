# 经济参数
struct FinanceParameters <: AbstractFinanceParameters
	Plow_cost::Float64       # 低温热泵价格 (¥/kWh)
	Phigh_cost::Float64      # 高温热泵价格 (¥/kWh)
	Pelec_cost::Float64      # 电锅炉价格 (¥/kWh)
	Storage_cost::Float64    # 蓄热设备价格 (¥/kWh)
	Life_years::Int64        # 设备寿命 (年)
	annual_days::Int64       # 年运行天数
	Discount_rate::Float64   # 折现率 
end
# 关键字构造函数
FinanceParameters(;
	Plow_cost::Float64,
	Phigh_cost::Float64,
	Pelec_cost::Float64,
	Storage_cost::Float64,
	Life_years::Int64,
	annual_days::Int64,
	Discount_rate::Float64,
) = FinanceParameters(
	Plow_cost,
	Phigh_cost,
	Pelec_cost,
	Storage_cost,
	Life_years,
	annual_days,
	Discount_rate,
)

# 总成本计算
function totalPresentWorth(
    ::PressedWaterOneStorageOneCompressor,
	financeParams::FinanceParameters,# 经济参数
	operationResult::Float64,          # 每日运行成本（元）
	heatPumpServiceCoff::Float64,    # 热泵服务系数
	heatStorageCapacity::Float64,    # 蓄热容量kwh
	maxheatStorageInputHour::Float64, # 蓄热电加热储满时长h
)
	# 初投资
	capitalCost = (
		heatPumpServiceCoff * (financeParams.Plow_cost + financeParams.Phigh_cost) +
		max(
            1- heatPumpServiceCoff + 
            heatStorageCapacity / maxheatStorageInputHour
        ,0.0) * financeParams.Pelec_cost +
		heatStorageCapacity * financeParams.Storage_cost
	)

	# 年运行费
	annualOperationCost = operationResult * financeParams.annual_days

	i = financeParams.Discount_rate
	n = financeParams.Life_years

	# 运行费折现系数
	annuity_pv_factor = (1 - (1 + i)^(-n)) / i

	# 总现值
	pw = capitalCost + annualOperationCost * annuity_pv_factor
	return pw,capitalCost,annuity_pv_factor
end

# 双层优化最低总成本
function get_bb_cost(
	::PressedWaterOneStorageOneCompressor,
	optimizeFunction::Function,
	# 给定经济参数结构体
	fp::FinanceParameters = FinanceParameters(
		Plow_cost = 1000.0,       # 低温热泵价格 (¥/kW)
		Phigh_cost = 1000.0,      # 高温热泵价格 (¥/kW)
		Pelec_cost = 300.0,       # 电锅炉价格 (¥/kW)
		Storage_cost = 400.0,     # 蓄热设备价格 (¥/kWh)
		Life_years = 20,          # 设备寿命 (年)
		annual_days = 300,        # 年运行天数
		Discount_rate = 0.032,      # 折现率 
	)
)
	function bb_cost(x::Vector{Float64})
		heatPumpServiceCoff = x[1]    # 热泵服务系数
		heatStorageCapacity = x[2]    # 蓄热容量
		maxheatStorageInputHour = x[3] # 蓄热电加热储满时长

		# 调用操作优化目标函数
		operationResult = optimizeFunction(
			heatPumpServiceCoff,
			heatStorageCapacity,
			maxheatStorageInputHour,
		)
		if operationResult[1] < 1.0
			@warn "目标函数值:$(operationResult[1]),小于1.0，请检查参数:$(round.(x,digits=3))"
		end


		# 计算总现值
		totalCost = totalPresentWorth(
			PressedWaterOneStorageOneCompressor(),
			fp,
			operationResult[1],   # 每日运行成本（元）
			heatPumpServiceCoff,
			heatStorageCapacity,
			maxheatStorageInputHour,
		)
		return totalCost[1]
	end
	return bb_cost
end