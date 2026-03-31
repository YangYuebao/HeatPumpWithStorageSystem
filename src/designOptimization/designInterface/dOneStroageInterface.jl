# 经济参数
struct OneStorageFinanceParameters <: AbstractFinanceParameters
	Pelec_cost::Float64      # 电锅炉价格 (¥/kWh)
	Storage_cost::Float64    # 蓄热设备价格 (¥/kWh)
	Life_years::Int64        # 设备寿命 (年)
	annual_days::Int64       # 年运行天数
	Discount_rate::Float64   # 折现率 
end
# 关键字构造函数
OneStorageFinanceParameters(;
	Pelec_cost::Float64,
	Storage_cost::Float64,
	Life_years::Int64,
	annual_days::Int64,
	Discount_rate::Float64,
) = OneStorageFinanceParameters(
	Pelec_cost,
	Storage_cost,
	Life_years,
	annual_days,
	Discount_rate,
)

# 总成本计算
function totalPresentWorth(
    ::OneStorage,
	financeParams::OneStorageFinanceParameters,# 经济参数
	operationResult::Float64,          # 每日运行成本（元）
	heatStorageCapacity::Float64,    # 蓄热容量
    PeMax::Float64    # 蓄热电加热功率
)
	# 初投资
	capitalCost = (
		PeMax * financeParams.Pelec_cost +
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
	return pw      # 单位：元
end

# 双层优化最低总成本
function get_bb_cost(::OneStorage,
	optimizeFunction::Function,
	# 给定经济参数结构体
	fp::OneStorageFinanceParameters = OneStorageFinanceParameters(
		Pelec_cost = 300.0,       # 电锅炉价格 (¥/kW)
		Storage_cost = 400.0,     # 蓄热设备价格 (¥/kWh)
		Life_years = 20,          # 设备寿命 (年)
		annual_days = 300,        # 年运行天数
		Discount_rate = 0.032,      # 折现率 
	)
)
	function bb_cost(x::Vector{Float64})
		heatStorageCapacity = x[1]    # 蓄热容量
		PeMax = x[2] # 蓄热电加热储满时长

		# 调用操作优化目标函数
		operationResult = optimizeFunction(
			heatStorageCapacity,
			PeMax,
		)
		if operationResult.objective < 1.0
			@warn "目标函数值:$(operationResult.objective),小于1.0，请检查参数:$(round.(x,digits=3))"
		end


		# 计算总现值
		totalCost = totalPresentWorth(
			OneStorage(),
			fp,
			operationResult.objective,   # 每日运行成本（元）
			heatStorageCapacity,
			PeMax,
		)
		return totalCost
	end
	return bb_cost
end