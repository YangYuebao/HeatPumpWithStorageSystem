"""
动态规划初始解生成参数结构体

字段说明:
- n: 状态系数量，即温度系的数量
- q: 单位热份数，用于计算温度变化量（每个温度系中相邻状态间放热恰好提供一份热量Q_unit）
- solver_type: 求解器类型，:Exhaustive（穷举搜索）或 :GoldenRatio（黄金分割搜索）
- dt: 动态规划每层间的时间间隔 (h)
"""
struct DP_INITIAL_PARAMS
    n::Int
    q::Int
    solver_type::Symbol
	dt::Float64
	y_s6::Float64	# 两次s6最小间隔 (h)，0表示无约束（仅 cop_mode=:continuous 生效）
end
# 兼容旧4参数构造（y_s6 默认0 = 无约束）
DP_INITIAL_PARAMS(n::Int, q::Int, solver_type::Symbol, dt::Float64) =
    DP_INITIAL_PARAMS(n, q, solver_type, dt, 0.0)

"""
系统变量结构体，封装各时段的热负荷和电价
"""
struct SysVariables
    load::Vector{Float64}
    price::Vector{Float64}
	epsilon::Float64	# DP平滑参数 
	SysVariables(load::Vector{Float64}, price::Vector{Float64}, epsilon::Float64)=new(load, price, epsilon)
end
SysVariables(load::Vector{Float64}, price::Vector{Float64})=SysVariables(load, price, 1e-4)

"""
获取温度状态空间（自适应探测模式）

参数:
- T_list: 初始温度列表，从这些温度开始探测
- q: 单位热份数，用于计算温度变化量
- params: MILPModelParameters 结构体，包含系统参数

返回:
- T_state_space: 字典，键是初始温度，值是从该温度出发可达的温度状态空间

算法逻辑:
1. 对每个初始温度T，从T开始向下探测（降温）
2. 根据COP分档和单位热份数计算温度变化量
3. 当温度变化跨越COP分档边界时记录温度点
4. 再从T开始向上探测（升温）
5. 返回每个初始温度对应的温度状态空间

参考: VaryLoadVaryArea.jl 中温度状态空间生成逻辑
"""
function get_temperature_state_space(
    #T_list::Vector{Float64},
	n::Int,	#状态系数量
    q::Int,
    params::MILPModelParameters,
)
    Q_unit = 1/q

	T_use = params.Tuse
	dTs = params.dTs
    Ts_min = params.Tsmin
    Ts_max = params.Tsmax
    COP2v = params.COP2v[:,1] |> Vector{Float64}
    T2g = params.T2g[1, :] |> Vector{Float64}
    cpm = params.cpm

    T_state_space = Dict{Float64, Vector{Float64}}()
	T_list = [Ts_min]

    idx = 1	# 记录序号，对于第一个温度，要计算它的下一个温度，在这之间划分状态系
	while idx <= n
    #for T in T_list
		T = T_list[idx]
        T_state_space[T] = [T]
        # 先往下探
        current_T = T
        next_idx = findfirst(x -> x > current_T, T2g)-1
        next_T = T
        while next_T >= Ts_min && next_idx >=1
            next_T = current_T - Q_unit/cpm*(COP2v[next_idx]-1)/(COP2v[next_idx])
            if next_T >= T2g[next_idx]
                next_T = min(next_T, T2g[next_idx+1])
                pushfirst!(T_state_space[T], next_T)
                current_T = next_T
                next_idx = findfirst(x -> x > current_T, T2g)-1
            else
                next_idx -= 1
            end
        end

        # 再往上探
        current_T = T
        current_idx = findfirst(x -> x > current_T, T2g)-1
        next_T = T
        while next_T <= Ts_max && current_idx <= length(T2g)-1
            next_T = current_T + Q_unit/cpm*(COP2v[current_idx]-1)/(COP2v[current_idx])
            if next_T <= Ts_max
                push!(T_state_space[T], next_T)
                current_T = next_T
                current_idx = findfirst(x -> x > current_T, T2g)-1
            else
                break
            end
        end
        
		if idx == 1
			#T_gap = (T_state_space[T][2]-T_state_space[T][1])/(n)
            sort!(T_state_space[T])
			T_list = collect(range(Ts_min, T_state_space[T][2], length=n+1))[1:end-1]
		end
		push!(T_state_space[T], Ts_max-dTs,T_use-dTs,T_use+dTs)
		sort!(T_state_space[T])
		idx +=1
    end

    return T_state_space
end

"""
单步变量结构体，封装计算单步转移成本所需的所有变量

字段说明:
- heat_load: 当前时段热负荷 (kW)
- COP_ca: 热泵直接供热模式 COP（设计工况）
- COP1v: 模式1 COP（直接供热，随温度变化）
- COP2v: 模式2 COP（蓄热供热，随温度变化）
- COP3v: 模式3 COP（热泵储热，随温度变化）
- COPwv: 水蒸气压缩机 COP（随温度变化）
- COP_lambda: 热泵模式切换时的COP（函数类型，输入温度返回COP）
- dt: 时间步长 (h)
- initialize: 是否为初始化阶段
- price: 当前时段电价 (元/kWh)
"""
struct SingleStepVariables
    heat_load::Float64           # 当前时段热负荷 (kW)
    COP_ca::Float64              # 热泵直接供热模式 COP
    COP1v::Float64               # 模式1 COP（直接供热）
    COP2v::Float64               # 模式2 COP（蓄热供热）
    COP3v::Float64               # 模式3 COP（热泵储热）
    COPwv::Float64               # 水蒸气压缩机 COP
    COP_lambda::Function         # 热泵模式切换时的COP（函数类型）
    dt::Float64                  # 时间步长 (h)
    initialize::Bool             # 是否为初始化阶段
    price::Float64               # 当前时段电价 (元/kWh)
end

"""
精确求解 s7 状态（热泵储热转热泵供热）的 JuMP 模型

参数:
- COP_lambda: 切换点处的 COP 值（Float64，已在切换温度处求值）
- COP_ca: 热泵直接供热模式 COP
- cpm: 蓄热罐热容质量乘积 (kJ/℃)
- dt: 时间步长 (h)
- Ts_start: 起始蓄热温度 (℃)
- T_lambda: 切换点温度 (℃)
- Ts_end: 结束蓄热温度 (℃)
- load: 热负荷 (kW)
- C_heatpump: 热泵额定容量 (kW)
- C_boiler: 电锅炉额定容量 (kW)

返回:
- obj: 目标函数值（总电功率），不可行时返回9999.0
- P3s_lambda, Pesl_lambda, P1e_lambda, Pees_lambda, Peel_lambda, lambda: 各功率变量值

状态描述:
s7 状态表示在前半段时间内热泵向蓄热罐储热，电锅炉承担热负荷；
后半段时间内热泵直接向工厂供热，电锅炉承担储热需求。

参考文档: 热泵运行调度算法_混合整数线性规划_SOS1.md 中 s7 状态定义
"""
function accurate_s7_model(
	COP_lambda::Float64,
	COP_ca::Float64,
	cpm::Float64,
	dt::Float64,
	Ts_start::Float64,
	T_lambda::Float64,
	Ts_end::Float64,
	load::Float64,
	C_heatpump::Float64,
	C_boiler::Float64,
)
    model = direct_model(COPT.Optimizer())
    set_silent(model)

    @variable(model, P3s_lambda >= 0.0)    # 前半段热泵储热功率
    @variable(model, Pesl_lambda >= 0.0)   # 前半段电锅炉供热功率
    @variable(model, Pess_lambda >= 0.0)   # 前半段电锅炉储热功率
    @variable(model, P1e_lambda >= 0.0)    # 后半段热泵直接供热功率
    @variable(model, Pees_lambda >= 0.0)   # 后半段电锅炉储热功率
    @variable(model, Peel_lambda >= 0.0)   # 后半段电锅炉供热功率
    @variable(model, 0 <= lambda <= 1.0)   # 切换时间比例

    @constraint(model, P3s_lambda * COP_lambda + Pess_lambda == cpm / dt * (T_lambda - Ts_start))
    @constraint(model, Pesl_lambda == load * lambda)
    @constraint(model, P1e_lambda * COP_ca + Peel_lambda == load * (1 - lambda))
    @constraint(model, Pees_lambda == cpm / dt * (Ts_end - T_lambda))
    @constraint(model, P1e_lambda <= (1 - lambda) * C_heatpump / COP_ca)
    @constraint(model, P3s_lambda <= lambda * C_heatpump / COP_lambda)
    @constraint(model, Pesl_lambda + Pess_lambda <= lambda * C_boiler)
    @constraint(model, Pees_lambda + Peel_lambda <= (1 - lambda) * C_boiler)

    @objective(model, Min, P3s_lambda + Pesl_lambda + Pess_lambda + P1e_lambda + Pees_lambda + Peel_lambda)

    optimize!(model)

    if primal_status(model) in [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]
        return objective_value(model), value.([P3s_lambda, Pesl_lambda, Pess_lambda, P1e_lambda, Pees_lambda, Peel_lambda, lambda])...
    else
        return (9999.0, fill(9999.0, 7)...)
    end
end

"""
精确求解 s8 状态（热泵同时供热和储热转热泵供热）的 JuMP 模型

参数:
- COP_lambda: 切换点处的 COP 值（Float64）
- COP_ca: 热泵直接供热模式 COP
- cpm: 蓄热罐热容质量乘积 (kJ/℃)
- dt: 时间步长 (h)
- Ts_start: 起始蓄热温度 (℃)
- T_lambda: 切换点温度 (℃)
- Ts_end: 结束蓄热温度 (℃)
- load: 热负荷 (kW)
- C_heatpump: 热泵额定容量 (kW)
- C_boiler: 电锅炉额定容量 (kW)

返回:
- obj: 目标函数值（总电功率），不可行时返回9999.0
- P3s_lambda, P1s_lambda, Pesl_lambda, P1e_lambda, Pees_lambda, Peel_lambda, lambda: 各功率变量值

状态描述:
s8 状态表示在前半段时间内热泵同时向蓄热罐储热和向工厂供热；
后半段时间内蓄热完成，热泵仅向工厂供热，电锅炉承担储热需求。

参考文档: 热泵运行调度算法_混合整数线性规划_SOS1.md 中 s8 状态定义
"""
function accurate_s8_model(
	COP_lambda::Float64,
	COP_ca::Float64,
	cpm::Float64,
	dt::Float64,
	Ts_start::Float64,
	T_lambda::Float64,
	Ts_end::Float64,
	load::Float64,
	C_heatpump::Float64,
	C_boiler::Float64,
)
    model = direct_model(COPT.Optimizer())
    set_silent(model)

    @variable(model, P3s_lambda >= 0.0)    # 前半段热泵储热功率
    @variable(model, P1s_lambda >= 0.0)    # 前半段热泵直接供热功率
    @variable(model, Pesl_lambda >= 0.0)   # 前半段电锅炉供热功率
    @variable(model, Pess_lambda >= 0.0)   # 前半段电锅炉储热功率
    @variable(model, P1e_lambda >= 0.0)    # 后半段热泵直接供热功率
    @variable(model, Pees_lambda >= 0.0)   # 后半段电锅炉储热功率
    @variable(model, Peel_lambda >= 0.0)   # 后半段电锅炉供热功率
    @variable(model, 0 <= lambda <= 1.0)   # 切换时间比例

    @constraint(model, P3s_lambda * COP_lambda + Pess_lambda == cpm / dt * (T_lambda - Ts_start))
    @constraint(model, P1s_lambda * COP_lambda + Pesl_lambda == load * lambda)
    @constraint(model, P1e_lambda * COP_ca + Peel_lambda == load * (1 - lambda))
    @constraint(model, Pees_lambda == cpm / dt * (Ts_end - T_lambda))
    @constraint(model, P1e_lambda <= (1 - lambda) * C_heatpump / COP_ca)
    @constraint(model, P3s_lambda + P1s_lambda <= lambda * C_heatpump / COP_lambda)
    @constraint(model, Pesl_lambda + Pess_lambda <= lambda * C_boiler)
    @constraint(model, Pees_lambda + Peel_lambda <= (1 - lambda) * C_boiler)

    @objective(model, Min, P3s_lambda + P1s_lambda + Pesl_lambda + P1e_lambda + Pees_lambda + Peel_lambda)

    optimize!(model)

    if primal_status(model) in [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]
        return objective_value(model), value.([P3s_lambda, P1s_lambda, Pesl_lambda,Pess_lambda, P1e_lambda, Pees_lambda, Peel_lambda, lambda])...
    else
        return (9999.0, fill(9999.0, 8)...)
    end
end

"""
单步计算的最小单元，计算从一个温度状态到另一个温度状态的最优功率分配

参数:
- Ts_start: 起始蓄热温度 (℃)
- Ts_end: 结束蓄热温度 (℃)
- variables: SingleStepVariables 结构体，包含当前时段的系统变量
- params: MILPModelParameters 结构体，包含系统参数

返回:
- best_state: 最优运行状态编号 (1-8)
- best_power: 最优总电功率消耗 (kW)
- P1s_best, P2s_best, P3s_best: 各模式起始功率 (kW)
- Pesl_best, Pess_best: 电锅炉供热/储热起始功率 (kW)
- P1e_best, P2e_best, P3e_best: 各模式结束功率 (kW)
- Peel_best, Pees_best: 电锅炉补充结束功率 (kW)
- lambda_best: 状态切换时间比例

状态定义:
1: 热泵直接供热（蓄热温度不变或升高，热泵承担热负荷，电锅炉补充）
2: 蓄热供热（水蒸气压缩机驱动，蓄热温度降低）
3: 热泵向蓄热储热（蓄热温度升高，热泵承担储热，电锅炉承担热负荷）
4: 热泵和蓄热同时供热（蓄热温度降低，热泵和水蒸气压缩机同时工作）
5: 热泵同时向产线供热和向蓄热罐储热（蓄热温度升高，热泵同时承担两种负荷）
6: 蓄热用完后切换为热泵供热（蓄热温度降低至低于用热温度，前半段蓄热供热，后半段热泵供热）
7: 储热完成后热泵向工厂供热（蓄热温度升高，前半段热泵储热，后半段热泵供热）
8: 热泵同时供热和储热，蓄热完成后热泵供热（蓄热温度升高，前半段同时供热储热，后半段仅供热）

参考文档: 热泵运行调度算法_混合整数线性规划_SOS1.md 中状态定义
"""
function single_step(
    Ts_start::Float64,
    Ts_end::Float64,
    variables::SingleStepVariables,
    params::MILPModelParameters
)
    heat_load = variables.heat_load
    COP_ca = variables.COP_ca
    COP1v = variables.COP1v
    COP2v = variables.COP2v
    COP3v = variables.COP3v
    COPwv = variables.COPwv
    initialize = variables.initialize
    dt = variables.dt
    COP_lambda = variables.COP_lambda

    C_heatpump = params.C_heatpump
    C_boiler = params.C_boiler
    Tcmax = params.Tcmax
    dTs = params.dTs
    Tuse = params.Tuse
    cpm = params.cpm

    best_power = 9999.0
    best_state = 0
    P1s_best = 9999.0
    P2s_best = 9999.0
    P3s_best = 9999.0
    Pesl_best = 9999.0
    Pess_best = 9999.0
    P1e_best = 9999.0
    P2e_best = 9999.0
    P3e_best = 9999.0
    Peel_best = 9999.0
    Pees_best = 9999.0
    lambda_best = 0.0

	current_power = 9999.0
	P1s = 9999.0
	P2s = 9999.0
	P3s = 9999.0
	Pesl = 9999.0
	Pess = 9999.0
	P1e = 9999.0
	P2e = 9999.0
	P3e = 9999.0
	Peel = 9999.0
	Pees = 9999.0
	lambda = 0.0	# 状态切换时间比例

	# 计算中间变量
	storage_load = cpm*(Ts_end-Ts_start)/dt

	# s1 热泵供热模式
	if  storage_load >= 0	# 忽略自然散热，这部分很小
		P1s = min(heat_load,C_heatpump)/COP_ca
		Pesl = heat_load - P1s * COP_ca
		Pess = storage_load
		current_power = P1s + Pesl + Pess
		if current_power < best_power && C_boiler >= Pess + Pesl >= 0
			best_power = current_power
			best_state = 1
			P1s_best = P1s
			P2s_best = 0.0
			P3s_best = 0.0
			Pesl_best = Pesl
			Pess_best = Pess
			P1e_best = 0.0
			P2e_best = 0.0
			P3e_best = 0.0
			Peel_best = 0.0
			Pees_best = 0.0
			lambda_best = 1.0
		end
	end
	# s2 蓄热供热模式
	if storage_load < 0
		P2s = -storage_load/(COP2v-1)
		Pesl = heat_load - P2s * COP2v
		current_power = P2s + Pesl
		if current_power < best_power && C_boiler >= Pesl >= 0
			best_power = current_power
			best_state = 2
			P1s_best = 0.0
			P2s_best = P2s
			P3s_best = 0.0
			Pesl_best = Pesl
			Pess_best = 0.0
			P1e_best = 0.0
			P2e_best = 0.0
			P3e_best = 0.0
			Peel_best = 0.0
			Pees_best = 0.0
			lambda_best = 1.0
		end
	end

	# s3
	if storage_load > 0 && Ts_end <= Tcmax - dTs
		P3s = min(C_heatpump,storage_load)/COP3v
		Pess = storage_load - P3s * COP3v
		Pesl = heat_load
		lambda = 1.0
		current_power = P3s + Pesl + Pess
		if current_power < best_power && C_boiler >= Pess + Pesl >= 0
			best_power = current_power
			best_state = 3
			P1s_best = 0.0
			P2s_best = 0.0
			P3s_best = P3s
			Pesl_best = Pesl
			Pess_best = Pess
			P1e_best = 0.0
			P2e_best = 0.0
			P3e_best = 0.0
			Peel_best = 0.0
			Pees_best = 0.0
			lambda_best = lambda
		end
	end

	
	# s4
	# s4状态的经济结果和s6完全不同，只是在操作上有所不同，所以计算合并到s6去

	# s5
	if storage_load > 0 && (heat_load < C_heatpump || storage_load < C_heatpump) && Ts_end <= Tcmax - dTs
		P3s = storage_load/COP3v	# 这个必须小于C_heatpump，否则退化到s1或者s3
		if P3s * COP3v < C_heatpump	# COP3v = COP1v
			P1s = min(heat_load,C_heatpump - P3s * COP3v)/COP_ca
			Pesl = heat_load - P1s * COP_ca 
			lambda = 1.0
			current_power = P1s + P3s + Pesl
			if current_power < best_power && C_boiler >= Pesl >= 0
				best_power = current_power
				best_state = 5
				P1s_best = P1s
				P2s_best = 0.0
				P3s_best = P3s
				Pesl_best = Pesl
				Pess_best = 0.0
				P1e_best = 0.0
				P2e_best = 0.0
				P3e_best = 0.0
				Peel_best = 0.0
				Pees_best = 0.0
				lambda_best = lambda
			end
		end
	end

	# s6
	if storage_load < 0 && !initialize
		P2s = -storage_load/(COP2v-1)
		flag = P2s * COP2v <= C_heatpump	# 蓄热供热功率不能超过水蒸气压缩机的额定供热功率，C_heatpump是最大供热量
		P1e = min(heat_load-P2s*COP2v,C_heatpump) / COP_ca
		lambda = P2s*COP2v/heat_load
		Pesl = lambda * heat_load - P2s * COP2v
		Peel = (1-lambda) * heat_load - P1e * COP_ca
		current_power = P1e + P2s + Pesl + Peel
		
		if (0<lambda <1 && 
			current_power < best_power &&
			C_boiler >= Pesl + Peel >= 0 && P1e >= 0 && flag)
			best_power = current_power
			current_state = (Ts_end >= Tuse + dTs)&& !initialize ? 4 : 6
			if current_state == 4
				# 模式4同时供热
				best_state = 4
				P1s_best = P1e
				P2s_best = P2s
				P3s_best = 0.0
				Pesl_best = Pesl+Peel
				Pess_best = 0.0
				P1e_best = 0.0
				P2e_best = 0.0
				P3e_best = 0.0
				Peel_best = 0.0
				Pees_best = 0.0
				lambda_best = 1.0
			elseif current_state == 6
				best_state = 6
				P1s_best = 0.0
				P2s_best = P2s
				P3s_best = 0.0
				Pesl_best = Pesl
				Pess_best = 0.0
				P1e_best = P1e	# 4和6计算结果的经济性是一样的，但数值上6的后半段取P1s的值
				P2e_best = 0.0
				P3e_best = 0.0
				Peel_best = Peel
				Pees_best = 0.0
				lambda_best = lambda
			end
		end
	end
    #=
	# s7 热泵蓄热转热泵供热
	if storage_load > 0 && (storage_load < C_heatpump || heat_load < C_heatpump) && !initialize
		if initialize
			if storage_load < C_heatpump && Ts_end <= Tcmax - dTs
				P3s = min(C_heatpump, storage_load)/COP3v
				P1e = min(heat_load, C_heatpump)/COP_ca
				Pesl = heat_load
				Pess = storage_load - P3s * COP3v
				Peel = heat_load - P1e * COP_ca	# 简化后没法算出来lambda，直接合并到一起
				current_power = P1e + P3s + Pesl + Pess + Peel
				if current_power < best_power && C_boiler >= Pesl + Peel + Pess >= 0
					best_power = current_power
					best_state = 7
					P1s_best = 0.0
					P2s_best = 0.0
					P3s_best = P3s
					Pesl_best = Pesl
					Pess_best = Pess
					P1e_best = P1e
					P2e_best = 0.0
					P3e_best = 0.0
					Peel_best = Peel
					Pees_best = 0.0
					lambda_best = 1.0
				end
			end
		else
			#=s7和x1+x3(s5)的区别是，x1+x3的COP与x3->x1的COP不同=#
			obj = 9999.0
			P3s_lambda = Pesl_lambda = P1e_lambda = Pees_lambda = Peel_lambda = 9999.0
			Pess_lambda = 0.0
			Ts_lambda = Ts_start
			lambda = 0.0

			# 遍历不同的Ts_lambda，找到最优解
			Ts_lambda_list = range(Ts_start, min(Ts_end, Tcmax - dTs), length=10)
			for T_lambda_temp in Ts_lambda_list
				COP_lambda_v = COP_lambda(T_lambda_temp)
				obj_temp, P3s_lambda_temp, Pesl_lambda_temp, P1e_lambda_temp, Pees_lambda_temp, Peel_lambda_temp, lambda_temp = accurate_s7_model(
					COP_lambda_v,
					COP_ca,
					cpm,
					dt,
					Ts_start,
					T_lambda_temp,
					Ts_end,
					heat_load,
					C_heatpump,
					C_boiler
				)
				if obj_temp < obj
					obj = obj_temp
					P3s_lambda = P3s_lambda_temp
					Pesl_lambda = Pesl_lambda_temp
					P1e_lambda = P1e_lambda_temp
					Pees_lambda = Pees_lambda_temp
					Peel_lambda = Peel_lambda_temp
					Ts_lambda = T_lambda_temp
					lambda = lambda_temp
				end
			end

			current_power = P3s_lambda + Pess_lambda + Pesl_lambda + P1e_lambda + Pees_lambda + Peel_lambda
			if current_power < best_power
				best_power = current_power
				best_state = 7
				P1s_best = 0.0
				P2s_best = 0.0
				P3s_best = P3s_lambda
				Pesl_best = Pesl_lambda
				Pess_best = Pess_lambda
				P1e_best = P1e_lambda
				P2e_best = 0.0
				P3e_best = 0.0
				Peel_best = Peel_lambda
				Pees_best = Pees_lambda
				lambda_best = lambda
			end
		end
	end

	# s8
	if storage_load > 0 && C_heatpump > heat_load && !initialize
		if !initialize
			#=
			if storage_load < C_heatpump && Ts_end <= Tcmax - dTs
				P3s = storage_load/COP3v
				P1s = min(heat_load,C_heatpump - P3s * COP3v)/COP_3v
				Pes = heat_load - P1s * COP_3v

				Pes = heat_load - P1e * COP_ca
				current_power = P1s + P2s + P3s + Pes
				if current_power < best_power && C_boiler >= Pes >= 0
					best_power = current_power
					best_state = 7
					P1s_best = 0.0
					P2s_best = 0.0
					P3s_best = P3s
					Pes_best = Pes
					P1e_best = P1e
					P2e_best = 0.0
					P3e_best = 0.0
					Pee_best = 0.0
					lambda_best = 1.0
				end
			end
			=#
			obj = 9999.0
			P3s_lambda = P1s_lambda = Pesl_lambda = P1e_lambda = Pees_lambda = Peel_lambda = 9999.0
			Pess_lambda = 0.0   # s8前半段电锅炉不储热，储热由热泵完成
			Ts_lambda = Ts_start
			lambda = 0.0

			# 遍历不同的Ts_lambda，找到最优解
			Ts_lambda_list = range(Ts_start, min(Ts_end, Tcmax - dTs), length=10)
			for T_lambda_temp in Ts_lambda_list
				COP_lambda_v = COP_lambda(T_lambda_temp)
				obj_temp, P3s_lambda_temp, P1s_lambda_temp, Pesl_lambda_temp, P1e_lambda_temp, Pees_lambda_temp, Peel_lambda_temp, lambda_temp = accurate_s8_model(
					COP_lambda_v,
					COP_ca,
					cpm,
					dt,
					Ts_start,
					T_lambda_temp,
					Ts_end,
					heat_load,
					C_heatpump,
					C_boiler
				)
				if obj_temp < obj
					obj = obj_temp
					P3s_lambda = P3s_lambda_temp
					P1s_lambda = P1s_lambda_temp
					Pesl_lambda = Pesl_lambda_temp
					P1e_lambda = P1e_lambda_temp
					Pees_lambda = Pees_lambda_temp
					Peel_lambda = Peel_lambda_temp
					Ts_lambda = T_lambda_temp
					lambda = lambda_temp
				end
			end

			current_power = P3s_lambda + P1s_lambda + Pess_lambda + Pesl_lambda + P1e_lambda + Pees_lambda + Peel_lambda
			if current_power < best_power
				best_power = current_power
				best_state = 8
				P1s_best = P1s_lambda
				P2s_best = 0.0
				P3s_best = P3s_lambda
				Pesl_best = Pesl_lambda
				Pess_best = Pess_lambda
				P1e_best = P1e_lambda
				P2e_best = 0.0
				P3e_best = 0.0
				Peel_best = Peel_lambda
				Pees_best = Pees_lambda
				lambda_best = lambda
			end
		end
	end
	=#
	return best_state, best_power, P1s_best, P2s_best, P3s_best, Pesl_best, Pess_best, P1e_best, P2e_best, P3e_best, Peel_best, Pees_best, lambda_best
end

"""
计算温度状态转移的运行成本（动态规划核心函数）

该函数是 single_step 的封装层，负责根据当前时段的温度获取对应的 COP 参数，
然后调用 single_step 计算最优功率分配。

参数:
- Ts_start: 起始蓄热温度 (℃)
- Ts_end: 结束蓄热温度 (℃)
- params: MILPModelParameters 结构体
- sysVariables: 当前时段的系统变量（热负荷、电价等）
- t: 当前时段索引
- dt: 时间步长 (h)

返回:
- C: 总电功率消耗 (kW)，不可行时返回9999.0（不乘电价和时间步长）
- state: 最优运行状态编号 (1-8)
- P1s, P2s, P3s: 各模式起始功率 (kW)
- Pesl, Pess: 电锅炉供热/储热起始功率 (kW)
- P1e, P2e, P3e: 各模式结束功率 (kW)
- Peel, Pees: 电锅炉供热/储热结束功率 (kW)
- lambda: 状态切换时间比例

注意: 返回的 C 是纯功率值（kW），不乘以电价和时间步长。
电价在 DP 求解前统一乘以成本张量，因为不同电价但相同功率和蒸发温度的时段，
功率的状态转移成本是相同的，这样可以避免重复计算。

参考: VaryLoadVaryArea.jl 中 getStateTransitionCost_SingleStep 函数
"""
function calculateTransitionCost(Ts_start::Float64, Ts_end::Float64,
    params::MILPModelParameters, sysVariables, t::Int, dt::Float64)

    heat_load = sysVariables.load[t]

    m1 = params.m1
    m2 = params.m2
    m3 = params.m3
    mw = params.mw

    Tg1 = params.T1g[t, :]
    Tg2 = params.T2g[t, :]
    Tg3 = params.T3g[t, :]
    Tgw = params.Twg[t, :]

    COP1v_t = params.COP1v[:, t]
    COP2v_t = params.COP2v[:, t]
    COP3v_t = params.COP3v[:, t]
    COPwv_t = params.COPwv[:, t]

    # 使用时段结束时温度对应的COP（隐式格式）
    # 这与MILP模型中COP取时段末值的处理方式一致
    seg1 = findCOPSegment(Ts_end, Tg1)
    seg2 = findCOPSegment(Ts_end, Tg2)
    seg3 = findCOPSegment(Ts_end, Tg3)
    segw = findCOPSegment(Ts_end, Tgw)

    COP_ca = params.COPca[t]
    COP1v = COP1v_t[seg1]
    COP2v = COP2v_t[seg2]
    COP3v = COP3v_t[seg3]
    COPwv = COPwv_t[segw]

    function COP_lambda(T_lambda::Float64)
        seg_lambda = findCOPSegment(T_lambda, Tg1)
        return COP1v_t[seg_lambda]
    end

    variables = SingleStepVariables(
        heat_load,
        COP_ca,
        COP1v,
        COP2v,
        COP3v,
        COPwv,
        COP_lambda,
        dt,
        false,
        0.0
    )

    state, C, P1s, P2s, P3s, Pesl, Pess, P1e, P2e, P3e, Peel, Pees, lambda = single_step(
        Ts_start, Ts_end, variables, params
    )

    return C, state, P1s, P2s, P3s, Pesl, Pess, P1e, P2e, P3e, Peel, Pees, lambda
end

"""
构建状态转移成本张量（动态规划核心函数）

参数:
- Ts_list: 温度状态列表 [nT]
- nt: 时间段数
- dt: 时间步长 (h)
- params: MILPModelParameters 结构体
- sysVariables: 各时段的系统变量（包含热负荷、电价等）
- use_cache: 是否使用工况缓存，默认为 true（仅在计算初始解时使用）

返回:
- C: 状态转移功率成本张量 [nt, nT, nT]，C[t, i, j] 表示第t时段从温度i转移到温度j的电功率(kW)
- stateMatrix: 最优状态编号矩阵 [nt, nT, nT]
- P1sMatrix, P2sMatrix, P3sMatrix: 各模式起始功率矩阵 [nt, nT, nT]
- P1eMatrix, P2eMatrix, P3eMatrix: 各模式结束功率矩阵 [nt, nT, nT]
- PeslMatrix, PessMatrix: 电锅炉供热/储热起始功率矩阵 [nt, nT, nT]
- PeelMatrix, PeesMatrix: 电锅炉供热/储热结束功率矩阵 [nt, nT, nT]
- lambdaMatrix: 状态切换时间比例矩阵 [nt, nT, nT]

算法逻辑:
1. 初始化成本张量为 Inf（不可行）
2. 遍历每个时段
3. 如果启用缓存，先检查 (Tair, load) 是否在缓存中
   - 命中缓存：直接复制缓存的矩阵数据
   - 未命中缓存：计算并存入缓存
4. 未启用缓存时，直接计算每个状态转移的成本

缓存机制:
- 以 (Tair, load) 为 key 缓存整个时段的功率矩阵
- 当多个时段的蒸发温度和热负荷相同时，避免重复计算
- 适用于温度状态空间一致的场景（如初始化计算）

注意: 返回的 C 是纯功率值（kW），不乘以电价和时间步长。
电价在 DP 求解前统一乘以成本张量，因为不同电价但相同功率和蒸发温度的时段，
功率的状态转移成本是相同的，这样可以避免重复计算。

参考: VaryLoadVaryArea.jl 中 getStateTransitionCost 函数
"""
function buildCostTensor(Ts_list::Vector{Float64}, nt::Int, dt::Float64,
    params::MILPModelParameters, sysVariables; use_cache::Bool=true)

    nT = length(Ts_list)

    C = fill(9999.0, nt, nT, nT)

    stateMatrix = zeros(Int, nt, nT, nT)

    P1sMatrix = zeros(nt, nT, nT)
    P2sMatrix = zeros(nt, nT, nT)
    P3sMatrix = zeros(nt, nT, nT)
    P1eMatrix = zeros(nt, nT, nT)
    P2eMatrix = zeros(nt, nT, nT)
    P3eMatrix = zeros(nt, nT, nT)
    PeslMatrix = zeros(nt, nT, nT)
    PessMatrix = zeros(nt, nT, nT)
    PeelMatrix = zeros(nt, nT, nT)
    PeesMatrix = zeros(nt, nT, nT)
    lambdaMatrix = zeros(nt, nT, nT)

    # 工况缓存：key = (Tair, load)，value = 该工况下所有矩阵数据
    cache = Dict{Tuple{Float64, Float64}, NTuple{12, Matrix{Float64}}}()
    state_cache = Dict{Tuple{Float64, Float64}, Matrix{Int}}()

    for t in 1:nt
        Tair_t = params.Tair[t]
        load_t = sysVariables.load[t]
        cache_key = (Tair_t, load_t)

        if use_cache && haskey(cache, cache_key)
            # 命中缓存：直接复制
            cached_C, cached_P1s, cached_P2s, cached_P3s,
            cached_Pesl, cached_Pess, cached_P1e, cached_P2e,
            cached_P3e, cached_Peel, cached_Pees, cached_lambda = cache[cache_key]

            C[t, :, :] .= cached_C
            stateMatrix[t, :, :] .= state_cache[cache_key]
            P1sMatrix[t, :, :] .= cached_P1s
            P2sMatrix[t, :, :] .= cached_P2s
            P3sMatrix[t, :, :] .= cached_P3s
            P1eMatrix[t, :, :] .= cached_P1e
            P2eMatrix[t, :, :] .= cached_P2e
            P3eMatrix[t, :, :] .= cached_P3e
            PeslMatrix[t, :, :] .= cached_Pesl
            PessMatrix[t, :, :] .= cached_Pess
            PeelMatrix[t, :, :] .= cached_Peel
            PeesMatrix[t, :, :] .= cached_Pees
            lambdaMatrix[t, :, :] .= cached_lambda
        else
            # 未命中缓存：计算该时段所有状态转移
            for i in 1:nT
                Ts_start = Ts_list[i]

                for j in 1:nT
                    Ts_end = Ts_list[j]

                    power_cost, state, P1s, P2s, P3s, Pesl, Pess, P1e, P2e, P3e, Peel, Pees, lambda = calculateTransitionCost(
                        Ts_start, Ts_end, params, sysVariables, t, dt
                    )

                    if power_cost < 9999.0
                        C[t, i, j] = power_cost
                        stateMatrix[t, i, j] = state
                        P1sMatrix[t, i, j] = P1s
                        P2sMatrix[t, i, j] = P2s
                        P3sMatrix[t, i, j] = P3s
                        P1eMatrix[t, i, j] = P1e
                        P2eMatrix[t, i, j] = P2e
                        P3eMatrix[t, i, j] = P3e
                        PeslMatrix[t, i, j] = Pesl
                        PessMatrix[t, i, j] = Pess
                        PeelMatrix[t, i, j] = Peel
                        PeesMatrix[t, i, j] = Pees
                        lambdaMatrix[t, i, j] = lambda
                    end
                end
            end

            # 存入缓存（如果启用缓存）
            if use_cache
                cache[cache_key] = (
                    C[t, :, :],
                    P1sMatrix[t, :, :],
                    P2sMatrix[t, :, :],
                    P3sMatrix[t, :, :],
                    PeslMatrix[t, :, :],
                    PessMatrix[t, :, :],
                    P1eMatrix[t, :, :],
                    P2eMatrix[t, :, :],
                    P3eMatrix[t, :, :],
                    PeelMatrix[t, :, :],
                    PeesMatrix[t, :, :],
                    lambdaMatrix[t, :, :]
                )
                state_cache[cache_key] = stateMatrix[t, :, :]
            end
        end
    end

    return C, stateMatrix,
           P1sMatrix, P2sMatrix, P3sMatrix,
           P1eMatrix, P2eMatrix, P3eMatrix,
           PeslMatrix, PessMatrix, PeelMatrix, PeesMatrix,
           lambdaMatrix
end

"""
使用动态规划生成 MILP 初始解

参数:
- dp_params: DP_INITIAL_PARAMS 结构体，包含动态规划求解参数
- params: MILPModelParameters 结构体，包含系统参数
- sysVariables: 各时段的系统变量列表

返回:
- initial_solution: InitialSolution 结构体，包含初始解的所有变量

算法逻辑:
1. 生成温度状态空间（自适应探测模式，返回多个温度系）
2. 对每个温度系分别求解动态规划，选择运行成本最低的温度系
3. 根据最优温度序列提取功率值和状态信息
4. 构建初始解结构体

关键设计:
- 温度状态空间按照"相邻状态间放热恰好提供一份热量Q_unit"的原则生成
- 每个温度系中相邻温度的间隔由COP2v（随温度非线性变化）决定
- 对多个温度系分别求解DP，选择最优结果
- 功率成本张量 C[t, i, j] 存储的是纯电功率(kW)，不乘电价
- 在 DP 求解前，根据各时段电价 price[t] 和时间步长 dt 计算经济成本
- 这样设计的好处是：不同电价但相同功率和蒸发温度的时段，功率的状态转移成本是相同的，
  可以避免重复计算，提高效率

参考:
- VaryLoadVaryArea.jl 中 generateAndSolve 函数
- ConstantLoadConstantArea.jl 中 dpSolve 函数
- DPSolverCore.jl 中 ExhaustiveSolver 和 GoldenRatioSolver 函数
"""
function generateInitialSolution_DP(dp_params::DP_INITIAL_PARAMS,
    params::MILPModelParameters,
    sysVariables;
    cop_mode::Symbol = :piecewise,      # :piecewise 分档COP（默认，与MILP一致）| :continuous 连续COP函数（与DP精确解一致）
    designParameters = nothing,         # cop_mode=:continuous 时必须传入 designParameters
)

    n = params.n_segments
    m1 = params.m1
    m2 = params.m2
    m3 = params.m3
    mw = params.mw
    dt = dp_params.dt

    # 是否启用s6最小间隔约束（仅连续COP模式支持，分档模式s4/s6共用同一计算无法分离）
    use_s6 = (cop_mode == :continuous) && dp_params.y_s6 > 0.0

    # 连续COP模式：构造各时刻环境温度列表（由 params.Tair 提供，不足nt+1时补末值）
    Tair_list = collect(params.Tair)
    while length(Tair_list) < n + 1
        push!(Tair_list, Tair_list[end])
    end

    # 生成温度状态空间字典（多个温度系）
    T_state_space = get_temperature_state_space(dp_params.n, dp_params.q, params)

    @info "动态规划求解初始化"
    @info "温度系数量: $(length(T_state_space)), 时段数: $n"

    # 记录最优结果
    best_cost = 9999.0
    best_T_key = -1
    best_TsIndexList = nothing
    best_Ts_list = nothing
    best_stateMatrix = nothing
    best_P1sMatrix = nothing
    best_P2sMatrix = nothing
    best_P3sMatrix = nothing
    best_P1eMatrix = nothing
    best_P2eMatrix = nothing
    best_P3eMatrix = nothing
    best_PeslMatrix = nothing
    best_PessMatrix = nothing
    best_PeelMatrix = nothing
    best_PeesMatrix = nothing
    best_lambdaMatrix = nothing
    # 连续COP模式的s6专用最优矩阵（分档模式不使用）
    best_WaitIndexList = Int[]
    best_P2s_s6_Matrix = nothing
    best_Pesl_s6_Matrix = nothing
    best_P1e_s6_Matrix = nothing
    best_Peel_s6_Matrix = nothing
    best_lambda_s6_Matrix = nothing

    # 对每个温度系分别求解DP，选择成本最低的
    for (T_key, Ts_list) in T_state_space
        nT = length(Ts_list)
        nt = n

        @info "求解温度系: 起始温度=$(T_key), 状态数=$(nT)"

        if cop_mode == :continuous
            @assert designParameters !== nothing "cop_mode=:continuous 需要传入 designParameters"

            # 连续COP模式：复用精确解的张量构建（各时刻状态空间相同，首列重复）
            TsMatrix = repeat(reshape(Ts_list, :, 1), 1, nt + 1)
            Cn_power, Cs_power, stateMatrix,
            P1sMatrix, P2sMatrix, P3sMatrix,
            P1eMatrix, P2eMatrix, P3eMatrix,
            PeslMatrix, PessMatrix, PeelMatrix, PeesMatrix,
            lambdaMatrix, P2s_s6_Matrix, Pesl_s6_Matrix, P1e_s6_Matrix, Peel_s6_Matrix, lambda_s6_Matrix =
                buildCostTensor_Precise(
                    TsMatrix, nt, dt, designParameters, params, sysVariables, Tair_list;
                    y_s6 = dp_params.y_s6,
                    extra_forbidden = Set{Int}([7, 8]),   # 初始解禁用s7/s8
                )

            # 根据各时段电价将功率成本转换为经济成本
            if use_s6
                Cn_econ = fill(9999.0, nt, nT, nT)
                Cs_econ = fill(9999.0, nt, nT, nT)
                for t in 1:nt
                    price_t = sysVariables.price[t]
                    for i in 1:nT
                        for j in 1:nT
                            if Cn_power[t, i, j] < 9999.0
                                Cn_econ[t, i, j] = Cn_power[t, i, j] * price_t * dt
                            end
                            if Cs_power[t, i, j] < 9999.0
                                Cs_econ[t, i, j] = Cs_power[t, i, j] * price_t * dt
                            end
                        end
                    end
                end

                # 图DP求解（带s6间隔约束）
                cost, TsIndexList, WaitIndexList, best_start_idx, best_wait_idx =
                    DC.GraphLayerSolver(Cn_econ, Cs_econ, nT, nt, dp_params.y_s6, dt)
            else
                C_economic = fill(9999.0, nt, nT, nT)
                for t in 1:nt
                    price_t = sysVariables.price[t]
                    for i in 1:nT
                        for j in 1:nT
                            if Cn_power[t, i, j] < 9999.0
                                C_economic[t, i, j] = Cn_power[t, i, j] * price_t * dt
                            end
                        end
                    end
                end

                # 标准DP求解（无s6约束）
                if dp_params.solver_type == :Exhaustive
                    cost, TsIndexList, best_start_idx = DC.ExhaustiveSolver(C_economic, nT, nt)
                else
                    cost, TsIndexList, best_start_idx = DC.GoldenRatioSolver(C_economic, nT, nt)
                end
                WaitIndexList = Int[]
            end
        else
            # 分档COP模式：与MILP一致，取时段末温度对应档位COP
            C_power, stateMatrix,
            P1sMatrix, P2sMatrix, P3sMatrix,
            P1eMatrix, P2eMatrix, P3eMatrix,
            PeslMatrix, PessMatrix, PeelMatrix, PeesMatrix,
            lambdaMatrix = buildCostTensor(
                Ts_list, nt, dt, params, sysVariables
            )

            # 根据各时段电价将功率成本转换为经济成本
            C_economic = fill(9999.0, nt, nT, nT)
            for t in 1:nt
                price_t = sysVariables.price[t]
                for i in 1:nT
                    for j in 1:nT
                        if C_power[t, i, j] < 9999.0
                            C_economic[t, i, j] = C_power[t, i, j] * price_t * dt
                        end
                    end
                end
            end

            # 求解DP
            cost, TsIndexList, best_start_idx = if dp_params.solver_type == :Exhaustive
                DC.ExhaustiveSolver(C_economic, nT, nt)
            else
                DC.GoldenRatioSolver(C_economic, nT, nt)
            end
            WaitIndexList = Int[]
        end

        @info "温度系 起始温度=$(T_key) 求解完成，成本: $(cost)"

        # 记录最优结果
        if cost < best_cost
            best_cost = cost
            best_T_key = T_key
            best_TsIndexList = TsIndexList
            best_Ts_list = Ts_list
            if cop_mode == :continuous
                best_WaitIndexList = copy(WaitIndexList)
                best_P2s_s6_Matrix = copy(P2s_s6_Matrix)
                best_Pesl_s6_Matrix = copy(Pesl_s6_Matrix)
                best_P1e_s6_Matrix = copy(P1e_s6_Matrix)
                best_Peel_s6_Matrix = copy(Peel_s6_Matrix)
                best_lambda_s6_Matrix = copy(lambda_s6_Matrix)
            end
            best_stateMatrix = stateMatrix
            best_P1sMatrix = P1sMatrix
            best_P2sMatrix = P2sMatrix
            best_P3sMatrix = P3sMatrix
            best_P1eMatrix = P1eMatrix
            best_P2eMatrix = P2eMatrix
            best_P3eMatrix = P3eMatrix
            best_PeslMatrix = PeslMatrix
            best_PessMatrix = PessMatrix
            best_PeelMatrix = PeelMatrix
            best_PeesMatrix = PeesMatrix
            best_lambdaMatrix = lambdaMatrix
        end
    end

    @info "动态规划求解完成，最优温度系起始温度: $(best_T_key), 最优成本: $(best_cost)"

    # 根据最优温度系提取结果
    Ts_list = best_Ts_list
    TsIndexList = best_TsIndexList
    stateMatrix = best_stateMatrix
    P1sMatrix = best_P1sMatrix
    P2sMatrix = best_P2sMatrix
    P3sMatrix = best_P3sMatrix
    P1eMatrix = best_P1eMatrix
    P2eMatrix = best_P2eMatrix
    P3eMatrix = best_P3eMatrix
    PeslMatrix = best_PeslMatrix
    PessMatrix = best_PessMatrix
    PeelMatrix = best_PeelMatrix
    PeesMatrix = best_PeesMatrix
    lambdaMatrix = best_lambdaMatrix

    Ts_result = Ts_list[TsIndexList]

    if best_cost >= 9999.0
        @warn "动态规划未找到可行解，使用默认初始温度"
        Ts_result = fill(params.Tsmin, n+1)
        TsIndexList = fill(1, n+1)
    end

    states = zeros(Int, n)
    P1_list = zeros(n)
    P2_list = zeros(n)
    P3_list = zeros(n)
    Pe_l_list = zeros(n)
    Pe_s_list = zeros(n)
    P1s_list = zeros(n)
    P2s_list = zeros(n)
    P3s_list = zeros(n)
    P1e_list = zeros(n)
    P2e_list = zeros(n)
    P3e_list = zeros(n)
    lambda_list = zeros(n)

    for i in 1:n
        idx_start = TsIndexList[i]
        idx_end = TsIndexList[i+1]

        # s6判断：连续COP模式启用y_s6时，wait被重置为0表示该时段使用了s6
        s6_used = use_s6 && (best_WaitIndexList[i+1] == 0)

        if s6_used
            # s6状态：从s6专用矩阵提取（状态6没有P1s, P3s, Pess, Pees, P2e, P3e）
            states[i] = 6
            P1s_list[i] = 0.0
            P2s_list[i] = best_P2s_s6_Matrix[i, idx_start, idx_end]
            P3s_list[i] = 0.0
            P1e_list[i] = best_P1e_s6_Matrix[i, idx_start, idx_end]
            P2e_list[i] = 0.0
            P3e_list[i] = 0.0
            Pe_l_list[i] = best_Pesl_s6_Matrix[i, idx_start, idx_end] + best_Peel_s6_Matrix[i, idx_start, idx_end]
            Pe_s_list[i] = 0.0
            lambda_list[i] = best_lambda_s6_Matrix[i, idx_start, idx_end]
        else
            states[i] = stateMatrix[i, idx_start, idx_end]
            P1s_list[i] = P1sMatrix[i, idx_start, idx_end]
            P2s_list[i] = P2sMatrix[i, idx_start, idx_end]
            P3s_list[i] = P3sMatrix[i, idx_start, idx_end]
            P1e_list[i] = P1eMatrix[i, idx_start, idx_end]
            P2e_list[i] = P2eMatrix[i, idx_start, idx_end]
            P3e_list[i] = P3eMatrix[i, idx_start, idx_end]
            Pe_l_list[i] = PeslMatrix[i, idx_start, idx_end] + PeelMatrix[i, idx_start, idx_end]
            Pe_s_list[i] = PessMatrix[i, idx_start, idx_end] + PeesMatrix[i, idx_start, idx_end]
            lambda_list[i] = lambdaMatrix[i, idx_start, idx_end]
        end
        P1_list[i] = P1s_list[i] + P1e_list[i]
        P2_list[i] = P2s_list[i] + P2e_list[i]
        P3_list[i] = P3s_list[i] + P3e_list[i]
    end

    # 从 best_cost 中扣除功率正则项（仅连续COP模式带epsilon平滑，分档模式无正则项）
    if cop_mode == :continuous
        regularization_sum = 0.0
        epsilon = sysVariables.epsilon
        for t in 1:n
            i = TsIndexList[t]
            j = TsIndexList[t+1]
            price_t = sysVariables.price[t]
            s6_used = use_s6 && (best_WaitIndexList[t+1] == 0)
            if s6_used
                P2s_s6 = best_P2s_s6_Matrix[t, i, j]
                P1e_s6 = best_P1e_s6_Matrix[t, i, j]
                Pesl_s6 = best_Pesl_s6_Matrix[t, i, j]
                Peel_s6 = best_Peel_s6_Matrix[t, i, j]
                reg_t = epsilon * ((P2s_s6 + P1e_s6)^2 + (Pesl_s6 + Peel_s6)^2)
            else
                P1s = best_P1sMatrix[t, i, j]
                P1e = best_P1eMatrix[t, i, j]
                P2s = best_P2sMatrix[t, i, j]
                P2e = best_P2eMatrix[t, i, j]
                P3s = best_P3sMatrix[t, i, j]
                P3e = best_P3eMatrix[t, i, j]
                Pesl = best_PeslMatrix[t, i, j]
                Pess = best_PessMatrix[t, i, j]
                Peel = best_PeelMatrix[t, i, j]
                Pees = best_PeesMatrix[t, i, j]
                reg_t = epsilon * ((P1s + P1e)^2 + (P2s + P2e)^2 + (P3s + P3e)^2 + (Pesl + Pess + Peel + Pees)^2)
            end
            regularization_sum += reg_t * price_t * dt
        end
        best_cost -= regularization_sum
    end

    s = zeros(8, n)
    for i in 1:n
        state = states[i]
        if 1 <= state <= 8
            s[state, i] = 1.0
        else
            s[1, i] = 1.0
        end
    end

    delta_su = [Ts_result[i] >= params.Tuse + params.dTs ? 1.0 : 0.0 for i in 1:n]
    delta_tilde = zeros(n)

    z1 = zeros(n, m1)
    z2 = zeros(n, m2)
    z3 = zeros(n, m3)
    zw = zeros(n, mw)

    for i in 1:n
        Tg1 = params.T1g[i, :]
        Tg2 = params.T2g[i, :]
        Tg3 = params.T3g[i, :]
        Tgw = params.Twg[i, :]

        seg1 = findCOPSegment(Ts_result[i], Tg1)
        seg2 = findCOPSegment(Ts_result[i], Tg2)
        seg3 = findCOPSegment(Ts_result[i], Tg3)
        segw = findCOPSegment(Ts_result[i], Tgw)

        z1[i, seg1] = 1.0
        z2[i, seg2] = 1.0
        z3[i, seg3] = 1.0
        zw[i, segw] = 1.0
    end

    y1 = zeros(n)
    w1 = zeros(n, m1)
    w3 = zeros(n, m3)
    for i in 1:n
        w3[i, :] .= z3[i, :]
    end

    P_k_s = zeros(8, n, 3)
    P_k_e = zeros(8, n, 3)

    for i in 1:n
        state = states[i]
        COP1_i = sum(params.COP1v[j, i] * z1[i, j] for j ∈ 1:m1)

        if state == 1
            P_k_s[1, i, 1] = P1s_list[i]
            P_k_e[1, i, 1] = P1e_list[i]
        elseif state == 2
            P_k_s[2, i, 2] = P2s_list[i]
            P_k_e[2, i, 2] = P2e_list[i]
        elseif state == 3
            P_k_s[3, i, 3] = P3s_list[i]
            P_k_e[3, i, 3] = P3e_list[i]
        elseif state == 4
            P_k_s[4, i, 1] = P1s_list[i]
            P_k_s[4, i, 2] = P2s_list[i]
            P_k_e[4, i, 1] = P1e_list[i]
            P_k_e[4, i, 2] = P2e_list[i]
        elseif state == 5
            P_k_s[5, i, 1] = P1s_list[i]
            P_k_s[5, i, 3] = P3s_list[i]
            P_k_e[5, i, 1] = P1e_list[i]
            P_k_e[5, i, 3] = P3e_list[i]
        elseif state == 6
            P_k_s[6, i, 2] = P2s_list[i]
            P_k_e[6, i, 1] = P1e_list[i]
        elseif state == 7
            P_k_s[7, i, 3] = P3s_list[i]
            P_k_e[7, i, 1] = P1e_list[i]
        elseif state == 8
            P_k_s[8, i, 1] = P1s_list[i]
            P_k_s[8, i, 3] = P3s_list[i]
            P_k_e[8, i, 1] = P1e_list[i]
        end
    end

    P_el = Pe_l_list
    P_es = Pe_s_list

    u1 = zeros(8, n)
    u2 = zeros(8, n)
    u3 = zeros(n)
    for k in 1:8, i in 1:n
        u1[k, i] = s[k, i] * (1 - y1[i])
        if k in [1, 2, 3, 4, 5, 8]
            u2[k, i] = u1[k, i] * P_k_s[k, i, 1]
        else
            u2[k, i] = u1[k, i] * P_k_e[k, i, 1]
        end
    end

    v1 = zeros(8, n, m1)
    v2 = zeros(8, n, m1)
    v3 = zeros(8, n, m2)
    v4 = zeros(8, n, m2)
    v5 = zeros(8, n)
    v6 = zeros(8, n, m3)
    v7 = zeros(8, n, m1)
    v8 = zeros(8, n)
    v9 = zeros(8, n, m2)

    for k in 1:8, i in 1:n
        v5[k, i] = s[k, i] * P_k_s[k, i, 3]
        v8[k, i] = s[k, i] * P_k_s[k, i, 2]

        for j in 1:m1
            v1[k, i, j] = w1[i, j] * s[k, i]
            v2[k, i, j] = v1[k, i, j] * P_k_s[k, i, 1]
            v7[k, i, j] = v5[k, i] * w1[i, j]
        end

        for j in 1:m2
            v3[k, i, j] = z2[i, j] * s[k, i]
            v4[k, i, j] = v3[k, i, j] * P_k_s[k, i, 2]
            v9[k, i, j] = v8[k, i] * z2[i, j]
        end

        for j in 1:m3
            v6[k, i, j] = v5[k, i] * w3[i, j]
        end
    end

    COP1e = zeros(n)
    COP2e = zeros(n)
    COP3e = zeros(n)
    COPwe = zeros(n)
    COP1 = zeros(n)
    COP2 = zeros(n)
    COP3 = zeros(n)
    COPw = zeros(n)

    for i in 1:n
        COP1e[i] = sum(params.COP1v[j, i] * z1[i, j] for j ∈ 1:m1)
        COP2e[i] = sum(params.COP2v[j, i] * z2[i, j] for j ∈ 1:m2)
        COP3e[i] = sum(params.COP3v[j, i] * z3[i, j] for j ∈ 1:m3)
        COPwe[i] = sum(params.COPwv[j, i] * zw[i, j] for j ∈ 1:mw)

        COP1[i] = COP1e[i]
        COP2[i] = COP2e[i]
        COP3[i] = COP3e[i]
        COPw[i] = COPwe[i]
    end

    P_total = zeros(4, n)
    for i in 1:n
        P_total[1, i] = sum(P_k_s[k, i, 1] + P_k_e[k, i, 1] for k=1:8)
        P_total[2, i] = sum(P_k_s[k, i, 2] + P_k_e[k, i, 2] for k=1:8)
        P_total[3, i] = sum(P_k_s[k, i, 3] + P_k_e[k, i, 3] for k=1:8)
        P_total[4, i] = P_el[i] + P_es[i]
    end

    dP_pos = zeros(4, n)
    dP_neg = zeros(4, n)
    for j in 1:4
        for i in 1:n
            i_next = mod1(i+1, n)
            diff = P_total[j, i_next] - P_total[j, i]
            dP_pos[j, i] = max(diff, 0)
            dP_neg[j, i] = max(-diff, 0)
        end
    end

    initial_solution = InitialSolution(
        Ts_result,
        s, delta_su, delta_tilde,
        z1, z2, z3, zw,
        w1, w3,
        y1,
        P_k_s, P_k_e,
        P_el, P_es,
        u1, u2, u3,
        v1, v2, v3, v4, v5, v6, v7, v8, v9,
        COP1e, COP2e, COP3e, COPwe,
        COP1, COP2, COP3, COPw,
        params.C_heatpump,
        params.C_boiler,
        dP_pos,
        dP_neg,
    )

    return initial_solution, best_cost
end

export DP_INITIAL_PARAMS, generateInitialSolution_DP, get_temperature_state_space