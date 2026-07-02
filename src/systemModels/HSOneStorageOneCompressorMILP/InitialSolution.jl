"""
    InitialSolution

MILP模型的初值数据结构，用于为求解器提供热启动点。
"""
struct InitialSolution
    # 温度变量
    Ts::Vector{Float64}            # 蓄热温度 [n+1]
    
    # 状态变量
    s::Matrix{Float64}             # 运行状态 [8, n]
    delta_su::Vector{Float64}      # 温度标志 [n]
    delta_tilde::Vector{Float64}   # 过渡标志 [n]
    
    # COP分档变量
    z1::Matrix{Float64}            # COP1分档 [n, m1]
    z2::Matrix{Float64}            # COP2分档 [n, m2]
    z3::Matrix{Float64}            # COP3分档 [n, m3]
    zw::Matrix{Float64}            # COPw分档 [n, mw]
    
    # 功率变量
    P_k_s::Array{Float64, 3}       # 起始功率 [8, n, 3]
    P_k_e::Array{Float64, 3}       # 结束功率 [8, n, 3]
    
    # 电加热功率
    P_el::Vector{Float64}          # 电加热补热功率 [n]
    P_es::Vector{Float64}          # 电加热蓄热功率 [n]
    
    # 容量变量初值（当 optimizeCapacity=true 时使用）
    C_heatpump::Float64            # 热泵容量初值 kW
    C_boiler::Float64              # 电锅炉容量初值 kW
end

"""
    findCOPSegment(T::Float64, Tg::Vector{Float64}) -> Int

根据温度 T 找到对应的 COP 档位。

# 参数
- `T`: 蓄热温度 (℃)
- `Tg`: 温度分档界限 [m-1]，例如 [130, 140, 150, 160, 170, 180]

# 返回
- 档位索引 j (1 到 m)

# 分档规则
- 档位1: T <= Tg[1]
- 档位2: Tg[1] < T <= Tg[2]
- ...
- 档位m: T > Tg[m-1]
"""
function findCOPSegment(T::Float64, Tg::Vector{Float64})
    m = length(Tg) + 1
    for j in 1:m-1
        if T <= Tg[j]
            return j
        end
    end
    return m  # 最高档位
end

"""
    generateInitialSolution_HeatPumpOnly(params::MILPModelParameters)

生成"全程热泵供热"的初值方案：
- 蓄热温度恒为 Tuse
- 运行状态恒为 s1（全程热泵供热）

# 参数
- `params`: MILP模型参数

# 返回
- `InitialSolution`: 初值数据结构

# 初值设置说明
1. **蓄热温度**: Ts[i] = Tuse（用热温度），所有时段温度相同
2. **运行状态**: s[1, i] = 1，其他状态 = 0（全程热泵供热）
3. **温度标志**: delta_su = 0，delta_tilde = 0（温度低于 Tu+ΔTs）
4. **COP分档**: 根据 Tuse 选择对应的温度档位
5. **功率变量**: 根据热负荷计算热泵功率
6. **电加热功率**: 为 0（热泵满足全部负荷）
"""
function generateInitialSolution_HeatPumpOnly(params::MILPModelParameters)
    n = params.n_segments
    m1, m2, m3, mw = params.m1, params.m2, params.m3, params.mw
    
    # 1. 蓄热温度：恒为 Tuse
    Ts = fill(params.Tuse, n + 1)
    
    # 2. 运行状态：s1=1，其他=0
    s = zeros(8, n)
    s[1, :] .= 1.0
    
    # 3. 温度标志：温度低于 Tu+ΔTs，所以 delta_su=0
    delta_su = zeros(n)
    delta_tilde = zeros(n)
    
    # 4. COP分档：根据 Ts 的温度选择对应档位
    z1 = zeros(n, m1)
    z2 = zeros(n, m2)
    z3 = zeros(n, m3)
    zw = zeros(n, mw)
    
    # 获取温度分档界限（需要从参数中提取）
    # 注意：T1g 的维度是 [n, m1-1]，这里取第一个时段的分档界限
    for i in 1:n
        # 找到 Ts[i] 所在的档位
        # T1g 的维度是 [n, m1+1]，其中第1列和第m1+1列是端点
        # 分界点是 T1g[i, 2:m1]，共 m1-1 个
        Tg1 = params.T1g[i, 2:m1]  # COP1温度分档分界点 [m1-1]
        Tg2 = params.T2g[i, 2:m2]
        Tg3 = params.T3g[i, 2:m3]
        Tgw = params.Twg[i, 2:mw]
        
        seg1 = findCOPSegment(params.Tuse, Tg1)
        seg2 = findCOPSegment(params.Tuse, Tg2)
        seg3 = findCOPSegment(params.Tuse, Tg3)
        segw = findCOPSegment(params.Tuse, Tgw)
        
        #=
        z1[i, seg1] = 1.0
        z2[i, seg2] = 1.0
        z3[i, seg3] = 1.0
        zw[i, segw] = 1.0
        =#
        z1[i, 1] = 1.0
        z2[i, 1] = 1.0
        z3[i, 1] = 1.0
        zw[i, 1] = 1.0
    end
    
    # 5. 功率变量：根据热负荷计算
    P_k_s = zeros(8, n, 3)
    P_k_e = zeros(8, n, 3)

    for i in 1:n
        # 状态1：热泵直接供热，功率 = 热负荷 / COP
        COP1_i = sum(params.COP1v[j, i] * z1[i, j] for j=1:m1)
        #println("调试: 第$(i)时段 COP1_i=$COP1_i, heatLoad=$(params.heatLoad[i])")
        
        P_k_s[1, i, 1] = params.heatLoad[i] / COP1_i
        
        #=
        if params.heatLoad[i]!=0
            println("调试: 第$(i)时段 heatLoad=$(params.heatLoad[1]), COP1_i=$COP1_i, P_k_s[1,1,1]=$(P_k_s[1,i,1])")
            #println("调试: COP1v[:,1] = $(params.COP1v[:,1])")
            #println("调试: z1[1,:] = $(z1[1,:])")
        end
        =#
    end
    
    # 6. 电加热功率：为0（热泵满足全部负荷）
    P_el = zeros(n)
    P_es = zeros(n)
    
    # 7. 容量变量初值：当 optimizeCapacity=true 时使用
    # 默认初值 = 1.0 * maximum(heatLoad)，确保能满足最大热负荷需求
    max_heat_load = maximum(params.heatLoad)
    C_heatpump_init = 1.0 * max_heat_load
    C_boiler_init = 0.0
    
    return InitialSolution(
        Ts, s, delta_su, delta_tilde,
        z1, z2, z3, zw,
        P_k_s, P_k_e,
        P_el, P_es,
        C_heatpump_init, C_boiler_init
    )
end

"""
    setInitialSolution(model::Model, initial::InitialSolution)

将初值设置到 JuMP 模型中，为求解器提供热启动点。

# 参数
- `model`: JuMP 模型
- `initial`: 初值数据结构

# 注意
此函数使用 JuMP 的 `set_start_value` 方法设置变量初值。
求解器会使用这些初值作为搜索起点，可能加速求解过程。
"""
function setInitialSolution(model::Model, initial::InitialSolution)
    # 设置温度初值
    if haskey(model, :Ts)
        set_start_value.(model[:Ts], initial.Ts)
    end
    
    # 设置状态初值
    if haskey(model, :s)
        set_start_value.(model[:s], initial.s)
    end
    
    # 设置温度标志初值
    if haskey(model, :delta_su)
        set_start_value.(model[:delta_su], initial.delta_su)
    end
    if haskey(model, :delta_tilde)
        set_start_value.(model[:delta_tilde], initial.delta_tilde)
    end
    
    # 设置COP分档初值
    if haskey(model, :z1)
        set_start_value.(model[:z1], initial.z1)
    end
    if haskey(model, :z2)
        set_start_value.(model[:z2], initial.z2)
    end
    if haskey(model, :z3)
        set_start_value.(model[:z3], initial.z3)
    end
    if haskey(model, :zw)
        set_start_value.(model[:zw], initial.zw)
    end
    
    # 设置功率初值
    if haskey(model, :P_k_s)
        set_start_value.(model[:P_k_s], initial.P_k_s)
    end
    if haskey(model, :P_k_e)
        set_start_value.(model[:P_k_e], initial.P_k_e)
    end
    
    # 设置电加热初值
    if haskey(model, :P_el)
        set_start_value.(model[:P_el], initial.P_el)
    end
    if haskey(model, :P_es)
        set_start_value.(model[:P_es], initial.P_es)
    end
    
    # 设置容量变量初值（当 optimizeCapacity=true 时使用）
    if haskey(model, :C_heatpump)
        set_start_value(model[:C_heatpump], initial.C_heatpump)
    end
    if haskey(model, :C_boiler)
        set_start_value(model[:C_boiler], initial.C_boiler)
    end
end

export findCOPSegment, generateInitialSolution_HeatPumpOnly
setInitialSolution