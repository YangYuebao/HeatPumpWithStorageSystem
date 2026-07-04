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
    
    # COP辅助变量
    w1::Matrix{Float64}            # w1 = y1 * z1 [n, m1]
    w3::Matrix{Float64}            # w3 = (1-y1) * z3 [n, m3]
    
    # 状态组合变量
    y1::Vector{Float64}            # y1 = s5 || s8 [n]
    
    # 功率变量
    P_k_s::Array{Float64, 3}       # 起始功率 [8, n, 3]
    P_k_e::Array{Float64, 3}       # 结束功率 [8, n, 3]
    
    # 电加热功率
    P_el::Vector{Float64}          # 电加热补热功率 [n]
    P_es::Vector{Float64}          # 电加热蓄热功率 [n]
    
    # 热负荷辅助变量
    u1::Matrix{Float64}            # u_k,i^1 = s_k,i * (1-y1,i) [8, n]
    u2::Matrix{Float64}            # u_k,i^2 = u_k,i^1 * P_k,i,1^s [8, n]
    u3::Vector{Float64}            # u_8,i^3 = s_8,i * P_8,i,1^e [n]
    
    # 功率乘积辅助变量
    v1::Array{Float64, 3}          # v_k,i,j^1 = w1,i,j * s_k,i [8, n, m1]
    v2::Array{Float64, 3}          # v_k,i,j^2 = v_k,i,j^1 * P_k,i,1^s [8, n, m1]
    v3::Array{Float64, 3}          # v_k,i,j^3 = z2,i,j * s_k,i [8, n, m2]
    v4::Array{Float64, 3}          # v_k,i,j^4 = v_k,i,j^3 * P_k,i,2^s [8, n, m2]
    v5::Matrix{Float64}            # v_k,i^5 = s_k,i * P_k,i,3^s [8, n]
    v6::Array{Float64, 3}          # v_k,i,j^6 = v_k,i^5 * w3,i,j [8, n, m3]
    v7::Array{Float64, 3}          # v_k,i,j^7 = v_k,i^5 * w1,i,j [8, n, m1]
    v8::Matrix{Float64}            # v_k,i^8 = s_k,i * P_k,i,2^s [8, n]
    v9::Array{Float64, 3}          # v_k,i,j^9 = v_k,i^8 * z2,i,j [8, n, m2]
    
    # COP估计值和实际值
    COP1e::Vector{Float64}         # COP1估计值 [n]
    COP2e::Vector{Float64}         # COP2估计值 [n]
    COP3e::Vector{Float64}         # COP3估计值 [n]
    COPwe::Vector{Float64}         # COPw估计值 [n]
    COP1::Vector{Float64}          # COP1实际值 [n]
    COP2::Vector{Float64}          # COP2实际值 [n]
    COP3::Vector{Float64}          # COP3实际值 [n]
    COPw::Vector{Float64}          # COPw实际值 [n]
    
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
    
    for i in 1:n
        Tg1 = params.T1g[i, 2:m1]
        Tg2 = params.T2g[i, 2:m2]
        Tg3 = params.T3g[i, 2:m3]
        Tgw = params.Twg[i, 2:mw]
        
        seg1 = findCOPSegment(params.Tuse, Tg1)
        seg2 = findCOPSegment(params.Tuse, Tg2)
        seg3 = findCOPSegment(params.Tuse, Tg3)
        segw = findCOPSegment(params.Tuse, Tgw)
        
        z1[i, seg1] = 1.0
        z2[i, seg2] = 1.0
        z3[i, seg3] = 1.0
        zw[i, segw] = 1.0
    end
    
    # 5. COP辅助变量
    # w1 = y1 * z1, w3 = (1-y1) * z3
    y1 = zeros(n)
    w1 = zeros(n, m1)
    w3 = zeros(n, m3)
    for i in 1:n
        w3[i, :] .= z3[i, :]
    end
    
    # 6. 功率变量：根据热负荷计算
    P_k_s = zeros(8, n, 3)
    P_k_e = zeros(8, n, 3)
    
    for i in 1:n
        COP1_i = sum(params.COP1v[j, i] * z1[i, j] for j=1:m1)
        P_k_s[1, i, 1] = min(params.heatLoad[i], params.C_heatpump) / COP1_i
    end
    
    # 7. 电加热功率
    P_el = max.(params.heatLoad .- params.C_heatpump, 0)
    P_es = zeros(n)
    
    # 8. 热负荷辅助变量
    u1 = zeros(8, n)
    u2 = zeros(8, n)
    u3 = zeros(n)
    for k in 1:8, i in 1:n
        u1[k, i] = s[k, i] * (1 - y1[i])
        if k in [1,2,3,4,5,8]
            u2[k, i] = u1[k, i] * P_k_s[k, i, 1]
        else
            u2[k, i] = u1[k, i] * P_k_e[k, i, 1]
        end
    end
    
    # 9. 功率乘积辅助变量
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
        for j in 1:m1
            v1[k, i, j] = w1[i, j] * s[k, i]
            if k in [1,2,3,4,5,8]
                v2[k, i, j] = v1[k, i, j] * P_k_s[k, i, 1]
            else
                v2[k, i, j] = v1[k, i, j] * P_k_e[k, i, 1]
            end
        end
        for j in 1:m2
            v3[k, i, j] = z2[i, j] * s[k, i]
            v4[k, i, j] = v3[k, i, j] * P_k_s[k, i, 2]
        end
        v5[k, i] = s[k, i] * P_k_s[k, i, 3]
        for j in 1:m3
            v6[k, i, j] = v5[k, i] * w3[i, j]
        end
        for j in 1:m1
            v7[k, i, j] = v5[k, i] * w1[i, j]
        end
        v8[k, i] = s[k, i] * P_k_s[k, i, 2]
        for j in 1:m2
            v9[k, i, j] = v8[k, i] * z2[i, j]
        end
    end
    
    # 10. COP估计值和实际值
    COP1e = zeros(n)
    COP2e = zeros(n)
    COP3e = zeros(n)
    COPwe = zeros(n)
    COP1 = zeros(n)
    COP2 = zeros(n)
    COP3 = zeros(n)
    COPw = zeros(n)
    
    for i in 1:n
        COP2e[i] = sum(params.COP2v[j, i] * z2[i, j] for j=1:m2)
        COPwe[i] = sum(params.COPwv[j, i] * zw[i, j] for j=1:mw)
        COP2[i] = COP2e[i]
        COPw[i] = COPwe[i]
        COP1[i] = (1 - y1[i]) * params.COPca[i] + sum(params.COP1v[j, i] * w1[i, j] for j=1:m1)
        COP3[i] = sum(params.COP3v[j, i] * w3[i, j] for j=1:m3) + sum(params.COP1v[j, i] * w1[i, j] for j=1:m1)
    end
    
    # 11. 容量变量初值
    C_heatpump_init = params.C_heatpump
    C_boiler_init = params.C_boiler
    
    return InitialSolution(
        Ts, s, delta_su, delta_tilde,
        z1, z2, z3, zw,
        w1, w3, y1,
        P_k_s, P_k_e,
        P_el, P_es,
        u1, u2, u3,
        v1, v2, v3, v4, v5, v6, v7, v8, v9,
        COP1e, COP2e, COP3e, COPwe,
        COP1, COP2, COP3, COPw,
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
function setInitialSolution(model::Model, initial::InitialSolution, params::MILPModelParameters)
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
    
    # 设置COP辅助变量初值
    if haskey(model, :w1)
        set_start_value.(model[:w1], initial.w1)
    end
    if haskey(model, :w3)
        set_start_value.(model[:w3], initial.w3)
    end
    
    # 设置状态组合变量初值
    if haskey(model, :y1)
        set_start_value.(model[:y1], initial.y1)
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
    
    # 设置热负荷辅助变量初值
    if haskey(model, :u1)
        set_start_value.(model[:u1], initial.u1)
    end
    if haskey(model, :u2)
        set_start_value.(model[:u2], initial.u2)
    end
    if haskey(model, :u3)
        set_start_value.(model[:u3], initial.u3)
    end
    
    # 设置功率乘积辅助变量初值
    if haskey(model, :v1)
        set_start_value.(model[:v1], initial.v1)
    end
    if haskey(model, :v2)
        set_start_value.(model[:v2], initial.v2)
    end
    if haskey(model, :v3)
        set_start_value.(model[:v3], initial.v3)
    end
    if haskey(model, :v4)
        set_start_value.(model[:v4], initial.v4)
    end
    if haskey(model, :v5)
        set_start_value.(model[:v5], initial.v5)
    end
    if haskey(model, :v6)
        set_start_value.(model[:v6], initial.v6)
    end
    if haskey(model, :v7)
        set_start_value.(model[:v7], initial.v7)
    end
    if haskey(model, :v8)
        set_start_value.(model[:v8], initial.v8)
    end
    if haskey(model, :v9)
        set_start_value.(model[:v9], initial.v9)
    end
    
    # 设置COP估计值和实际值初值
    if haskey(model, :COP1e)
        set_start_value.(model[:COP1e], initial.COP1e)
    end
    if haskey(model, :COP2e)
        set_start_value.(model[:COP2e], initial.COP2e)
    end
    if haskey(model, :COP3e)
        set_start_value.(model[:COP3e], initial.COP3e)
    end
    if haskey(model, :COPwe)
        set_start_value.(model[:COPwe], initial.COPwe)
    end
    if haskey(model, :COP1)
        set_start_value.(model[:COP1], initial.COP1)
    end
    if haskey(model, :COP2)
        set_start_value.(model[:COP2], initial.COP2)
    end
    if haskey(model, :COP3)
        set_start_value.(model[:COP3], initial.COP3)
    end
    if haskey(model, :COPw)
        set_start_value.(model[:COPw], initial.COPw)
    end
    
    # 设置容量变量初值（当 optimizeCapacity=true 时使用）
    if haskey(model, :C_heatpump) && params.optimizeCapacity
        set_start_value(model[:C_heatpump], initial.C_heatpump)
    end
    if haskey(model, :C_boiler) && params.optimizeCapacity
        set_start_value(model[:C_boiler], initial.C_boiler)
    end
end

export findCOPSegment, generateInitialSolution_HeatPumpOnly
setInitialSolution