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
        Tg1 = params.T1g[i, :]
        Tg2 = params.T2g[i, :]
        Tg3 = params.T3g[i, :]
        Tgw = params.Twg[i, :]
        
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
    
    # 12. 功率变差变量初值（dP_pos, dP_neg）
    # 根据实际功率变化计算
    # j=1: 热泵直接供热功率
    # j=2: 蓄热供热功率
    # j=3: 热泵向蓄热储热功率
    # j=4: 电加热功率
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
        C_heatpump_init, C_boiler_init,
        dP_pos, dP_neg
    )
end

export generateInitialSolution_HeatPumpOnly