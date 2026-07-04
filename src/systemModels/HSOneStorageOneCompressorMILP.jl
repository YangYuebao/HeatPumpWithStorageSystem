"""
    PressedWaterOneStorageOneCompressor_MILP 系统模型
    基于混合整数线性规划(SOS1/二进制)的热泵承压水蓄热系统运行调度优化

    参考文档: algorithms/热泵运行调度算法_混合整数线性规划_SOS1.md

    支持求解器:
    - HiGHS (开源求解器)
    - COPT (商业求解器,需要安装COPT.jl)
"""

# MOI 是 JuMP 重新导出的 MathOptInterface 接口
# 用于定义 SOS1 等特殊约束
const MOI = JuMP.MOI

# COPT求解器支持 (可选加载)
# 如果未安装COPT.jl，使用HiGHS作为备选求解器

# =============================================================
# 1. 参数结构体
# =============================================================

"""
    MILPModelParameters
MILP模型参数，包含系统设计参数、外部输入曲线和离散化COP数据

文档引用: 第1.2.1节 预设参数与常量
"""
@kwdef mutable struct MILPModelParameters
    # === 系统设计参数 (文档1.2.1节) ===
    Tuse::Float64              # 用热温度 T_u ℃
    Tsmin::Float64             # 蓄热最小温度 T_s,min ℃
    Tsmax::Float64             # 蓄热最大温度 T_s,max ℃
    Tcmax::Float64             # 压缩机最高排气饱和温度 T_c,max ℃
    dTs::Float64               # 传热温差 ΔT_s ℃
    Te::Float64                # 低温热泵吸气饱和温度 T_e ℃
    cpm::Float64               # 蓄热罐热容 c_p·m kJ/K
    dt::Float64                # 时间步长 Δt_i h
    COPmax::Float64            # 最大制热系数 COP_max
    COPdesign::Vector{Float64}         # 设计工况COP COP_design
    lambdaNorm::Float64        # 正则化参数 λ_norm
    smoother::Float64          # 平滑系数
    
    # === 设备容量参数 (文档3.2节) ===
    C_heatpump::Float64        # 热泵容量 kW (决策变量时为0，表示需优化)
    C_boiler::Float64          # 电锅炉容量 kW (决策变量时为0，表示需优化)
    C_storage::Float64         # 蓄热系统容量 kWh (常数，与cpm成正比)
    optimizeCapacity::Bool     # 是否优化设备容量

    # === 设备成本系数 (文档1.2.1节) ===
    c_heatpump::Float64        # 单位热泵投资成本 元/kW
    c_storage::Float64         # 单位蓄热系统成本 元/kWh
    c_boiler::Float64         # 单位电锅炉成本 元/kW

    # === 经济性参数 (文档1.2.1节) ===
    D::Float64                 # 年运行天数 天/年
    r::Float64                 # 折现率
    N::Float64                 # 运行年限 年
    
    # === 外部输入曲线 ===
    hourlyTariff::Vector{Float64}    # 电价向量 c_g,i 元/kWh
    heatLoad::Vector{Float64}        # 热负荷向量 L_i kW
    Tair::Vector{Float64}            # 环境温度向量 ℃
    dt_list::Vector{Float64}         # 各时段时长 Δt_i h
    
    # === 时间分段参数 ===
    n_segments::Int                  # 合并后的时段数 (length(dt_list))
    segmentTariff::Vector{Float64}   # 合并后的时段电价 [n_segments]
    
    # === COP分档数据 (文档1.2.1节) ===
    m1::Int                                   # COP1分档数
    COP1v::Matrix{Float64}                    # COP1分档值 [m1, n_segments]
    T1g::Matrix{Float64}                      # COP1温度分档界限 [n_segments, m1]
    
    m2::Int                                   # COP2分档数
    COP2v::Matrix{Float64}                    # COP2分档值 [m2, n_segments]
    T2g::Matrix{Float64}                      # COP2温度分档界限 [n_segments, m2]
    
    m3::Int                                   # COP3分档数
    COP3v::Matrix{Float64}                    # COP3分档值 [m3, n_segments]
    T3g::Matrix{Float64}                      # COP3温度分档界限 [n_segments, m3]
    
    mw::Int                                   # COPw分档数
    COPwv::Matrix{Float64}                    # COPw分档值 [mw, n_segments]
    Twg::Matrix{Float64}                      # COPw温度分档界限 [n_segments, mw]
    
    # === 定参数COP ===
    COPca::Vector{Float64}                    # COP_ca(T_e, T_u)

    # === 大M参数 (文档1.2.1节) ===
    M1::Float64               # 蓄热温度范围 M1 = Tsmax - Tsmin ℃
    M2::Float64               # 温度约束大M M2 = max(Tcmax-dTs-Tsmin, Tsmax-Tcmax+dTs) ℃
    M3::Float64               # 温度约束大M M3 = max(Tuse+dTs-Tsmin, Tsmax-Tuse-dTs) ℃
    M4::Float64               # 档位温差最大值 ℃
    M5::Float64               # 功率最大值 kW

    # === 求解器选择 ===
    solver::Symbol            # 求解器类型: :HiGHS 或 :COPT
end

# 引入初值相关模块（必须在 MILPModelParameters 定义之后）
include("HSOneStorageOneCompressorMILP/InitialSolution.jl")

"""
    MILPModelResult
MILP模型求解结果
"""
struct MILPModelResult
    objective::Float64
    C_operation::Float64      # 运行成本 元/天
    C_norm::Float64           # 正则项
    C_initial::Float64        # 初投资成本 元
    C_total::Float64          # 总现值 元
    P1::Vector{Float64}       # 热泵直接供热功率 kW
    P2::Vector{Float64}       # 蓄热供热功率 kW
    P3::Vector{Float64}       # 热泵向蓄热储热功率 kW
    P_el::Vector{Float64}     # 电加热补热功率 kW
    P_es::Vector{Float64}     # 电加热储热功率 kW
    Ts::Vector{Float64}       # 蓄热温度 ℃
    states::Vector{Int}       # 各时段状态编号 1..8
    COP1::Vector{Float64}     # 实际COP1
    COP2::Vector{Float64}     # 实际COP2
    COP3::Vector{Float64}     # 实际COP3
    COPw::Vector{Float64}     # 实际COPw
    C_heatpump_opt::Float64   # 优化后的热泵容量
    C_boiler_opt::Float64     # 优化后的电锅炉容量
    isFeasible::Bool
    gap::Float64              # MIP gap
    best_bound::Float64       # Best bound
end

"""
计算构建分段COP约束的数据

文档引用: 第2.2节 制热系数约束
"""
function getCOP_piecewise_data(
    designParameters,
    T_s_1g_list::Vector{T},
    T_s_2g_list::Vector{T},
    T_s_3g_list::Vector{T},
    T_s_wg_list::Vector{T},
    t_list::Vector{T}
) where {T<:Real}
    n = length(t_list)-1
    
    # COP函数定义
    function COP1(Te, Ts)
        if Ts + designParameters.dT_EvaporationStandard > designParameters.ThMax
            return 1.0
        end
        return designParameters.COPOverlapFunction(Te, max(designParameters.Tuse, Ts + designParameters.dT_EvaporationStandard))
    end
    COP2(Ts) = designParameters.COPWater(Ts - designParameters.dT_EvaporationStandard, designParameters.Tuse)
    function COP3(Te, Ts)
        if Ts + designParameters.dT_EvaporationStandard > designParameters.ThMax
            return 1.0
        end
        return designParameters.COPOverlapFunction(Te, Ts + designParameters.dT_EvaporationStandard)
    end
    function COPw(Ts)
        if Ts > designParameters.ThMax
            return 1.0
        end
        return designParameters.COPWater(designParameters.TCompressorIn, Ts)
    end
    
    # 分档数量
    m1 = length(T_s_1g_list) - 1
    m2 = length(T_s_2g_list) - 1
    m3 = length(T_s_3g_list) - 1
    mw = length(T_s_wg_list) - 1
    
    # COP分档值矩阵 [档位数, 时段数]
    COP1v = zeros(m1, n)
    COP2v = zeros(m2, n)
    COP3v = zeros(m3, n)
    COPwv = zeros(mw, n)
    
    # COP1取后一个温度（蓄热温度升高）
    for i in 1:n
        for j in 1:m1
            COP1v[j, i] = COP1(designParameters.TairFunction(t_list[i])-designParameters.dT_EvaporationStandard, T_s_1g_list[j])
        end
    end
    
    # COP2取前一个温度（蓄热温度降低）
    for i in 1:n
        for j in 1:m2
            COP2v[j, i] = COP2(T_s_2g_list[j])
        end
    end
    
    # COP3取后一个温度（蓄热温度升高）
    for i in 1:n
        for j in 1:m3
            COP3v[j, i] = COP3(designParameters.TairFunction(t_list[i])-designParameters.dT_EvaporationStandard, T_s_3g_list[j])
        end
    end
    
    # COPw直接在温度上离散
    for i in 1:n
        for j in 1:mw
            COPwv[j, i] = COPw(T_s_wg_list[j])
        end
    end
    
    # 温度分档界限矩阵 [时段数, 档位数-1]
    T1g = repeat(T_s_1g_list', n, 1)
    T2g = repeat(T_s_2g_list', n, 1)
    T3g = repeat(T_s_3g_list', n, 1)
    Twg = repeat(T_s_wg_list', n, 1)
    
    return m1, m2, m3, mw, COP1v, COP2v, COP3v, COPwv, T1g, T2g, T3g, Twg
end

"""
    generateAndSolve(::PressedWaterOneStorageOneCompressor_MILP, params::MILPModelParameters)

构建并求解MILP模型，严格按照文档实现所有约束和目标函数。

文档引用: 全文所有约束
"""
function generate_model(
    ::PressedWaterOneStorageOneCompressor_MILP,
    params::MILPModelParameters;
    callback::Union{Function, Nothing} = nothing
)
    n = params.n_segments
    
    # =============================================================
    # 大M参数 (文档1.2.1节) - 从外部传入
    # =============================================================
    M1 = params.M1  # 蓄热温度范围
    M2 = params.M2  # 温度约束大M
    M3 = params.M3  # 温度约束大M
    M4 = params.M4  # 档位温差最大值
    M5 = params.M5  # 功率大M
    
    m1, m2, m3, mw = params.m1, params.m2, params.m3, params.mw

    max_load = maximum(params.heatLoad)
    # =============================================================
    # 创建模型 - 根据求解器类型选择
    # =============================================================
    if params.solver == :COPT
        # 检查COPT是否可用
        copt_available = isdefined(@__MODULE__, :COPT)
        if !copt_available
            error("COPT求解器未安装，请先安装COPT.jl包，或将solver设置为:HiGHS")
        end
        model = direct_model(COPT.Optimizer())
        # COPT参数设置
        # set_silent(model)
        set_attribute(model, "TimeLimit", 60*10)
        set_attribute(model, "Presolve", 3)
        set_attribute(model, "Threads", 24)

        set_attribute(model, "CutLevel", 1)            # 增强割平面
        #set_attribute(model, "RootCutRounds", 20)      # 根节点多轮割
        #set_attribute(model, "StrongBranching", 1)     # 启用强分支
        set_attribute(model, "RelGap", 0.05)          # 设置一个合理的最优间隙，避免过度证明

        # 2. 根节点割平面强度与轮数
        # 作用: 专门控制根节点上的割。根节点割得好，能大幅提升初始下界。
        # 建议: 将强度设为 2，并将轮数从默认的少轮增加到 10 轮。
        set_attribute(model, "RootCutLevel", 3)
        #set_attribute(model, "RootCutRounds", 10)
        
        # 3. 搜索树中的割平面策略
        # 作用: 控制搜索树节点上的割。
        # 建议: 如果单节点 LP 求解变慢 (LPit/n 升高)，可以降低此参数。
        # 如果希望更激进地剪枝，可以设为 1 或 2。
        #set_attribute(model, "TreeCutLevel", 1)

        # 4. 节点割平面轮数 (可选)
        # 作用: 限制在搜索树节点上生成割的轮数。
        # 建议: 保持默认，或设为 2~3 以限制单节点开销。
        #set_attribute(model, "NodeCutRounds", 2)

        # 5. 启发式算法强度
        set_attribute(model, "PreRootHeurLevel", 3)
        set_attribute(model, "DivingHeurLevel", 3)
        set_attribute(model, "SubMipHeurLevel", 3)
        set_attribute(model, "FAPHeurLevel", 3)

    elseif params.solver == :HiGHS
        model = Model(HiGHS.Optimizer)
        set_silent(model)
        set_attribute(model, "time_limit", 300.0)
        set_attribute(model, "mip_rel_gap", 0.01)
    else
        error("不支持的求解器类型: $(params.solver)")
    end
    
    # =============================================================
    # 注册回调函数 (只对COPT求解器有效)
    # =============================================================
    if params.solver == :COPT && callback !== nothing
        # 创建wrapper回调函数，将model作为第三个参数传递给用户回调
        MOI.set(model, COPT.CallbackFunction(), (cb_data, cb_context) -> callback(cb_data, cb_context, model))
    end

    # 3.1 状态变量 s_{k,i} (文档2.1.1节 系统运行模式组合)
    # =============================================================
    """
    状态定义:
    - s1: 热泵直接供热
    - s2: 蓄热供热
    - s3: 热泵向蓄热储热
    - s4: 热泵和蓄热同时供热
    - s5: 热泵同时供热和储热
    - s6: 蓄热供热切换为热泵供热
    - s7: 热泵储热切换为热泵供热
    - s8: 热泵同时供热储热切换为热泵供热

    文档约束(2.1.1节) - SOS1形式:
    {s_1,i, s_2,i, ..., s_8,i} ∈ SOS1
    ∑_{k=1}^{8} s_k,i = 1
    """
    # 定义状态变量 (使用 SOS1 约束)
    @variable(model, 0 <= s[1:8, 1:n] <= 1)
    
    # 文档2.1.1节: 状态唯一性约束
    @constraint(model, state_unique[i=1:n], sum(s[k, i] for k=1:8) == 1)
    
    # SOS1约束: 每个时段的状态变量组
    # 使用 MOI.SOS1 约束，权重为 [1.0, 2.0, ..., 8.0]
    for i in 1:n
        # 创建该时段的SOS1约束 (匿名约束，不指定名称)
        @constraint(model, s[:, i] in MOI.SOS1([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]))
    end

    # --- 4.1.2 s3, s5温度约束 (文档2.1.2.2节) ---
    """
    T_s,i: 第i个时段初始时刻的蓄热温度
    文档2.6.2节: T_s,min ≤ T_s,i ≤ T_s,max
    """
    @variable(model, params.Tsmin <= Ts[1:n+1] <= params.Tsmax)
    """
    文档2.1.2.2节:
    T_s,i - T_s,i+1 ≤ M1·(1-s3,i)
    T_s,i+1 ≤ T_c,max - ΔT_s + M2·(1-s3,i)
    T_s,i - T_s,i+1 ≤ M1·(1-s5,i)
    T_s,i+1 ≤ T_c,max - ΔT_s + M2·(1-s5,i)
    """
    @constraint(model, s3_temp_drop[i=1:n], Ts[i] - Ts[i+1] <= M1 * (1 - s[3, i]))
    @constraint(model, s3_temp_upper[i=1:n], Ts[i+1] <= params.Tcmax - params.dTs + M2 * (1 - s[3, i]))
    @constraint(model, s5_temp_drop[i=1:n], Ts[i] - Ts[i+1] <= M1 * (1 - s[5, i]))
    @constraint(model, s5_temp_upper[i=1:n], Ts[i+1] <= params.Tcmax - params.dTs + M2 * (1 - s[5, i]))

    """
    文档2.1.2.3节:
    T_s,i+1 - T_s,i ≤ M1·(1-s4,i)
    T_u + ΔT_s - T_s,i+1 ≤ M3·(1-s4,i)
    """
    @constraint(model, s4_temp_rise[i=1:n], Ts[i+1] - Ts[i] <= M1 * (1 - s[4, i]))
    @constraint(model, s4_temp_lower[i=1:n], params.Tuse + params.dTs - Ts[i+1] <= M3 * (1 - s[4, i]))

    # --- δ_su: 蓄热温度高于T_u+ΔT_s标志 ---
    # --- δ̃: 蓄热温度从高到低过渡标志 ---
    """
    文档2.1.2.4节:
    δ_su,i = 1 when T_s,i ≥ T_u + ΔT_s
    δ̃_i = 1 when δ_su,i = 1 AND δ_su,i+1 = 0
    """
    @variable(model, delta_su[1:n], Bin)
    @variable(model, 0 <= delta_tilde[1:n] <= 1)

    """
    文档2.1.2.4节:
    T_s,i ≥ T_u + ΔT_s - M1·(1-δ_su,i)
    T_s,i ≤ T_u + ΔT_s + M1·δ_su,i
    """
    @constraint(model, delta_su_lower[i=1:n], Ts[i] >= params.Tuse + params.dTs - M1 * (1 - delta_su[i]))
    @constraint(model, delta_su_upper[i=1:n], Ts[i] <= params.Tuse + params.dTs + M1 * delta_su[i])
    
    # --- 4.1.5 δ̃约束 (文档2.1.2.4节) ---
    """
    文档2.1.2.4节:
    δ̃_i ≤ δ_su,i
    δ̃_i ≤ 1 - δ_su,i+1
    δ̃_i ≥ δ_su,i - δ_su,i+1
    δ̃_i ≥ 0
    """
    for i in 1:n-1
        @constraint(model, delta_tilde[i] <= delta_su[i])
        @constraint(model, delta_tilde[i] <= 1 - delta_su[i+1])
        @constraint(model, delta_tilde[i] >= delta_su[i] - delta_su[i+1])
    end
    
    """
    文档2.1.2.4节:
    ∑_{i=1}^{k} s6,i ≤ ∑_{i=1}^{k} δ̃_i, k=1,...,n-1
    ∑_{i=1}^{n} s6,i = ∑_{i=1}^{n} δ̃_i
    ∑_{i=k}^{n} s6,i ≤ ∑_{i=k}^{n} δ̃_i + 1, k=1,...,n
    """
    for k in 1:n-1
        @constraint(model, sum(s[6, i] for i=1:k) <= sum(delta_tilde[i] for i=1:k))
    end
    @constraint(model, sum(s[6, i] for i=1:n) == sum(delta_tilde[i] for i=1:n-1))
    for k in 1:n-1
        @constraint(model, sum(s[6, i] for i=k:n) <= sum(delta_tilde[i] for i=k:n) + 1)
    end

    """
    文档2.1.2.5节:
    T_s,i - T_s,i+1 ≤ M1·(1-s7,i)
    T_s,i - T_c,max + ΔT_s ≤ M2·(1-s7,i)
    T_s,i - T_s,i+1 ≤ M1·(1-s8,i)
    T_s,i - T_c,max + ΔT_s ≤ M2·(1-s8,i)
    T_{s,i}  >= T_u - M_1(1- s_{8,i})
    """
    @constraint(model, s7_temp_drop[i=1:n], Ts[i] - Ts[i+1] <= M1 * (1 - s[7, i]))
    @constraint(model, s7_temp_upper[i=1:n], Ts[i] - params.Tcmax + params.dTs <= M2 * (1 - s[7, i]))
    @constraint(model, s8_temp_drop[i=1:n], Ts[i] - Ts[i+1] <= M1 * (1 - s[8, i]))
    @constraint(model, s8_temp_upper[i=1:n], Ts[i] - params.Tcmax + params.dTs <= M2 * (1 - s[8, i]))
    @constraint(model, s8_temp_lower[i=1:n], Ts[i] >= params.Tuse - params.dTs - M1 * (1 - s[8, i]))


    """
    文档2.2节 - COP分档选择变量 (SOS1形式):
    {z_1,i, z_2,i, ..., z_m,i} ∈ SOS1
    ∑_{j=1}^{m} z_j,i = 1
    """
    # COP分档选择变量 z_{1,2,3,w,i,j} - 使用 SOS1 约束
    @variable(model, 0 <= z1[1:n, 1:m1] <= 1)
    @variable(model, 0 <= z2[1:n, 1:m2] <= 1)
    @variable(model, 0 <= z3[1:n, 1:m3] <= 1)
    @variable(model, 0 <= zw[1:n, 1:mw] <= 1)

    @variable(model, 0 <= w1[1:n, 1:m1] <= 1)
    @variable(model, 0 <= w3[1:n, 1:m3] <= 1)
    
    # 添加 SOS1 约束和唯一性约束
    for i in 1:n
        # COP1 分档 SOS1 约束 (匿名约束)
        @constraint(model, sum(z1[i, j] for j=1:m1) == 1)
        @constraint(model, z1[i, :] in MOI.SOS1(collect(1.0:m1)))
        
        # COP2 分档 SOS1 约束
        @constraint(model, sum(z2[i, j] for j=1:m2) == 1)
        @constraint(model, z2[i, :] in MOI.SOS1(collect(1.0:m2)))
        
        # COP3 分档 SOS1 约束
        @constraint(model, sum(z3[i, j] for j=1:m3) == 1)
        @constraint(model, z3[i, :] in MOI.SOS1(collect(1.0:m3)))
        
        # COPw 分档 SOS1 约束
        @constraint(model, sum(zw[i, j] for j=1:mw) == 1)
        @constraint(model, zw[i, :] in MOI.SOS1(collect(1.0:mw)))
    end
    
    # COP估计值 COP_{1e,2e,3e,we,i}
    @variable(model, 0 <= COP1e[1:n] <= 21.0)
    @variable(model, 0 <= COP2e[1:n] <= 21.0)
    @variable(model, 0 <= COP3e[1:n] <= 21.0)
    @variable(model, 0 <= COPwe[1:n] <= 21.0)
    
    # COP实际值 COP_{1,2,3,w,i}
    @variable(model, 0 <= COP1[1:n] <= 21.0)
    @variable(model, 0 <= COP2[1:n] <= 21.0)
    @variable(model, 0 <= COP3[1:n] <= 21.0)
    @variable(model, 0 <= COPw[1:n] <= 21.0)

    """
    文档2.2.1节: y1,i = s5 || s8
    y1 ≤ 1
    y1 ≥ s5
    y1 ≥ s8
    y1 ≥ 0
    """
    @variable(model, 0 <= y1[1:n] <= 1)
    @constraint(model, y1_upper[i=1:n], y1[i] <= s[5, i] + s[8, i])
    @constraint(model, y1_s5[i=1:n], y1[i] >= s[5, i])
    @constraint(model, y1_s8[i=1:n], y1[i] >= s[8, i])

    
    """
    文档2.2.1节:
    ∑_{j=1}^{m1} z1,i,j = 1
    COP1e,i = ∑_{j=1}^{m1} COP1v,i,j · z1,i,j
    T_s,i ≤ T_s,1g,i,1 + M1·(1-z1,i,1)
    T_s,i ≥ T_s,1g,i,m1-1 - M1·(1-z1,i,m1)
    对于 j=2,...,m1-1:
        T_s,i ≥ T_s,1g,i,j-1 + M4·(1-z1,i,j)
        T_s,i ≤ T_s,1g,i,j + M4·(1-z1,i,j)
    """    
    for i in 1:n
        #@constraint(model, COP1e[i] == sum(params.COP1v[j, i] * z1[i, j] for j=1:m1))
        if m1 > 1
            @constraint(model, Ts[i] <= params.T1g[i, 2] + M1 * (1 - z1[i, 1]))
            @constraint(model, Ts[i] >= params.T1g[i, m1] - M1 * (1 - z1[i, m1]))            
            for j in 2:m1-1
                @constraint(model, Ts[i] >= params.T1g[i, j] - M1 * (1 - z1[i, j]))
                @constraint(model, Ts[i] <= params.T1g[i, j+1] + M1 * (1 - z1[i, j]))
            end
        end
    end

    """
    文档2.2.1节:w1 = y1 * z1
    w1,i,j ≤ y1,i
    w1,i,j ≤ z1,i,j
    w1,i,j ≥ y1,i + z1,i,j - 1
    w1,i,j ≥ 0
    """
    for i in 1:n, j in 1:m1
        @constraint(model, w1[i, j] <= y1[i])
        @constraint(model, w1[i, j] <= z1[i, j])
        @constraint(model, w1[i, j] >= y1[i] + z1[i, j] - 1)
    end

    """
    文档2.2.1节:
    COP1,i = (1-y1,i)·COP_ca(T_e,T_u) + ∑_{j=1}^{m1} COP1v,i,j·w1,i,j
    """
    for i in 1:n
        @constraint(model, COP1[i] == (1 - y1[i]) * params.COPca[i] + sum(params.COP1v[j, i] * w1[i, j] for j=1:m1))
    end


    """
    文档2.2.2节:
    ∑_{j=1}^{m2} z2,i,j = 1
    COP2e,i = ∑_{j=1}^{m2} COP2v,i,j · z2,i,j
    T_s,i ≤ T_s,2g,i,1 + M1·(1-z2,i,1)
    T_s,i ≥ T_s,2g,i,m2-1 - M1·(1-z2,i,m2)
    对于 j=2,...,m2-1:
        T_s,i ≥ T_s,2g,i,j-1 + M4·(1-z2,i,j)
        T_s,i ≤ T_s,2g,i,j + M4·(1-z2,i,j)
    COP2,i = COP2e,i
    """
    for i in 1:n
        @constraint(model, COP2e[i] == sum(params.COP2v[j, i] * z2[i, j] for j=1:m2))
        if m2 > 1
            @constraint(model, Ts[i] <= params.T2g[i, 2] + M1 * (1 - z2[i, 1]))
            @constraint(model, Ts[i] >= params.T2g[i, m2] - M1 * (1 - z2[i, m2]))
            for j in 2:m2-1
                @constraint(model, Ts[i] >= params.T2g[i, j] - M1 * (1 - z2[i, j]))
                @constraint(model, Ts[i] <= params.T2g[i, j+1] + M1 * (1 - z2[i, j]))
            end
        end
        @constraint(model, COP2[i] == COP2e[i])
    end

    """
    文档2.2.3节:
    ∑_{j=1}^{m3} z3,i,j = 1
    COP3e,i = ∑_{j=1}^{m3} COP3v,i,j · z3,i,j
    温度分档约束同COP1
    """
    for i in 1:n
        #@constraint(model, COP3e[i] == sum(params.COP3v[j, i] * z3[i, j] for j=1:m3))
        if m3 > 1
            @constraint(model, Ts[i] <= params.T3g[i, 2] + M1 * (1 - z3[i, 1]))
            @constraint(model, Ts[i] >= params.T3g[i, m3] - M1 * (1 - z3[i, m3]))
            for j in 2:m3-1
                @constraint(model, Ts[i] >= params.T3g[i, j] - M1 * (1 - z3[i, j]))
                @constraint(model, Ts[i] <= params.T3g[i, j+1] + M1 * (1 - z3[i, j]))
            end
        end
    end
    """
    文档2.2.3节:w3 = (1-y1) * z3
    w3,i,j ≤ 1-y1,i
    w3,i,j ≤ z3,i,j
    w3,i,j ≥ z3,i,j - y1,i
    w3,i,j ≥ 0
    """
    for i in 1:n, j in 1:m3
        @constraint(model, w3[i, j] <= 1 - y1[i])
        @constraint(model, w3[i, j] <= z3[i, j])
        @constraint(model, w3[i, j] >= z3[i, j] - y1[i])
    end
    
    """
    文档2.2.3节:
    COP3,i = ∑_{j=1}^{m3} COP3v,i,j·w3,i,j + ∑_{j=1}^{m1} COP1v,i,j·w1,i,j
    """
    for i in 1:n
        @constraint(model, COP3[i] == sum(params.COP3v[j, i] * w3[i, j] for j=1:m3) + sum(params.COP1v[j, i] * w1[i, j] for j=1:m1))
    end

    """
    文档2.2.4节:
    ∑_{j=1}^{mw} zw,i,j = 1
    COPwe,i = ∑_{j=1}^{mw} COPwv,i,j · zw,i,j
    温度分档约束同COP2
    COPw,i = COPwe,i
    """
    for i in 1:n
        @constraint(model, COPwe[i] == sum(params.COPwv[j, i] * zw[i, j] for j=1:mw))
        if mw > 1
            @constraint(model, Ts[i] <= params.Twg[i, 2] + M1 * (1 - zw[i, 1]))
            @constraint(model, Ts[i] >= params.Twg[i, mw] - M1 * (1 - zw[i, mw]))
            for j in 2:mw-1
                @constraint(model, Ts[i] >= params.Twg[i, j] - M1 * (1 - zw[i, j]))
                @constraint(model, Ts[i] <= params.Twg[i, j+1] + M1 * (1 - zw[i, j]))
            end
        end
        @constraint(model, COPw[i] == COPwe[i])
    end
    
    
    # =============================================================
    # 3.2 功率变量 P_{k,i,j}^s 和 P_{k,i,j}^e (文档2.3.1节)
    # =============================================================
    """
    文档2.3.1节变量描述:
    P_{k,i,j}^s: 第i个时段，状态k在模式j的起始功率
    P_{k,i,j}^e: 第i个时段，状态k在模式j的结束功率
    
    模式j定义:
    j=1: 热泵直接供热模式
    j=2: 蓄热供热模式
    j=3: 热泵向蓄热储热模式
    
    功率约束表(文档2.3.1节):
    """
    @variable(model, 0 <= P_k_s[1:8, 1:n, 1:3] <= 5)  # P_{k,i,j}^s
    @variable(model, 0 <= P_k_e[1:8, 1:n, 1:3] <= 5)  # P_{k,i,j}^e
    
    # 根据功率约束表(文档2.3.1节)，某些功率恒为0
    # k=1: 只有P_1^s非零
    @constraint(model, P1_zero_s[j=2:3, i=1:n], P_k_s[1, i, j] == 0)
    @constraint(model, P1_zero_e[j=1:3, i=1:n], P_k_e[1, i, j] == 0)
    
    # k=2: 只有P_2^s非零
    @constraint(model, P2_zero_s[i=1:n], P_k_s[2, i, 1] == 0)
    @constraint(model, P2_zero_e[i=1:n], P_k_e[2, i, 2] == 0)
    @constraint(model, P2_zero_e2[i=1:n], P_k_e[2, i, 1] == 0)
    @constraint(model, P2_zero_e3[i=1:n], P_k_e[2, i, 3] == 0)
    @constraint(model, P2_zero_s3[i=1:n], P_k_s[2, i, 3] == 0)
    
    # k=3: 只有P_3^s非零
    @constraint(model, P3_zero_s[i=1:n], P_k_s[3, i, 1] == 0)
    @constraint(model, P3_zero_s2[i=1:n], P_k_s[3, i, 2] == 0)
    @constraint(model, P3_zero_e[i=1:n], P_k_e[3, i, 3] == 0)
    @constraint(model, P3_zero_e1[i=1:n], P_k_e[3, i, 1] == 0)
    @constraint(model, P3_zero_e2[i=1:n], P_k_e[3, i, 2] == 0)
    
    # k=4: 只有P_1^s和P_2^s非零
    @constraint(model, P4_zero_s3[i=1:n], P_k_s[4, i, 3] == 0)
    @constraint(model, P4_zero_e[i=1:n, j=1:3], P_k_e[4, i, j] == 0)
    
    # k=5: 只有P_1^s和P_3^s非零
    @constraint(model, P5_zero_s2[i=1:n], P_k_s[5, i, 2] == 0)
    @constraint(model, P5_zero_e[i=1:n, j=1:3], P_k_e[5, i, j] == 0)
    
    # k=6: P_2^s和P_1^e非零
    @constraint(model, P6_zero_s1[i=1:n], P_k_s[6, i, 1] == 0)
    @constraint(model, P6_zero_s3[i=1:n], P_k_s[6, i, 3] == 0)
    @constraint(model, P6_zero_e2[i=1:n], P_k_e[6, i, 2] == 0)
    @constraint(model, P6_zero_e3[i=1:n], P_k_e[6, i, 3] == 0)
    
    # k=7: P_3^s和P_1^e非零
    @constraint(model, P7_zero_s1[i=1:n], P_k_s[7, i, 1] == 0)
    @constraint(model, P7_zero_s2[i=1:n], P_k_s[7, i, 2] == 0)
    @constraint(model, P7_zero_e3[i=1:n], P_k_e[7, i, 3] == 0)
    @constraint(model, P7_zero_e2[i=1:n], P_k_e[7, i, 2] == 0)
    
    # k=8: P_1^s, P_3^s, P_1^e非零
    @constraint(model, P8_zero_s2[i=1:n], P_k_s[8, i, 2] == 0)
    @constraint(model, P8_zero_e2[i=1:n], P_k_e[8, i, 2] == 0)
    @constraint(model, P8_zero_e3[i=1:n], P_k_e[8, i, 3] == 0)
    
    # =============================================================
    # 3.3 电加热功率变量 (文档1.2.3节 设备运行相关)
    # =============================================================
    """
    P_e^l: 电加热为用热负荷补热的功率
    P_e^s: 电加热为蓄热储热的功率
    """
    @variable(model, 0 <= P_el[1:n] <= 5)  # P_{e,i}^l
    @variable(model, 0 <= P_es[1:n] <= 5)  # P_{e,i}^s

    # 计算总功率 P_{i,j}
    @expression(model, P_total[i=1:n, j=1:3], sum(P_k_s[k, i, j] + P_k_e[k, i, j] for k=1:8))
    @expression(model, P_total_4[i=1:n], P_el[i] + P_es[i])  # P_{i,4}

    """
    文档2.3.2节:
    u_k,i^1 = s_k,i * (1-y1,i)
    u_k,i^2 = u_k,i^1 * P_k,i,1^s
    u_8,i^3 = s_8,i * P_8,i,1^e
    """
    @variable(model, 0 <= u1[1:8, 1:n] <= 1)      # u_k,i^1
    @variable(model, 0 <= u2[1:8, 1:n] <= 5)      # u_k,i^2
    @variable(model, 0 <= u3[1:n] <= 5)           # u_8,i^3 (仅状态8)

    """
    文档2.3.2节和2.4节:
    v_k,i,j^1 = w1,i,j * s_k,i
    v_k,i,j^2 = v_k,i,j^1 * P_k,i,1^s
    v_k,i,j^3 = z2,i,j * s_k,i
    v_k,i,j^4 = v_k,i,j^3 * P_k,i,2^s
    v_k,i,j^5 = s_k,i * P_k,i,3^s
    v_k,i,j^6 = v_k,i,j^5 * w3,i,j
    v_k,i,j^7 = v_k,i,j^5 * w1,i,j
    v_k,i,j^8 = s_k,i * P_k,i,2^s
    v_k,i,j^9 = v_k,i,j^8 * z2,i,j
    """
    @variable(model, 0 <= v1[1:8, 1:n, 1:m1] <= 1)  # v_k,i,j^1
    @variable(model, 0 <= v2[1:8, 1:n, 1:m1] <= 5)  # v_k,i,j^2
    @variable(model, 0 <= v3[1:8, 1:n, 1:m2] <= 1)  # v_k,i,j^3
    @variable(model, 0 <= v4[1:8, 1:n, 1:m2] <= 5)  # v_k,i,j^4
    @variable(model, 0 <= v5[1:8, 1:n] <= 5)        # v_k,i^5 (简化为不依赖j)
    @variable(model, 0 <= v6[1:8, 1:n, 1:m3] <= 5)  # v_k,i,j^6
    @variable(model, 0 <= v7[1:8, 1:n, 1:m1] <= 5)  # v_k,i,j^7
    @variable(model, 0 <= v8[1:8, 1:n] <= 5)        # v_k,i^8
    @variable(model, 0 <= v9[1:8, 1:n, 1:m2] <= 5)  # v_k,i,j^9

    # u1 = s * (1-y1)
    for k in 1:8, i in 1:n
        @constraint(model, u1[k, i] <= s[k, i])
        @constraint(model, u1[k, i] <= 1 - y1[mod1(i+1,n)])
        @constraint(model, u1[k, i] >= s[k, i] - y1[mod1(i+1,n)])
    end

    # v1 = w1 * s
    for k in 1:8, i in 1:n, j in 1:m1
        @constraint(model, v1[k, i, j] <= s[k, i])
        @constraint(model, v1[k, i, j] <= w1[mod1(i+1,n), j])
        @constraint(model, v1[k, i, j] >= s[k, i] + w1[mod1(i+1,n), j] - 1)
    end

    # --- 4.3.2 u2约束: u2 = u1 * P^s (文档2.3.2节) ---
    """
    文档2.3.2节:
    u_k,i^2 ≤ M5·u_k,i^1
    u_k,i^2 ≥ P_k,i,1^s - M5·(1-u_k,i^1)
    u_k,i^2 ≤ P_k,i,1^s
    u_k,i^2 ≥ 0
    """
    for k in vcat(1:5,8), i in 1:n
        @constraint(model, u2[k, i] <= M5 * u1[k, i])
        @constraint(model, u2[k, i] >= P_k_s[k, i, 1] - M5 * (1 - u1[k, i]))
        @constraint(model, u2[k, i] <= P_k_s[k, i, 1])
    end

    # k=6
    for k in 6:7, i in 1:n
        @constraint(model, u2[k, i] <= M5 * u1[k, i])
        @constraint(model, u2[k, i] >= P_k_e[k, i, 1] - M5 * (1 - u1[k, i]))
        @constraint(model, u2[k, i] <= P_k_e[k, i, 1])
    end

    # --- 4.3.4 v2约束: v2 = v1 * P^s (文档2.3.2节) ---
    """
    文档2.3.2节:
    v_k,i,j^2 ≤ M5·v_k,i,j^1
    v_k,i,j^2 ≥ P_k,i,1^s - M5·(1-v_k,i,j^1)
    v_k,i,j^2 ≤ P_k,i,1^s
    v_k,i,j^2 ≥ 0
    """
    for k in vcat(1:5,8), i in 1:n, j in 1:m1
        @constraint(model, v2[k, i, j] <= M5 * v1[k, i, j])
        @constraint(model, v2[k, i, j] >= P_k_s[k, i, 1] - M5 * (1 - v1[k, i, j]))
        @constraint(model, v2[k, i, j] <= P_k_s[k, i, 1])
    end

    for k in 6:7, i in 1:n, j in 1:m1
        @constraint(model, v2[k, i, j] <= M5 * v1[k, i, j])
        @constraint(model, v2[k, i, j] >= P_k_e[k, i, 1] - M5 * (1 - v1[k, i, j]))
        @constraint(model, v2[k, i, j] <= P_k_e[k, i, 1])
    end

    """
    文档2.3.2节:v3 = z2 * s
    v_k,i,j^3 ≤ s_k,i
    v_k,i,j^3 ≤ z2,i,j
    v_k,i,j^3 ≥ s_k,i + z2,i,j - 1
    v_k,i,j^3 ≥ 0
    """
    for k in 1:8, i in 1:n, j in 1:m2
        @constraint(model, v3[k, i, j] <= s[k, i])
        @constraint(model, v3[k, i, j] <= z2[mod1(i+1,n), j])
        @constraint(model, v3[k, i, j] >= s[k, i] + z2[mod1(i+1,n), j] - 1)
    end
    
    """
    文档2.3.2节:v4 = v3 * P^s
    v_k,i,j^4 ≤ M5·v_k,i,j^3
    v_k,i,j^4 ≥ P_k,i,2^s - M5·(1-v_k,i,j^3)
    v_k,i,j^4 ≤ P_k,i,2^s
    v_k,i,j^4 ≥ 0
    """
    for k in 1:8, i in 1:n, j in 1:m2
        @constraint(model, v4[k, i, j] <= M5 * v3[k, i, j])
        @constraint(model, v4[k, i, j] >= P_k_s[k, i, 2] - M5 * (1 - v3[k, i, j]))
        @constraint(model, v4[k, i, j] <= P_k_s[k, i, 2])
    end

    """
    文档2.3.5节:u3 = s8 * P_e
    u_8,i^3 ≤ M5·s_8,i
    u_8,i^3 ≥ P_8,i,1^e - M5·(1-s_8,i)
    u_8,i^3 ≤ P_8,i,1^e
    u_8,i^3 ≥ 0
    """
    for i in 1:n
        @constraint(model, u3[i] <= M5 * s[8, i])
        @constraint(model, u3[i] >= P_k_e[8, i, 1] - M5 * (1 - s[8, i]))
        @constraint(model, u3[i] <= P_k_e[8, i, 1])
    end

    """
    文档2.3.5节:
    L_i ≤ ∑_{k=1}^{7} (COP_ca·u_k,i^2 + ∑_{j=1}^{m1} COP1v,j·v_k,i,j^2 + ∑_{j=1}^{m2} COP2v,j·v_k,i,j^4)
         + ∑_{j=1}^{m1} COP3v,j·v_{8,i,j}^2 + COP_ca·u_{8,i}^3 + P_e,i^l

    注意: 对于状态8, v_{8,i,j}^2的定义与前7个状态不同
    - k=1..7: v_{k,i,j}^2 = v_{k,i,j}^1 · P_{k,i,1}^s (起始功率)
    - k=8: v_{8,i,j}^2 = v_{8,i,j}^1 · P_{8,i,1}^s (前段供热功率,也是起始功率)
    """
    heat_k17_expr = @expression(model, [i=1:n],
        sum(
            params.COPca[i] * u2[k, i] +
            sum(params.COP1v[j, mod1(i+1,n)] * v2[k, i, j] for j=1:m1) +
            sum(params.COP2v[j, mod1(i+1,n)] * v4[k, i, j] for j=1:m2)
            for k=1:7
        )
    )
    heat_k8_expr = @expression(model, [i=1:n],
        sum(params.COP3v[j, mod1(i+1,n)] * v2[8, i, j] for j=1:m1) + params.COPca[i] * u3[i]
    )

    for i in 1:n
        @constraint(model, params.heatLoad[i] <= heat_k17_expr[i] + heat_k8_expr[i] + P_el[i])
    end

    """
    文档2.4节:v5 = s * P_3^s
    v_k,i^5 ≤ M5·s_k,i
    v_k,i^5 ≥ P_k,i,3^s - M5·(1-s_k,i)
    v_k,i^5 ≤ P_k,i,3^s
    v_k,i^5 ≥ 0
    """
    for k in 1:8, i in 1:n
        @constraint(model, v5[k, i] <= M5 * s[k, i])
        @constraint(model, v5[k, i] >= P_k_s[k, i, 3] - M5 * (1 - s[k, i]))
        @constraint(model, v5[k, i] <= P_k_s[k, i, 3])
    end
    
    """
    文档2.4节:v6 = v5 * w3
    v_k,i,j^6 ≤ M5·w3,i,j
    v_k,i,j^6 ≥ v_k,i^5 - M5·(1-w3,i,j)
    v_k,i,j^6 ≤ v_k,i^5
    v_k,i,j^6 ≥ 0
    """
    for k in 1:8, i in 1:n, j in 1:m3
        @constraint(model, v6[k, i, j] <= M5 * w3[mod1(i+1,n), j])
        @constraint(model, v6[k, i, j] >= v5[k, i] - M5 * (1 - w3[mod1(i+1,n), j]))
        @constraint(model, v6[k, i, j] <= v5[k, i])
    end

    """
    文档2.4节:v7 = v5 * w1
    v_k,i,j^7 ≤ M5·w1,i,j
    v_k,i,j^7 ≥ v_k,i^5 - M5·(1-w1,i,j)
    v_k,i,j^7 ≤ v_k,i^5
    v_k,i,j^7 ≥ 0
    """
    for k in 1:8, i in 1:n, j in 1:m1
        @constraint(model, v7[k, i, j] <= M5 * w1[mod1(i+1,n), j])
        @constraint(model, v7[k, i, j] >= v5[k, i] - M5 * (1 - w1[mod1(i+1,n), j]))
        @constraint(model, v7[k, i, j] <= v5[k, i])
    end
    
    """
    文档2.4节:v8 = s * P_2^s
    v_k,i^8 ≤ M5·s_k,i
    v_k,i^8 ≥ P_k,i,2^s - M5·(1-s_k,i)
    v_k,i^8 ≤ P_k,i,2^s
    v_k,i^8 ≥ 0
    """
    for k in 1:8, i in 1:n
        @constraint(model, v8[k, i] <= M5 * s[k, i])
        @constraint(model, v8[k, i] >= P_k_s[k, i, 2] - M5 * (1 - s[k, i]))
        @constraint(model, v8[k, i] <= P_k_s[k, i, 2])
    end
    

    """
    文档2.4节: v9 = v8 * z2
    v_k,i,j^9 ≤ M5·z2,i,j
    v_k,i,j^9 ≥ v_k,i^8 - M5·(1-z2,i,j)
    v_k,i,j^9 ≤ v_k,i^8
    v_k,i,j^9 ≥ 0
    """
    for k in 1:8, i in 1:n, j in 1:m2
        @constraint(model, v9[k, i, j] <= M5 * z2[i, j])
        @constraint(model, v9[k, i, j] >= v8[k, i] - M5 * (1 - z2[i, j]))
        @constraint(model, v9[k, i, j] <= v8[k, i])
    end

    """
    文档2.4节:
    c_p·m · (T_s,i+1 - T_s,i) / Δt_i = ∑_{k=1}^{8} (
        ∑_{j=1}^{m3} COP3v,j·v_k,i,j^6
        + ∑_{j=1}^{m1} COP1v,j·v_k,i,j^7
        - ∑_{j=1}^{m2} (COP2v,j-1)·v_k,i,j^9
    )
    + P_es,i^l
    """
    for i in 1:n
        heat_change = sum(
            sum(params.COP3v[j, mod1(i+1,n)] * v6[k, i, j] for j=1:m3) +
            sum(params.COP1v[j, mod1(i+1,n)] * v7[k, i, j] for j=1:m1) -
            sum((params.COP2v[j, mod1(i+1,n)] - 1) * v9[k, i, j] for j=1:m2)
            for k=1:8
        )
        @constraint(model, params.cpm * (Ts[i+1] - Ts[i]) / params.dt_list[i] == heat_change + P_es[i])
    end

    """
    文档2.5节:
    T_s,1 = T_s,n+1
    """
    @constraint(model, periodic_boundary, Ts[1] == Ts[n+1])

    """
    文档3.2节:
    C_heatpump: 热泵容量 kW
    C_boiler: 电锅炉容量 kW
    """
    if params.optimizeCapacity
        @variable(model, 0.0 <= C_heatpump <= 2.0*max_load)
        @variable(model, 0.0 <= C_boiler <= 5.5*max_load)
    else
        # 固定容量
        C_heatpump = params.C_heatpump
        C_boiler = params.C_boiler
    end

    """
    文档2.6.1节:
    P_e,i^l + P_e,i^s ≥ 0
    P_e,i^l + P_e,i^s ≤ C_boiler
    """
    for i in 1:n
        @constraint(model, P_el[i] + P_es[i] >= 0)
        if params.optimizeCapacity
            @constraint(model, P_el[i] + P_es[i] <= C_boiler)
        else
            @constraint(model, P_el[i] + P_es[i] <= params.C_boiler)
        end
    end

    """
    文档2.6.3节:
    P_i,1 + P_i,3 ≤ C_heatpump / COP_design
    
    注意: P_i,j = ∑_{k=1}^{8} (P_k,i,j^s + P_k,i,j^e)
    """
    for i in 1:n
        P_i1 = sum(P_k_s[k, i, 1] + P_k_e[k, i, 1] for k=1:8)
        P_i3 = sum(P_k_s[k, i, 3] + P_k_e[k, i, 3] for k=1:8)
        if params.optimizeCapacity
            @constraint(model, P_total[i, 1] + P_total[i, 3] <= C_heatpump / params.COPdesign[i])
        else
            @constraint(model, P_total[i, 1] + P_total[i, 3] <= params.C_heatpump / params.COPdesign[i])
        end
    end

    """
    文档3.1.1节:
    C_operation = ∑_{i=1}^{n} c_g,i · ∑_{j=1}^{4} P_{i,j} · Δt_i
    
    其中 P_{i,j} 定义如下:
    - j=1,2,3: 热泵模式功率
    - j=4: 电加热功率
    """
    @expression(model, C_operation,
        sum(params.segmentTariff[i] * sum(P_total[i, j] for j=1:3) * params.dt_list[i] for i=1:n) +
        sum(params.segmentTariff[i] * P_total_4[i] * params.dt_list[i] for i=1:n)
    )

     """
    文档3.1.2节:
    ΔP_{i,j}^positive: 功率正向变化量 (j=1,2,3,4)
    ΔP_{i,j}^negative: 功率负向变化量 (j=1,2,3,4)
    
    其中 j 对应四类功率:
    - j=1: 热泵直接供热功率 P_{i,1}
    - j=2: 蓄热供热功率 P_{i,2}
    - j=3: 热泵向蓄热储热功率 P_{i,3}
    - j=4: 电加热功率 P_{i,4}
    """
    @variable(model, 0 <= dP_pos[1:4, 1:n] <= 5)  # ΔP_{i,j}^positive
    @variable(model, 0 <= dP_neg[1:4, 1:n] <= 5)  # ΔP_{i,j}^negative

    """
    文档3.1.2节:
    C_norm = ∑_{i=1}^{n} ∑_{j=1}^{4}|P_{i+1,j} - P_{i,j}|
    
    线性化约束:
    P_{i+1,j} - P_{i,j} = ΔP_{i,j}^positive - ΔP_{i,j}^negative
    ΔP_{i,j}^positive ≥ 0
    ΔP_{i,j}^negative ≥ 0
    
    其中 P_{i,j} 定义如下 (文档2.3.1节):
    - j=1: 热泵直接供热功率 P_{i,1} = ∑_{k=1}^{8} (P_{k,i,1}^s + P_{k,i,1}^e)
    - j=2: 蓄热供热功率 P_{i,2} = ∑_{k=1}^{8} (P_{k,i,2}^s + P_{k,i,2}^e)
    - j=3: 热泵向蓄热储热功率 P_{i,3} = ∑_{k=1}^{8} (P_{k,i,3}^s + P_{k,i,3}^e)
    - j=4: 电加热功率 P_{i,4} = P_{e,i}^l + P_{e,i}^s
    """
    # 变差约束: P_{i+1,j} - P_{i,j} = ΔP_{i,j}^positive - ΔP_{i,j}^negative
    for j in 1:4
        for i in 1:n
            i_next = mod1(i+1,n)
            if j <= 3
                @constraint(model, P_total[i_next, j] - P_total[i, j] == dP_pos[j, i] - dP_neg[j, i])
            else
                @constraint(model, P_total_4[i_next] - P_total_4[i] == dP_pos[4, i] - dP_neg[4, i])
            end
        end
    end

    """
    文档3.1.2节:
    C_norm = ∑_{i=1}^{n} ∑_{j=1}^{4} (ΔP_{i,j}^positive + ΔP_{i,j}^negative)
    
    其中 j=1,2,3,4 对应四类功率:
    - j=1: 热泵直接供热功率
    - j=2: 蓄热供热功率
    - j=3: 热泵向蓄热储热功率
    - j=4: 电加热功率
    """
    @expression(model, C_norm,
        sum(dP_pos[j, i] + dP_neg[j, i] for j=1:4, i=1:n)
    )

    @expression(model, C_norm_2,
        sum(s[k,i]*k for k=1:8, i=1:n)
    )

    """
    文档3.2节:
    C_initial = c_heatpump·C_heatpump + c_storage·C_storage + c_boiler·C_boiler

    注意:
    - C_storage是常数,与cpm成正比: C_storage = cpm·(T_s,max - T_u)/3600
    - 当optimizeCapacity=false时,C_heatpump和C_boiler也是常数
    """
    if params.optimizeCapacity
        @expression(model, C_initial,
            params.c_heatpump * C_heatpump +
            params.c_storage * params.C_storage +
            params.c_boiler * C_boiler
        )
    else
        # 固定容量时初投资为常数，不影响优化
        @expression(model, C_initial,
            params.c_heatpump * params.C_heatpump +
            params.c_storage * params.C_storage +
            params.c_boiler * params.C_boiler
        )
    end

    """
    文档3.3节:
    C_total = C_initial + C_operation·D·(1-(1+r)^{-N})/r + λ_norm'·C_norm

    注意:
    - D: 年运行天数
    - r: 折现率
    - N: 运行年限
    - λ_norm' = λ_norm·D·(1-(1+r)^{-N})/r
    """
    annuity_factor = params.D * (1 - (1+params.r)^(-params.N)) / params.r
    lambda_norm_pv = params.lambdaNorm * annuity_factor

    @objective(model, Min, C_initial + C_operation * annuity_factor + lambda_norm_pv * (C_norm+C_norm_2))

    return model

end
function solve_model(::PressedWaterOneStorageOneCompressor_MILP,model,params; initial_solution = nothing)

    n = params.n_segments
    # =============================================================
    # 6. 设置初值（如果提供）
    # =============================================================
    
    if initial_solution !== nothing
        setInitialSolution(model, initial_solution, params)
    end

    # =============================================================
    # 7. 求解
    # =============================================================
    
    # 打印模型统计信息
    num_vars = num_variables(model)
    num_bin = count(is_binary, all_variables(model))        # 二进制变量
    num_cont = num_vars - num_bin      # 连续变量（包含SOS1中的变量）
    n_constraints = num_constraints(model; count_variable_in_set_constraints = false)
    num_sos1 = num_constraints(model; count_variable_in_set_constraints = true) - n_constraints
    println("\n模型统计:")
    println("  变量数量: $num_vars")
    println("    - 连续变量: $num_cont")
    println("    - 二进制变量: $num_bin")
    println("  约束数量: $n_constraints")
    println("  SOS1约束: $num_sos1")

    optimize!(model)
    
    isFeasible = primal_status(model) in [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]
    
    if !isFeasible
        return MILPModelResult(
            9999.0, 9999.0, 9999.0, 9999.0, 9999.0,
            fill(9999.0, n), fill(9999.0, n), fill(9999.0, n),
            fill(9999.0, n), fill(9999.0, n),
            fill(9999.0, n+1), fill(-1, n),
            fill(9999.0, n), fill(9999.0, n), fill(9999.0, n), fill(9999.0, n),
            9999.0, 9999.0,
            false,
            9999.0, 9999.0
        ), model
    end
    
    # 提取结果
    states = [argmax(value.(model[:s][:, i])) for i in 1:n]

    P1_result = [value(sum(model[:P_k_s][k, i, 1] + model[:P_k_e][k, i, 1] for k=1:8)) for i in 1:n]
    P2_result = [value(sum(model[:P_k_s][k, i, 2] + model[:P_k_e][k, i, 2] for k=1:8)) for i in 1:n]
    P3_result = [value(sum(model[:P_k_s][k, i, 3] + model[:P_k_e][k, i, 3] for k=1:8)) for i in 1:n]

    C_heatpump_opt = params.optimizeCapacity ? value(model[:C_heatpump]) : params.C_heatpump
    C_boiler_opt = params.optimizeCapacity ? value(model[:C_boiler]) : params.C_boiler

    # 计算总现值
    annuity_factor = params.D * (1 - (1+params.r)^(-params.N)) / params.r
    C_total_value = value(model[:C_initial]) + value(model[:C_operation]) * annuity_factor + params.lambdaNorm * annuity_factor * value(model[:C_norm])

    return MILPModelResult(
        objective_value(model),
        value(model[:C_operation]),
        value(model[:C_norm]),
        value(model[:C_initial]),
        C_total_value,
        P1_result, P2_result, P3_result,
        value.(model[:P_el]), value.(model[:P_es]),
        value.(model[:Ts]),
        states,
        value.(model[:COP1]), value.(model[:COP2]), value.(model[:COP3]), value.(model[:COPw]),
        C_heatpump_opt, C_boiler_opt,
        true,
        relative_gap(model),
        objective_bound(model)
    ), model
end

# 导出结构体和函数
export MILPModelParameters, MILPModelResult, getCOP_piecewise_data, generate_model, solve_model
export PressedWaterOneStorageOneCompressor_MILP
