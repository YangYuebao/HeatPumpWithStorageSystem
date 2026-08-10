# =============================================================================
# situation32 结果汇总绘图脚本
# 使用 situation32 的经济性结果绘制最优运行成本曲线图：
#   横轴 = 蓄热容量
#   纵轴 = 最小运行成本（对电锅炉容量维度取最小）
#   每条曲线 = 固定热泵容量
# 运行成本 = 每日运行成本 × 目标基准电价（dailyOperationCost × base_price）
# =============================================================================
using Plots
using JSON3

include(joinpath(pwd(), "tools", "plottool.jl"))

situation = "situation32"
file_path0 = joinpath(pwd(), "calculations", situation)

# 最优运行成本曲线图
println("\n========== situation32 结果汇总绘图 ==========")
num_oc = plot_optimal_operating_cost(
    file_path0;
    params_key = "",                   # 参数在JSON顶层
    econ_key = "economicResults",      # 经济数据在 economicResults
    cost_field = "dailyOperationCost", # 每日运行成本
    base_price = 1.0,                  # 目标基准电价（运行成本 = 日运行成本 × 基准电价）
)
println("最优运行成本绘图完成！共 $(num_oc) 张图")
