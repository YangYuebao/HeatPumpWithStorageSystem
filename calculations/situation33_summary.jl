# =============================================================================
# situation33 数值实验汇总脚本
# 读取 calculations/situation33 下所有成功算例，统计:
#   1. dp_precise（有初始化）与 dp_nostart（无初始化）的迭代次数、成本、收敛状态
#   2. 输出对比CSV + 终端汇总表 + 迭代过程绘图（可选）
# =============================================================================
using JSON3
using DataFrames, CSV
using Statistics
using Plots

situation = "situation33"
file_path0 = joinpath(pwd(), "calculations", situation)

storage_folders = filter(x -> isdir(joinpath(file_path0, x)), readdir(file_path0))

rows = []
for folder in storage_folders
    folder_path = joinpath(file_path0, folder)
    json_files = filter(x -> endswith(x, ".json"), readdir(folder_path))
    for json_file in json_files
        json_path = joinpath(folder_path, json_file)
        try
            data = JSON3.read(read(json_path, String))
            if data.status != "success"
                continue
            end

            hp = data.heatPumpServiceCoff
            sc = data.heatStorageCapacity
            mh = data.maxheatStorageInputHour

            # 有初始化组 (dp_precise)
            dp = haskey(data, :dp_precise) ? data.dp_precise : nothing
            # 无初始化组 (dp_nostart)
            dn = haskey(data, :dp_nostart) ? data.dp_nostart : nothing

            # 平均每次迭代最优解下降量（用 best_cost 差分，只统计 e>0 的改善轮）
            function avg_drop(hist)
                hist === nothing && return missing
                isempty(hist) && return missing
                b = [get(h, "best_cost", h.cost) for h in hist]  # 兼容旧JSON无best_cost
                n = length(b)
                n < 2 && return missing
                e = b[1:n-1] .- b[2:n]
                pos = e[e .> 0]
                isempty(pos) && return missing
                return mean(pos)
            end

            push!(rows, (
                heatPumpServiceCoff = hp,
                heatStorageCapacity = sc,
                maxheatStorageInputHour = mh,
                # 初始解成本（阶段1 DP初始解，用于与最优解对比）
                init_initial_cost = (haskey(data, :initial) && data.initial !== nothing) ? data.initial.cost : missing,
                # 有初始化
                init_converged = dp !== nothing ? dp.converged : missing,
                init_cost = dp !== nothing ? dp.cost : missing,
                init_n_iter = dp !== nothing ? length(dp.history) : missing,
                init_avg_drop = avg_drop(dp !== nothing ? dp.history : nothing),
                # 无初始化
                nostart_converged = dn !== nothing ? dn.converged : missing,
                nostart_cost = dn !== nothing ? dn.cost : missing,
                nostart_n_iter = dn !== nothing ? length(dn.history) : missing,
                nostart_avg_drop = avg_drop(dn !== nothing ? dn.history : nothing),
            ))
        catch e
            @warn "读取失败: $(json_file)" exception = e
        end
    end
end

if isempty(rows)
    println("没有找到有效的算例结果: ", file_path0)
    return
end

df = DataFrame(rows)
sort!(df, [:heatStorageCapacity, :heatPumpServiceCoff, :maxheatStorageInputHour])

# 迭代次数差值/比例
df.init_minus_nostart = df.init_n_iter .- df.nostart_n_iter
df.cost_diff = df.init_cost .- df.nostart_cost

# 保存CSV
csv_path = joinpath(file_path0, "experiment_summary.csv")
CSV.write(csv_path, df)
println("汇总CSV已保存: ", csv_path)

# ===== 终端汇总表 =====
println("\n" * "="^80)
println("situation33 数值实验汇总（有初始化 dp_precise vs 无初始化 dp_nostart）")
println("="^80)
println(DataFrames.select(df, Not([:init_converged, :nostart_converged])))

# ===== 统计 =====
println("\n----- 统计 -----")
valid = df[.!ismissing.(df.init_n_iter) .& .!ismissing.(df.nostart_n_iter), :]
if nrow(valid) > 0
    println("有效对比算例数: ", nrow(valid))
    println("有初始化平均迭代次数: ", round(mean(valid.init_n_iter), digits=2))
    println("无初始化平均迭代次数: ", round(mean(valid.nostart_n_iter), digits=2))
    println("平均减少迭代次数: ", round(mean(valid.init_minus_nostart), digits=2))
    println("平均减少比例: ", round(mean(-valid.init_minus_nostart ./ valid.nostart_n_iter) * 100, digits=1), "%")
    println("有初始化收敛率: ", round(sum(df.init_converged .== true) / count(!ismissing, df.init_converged) * 100, digits=1), "%")
    println("无初始化收敛率: ", round(sum(df.nostart_converged .== true) / count(!ismissing, df.nostart_converged) * 100, digits=1), "%")
    println("成本差异(有-无): 平均 ", round(mean(valid.cost_diff), digits=4))
end

# ===== 初始解价值分析三图 =====
# 横轴统一使用 CSV 排序后的算例编号（1..nrow，与 experiment_summary.csv 行号一致）
let
    out_dir = joinpath(file_path0, "experiment_plots")
    mkpath(out_dir)

    idx1 = Int[]; y1 = Float64[]   # 图1: 初始解成本 - 有初始化最优解成本
    idx2 = Int[]; y2 = Float64[]   # 图2: 无初始化迭代次数 - 有初始化迭代次数
    idx3 = Int[]; y3 = Float64[]   # 图3: 无初始化最优解成本 - 有初始化最优解成本

    for (i, r) in enumerate(eachrow(df))
        # 图1: 所有启用初始解的算例（需 initial.cost 与 dp_precise.cost 均存在）
        if !ismissing(r.init_initial_cost) && !ismissing(r.init_cost)
            push!(idx1, i)
            push!(y1, r.init_initial_cost - r.init_cost)
        end
        # 图2: 有效对比算例（两组迭代次数均存在）
        if !ismissing(r.init_n_iter) && !ismissing(r.nostart_n_iter)
            push!(idx2, i)
            push!(y2, r.nostart_n_iter - r.init_n_iter)
        end
        # 图3: 有效对比算例（两组最优解成本均存在）
        if !ismissing(r.init_cost) && !ismissing(r.nostart_cost)
            push!(idx3, i)
            push!(y3, r.nostart_cost - r.init_cost)
        end
    end

    # 图1: 初始解相对最优解的差距
    if !isempty(idx1)
        p1 = bar(idx1, y1,
            size = (900, 400),
            xlabel = "Case Index",
            ylabel = "Initial Cost - Optimal Cost (CNY)",
            title = "Initial Solution Gap to Optimal (With Init)",
            legend = false, color = :blue,
            right_margin = 18*Plots.mm,
            bottom_margin = 18*Plots.mm,
            left_margin = 5*Plots.mm,
            #top_margin = 18*Plots.mm,
        )
        savefig(p1, joinpath(out_dir, "gap_initial_vs_optimal.png"))
        println("\n图1 初始解-最优解差值已生成: ", length(idx1), " 个算例 → ", out_dir)
    else
        println("\n图1 跳过: 没有同时含 initial 与 dp_precise 成本的算例")
    end

    # 图2: 初始解减少的迭代次数
    if !isempty(idx2)
        p2 = bar(idx2, y2,
            size = (900, 400),
            xlabel = "Case Index",
            ylabel = "No-Init Iterations - With-Init Iterations",
            title = "Iteration Reduction by Initial Solution",
            legend = false, color = :green,
            right_margin = 18*Plots.mm,
            bottom_margin = 18*Plots.mm,
            left_margin = 5*Plots.mm,
            #top_margin = 18*Plots.mm,
        )
        savefig(p2, joinpath(out_dir, "iter_reduction.png"))
        println("图2 迭代次数减少已生成: ", length(idx2), " 个算例 → ", out_dir)
    else
        println("\n图2 跳过: 没有同时含两组迭代次数的算例")
    end

    # 图3: 初始解改善的最优性（避免过早陷入局部最优）
    if !isempty(idx3)
        p3 = bar(idx3, y3,
            size = (900, 400),
            xlabel = "Case Index",
            ylabel = "No-Init Optimal Cost - With-Init Optimal Cost (CNY)",
            title = "Optimality Improvement by Initial Solution",
            legend = false, color = :orange,
            right_margin = 18*Plots.mm,
            bottom_margin = 18*Plots.mm,
            left_margin = 5*Plots.mm,
            #top_margin = 18*Plots.mm,
        )
        savefig(p3, joinpath(out_dir, "optimality_improvement.png"))
        println("图3 最优性改善已生成: ", length(idx3), " 个算例 → ", out_dir)
    else
        println("\n图3 跳过: 没有同时含两组最优解成本的算例")
    end
end

# ===== 迭代过程绘图（可选） =====
# 对每个算例绘制两组成本下降曲线
let
    out_dir = joinpath(file_path0, "experiment_plots")
    mkpath(out_dir)

    # 对所有算例绘制成本下降对比图
    plot_count = 0
    for folder in storage_folders
        folder_path = joinpath(file_path0, folder)
        json_files = filter(x -> endswith(x, ".json"), readdir(folder_path))
        for json_file in json_files
            json_path = joinpath(folder_path, json_file)
            data = JSON3.read(read(json_path, String))
            if data.status != "success" || data.dp_precise === nothing || data.dp_nostart === nothing
                continue
            end
            hp = data.heatPumpServiceCoff
            sc = data.heatStorageCapacity
            mh = data.maxheatStorageInputHour

            # 取有history的记录
            hist_init = collect(data.dp_precise.history)
            hist_nostart = collect(data.dp_nostart.history)
            isempty(hist_init) && isempty(hist_nostart) && continue

            # 迭代成本曲线
            plt = plot(size=(900, 400), xlabel="Iteration", ylabel="Cost (CNY)",
                title="Cost Convergence: HP=$(hp), SC=$(sc), MH=$(mh)",
                legend=:topright)
            if !isempty(hist_init)
                plot!(plt, [h.cost for h in hist_init], label="With Init", lw=2, color=:blue)
            end
            if !isempty(hist_nostart)
                plot!(plt, [h.cost for h in hist_nostart], label="No Init", lw=2, color=:red, linestyle=:dash)
            end
            png = joinpath(out_dir, "convergence_hp$(hp)_sc$(sc)_mh$(mh).png")
            savefig(plt, png)
            plot_count += 1
        end
    end
    println("\n迭代过程图已生成: ", plot_count, " 张 → ", out_dir)
end

# ===== 最优解下降量差分对数图 =====
#= e[k] = best_cost[k] - best_cost[k+1]（只保留 e>0 的改善轮），log10 坐标
let
    out_dir = joinpath(file_path0, "experiment_plots")
    mkpath(out_dir)

    improvement_count = 0
    for folder in storage_folders
        folder_path = joinpath(file_path0, folder)
        json_files = filter(x -> endswith(x, ".json"), readdir(folder_path))
        for json_file in json_files
            json_path = joinpath(folder_path, json_file)
            data = JSON3.read(read(json_path, String))
            if data.status != "success" || data.dp_precise === nothing || data.dp_nostart === nothing
                continue
            end
            hp = data.heatPumpServiceCoff
            sc = data.heatStorageCapacity
            mh = data.maxheatStorageInputHour

            hist_init = collect(data.dp_precise.history)
            hist_nostart = collect(data.dp_nostart.history)
            isempty(hist_init) && isempty(hist_nostart) && continue

            # 最优解下降量: e[k] = b[k] - b[k+1]，只保留 e>0（log10 坐标要求正值）
            b_init  = [get(h, "best_cost", h.cost) for h in hist_init]
            b_nostart = [get(h, "best_cost", h.cost) for h in hist_nostart]
            e_init  = b_init[1:end-1]  .- b_init[2:end]
            e_nostart = b_nostart[1:end-1] .- b_nostart[2:end]
            pos_init  = findall(>(0), e_init)
            pos_nostart = findall(>(0), e_nostart)

            if isempty(pos_init) && isempty(pos_nostart)
                continue
            end

            plt = plot(size=(900, 400), xlabel="Iteration k (1st→2nd as k=1)",
                ylabel="Best cost drop (CNY, log10)",
                title="Improvement (log): HP=$(hp), SC=$(sc), MH=$(mh)",
                yscale=:log10, legend=:topright)
            if !isempty(pos_init)
                plot!(plt, pos_init, e_init[pos_init], label="With Init",
                    lw=2, color=:blue, marker=:circle)
            end
            if !isempty(pos_nostart)
                plot!(plt, pos_nostart, e_nostart[pos_nostart], label="No Init",
                    lw=2, color=:red, linestyle=:dash, marker=:circle)
            end
            png = joinpath(out_dir, "improvement_hp$(hp)_sc$(sc)_mh$(mh).png")
            savefig(plt, png)
            improvement_count += 1
        end
    end
    println("\n下降量对数图已生成: ", improvement_count, " 张 → ", out_dir)
end
=#

# ===== log-log 自回归散点图 =====
# e[k] = best_cost[k] - best_cost[k+1]（仅保留 e>0 的改善轮）
# 散点: x = log(e[k-1]), y = log(e[k])，检验 e 是否近似几何衰减
# 拟合: y = a*x + b（最小二乘）
# 分算例与合并分开画

"""
    线性拟合 y = a*x + b（最小二乘）

返回 (a, b, R²)，点数不足或方差为0时返回 (NaN, NaN, NaN)。
"""
function linfit(x, y)
    n = length(x)
    n < 2 && return (NaN, NaN, NaN)
    xbar = sum(x) / n
    ybar = sum(y) / n
    sxx = sum((x .- xbar) .^ 2)
    sxy = sum((x .- xbar) .* (y .- ybar))
    sxx <= 0 && return (NaN, NaN, NaN)
    a = sxy / sxx
    b = ybar - a * xbar
    ss_res = sum((y .- (a .* x .+ b)) .^ 2)
    ss_tot = sum((y .- ybar) .^ 2)
    r2 = ss_tot > 0 ? 1 - ss_res / ss_tot : NaN
    return (a, b, r2)
end

# 提取 (log(e[k-1]), log(e[k])) 点对的辅助函数
function loglog_pairs(hist)
    b = [get(h, "best_cost", h.cost) for h in hist]
    n = length(b)
    n < 3 && return (Float64[], Float64[])
    e = b[1:n-1] .- b[2:n]
    e = e[e .> 0]                     # 只保留正下降量
    length(e) < 2 && return (Float64[], Float64[])
    x = log.(e[1:end-1])              # log(e[k-1])
    y = log.(e[2:end])                # log(e[k])
    return (x, y)
end

# 分算例 + 合并收集
let
    out_dir = joinpath(file_path0, "experiment_plots")
    mkpath(out_dir)

    # 合并收集（所有算例）
    all_x_init = Float64[]
    all_y_init = Float64[]
    all_x_nostart = Float64[]
    all_y_nostart = Float64[]

    per_case_count = 0
    for folder in storage_folders
        folder_path = joinpath(file_path0, folder)
        json_files = filter(x -> endswith(x, ".json"), readdir(folder_path))
        for json_file in json_files
            json_path = joinpath(folder_path, json_file)
            data = JSON3.read(read(json_path, String))
            if data.status != "success" || data.dp_precise === nothing || data.dp_nostart === nothing
                continue
            end
            hp = data.heatPumpServiceCoff
            sc = data.heatStorageCapacity
            mh = data.maxheatStorageInputHour

            hist_init = collect(data.dp_precise.history)
            hist_nostart = collect(data.dp_nostart.history)

            # 分算例图
            x_init, y_init = loglog_pairs(hist_init)
            x_nostart, y_nostart = loglog_pairs(hist_nostart)
            if isempty(x_init) && isempty(x_nostart)
                continue
            end

            plt = plot(size=(900, 400),
                xlabel="log(e[k-1])", ylabel="log(e[k])",
                title="Log-log Autoregression: HP=$(hp), SC=$(sc), MH=$(mh)",
                legend=:topleft)
            if !isempty(x_init)
                scatter!(plt, x_init, y_init, label="With Init", color=:blue, ms=4)
                a, b, r2 = linfit(x_init, y_init)
                !isnan(a) && plot!(plt, [minimum(x_init), maximum(x_init)],
                    a .* [minimum(x_init), maximum(x_init)] .+ b,
                    label="Fit (a=$(round(a, digits=3)), R²=$(round(r2, digits=3)))",
                    color=:blue, lw=2)
            end
            if !isempty(x_nostart)
                scatter!(plt, x_nostart, y_nostart, label="No Init", color=:red, ms=4)
                a, b, r2 = linfit(x_nostart, y_nostart)
                !isnan(a) && plot!(plt, [minimum(x_nostart), maximum(x_nostart)],
                    a .* [minimum(x_nostart), maximum(x_nostart)] .+ b,
                    label="Fit (a=$(round(a, digits=3)), R²=$(round(r2, digits=3)))",
                    color=:red, lw=2, linestyle=:dash)
            end
            png = joinpath(out_dir, "loglog_hp$(hp)_sc$(sc)_mh$(mh).png")
            savefig(plt, png)
            per_case_count += 1

            # 收集到合并集合
            append!(all_x_init, x_init)
            append!(all_y_init, y_init)
            append!(all_x_nostart, x_nostart)
            append!(all_y_nostart, y_nostart)
        end
    end
    println("\nlog-log分算例图已生成: ", per_case_count, " 张 → ", out_dir)

    # 合并图（所有算例点集合在一张图）
    if !isempty(all_x_init) || !isempty(all_x_nostart)
        plt_all = plot(size=(900, 500),
            xlabel="log(e[k-1])", ylabel="log(e[k])",
            title="Log-log Autoregression (all cases combined)",
            legend=:topleft)
        if !isempty(all_x_init)
            scatter!(plt_all, all_x_init, all_y_init, label="With Init", color=:blue, ms=3, alpha=0.6)
            a, b, r2 = linfit(all_x_init, all_y_init)
            !isnan(a) && plot!(plt_all, [minimum(all_x_init), maximum(all_x_init)],
                a .* [minimum(all_x_init), maximum(all_x_init)] .+ b,
                label="Fit (a=$(round(a, digits=3)), b=$(round(b, digits=3)), R²=$(round(r2, digits=3)))",
                color=:blue, lw=2)
        end
        if !isempty(all_x_nostart)
            scatter!(plt_all, all_x_nostart, all_y_nostart, label="No Init", color=:red, ms=3, alpha=0.6)
            a, b, r2 = linfit(all_x_nostart, all_y_nostart)
            !isnan(a) && plot!(plt_all, [minimum(all_x_nostart), maximum(all_x_nostart)],
                a .* [minimum(all_x_nostart), maximum(all_x_nostart)] .+ b,
                label="Fit (a=$(round(a, digits=3)), b=$(round(b, digits=3)), R²=$(round(r2, digits=3)))",
                color=:red, lw=2, linestyle=:dash)
        end
        png_all = joinpath(out_dir, "loglog_all.png")
        savefig(plt_all, png_all)
        println("log-log合并图已生成: ", png_all)
    end
end

# ===== 层间 log-log 自回归散点图（合并版） =====
# 对每个算例按 dT_local 分层，每层取最小当轮 cost 得序列 c_min[1..m]（dT 从大到小）
# 层间下降量 e[k] = c_min[k] - c_min[k+1]（仅保留 e>0）
# 散点: x = log(e[k-1]), y = log(e[k])；所有算例点合并，检验层间是否几何衰减

"""
    每层最小当轮 cost 序列（按 dT_local 从大到小，即迭代顺序）

返回 (dT列表, 每层最小cost序列)。若无数据返回空。
"""
function layer_min_costs(hist)
    isempty(hist) && return (Float64[], Float64[])
    groups = Dict{Float64, Vector{Float64}}()
    for rec in hist
        push!(get!(groups, rec.dT_local, Float64[]), rec.cost)
    end
    dts = sort!(collect(keys(groups)); rev=true)   # dT 从大到小 = 网格细化顺序
    c_min = [minimum(groups[dt]) for dt in dts]
    return (dts, c_min)
end

# 层间点对: x = log(e[k-1]), y = log(e[k])，e = 层间最小cost差分（>0）
function layer_loglog_pairs(hist)
    _, c = layer_min_costs(hist)
    m = length(c)
    m < 3 && return (Float64[], Float64[])
    e = c[1:m-1] .- c[2:m]          # 层间下降量
    e = e[e .> 0]                   # 只保留正下降
    length(e) < 2 && return (Float64[], Float64[])
    return (log.(e[1:end-1]), log.(e[2:end]))
end

let
    out_dir = joinpath(file_path0, "experiment_plots")
    mkpath(out_dir)

    # 合并收集所有算例的层间点对
    layer_x_init = Float64[]
    layer_y_init = Float64[]
    layer_x_nostart = Float64[]
    layer_y_nostart = Float64[]

    for folder in storage_folders
        folder_path = joinpath(file_path0, folder)
        json_files = filter(x -> endswith(x, ".json"), readdir(folder_path))
        for json_file in json_files
            json_path = joinpath(folder_path, json_file)
            data = JSON3.read(read(json_path, String))
            if data.status != "success" || data.dp_precise === nothing || data.dp_nostart === nothing
                continue
            end
            hist_init = collect(data.dp_precise.history)
            hist_nostart = collect(data.dp_nostart.history)

            x_init, y_init = layer_loglog_pairs(hist_init)
            x_nostart, y_nostart = layer_loglog_pairs(hist_nostart)
            append!(layer_x_init, x_init)
            append!(layer_y_init, y_init)
            append!(layer_x_nostart, x_nostart)
            append!(layer_y_nostart, y_nostart)
        end
    end

    if !isempty(layer_x_init) || !isempty(layer_x_nostart)
        plt_layer = plot(size=(900, 500),
            xlabel="log(e_layer[k-1])", ylabel="log(e_layer[k])",
            title="Layer-wise Log-log Autoregression (all cases, min cost per dT layer)",
            legend=:topleft)
        if !isempty(layer_x_init)
            scatter!(plt_layer, layer_x_init, layer_y_init, label="With Init", color=:blue, ms=5, alpha=0.7)
            a, b, r2 = linfit(layer_x_init, layer_y_init)
            !isnan(a) && plot!(plt_layer, [minimum(layer_x_init), maximum(layer_x_init)],
                a .* [minimum(layer_x_init), maximum(layer_x_init)] .+ b,
                label="Fit (a=$(round(a, digits=3)), b=$(round(b, digits=3)), R²=$(round(r2, digits=3)), n=$(length(layer_x_init)))",
                color=:blue, lw=2)
        end
        if !isempty(layer_x_nostart)
            scatter!(plt_layer, layer_x_nostart, layer_y_nostart, label="No Init", color=:red, ms=5, alpha=0.7)
            a, b, r2 = linfit(layer_x_nostart, layer_y_nostart)
            !isnan(a) && plot!(plt_layer, [minimum(layer_x_nostart), maximum(layer_x_nostart)],
                a .* [minimum(layer_x_nostart), maximum(layer_x_nostart)] .+ b,
                label="Fit (a=$(round(a, digits=3)), b=$(round(b, digits=3)), R²=$(round(r2, digits=3)), n=$(length(layer_x_nostart)))",
                color=:red, lw=2, linestyle=:dash)
        end
        png_layer = joinpath(out_dir, "loglog_layer_all.png")
        savefig(plt_layer, png_layer)
        println("\n层间log-log合并图已生成: ", png_layer)
        println("  有初始化层间点数: ", length(layer_x_init), ", 无初始化层间点数: ", length(layer_x_nostart))
    else
        println("\n层间log-log: 无有效点对（各算例层数不足），跳过")
    end
end
