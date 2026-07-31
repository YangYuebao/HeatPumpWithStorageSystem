using Plots
using JSON3

function operation_result_plot(
        time_list,
        hourly_tariff_ori,
        result;
        w=0.45
    )
    items = [
        "Grid Peak",
        "Grid Shoulder",
        "Grid Valley",
        "HP Supply",
        "HS Supply",
        "HP Store",
        "Elec On"
    ]

    # 自定义时间段
    time_periods = [(time_list[i],time_list[i+1]) for i in 1:length(time_list)-1]

    peak_list = hourly_tariff_ori .> 1
    shoulder_list = hourly_tariff_ori .== 1.0
    valley_list = hourly_tariff_ori .< 1
    cost = result[1]
    Ts_list = result[2]
    P1_list = result[3]
    P2_list = result[4]
    P3_list = result[5]
    Pe_list = result[6]


    # 转换为0-1表示
    status = [peak_list shoulder_list valley_list]
    for list in [P1_list, P2_list, P3_list, Pe_list]
        status = hcat(status, list .> 0)
    end
    status=transpose(status)

    # 绘图
    #legend=:outertopright
    p = plot(legend=:none, size=(900, 300), ylim=(0.5, 7.5), xlim=(0, 24))

    # 为每个项目绘制矩形
    for i in 1:7  # 3个智能体
        y_pos = 7 - i + 1  # Agent A在顶部(y=3), C在底部(y=1)
        
        for j in 1:length(time_periods)
            start_time, end_time = time_periods[j]
            width = end_time - start_time
            
            # 设置颜色
            color = status[i, j] == 1 ? :green : :red
            
            # 绘制矩形: [x_start, x_end, y_low, y_high]
            plot!(p,[start_time, end_time, end_time, start_time, start_time], 
                [y_pos-w, y_pos-w, y_pos+w, y_pos+w, y_pos-w],
                fill=(0, color), linecolor=:black, linewidth=0.5,
                label=false)
        end
    end
    # 添加时间标签（每隔2小时显示）
    xticks!(p,0:2:24, string.(0:2:24, "h"))

    # 设置y轴
    yticks!(p,collect(1:7), reverse(items))
    #xlabel!("hour H")
    title!(p,"Operation chedule")


    p_t=plot(time_list,Ts_list,legend=:none, size=(900, 200), ylim=(115,225), xlim=(0, 24),title="Storage Temperature")
    xticks!(p_t,0:2:24, string.(0:2:24, "h"))
    ylabel!(p_t,"℃")

    plt_combined = plot(
        p_t, p,
        layout = grid(2, 1, heights=[0.33, 0.67]),
        size=(900, 800),
        link=:x
    )

    return plt_combined
end

# =============================================================================
# 批量绘制运行结果
# =============================================================================

"""
    batch_plot_operation_results(file_path0, t_list, tariff_list;
        result_key="operationResults", cost_field="objective",
        pe_l_field="P_el", pe_s_field="P_es", output_suffix="", output_dir="")

遍历指定目录下所有 storage_*/xxx.json 文件，读取运行结果并批量生成运行调度图。

# 参数
- `file_path0`: 算例根目录（如 `calculations/situation30/`）
- `t_list`: 时间边界列表
- `tariff_list`: 分时电价列表（与时间段一一对应）
- `result_key`: JSON中运行结果所在的键名，默认 "operationResults"。
  对于situation32的MILP结果用 "milp"，DP精确解用 "dp_precise"
- `cost_field`: 目标函数值字段名，默认 "objective"
- `pe_l_field`: 电加热补热功率字段名，默认 "P_el"
- `pe_s_field`: 电加热蓄热功率字段名，默认 "P_es"
- `output_suffix`: 输出PNG文件名后缀（用于区分MILP/DP结果）
- `output_dir`: 输出子文件夹名，空字符串表示直接放在算例子目录

# 返回值
- `(total_plotted, total_skipped, total_failed)`: 成功绘制数、跳过数、失败数
"""
function batch_plot_operation_results(
    file_path0::String,
    t_list::Vector{Float64},
    tariff_list::Vector{Float64};
    result_key::String="operationResults",
    cost_field::String="objective",
    pe_l_field::String="P_el",
    pe_s_field::String="P_es",
    output_suffix::String="",
    output_dir::String=""
)
    storage_folders = filter(x -> isdir(joinpath(file_path0, x)), readdir(file_path0))

    total_plotted = 0
    total_skipped = 0
    total_failed_plot = 0

    for folder in storage_folders
        folder_path = joinpath(file_path0, folder)

        # 根据output_dir决定PNG保存位置
        if !isempty(output_dir)
            target_dir = joinpath(folder_path, output_dir)
            mkpath(target_dir)
        else
            target_dir = folder_path
        end

        json_files = filter(x -> endswith(x, ".json"), readdir(folder_path))

        for json_file in json_files
            json_path = joinpath(folder_path, json_file)

            try
                data = JSON3.read(read(json_path, String))

                # 检查状态
                if data.status == "infeasible"
                    total_skipped += 1
                    continue
                end

                # 检查运行结果是否存在
                if !haskey(data, Symbol(result_key)) || data[Symbol(result_key)] === nothing
                    total_skipped += 1
                    continue
                end

                op_result = data[Symbol(result_key)]

                # 构造 result 数组：[cost, Ts, P1, P2, P3, Pe_total, unused]
                result = [
                    op_result[Symbol(cost_field)],
                    op_result.Ts,
                    op_result.P1,
                    op_result.P2,
                    op_result.P3,
                    op_result[Symbol(pe_l_field)] + op_result[Symbol(pe_s_field)],
                    zeros(length(op_result.P1))
                ]

                # 生成PNG文件名：原JSON名去掉.json，加上后缀
                base_name = replace(json_file, ".json" => "")
                png_name = base_name * output_suffix * ".png"
                png_path = joinpath(target_dir, png_name)

                plt = operation_result_plot(
                    t_list,
                    tariff_list,
                    result;
                    w=0.45
                )
                savefig(plt, png_path)

                total_plotted += 1

            catch e
                @warn "绘图失败: $(json_file)" exception=e
                total_failed_plot += 1
            end
        end
    end

    return total_plotted, total_skipped, total_failed_plot
end

# =============================================================================
# 经济性分析绘图
# =============================================================================

"""
    plot_economic_analysis(file_path0; params_key, econ_key, capital_field, total_pw_field, plot_modes, use_subfolders)

读取指定目录下所有成功算例的JSON结果，生成经济性分析图（堆叠柱状图 + 总现值对比图）。

# 参数
- `file_path0`: 算例根目录（如 `calculations/situation30/`）
- `params_key`: 案例参数所在的JSON键名。situation30用 `"caseParameters"`，situation32用 `""`（顶层）
- `econ_key`: 经济数据所在的JSON键名。situation30用 `"economicResults"`，situation32用 `"milp"`
- `capital_field`: 初投资字段名，situation30默认 `"capitalCost"`，situation32用 `"C_initial"`
- `total_pw_field`: 总现值字段名，situation30默认 `"totalPresentWorth"`，situation32用 `"C_total"`
- `plot_modes`: 要绘制的模式，默认 `[:all]` 表示全部绘制。可选值：
  - `:hp_vs_capacity`  - 堆叠图: 固定(蓄热,时长), x=热泵容量
  - `:hour_comparison`  - 对比图: 固定蓄热, 不同时长, x=热泵容量
  - `:storage_vs_capacity` - 堆叠图: 固定(热泵,时长), x=蓄热容量
  - `:storage_comparison`  - 对比图: 固定热泵, 不同时长, x=蓄热容量
  - `:hour_vs_capacity`    - 堆叠图: 固定(蓄热,热泵), x=储满时长
  - `:hp_comparison`       - 对比图: 固定蓄热, 不同热泵, x=储满时长
- `use_subfolders`: 是否使用子文件夹组织图片，默认 `true`

# 输出
根据 `use_subfolders` 决定输出位置：
- `true`: 输出到 `vary_heatpump/`, `vary_storage/`, `vary_hour/` 三个子文件夹
- `false`: 直接输出到 `file_path0`

# 返回值
- `(num_stacked, num_comparison)`: 绘制的堆叠图数量和对比图数量
"""
function plot_economic_analysis(
    file_path0::String;
    params_key::String="caseParameters",
    econ_key::String="economicResults",
    capital_field::String="capitalCost",
    total_pw_field::String="totalPresentWorth",
    plot_modes::Vector{Symbol}=[:all],
    use_subfolders::Bool=true
)
    all_data = []

    storage_folders = filter(x -> isdir(joinpath(file_path0, x)), readdir(file_path0))

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

                # 提取案例参数：params_key为空时从顶层读，否则从子对象读
                if isempty(params_key)
                    hp = data.heatPumpServiceCoff
                    sc = data.heatStorageCapacity
                    mh = data.maxheatStorageInputHour
                else
                    params = data[Symbol(params_key)]
                    hp = params.heatPumpServiceCoff
                    sc = params.heatStorageCapacity
                    mh = params.maxHeatStorageInputHour
                end

                # 提取经济数据
                if isempty(econ_key)
                    econ = data
                else
                    if !haskey(data, Symbol(econ_key)) || data[Symbol(econ_key)] === nothing
                        continue
                    end
                    econ = data[Symbol(econ_key)]
                end

                cap = econ[Symbol(capital_field)]
                tpw = econ[Symbol(total_pw_field)]

                if cap === nothing || tpw === nothing
                    continue
                end

                push!(all_data, (
                    heatPumpCapacity = hp,
                    storageCapacity = sc,
                    maxInputHour = mh,
                    capitalCost = cap,
                    totalPresentWorth = tpw
                ))

            catch e
                @warn "读取JSON失败: $(json_file)" exception=e
            end
        end
    end

    if isempty(all_data)
        @warn "没有找到有效的经济数据，跳过经济性分析绘图"
        return 0, 0
    end

    println("\n共收集 $(length(all_data)) 个成功算例")

    # 处理 :all 模式
    all_mode_list = [:hp_vs_capacity, :hour_comparison,
                     :storage_vs_capacity, :storage_comparison,
                     :hour_vs_capacity, :hp_comparison]
    if :all in plot_modes
        active_modes = all_mode_list
    else
        active_modes = plot_modes
    end

    max_total_PW = maximum([item.totalPresentWorth for item in all_data])

    # 创建子文件夹
    if use_subfolders
        vary_hp_dir = joinpath(file_path0, "vary_heatpump")
        vary_sc_dir = joinpath(file_path0, "vary_storage")
        vary_hour_dir = joinpath(file_path0, "vary_hour")
        mkpath(vary_hp_dir)
        mkpath(vary_sc_dir)
        mkpath(vary_hour_dir)
    else
        vary_hp_dir = file_path0
        vary_sc_dir = file_path0
        vary_hour_dir = file_path0
    end

    num_stacked = 0
    num_comparison = 0

    # 颜色调色板
    color_palette = [:blue, :red, :green, :orange, :purple, :brown, :pink, :gray]

    # ================================================================
    # 第一组: 变热泵容量 (vary_heatpump)
    # ================================================================

    # --- 堆叠图: 固定(蓄热,时长), x=热泵容量 ---
    if :hp_vs_capacity in active_modes
        println("\n[变热泵容量] 生成堆叠图...")

        grouped_by_config = Dict{Tuple{Float64, Float64}, Vector{Any}}()
        for item in all_data
            key = (item.storageCapacity, item.maxInputHour)
            if !haskey(grouped_by_config, key)
                grouped_by_config[key] = []
            end
            push!(grouped_by_config[key], item)
        end

        for ((storageCap, maxHour), items) in sort(collect(grouped_by_config), by=x -> (x[1][1], x[1][2]))
            sort!(items, by=x -> x.heatPumpCapacity)

            heatPumpCapacities = [item.heatPumpCapacity for item in items]
            capitalCosts = [item.capitalCost for item in items]
            totalPWs = [item.totalPresentWorth for item in items]

            plt = plot(
                size=(800, 600),
                grid=true,
                xlabel="Heat Pump Capacity",
                ylabel="Cost (CNY)",
                title="Economic Analysis\nStorage = $(round(storageCap, digits=1)) kWh, Max Input Hour = $(round(maxHour, digits=1)) h",
                legend=:topleft,
                ylim=(0, max_total_PW * 1.05)
            )

            plot!(plt, heatPumpCapacities, capitalCosts,
                label="Capital Cost", lw=2, color=:blue,
                fillrange=zeros(length(capitalCosts)), fillalpha=0.5, fillcolor=:blue)
            plot!(plt, heatPumpCapacities, totalPWs,
                label="Operating Cost PW", lw=2, color=:red,
                fillrange=capitalCosts, fillalpha=0.5, fillcolor=:red)

            storage_str = string(round(storageCap, digits=1))
            hour_str = string(round(maxHour, digits=1))
            png_path = joinpath(vary_hp_dir, "stacked_storage_$(storage_str)_hour_$(hour_str).png")
            savefig(plt, png_path)
            num_stacked += 1
        end
        println("  已生成 $(num_stacked) 张堆叠图")
    end

    # --- 对比图: 固定蓄热, 不同时长, x=热泵容量 ---
    if :hour_comparison in active_modes
        println("\n[变热泵容量] 生成时长对比图...")

        grouped_by_storage = Dict{Float64, Vector{Any}}()
        for item in all_data
            sc = item.storageCapacity
            if !haskey(grouped_by_storage, sc)
                grouped_by_storage[sc] = []
            end
            push!(grouped_by_storage[sc], item)
        end

        count = 0
        for (storageCap, items) in sort(collect(grouped_by_storage), by=x -> x[1])
            hour_groups = Dict{Float64, Vector{Any}}()
            for item in items
                mh = item.maxInputHour
                if !haskey(hour_groups, mh)
                    hour_groups[mh] = []
                end
                push!(hour_groups[mh], item)
            end

            plt = plot(
                size=(800, 600),
                grid=true,
                xlabel="Heat Pump Capacity",
                ylabel="Total Present Worth (CNY)",
                title="Total Present Worth Comparison\nStorage Capacity = $(round(storageCap, digits=1)) kWh",
                legend=:topleft,
                ylim=(0, max_total_PW * 1.05)
            )

            color_idx = 1
            for (hour, hour_items) in sort(collect(hour_groups), by=x -> x[1])
                sort!(hour_items, by=x -> x.heatPumpCapacity)
                hp_caps = [item.heatPumpCapacity for item in hour_items]
                tpws = [item.totalPresentWorth for item in hour_items]

                plot!(plt, hp_caps, tpws,
                    label="$(round(hour, digits=1)) h", lw=2,
                    marker=:circle, markersize=4,
                    color=color_palette[(color_idx - 1) % length(color_palette) + 1])
                color_idx += 1
            end

            storage_str = string(round(storageCap, digits=1))
            png_path = joinpath(vary_hp_dir, "comparison_storage_$(storage_str).png")
            savefig(plt, png_path)
            count += 1
        end
        num_comparison += count
        println("  已生成 $(count) 张对比图")
    end

    # ================================================================
    # 第二组: 变蓄热容量 (vary_storage)
    # ================================================================

    # --- 堆叠图: 固定(热泵,时长), x=蓄热容量 ---
    if :storage_vs_capacity in active_modes
        println("\n[变蓄热容量] 生成堆叠图...")

        grouped_by_hp_hour = Dict{Tuple{Float64, Float64}, Vector{Any}}()
        for item in all_data
            key = (item.heatPumpCapacity, item.maxInputHour)
            if !haskey(grouped_by_hp_hour, key)
                grouped_by_hp_hour[key] = []
            end
            push!(grouped_by_hp_hour[key], item)
        end

        count = 0
        for ((hpCap, maxHour), items) in sort(collect(grouped_by_hp_hour), by=x -> (x[1][1], x[1][2]))
            sort!(items, by=x -> x.storageCapacity)

            storageCaps = [item.storageCapacity for item in items]
            capitalCosts = [item.capitalCost for item in items]
            totalPWs = [item.totalPresentWorth for item in items]

            plt = plot(
                size=(800, 600),
                grid=true,
                xlabel="Heat Storage Capacity (kWh)",
                ylabel="Total Present Worth (CNY)",
                title="Total Present Worth vs Storage Capacity\nHP = $(round(hpCap, digits=1)), Max Input Hour = $(round(maxHour, digits=1)) h",
                legend=:topleft,
                ylim=(0, max_total_PW * 1.05)
            )

            plot!(plt, storageCaps, capitalCosts,
                label="Capital Cost", lw=2, color=:blue,
                fillrange=zeros(length(capitalCosts)), fillalpha=0.5, fillcolor=:blue)
            plot!(plt, storageCaps, totalPWs,
                label="Total Present Worth", lw=2, color=:red,
                fillrange=capitalCosts, fillalpha=0.5, fillcolor=:red)

            hp_str = string(round(hpCap, digits=1))
            hour_str = string(round(maxHour, digits=1))
            png_path = joinpath(vary_sc_dir, "stacked_hp_$(hp_str)_hour_$(hour_str).png")
            savefig(plt, png_path)
            count += 1
        end
        num_stacked += count
        println("  已生成 $(count) 张堆叠图")
    end

    # --- 对比图: 固定热泵, 不同时长, x=蓄热容量 ---
    if :storage_comparison in active_modes
        println("\n[变蓄热容量] 生成时长对比图...")

        grouped_by_hp = Dict{Float64, Vector{Any}}()
        for item in all_data
            hp = item.heatPumpCapacity
            if !haskey(grouped_by_hp, hp)
                grouped_by_hp[hp] = []
            end
            push!(grouped_by_hp[hp], item)
        end

        count = 0
        for (hpCap, items) in sort(collect(grouped_by_hp), by=x -> x[1])
            hour_groups = Dict{Float64, Vector{Any}}()
            for item in items
                mh = item.maxInputHour
                if !haskey(hour_groups, mh)
                    hour_groups[mh] = []
                end
                push!(hour_groups[mh], item)
            end

            plt = plot(
                size=(800, 600),
                grid=true,
                xlabel="Heat Storage Capacity (kWh)",
                ylabel="Total Present Worth (CNY)",
                title="Total Present Worth Comparison\nHeat Pump Capacity = $(round(hpCap, digits=1))",
                legend=:topleft,
                ylim=(0, max_total_PW * 1.05)
            )

            color_idx = 1
            for (hour, hour_items) in sort(collect(hour_groups), by=x -> x[1])
                sort!(hour_items, by=x -> x.storageCapacity)
                sc_caps = [item.storageCapacity for item in hour_items]
                tpws = [item.totalPresentWorth for item in hour_items]

                plot!(plt, sc_caps, tpws,
                    label="$(round(hour, digits=1)) h", lw=2,
                    marker=:circle, markersize=4,
                    color=color_palette[(color_idx - 1) % length(color_palette) + 1])
                color_idx += 1
            end

            hp_str = string(round(hpCap, digits=1))
            png_path = joinpath(vary_sc_dir, "comparison_hp_$(hp_str).png")
            savefig(plt, png_path)
            count += 1
        end
        num_comparison += count
        println("  已生成 $(count) 张对比图")
    end

    # ================================================================
    # 第三组: 变储满时长 (vary_hour)
    # ================================================================

    # --- 堆叠图: 固定(蓄热,热泵), x=储满时长 ---
    if :hour_vs_capacity in active_modes
        println("\n[变储满时长] 生成堆叠图...")

        grouped_by_sc_hp = Dict{Tuple{Float64, Float64}, Vector{Any}}()
        for item in all_data
            key = (item.storageCapacity, item.heatPumpCapacity)
            if !haskey(grouped_by_sc_hp, key)
                grouped_by_sc_hp[key] = []
            end
            push!(grouped_by_sc_hp[key], item)
        end

        count = 0
        for ((storageCap, hpCap), items) in sort(collect(grouped_by_sc_hp), by=x -> (x[1][1], x[1][2]))
            sort!(items, by=x -> x.maxInputHour)

            maxHours = [item.maxInputHour for item in items]
            capitalCosts = [item.capitalCost for item in items]
            totalPWs = [item.totalPresentWorth for item in items]

            plt = plot(
                size=(800, 600),
                grid=true,
                xlabel="Max Heat Storage Input Hour (h)",
                ylabel="Total Present Worth (CNY)",
                title="Total Present Worth vs Max Input Hour\nStorage = $(round(storageCap, digits=1)) kWh, HP = $(round(hpCap, digits=1))",
                legend=:topleft,
                ylim=(0, max_total_PW * 1.05)
            )

            plot!(plt, maxHours, capitalCosts,
                label="Capital Cost", lw=2, color=:blue,
                fillrange=zeros(length(capitalCosts)), fillalpha=0.5, fillcolor=:blue)
            plot!(plt, maxHours, totalPWs,
                label="Total Present Worth", lw=2, color=:red,
                fillrange=capitalCosts, fillalpha=0.5, fillcolor=:red)

            hp_str = string(round(hpCap, digits=1))
            storage_str = string(round(storageCap, digits=1))
            png_path = joinpath(vary_hour_dir, "stacked_storage_$(storage_str)_hp_$(hp_str).png")
            savefig(plt, png_path)
            count += 1
        end
        num_stacked += count
        println("  已生成 $(count) 张堆叠图")
    end

    # --- 对比图: 固定蓄热, 不同热泵, x=储满时长 ---
    if :hp_comparison in active_modes
        println("\n[变储满时长] 生成热泵对比图...")

        grouped_by_storage = Dict{Float64, Vector{Any}}()
        for item in all_data
            sc = item.storageCapacity
            if !haskey(grouped_by_storage, sc)
                grouped_by_storage[sc] = []
            end
            push!(grouped_by_storage[sc], item)
        end

        count = 0
        for (storageCap, items) in sort(collect(grouped_by_storage), by=x -> x[1])
            hp_groups = Dict{Float64, Vector{Any}}()
            for item in items
                hp = item.heatPumpCapacity
                if !haskey(hp_groups, hp)
                    hp_groups[hp] = []
                end
                push!(hp_groups[hp], item)
            end

            plt = plot(
                size=(800, 600),
                grid=true,
                xlabel="Max Heat Storage Input Hour (h)",
                ylabel="Total Present Worth (CNY)",
                title="Total Present Worth Comparison\nStorage Capacity = $(round(storageCap, digits=1)) kWh",
                legend=:topleft,
                ylim=(0, max_total_PW * 1.05)
            )

            color_idx = 1
            for (hp, hp_items) in sort(collect(hp_groups), by=x -> x[1])
                sort!(hp_items, by=x -> x.maxInputHour)
                hours = [item.maxInputHour for item in hp_items]
                tpws = [item.totalPresentWorth for item in hp_items]

                plot!(plt, hours, tpws,
                    label="HP $(round(hp, digits=1))", lw=2,
                    marker=:circle, markersize=4,
                    color=color_palette[(color_idx - 1) % length(color_palette) + 1])
                color_idx += 1
            end

            storage_str = string(round(storageCap, digits=1))
            png_path = joinpath(vary_hour_dir, "comparison_storage_$(storage_str).png")
            savefig(plt, png_path)
            count += 1
        end
        num_comparison += count
        println("  已生成 $(count) 张对比图")
    end

    println("\n经济性分析绘图完成！堆叠图 $(num_stacked) 张，对比图 $(num_comparison) 张")
    return num_stacked, num_comparison
end

#=
plt = operation_result_plot(
    0:dt:24,
    hourly_tariff_ori,
    result;
    w=0.45
)

savefig(plt, "plots/plt_2.png")
=#
