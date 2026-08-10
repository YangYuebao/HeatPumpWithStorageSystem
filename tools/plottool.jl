using Plots
using JSON3

function operation_result_plot(
        time_list,
        hourly_tariff_ori,
        result;
        w = 0.45,
        pe_l_list = nothing,   # 电锅炉补热功率（堆叠柱状图分层用），nothing时回退为总Pe
        pe_s_list = nothing,   # 电锅炉储热功率（堆叠柱状图分层用），nothing时回退为0
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
    Pe_list = result[6]    # 合并电锅炉功率（补热+储热），用于状态时序图"电加热是否启用"判断

    # 堆叠柱状图用的拆分电锅炉功率（Pe_l补热 / Pe_s储热）
    # 未提供时回退为 Pe_l=总Pe、Pe_s=0，保证旧调用点不报错
    Pe_l = pe_l_list === nothing ? Pe_list : pe_l_list
    Pe_s = pe_s_list === nothing ? zeros(length(Pe_list)) : pe_s_list


    # 转换为0-1表示
    status = [peak_list shoulder_list valley_list]
    for list in [P1_list, P2_list, P3_list, Pe_list]
        status = hcat(status, list .> 0)
    end
    status=transpose(status)

    # 绘图
    #legend=:outertopright
    p = plot(legend=:none, size=(900, 200), ylim=(0.5, 7.5), xlim=(0, 24))

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
    #title!(p,"Operation chedule")


    p_t = mode_energy_stacked_plot(
        Float64.(collect(time_list)),
        Float64.(collect(Ts_list)),
        Float64.(collect(P1_list)),
        Float64.(collect(P2_list)),
        Float64.(collect(P3_list)),
        Float64.(collect(Pe_l)),
        Float64.(collect(Pe_s)),
        labels = ["P1 (HP direct)" "P2 (HS supply)" "P3 (HP store)" "Pe_l (EH heat)" "Pe_s (EH store)"],
        palette = [:steelblue, :tomato, :orange, :greenyellow, :purple],
        size = (900, 400),
        ylabel = "Power (kW)",
    )
    title!(p_t,"Operation chedule")
    #=
    p_t=plot(time_list,Ts_list,legend=:none, size=(900, 200), ylim=(115,225), xlim=(0, 24),title="Storage Temperature")
    xticks!(p_t,0:2:24, string.(0:2:24, "h"))
    ylabel!(p_t,"℃")
    =#
    plt_combined = plot(
        p_t, p,
        layout = grid(2, 1, heights=[0.7, 0.3]),
        size=(900, 800),
        link=:x,
        #title=["Storage Temperature" "Operation chedule"]
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
                    w = 0.45,
                    pe_l_list = op_result[Symbol(pe_l_field)],
                    pe_s_list = op_result[Symbol(pe_s_field)],
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
    plot_economic_analysis(file_path0; params_key, econ_key, capital_field, total_pw_field, plot_modes, use_subfolders, base_price)

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
- `base_price`: 目标基准电价（¥/kWh），默认 `1.0`。
  算例按基准电价 `economicResults.basePrice` 计算，总现值 = 初投资 + 运行成本现值 × (base_price / basePrice)。
  用于研究不同基准电价下的最优配置：算例按基准电价1计算，绘图时传实际电价即可。

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
    use_subfolders::Bool=true,
    base_price::Float64=1.0
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

                # 按基准电价折算总现值：
                # 总现值 = 初投资 + 运行成本现值，运行成本现值 ∝ 电价。
                # JSON中记录计算时的基准电价 basePrice，此处折算到目标基准电价 base_price。
                json_base_price = haskey(econ, :basePrice) ? econ.basePrice : 1.0
                operatingPV = haskey(econ, :operatingPV) ? econ.operatingPV : (tpw - cap)
                tpw_adj = cap + operatingPV * (base_price / json_base_price)

                push!(all_data, (
                    heatPumpCapacity = hp,
                    storageCapacity = sc,
                    maxInputHour = mh,
                    capitalCost = cap,
                    totalPresentWorth = tpw_adj
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
    all_mode_list = [:hp_vs_capacity, :hour_comparison, :storage_comparison_hp,
                     :storage_vs_capacity, :storage_comparison, :hp_comparison_sc,
                     :hour_vs_capacity, :hp_comparison, :storage_comparison_hr]
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

    # --- 对比图: 固定时长, 不同蓄热, x=热泵容量 ---
    if :storage_comparison_hp in active_modes
        println("\n[变热泵容量] 生成蓄热对比图...")

        grouped_by_hour = Dict{Float64, Vector{Any}}()
        for item in all_data
            mh = item.maxInputHour
            if !haskey(grouped_by_hour, mh)
                grouped_by_hour[mh] = []
            end
            push!(grouped_by_hour[mh], item)
        end

        count = 0
        for (maxHour, items) in sort(collect(grouped_by_hour), by=x -> x[1])
            storage_groups = Dict{Float64, Vector{Any}}()
            for item in items
                sc = item.storageCapacity
                if !haskey(storage_groups, sc)
                    storage_groups[sc] = []
                end
                push!(storage_groups[sc], item)
            end

            plt = plot(
                size=(800, 600),
                grid=true,
                xlabel="Heat Pump Capacity",
                ylabel="Total Present Worth (CNY)",
                title="Total Present Worth Comparison\nMax Input Hour = $(round(maxHour, digits=1)) h",
                legend=:topleft,
                ylim=(0, max_total_PW * 1.05)
            )

            color_idx = 1
            for (storage, storage_items) in sort(collect(storage_groups), by=x -> x[1])
                sort!(storage_items, by=x -> x.heatPumpCapacity)
                hp_caps = [item.heatPumpCapacity for item in storage_items]
                tpws = [item.totalPresentWorth for item in storage_items]

                plot!(plt, hp_caps, tpws,
                    label="Storage $(round(storage, digits=1)) kWh", lw=2,
                    marker=:circle, markersize=4,
                    color=color_palette[(color_idx - 1) % length(color_palette) + 1])
                color_idx += 1
            end

            hour_str = string(round(maxHour, digits=1))
            png_path = joinpath(vary_hp_dir, "comparison_hour_$(hour_str).png")
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

    # --- 对比图: 固定时长, 不同热泵, x=蓄热容量 ---
    if :hp_comparison_sc in active_modes
        println("\n[变蓄热容量] 生成热泵对比图...")

        grouped_by_hour = Dict{Float64, Vector{Any}}()
        for item in all_data
            mh = item.maxInputHour
            if !haskey(grouped_by_hour, mh)
                grouped_by_hour[mh] = []
            end
            push!(grouped_by_hour[mh], item)
        end

        count = 0
        for (maxHour, items) in sort(collect(grouped_by_hour), by=x -> x[1])
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
                xlabel="Heat Storage Capacity (kWh)",
                ylabel="Total Present Worth (CNY)",
                title="Total Present Worth Comparison\nMax Input Hour = $(round(maxHour, digits=1)) h",
                legend=:topleft,
                ylim=(0, max_total_PW * 1.05)
            )

            color_idx = 1
            for (hp, hp_items) in sort(collect(hp_groups), by=x -> x[1])
                sort!(hp_items, by=x -> x.storageCapacity)
                sc_caps = [item.storageCapacity for item in hp_items]
                tpws = [item.totalPresentWorth for item in hp_items]

                plot!(plt, sc_caps, tpws,
                    label="HP $(round(hp, digits=1))", lw=2,
                    marker=:circle, markersize=4,
                    color=color_palette[(color_idx - 1) % length(color_palette) + 1])
                color_idx += 1
            end

            hour_str = string(round(maxHour, digits=1))
            png_path = joinpath(vary_sc_dir, "comparison_hour_$(hour_str).png")
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

    # --- 对比图: 固定热泵, 不同蓄热, x=储满时长 ---
    if :storage_comparison_hr in active_modes
        println("\n[变储满时长] 生成蓄热对比图...")

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
            storage_groups = Dict{Float64, Vector{Any}}()
            for item in items
                sc = item.storageCapacity
                if !haskey(storage_groups, sc)
                    storage_groups[sc] = []
                end
                push!(storage_groups[sc], item)
            end

            plt = plot(
                size=(800, 600),
                grid=true,
                xlabel="Max Heat Storage Input Hour (h)",
                ylabel="Total Present Worth (CNY)",
                title="Total Present Worth Comparison\nHeat Pump Capacity = $(round(hpCap, digits=1))",
                legend=:topleft,
                ylim=(0, max_total_PW * 1.05)
            )

            color_idx = 1
            for (storage, storage_items) in sort(collect(storage_groups), by=x -> x[1])
                sort!(storage_items, by=x -> x.maxInputHour)
                hours = [item.maxInputHour for item in storage_items]
                tpws = [item.totalPresentWorth for item in storage_items]

                plot!(plt, hours, tpws,
                    label="Storage $(round(storage, digits=1)) kWh", lw=2,
                    marker=:circle, markersize=4,
                    color=color_palette[(color_idx - 1) % length(color_palette) + 1])
                color_idx += 1
            end

            hp_str = string(round(hpCap, digits=1))
            png_path = joinpath(vary_hour_dir, "comparison_hp_$(hp_str).png")
            savefig(plt, png_path)
            count += 1
        end
        num_comparison += count
        println("  已生成 $(count) 张对比图")
    end

    println("\n经济性分析绘图完成！堆叠图 $(num_stacked) 张，对比图 $(num_comparison) 张")
    return num_stacked, num_comparison
end

# ================================================================
# 最小总现值包络分析绘图
# ================================================================
"""
    plot_optimal_envelope(
        file_path0::String;
        params_key::String="caseParameters",
        econ_key::String="economicResults",
        capital_field::String="capitalCost",
        total_pw_field::String="totalPresentWorth",
        output_dir::String="optimal_envelope"
    )

绘制"最小总现值包络分析"三类图。对电锅炉容量（maxInputHour）维度取最小值，
展示其余两个设计参数与最小总现值的关系：

1. `optimal_hp_vs_tpw.png` : 横轴=热泵容量, 纵轴=最小总现值,
   每条曲线=固定蓄热容量（对每个热泵容量取所有电锅炉容量中的最小总现值）
2. `optimal_sc_vs_tpw.png` : 横轴=蓄热容量, 纵轴=最小总现值,
   每条曲线=固定热泵容量（对每个蓄热容量取所有电锅炉容量中的最小总现值）
3. `optimal_hp_sc_heatmap.png` : 二维云图, 横轴=热泵容量, 纵轴=蓄热容量,
   颜色=该 (热泵, 蓄热) 组合下所有电锅炉容量中的最小总现值

# 参数
- `file_path0`: 结果根目录，遍历其下的 `storage_*` 子文件夹中的 JSON
- `params_key`: 参数所在JSON键名，空字符串表示从顶层读取
- `econ_key`: 经济数据所在JSON键名
- `capital_field` / `total_pw_field`: 初投资 / 总现值字段名
- `output_dir`: 输出子文件夹名（相对于 file_path0）

# 返回值
- `Int`: 生成的图数量（3）
"""
function plot_optimal_envelope(
    file_path0::String;
    params_key::String="caseParameters",
    econ_key::String="economicResults",
    capital_field::String="capitalCost",
    total_pw_field::String="totalPresentWorth",
    output_dir::String="optimal_envelope",
    base_price::Float64=1.0
)
    # ---- 数据收集（与 plot_economic_analysis 相同的读取逻辑） ----
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

                tpw = econ[Symbol(total_pw_field)]

                if tpw === nothing
                    continue
                end

                # 按基准电价折算总现值（与 plot_economic_analysis 相同逻辑）
                cap = haskey(econ, Symbol(capital_field)) ? econ[Symbol(capital_field)] : 0.0
                json_base_price = haskey(econ, :basePrice) ? econ.basePrice : 1.0
                operatingPV = haskey(econ, :operatingPV) ? econ.operatingPV : (tpw - cap)
                tpw_adj = cap + operatingPV * (base_price / json_base_price)

                push!(all_data, (
                    heatPumpCapacity = hp,
                    storageCapacity = sc,
                    maxInputHour = mh,
                    capitalCost = cap,
                    totalPresentWorth = tpw_adj
                ))

            catch e
                @warn "读取JSON失败: $(json_file)" exception=e
            end
        end
    end

    if isempty(all_data)
        @warn "没有找到有效的经济数据，跳过包络分析绘图"
        return 0
    end

    println("\n共收集 $(length(all_data)) 个成功算例")

    # ---- 降维: 对每个 (热泵, 蓄热) 组合取所有电锅炉容量中的最小总现值 ----
    # 键: (heatPumpCapacity, storageCapacity), 值: 最小 totalPresentWorth
    min_map = Dict{Tuple{Float64, Float64}, Float64}()
    for item in all_data
        key = (item.heatPumpCapacity, item.storageCapacity)
        if !haskey(min_map, key) || item.totalPresentWorth < min_map[key]
            min_map[key] = item.totalPresentWorth
        end
    end

    println("降维后 (热泵, 蓄热) 组合数: $(length(min_map))")

    # 创建输出文件夹
    out_path = joinpath(file_path0, output_dir)
    mkpath(out_path)

    # 颜色调色板
    color_palette = [:blue, :red, :green, :orange, :purple, :brown, :pink, :gray]

    num_plots = 0

    # ================================================================
    # 图1: 横轴=热泵容量, 每条曲线=固定蓄热容量
    # ================================================================
    println("\n[包络分析] 生成 热泵容量-最小总现值 曲线图...")

    # 分组: sc => Dict(hp => min_tpw)
    grouped_by_sc = Dict{Float64, Dict{Float64, Float64}}()
    for ((hp, sc), tpw) in min_map
        if !haskey(grouped_by_sc, sc)
            grouped_by_sc[sc] = Dict{Float64, Float64}()
        end
        grouped_by_sc[sc][hp] = tpw
    end

    max_tpw = maximum(values(min_map))
    min_tpw = minimum(values(min_map))

    plt1 = plot(
        size=(800, 600),
        grid=true,
        xlabel="Heat Pump Capacity",
        ylabel="Min Total Present Worth (CNY)",
        title="Min Total Present Worth vs Heat Pump Capacity\n(over all elec boiler capacities)",
        legend=:topleft,
        ylim=(min_tpw*0.95, max_tpw * 1.05)
    )

    color_idx = 1
    for (sc, hp_map) in sort(collect(grouped_by_sc), by=x -> x[1])
        hp_vec = sort(collect(keys(hp_map)))
        tpw_vec = [hp_map[hp] for hp in hp_vec]

        plot!(plt1, hp_vec, tpw_vec,
            label="Storage $(round(sc, digits=1)) kWh", lw=2,
            marker=:circle, markersize=4,
            color=color_palette[(color_idx - 1) % length(color_palette) + 1])
        color_idx += 1
    end

    png1 = joinpath(out_path, "optimal_hp_vs_tpw.png")
    savefig(plt1, png1)
    num_plots += 1
    println("  已生成: $(png1)")

    # ================================================================
    # 图2: 横轴=蓄热容量, 每条曲线=固定热泵容量
    # ================================================================
    println("\n[包络分析] 生成 蓄热容量-最小总现值 曲线图...")

    # 分组: hp => Dict(sc => min_tpw)
    grouped_by_hp = Dict{Float64, Dict{Float64, Float64}}()
    for ((hp, sc), tpw) in min_map
        if !haskey(grouped_by_hp, hp)
            grouped_by_hp[hp] = Dict{Float64, Float64}()
        end
        grouped_by_hp[hp][sc] = tpw
    end

    plt2 = plot(
        size=(800, 600),
        grid=true,
        xlabel="Heat Storage Capacity (kWh)",
        ylabel="Min Total Present Worth (CNY)",
        title="Min Total Present Worth vs Heat Storage Capacity\n(over all elec boiler capacities)",
        legend=:topleft,
        ylim=(min_tpw*0.95, max_tpw * 1.05)
    )

    color_idx = 1
    for (hp, sc_map) in sort(collect(grouped_by_hp), by=x -> x[1])
        sc_vec = sort(collect(keys(sc_map)))
        tpw_vec = [sc_map[sc] for sc in sc_vec]

        plot!(plt2, sc_vec, tpw_vec,
            label="HP $(round(hp, digits=1))", lw=2,
            marker=:circle, markersize=4,
            color=color_palette[(color_idx - 1) % length(color_palette) + 1])
        color_idx += 1
    end

    png2 = joinpath(out_path, "optimal_sc_vs_tpw.png")
    savefig(plt2, png2)
    num_plots += 1
    println("  已生成: $(png2)")

    # ================================================================
    # 图3: 二维云图, 横轴=热泵容量, 纵轴=蓄热容量, 颜色=最小总现值
    # ================================================================
    println("\n[包络分析] 生成 热泵-蓄热 最小总现值云图...")

    hp_vec = sort(unique([k[1] for k in keys(min_map)]))
    sc_vec = sort(unique([k[2] for k in keys(min_map)]))

    # heatmap(x, y, Z) 要求 Z 维度为 (length(y), length(x))
    Z = fill(NaN, length(sc_vec), length(hp_vec))
    for ((hp, sc), tpw) in min_map
        i = findfirst(==(sc), sc_vec)
        j = findfirst(==(hp), hp_vec)
        Z[i, j] = tpw
    end

    plt3 = heatmap(
        hp_vec, sc_vec, Z;
        size=(800, 600),
        xlabel="Heat Pump Capacity",
        ylabel="Heat Storage Capacity (kWh)",
        title="Min Total Present Worth (CNY)\n(over all elec boiler capacities)",
        color=:viridis,
        colorbar_title="TPW (CNY)",
    )

    png3 = joinpath(out_path, "optimal_hp_sc_heatmap.png")
    savefig(plt3, png3)
    num_plots += 1
    println("  已生成: $(png3)")

    # ================================================================
    # 图4: 3D曲面图, x=热泵容量, y=蓄热容量, z=最小总现值
    # ================================================================
    println("\n[包络分析] 生成 热泵-蓄热 最小总现值3D曲面图...")

    if length(hp_vec) < 2 || length(sc_vec) < 2
        # GR 的 gridit 需要至少 2x2 的网格，点数不足时跳过
        println("  点数不足（hp=$(length(hp_vec)), sc=$(length(sc_vec))），跳过3D曲面图")
    else
        # surface(x, y, Z) 要求 Z 维度为 (length(y), length(x))，与 heatmap 一致，可直接复用 Z
        plt4 = surface(
            hp_vec, sc_vec, Z;
            size=(800, 600),
            camera=(120, 30),              # 沿z轴顺时针转90°（默认(30,30)）
            xlabel="Heat Pump Capacity",
            ylabel="Heat Storage Capacity (kWh)",
            zlabel="Min TPW (CNY)",
            title="Min Total Present Worth (CNY)\n(over all elec boiler capacities)",
            color=:viridis,
            colorbar_title="TPW (CNY)",
        )

        png4 = joinpath(out_path, "optimal_hp_sc_3d.png")
        savefig(plt4, png4)
        num_plots += 1
        println("  已生成: $(png4)")
    end

    println("\n包络分析绘图完成！共 $(num_plots) 张图")
    return num_plots
end

# =============================================================================
# 各时段能耗模式构成堆叠柱状图
# =============================================================================

"""
    mode_energy_stacked_plot(time_list, Ts_list, P1_list, P2_list, P3_list, Pe_l_list, Pe_s_list;
        labels, palette, size, ylabel)

绘制各时段能耗的模式构成堆叠柱状图（5层：P1热泵供热 / P2蓄热供热 / P3热泵储热 / Pe_l电锅炉补热 / Pe_s电锅炉储热），
并叠加蓄热温度曲线（右侧y轴）。

# 参数
- `time_list`: 时间边界列表 [nt+1]
- `Ts_list`: 蓄热温度轨线 [nt+1]，绘制在右侧y轴
- `P1_list`: 模式1（热泵直接供热）功率 [nt]
- `P2_list`: 模式2（蓄热供热）功率 [nt]
- `P3_list`: 模式3（热泵储热）功率 [nt]
- `Pe_l_list`: 电锅炉补热功率 [nt]
- `Pe_s_list`: 电锅炉储热功率 [nt]
- `labels`: 图例标签（默认5个模式）
- `palette`: 各模式颜色
- `size`: 图幅
- `ylabel`: 左y轴标签

# 返回值
- Plots.jl 绘图对象
"""
function mode_energy_stacked_plot(
    time_list::Vector{Float64},
    Ts_list::Vector{Float64},
    P1_list::Vector{Float64},
    P2_list::Vector{Float64},
    P3_list::Vector{Float64},
    Pe_l_list::Vector{Float64},
    Pe_s_list::Vector{Float64};
    labels = ["P1 (HP direct)" "P2 (HS supply)" "P3 (HP store)" "Pe_l (EH heat)" "Pe_s (EH store)"],
    palette = [:steelblue, :tomato, :orange, :greenyellow, :purple],
    size = (900, 600),
    ylabel = "Power (kW)",
)
    n_seg = length(P1_list)

    # 5 层堆叠数据矩阵（累积值，配合 stacked=false，避免 Plots 二次累加）
    Pe_sv = Pe_s_list
    Pe_lv = Pe_sv + Pe_l_list
    P3v = Pe_lv + P3_list
    P2v = P3v + P2_list
    P1v = P2v + P1_list

    # x 坐标取各时间段中心点，柱宽取各时段实际长度（不等宽柱）
    x_centers = [(time_list[i] + time_list[i+1]) / 2 for i in 1:n_seg]
    widths = [time_list[i+1] - time_list[i] for i in 1:n_seg]

    data = hcat(P1v, P2v, P3v, Pe_lv, Pe_sv)

    p = bar(
        x_centers, data,
        stacked = false,
        bar_width = widths,          # 柱宽按实际时段长度
        labels = labels,
        palette = palette,
        legend = :topright,          # 图例放图内（放图外会导致右轴刻度被GR挤出画布）
        right_margin = 18Plots.mm,   # 为右轴刻度/标签预留空间，避免被裁剪
        left_margin = 5Plots.mm,     # 为左轴刻度/标签预留空间，避免被裁剪
        size = size,
        xlabel = "Time (h)",
        ylabel = ylabel,
        xlim = (0, time_list[end]),
    )
    xticks!(p, 0:2:round(Int, time_list[end]), string.(0:2:round(Int, time_list[end]), "h"))

    # 添加右侧y轴绘制蓄热温度曲线
    # xlabel="" 清除twin图共享x轴的重复标签，避免与主图"Time (h)"重影
    pt = twinx(p)
    plot!(pt, time_list, Ts_list, label = "Ts", ylabel = "Ts ℃",
        ylims = (50.0, 250.0), yticks = 50:50:250,
        color = :red, linestyle = :solid,
        legend = :none,              # 关闭twin图的图例，避免重复
        right_margin = 18Plots.mm,   # 与主图保持一致，右轴才不被裁剪
        left_margin = 5Plots.mm,     # 与主图保持一致，左轴才不被裁剪
        xlim = (0, time_list[end]),
        xlabel = "")

    # 合并图例: 在主图添加NaN副本系列（不画线，只贡献图例条目 "Ts"）
    plot!(p, time_list, fill(NaN, length(time_list)),
        label = "Ts", color = :red, linestyle = :solid, legend = :topright)

    return p
end

"""
    batch_plot_mode_energy_stacked(file_path0, t_list;
        result_key, pe_l_field, pe_s_field, output_suffix, output_dir)

批量绘制各时段能耗模式构成堆叠柱状图，输出 PNG 到各算例子目录。

# 参数
- `file_path0`: 算例根目录（如 `calculations/situation32/`）
- `t_list`: 时间边界列表
- `result_key`: JSON中运行结果所在的键名，如 "milp" / "dp_precise" / "initial"
- `pe_l_field`: 电锅炉补热功率字段名，MILP/初始解用 "P_el"，DP精确解用 "Pe_l"
- `pe_s_field`: 电锅炉储热功率字段名，MILP/初始解用 "P_es"，DP精确解用 "Pe_s"
- `output_suffix`: 输出PNG文件名后缀
- `output_dir`: 输出子文件夹名，空字符串表示直接放在算例子目录

# 返回值
- `(total_plotted, total_skipped, total_failed)`: 成功绘制数、跳过数、失败数
"""
function batch_plot_mode_energy_stacked(
    file_path0::String,
    t_list::Vector{Float64};
    result_key::String = "operationResults",
    pe_l_field::String = "P_el",
    pe_s_field::String = "P_es",
    output_suffix::String = "_energy",
    output_dir::String = "energy_stacked"
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

                # 检查所需的6个字段是否齐全（5个功率 + Ts温度轨线）
                required_fields = [:P1, :P2, :P3, Symbol(pe_l_field), Symbol(pe_s_field), :Ts]
                if !all(haskey(op_result, f) for f in required_fields)
                    total_skipped += 1
                    continue
                end

                plt = mode_energy_stacked_plot(
                    t_list,
                    collect(op_result.Ts) .|> Float64,
                    collect(op_result.P1) .|> Float64,
                    collect(op_result.P2) .|> Float64,
                    collect(op_result.P3) .|> Float64,
                    collect(op_result[Symbol(pe_l_field)]) .|> Float64,
                    collect(op_result[Symbol(pe_s_field)]) .|> Float64,
                )

                # 生成PNG文件名：原JSON名去掉.json，加上后缀
                base_name = replace(json_file, ".json" => "")
                png_name = base_name * output_suffix * ".png"
                png_path = joinpath(target_dir, png_name)

                savefig(plt, png_path)

                total_plotted += 1

            catch e
                @warn "堆叠能耗图绘制失败: $(json_file)" exception=e
                total_failed_plot += 1
            end
        end
    end

    return total_plotted, total_skipped, total_failed_plot
end

# =============================================================================
# 基准电价对最优配置影响分析绘图
# =============================================================================

"""
    plot_optimal_capacity_vs_price(
        file_path0::String;
        params_key::String="caseParameters",
        econ_key::String="economicResults",
        capital_field::String="capitalCost",
        total_pw_field::String="totalPresentWorth",
        base_prices::AbstractVector{<:Real}=0.6:0.05:1.0,
        output_dir::String="optimal_capacity_vs_price"
    )

绘制"不同基准电价下最优配置"单张双y轴图。总现值不读 JSON 中的 totalPresentWorth，
而是按 `tpw_adj = capitalCost + operatingPV × (bp / jsonBasePrice)` 对每个基准电价实时计算
（operatingPV 为计算基准电价 jsonBasePrice 下的运行成本现值，与基准电价成正比），
对每个基准电价求使总现值最小的全局最优配置：

- 横轴: 基准电价
- 左y轴: 最优热泵容量（蓝色）
- 右y轴: 最优蓄热容量（红色）

输出: `optimal_capacity_vs_price.png`

# 参数
- `file_path0`: 结果根目录，遍历其下的 `storage_*` 子文件夹中的 JSON
- `params_key` / `econ_key` / `capital_field`: 同 plot_optimal_envelope
- `operating_pv_field` / `base_price_field`: 运行成本现值 / 计算基准电价 字段名
- `base_prices`: 基准电价范围，按此实时缩放运行成本现值折算总现值
- `output_dir`: 输出子文件夹名（相对于 file_path0）

# 返回值
- `Int`: 生成的图数量（1）
"""
function plot_optimal_capacity_vs_price(
    file_path0::String;
    params_key::String = "caseParameters",
    econ_key::String = "economicResults",
    capital_field::String = "capitalCost",
    operating_pv_field::String = "operatingPV",
    base_price_field::String = "basePrice",
    base_prices::AbstractVector{<:Real} = 0.6:0.05:1.0,
    output_dir::String = "optimal_capacity_vs_price",
)
    # ---- 数据收集（与 plot_optimal_envelope 相同的读取逻辑） ----
    # 元组: (heatPumpCapacity, storageCapacity, maxInputHour, capitalCost, operatingPV, jsonBasePrice)
    all_data = Tuple{Float64, Float64, Float64, Float64, Float64, Float64}[]

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

                # 运行成本现值 operatingPV（按计算基准电价 basePrice 计算，与基准电价成正比）
                if !haskey(econ, Symbol(operating_pv_field)) || !haskey(econ, Symbol(base_price_field))
                    continue
                end
                operatingPV = econ[Symbol(operating_pv_field)]
                json_base_price = econ[Symbol(base_price_field)]

                cap = haskey(econ, Symbol(capital_field)) ? econ[Symbol(capital_field)] : 0.0

                push!(all_data, (
                    Float64(hp), Float64(sc), Float64(mh),
                    Float64(cap), Float64(operatingPV), Float64(json_base_price),
                ))

            catch e
                @warn "读取JSON失败: $(json_file)" exception=e
            end
        end
    end

    if isempty(all_data)
        @warn "没有找到有效的经济数据，跳过基准电价影响分析绘图"
        return 0
    end

    println("\n共收集 $(length(all_data)) 个成功算例")

    prices = sort(collect(Float64.(base_prices)))

    # 折算总现值: tpw_adj = capitalCost + operatingPV × (bp / jsonBasePrice)
    tpw_adj = (d, bp) -> d[4] + d[5] * (bp / d[6])

    # 对每个基准电价求全局最优配置（遍历所有 (hp, sc, mh) 组合，电锅炉容量维度隐式取最小）
    hp_opt = Float64[]   # 每个电价的最优热泵容量
    sc_opt = Float64[]   # 每个电价的最优蓄热容量
    for bp in prices
        best_hp = NaN
        best_sc = NaN
        best_cost = Inf
        for d in all_data
            c = tpw_adj(d, bp)
            if c < best_cost
                best_cost = c
                best_hp = d[1]
                best_sc = d[2]
            end
        end
        push!(hp_opt, best_hp)
        push!(sc_opt, best_sc)
    end

    # 创建输出文件夹
    out_path = joinpath(file_path0, output_dir)
    mkpath(out_path)

    num_plots = 0

    # 单张双y轴图：横轴=基准电价，左y轴=最优热泵容量，右y轴=最优蓄热容量
    # 边距参考 mode_energy_stacked_plot 的 twinx 经验：为右轴刻度/标签预留空间，避免被GR裁剪
    println("\n[基准电价影响] 生成 基准电价-最优配置 双y轴图...")

    plt = plot(
        prices, hp_opt,
        size = (800, 600),
        grid = true,
        xlabel = "Base Electricity Price",
        ylabel = "Optimal Heat Pump Capacity",
        title = "Optimal Configuration vs Base Electricity Price",
        legend = :topright,
        label = "Optimal Heat Pump Capacity",
        color = :blue,
        marker = :circle, markersize = 4, linewidth = 2,
        right_margin = 18Plots.mm,   # 为右轴刻度/标签预留空间，避免被裁剪
        left_margin = 5Plots.mm,     # 为左轴刻度/标签预留空间，避免被裁剪
    )

    # 右y轴绘制最优蓄热容量
    pt = twinx(plt)
    plot!(pt, prices, sc_opt,
        label = "Optimal Storage Capacity",
        ylabel = "Optimal Storage Capacity",
        color = :red,
        marker = :circle, markersize = 4, linewidth = 2,
        legend = :none,              # 关闭twin图的图例，避免重复
        right_margin = 18Plots.mm,   # 与主图保持一致，右轴才不被裁剪
        left_margin = 5Plots.mm,     # 与主图保持一致，左轴才不被裁剪
    )

    # 合并图例: 在主图添加NaN副本系列（不画线，只贡献图例条目）
    plot!(plt, prices, fill(NaN, length(prices)),
        label = "Optimal Storage Capacity", color = :red,
        marker = :circle, markersize = 4, linewidth = 2, legend = :topright)

    png = joinpath(out_path, "optimal_capacity_vs_price.png")
    savefig(plt, png)
    num_plots += 1
    println("  已生成: $(png)")

    println("\n基准电价影响分析绘图完成！共 $(num_plots) 张图")
    return num_plots
end

# =============================================================================
# 最优运行成本绘图
# =============================================================================

"""
    plot_optimal_operating_cost(
        file_path0::String;
        params_key::String="caseParameters",
        econ_key::String="economicResults",
        cost_field::String="operatingPV",
        base_price::Float64=1.0,
        output_dir::String="optimal_operating_cost"
    )

绘制"最优运行成本"曲线图。对每个 (热泵容量, 蓄热容量) 组合，取所有电锅炉容量
（maxInputHour）中最小的运行成本，按固定热泵容量分组绘制：

- 横轴: 蓄热容量
- 纵轴: 最小运行成本
- 每条曲线: 固定热泵容量

运行成本按目标基准电价缩放：`cost_adj = cost × base_price`
（cost 取每日运行成本 dailyOperationCost，乘以目标基准电价得该电价下的日运行成本）。

# 参数
- `file_path0`: 结果根目录，遍历其下的 `storage_*` 子文件夹中的 JSON
- `params_key` / `econ_key`: 同 plot_optimal_envelope
- `cost_field`: 运行成本字段名，默认 "dailyOperationCost"（每日运行成本）
- `base_price`: 目标基准电价，运行成本 = cost × base_price
- `output_dir`: 输出子文件夹名（相对于 file_path0）

# 返回值
- `Int`: 生成的图数量（1）
"""
function plot_optimal_operating_cost(
    file_path0::String;
    params_key::String = "caseParameters",
    econ_key::String = "economicResults",
    cost_field::String = "dailyOperationCost",
    base_price::Float64 = 1.0,
    output_dir::String = "optimal_operating_cost",
)
    # ---- 数据收集 ----
    # 元组: (heatPumpCapacity, storageCapacity, maxInputHour, cost_adj)
    all_data = Tuple{Float64, Float64, Float64, Float64}[]

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

                if isempty(econ_key)
                    econ = data
                else
                    if !haskey(data, Symbol(econ_key)) || data[Symbol(econ_key)] === nothing
                        continue
                    end
                    econ = data[Symbol(econ_key)]
                end

                if !haskey(econ, Symbol(cost_field))
                    continue
                end
                cost = econ[Symbol(cost_field)]

                # 运行成本 = 每日运行成本 × 目标基准电价
                cost_adj = Float64(cost) * base_price

                push!(all_data, (Float64(hp), Float64(sc), Float64(mh), cost_adj))

            catch e
                @warn "读取JSON失败: $(json_file)" exception = e
            end
        end
    end

    if isempty(all_data)
        @warn "没有找到有效的经济数据，跳过最优运行成本绘图"
        return 0
    end

    println("\n共收集 $(length(all_data)) 个成功算例")

    # 对每个 (hp, sc) 组合取所有电锅炉容量中最小的运行成本
    min_map = Dict{Tuple{Float64, Float64}, Float64}()
    for d in all_data
        key = (d[1], d[2])
        if !haskey(min_map, key) || d[4] < min_map[key]
            min_map[key] = d[4]
        end
    end

    sc_vals = sort(unique([k[2] for k in keys(min_map)]))
    hp_vals = sort(unique([k[1] for k in keys(min_map)]))

    out_path = joinpath(file_path0, output_dir)
    mkpath(out_path)

    color_palette = [:blue, :red, :green, :orange, :purple, :brown, :pink, :gray, :black]

    # 图: 横轴=蓄热容量, 每条曲线=固定热泵容量, 纵轴=最小运行成本
    println("\n[最优运行成本] 生成 蓄热容量-最小运行成本 图...")

    plt = plot(
        size = (800, 600),
        grid = true,
        xlabel = "Storage Capacity",
        ylabel = "Min Operating Cost (CNY/day)",
        title = "Optimal Operating Cost vs Storage Capacity",
        legend = :topright,
    )
    for (idx, hp) in enumerate(hp_vals)
        xs = Float64[]
        ys = Float64[]
        for sc in sc_vals
            if haskey(min_map, (hp, sc))
                push!(xs, sc)
                push!(ys, min_map[(hp, sc)])
            end
        end
        plot!(plt, xs, ys,
            label = string("HP = ", hp),
            color = color_palette[mod1(idx, length(color_palette))],
            marker = :circle, markersize = 4, linewidth = 2,
        )
    end

    png = joinpath(out_path, "optimal_operating_cost.png")
    savefig(plt, png)
    println("  已生成: $(png)")

    println("\n最优运行成本绘图完成！共 1 张图")
    return 1
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
