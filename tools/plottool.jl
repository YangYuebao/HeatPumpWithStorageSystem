using Plots

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

#=
plt = operation_result_plot(
    0:dt:24,
    hourly_tariff_ori,
    result;
    w=0.45
)

savefig(plt, "plots/plt_2.png")
=#
