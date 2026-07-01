"""
参数文档生成工具
用于自动生成热泵与蓄热系统设计参数的 Markdown 文档
"""

using DataFrames

"""
合并连续相同的时段
返回: (起始索引, 结束索引, 值) 的列表
"""
function mergeContinuousPeriods(values::Vector)
    if isempty(values)
        return Tuple{Int, Int, Any}[]
    end
    
    periods = Tuple{Int, Int, Any}[]
    start_idx = 1
    current_value = values[1]
    
    for i in 2:length(values)
        if values[i] != current_value
            push!(periods, (start_idx, i-1, current_value))
            start_idx = i
            current_value = values[i]
        end
    end
    push!(periods, (start_idx, length(values), current_value))
    
    return periods
end

"""
生成电价曲线表格
"""
function generateTariffTable(hourlyTariff::Vector)
    n = length(hourlyTariff)
    maxTariff = maximum(hourlyTariff)
    
    periods = mergeContinuousPeriods(hourlyTariff)
    
    df = DataFrame(
        时段序号 = String[],
        时间区间 = String[],
        电价系数 = Float64[],
        电价_元千瓦时 = Float64[]
    )
    
    for (start_idx, end_idx, value) in periods
        # 时段序号
        period_str = if start_idx == end_idx
            "$start_idx"
        else
            "$start_idx-$end_idx"
        end
        
        # 时间区间
        start_hour = start_idx - 1
        end_hour = end_idx
        time_str = "$(lpad(start_hour, 2, '0')):00-$(lpad(end_hour, 2, '0')):00"
        
        # 电价系数
        coefficient = round(value / maxTariff, digits=2)
        
        push!(df, (period_str, time_str, coefficient, value))
    end
    
    return df
end

"""
生成负荷曲线表格
"""
function generateLoadTable(heatConsumptionPower::Vector)
    periods = mergeContinuousPeriods(heatConsumptionPower)
    
    df = DataFrame(
        时段序号 = String[],
        时间区间 = String[],
        用热负荷_千瓦 = Float64[]
    )
    
    for (start_idx, end_idx, value) in periods
        period_str = if start_idx == end_idx
            "$start_idx"
        else
            "$start_idx-$end_idx"
        end
        
        start_hour = start_idx - 1
        end_hour = end_idx
        time_str = "$(lpad(start_hour, 2, '0')):00-$(lpad(end_hour, 2, '0')):00"
        
        push!(df, (period_str, time_str, value))
    end
    
    return df
end

"""
生成环境温度表格
"""
function generateTemperatureTable(Tair::Vector)
    periods = mergeContinuousPeriods(Tair)
    
    df = DataFrame(
        时段序号 = String[],
        时间区间 = String[],
        环境温度_°C = Float64[]
    )
    
    for (start_idx, end_idx, value) in periods
        period_str = if start_idx == end_idx
            "$start_idx"
        else
            "$start_idx-$end_idx"
        end
        
        start_hour = start_idx - 1
        end_hour = end_idx
        time_str = "$(lpad(start_hour, 2, '0')):00-$(lpad(end_hour, 2, '0')):00"
        
        push!(df, (period_str, time_str, value))
    end
    
    return df
end

"""
从 OneStorageFinanceParameters 提取经济性参数
"""
function extractOneStorageParams(
    financeParams::OneStorageFinanceParameters
)
    df = DataFrame(
        参数 = String[],
        单位 = String[],
        数值 = Float64[],
        备注 = String[]
    )
    
    push!(df, ("电锅炉单位功率成本", "元/kW", financeParams.Pelec_cost, ""))
    push!(df, ("蓄热单位小时成本", "元/h", financeParams.Storage_cost, ""))
    push!(df, ("项目年限", "年", Float64(financeParams.Life_years), ""))
    push!(df, ("年运行天数", "天", Float64(financeParams.annual_days), ""))
    push!(df, ("折现率", "%", financeParams.Discount_rate * 100, ""))
    
    return df
end

"""
从 FinanceParameters 提取经济性参数（包含热泵）
"""
function extractPressedWaterParams(
    financeParams::FinanceParameters
)
    df = DataFrame(
        参数 = String[],
        单位 = String[],
        数值 = Float64[],
        备注 = String[]
    )
    
    # 计算热泵综合成本
    heatPumpCost = financeParams.Plow_cost + financeParams.Phigh_cost
    
    push!(df, ("热泵单位供热功率成本", "元/kW", heatPumpCost, "低温热泵与水蒸气压缩机统一折算"))
    push!(df, ("电锅炉单位功率成本", "元/kW", financeParams.Pelec_cost, ""))
    push!(df, ("蓄热单位小时成本", "元/h", financeParams.Storage_cost, ""))
    push!(df, ("项目年限", "年", Float64(financeParams.Life_years), ""))
    push!(df, ("年运行天数", "天", Float64(financeParams.annual_days), ""))
    push!(df, ("折现率", "%", financeParams.Discount_rate * 100, ""))
    
    return df
end

"""
生成运行条件限制表格（热泵系统）
"""
function generateOperatingConditionsTable(
    designInput::DesignOptimizeInput
)
    df = DataFrame(
        参数 = String[],
        单位 = String[],
        数值 = Float64[]
    )
    
    push!(df, ("绝热效率", "-", designInput.eta_s))
    push!(df, ("低温热泵蒸发温度", "°C", designInput.TWaste))
    push!(df, ("水蒸气压缩机吸气温度", "°C", designInput.TCompressorIn))
    push!(df, ("水蒸气压缩机最大饱和温度", "°C", designInput.maxTcHigh))
    push!(df, ("传热温差", "°C", designInput.dT_EvaporationStandard))
    push!(df, ("蓄热最低温度", "°C", designInput.Tsmin))
    push!(df, ("蓄热最高温度", "°C", designInput.Tsmax))
    
    return df
end

"""
生成运行条件限制表格（简单蓄热系统）
"""
function generateSimpleOperatingConditionsTable(
    Tsmin::Float64,
    Tsmax::Float64
)
    df = DataFrame(
        参数 = String[],
        单位 = String[],
        数值 = Float64[]
    )
    
    push!(df, ("蓄热最低温度", "°C", Tsmin))
    push!(df, ("蓄热最高温度", "°C", Tsmax))
    
    return df
end

"""
生成 Markdown 表格字符串
"""
function dataframeToMarkdown(df::DataFrame)
    if nrow(df) == 0
        return ""
    end
    
    # 表头
    headers = names(df)
    header_line = "| " * join(headers, " | ") * " |"
    
    # 分隔线
    separator = "|" * join([":--------:" for _ in 1:ncol(df)], "|") * "|"
    
    # 数据行
    rows = String[]
    for row in eachrow(df)
        row_str = "| " * join([string(val) for val in row], " | ") * " |"
        push!(rows, row_str)
    end
    
    return header_line * "\n" * separator * "\n" * join(rows, "\n")
end

"""
确保目录存在，如果不存在则创建
"""
function ensureDirectoryExists(filepath::String)
    dir_path = dirname(filepath)
    if !isdir(dir_path)
        mkpath(dir_path)
        @info "创建目录: $dir_path"
    end
end

"""
主函数：生成设计参数文档（热泵系统版本）
"""
function generateDesignDoc(
    scenario_name::String,
    designInput::DesignOptimizeInput,
    designParameters::DesignOptimizeParameters,
    financeParams::FinanceParameters,
    output_path::String
)
    @info "正在生成 $scenario_name 的设计参数文档..."
    
    # 生成各部分表格
    tariff_df = generateTariffTable(designInput.hourlyTariff)
    load_df = generateLoadTable(designInput.heatConsumptionPower)
    temp_df = generateTemperatureTable(designInput.Tair)
    operating_df = generateOperatingConditionsTable(designInput)
    econ_df = extractPressedWaterParams(financeParams)
    
    # 构建 Markdown 文档
    md_content = "# 案例参数\n\n"
    
    # 电价曲线
    md_content *= "## 电价曲线：\n"
    md_content *= dataframeToMarkdown(tariff_df)
    md_content *= "\n\n"
    
    # 负荷曲线
    md_content *= "## 负荷曲线：\n\n"
    md_content *= dataframeToMarkdown(load_df)
    md_content *= "\n\n采用单位化负荷曲线，最大用热负荷为$(maximum(designInput.heatConsumptionPower))。\n\n\n"
    
    # 环境温度
    md_content *= "## 环境温度变化：\n\n"
    md_content *= dataframeToMarkdown(temp_df)
    md_content *= "\n\n"
    
    # 运行条件
    md_content *= "## 运行条件限制：\n\n"
    md_content *= dataframeToMarkdown(operating_df)
    md_content *= "\n\n"
    
    # 经济性参数
    md_content *= "## 经济性参数\n\n"
    md_content *= dataframeToMarkdown(econ_df)
    md_content *= "\n\n"
    
    # 确保目录存在
    ensureDirectoryExists(output_path)
    
    # 写入文件
    open(output_path, "w") do f
        write(f, md_content)
    end
    
    @info "设计参数文档已保存至: $output_path"
end

"""
主函数：生成设计参数文档（简单蓄热系统版本 - OneStorage）
"""
function generateDesignDoc(
    scenario_name::String,
    hourlyTariff::Vector,
    heatConsumptionPower::Vector,
    Tair::Vector,
    Tsmin::Float64,
    Tsmax::Float64,
    financeParams::OneStorageFinanceParameters,
    output_path::String
)
    @info "正在生成 $scenario_name 的设计参数文档..."
    
    # 生成各部分表格
    tariff_df = generateTariffTable(hourlyTariff)
    load_df = generateLoadTable(heatConsumptionPower)
    temp_df = generateTemperatureTable(Tair)
    operating_df = generateSimpleOperatingConditionsTable(Tsmin, Tsmax)
    econ_df = extractOneStorageParams(financeParams)
    
    # 构建 Markdown 文档
    md_content = "# 案例参数\n\n"
    
    # 电价曲线
    md_content *= "## 电价曲线：\n"
    md_content *= dataframeToMarkdown(tariff_df)
    md_content *= "\n\n"
    
    # 负荷曲线
    md_content *= "## 负荷曲线：\n\n"
    md_content *= dataframeToMarkdown(load_df)
    md_content *= "\n\n采用单位化负荷曲线，最大用热负荷为$(maximum(heatConsumptionPower))。\n\n\n"
    
    # 环境温度
    md_content *= "## 环境温度变化：\n\n"
    md_content *= dataframeToMarkdown(temp_df)
    md_content *= "\n\n"
    
    # 运行条件
    md_content *= "## 运行条件限制：\n\n"
    md_content *= dataframeToMarkdown(operating_df)
    md_content *= "\n\n"
    
    # 经济性参数
    md_content *= "## 经济性参数\n\n"
    md_content *= dataframeToMarkdown(econ_df)
    md_content *= "\n\n"
    
    # 确保目录存在
    ensureDirectoryExists(output_path)
    
    # 写入文件
    open(output_path, "w") do f
        write(f, md_content)
    end
    
    @info "设计参数文档已保存至: $output_path"
end
