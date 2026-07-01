"""
分段COP数据结构。
"""
struct COPPiecewiseData
    minTe::Real                  # 蒸发温度下限
    maxTe::Real                  # 蒸发温度上限
    minTc::Real                  # 冷凝温度下限
    maxTc::Real                  # 冷凝温度上限
    refrigerant::Refrigerant     # 工质
    maxCOP::Real                 # 最大COP
    eta_s::Real                  # 绝热效率
    dT::Real                     # 插值步长
    Te_points::Vector{Real}      # 蒸发温度节点
    Tc_points::Vector{Real}      # 冷凝温度节点
    COP_points::Matrix{Real}     # COP值矩阵
    Te_pieces::Int               # 蒸发温度段数
    Tc_pieces::Int               # 冷凝温度段数
end


