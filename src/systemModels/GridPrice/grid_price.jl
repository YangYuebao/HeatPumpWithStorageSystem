struct GridPrice
    name::String
    price::Float64
    time_gap::Float64
    n::Int
end

JiangSu_hourly_tariff_ori = ones(48)
p = 1.7
pp = p * 1.2
v = 0.35
vv = v * 0.8
JiangSu_hourly_tariff_ori[1:12] .*= v
JiangSu_hourly_tariff_ori[23:26] .*= v
JiangSu_hourly_tariff_ori[29:30] .*= pp
JiangSu_hourly_tariff_ori[31:39] .*= p
JiangSu_hourly_tariff_ori[40:43] .*= pp

JiangSu_grid_price = GridPrice(
    "Jiangsu",
    JiangSu_hourly_tariff_ori,
    0.5,
    48
)

pp=2.05
p = 1.85
v=0.4
vv=0.2
ZheJiang_hourly_tariff_ori = ones(24)
ZheJiang_hourly_tariff_ori[1:7] *= v
ZheJiang_hourly_tariff_ori[12:14] *=v
ZheJiang_hourly_tariff_ori[17:18] *=p
ZheJiang_hourly_tariff_ori[19:22] *=pp
ZheJiang_hourly_tariff_ori[23] *= p
ZheJiang_grid_price = GridPrice(
    "ZheJiang",
    ZheJiang_hourly_tariff_ori,
    1.0,
    24
)