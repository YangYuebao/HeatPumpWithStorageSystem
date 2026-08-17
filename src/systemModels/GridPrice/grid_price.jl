struct GridPrice
    name::String
    city::String
    season::String
    price::Vector{Float64}   # 分时电价向量（相对基准电价的倍数，[n]个时段）
    time_gap::Float64
    n::Int
end

# =============================================================================
# 江苏省工商业分时电价（苏发改价格发〔2025〕426号，2025-06-01起执行）
# 参考: https://fzggw.jiangsu.gov.cn/art/2025/4/30/art_284_11557001.html
# 浮动比例（以用户购电价格作为平段电价，此处采用单一制100千伏安及以上档）：
#   峰段上浮70%   → 峰   p  = 1.7
#   谷段下浮65%   → 谷   v  = 0.35
#   尖峰 = 峰段上浮20%     → 尖峰 pp = 1.7 × 1.2 = 2.04
#   深谷 = 谷段下浮20%     → 深谷 vv = 0.35 × 0.8 = 0.28
# 注: 尖峰（仅315千伏安及以上工业用电执行）因季节而异；
#     深谷时段为节假日10:00-14:00（315千伏安及以上），此处不单独建模
# 时段: 夏冬季 峰14:00-22:00，谷0:00-6:00、11:00-13:00，平6:00-11:00、13:00-14:00、22:00-24:00
#       春秋季 峰15:00-22:00，谷2:00-6:00、10:00-14:00，平6:00-10:00、14:00-15:00、22:00-次日2:00
# 半小时精度（48点，time_gap=0.5）：尖峰时段含19:30/21:30边界
# =============================================================================
p = 1.7     # 峰: 上浮70%
pp = 2.04   # 尖峰: 峰段上浮20% = 1.7 × 1.2
v = 0.35    # 谷: 下浮65%
vv = 0.28   # 深谷: 谷段下浮20% = 0.35 × 0.8

# 夏、冬两季（6-8月、12月-次年2月）基础时段
JiangSu_6_8_12_2_hourly_tariff_ori = ones(48)
JiangSu_6_8_12_2_hourly_tariff_ori[1:12] .*= v      # 谷 0:00-6:00
JiangSu_6_8_12_2_hourly_tariff_ori[23:26] .*= v     # 谷 11:00-13:00
JiangSu_6_8_12_2_hourly_tariff_ori[29:44] .*= p     # 峰 14:00-22:00

# 夏冬季 + 尖峰（7-8月: 14:00-15:00、19:30-21:30）
JiangSu_7_8_hourly_tariff_ori = copy(JiangSu_6_8_12_2_hourly_tariff_ori)
JiangSu_7_8_hourly_tariff_ori[29:30] .= pp          # 尖峰 14:00-15:00
JiangSu_7_8_hourly_tariff_ori[40:43] .= pp          # 尖峰 19:30-21:30
JiangSu_7_8_grid_price = GridPrice(
    "JiangSu_7_8",          # 名称
    "JiangSu",              # 城市
    "summer",               # 季节（覆盖季节，跨季用逗号分隔，首个=调度表计算轮次）
    JiangSu_7_8_hourly_tariff_ori,
    0.5,
    48
)

# 夏冬季 + 尖峰（12月-次年1月: 18:00-20:00）
JiangSu_12_1_hourly_tariff_ori = copy(JiangSu_6_8_12_2_hourly_tariff_ori)
JiangSu_12_1_hourly_tariff_ori[37:40] .= pp         # 尖峰 18:00-20:00
JiangSu_12_1_grid_price = GridPrice(
    "JiangSu_12_1",         # 名称
    "JiangSu",              # 城市
    "winter",               # 季节
    JiangSu_12_1_hourly_tariff_ori,
    0.5,
    48
)

# 春、秋两季（3-5月、9-11月），无尖峰
JiangSu_3_5_9_11_hourly_tariff_ori = ones(48)
JiangSu_3_5_9_11_hourly_tariff_ori[5:8] .*= v       # 谷 2:00-6:00
JiangSu_3_5_9_11_hourly_tariff_ori[21:28] .*= v     # 谷 10:00-14:00
JiangSu_3_5_9_11_hourly_tariff_ori[31:44] .*= p     # 峰 15:00-22:00
JiangSu_3_5_9_11_grid_price = GridPrice(
    "JiangSu_3_5_9_11",     # 名称
    "JiangSu",              # 城市
    "spring,autumn",        # 季节（春+秋）
    JiangSu_3_5_9_11_hourly_tariff_ori,
    0.5,
    48
)

# =============================================================================
# 浙江省工商业分时电价（浙发改价格〔2026〕112号，2026-07-01起执行）
# 参考: https://www.zj.gov.cn/col/col1229203589/art/2026/art_db07d88b6c68f1412561a0e8640cc135.html
# 浮动比例: 尖峰:高峰:平段:低谷:深谷 = 2.05 : 1.85 : 1 : 0.4 : 0.2
#   → 尖峰 pp=2.05，峰 p=1.85，谷 v=0.4，深谷 vv=0.2
# 时段划分（hourly_tariff_ori[i] 对应 (i-1):00 至 i:00）：
#   春秋季（2-6月、9-11月）:
#     高峰 16:00-23:00；平段 7:00-11:00、14:00-16:00、23:00-24:00；低谷 0:00-7:00、11:00-14:00
#   夏冬季（1月、7月、8月、12月）:
#     尖峰 18:00-22:00；高峰 16:00-18:00、22:00-23:00；平/谷同春秋季
# 注: 节假日（劳动节/国庆节前三天、春节）设深谷 9:00-15:00（vv=0.2），此处不单独建模
# =============================================================================
pp = 2.05   # 尖峰: 2.05
p = 1.85    # 高峰: 1.85
v = 0.4     # 低谷: 0.4
vv = 0.2    # 深谷: 0.2

# 春秋季（2-6月、9-11月），无尖峰
ZheJiang_2_6_9_11_hourly_tariff_ori = ones(24)
ZheJiang_2_6_9_11_hourly_tariff_ori[1:7] *= v        # 谷 0:00-7:00
ZheJiang_2_6_9_11_hourly_tariff_ori[12:14] *= v      # 谷 11:00-14:00
ZheJiang_2_6_9_11_hourly_tariff_ori[17:23] *= p      # 峰 16:00-23:00
ZheJiang_2_6_9_11_grid_price = GridPrice(
    "ZheJiang_2_6_9_11",    # 名称
    "ZheJiang",             # 城市
    "spring,autumn",        # 季节（春+秋）
    ZheJiang_2_6_9_11_hourly_tariff_ori,
    1.0,
    24
)

# 夏冬季（1月、7月、8月、12月），含尖峰
ZheJiang_1_7_8_12_hourly_tariff_ori = ones(24)
ZheJiang_1_7_8_12_hourly_tariff_ori[1:7] *= v        # 谷 0:00-7:00
ZheJiang_1_7_8_12_hourly_tariff_ori[12:14] *= v      # 谷 11:00-14:00
ZheJiang_1_7_8_12_hourly_tariff_ori[17:18] *= p      # 峰 16:00-18:00
ZheJiang_1_7_8_12_hourly_tariff_ori[19:22] *= pp     # 尖峰 18:00-22:00
ZheJiang_1_7_8_12_hourly_tariff_ori[23] *= p         # 峰 22:00-23:00
ZheJiang_1_7_8_12_grid_price = GridPrice(
    "ZheJiang_1_7_8_12",    # 名称
    "ZheJiang",             # 城市
    "summer,winter",        # 季节（夏+冬）
    ZheJiang_1_7_8_12_hourly_tariff_ori,
    1.0,
    24
)

# =============================================================================
# 山东省工商业分时电价（2026年，依据济南市发改委公告）
# 参考: https://jndpc.jinan.gov.cn/col2191/art/2025/art_2191_4789557.html
# 浮动比例（相对平时段 = 1.0）：
#   高峰上浮70%   → 峰   p  = 1.7
#   尖峰上浮100%  → 尖峰 pp = 2.0
#   低谷下浮70%   → 谷   v  = 0.3
#   深谷下浮90%   → 深谷 vv = 0.1
# 时段划分（hourly_tariff_ori[i] 对应 (i-1):00 至 i:00）：
#   峰/尖峰先按峰设置，再按尖峰覆盖；谷/深谷同理
# 注：1-2月时段与12月相同
# =============================================================================
p = 1.7    # 峰: 上浮70%
pp = 2.0   # 尖峰: 上浮100%
v = 0.3    # 谷: 下浮70%
vv = 0.1   # 深谷: 下浮90%

# 1月至2月、12月
# 谷: 2:00-6:00、10:00-15:00；深谷: 11:00-14:00
# 峰: 7:00-9:00、16:00-21:00；尖峰: 16:00-19:00
ShanDong_1_2_hourly_tariff_ori = ones(24)
ShanDong_1_2_hourly_tariff_ori[3:6] *= v        # 谷 2:00-6:00
ShanDong_1_2_hourly_tariff_ori[11:15] *= v      # 谷 10:00-15:00
ShanDong_1_2_hourly_tariff_ori[12:14] *= vv     # 深谷 11:00-14:00
ShanDong_1_2_hourly_tariff_ori[8:9] *= p        # 峰 7:00-9:00
ShanDong_1_2_hourly_tariff_ori[17:21] *= p      # 峰 16:00-21:00
ShanDong_1_2_hourly_tariff_ori[17:19] *= pp     # 尖峰 16:00-19:00
ShanDong_1_2_grid_price = GridPrice(
    "ShanDong_1_2",         # 名称
    "ShanDong",             # 城市
    "winter",               # 季节
    ShanDong_1_2_hourly_tariff_ori,
    1.0,
    24
)

# 3月至5月
# 谷: 10:00-15:00；深谷: 11:00-14:00
# 峰: 17:00-22:00；尖峰: 17:00-20:00
ShanDong_3_5_hourly_tariff_ori = ones(24)
ShanDong_3_5_hourly_tariff_ori[11:15] *= v      # 谷 10:00-15:00
ShanDong_3_5_hourly_tariff_ori[12:14] *= vv     # 深谷 11:00-14:00
ShanDong_3_5_hourly_tariff_ori[18:22] *= p      # 峰 17:00-22:00
ShanDong_3_5_hourly_tariff_ori[18:20] *= pp     # 尖峰 17:00-20:00
ShanDong_3_5_grid_price = GridPrice(
    "ShanDong_3_5",         # 名称
    "ShanDong",             # 城市
    "spring",               # 季节
    ShanDong_3_5_hourly_tariff_ori,
    1.0,
    24
)

# 6月
# 谷: 7:00-12:00
# 峰: 16:00-23:00；尖峰: 17:00-22:00（无深谷）
ShanDong_6_hourly_tariff_ori = ones(24)
ShanDong_6_hourly_tariff_ori[8:12] *= v         # 谷 7:00-12:00
ShanDong_6_hourly_tariff_ori[17:23] *= p        # 峰 16:00-23:00
ShanDong_6_hourly_tariff_ori[18:22] *= pp       # 尖峰 17:00-22:00
ShanDong_6_grid_price = GridPrice(
    "ShanDong_6",           # 名称
    "ShanDong",             # 城市
    "summer",               # 季节
    ShanDong_6_hourly_tariff_ori,
    1.0,
    24
)

# 7月至8月
# 谷: 1:00-6:00
# 峰: 16:00-23:00；尖峰: 17:00-22:00（无深谷）
ShanDong_7_8_hourly_tariff_ori = ones(24)
ShanDong_7_8_hourly_tariff_ori[2:6] *= v        # 谷 1:00-6:00
ShanDong_7_8_hourly_tariff_ori[17:23] *= p      # 峰 16:00-23:00
ShanDong_7_8_hourly_tariff_ori[18:22] *= pp     # 尖峰 17:00-22:00
ShanDong_7_8_grid_price = GridPrice(
    "ShanDong_7_8",         # 名称
    "ShanDong",             # 城市
    "summer",               # 季节
    ShanDong_7_8_hourly_tariff_ori,
    1.0,
    24
)

# 9月至11月
# 谷: 10:00-15:00；深谷: 11:00-14:00
# 峰: 16:00-21:00；尖峰: 17:00-19:00
ShanDong_9_11_hourly_tariff_ori = ones(24)
ShanDong_9_11_hourly_tariff_ori[11:15] *= v     # 谷 10:00-15:00
ShanDong_9_11_hourly_tariff_ori[12:14] *= vv    # 深谷 11:00-14:00
ShanDong_9_11_hourly_tariff_ori[17:21] *= p     # 峰 16:00-21:00
ShanDong_9_11_hourly_tariff_ori[18:19] *= pp    # 尖峰 17:00-19:00
ShanDong_9_11_grid_price = GridPrice(
    "ShanDong_9_11",        # 名称
    "ShanDong",             # 城市
    "autumn",               # 季节
    ShanDong_9_11_hourly_tariff_ori,
    1.0,
    24
)

# =============================================================================
# 广东省工商业峰谷分时电价（粤发改价格〔2021〕331号）
# 参考: https://drc.gd.gov.cn/ywtz/content/post_3500421.html
# 时段划分（全年统一）：
#   高峰: 10-12点、14-19点
#   低谷: 0-8点
#   其余: 平时段
# 峰平谷比价 = 1.7 : 1 : 0.38（峰 p=1.7，谷 v=0.38）
# 尖峰（仅7-9月及高温日，执行时段 11-12时、15-17时共3小时）：
#   在峰段基础上上浮25% → pp = 1.7 × 1.25 = 2.125
# =============================================================================
v = 0.38    # 谷: 峰平谷比价 1.7:1:0.38
pp = 2.125  # 尖峰: 峰段上浮25% = 1.7 × 1.25

# 全年基础（无尖峰）
GuangDong_hourly_tariff_ori = ones(24)
GuangDong_hourly_tariff_ori[1:8] *= v          # 低谷 0:00-8:00
GuangDong_hourly_tariff_ori[11:12] *= p        # 高峰 10:00-12:00
GuangDong_hourly_tariff_ori[15:19] *= p        # 高峰 14:00-19:00
GuangDong_grid_price = GridPrice(
    "GuangDong",            # 名称
    "GuangDong",            # 城市
    "other",                # 季节（全年统一）
    GuangDong_hourly_tariff_ori,
    1.0,
    24
)

# 7-9月（含尖峰：11-12时、15-17时）
GuangDong_7_9_hourly_tariff_ori = ones(24)
GuangDong_7_9_hourly_tariff_ori[1:8] *= v      # 低谷 0:00-8:00
GuangDong_7_9_hourly_tariff_ori[11:12] *= p    # 高峰 10:00-12:00
GuangDong_7_9_hourly_tariff_ori[15:19] *= p    # 高峰 14:00-19:00
GuangDong_7_9_hourly_tariff_ori[12] *= pp      # 尖峰 11:00-12:00
GuangDong_7_9_hourly_tariff_ori[16:17] *= pp   # 尖峰 15:00-17:00
GuangDong_7_9_grid_price = GridPrice(
    "GuangDong_7_9",        # 名称
    "GuangDong",            # 城市
    "summer",               # 季节
    GuangDong_7_9_hourly_tariff_ori,
    1.0,
    24
)

# =============================================================================
# 山西省工商业分时电价（晋发改商品发〔2026〕15号，2026-05-01起执行）
# 参考: https://fgw.shanxi.gov.cn/sxfgwzwgk/sxsfgwxxgk/xxgkml/tz/202602/t20260228_10069158.shtml
# 浮动比例（相对平时段 = 1.0）：
#   高峰上浮60%   → 峰   p  = 1.6
#   低谷下浮55%   → 谷   v  = 0.45
#   尖峰 = 高峰基础上上浮20%  → 尖峰 pp = 1.6 × 1.2 = 1.92
#   深谷 = 低谷基础上下浮20%  → 深谷 vv = 0.45 × 0.8 = 0.36
# 时段划分（hourly_tariff_ori[i] 对应 (i-1):00 至 i:00）：
#   春/秋无尖峰；夏尖峰 18:00-21:00、冬尖峰 17:00-20:00
# 注: 节假日（春节/劳动节/国庆节）设深谷 13:00-15:00（vv=0.36），此处不单独建模
# =============================================================================
p = 1.6     # 峰: 上浮60%
pp = 1.92   # 尖峰: 高峰上浮20% = 1.6 × 1.2
v = 0.45    # 谷: 下浮55%
vv = 0.36   # 深谷: 低谷下浮20% = 0.45 × 0.8

# 春季（3-5月），无尖峰
# 峰: 6:00-8:00、17:00-24:00；谷: 2:00-5:00、9:00-15:00
ShanXi_3_5_hourly_tariff_ori = ones(24)
ShanXi_3_5_hourly_tariff_ori[3:5] *= v          # 谷 2:00-5:00
ShanXi_3_5_hourly_tariff_ori[10:15] *= v        # 谷 9:00-15:00
ShanXi_3_5_hourly_tariff_ori[7:8] *= p          # 峰 6:00-8:00
ShanXi_3_5_hourly_tariff_ori[18:24] *= p        # 峰 17:00-24:00
ShanXi_3_5_grid_price = GridPrice(
    "ShanXi_3_5",           # 名称
    "ShanXi",               # 城市
    "spring",               # 季节
    ShanXi_3_5_hourly_tariff_ori,
    1.0,
    24
)

# 夏季（6-8月），尖峰 18:00-21:00
# 峰: 6:00-8:00、18:00-24:00；谷: 2:00-5:00、10:00-15:00
ShanXi_6_8_hourly_tariff_ori = ones(24)
ShanXi_6_8_hourly_tariff_ori[3:5] *= v          # 谷 2:00-5:00
ShanXi_6_8_hourly_tariff_ori[11:15] *= v        # 谷 10:00-15:00
ShanXi_6_8_hourly_tariff_ori[7:8] *= p          # 峰 6:00-8:00
ShanXi_6_8_hourly_tariff_ori[19:24] *= p        # 峰 18:00-24:00
ShanXi_6_8_hourly_tariff_ori[19:21] .= pp       # 尖峰 18:00-21:00（覆盖峰）
ShanXi_6_8_grid_price = GridPrice(
    "ShanXi_6_8",           # 名称
    "ShanXi",               # 城市
    "summer",               # 季节
    ShanXi_6_8_hourly_tariff_ori,
    1.0,
    24
)

# 秋季（9-11月），无尖峰
# 峰: 6:00-8:00、17:00-24:00；谷: 2:00-5:00、11:00-15:00
ShanXi_9_11_hourly_tariff_ori = ones(24)
ShanXi_9_11_hourly_tariff_ori[3:5] *= v          # 谷 2:00-5:00
ShanXi_9_11_hourly_tariff_ori[12:15] *= v        # 谷 11:00-15:00
ShanXi_9_11_hourly_tariff_ori[7:8] *= p          # 峰 6:00-8:00
ShanXi_9_11_hourly_tariff_ori[18:24] *= p        # 峰 17:00-24:00
ShanXi_9_11_grid_price = GridPrice(
    "ShanXi_9_11",          # 名称
    "ShanXi",               # 城市
    "autumn",               # 季节
    ShanXi_9_11_hourly_tariff_ori,
    1.0,
    24
)

# 冬季（12月、1-2月），尖峰 17:00-20:00
# 峰: 6:00-8:00、16:00-23:00；谷: 2:00-5:00、10:00-15:00
ShanXi_12_1_2_hourly_tariff_ori = ones(24)
ShanXi_12_1_2_hourly_tariff_ori[3:5] *= v        # 谷 2:00-5:00
ShanXi_12_1_2_hourly_tariff_ori[11:15] *= v      # 谷 10:00-15:00
ShanXi_12_1_2_hourly_tariff_ori[7:8] *= p        # 峰 6:00-8:00
ShanXi_12_1_2_hourly_tariff_ori[17:23] *= p      # 峰 16:00-23:00
ShanXi_12_1_2_hourly_tariff_ori[18:20] .= pp     # 尖峰 17:00-20:00（覆盖峰）
ShanXi_12_1_2_grid_price = GridPrice(
    "ShanXi_12_1_2",        # 名称
    "ShanXi",               # 城市
    "winter",               # 季节
    ShanXi_12_1_2_hourly_tariff_ori,
    1.0,
    24
)

# =============================================================================
# 安徽省工商业分时电价（2025-07-01起执行）
# 参考: https://fzggw.ah.gov.cn/public/7011/149979901.html
# 浮动比例：
#   低谷下浮61.8%           → 谷   v        = 0.382
#   春秋季高峰上浮74%       → 峰   p_spring = 1.74
#   夏冬季高峰上浮84.3%     → 峰   p_sw     = 1.843
#   尖峰 = 高峰上浮20%      → 尖峰 pp       = 1.843 × 1.2 = 2.2116
#   深谷 = 低谷下浮20%      → 深谷 vv       = 0.382 × 0.8 = 0.3056
# 尖峰执行期（仅315千伏安及以上两部制工业）：
#   每年7月15日-8月31日 20:00-22:00；12月15日-次年1月31日 19:00-21:00
# 注: 节假日深谷 11:00-15:00（vv=0.3056），此处不单独建模
# =============================================================================
v = 0.382        # 谷: 下浮61.8%
p_spring = 1.74  # 春秋季峰: 上浮74%
p_sw = 1.843     # 夏冬季峰: 上浮84.3%
pp = 2.2116      # 尖峰: 高峰上浮20% = 1.843 × 1.2
vv = 0.3056      # 深谷: 低谷下浮20% = 0.382 × 0.8

# 春秋季（2-6月、10月、11月），无尖峰
# 峰: 6:00-8:00、16:00-22:00；谷: 11:00-14:00、23:00-次日6:00
AnHui_2_6_10_11_hourly_tariff_ori = ones(24)
AnHui_2_6_10_11_hourly_tariff_ori[1:6] *= v            # 谷 0:00-6:00
AnHui_2_6_10_11_hourly_tariff_ori[12:14] *= v          # 谷 11:00-14:00
AnHui_2_6_10_11_hourly_tariff_ori[24] *= v             # 谷 23:00-24:00
AnHui_2_6_10_11_hourly_tariff_ori[7:8] *= p_spring     # 峰 6:00-8:00
AnHui_2_6_10_11_hourly_tariff_ori[17:22] *= p_spring   # 峰 16:00-22:00
AnHui_2_6_10_11_grid_price = GridPrice(
    "AnHui_2_6_10_11",      # 名称
    "AnHui",                # 城市
    "spring,autumn",        # 季节（春+秋）
    AnHui_2_6_10_11_hourly_tariff_ori,
    1.0,
    24
)

# 夏季基础（7-9月），无尖峰
# 峰: 16:00-24:00；谷: 2:00-9:00、11:00-13:00
AnHui_7_9_hourly_tariff_ori = ones(24)
AnHui_7_9_hourly_tariff_ori[3:9] *= v          # 谷 2:00-9:00
AnHui_7_9_hourly_tariff_ori[12:13] *= v        # 谷 11:00-13:00
AnHui_7_9_hourly_tariff_ori[17:24] *= p_sw     # 峰 16:00-24:00
AnHui_7_9_grid_price = GridPrice(
    "AnHui_7_9",            # 名称
    "AnHui",                # 城市
    "summer",               # 季节
    AnHui_7_9_hourly_tariff_ori,
    1.0,
    24
)

# 夏季 + 尖峰（7月15日-8月31日: 20:00-22:00）
AnHui_7_8_hourly_tariff_ori = copy(AnHui_7_9_hourly_tariff_ori)
AnHui_7_8_hourly_tariff_ori[21:22] .= pp      # 尖峰 20:00-22:00（覆盖峰）
AnHui_7_8_grid_price = GridPrice(
    "AnHui_7_8",            # 名称
    "AnHui",                # 城市
    "summer",               # 季节
    AnHui_7_8_hourly_tariff_ori,
    1.0,
    24
)

# 冬季基础（1月、12月），无尖峰
# 峰: 15:00-23:00；谷: 12:00-14:00、23:00-次日6:00
AnHui_1_12_hourly_tariff_ori = ones(24)
AnHui_1_12_hourly_tariff_ori[1:6] *= v          # 谷 0:00-6:00
AnHui_1_12_hourly_tariff_ori[13:14] *= v        # 谷 12:00-14:00
AnHui_1_12_hourly_tariff_ori[24] *= v           # 谷 23:00-24:00
AnHui_1_12_hourly_tariff_ori[16:23] *= p_sw     # 峰 15:00-23:00
AnHui_1_12_grid_price = GridPrice(
    "AnHui_1_12",           # 名称
    "AnHui",                # 城市
    "winter",               # 季节
    AnHui_1_12_hourly_tariff_ori,
    1.0,
    24
)

# 冬季 + 尖峰（12月15日-次年1月31日: 19:00-21:00）
AnHui_12_1_hourly_tariff_ori = copy(AnHui_1_12_hourly_tariff_ori)
AnHui_12_1_hourly_tariff_ori[20:21] .= pp     # 尖峰 19:00-21:00（覆盖峰）
AnHui_12_1_grid_price = GridPrice(
    "AnHui_12_1",           # 名称
    "AnHui",                # 城市
    "winter",               # 季节
    AnHui_12_1_hourly_tariff_ori,
    1.0,
    24
)

# =============================================================================
# 云南省工商业分时电价（云发改价格〔2026〕51号，2026-03-01起执行）
# 参考: https://yndrc.yn.gov.cn/html/2026/tongzhigonggao_0120/25211.html
# 浮动比例: 高峰:平时:低谷 = 1.5 : 1 : 0.5（峰 p=1.5、谷 v=0.5）
# 时段划分（全年统一，无季节差异）：
#   高峰: 7:00-9:00、18:00-24:00
#   低谷: 2:00-6:00、12:00-16:00
#   平时: 0:00-2:00、6:00-7:00、9:00-12:00、16:00-18:00
# 注: 尖峰电价暂缓执行；无深谷时段
# =============================================================================
v = 0.5     # 谷: 下浮50%
p = 1.5     # 峰: 上浮50%

YunNan_hourly_tariff_ori = ones(24)
YunNan_hourly_tariff_ori[3:6] *= v          # 谷 2:00-6:00
YunNan_hourly_tariff_ori[13:16] *= v        # 谷 12:00-16:00
YunNan_hourly_tariff_ori[8:9] *= p          # 峰 7:00-9:00
YunNan_hourly_tariff_ori[19:24] *= p        # 峰 18:00-24:00
YunNan_grid_price = GridPrice(
    "YunNan",               # 名称
    "YunNan",               # 城市
    "other",                # 季节（全年统一）
    YunNan_hourly_tariff_ori,
    1.0,
    24
)

# =============================================================================
# 广西壮族自治区工商业分时电价（桂发改价格规〔2021〕1029号，2021-12-20起执行）
# 参考: http://fgw.gxzf.gov.cn/zfxxgkzl/wjzx/zyzc/sfjg/t10811420.shtml
# 浮动比例: 峰平谷比价 = 1.5 : 1 : 0.5（峰 p=1.5、谷 v=0.5）
#   尖峰 = 峰段上浮20% → 尖峰 pp = 1.5 × 1.2 = 1.8
# 时段划分（全年统一）：
#   高峰: 10:00-12:00、16:00-22:00
#   平段: 7:00-10:00、12:00-16:00、22:00-23:00
#   低谷: 23:00-24:00、0:00-7:00
# 尖峰（仅7、8、9、12月四个整月）: 11:00-12:00、17:00-18:00
# =============================================================================
v = 0.5     # 谷: 下浮50%
p = 1.5     # 峰: 上浮50%
pp = 1.8    # 尖峰: 峰段上浮20% = 1.5 × 1.2

# 基础（无尖峰）
GuangXi_hourly_tariff_ori = ones(24)
GuangXi_hourly_tariff_ori[1:7] *= v          # 谷 0:00-7:00
GuangXi_hourly_tariff_ori[24] *= v           # 谷 23:00-24:00
GuangXi_hourly_tariff_ori[11:12] *= p        # 峰 10:00-12:00
GuangXi_hourly_tariff_ori[17:22] *= p        # 峰 16:00-22:00
GuangXi_grid_price = GridPrice(
    "GuangXi",              # 名称
    "GuangXi",              # 城市
    "other",                # 季节（全年统一）
    GuangXi_hourly_tariff_ori,
    1.0,
    24
)

# 尖峰月（7、8、9、12月）: 尖峰 11:00-12:00、17:00-18:00
GuangXi_7_8_9_12_hourly_tariff_ori = copy(GuangXi_hourly_tariff_ori)
GuangXi_7_8_9_12_hourly_tariff_ori[12] = pp     # 尖峰 11:00-12:00（覆盖峰）
GuangXi_7_8_9_12_hourly_tariff_ori[18] = pp     # 尖峰 17:00-18:00（覆盖峰）
GuangXi_7_8_9_12_grid_price = GridPrice(
    "GuangXi_7_8_9_12",     # 名称
    "GuangXi",              # 城市
    "summer,winter",        # 季节（夏+冬）
    GuangXi_7_8_9_12_hourly_tariff_ori,
    1.0,
    24
)

# =============================================================================
# 上海市分时电价（沪发改价管〔2022〕50号，2023-01-01起执行）
# 参考: https://fgw.sh.gov.cn/fgw_jggl/20221216/e2652e3ab7ee49438d6e82af8880b160.html
# 两类用户、浮动比例分夏冬/其他两档：
#   两部制/大工业:
#     夏冬(1、7、8、9、12月): 峰上浮80% p=1.8，谷下浮60% v=0.4，尖峰=峰×1.25=2.25
#     其他月:                  峰上浮60% p=1.6，谷下浮50% v=0.5
#     尖峰: 7-8月 12:00-14:00；冬季(1、12月) 19:00-21:00
#   单一制:
#     峰 6:00-22:00、谷 22:00-次日6:00
#     夏冬峰上浮20% p=1.2，其他月峰上浮17% p=1.17，谷均下浮45% v=0.55
# =============================================================================

# ---------- 两部制/大工业 ----------
p_two_sw = 1.8      # 夏冬季峰: 上浮80%
p_two_o = 1.6       # 其他月峰: 上浮60%
v_two_sw = 0.4      # 夏冬季谷: 下浮60%
v_two_o = 0.5       # 其他月谷: 下浮50%
pp_two = 2.25       # 尖峰: 峰段上浮25% = 1.8 × 1.25

# 夏季基础（7-9月），无尖峰
# 峰: 8:00-15:00、18:00-21:00；平: 6:00-8:00、15:00-18:00、21:00-22:00；谷: 22:00-次日6:00
ShangHai_7_9_hourly_tariff_ori = ones(24)
ShangHai_7_9_hourly_tariff_ori[1:6] *= v_two_sw       # 谷 0:00-6:00
ShangHai_7_9_hourly_tariff_ori[23:24] *= v_two_sw     # 谷 22:00-24:00
ShangHai_7_9_hourly_tariff_ori[9:15] *= p_two_sw      # 峰 8:00-15:00
ShangHai_7_9_hourly_tariff_ori[19:21] *= p_two_sw     # 峰 18:00-21:00
ShangHai_7_9_grid_price = GridPrice(
    "ShangHai_7_9",         # 名称
    "ShangHai",             # 城市
    "summer",               # 季节
    ShangHai_7_9_hourly_tariff_ori,
    1.0,
    24
)

# 夏季 + 尖峰（7、8月: 12:00-14:00）
ShangHai_7_8_hourly_tariff_ori = copy(ShangHai_7_9_hourly_tariff_ori)
ShangHai_7_8_hourly_tariff_ori[13:14] .= pp_two      # 尖峰 12:00-14:00（覆盖峰）
ShangHai_7_8_grid_price = GridPrice(
    "ShangHai_7_8",         # 名称
    "ShangHai",             # 城市
    "summer",               # 季节
    ShangHai_7_8_hourly_tariff_ori,
    1.0,
    24
)

# 冬季（1、12月），尖峰 19:00-21:00
# 峰: 8:00-11:00、18:00-21:00；平: 6:00-8:00、11:00-18:00、21:00-22:00；谷: 22:00-次日6:00
ShangHai_1_12_hourly_tariff_ori = ones(24)
ShangHai_1_12_hourly_tariff_ori[1:6] *= v_two_sw       # 谷 0:00-6:00
ShangHai_1_12_hourly_tariff_ori[23:24] *= v_two_sw     # 谷 22:00-24:00
ShangHai_1_12_hourly_tariff_ori[9:11] *= p_two_sw      # 峰 8:00-11:00
ShangHai_1_12_hourly_tariff_ori[19:21] *= p_two_sw     # 峰 18:00-21:00
ShangHai_1_12_hourly_tariff_ori[20:21] .= pp_two       # 尖峰 19:00-21:00（覆盖峰）
ShangHai_1_12_grid_price = GridPrice(
    "ShangHai_1_12",        # 名称
    "ShangHai",             # 城市
    "winter",               # 季节
    ShangHai_1_12_hourly_tariff_ori,
    1.0,
    24
)

# 其他月份（2-6月、10-11月），无尖峰
ShangHai_2_6_10_11_hourly_tariff_ori = ones(24)
ShangHai_2_6_10_11_hourly_tariff_ori[1:6] *= v_two_o       # 谷 0:00-6:00
ShangHai_2_6_10_11_hourly_tariff_ori[23:24] *= v_two_o     # 谷 22:00-24:00
ShangHai_2_6_10_11_hourly_tariff_ori[9:11] *= p_two_o      # 峰 8:00-11:00
ShangHai_2_6_10_11_hourly_tariff_ori[19:21] *= p_two_o     # 峰 18:00-21:00
ShangHai_2_6_10_11_grid_price = GridPrice(
    "ShangHai_2_6_10_11",   # 名称
    "ShangHai",             # 城市
    "spring,autumn",        # 季节（春+秋）
    ShangHai_2_6_10_11_hourly_tariff_ori,
    1.0,
    24
)

# ---------- 单一制 ----------
# 峰: 6:00-22:00；谷: 22:00-次日6:00
p_single_sw = 1.2     # 夏冬峰: 上浮20%
p_single_o = 1.17     # 其他月峰: 上浮17%
v_single = 0.55       # 谷: 下浮45%（全年统一）

# 夏冬季（1、7、8、9、12月）
ShangHai_single_1_7_8_9_12_hourly_tariff_ori = ones(24)
ShangHai_single_1_7_8_9_12_hourly_tariff_ori[1:6] *= v_single       # 谷 0:00-6:00
ShangHai_single_1_7_8_9_12_hourly_tariff_ori[23:24] *= v_single     # 谷 22:00-24:00
ShangHai_single_1_7_8_9_12_hourly_tariff_ori[7:22] *= p_single_sw   # 峰 6:00-22:00
ShangHai_single_1_7_8_9_12_grid_price = GridPrice(
    "ShangHai_single_1_7_8_9_12",   # 名称
    "ShangHai",                     # 城市
    "summer,winter",                # 季节（夏+冬）
    ShangHai_single_1_7_8_9_12_hourly_tariff_ori,
    1.0,
    24
)

# 其他月份（2-6月、10-11月）
ShangHai_single_2_6_10_11_hourly_tariff_ori = ones(24)
ShangHai_single_2_6_10_11_hourly_tariff_ori[1:6] *= v_single       # 谷 0:00-6:00
ShangHai_single_2_6_10_11_hourly_tariff_ori[23:24] *= v_single     # 谷 22:00-24:00
ShangHai_single_2_6_10_11_hourly_tariff_ori[7:22] *= p_single_o    # 峰 6:00-22:00
ShangHai_single_2_6_10_11_grid_price = GridPrice(
    "ShangHai_single_2_6_10_11",    # 名称
    "ShangHai",                     # 城市
    "spring,autumn",                # 季节（春+秋）
    ShangHai_single_2_6_10_11_hourly_tariff_ori,
    1.0,
    24
)

# =============================================================================
# 陕西省分时电价（陕发改价格〔2025〕1034号，2025-08-01起执行）
# 参考: https://sndrc.shaanxi.gov.cn/zfxxgk/zc/fgwj/sfzggwwj/2025/202507/t20250723_3546317.html
# 浮动比例（代理购电工商业用电，以平段电价为基准）：
#   高峰上浮70% → 峰 p = 1.7
#   低谷下浮70% → 谷 v = 0.3
#   尖峰在平段基础上上浮90% → 尖峰 pp = 1.9
# 时段划分（全年统一）：
#   高峰: 16:00-23:00；低谷: 0:00-6:00、11:00-14:00；平段: 6:00-11:00、14:00-16:00、23:00-24:00
# 尖峰（迎峰度夏/冬）:
#   夏季(7、8月): 19:00-21:00；冬季(1、12月): 18:00-20:00
# 注: 农业生产用电浮动比例不同（峰上浮50%、谷下浮50%），此处不单独建模
# =============================================================================
p = 1.7     # 峰: 上浮70%
v = 0.3     # 谷: 下浮70%
pp = 1.9    # 尖峰: 平段基础上上浮90%

# 基础（其他月份，无尖峰）
ShaanXi_hourly_tariff_ori = ones(24)
ShaanXi_hourly_tariff_ori[1:6] *= v          # 谷 0:00-6:00
ShaanXi_hourly_tariff_ori[12:14] *= v        # 谷 11:00-14:00
ShaanXi_hourly_tariff_ori[17:23] *= p        # 峰 16:00-23:00
ShaanXi_grid_price = GridPrice(
    "ShaanXi",              # 名称
    "ShaanXi",              # 城市
    "other",                # 季节（全年统一）
    ShaanXi_hourly_tariff_ori,
    1.0,
    24
)

# 夏季（7、8月），尖峰 19:00-21:00
ShaanXi_7_8_hourly_tariff_ori = copy(ShaanXi_hourly_tariff_ori)
ShaanXi_7_8_hourly_tariff_ori[20:21] .= pp   # 尖峰 19:00-21:00（覆盖峰）
ShaanXi_7_8_grid_price = GridPrice(
    "ShaanXi_7_8",          # 名称
    "ShaanXi",              # 城市
    "summer",               # 季节
    ShaanXi_7_8_hourly_tariff_ori,
    1.0,
    24
)

# 冬季（1、12月），尖峰 18:00-20:00
ShaanXi_1_12_hourly_tariff_ori = copy(ShaanXi_hourly_tariff_ori)
ShaanXi_1_12_hourly_tariff_ori[19:20] .= pp  # 尖峰 18:00-20:00（覆盖峰）
ShaanXi_1_12_grid_price = GridPrice(
    "ShaanXi_1_12",         # 名称
    "ShaanXi",              # 城市
    "winter",               # 季节
    ShaanXi_1_12_hourly_tariff_ori,
    1.0,
    24
)

# =============================================================================
# 分时电价对象清单（共33个，全部对象）
# 用途: 作为 GRID_PRICE_BY_NAME 的数据源；实际计算顺序由
#       build_grid_price_order() 依据 GRID_PRICE_SCHEDULE 决定（季节轮次制）
# =============================================================================
GRID_PRICE_ORDER = [
    # ---- 第1轮: 夏/冬/春秋/全年 ----
    ZheJiang_1_7_8_12_grid_price,           # 夏
    ShanDong_1_2_grid_price,                # 冬
    ZheJiang_2_6_9_11_grid_price,           # 春秋
    GuangDong_grid_price,                   # 全年
    # ---- 第2轮 ----
    ShanDong_7_8_grid_price,                # 夏
    JiangSu_12_1_grid_price,                # 冬
    ShanDong_3_5_grid_price,                # 春秋
    YunNan_grid_price,                      # 全年
    # ---- 第3轮 ----
    ShanDong_6_grid_price,                  # 夏
    AnHui_1_12_grid_price,                  # 冬
    ShanDong_9_11_grid_price,               # 春秋
    GuangXi_grid_price,                     # 全年
    # ---- 第4轮 ----
    JiangSu_7_8_grid_price,                 # 夏
    AnHui_12_1_grid_price,                  # 冬
    JiangSu_3_5_9_11_grid_price,            # 春秋
    ShaanXi_grid_price,                     # 全年
    # ---- 第5轮（全年桶取完） ----
    GuangDong_7_9_grid_price,               # 夏
    ShanXi_12_1_2_grid_price,               # 冬
    AnHui_2_6_10_11_grid_price,             # 春秋
    # ---- 第6轮 ----
    AnHui_7_8_grid_price,                   # 夏
    ShangHai_1_12_grid_price,               # 冬
    ShanXi_3_5_grid_price,                  # 春秋
    # ---- 第7轮 ----
    AnHui_7_9_grid_price,                   # 夏
    ShaanXi_1_12_grid_price,                # 冬
    ShanXi_9_11_grid_price,                 # 春秋
    # ---- 第8轮 ----
    ShanXi_6_8_grid_price,                  # 夏
    ShangHai_single_1_7_8_9_12_grid_price,  # 冬
    ShangHai_2_6_10_11_grid_price,          # 春秋
    # ---- 第9轮（冬桶取完） ----
    ShangHai_7_8_grid_price,                # 夏
    ShangHai_single_2_6_10_11_grid_price,   # 春秋
    # ---- 剩余（仅夏桶） ----
    ShangHai_7_9_grid_price,                # 夏
    GuangXi_7_8_9_12_grid_price,            # 夏
    ShaanXi_7_8_grid_price,                 # 夏
]

# =============================================================================
# 分时电价季节调度表（手动编辑）
# - 每个键为一个计算轮次（夏/冬/春/秋/其他），列表内对象按顺序计算
# - 每城市每季节只放1条代表曲线（优先尖峰），其余对象放入 :other 最后计算
# - 对象名跨轮重复时，build_grid_price_order() 自动去重（只计算首次出现的轮次）
# 手动调节方式: 直接编辑本表（调整轮内顺序 / 改变对象归属轮次 / 增删对象）
# =============================================================================
const GRID_PRICE_SCHEDULE = Dict(
    :summer => [
        "ZheJiang_1_7_8_12",    # 浙江 夏+冬（含尖峰）
        "ShanDong_7_8",         # 山东 夏（尖峰）
        "JiangSu_7_8",          # 江苏 夏（尖峰）
        "GuangDong_7_9",        # 广东 夏（尖峰）
        "AnHui_7_8",            # 安徽 夏（尖峰）
        "ShanXi_6_8",           # 山西 夏（尖峰）
        "ShangHai_7_8",         # 上海 夏（尖峰）
        "GuangXi_7_8_9_12",     # 广西 夏+冬（尖峰）
        "ShaanXi_7_8",          # 陕西 夏（尖峰）
    ],
    :winter => [
        "ShanDong_1_2",         # 山东 冬（尖峰）
        "JiangSu_12_1",         # 江苏 冬（尖峰）
        "AnHui_12_1",           # 安徽 冬（尖峰）
        "ShanXi_12_1_2",        # 山西 冬（尖峰）
        "ShangHai_1_12",        # 上海 冬（尖峰）
        "ShaanXi_1_12",         # 陕西 冬（尖峰）
    ],
    :spring => [
        "ZheJiang_2_6_9_11",    # 浙江 春+秋
        "ShanDong_3_5",         # 山东 春
        "JiangSu_3_5_9_11",     # 江苏 春+秋
        "AnHui_2_6_10_11",      # 安徽 春+秋
        "ShanXi_3_5",           # 山西 春
        "ShangHai_2_6_10_11",   # 上海 春+秋
    ],
    :autumn => [
        "ShanDong_9_11",        # 山东 秋
        "ShanXi_9_11",          # 山西 秋
    ],
    :other => [   # 剩余对象（同城市同季节次选 + 全年基础），顺序任意
        "ShanDong_6",                  # 山东 夏（次选）
        "AnHui_7_9",                   # 安徽 夏基础（次选）
        "AnHui_1_12",                  # 安徽 冬基础（次选）
        "ShangHai_7_9",                # 上海 夏基础（次选）
        "ShangHai_single_1_7_8_9_12",  # 上海 单一制夏冬
        "ShangHai_single_2_6_10_11",   # 上海 单一制春秋
        "GuangDong",                   # 广东 全年
        "YunNan",                      # 云南 全年
        "GuangXi",                     # 广西 全年
        "ShaanXi",                     # 陕西 全年
    ],
)

# 电价对象名 → GridPrice 映射（数据源: GRID_PRICE_ORDER）
const GRID_PRICE_BY_NAME = Dict(gp.name => gp for gp in GRID_PRICE_ORDER)

"""
    获取电价对象所属的季节轮（与计算顺序一致，取首次出现的轮次）

参数:
- name: 电价曲线名称（如 "ZheJiang_1_7_8_12"）

返回:
- 季节轮 Symbol（:summer / :winter / :spring / :autumn / :other）
"""
function get_gp_season(name::String)
    for season in [:summer, :winter, :spring, :autumn, :other]
        if name in get(GRID_PRICE_SCHEDULE, season, String[])
            return season
        end
    end
    return :other
end

"""
    生成分时电价计算顺序（按季节轮次 + 自动去重）

排序规则:
1. 按季节轮: 夏 → 冬 → 春 → 秋 → 其他
2. 每轮内按 GRID_PRICE_SCHEDULE 声明的顺序（优先尖峰）
3. 对象名跨轮重复时只计算第一次出现的轮次（computed 集合去重）

返回:
- 按计算顺序排列的电价对象数组 Vector{GridPrice}
"""
function build_grid_price_order()
    result = GridPrice[]
    computed = Set{String}()
    for season in [:summer, :winter, :spring, :autumn, :other]
        for name in get(GRID_PRICE_SCHEDULE, season, String[])
            if haskey(GRID_PRICE_BY_NAME, name) && !(name in computed)
                push!(result, GRID_PRICE_BY_NAME[name])
                push!(computed, name)
            else
                @warn "跳过重复/未知电价对象: $name"
            end
        end
    end
    return result
end

export build_grid_price_order