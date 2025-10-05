using Revise
using Pkg
using HeatPumpWithStorageSystem
using DataFrames, CSV

#=
测试混合整数线性规划求解器好使不好使
=#
situation = "situation14"

cpm=8/(220-150)
TCompressorIn=115.0
Tuse=150.0
TWaste = 30.0
maxTcHigh=180.0
Tair=25.0

heatLoad=1.0

maxCOP=21.0
eta_s=0.7
dT=1.0
dT_heatTransfer = 5.0

heatPumpServiceCoff=1.0
PeMax=1.0
or = NH3_Water
refLow = or.refrigerantLow


COPLowFunction = getCOP(
	refLow.minTe,# 蒸发温度下限,这里是实际设计中的蒸发冷凝温度界限
	refLow.maxTe,# 蒸发温度上限
	refLow.minTc,# 冷凝温度下限
	refLow.maxTc,# 冷凝温度上限
	refLow,# 工质
	maxCOP,# 最大COP
	eta_s,# 绝热效率
	dT,# 插值步长
)
COPOverlap = getOverlapCOP_fixMidTemperature(
	or,
	TCompressorIn + or.midTDifference / 2;
	maxCOP = maxCOP,# 最大COP
	eta_s = eta_s,# 绝热效率
	dT = 1.0,# 插值步长
)
COPWater = getCOP(
	TCompressorIn,# 蒸发温度下限,这里是实际设计中的蒸发冷凝温度界限
	maxTcHigh,# 蒸发温度上限
	TCompressorIn,# 冷凝温度下限
	maxTcHigh,# 冷凝温度上限
	refWater,# 工质
	maxCOP,# 最大COP
	eta_s,# 绝热效率
	dT,# 插值步长
)

# 定义参数结构体
params =  SystemParameters(
	ThMax=maxTcHigh,
	Tuse=Tuse,
	load=heatLoad,
	dT=dT_heatTransfer,
	TCompressorIn=TCompressorIn,
	cpm=cpm,
	COPWater=COPWater,
	COPl = COPOverlap(TWaste,Tuse),
	PhMax = heatPumpServiceCoff*heatLoad/COPWater(TCompressorIn,Tuse),
	PeMax =PeMax,
	Tair=Tair,
	TWaste=TWaste
)

TsStart=128.0
TsEnd=130.0
dt=0.5

@time getMinimumCost(TsStart,TsEnd,dt,params)
"""
COPl,coph1,coph2,coph3,coph4,
"""
result=[getMinimumCost(TsStart,TsEnd,dt,params)...]
