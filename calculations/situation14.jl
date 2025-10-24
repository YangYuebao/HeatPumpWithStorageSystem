using Revise
using Pkg
using HeatPumpWithStorageSystem
using DataFrames, CSV
using JuMP, HiGHS

#=
测试线性规划求解器好使不好使
=#
situation = "situation14"

cpm=8/(220-150)
TCompressorIn=115.0
Tuse=170.0
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

heatStorageCapacity=3.0
Tsmin = 120.0
Tsmax = 220.0
TcChangeToElec=180.0
dT_EvaporationStandard = 5.0
cpm_h = heatStorageCapacity/(Tsmax-Tuse)

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

PWaterCompressorMax = heatPumpServiceCoff * heatStorageCapacity / COPWater(TCompressorIn,Tuse)
PelecHeatMax = heatStorageCapacity / 4.0

# 定义参数结构体
params = SystemParameters(
				ThMax = TcChangeToElec,
				Tuse = Tuse,
				dT = dT_EvaporationStandard,
				TCompressorIn = TCompressorIn,
				cpm = cpm_h,
				COPWater = COPWater,
				PhMax = PWaterCompressorMax,
				PeMax = PelecHeatMax,
				cp_cw = 4.275,
				latentHeat = 2150.0,
				Tsmin = Tsmin,
				Tsmax = Tsmax,
				dTRecycleSupply = 5.0,
				dTRecycleBackward = 5.0,
				sysStruct = RecycleStruct(1,0,1)
			)

sysVariables = HeatPumpWithStorageSystem.SystemVariables(
    heatLoad,
    COPLowFunction(TWaste,TCompressorIn+5.0),
    26,
    30
)

160.23854
168.25839

TsStart=160.23854
TsEnd=168.25839
dt=0.5

COPl=sysVariables.COPl


@time result = getMinimumCost(TsStart,TsEnd,dt,params,sysVariables)
@time resultMILP = HeatPumpWithStorageSystem.getMinimumCost_MILP(TsStart,TsEnd,dt,params,sysVariables)
@show result
@show resultMILP

tags = [
	"cost",
	"flag",
	"P1",
	"P2",
	"P3",
	"Pe"
]

differenct = collect(Float64.(result).-Float64.(resultMILP))
df = DataFrame(:tags => tags,:result =>collect(result),:resultMILP => collect(resultMILP), :differenct => differenct)

vscodedisplay(df)
