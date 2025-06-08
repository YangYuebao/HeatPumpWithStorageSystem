
using Revise
using HeatPumpWithStorageSystem
using Pkg

#=
无关性检验
=#
KList=1:1:30
n=length(KList)
costList=zeros(n)

for i in 1:n
begin
	hourlyTariff = zeros(24)
	hourlyTariff[1:7] .= 0.3340
	hourlyTariff[8:11] .= 0.7393
	hourlyTariff[12:14] .= 1.2360
	hourlyTariff[15:18] .= 0.7393
	hourlyTariff[19:23] .= 1.2360
	hourlyTariff[24] = 0.3340

	or = NH3_Water
	heatStorageCapacity = 10.0
	Tuse = 180.0

	# 系数
	heatPumpServiceCoff = 1.0
	maxCOP = 21# 最大COP
	eta_s = 0.7# 绝热效率
	workingStartHour = 0                # 生产开始时间
	workingHours = 24                   # 每日工作小时数
	TWaste = 30.0                     # 废热源温度
	Tair = fill(25.0,25)                        # 外部环境温度
	TCompressorIn = 115.0
	maxTcHigh = 180.0
	dT_EvaporationStandard = 5.0
	heatConsumptionPower = heatConsumptionPower = fill(1.0, 25)
	# 计算参数
	dT = 1.0
	dt = 1.0# 时间步长过小会导致初始温度优化的目标不是一个单峰函数
	K=KList[i]

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

	COPOverlap = getOverlapCOP_fixMidTemperature(
		or,
		TCompressorIn + or.midTDifference / 2;
		maxCOP = maxCOP,# 最大COP
		eta_s = eta_s,# 绝热效率
		dT = dT,# 插值步长
	)
	#COP1_design=COPOverlap(TWaste, Tuse)
	COP1_design=1.0
end

#=
			GoldenRatioMethod	
步长1/1	cost 6.040405104886531
步长1/6 cost 6.211648136391782
步长1/10 cost 6.13
步长1/15 cost 6.359773767716163
=#

hourlyTariffFunction, heatConsumptionPowerFunction, TairFunction,
Tuse, TCompressorIn,
dT_EvaporationStandard,
latentHeat, cp_cw, cp_cs,
TcChangeToElec, TWaste,
cpm_h,
TstorageTankMax, PheatPumpMax, PelecHeatMax = generateSystemCoff(HeatPumpWithStorageSystem.PressedWaterOneStorageOneCompressor();
	overlapRefrigerant = or,    # 复叠工质
	COP1_design=COP1_design,
	maxTcHigh = maxTcHigh,                  # 高温热泵冷凝器温度上限
	TCompressorIn = TCompressorIn,              # 中间温度
	TWaste = TWaste,                      # 废热源温度
	Tair = Tair,                        # 外部环境温度
	Tuse = Tuse,                       # 工厂使用温度
	heatStorageCapacity = heatStorageCapacity,  # 蓄热量kWh(承压水蓄热)
	TstorageTankMax = 220.0,            # 蓄热罐的最高温度
	maxheatStorageInputHour = 4,        # 蓄热充满时长
	dT_EvaporationStandard = dT_EvaporationStandard,       # 全蒸温差
	heatConsumptionPower = heatConsumptionPower,         # 每小时用热功率kW
	workingStartHour = workingStartHour,# 生产开始时间
	workingHours = workingHours,                  # 每日工作小时数
	heatPumpServiceCoff = heatPumpServiceCoff,    # 热泵功率服务系数
	hourlyTariff = hourlyTariff,       # 电价向量
)

minCostTest, minTsListTest, P1ListTest, P2ListTest, P3ListTest, PeListTest = generateAndSolve(PressedWaterOneStorageOneCompressor(), MinimizeCost(), ConstloadandArea(), GoldenRatioMethod();
	COPOverlap = COPOverlap,
	COPWater = COPWater,
	hourlyTariffFunction = hourlyTariffFunction,
	heatConsumptionPowerFunction = heatConsumptionPowerFunction,
	TairFunction = TairFunction,
	Tuse = Tuse,
	TCompressorIn = TCompressorIn,
	dT_EvaporationStandard = dT_EvaporationStandard,
	latentHeat = latentHeat,
	cp_cw = cp_cw,
	cp_cs = cp_cs,
	TcChangeToElec = TcChangeToElec,
	TWaste = TWaste,
	cpm_h = cpm_h,
	TstorageTankMax = TstorageTankMax,
	PheatPumpMax = PheatPumpMax,
	PelecHeatMax = PelecHeatMax,
	# 求解参数
	dT = dT,# 状态参数高温蓄热温度离散步长
	dt = dt,# 时间步长
	K=K
)
costList[i]=minCostTest
println("K=$K,cost=$minCostTest")
end


using CSV,DataFrames

df = DataFrame(
	"dTs"=>1.0./KList,
	"180℃"=>costList[:,3],
)

CSV.write(joinpath(pwd(),"calculations","situation9","dt_1h.csv"), df)