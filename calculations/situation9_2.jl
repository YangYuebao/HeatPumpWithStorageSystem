using Pkg
using Revise
using HeatPumpWithStorageSystem


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
	Tuse = 150.0

	# 系数
	heatPumpServiceCoff = 1.0
	maxCOP = 21# 最大COP
	eta_s = 0.7# 绝热效率
	workingStartHour = 0                # 生产开始时间
	workingHours = 24                   # 每日工作小时数
	TWaste = 30.0                     # 废热源温度
	Tair = 25.0                        # 外部环境温度
	TCompressorIn = 115.0
	maxTcHigh = 180.0
	dT_EvaporationStandard = 5.0
	heatConsumptionPower = 1.0
	# 计算参数
	dT = 0.025
	temp=1.0
	dt = 1/temp# 时间步长过小会导致初始温度优化的目标不是一个单峰函数
	K=temp

	COPWater = getCOP(
		TCompressorIn,# 蒸发温度下限,这里是实际设计中的蒸发冷凝温度界限
		maxTcHigh,# 蒸发温度上限
		TCompressorIn,# 冷凝温度下限
		maxTcHigh,# 冷凝温度上限
		refWater,# 工质
		maxCOP,# 最大COP
		eta_s,# 绝热效率
		min(dT,1.0),# 插值步长
	)

	COPLow = getCOP(
		-10.0,# 蒸发温度下限,这里是实际设计中的蒸发冷凝温度界限
		125.0,# 蒸发温度上限
		20.0,# 冷凝温度下限
		125,# 冷凝温度上限
		or.refrigerantLow,# 工质
		maxCOP,# 最大COP
		eta_s,# 绝热效率
		min(dT,1.0),# 插值步长
	)

	COPOverlap = getOverlapCOP_fixMidTemperature(
		or,
		TCompressorIn + or.midTDifference / 2;
		maxCOP = maxCOP,# 最大COP
		eta_s = eta_s,# 绝热效率
		dT = min(dT,1.0),# 插值步长
	)
	COP1_design=COPOverlap(TWaste, Tuse)
	#COP1_design=1.0
	COPWater_design = COPWater(TCompressorIn,Tuse)
	COPLow_design = COPLow(TWaste,TCompressorIn+dT_EvaporationStandard)
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
TstorageTankMax, PheatPumpMax, PelecHeatMax,
PWaterCompressorMax = generateSystemCoff(HeatPumpWithStorageSystem.PressedWaterOneStorageOneCompressor();
	overlapRefrigerant = or,    # 复叠工质
	COP1_design=COP1_design,
	COPWater_design=COPWater_design,
	maxTcHigh = maxTcHigh,                  # 高温热泵冷凝器温度上限
	TCompressorIn = TCompressorIn,              # 中间温度
	TWaste = TWaste,                      # 废热源温度
	Tair = fill(Tair,25),                        # 外部环境温度
	Tuse = Tuse,                       # 工厂使用温度
	heatStorageCapacity = heatStorageCapacity,  # 蓄热量kWh(承压水蓄热)
	TstorageTankMax = 220.0,            # 蓄热罐的最高温度
	maxheatStorageInputHour = 4,        # 蓄热充满时长
	dT_EvaporationStandard = dT_EvaporationStandard,       # 全蒸温差
	heatConsumptionPower = fill(heatConsumptionPower,25),         # 每小时用热功率kW
	workingStartHour = workingStartHour,# 生产开始时间
	workingHours = workingHours,                  # 每日工作小时数
	heatPumpServiceCoff = heatPumpServiceCoff,    # 热泵功率服务系数
	hourlyTariff = hourlyTariff,       # 电价向量
)


@time minCostTest, minTsListTest, P1ListTest, P2ListTest, P3ListTest, PeListTest, C = generateAndSolve(PressedWaterOneStorageOneCompressor(), MinimizeCost(), ConstloadandArea(), GoldenRatioMethod();
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


using Plots
using DataFrames, CSV

TsList=collect(120:dT:220)
df=DataFrame(hcat(TsList,C),vcat(["Ts"],string.(TsList).*"℃"))
CSV.write(joinpath(pwd(),"calculations","situation7","C.csv"), df)


tList = 0:dt:(24-dt)

# 看看不同的初始温度下运行费用的变化

# 看看最优运行费用下蓄热温度的变化
plt = plot(tList .+ 0.5, [P1ListTest P2ListTest P3ListTest PeListTest], label = ["P1" "P2" "P3" "Pe"],ylims=(0,2.5))
plot!(plt, tList, minTsListTest[1:end-1] / 220, label = "Ts")



#=
dt	K	Tuse	Storage		Cost		Ts0
1/3 12	170		3.0			6.68052		132.18
1.0 12	170		3.0			6.68007		137.83
1/2	6	170		3.0			6.68341		140.5
1.0	1	170		3.0			6.69800		130.3
1/2	2	170		3.0			6.70000		139.35
1/2 8	170		3.0			6.68623		134.68
=#

