
using Revise
using HeatPumpWithStorageSystem

begin
	# 夏季7-8月,315千伏安以上，有尖峰
	hourlyTariff = ones(48)
	p=1.7
	pp=p*1.2
	v=0.35
	vv=v*0.8

	hourlyTariff[1:12] .= v
	hourlyTariff[23:26] .= v
	hourlyTariff[29:30] .= pp
	hourlyTariff[31:39] .= p
	hourlyTariff[40:43] .= pp
	hourlyTariff[44] = p


	Tair = vcat(
		fill(26.0, 7),
		fill(27.0, 7),
		fill(26.0,11)
	)

	#=
	hourlyTariff = zeros(24)
	hourlyTariff[1:7] .= 0.3340
	hourlyTariff[8:11] .= 0.7393
	hourlyTariff[12:14] .= 1.2360
	hourlyTariff[15:18] .= 0.7393
	hourlyTariff[19:23] .= 1.2360
	hourlyTariff[24] = 0.3340

	Tair = fill(25.0, 25)
	=#

	heatConsumptionPower = fill(1.0, 25)

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
	#Tair = 25.0                        # 外部环境温度
	TCompressorIn = 115.0
	maxTcHigh = 180.0
	dT_EvaporationStandard = 5.0
	#heatConsumptionPower = 1.0
	# 计算参数
	dT = 0.1
	dt = 1/2# 时间步长过小会导致初始温度优化的目标不是一个单峰函数
	smoother=0.001

	COPWater = getCOP(
		TCompressorIn,# 蒸发温度下限,这里是实际设计中的蒸发冷凝温度界限
		maxTcHigh,# 蒸发温度上限
		TCompressorIn,# 冷凝温度下限
		maxTcHigh,# 冷凝温度上限
		refWater,# 工质
		maxCOP,# 最大COP
		eta_s,# 绝热效率
		1.0,# 插值步长
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
		dT = 1.0,# 插值步长
	)
	COP1_design=COPOverlap(TWaste, Tuse)
	COPWater_design = COPWater(TCompressorIn,Tuse)
	COPLow_design = COPLow(TWaste,TCompressorIn+dT_EvaporationStandard)
	#COP1_design=1.0
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

@time minCostTest2, minTsListTest2, P1ListTest2, P2ListTest2, P3ListTest2, PeListTest2 = generateAndSolve(PressedWaterOneStorageOneCompressor(), MinimizeCost(), VaryLoadVaryArea(), GoldenRatioMethod();
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
	PWaterCompressorMax = PWaterCompressorMax,#水蒸气压缩机最大功率
	Tsmin=120.0,
	# 求解参数
	dT = dT/10,# 状态参数高温蓄热温度离散步长
	dt = dt,# 时间步长
	smoother=0.001
)

minCostTest2-=0.001*(sum(abs2.(P1ListTest2))+sum(abs2.(P2ListTest2))+sum(abs2.(P3ListTest2))+sum(abs2.(PeListTest2)))

using Plots, DataFrames, CSV

tList = 0:dt:(24-dt)

# 看看不同的初始温度下运行费用的变化

# 看看最优运行费用下蓄热温度的变化
plt = plot(tList .+ 0.5, [P1ListTest P2ListTest P3ListTest PeListTest], label = ["P1" "P2" "P3" "Pe"],ylims=(0,2.5))
plot!(plt, tList, minTsListTest[1:end-1] / 220, label = "Ts")

plt2 = plot(tList .+ 0.5, [P1ListTest2 P2ListTest2 P3ListTest2 PeListTest2], label = ["P1" "P2" "P3" "Pe"],ylims=(0,2.5))
plot!(plt2, tList, minTsListTest2[1:end-1] / 220, label = "Ts")
