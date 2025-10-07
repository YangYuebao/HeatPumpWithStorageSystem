using Revise
using Pkg
using HeatPumpWithStorageSystem
using DataFrames, CSV
using JuMP, HiGHS

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
	#load=heatLoad,
    load=0,
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


function getMinimumCost_test(TsStart,TsEnd,dt,params)
    @time begin
    ThMax = params.ThMax
    Tuse = params.Tuse

    if (TsStart-Tuse)*(TsEnd-Tuse)<0
        dt1=dt*(Tuse-TsStart)/(TsEnd-TsStart)
        dt2=dt-dt1
        COPl,coph1,coph2,coph3,coph4,C1, flag1, P11, P21, P31,Pe1=getMinimumCost(TsStart,Tuse,dt1,params)
        COPl,coph1,coph2,coph3,coph4,C2, flag2, P12, P22, P32,Pe2=getMinimumCost(Tuse,TsEnd,dt2,params)
        return COPl,coph1,coph2,coph3,coph4,C1+C2, flag1&flag2, (P11*dt1+P12*dt2)/dt, (P21*dt1+P22*dt2)/dt, (P31*dt1+P32*dt2)/dt, (Pe1*dt1+Pe2*dt2)/dt 
    end

    if (TsStart-ThMax)*(TsEnd-ThMax)<0
        dt1=dt*(ThMax-TsStart)/(TsEnd-TsStart)
        dt2=dt-dt1
        COPl,coph1,coph2,coph3,coph4,C1, flag1, P11, P21, P31,Pe1=getMinimumCost(TsStart,ThMax,dt1,params)
        COPl,coph1,coph2,coph3,coph4,C2, flag2, P12, P22, P32,Pe2=getMinimumCost(ThMax,TsEnd,dt2,params)
        return COPl,coph1,coph2,coph3,coph4,C1+C2, flag1&flag2, (P11*dt1+P12*dt2)/dt, (P21*dt1+P22*dt2)/dt, (P31*dt1+P32*dt2)/dt, (Pe1*dt1+Pe2*dt2)/dt 
    end

     # 解构参数以便使用
    
    load = params.load
    dT = params.dT
    TCompressorIn = params.TCompressorIn
    cpm = params.cpm
    COPWater = params.COPWater
    COPl = params.COPl
    PhMax = params.PhMax
    PeMax = params.PeMax

    Ts=0.5*(TsStart+TsEnd)

    coph1=COPWater(TCompressorIn,Tuse)
    coph2=COPWater(Ts-dT,Tuse)
    coph3=COPWater(Tuse,Ts+dT)
    coph4=COPWater(TCompressorIn,max(Tuse,Ts+dT))

    coph1=coph3=min(coph1,coph3)

    Ts=0.5*(TsStart+TsEnd)
    Mp=2.0
    Mt=150.0

    
    model = Model(HiGHS.Optimizer)

    # 取消显示
    set_silent(model)

    @variable(model, x[1:3], Bin)#三个模式是否开启

    #=
    set_start_value(x[1], 1)
    set_start_value(x[2], 0)
    set_start_value(x[3], 0)
    =#

    @variable(model, P[1:3]>=0)#各个模式的功率
    @variable(model, Pe[1:2]>=0)#供热电加热功率与储热电加热功率

    #=
    @variable(model, delta[1:2], Bin)#表示蓄热温度与用热温度、切换温度的关系

    set_start_value(delta[1], Ts>=Tuse)
    set_start_value(delta[2], Ts>=ThMax)
    =#

    delta=[Ts>=Tuse,Ts>=ThMax]

    @variable(model, y[1:3]>=0)#辅助变量

    # 状态约束
    @constraint(model, x[1]+x[2]<=1+delta[1])
    @constraint(model, x[2]+x[3]<=1+delta[1])
    #@constraint(model, Mt*delta[1]>=Ts-Tuse-dT)#蓄热温度与用热温度关系
    #@constraint(model, Mt*delta[1]<=Tuse-Ts-dT+Mt)
    #@constraint(model, Mt*delta[2]>=Ts-ThMax)#蓄热温度与热泵最高温度关系
    #@constraint(model, Mt*delta[2]<=Ts-ThMax+Mt)
    @constraint(model, x[3]<=1-delta[2])#蓄热温度超过热泵最大温度时不用热泵储热

    # 松弛变量y[i]=x[i]*P[i]
    @constraint(model,[i=1:3], y[i]<=P[i]+Mp*(1-x[i]))
    @constraint(model,[i=1:3], y[i]>=P[i]-Mp*(1-x[i]))
    @constraint(model,[i=1:3], y[i]<=Mp*x[i])
    @constraint(model,[i=1:3], y[i]>=-Mp*x[i])

    # 测试约束
    @constraint(model,y[1]*coph1+y[2]*coph2+Pe[1]>=load)
    @constraint(model,y[3]*coph3-y[2]*(coph2-1)+Pe[2]==cpm/dt*(TsEnd-TsStart))

    # 变量范围约束
    @constraint(model,sum(Pe[i] for i in 1:2)<=PeMax)
    @constraint(model,P[1]+P[3]<=PhMax)

    # 目标函数
    @objective(model, Min, P[1]*coph1/COPl+P[2]+P[3]*coph3/COPl+Pe[1]+Pe[2])
    end
    # 求解模型
    println("opt!")
    @time optimize!(model)

    # 查询求解结果
    isFeasible = primal_status(model)
    
    flag = isFeasible in [FEASIBLE_POINT,NEARLY_FEASIBLE_POINT]

    if !flag
        println("求解失败,$isFeasible,TsStart=$TsStart,TsEnd=$TsEnd,dt=$dt")
        return COPl,coph1,coph2,coph3,coph4,9999.0,flag,9999.0,9999.0,9999.0,9999.0
    end

    return COPl,coph1,coph2,coph3,coph4,objective_value(model)*dt,flag , value(P[1])*coph1/COPl, value(P[2]), value(P[3])*coph3/COPl,value(Pe[1]+Pe[2])
end


TsStart=120.0
TsEnd=120.0
dt=2.0

COPl=params.COPl


@time result = getMinimumCost_test(TsStart,TsEnd,dt,params)


"""
COPl,coph1,coph2,coph3,coph4,
"""
title=[
	"COPl","coph1","coph2","coph3","coph4","power_consumption","flag" , "p1", "p2", "p3","pe","cpm","dT","dt"
]

df=DataFrame()
for t in title
	df[!,t]=[]
end

push!(df,vcat(result...,cpm,TsEnd-TsStart,dt))


CSV.write(joinpath(pwd(),"calculations","situation14","result.csv"),df)

