
# 定义参数结构体
struct SystemParameters
    ThMax::Real          # 热泵最高温度
    Tuse::Real           # 使用温度
    load::Real           # 负载
    dT::Real             # 温度差
    TCompressorIn::Real  # 压缩机入口温度
    cpm::Real            # 热容质量乘积
    COPWater::Function      # 水COP计算函数
    COPl::Real           # 其他COP值
    PhMax::Real          # 热功率最大值
    PeMax::Real          # 电功率最大值,
    Tair::Real           # 环境温度
    TWaste::Real         # 废热温度
    SystemParameters(;
        ThMax::Real,Tuse::Real,load::Real,dT::Real,TCompressorIn::Real,cpm::Real,COPWater::Function,COPl::Real,PhMax::Real,PeMax::Real,Tair::Real,TWaste::Real
    )=new(ThMax,Tuse,load,dT,TCompressorIn,cpm,COPWater,COPl,PhMax,PeMax,Tair,TWaste)
end

"""
Is status valid by temperature
"""
function isValidTemperature(x1,x2,x3,)

end

"""
Calculate the COP for different mode
Return flag and COPh1,COPh2,COPh3
"""
function getCOPbyMode(x1::Int,x2::Int,x3::Int,TsStart::Real,TsEnd::Real,params::SystemParameters)
    TsMid = 0.5*(TsStart+TsEnd)
    # delta[1] mode 2 can work with mode 1 and 3
    # delta[2] mode 3 can work
    delta=[TsMid+params.dT>=params.Tuse,TsMid>=params.ThMax]
    # Check if status valid by temperature
    if !(x1+x2<=1+delta[1] &&
        x2+x3<=1 &&
        x3<=1-delta[2])

        return false,1.0,1.0,1.0
    end

    #=
    status=[
        (0,0,0),
        (1,0,0),
        (0,1,0),
        (1,1,0),
        (0,0,1),
        (1,0,1)
    ]
    =#
    if x1==0 && x2==0 && x3==0
        return true,1.0,1.0,1.0
    elseif x1==1 && x2==0 && x3==0
        return true,
        params.COPWater(params.TCompressorIn,params.Tuse),
        1.0,
        1.0
    elseif x1==0 && x2==1 && x3==0
        return true,
        1.0,
        params.COPWater(TsMid-params.dT,params.Tuse),
        1.0
    elseif x1==1 && x2==1 && x3==0
        return true,
        params.COPWater(params.TCompressorIn,params.Tuse),
        params.COPWater(TsMid-params.dT,params.Tuse),
        1.0
    elseif x1==0 && x2==0 && x3==1
        return true,
        1.0,
        1.0,
        params.COPWater(params.TCompressorIn,TsMid+params.dT)
    elseif x1==1 && x2==0 && x3==1
        tempCOP=params.COPWater(params.TCompressorIn,max(params.Tuse,TsMid+params.dT))
        return true,
        tempCOP,
        1.0,
        tempCOP
    elseif x1==0 && x2==1 && x3==1
        # invalid
        return false,1.0,1.0,1.0
    elseif x1==1 && x2==1 && x3==1
        # invalid
        return false,1.0,1.0,1.0
    else
        @warn "error status x=[$x1,$x2,$x3]"
        return false,1.0,1.0,1.0
    end
end

const allowedStatus = [
    (0,0,0),
    (1,0,0),
    (0,1,0),
    (1,1,0),
    (0,0,1),
    (1,0,1)
]

function getMinimumCost(TsStart::Real,TsEnd::Real,dt::Real,params::SystemParameters)
    
    ThMax = params.ThMax
    Tuse = params.Tuse

    if (TsStart-Tuse)*(TsEnd-Tuse)<0
        dt1=dt*(Tuse-TsStart)/(TsEnd-TsStart)
        dt2=dt-dt1
        C1, flag1, P11, P21, P31,Pe1=getMinimumCost(TsStart,Tuse,dt1,params)
        C2, flag2, P12, P22, P32,Pe2=getMinimumCost(Tuse,TsEnd,dt2,params)
        return C1+C2, flag1&flag2, (P11*dt1+P12*dt2)/dt, (P21*dt1+P22*dt2)/dt, (P31*dt1+P32*dt2)/dt, (Pe1*dt1+Pe2*dt2)/dt 
    end

    if (TsStart-ThMax)*(TsEnd-ThMax)<0
        dt1=dt*(ThMax-TsStart)/(TsEnd-TsStart)
        dt2=dt-dt1
        C1, flag1, P11, P21, P31,Pe1=getMinimumCost(TsStart,ThMax,dt1,params)
        C2, flag2, P12, P22, P32,Pe2=getMinimumCost(ThMax,TsEnd,dt2,params)
        return C1+C2, flag1&flag2, (P11*dt1+P12*dt2)/dt, (P21*dt1+P22*dt2)/dt, (P31*dt1+P32*dt2)/dt, (Pe1*dt1+Pe2*dt2)/dt 
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

    for (x1,x2,x3) in allowedStatus
        flag,coph1,coph2,coph3 = getCOPbyMode(x1,x2,x3,TsStart,TsEnd,params)
    
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    @variable(model, P[1:3]>=0)#各个模式的功率
    @variable(model, Pe[1:2]>=0)#供热电加热功率与储热电加热功率

    # 测试约束
    @constraint(model,P[1]*x1*coph1+P[2]*x2*coph2+Pe[1]>=load)
    @constraint(model,P[3]*x3*coph3-P[2]*x2*(coph2-1)+Pe[2]==cpm/dt*(TsEnd-TsStart))

    # 变量范围约束
    @constraint(model,sum(Pe[i] for i in 1:2)<=PeMax)
    @constraint(model,P[1]+P[3]<=PhMax)

    # 目标函数
    @objective(model, Min, P[1]*coph1/COPl+P[2]+P[3]*coph3/COPl+Pe[1]+Pe[2])

    # 求解模型
    optimize!(model)

    # 查询求解结果
    isFeasible = primal_status(model)
    
    flag = isFeasible in [FEASIBLE_POINT,NEARLY_FEASIBLE_POINT]

    #=
    if params.PhMax==0
        println("params.PhMax==0")
    end
    =#
    if !flag
        println("求解失败,$isFeasible,TsStart=$TsStart,TsEnd=$TsEnd,dt=$dt")
        #=
        ThMax
        Tuse::Real           # 使用温度
        load::Real           # 负载
        dT::Real             # 温度差
        TCompressorIn::Real  # 压缩机入口温度
        cpm::Real            # 热容质量乘积
        COPWater::Function      # 水COP计算函数
        COPl::Real           # 其他COP值
        PhMax::Real          # 热功率最大值
        PeMax::Real          # 电功率最大值,
        Tair::Real           # 环境温度
        TWaste::Real         # 废热温度
        =#
        if TsStart==TsEnd
            println("TsStart==TsEnd. Parametets:")
            println("""
                ThMax=$ThMax,
                Tuse=$Tuse,
                load=$load,
                dT=$dT,
                TCompressorIn=$TCompressorIn,
                cpm=$cpm,
                COPl=$COPl,
                PhMax=$PhMax,
                PeMax=$PeMax,
                Tair=$(params.Tair),
                TWaste=$(params.TWaste)
            """)
        end
        return 9999.0,flag,9999.0,9999.0,9999.0,9999.0
    end

    return objective_value(model)*dt,flag , value(P[1])*coph1/COPl, value(P[2]), value(P[3])*coph3/COPl,value(Pe[1]+Pe[2])
end

function getMinimumCost_MILP(TsStart::Real,TsEnd::Real,dt::Real,params::SystemParameters)
    
    ThMax = params.ThMax
    Tuse = params.Tuse

    if (TsStart-Tuse)*(TsEnd-Tuse)<0
        dt1=dt*(Tuse-TsStart)/(TsEnd-TsStart)
        dt2=dt-dt1
        C1, flag1, P11, P21, P31,Pe1=getMinimumCost_MILP(TsStart,Tuse,dt1,params)
        C2, flag2, P12, P22, P32,Pe2=getMinimumCost_MILP(Tuse,TsEnd,dt2,params)
        return C1+C2, flag1&flag2, (P11*dt1+P12*dt2)/dt, (P21*dt1+P22*dt2)/dt, (P31*dt1+P32*dt2)/dt, (Pe1*dt1+Pe2*dt2)/dt 
    end

    if (TsStart-ThMax)*(TsEnd-ThMax)<0
        dt1=dt*(ThMax-TsStart)/(TsEnd-TsStart)
        dt2=dt-dt1
        C1, flag1, P11, P21, P31,Pe1=getMinimumCost_MILP(TsStart,ThMax,dt1,params)
        C2, flag2, P12, P22, P32,Pe2=getMinimumCost_MILP(ThMax,TsEnd,dt2,params)
        return C1+C2, flag1&flag2, (P11*dt1+P12*dt2)/dt, (P21*dt1+P22*dt2)/dt, (P31*dt1+P32*dt2)/dt, (Pe1*dt1+Pe2*dt2)/dt 
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

    coph1=COPWater(TCompressorIn,Tuse)# mode 1
    coph2=COPWater(Ts-dT,Tuse)# mode 2
    coph3=COPWater(TCompressorIn,Ts+dT)# mode 3
    #coph4=COPWater(TCompressorIn,max(Tuse,Ts+dT))

    println("COPh1=$coph1 COPh3=$coph3")
    # set condensor temperature as the higher one in Tuse and Ts+dT 
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
    @constraint(model, x[2]+x[3]<=1)
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

    # 求解模型
    optimize!(model)

    # 查询求解结果
    isFeasible = primal_status(model)
    
    flag = isFeasible in [FEASIBLE_POINT,NEARLY_FEASIBLE_POINT]

    #=
    if params.PhMax==0
        println("params.PhMax==0")
    end
    =#
    if !flag
        println("求解失败,$isFeasible,TsStart=$TsStart,TsEnd=$TsEnd,dt=$dt")
        #=
        ThMax
        Tuse::Real           # 使用温度
        load::Real           # 负载
        dT::Real             # 温度差
        TCompressorIn::Real  # 压缩机入口温度
        cpm::Real            # 热容质量乘积
        COPWater::Function      # 水COP计算函数
        COPl::Real           # 其他COP值
        PhMax::Real          # 热功率最大值
        PeMax::Real          # 电功率最大值,
        Tair::Real           # 环境温度
        TWaste::Real         # 废热温度
        =#
        if TsStart==TsEnd
            println("TsStart==TsEnd. Parametets:")
            println("""
                ThMax=$ThMax,
                Tuse=$Tuse,
                load=$load,
                dT=$dT,
                TCompressorIn=$TCompressorIn,
                cpm=$cpm,
                COPl=$COPl,
                PhMax=$PhMax,
                PeMax=$PeMax,
                Tair=$(params.Tair),
                TWaste=$(params.TWaste)
            """)
        end
        return 9999.0,flag,9999.0,9999.0,9999.0,9999.0
    end

    return objective_value(model)*dt,flag , value(P[1])*coph1/COPl, value(P[2]), value(P[3])*coph3/COPl,value(Pe[1]+Pe[2])
end
