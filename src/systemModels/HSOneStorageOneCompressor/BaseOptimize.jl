

# 定义参数结构体
struct SystemParameters	# 系统常量
    ThMax::Real          # 热泵最高温度
    Tuse::Real           # 使用温度
    dT::Real             # 温度差
    TCompressorIn::Real  # 压缩机入口温度
    cpm::Real            # 热容质量乘积
    COPWater::Function      # 水COP计算函数
    PhMax::Real          # 热功率最大值
    PeMax::Real          # 电功率最大值
    #cp_cw
    
    SystemParameters(;
        ThMax::Real,Tuse::Real,dT::Real,TCompressorIn::Real,cpm::Real,COPWater::Function,PhMax::Real,PeMax::Real
    )=new(ThMax,Tuse,dT,TCompressorIn,cpm,COPWater,PhMax,PeMax)
end

macro unpackParameters(param_name)
    # 生成展开语句
    ex = quote
        dT = $(param_name).dT
        TCompressorIn = $(param_name).TCompressorIn
        cpm = $(param_name).cpm
        COPWater = $(param_name).COPWater
        PhMax = $(param_name).PhMax
        PeMax = $(param_name).PeMax
    end
    # 取消卫生宏。卫生宏为防止与空间中变量名冲突会进行重命名。esc可以取消这一过程
    return esc(ex)
end

struct SystemVariables  #系统变量
    load::Real           # 负载
    COPl::Real           # 低温COP，这一项临时有用
    Tair::Real           # 环境温度
    TWaste::Real         # 废热温度
end

macro unpackVariables(param_name)
    ex = quote
        load = $(param_name).load
        COPl = $(param_name).COPl
        Tair = $(param_name).Tair
        TWaste = $(param_name).TWaste
    end
    return esc(ex)
end

"""
Calculate the COP for different mode
Return flag and COPh1,COPh2,COPh3 and COPOverlap(leap)
"""
function getCOPbyMode(x1::Union{Int,Bool},x2::Union{Int,Bool},x3::Union{Int,Bool},TsStart::Real,TsEnd::Real,params::SystemParameters,sysVariables::SystemVariables)
    TsMid = 0.5*(TsStart+TsEnd)
    # delta[1] mode 2 can work with mode 1 and 3
    # delta[2] mode 3 can work
    delta=[TsMid+params.dT>=params.Tuse,TsMid>=params.ThMax]
    # Check if status valid by temperature
    if !(x1+x2<=1+delta[1] &&
        x2+x3<=1 &&
        x3<=1-delta[2])

        return false,1.0,1.0,1.0,1.0
    end

    if x1==0 && x2==0 && x3==0
        if TsStart == TsEnd
            return true,1.0,1.0,1.0,1.0
        else
            return false,1.0,1.0,1.0,1.0
        end
    elseif x1==1 && x2==0 && x3==0
        coph1 = params.COPWater(params.TCompressorIn,params.Tuse)
        return true,
        coph1,
        1.0,
        1.0,
        sysVariables.COPl*coph1/(coph1+sysVariables.COPl-1)
    elseif x1==0 && x2==1 && x3==0
        if TsStart > TsEnd
            return true,
            1.0,
            params.COPWater(TsMid-params.dT,params.Tuse),
            1.0,
            1.0
        else
            return false,1.0,1.0,1.0,1.0
        end
    elseif x1==1 && x2==1 && x3==0
        if TsStart > TsEnd >= params.Tuse+params.dT
            coph1 = params.COPWater(params.TCompressorIn,params.Tuse)
            return true,
            coph1,
            params.COPWater(TsMid-params.dT,params.Tuse),
            1.0,
            sysVariables.COPl*coph1/(coph1+sysVariables.COPl-1)
        else
            return false,1.0,1.0,1.0,1.0
        end
    elseif x1==0 && x2==0 && x3==1
        if TsStart < TsEnd
            coph3=params.COPWater(params.TCompressorIn,TsMid+params.dT)
            return true,
            1.0,
            1.0,
            coph3,
            sysVariables.COPl*coph3/(coph3+sysVariables.COPl-1)
        else
            return false,1.0,1.0,1.0,1.0
        end
        
    elseif x1==1 && x2==0 && x3==1
        if TsStart < TsEnd
            tempCOP=params.COPWater(params.TCompressorIn,max(params.Tuse,TsMid+params.dT))
            return true,
            tempCOP,
            1.0,
            tempCOP,
            sysVariables.COPl*tempCOP/(tempCOP+sysVariables.COPl-1)
        else
            return false,1.0,1.0,1.0,1.0
        end
    else
        @warn "error status x=[$x1,$x2,$x3]"
        return false,1.0,1.0,1.0,1.0
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

function getMinimumCost(TsStart::Real,TsEnd::Real,dt::Real,params::SystemParameters,sysVariables::SystemVariables)
    
    ThMax = params.ThMax
    Tuse = params.Tuse

    if (TsStart-Tuse)*(TsEnd-Tuse)<0
        dt1=dt*(Tuse-TsStart)/(TsEnd-TsStart)
        dt2=dt-dt1
        C1, flag1, P11, P21, P31,Pe1=getMinimumCost(TsStart,Tuse,dt1,params,sysVariables)
        C2, flag2, P12, P22, P32,Pe2=getMinimumCost(Tuse,TsEnd,dt2,params,sysVariables)
        return C1+C2, flag1&flag2, (P11*dt1+P12*dt2)/dt, (P21*dt1+P22*dt2)/dt, (P31*dt1+P32*dt2)/dt, (Pe1*dt1+Pe2*dt2)/dt 
    end

    if (TsStart-ThMax)*(TsEnd-ThMax)<0
        dt1=dt*(ThMax-TsStart)/(TsEnd-TsStart)
        dt2=dt-dt1
        C1, flag1, P11, P21, P31,Pe1=getMinimumCost(TsStart,ThMax,dt1,params,sysVariables)
        C2, flag2, P12, P22, P32,Pe2=getMinimumCost(ThMax,TsEnd,dt2,params,sysVariables)
        return C1+C2, flag1&flag2, (P11*dt1+P12*dt2)/dt, (P21*dt1+P22*dt2)/dt, (P31*dt1+P32*dt2)/dt, (Pe1*dt1+Pe2*dt2)/dt 
    end

    # 解构常量
    @unpackParameters(params)
    @unpackVariables(sysVariables)

    # 解构变量
    cost=9999.0
    flagAll=false
    P1Value=0.0
    P2Value=0.0
    P3Value=0.0
    PeValue=0.0

    PsView=cpm/dt*(TsEnd-TsStart)# 蓄热罐温度变化示数功率（不是真实的）
    for (x1,x2,x3) in allowedStatus
        flag,coph1,coph2,coph3,copoverlap = getCOPbyMode(x1,x2,x3,TsStart,TsEnd,params,sysVariables)

        if !flag
            continue
        end
        
        model = Model(HiGHS.Optimizer)
        set_silent(model)

        @variable(model, P[1:3]>=0)#各个模式的功率
        @variable(model, Pe[1:2]>=0)#供热电加热功率与储热电加热功率

        set_start_value(P[1], PsView*x1*coph1/copoverlap)
        #=
        if x1 == 0
            @constraint(model, P[1] == 0.0)  # 显式固定为0，避免任意值
        end
        if x2 == 0
            @constraint(model, P[2] == 0.0)
        end
        if x3 == 0
            @constraint(model, P[3] == 0.0)
        end
        =#

        # 测试约束
        @constraint(model,P[1]*x1*coph1+P[2]*x2*coph2+Pe[1]>=load)
        @constraint(model,P[3]*x3*coph3-P[2]*x2*(coph2-1)+Pe[2]==PsView)

        # 变量范围约束
        @constraint(model,sum(Pe[i] for i in 1:2)<=PeMax)
        @constraint(model,P[1]+P[3]<=PhMax)

        # 目标函数
        @objective(model, Min, P[1]*coph1/copoverlap+P[2]+P[3]*coph3/copoverlap+Pe[1]+Pe[2])

        # 求解模型
        optimize!(model)

        # 查询求解结果，记录更优值
        if primal_status(model) in [FEASIBLE_POINT,NEARLY_FEASIBLE_POINT]
            if objective_value(model)*dt < cost
                cost = objective_value(model)*dt
                flagAll = flag
                P1Value = value(P[1])*coph1/copoverlap
                P2Value = value(P[2])
                P3Value = value(P[3])*coph3/copoverlap
                PeValue = value(Pe[1]+Pe[2])
            end
        else
            continue
        end
    end
    if !flagAll
        #println("求解失败,TsStart=$TsStart,TsEnd=$TsEnd,dt=$dt")
        return 9999.0, false, 9999.0, 9999.0, 9999.0,9999.0
    end

    return cost, flagAll, P1Value, P2Value, P3Value, PeValue
end

"""
混合整数线性规划求解器，有问题，COP约束不对 2025.10.09
"""
function getMinimumCost_MILP(TsStart::Real,TsEnd::Real,dt::Real,params::SystemParameters,sysVariables::SystemVariables)
    
    ThMax = params.ThMax
    Tuse = params.Tuse

    if (TsStart-Tuse)*(TsEnd-Tuse)<0
        dt1=dt*(Tuse-TsStart)/(TsEnd-TsStart)
        dt2=dt-dt1
        C1, flag1, P11, P21, P31,Pe1=getMinimumCost_MILP(TsStart,Tuse,dt1,params,sysVariables::SystemVariables)
        C2, flag2, P12, P22, P32,Pe2=getMinimumCost_MILP(Tuse,TsEnd,dt2,params,sysVariables::SystemVariables)
        return C1+C2, flag1&flag2, (P11*dt1+P12*dt2)/dt, (P21*dt1+P22*dt2)/dt, (P31*dt1+P32*dt2)/dt, (Pe1*dt1+Pe2*dt2)/dt 
    end

    if (TsStart-ThMax)*(TsEnd-ThMax)<0
        dt1=dt*(ThMax-TsStart)/(TsEnd-TsStart)
        dt2=dt-dt1
        C1, flag1, P11, P21, P31,Pe1=getMinimumCost_MILP(TsStart,ThMax,dt1,params,sysVariables::SystemVariables)
        C2, flag2, P12, P22, P32,Pe2=getMinimumCost_MILP(ThMax,TsEnd,dt2,params,sysVariables::SystemVariables)
        return C1+C2, flag1&flag2, (P11*dt1+P12*dt2)/dt, (P21*dt1+P22*dt2)/dt, (P31*dt1+P32*dt2)/dt, (Pe1*dt1+Pe2*dt2)/dt 
    end

    # 解构参数以便使用
    @unpackParameters(params)
    @unpackVariables(sysVariables)
    #=
    load = sysVariables.load
    dT = params.dT
    TCompressorIn = params.TCompressorIn
    cpm = params.cpm
    COPWater = params.COPWater
    COPl = sysVariables.COPl
    PhMax = params.PhMax
    PeMax = params.PeMax
    =#

    Ts=0.5*(TsStart+TsEnd)

    coph1=COPWater(TCompressorIn,Tuse)# mode 1
    coph2=COPWater(Ts-dT,Tuse)# mode 2
    coph3=COPWater(TCompressorIn,Ts+dT)# mode 3
    #coph4=COPWater(TCompressorIn,max(Tuse,Ts+dT))
    # set condensor temperature as the higher one in Tuse and Ts+dT 
    coph1=coph3=min(coph1,coph3)
    copoverlap=coph1*COPl/(coph1+COPl-1)

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
    # 负载等号约束需要处理
    @constraint(model,y[1]*coph1+y[2]*coph2+Pe[1]==load)
    @constraint(model,y[3]*coph3-y[2]*(coph2-1)+Pe[2]==cpm/dt*(TsEnd-TsStart))

    # 变量范围约束
    @constraint(model,sum(Pe[i] for i in 1:2)<=PeMax)
    @constraint(model,P[1]+P[3]<=PhMax)

    # 目标函数
    @objective(model, Min, P[1]*coph1/copoverlap+P[2]+P[3]*coph3/copoverlap+Pe[1]+Pe[2])

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
        #println("求解失败,$isFeasible,TsStart=$TsStart,TsEnd=$TsEnd,dt=$dt")
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

    return objective_value(model)*dt,flag , value(y[1])*coph1/COPl, value(y[2]), value(y[3])*coph3/COPl,value(Pe[1]+Pe[2])
end
