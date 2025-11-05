

# 定义参数结构体
struct SystemParameters	# 系统常量
    ThMax::Real          # 热泵最高温度
    Tuse::Real           # 使用温度
    dT::Real             # 温度差
    TCompressorIn::Real  # 压缩机入口温度
    cpm::Real            # 热容质量乘积,提出来
    COPWater::Function      # 水COP计算函数
    PhMax::Real          # 热功率最大值,提出来
    PeMax::Real          # 电功率最大值,提出来
    cp_cw::Real          # 液态水比热容
    latentHeat::Real     # 汽化潜热
    Tsmin::Real          # 蓄热最低温度
    Tsmax::Real          # 蓄热最高温度
    dTRecycleSupply::Real   # 蓄热到用热热回收最小温差
    dTRecycleBackward::Real # 用热到蓄热热回收最小温差
    sysStruct::SystemStructure  # 系统结构
    
    SystemParameters(;
        ThMax::Real,Tuse::Real,dT::Real,TCompressorIn::Real,cpm::Real,COPWater::Function,PhMax::Real,PeMax::Real,cp_cw::Real,latentHeat::Real,Tsmin::Real,Tsmax::Real,dTRecycleSupply::Real,dTRecycleBackward::Real,sysStruct::SystemStructure
    )=new(ThMax,Tuse,dT,TCompressorIn,cpm,COPWater,PhMax,PeMax,cp_cw,latentHeat,Tsmin,Tsmax,dTRecycleSupply,dTRecycleBackward,sysStruct)
end

macro unpackParameters(param_name)
    # 生成展开语句
    ex = quote
        ThMax = $(param_name).ThMax
        Tuse = $(param_name).Tuse
        dT = $(param_name).dT
        TCompressorIn = $(param_name).TCompressorIn
        cpm = $(param_name).cpm
        COPWater = $(param_name).COPWater
        PhMax = $(param_name).PhMax
        PeMax = $(param_name).PeMax
        cp_cw = $(param_name).cp_cw
        latentHeat = $(param_name).latentHeat
        Tsmin = $(param_name).Tsmin
        Tsmax = $(param_name).Tsmax
        dTRecycleSupply = $(param_name).dTRecycleSupply
        dTRecycleBackward = $(param_name).dTRecycleBackward
        sysStruct = $(param_name).sysStruct
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
    # delta[2] mode 3 can't work
    delta=[TsMid+params.dT>=params.Tuse,TsMid+params.dT>=params.ThMax]
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

function getLPModel()
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    #各个模式的功率是P[1:3],各环节电加热功率是P[4:5]
    @variable(model, P[1:5]>=0)
    return model
end

const model_lp = getLPModel()

function getMinimumCost(TsStart::Real,TsEnd::Real,dt::Real,params::SystemParameters,sysVariables::SystemVariables)
    
    ThMax = params.ThMax
    Tuse = params.Tuse
        
    if TsStart<Tuse-params.dT<TsEnd
        dt1=dt*(Tuse-params.dT-TsStart)/(TsEnd-TsStart)
        dt2=dt-dt1
        C1, flag1, P11, P21, P31,Pe1=getMinimumCost(TsStart,Tuse-params.dT,dt1,params,sysVariables)
        C2, flag2, P12, P22, P32,Pe2=getMinimumCost(Tuse-params.dT,TsEnd,dt2,params,sysVariables)
        return C1+C2, flag1&flag2, (P11*dt1+P12*dt2)/dt, (P21*dt1+P22*dt2)/dt, (P31*dt1+P32*dt2)/dt, (Pe1*dt1+Pe2*dt2)/dt 
    end

    if TsStart>Tuse+params.dT>TsEnd
        dt1=dt*(Tuse+params.dT-TsStart)/(TsEnd-TsStart)
        dt2=dt-dt1
        C1, flag1, P11, P21, P31,Pe1=getMinimumCost(TsStart,Tuse+params.dT,dt1,params,sysVariables)
        C2, flag2, P12, P22, P32,Pe2=getMinimumCost(Tuse+params.dT,TsEnd,dt2,params,sysVariables)
        return C1+C2, flag1&flag2, (P11*dt1+P12*dt2)/dt, (P21*dt1+P22*dt2)/dt, (P31*dt1+P32*dt2)/dt, (Pe1*dt1+Pe2*dt2)/dt 
    end

    if (TsStart+params.dT-ThMax)*(TsEnd+params.dT-ThMax)<0
        dt1=dt*(ThMax-params.dT-TsStart)/(TsEnd-TsStart)
        dt2=dt-dt1
        C1, flag1, P11, P21, P31,Pe1=getMinimumCost(TsStart,ThMax-params.dT,dt1,params,sysVariables)
        C2, flag2, P12, P22, P32,Pe2=getMinimumCost(ThMax-params.dT,TsEnd,dt2,params,sysVariables)
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

    recycle=[
        sysStruct.heatpumpWithStorage,
	    sysStruct.toUserRecycle,
	    sysStruct.toStorageRecycle
    ]
    Ts=0.5*(TsStart+TsEnd)
    recycleValid = [Ts + dT - Tuse >= dTRecycleSupply, Tuse - Ts - dT >= dTRecycleBackward]

    PsView=cpm/dt*(TsEnd-TsStart)# 蓄热罐温度变化示数功率（不是真实的）
    for (x1,x2,x3) in allowedStatus
        if x1 == 1 && x2 == 1 && recycle[1] ==0
            #不允许同时运行
            continue
        end
        flag,coph1,coph2,coph3,copoverlap = getCOPbyMode(x1,x2,x3,TsStart,TsEnd,params,sysVariables)

        if !flag
            continue
        end
        
        #=
        model = Model(HiGHS.Optimizer)
        set_silent(model)

        @variable(model, P[1:3]>=0)#各个模式的功率
        @variable(model, Pe[1:2]>=0)#供热电加热功率与储热电加热功率

        set_start_value(P[1], PsView*x1*coph1/copoverlap)
        
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

        model = copy(model_lp)
        set_optimizer(model, () -> HiGHS.Optimizer())
        set_silent(model)

        if x1 == 0
            @constraint(model, model[:P][1] == 0.0)  # 显式固定为0，避免任意值
        end
        if x2 == 0
            @constraint(model, model[:P][2] == 0.0)
        end
        if x3 == 0
            @constraint(model, model[:P][3] == 0.0)
        end

        A=[
            x1*coph1 x2*coph2 recycle[2]*recycleValid[1]*cp_cw/latentHeat*(Ts+dT-Tuse) 1.0 0.0;
            0.0 -x2*(coph2-1) x3*coph3 0.0 1.0;
            0.0 0.0 0.0 -1.0 -1.0;
            -1.0 0.0 -1.0 0.0 0.0
        ]
        b = [
            load,
            PsView - recycle[3]*recycleValid[2]*cp_cw/latentHeat*load*(Tuse - Ts - dT),
            -PeMax,
            -PhMax
        ]
        c = [
            coph1/copoverlap,
            1.0,
            coph3/copoverlap,
            1.0,
            1.0
        ]

        #=
        # 测试约束
        @constraint(model,P[1]*x1*coph1+P[2]*x2*coph2+Pe[1]+recycle[2]*recycleValid[1]*cp_cw/latentHeat*(Ts+dT-Tuse)*P[3]*x3*coph3>=load)
        @constraint(model,P[3]*x3*coph3-P[2]*x2*(coph2-1)+Pe[2]+recycle[3]*recycleValid[2]*cp_cw/latentHeat*load*(Tuse - Ts - dT)>=PsView)

        # 变量范围约束
        @constraint(model,sum(Pe[i] for i in 1:2)<=PeMax)
        @constraint(model,P[1]+P[3]<=PhMax)
        =#

        # 向量化约束
        @constraint(model,A*model[:P] - b .>= 0)
        # 目标函数
        #@objective(model, Min, P[1]*coph1/copoverlap+P[2]+P[3]*coph3/copoverlap+Pe[1]+Pe[2])
        @objective(model, Min,sum(c.*model[:P]))

        # 求解模型
        optimize!(model)

        # 查询求解结果，记录更优值
        if primal_status(model) in [FEASIBLE_POINT,NEARLY_FEASIBLE_POINT]
            if objective_value(model)*dt < cost
                cost = objective_value(model)*dt
                flagAll = flag
                P1Value = value(model[:P][1])*coph1/copoverlap
                P2Value = value(model[:P][2])
                P3Value = value(model[:P][3])*coph3/copoverlap
                PeValue = value(model[:P][4]+model[:P][5])
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
        
    if TsStart<Tuse-params.dT<TsEnd
        dt1=dt*(Tuse-params.dT-TsStart)/(TsEnd-TsStart)
        dt2=dt-dt1
        C1, flag1, P11, P21, P31,Pe1=getMinimumCost_MILP(TsStart,Tuse-params.dT,dt1,params,sysVariables)
        C2, flag2, P12, P22, P32,Pe2=getMinimumCost_MILP(Tuse-params.dT,TsEnd,dt2,params,sysVariables)
        return C1+C2, flag1&flag2, (P11*dt1+P12*dt2)/dt, (P21*dt1+P22*dt2)/dt, (P31*dt1+P32*dt2)/dt, (Pe1*dt1+Pe2*dt2)/dt 
    end

    if TsStart>Tuse+params.dT>TsEnd
        dt1=dt*(Tuse+params.dT-TsStart)/(TsEnd-TsStart)
        dt2=dt-dt1
        C1, flag1, P11, P21, P31,Pe1=getMinimumCost_MILP(TsStart,Tuse+params.dT,dt1,params,sysVariables)
        C2, flag2, P12, P22, P32,Pe2=getMinimumCost_MILP(Tuse+params.dT,TsEnd,dt2,params,sysVariables)
        return C1+C2, flag1&flag2, (P11*dt1+P12*dt2)/dt, (P21*dt1+P22*dt2)/dt, (P31*dt1+P32*dt2)/dt, (Pe1*dt1+Pe2*dt2)/dt 
    end

    if (TsStart+params.dT-ThMax)*(TsEnd+params.dT-ThMax)<0
        dt1=dt*(ThMax-params.dT-TsStart)/(TsEnd-TsStart)
        dt2=dt-dt1
        C1, flag1, P11, P21, P31,Pe1=getMinimumCost_MILP(TsStart,ThMax-params.dT,dt1,params,sysVariables)
        C2, flag2, P12, P22, P32,Pe2=getMinimumCost_MILP(ThMax-params.dT,TsEnd,dt2,params,sysVariables)
        return C1+C2, flag1&flag2, (P11*dt1+P12*dt2)/dt, (P21*dt1+P22*dt2)/dt, (P31*dt1+P32*dt2)/dt, (Pe1*dt1+Pe2*dt2)/dt 
    end

    # 解构常量
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
    coph4=min(coph1,coph3)
    # set condensor temperature as the higher one in Tuse and Ts+dT 
    Ts=0.5*(TsStart+TsEnd)
    Mp=2.0


    model = Model(HiGHS.Optimizer)
    set_silent(model)   # 取消显示
    # 热回收模式
    recycle=[
        sysStruct.heatpumpWithStorage,
	    sysStruct.toUserRecycle,
	    sysStruct.toStorageRecycle
    ]

    @variable(model, x[1:3], Bin)   # 三个模式是否开启
    @variable(model, P[1:3]>=0)     # 各个模式的功率
    @variable(model, Pe[1:2]>=0)    # 供热电加热功率与储热电加热功率

    delta=[Ts>=Tuse,Ts>=ThMax]      # 系统运行状态描述

    # 热回收状态描述
    # 从蓄热到用热端，温度为Ts+dT,目标温度为Tuse,热回收温差为 dTRecycleSupply
    # 从用热到蓄热端，温度为Tuse,目标温度为Ts+dT,热回收温差为 dTRecycleBackward
    recycleValid = [Ts + dT - Tuse >= dTRecycleSupply, Tuse - Ts - dT >= dTRecycleBackward]

    @variable(model, y[1:3]>=0) #辅助变量 y[i]=x[i]*P[i]
    @variable(model,y1x3>=0)    #y1x3=x[1]*x[3]*P[1]
    @variable(model,y3x1>=0)    #y3x1=x[1]*x[3]*P[3]
    @variable(model,y1cop13>=0) #真正的模式1功率
    @variable(model,y3cop13>=0) #真正的模式3功率

    #模式1的总功率，水蒸气压缩机加低温热泵
    @variable(model,mode1Power>=0)
    @variable(model,mode3Power>=0)  #模式3的总功率

    # 状态约束
    @constraint(model, x[1]+x[2]<=1+recycle[1]*delta[1])
    @constraint(model, x[2]+x[3]<=1)
    @constraint(model, x[3]<=1-delta[2])#蓄热温度超过热泵最大温度时不用热泵储热

    # 松弛变量y[i]=x[i]*P[i]
    @constraint(model,[i=1:3], y[i]<=P[i]+Mp*(1-x[i]))
    @constraint(model,[i=1:3], y[i]>=P[i]-Mp*(1-x[i]))
    @constraint(model,[i=1:3], y[i]<=Mp*x[i])
    @constraint(model,[i=1:3], y[i]>=-Mp*x[i])

    @constraint(model, y1x3<=y[1]+Mp*(1-x[3]))
    @constraint(model, y1x3>=y[1]-Mp*(1-x[3]))
    @constraint(model, y1x3<=Mp*x[3])
    @constraint(model, y1x3>=-Mp*x[3])

    @constraint(model, y3x1<=y[3]+Mp*(1-x[1]))
    @constraint(model, y3x1>=y[3]-Mp*(1-x[1]))
    @constraint(model, y3x1<=Mp*x[1])
    @constraint(model, y3x1>=-Mp*x[1])

    @constraint(model, y1cop13==y[1]*coph1+(coph4-coph1)*y1x3)
    @constraint(model, y3cop13==y[3]*coph3+(coph4-coph3)*y3x1)

    @constraint(model, mode1Power == y1cop13/COPl+y[1]*(COPl-1)/COPl)
    @constraint(model, mode3Power == y3cop13/COPl+y[3]*(COPl-1)/COPl)

    # 测试约束
    # 负载等号约束需要处理
    @constraint(model, y1cop13+y[2]*coph2+Pe[1]+recycle[2]*recycleValid[1]*cp_cw/latentHeat*(Ts+dT-Tuse)*y3cop13>=load)
    @constraint(model, y3cop13-y[2]*(coph2-1)+Pe[2]+recycle[3]*recycleValid[2]*cp_cw/latentHeat*load*(Tuse - Ts - dT)>=cpm/dt*(TsEnd-TsStart))

    # 变量范围约束
    @constraint(model,sum(Pe[i] for i in 1:2)<=PeMax)
    @constraint(model,P[1]+P[3]<=PhMax)

    # 目标函数
    @objective(model, Min, mode1Power+P[2]+mode3Power+Pe[1]+Pe[2])

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

    return objective_value(model)*dt,flag , value(mode1Power), value(y[2]), value(mode3Power),value(Pe[1]+Pe[2])
end
