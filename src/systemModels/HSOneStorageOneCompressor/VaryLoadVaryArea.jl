#=
变工况与定工况的区别为：热泵配置可能小于标准配置，原先热泵的配置大小足够，蓄热在放热时因为效率总是高于标准工况所以不用考虑需要电热补热；现在需要考虑电热补热。因此：
1.需要考虑热泵供热不足时用电热补热
2.需要给定压缩机的运行功率约束
=#

"""
单个时间层的状态转移成本函数
"""
function getStateTransitionCost_SingleStep(::PressedWaterOneStorageOneCompressor;
	COPOverlap::Function,
	COPWater::Function,
	heatLoad::Real,#热负荷
	Tair::Real,# 外部环境温度
	# 总循环参数
	Tuse::Real,# 供热蒸汽温度
	dT_EvaporationStandard::Real,#全蒸温差
	latentHeat::Real,# 汽化潜热
	cp_cw::Real,# 循环水定压热容
	TcChangeToElec::Real,
	TWaste::Real,# 废热回收蒸发器温度

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容

	# 设备运行约束
	PheatPumpMax::Real,# 热泵最大功率
    PelecHeatMax::Real,# 电锅炉最大功率
    PWaterCompressorMax::Real,# 水蒸气压缩机的最大功率
    Tsmin::Real,# 蓄热的最小温度
    
	# 求解参数
	TsListStart::Matrix,# 状态参数高温蓄热温度起始值列表
	TsListEnd::Matrix,# 状态参数高温蓄热温度结束值列表
	dt::Real = 0.1,# 时间步长
)
	COP1 = COPOverlap(TWaste, Tuse)
	COP2 = COPWater
	COP3 = COPOverlap

	nT = length(TsListStart)# 温度步数

	#C[i,j]表示温度从TsList[i]到TsList[j]时的最低功率；在负载、环境不变的情况下，C[i,j]是不变的
	C = fill(99.0, nT, nT)# 状态转移矩阵
	P1Matrix = zeros(nT, nT)# 状态转移功率1参数
	P2Matrix = zeros(nT, nT)# 状态转移功率2参数
	P3Matrix = zeros(nT, nT)# 状态转移功率3参数
	PeMatrix = zeros(nT, nT)# 状态转移功率电加热参数

	# 只用热泵供热时的功率
	P1Only = heatLoad / COP1# 只用热泵供热时的功率

	"""计算不同蓄热温度下蓄热温度降低的电度"""
	function powerCalculate_disCharge(Tsaim, Tsstart, dt)
		COP2value=COP2((TsList[i] + TsList[j]) / 2 - dT_EvaporationStandard, Tuse)
		Psout=cpm_h*(Tsstart-Tsaim)/dt
		P2=Psout / (COP2value-1)
		if (P2+Psout>heatLoad) || (P2 > PWaterCompressorMax) # 降不到这个温度
			return 9999.0,false,9999.0,9999.0,9999.0
		end
        #先看是不是温度到下限了
        if Tsaim == Tsmin
            P1 = min((heatLoad-P2-Psout) / COP1,PWaterCompressorMax)
            Pe = heatLoad-P2-Psout-P1 * COP1
            P21=min(PheatPumpMax,heatLoad/COP2value)  # 前半段蓄热供热最大P2功率
            Pe1=heatLoad-P21*COP2value
            P12 = min(heatLoad / COP1,PWaterCompressorMax)   # 后半段热泵供热的功率
			Pe2 = heatLoad-P12 * COP1
            if Pe1 <= PheatPumpMax && Pe2 <= PheatPumpMax
			    return (P1+P2+Pe)*dt,true,P1,P2,Pe
            else
                return 9999.0,false,9999.0,9999.0,9999.0
            end
        end

        #蓄热温度的末态高于用热温度,放热全程热泵补热;或放热到蓄热的最低温度,先放热后热泵
		if (Tsaim-dT_EvaporationStandard>=Tuse)
			P1 = min((heatLoad-P2-Psout) / COP1,PWaterCompressorMax)
			Pe = heatLoad-P2-Psout-P1 * COP1
            if Pe <= PheatPumpMax
			    return (P1+P2+Pe)*dt,true,P1,P2,Pe
            else
                return 9999.0,false,9999.0,9999.0,9999.0
            end
		elseif Tsstart-dT_EvaporationStandard<=Tuse#蓄热温度的初态低于用热温度，放热全程电热补热
			P1 = 0
			Pe = heatLoad-P2-Psout
            if Pe <= PheatPumpMax
			    return (Pe+P2)*dt,true,P1,P2,Pe
            else
                return 9999.0,false,9999.0,9999.0,9999.0
            end
		else# 蓄热温度的初态在用热温度中间，放热部分热泵补热，部分电热补热
			tmid=dt*(Tsstart - Tuse - dT_EvaporationStandard)/(Tsstart-Tsaim)
			C1,flag1,P1,_,Pe1=powerCalculate_disCharge(Tuse+dT_EvaporationStandard, Tsstart, tmid)
			C2,flag2,_,_,Pe2=powerCalculate_disCharge(Tsaim, Tuse+dT_EvaporationStandard, dt-tmid)
			P1*=tmid/dt
			Pe=Pe1*tmid/dt+Pe2*(dt-tmid)/dt
			if flag1 && flag2
				return C1+C2,true,P1,P2,Pe
			else
				return 9999.0,false,9999.0,9999.0,9999.0
			end
		end
	end

	"""计算蓄热温度相对较低时升温的电度情况"""
	function powerCalculate_lowTs(Tsaim, Tsstart, dt)
		P3 = min(PheatPumpMax - P1Only, cpm_h / dt / COP1 * (Tsaim - Tsstart))
		Pe = cpm_h / dt * (Tsaim - Tsstart) - P3 * COP1
		flag = Pe <= PelecHeatMax
		return (P1Only + P3 + Pe) * dt, flag, P1Only, P3, Pe
	end

	"""计算蓄热温度相对较高时升温的电度情况(考虑蓄热温度超过热泵温度上限)"""
	function powerCalculate_highTs(Tsaim, Tsstart, dt)
		COP3value = COP3(TWaste, (Tsaim + Tsstart) / 2 + dT_EvaporationStandard)

		function mode1(dt1)# 模式1：电热加热
			Pe = cpm_h / dt1 * (Tsaim - Tsstart)
			if Pe <= PelecHeatMax
				return true, P1Only, 0.0, Pe
			else
				return false, 9999.0, 9999.0, 9999.0
			end
			
		end

		function mode2(Tsaim2, Tsstart2, dt2)
			P3 = min(
				PheatPumpMax,
				(PheatPumpMax - heatLoad / COP3value) / (1 - cp_cw / latentHeat * (Tsaim2 - Tsstart2)),
				cpm_h / dt2 / COP3value * (Tsaim2 - Tsstart2),
			)
			P1 = max(0, heatLoad / COP3value - cp_cw / latentHeat * (Tsaim2 - Tsstart2) * P3)
			Pe = cpm_h / dt2 * (Tsaim2 - Tsstart2) - P3 * COP3value
			if Pe <= PelecHeatMax
				return true, P1, P3, Pe
			else
				return false, 9999.0, 9999.0, 9999.0
			end
		end
		#COP3value = COP3(TWaste, (min(Tsaim,TcChangeToElec - dT_EvaporationStandard) + Tsstart) / 2 + dT_EvaporationStandard)
		
		Ptotal = 999.0
		P1res = 0.0
		P3res = 0.0
		Peres = 0.0
		flag = false
		# 模式1：电热加热
		flag, P1res, P3res, Peres= mode1(dt)
		Ptotal=P1res+P3res+Peres

		# 模式2：热泵加电热,又分两种模式：末态温度不大于电加热温度界限与末态温度大于电加热温度界线
		if (Tsaim <= TcChangeToElec - dT_EvaporationStandard)	# 该时间层的末态温度大于电加热温度界限
			flag_2, P1_2, P3_2, Pe_2= mode2(Tsaim,Tsstart,dt)
			Ptotal_2=P1_2+P3_2+Pe_2
		else
			tmid=dt*(TcChangeToElec - dT_EvaporationStandard-Tsstart)/(Tsaim - Tsstart)
			flag_21, P1_21, P3_21, Pe_21= mode2(TcChangeToElec - dT_EvaporationStandard,Tsstart,tmid)
			flag_22, P1_22, P3_22, Pe_22= mode2(Tsaim,TcChangeToElec - dT_EvaporationStandard,dt-tmid)
			flag_2=flag_21 && flag_22
			P1_2=P1_21*tmid/dt+P1_22*(dt-tmid)/dt
			P3_2=P3_21*tmid/dt+P3_22*(dt-tmid)/dt
			Pe_2=Pe_21*tmid/dt+Pe_22*(dt-tmid)/dt
			Ptotal_2=P1_2+P3_2+Pe_2
		end
		if (Ptotal_2 < Ptotal)
			flag = flag || flag_2
			P1res = P1_2
			P3res = P3_2
			Peres = Pe_2
			Ptotal = Ptotal_2
		end
		return Ptotal * dt, flag, P1res, P3res, Peres
	end

	for i ∈ 1:nT
		# 先计算温度下降的功率,存在一边电加热一边开2号模式的情况
		# 首先计算该模式在下一刻的最低温度
		TsNextMin = TsList[i]
		# 简单迭代求蓄热的最低温度
		for _ ∈ 1:3
			TsNextMin = TsList[i] - heatLoad * dt / cpm_h * (1 - 1 / COP2(0.5 * (TsNextMin + TsList[i]) - dT_EvaporationStandard, Tuse))
		end
		# 计算最低温度的下标
		ibias=findfirst(x -> x >= TsNextMin, TsList)
		TsDecreaseIndexList[i]=i-ibias
		TsNextMinGrid = TsList[ibias]
		elecP = (TsNextMinGrid - TsNextMin) * cpm_h

		# 如果温度在一个时间层内会跌到低于最低温度，那么将这个时间层分为两段，前段用蓄热，后段用热泵直供
		index = i - TsDecreaseIndexList[i]
		if index == 1
			j=1
			P21 = heatLoad / (COP2((TsList[i] + TsList[j]) / 2 - dT_EvaporationStandard, Tuse))
			elecP1 = (TsList[j] - TsNextMin) * cpm_h
			#补热的功率基本上超不了，因为这个工况下从蓄热取走的热量小于工厂的负载，而补热的功率大于等于工厂的负载

			P22 = cpm_h*(TsList[i] - TsList[j])/(COP2((TsList[i] + TsList[j]) / 2 - dT_EvaporationStandard, Tuse)-1)
			P12 = (heatLoad-P22-cpm_h*(TsList[i] - TsList[j]))/COP1
			if P21+elecP1 < P22+P12
				C[i, j] = P21 * dt + elecP1
				P2Matrix[i, j] = P21
				PeMatrix[i, j] = elecP1 / dt
			else
				C[i, j] = (P22+P12)*dt
				P2Matrix[i, j] = P22
				P1Matrix[i, j] = P12
			end
		end
		#=
		# 纯蓄热供热的工况下不考虑电加热的能量守恒会导致严重的计算错误,但在时间步长减小时由于温度降低幅度减小，又要依附到网格上，导致电加热的补齐效应一直开启。
		P2 = heatLoad / (COP2((TsList[i] + TsList[index]) / 2 - dT_EvaporationStandard, Tuse))
		C[i, index] = P2 * dt
		P2Matrix[i, index] = P2
		=#

		for j ∈ index+1:i-1
			P2 = heatLoad / (COP2((TsList[i] + TsList[j]) / 2 - dT_EvaporationStandard, Tuse))
			elecP = (TsList[j] - TsNextMin) * cpm_h
			#补热的功率基本上超不了，因为这个工况下从蓄热取走的热量小于工厂的负载，而补热的功率大于等于工厂的负载
			if elecP < PelecHeatMax * dt
				C[i, j] = P2 * dt + elecP
				P2Matrix[i, j] = P2
				PeMatrix[i, j] = elecP / dt
			end
		end

		# 计算温度不变的电度
		C[i, i] = P1Only * dt

		# 计算温度上升的功率，首先计算该模式的转换温度。升温模式有两种工作状态：1.热泵直供，电加热提高蓄热温度；2.热泵按照两者中高的供热
		# 需要一个计算温度低于Tuse和温度高于Tuse时功率的函数
		j = i + 1
		flag = true
		while flag && j < nT
			if TsList[j] <= Tuse - dT_EvaporationStandard
				C1, flag, P1, P3, Pe = powerCalculate_lowTs(TsList[j], TsList[i], dt)
				if flag
					C[i, j] = C1
					P1Matrix[i, j] = P1
					P3Matrix[i, j] = P3
					PeMatrix[i, j] = Pe
				end
				#C[i, j] = flag ? C1 : 99
			elseif TsList[i] < Tuse - dT_EvaporationStandard < TsList[j]
				dt1 = (Tuse - dT_EvaporationStandard - TsList[i]) / (TsList[j] - TsList[i]) * dt
				dt2 = dt - dt1
				C1, flag1, P11, P31, Pe1 = powerCalculate_lowTs(Tuse - dT_EvaporationStandard, TsList[i], dt1)
				C2, flag2, P12, P32, Pe2 = powerCalculate_highTs(TsList[j], Tuse - dT_EvaporationStandard, dt2)
				flag = flag1 && flag2
				if flag
					C[i, j] = C1 + C2
					P1Matrix[i, j] = (P11 * dt1 + P12 * dt2) / dt
					P3Matrix[i, j] = (P31 * dt1 + P32 * dt2) / dt
					PeMatrix[i, j] = (Pe1 * dt1 + Pe2 * dt2) / dt
				end
			elseif TsList[i] >= Tuse - dT_EvaporationStandard
				C1, flag, P1, P3, Pe = powerCalculate_highTs(TsList[j], TsList[i], dt)
				if flag
					C[i, j] = C1
					P1Matrix[i, j] = P1
					P3Matrix[i, j] = P3
					PeMatrix[i, j] = Pe
				end
			end
			j += 1
		end
		TsIncreaseIndexList[i] = flag ? j - 1 - i : j - 2 - i
	end

	return C, P1Matrix, P2Matrix, P3Matrix, PeMatrix
end

function getStateTransitionCost(::PressedWaterOneStorageOneCompressor, ::VaryLoadVaryArea;
	COPOverlap::Function,
	COPWater::Function,
	heatLoad::Vector,#热负荷
	Tair::Vector,# 外部环境温度
	# 总循环参数
	Tuse::Real,# 供热蒸汽温度
	dT_EvaporationStandard::Real,#全蒸温差
	latentHeat::Real,# 汽化潜热
	cp_cw::Real,# 循环水定压热容
	TcChangeToElec::Real,
	TWaste::Real,# 废热回收蒸发器温度

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容

	# 设备运行约束
	PheatPumpMax::Real,# 热泵最大功率
	PelecHeatMax::Real,# 电锅炉最大功率

	# 求解参数
	TsListStart::Matrix,# 状态参数高温蓄热温度起始值列表
	TsListEnd::Matrix,# 状态参数高温蓄热温度结束值列表
	dt::Real = 0.1,# 时间步长
)
	COP1 = COPOverlap(TWaste, Tuse)
	COP2 = COPWater
	COP3 = COPOverlap

	nT = length(TsListStart)# 温度步数

	#C[i,j]表示温度从TsList[i]到TsList[j]时的最低功率；在负载、环境不变的情况下，C[i,j]是不变的
	C = fill(99.0, nT, nT)# 状态转移矩阵
	P1Matrix = zeros(nT, nT)# 状态转移功率1参数
	P2Matrix = zeros(nT, nT)# 状态转移功率2参数
	P3Matrix = zeros(nT, nT)# 状态转移功率3参数
	PeMatrix = zeros(nT, nT)# 状态转移功率电加热参数
	TsDecreaseIndexList = zeros(Int, nT)# 记录温度下降最多偏移的index
	TsIncreaseIndexList = zeros(Int, nT)# 记录温度上升最多偏移的index

	# 只用热泵供热时的功率
	P1Only = heatLoad / COP1# 只用热泵供热时的功率

	"""计算蓄热温度相对较低时的电度情况"""
	function powerCalculate_lowTs(Tsaim, Tsstart, dt)
		P3 = min(PheatPumpMax - P1Only, cpm_h / dt / COP1 * (Tsaim - Tsstart))
		Pe = cpm_h / dt * (Tsaim - Tsstart) - P3 * COP1
		flag = Pe <= PelecHeatMax
		return (P1Only + P3 + Pe) * dt, flag, P1Only, P3, Pe
	end

	"""计算蓄热温度相对较高时的电度情况"""
	function powerCalculate_highTs(Tsaim, Tsstart, dt)
		COP3value = COP3(TWaste, (Tsaim + Tsstart) / 2 + dT_EvaporationStandard)

		function mode1(dt1)# 模式1：电热加热
			Pe = cpm_h / dt1 * (Tsaim - Tsstart)
			if Pe <= PelecHeatMax
				return true, P1Only, 0.0, Pe
			else
				return false, 9999.0, 9999.0, 9999.0
			end
			
		end

		function mode2(Tsaim2, Tsstart2, dt2)
			P3 = min(
				PheatPumpMax,
				(PheatPumpMax - heatLoad / COP3value) / (1 - cp_cw / latentHeat * (Tsaim2 - Tsstart2)),
				cpm_h / dt2 / COP3value * (Tsaim2 - Tsstart2),
			)
			P1 = max(0, heatLoad / COP3value - cp_cw / latentHeat * (Tsaim2 - Tsstart2) * P3)
			Pe = cpm_h / dt2 * (Tsaim2 - Tsstart2) - P3 * COP3value
			if Pe <= PelecHeatMax
				return true, P1, P3, Pe
			else
				return false, 9999.0, 9999.0, 9999.0
			end
		end
		#COP3value = COP3(TWaste, (min(Tsaim,TcChangeToElec - dT_EvaporationStandard) + Tsstart) / 2 + dT_EvaporationStandard)
		
		Ptotal = 999.0
		P1res = 0.0
		P3res = 0.0
		Peres = 0.0
		flag = false
		# 模式1：电热加热
		flag, P1res, P3res, Peres= mode1(dt)
		Ptotal=P1res+P3res+Peres

		# 模式2：热泵加电热,又分两种模式：末态温度不大于电加热温度界限与末态温度大于电加热温度界线
		if (Tsaim <= TcChangeToElec - dT_EvaporationStandard)	# 该时间层的末态温度大于电加热温度界限
			flag_2, P1_2, P3_2, Pe_2= mode2(Tsaim,Tsstart,dt)
			Ptotal_2=P1_2+P3_2+Pe_2
		else
			tmid=dt*(TcChangeToElec - dT_EvaporationStandard-Tsstart)/(Tsaim - Tsstart)
			flag_21, P1_21, P3_21, Pe_21= mode2(TcChangeToElec - dT_EvaporationStandard,Tsstart,tmid)
			flag_22, P1_22, P3_22, Pe_22= mode2(Tsaim,TcChangeToElec - dT_EvaporationStandard,dt-tmid)
			flag_2=flag_21 && flag_22
			P1_2=P1_21*tmid/dt+P1_22*(dt-tmid)/dt
			P3_2=P3_21*tmid/dt+P3_22*(dt-tmid)/dt
			Pe_2=Pe_21*tmid/dt+Pe_22*(dt-tmid)/dt
			Ptotal_2=P1_2+P3_2+Pe_2
		end
		if (Ptotal_2 < Ptotal)
			flag = flag || flag_2
			P1res = P1_2
			P3res = P3_2
			Peres = Pe_2
			Ptotal = Ptotal_2
		end
		return Ptotal * dt, flag, P1res, P3res, Peres
	end

	for i ∈ 1:nT
		# 先计算温度下降的功率,存在一边电加热一边开2号模式的情况
		# 首先计算该模式在下一刻的最低温度
		TsNextMin = TsList[i]
		# 简单迭代求蓄热的最低温度
		for _ ∈ 1:3
			TsNextMin = TsList[i] - heatLoad * dt / cpm_h * (1 - 1 / COP2(0.5 * (TsNextMin + TsList[i]) - dT_EvaporationStandard, Tuse))
		end
		# 计算最低温度的下标
		ibias=findfirst(x -> x >= TsNextMin, TsList)
		TsDecreaseIndexList[i]=i-ibias
		TsNextMinGrid = TsList[ibias]
		elecP = (TsNextMinGrid - TsNextMin) * cpm_h

		# 如果温度在一个时间层内会跌到低于最低温度，那么将这个时间层分为两段，前段用蓄热，后段用热泵直供
		index = i - TsDecreaseIndexList[i]
		if index == 1
			j=1
			P21 = heatLoad / (COP2((TsList[i] + TsList[j]) / 2 - dT_EvaporationStandard, Tuse))
			elecP1 = (TsList[j] - TsNextMin) * cpm_h
			#补热的功率基本上超不了，因为这个工况下从蓄热取走的热量小于工厂的负载，而补热的功率大于等于工厂的负载

			P22 = cpm_h*(TsList[i] - TsList[j])/(COP2((TsList[i] + TsList[j]) / 2 - dT_EvaporationStandard, Tuse)-1)
			P12 = (heatLoad-P22-cpm_h*(TsList[i] - TsList[j]))/COP1
			if P21+elecP1 < P22+P12
				C[i, j] = P21 * dt + elecP1
				P2Matrix[i, j] = P21
				PeMatrix[i, j] = elecP1 / dt
			else
				C[i, j] = (P22+P12)*dt
				P2Matrix[i, j] = P22
				P1Matrix[i, j] = P12
			end
		end
		#=
		# 纯蓄热供热的工况下不考虑电加热的能量守恒会导致严重的计算错误,但在时间步长减小时由于温度降低幅度减小，又要依附到网格上，导致电加热的补齐效应一直开启。
		P2 = heatLoad / (COP2((TsList[i] + TsList[index]) / 2 - dT_EvaporationStandard, Tuse))
		C[i, index] = P2 * dt
		P2Matrix[i, index] = P2
		=#

		for j ∈ index+1:i-1
			P2 = heatLoad / (COP2((TsList[i] + TsList[j]) / 2 - dT_EvaporationStandard, Tuse))
			elecP = (TsList[j] - TsNextMin) * cpm_h
			#补热的功率基本上超不了，因为这个工况下从蓄热取走的热量小于工厂的负载，而补热的功率大于等于工厂的负载
			if elecP < PelecHeatMax * dt
				C[i, j] = P2 * dt + elecP
				P2Matrix[i, j] = P2
				PeMatrix[i, j] = elecP / dt
			end
		end

		# 计算温度不变的电度
		C[i, i] = P1Only * dt

		# 计算温度上升的功率，首先计算该模式的转换温度。升温模式有两种工作状态：1.热泵直供，电加热提高蓄热温度；2.热泵按照两者中高的供热
		# 需要一个计算温度低于Tuse和温度高于Tuse时功率的函数
		j = i + 1
		flag = true
		while flag && j < nT
			if TsList[j] <= Tuse - dT_EvaporationStandard
				C1, flag, P1, P3, Pe = powerCalculate_lowTs(TsList[j], TsList[i], dt)
				if flag
					C[i, j] = C1
					P1Matrix[i, j] = P1
					P3Matrix[i, j] = P3
					PeMatrix[i, j] = Pe
				end
				#C[i, j] = flag ? C1 : 99
			elseif TsList[i] < Tuse - dT_EvaporationStandard < TsList[j]
				dt1 = (Tuse - dT_EvaporationStandard - TsList[i]) / (TsList[j] - TsList[i]) * dt
				dt2 = dt - dt1
				C1, flag1, P11, P31, Pe1 = powerCalculate_lowTs(Tuse - dT_EvaporationStandard, TsList[i], dt1)
				C2, flag2, P12, P32, Pe2 = powerCalculate_highTs(TsList[j], Tuse - dT_EvaporationStandard, dt2)
				flag = flag1 && flag2
				if flag
					C[i, j] = C1 + C2
					P1Matrix[i, j] = (P11 * dt1 + P12 * dt2) / dt
					P3Matrix[i, j] = (P31 * dt1 + P32 * dt2) / dt
					PeMatrix[i, j] = (Pe1 * dt1 + Pe2 * dt2) / dt
				end
			elseif TsList[i] >= Tuse - dT_EvaporationStandard
				C1, flag, P1, P3, Pe = powerCalculate_highTs(TsList[j], TsList[i], dt)
				if flag
					C[i, j] = C1
					P1Matrix[i, j] = P1
					P3Matrix[i, j] = P3
					PeMatrix[i, j] = Pe
				end
			end
			j += 1
		end
		TsIncreaseIndexList[i] = flag ? j - 1 - i : j - 2 - i
	end

	return C, P1Matrix, P2Matrix, P3Matrix, PeMatrix, TsDecreaseIndexList, TsIncreaseIndexList
end

function generateAndSolve(::PressedWaterOneStorageOneCompressor, ::MinimizeCost, ::VaryLoadVaryArea, ::GoldenRatioMethod;
	COPOverlap::Function,
	COPWater::Function,
	hourlyTariffFunction::Function,   # 电价函数
	heatConsumptionPowerFunction::Function,  # 用热负载函数
	TairFunction::Function,# 环境温度函数

	# 总循环参数
	Tuse::Real,# 供热蒸汽温度
	TCompressorIn::Real,# 中间级温度
	dT_EvaporationStandard::Real,#全蒸温差
	latentHeat::Real,# 汽化潜热
	cp_cw::Real,# 循环水定压热容
	cp_cs::Real,# 蒸汽定压热容 cp cycled steam
	TcChangeToElec::Real,
	TWaste::Real,# 废热回收蒸发器温度

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容

	# 设备运行约束
	TstorageTankMax::Real,# 蓄热罐的最高温度
	PheatPumpMax::Real,# 热泵最大功率
	PelecHeatMax::Real,# 电锅炉最大功率
	# 求解参数
	dT::Real = 0.01,# 状态参数高温蓄热温度离散步长
	dt::Real = 1 / 6,# 时间步长
)
	tList = 0:dt:24

	TsList = repeat(TCompressorIn+dT_EvaporationStandard:1.0:TstorageTankMax,1,length(tList))
	
	nT = size(TsList,1)# 温度步数
	nt = length(tList)# 时间步数
	TsMatrix = zeros(Int, nt, nT)# 存储状态参数：高温蓄热温度
	bestValueList = fill(99.0, nT)

	# 已经生成了C, TsDecreaseIndexList, TsIncreaseIndexList


	C, P1Matrix, P2Matrix, P3Matrix, PeMatrix, TsDecreaseIndexList, TsIncreaseIndexList = getStateTransitionCost(
		PressedWaterOneStorageOneCompressor(),
		ConstloadandArea();
		COPOverlap = COPOverlap,
		COPWater = COPWater,
		heatLoad = heatLoad,#热负荷
		Tair = Tair,# 外部环境温度
		# 总循环参数
		Tuse = Tuse,# 供热蒸汽温度
		dT_EvaporationStandard = dT_EvaporationStandard,#全蒸温差
		latentHeat = latentHeat,# 汽化潜热
		cp_cw = cp_cw,# 循环水定压热容
		TcChangeToElec = TcChangeToElec,
		TWaste = TWaste,# 废热回收蒸发器温度

		# 高温蓄热参数
		cpm_h = cpm_h,# 高温蓄热热容

		# 设备运行约束
		PheatPumpMax = PheatPumpMax,# 热泵最大功率
		PelecHeatMax = PelecHeatMax,# 电锅炉最大功率
		# 求解参数
		TsList = collect(TsList),# 状态参数高温蓄热温度列表
		dt = dt,# 时间步长
	)
	TsDecreaseIndex = maximum(TsDecreaseIndexList)
	TsIncreaseIndex = maximum(TsIncreaseIndexList)
	#最多降温TsDecreaseIndex,最多升温TsIncreaseIndex,所以从目标出发要多TsDecreaseIndex个，少TsIncreaseIndex个
	# 先直接正向计算
	# 正向计算步数
	count = 0
	# 初始化黄金分割法
	phi = 0.618
	jList = [1, nT - round(Int, phi * (nT - 1)), 1 + round(Int, phi * (nT - 1)), nT]
	valueList = zeros(4)
	TsIndexListMatirx = zeros(Int, nt, 4)
	for (i, j) in enumerate(jList)
		valueList[i], TsIndexListMatirx[:, i] = dpSolve(
			nT,
			nt,
			C,
			hourlyTariffFunction,
			dt,
			j,
			TsDecreaseIndex,# 温度下降最多偏移的index
			TsIncreaseIndex,# 温度上升最多偏移的index
			ConstloadandArea(),
		)
	end
	count = 0
	while (jList[4] - jList[1] >= 4) && count < 100
		if valueList[2] < valueList[3]# 最值在j1到j3之间
			jList[4] = jList[3]
			valueList[4] = valueList[3]
			TsIndexListMatirx[:, 4] = TsIndexListMatirx[:, 3]
			jList[3] = jList[2]
			valueList[3] = valueList[2]
			TsIndexListMatirx[:, 3] = TsIndexListMatirx[:, 2]
			jList[2] = min(max(floor(Int, jList[4] - phi * (jList[4] - jList[1])), jList[1] + 1), jList[3] - 1)
			valueList[2], TsIndexListMatirx[:, 2] = dpSolve(
				nT,
				nt,
				C,
				hourlyTariffFunction,
				dt,
				jList[2],
				TsDecreaseIndex,# 温度下降最多偏移的index
				TsIncreaseIndex,# 温度上升最多偏移的index
				ConstloadandArea(),
			)
		elseif valueList[2] >= valueList[3]
			jList[1] = jList[2]
			valueList[1] = valueList[2]
			TsIndexListMatirx[:, 1] = TsIndexListMatirx[:, 2]
			jList[2] = jList[3]
			valueList[2] = valueList[3]
			TsIndexListMatirx[:, 2] = TsIndexListMatirx[:, 3]
			jList[3] = max(min(ceil(Int, jList[1] + phi * (jList[4] - jList[1])), jList[4] - 1), jList[2] + 1)
			valueList[3], TsIndexListMatirx[:, 3] = dpSolve(
				nT,
				nt,
				C,
				hourlyTariffFunction,
				dt,
				jList[3],
				TsDecreaseIndex,# 温度下降最多偏移的index
				TsIncreaseIndex,# 温度上升最多偏移的index
				ConstloadandArea(),
			)
		end
		count += 1
		println("$count/$nT", "j1=$(jList[1]),j2=$(jList[2]),j3=$(jList[3]),j4=$(jList[4])")
	end

	minCost, index = findmin(valueList)
	minTsList = TsIndexListMatirx[:, index]

	P1List = [P1Matrix[minTsList[i], minTsList[i+1]] for i in 1:nt-1]
	P2List = [P2Matrix[minTsList[i], minTsList[i+1]] for i in 1:nt-1]
	P3List = [P3Matrix[minTsList[i], minTsList[i+1]] for i in 1:nt-1]
	PeList = [PeMatrix[minTsList[i], minTsList[i+1]] for i in 1:nt-1]

	return minCost, TsList[minTsList], P1List, P2List, P3List, PeList
end


