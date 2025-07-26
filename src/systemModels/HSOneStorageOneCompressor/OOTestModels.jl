"""000"""
function getStateTransitionCost_SingleStep(
	::OO000;
	#C::Matrix, P1Matrix::Matrix, P2Matrix::Matrix, P3Matrix::Matrix, PeMatrix::Matrix;
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
	TCompressorIn::Real,# 压缩机吸气温度

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容

	# 设备运行约束
	PheatPumpMax::Real,# 热泵最大功率
    PelecHeatMax::Real,# 电锅炉最大功率
    PWaterCompressorMax::Real,# 水蒸气压缩机的最大功率
	#PheatPumpLowMax::Real,# 低温热泵的最大功率,一般来说，水蒸气压缩机和低温热泵同时运行时，由于低温热泵按照最劣工况设计，所以水蒸气压缩机会先触碰运行边界，不需要考虑低温热泵的边界问题
    Tsmin::Real,# 蓄热的最小温度
    
	# 求解参数
	TsListStart::Vector,# 状态参数高温蓄热温度起始值列表
	TsListEnd::Vector,# 状态参数高温蓄热温度结束值列表
	dt::Real = 1.0,# 时间步长
	Tsmax::Real = 220.0
)
	COP1 = COPOverlap(TWaste, Tuse)
	COP2 = COPWater
	COP2_design = COP2(TCompressorIn,Tuse)
	COP3 = COPOverlap
	heatLoadPumpMax=PWaterCompressorMax*COP2_design	# 热泵设计的供热功率

	nT = length(TsListStart)# 温度步数

	#C[i,j]表示温度从TsList[i]到TsList[j]时的最低功率；在负载、环境不变的情况下，C[i,j]是不变的
	C = fill(9999.0, nT, nT)# 状态转移矩阵
	P1Matrix = zeros(nT, nT)# 状态转移功率1参数
	P2Matrix = zeros(nT, nT)# 状态转移功率2参数
	P3Matrix = zeros(nT, nT)# 状态转移功率3参数
	PeMatrix = zeros(nT, nT)# 状态转移功率电加热参数

	
	# 只用热泵供热时的功率
	P1Only = min(heatLoad, heatLoadPumpMax)/ COP1# 只用热泵供热时的功率
	P1hOnly = P1Only*COP1/COP2_design	#只用热泵供热时的压缩机功率

	"""计算不同蓄热温度下蓄热温度降低的电度"""
	function powerCalculate_disCharge(Tsaim, Tsstart, dt)
		COP2value=COP2((Tsaim + Tsstart) / 2 - dT_EvaporationStandard, Tuse)
		Psout=cpm_h*(Tsstart-Tsaim)/dt
		P2=Psout / (COP2value-1)
		if (P2+Psout>heatLoad) || (P2 > PWaterCompressorMax) # 降不到这个温度
			return 9999.0,false,9999.0,9999.0,9999.0
		end
        #先看是不是温度到下限了
        if Tsaim == Tsmin
            P1 = min((heatLoad-P2-Psout) / COP1,PheatPumpMax)
            Pe = heatLoad-P2-Psout-P1 * COP1
            P21=min(PWaterCompressorMax,heatLoad/COP2value)  # 前半段蓄热供热最大P2功率
            Pe1=heatLoad-P21*COP2value
            P12 = min(heatLoad / COP1,P1Only)   # 后半段热泵供热的功率,P1Only保证了压缩机功率不超过设计功率，所以后面不用比较压缩机功率了
			Pe2 = heatLoad-P12 * COP1
			#P1h=P1*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P1l=P1-P1h
			#P12h=P12*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P12l=P12-P12h
            if (Pe1 <= PelecHeatMax) && (Pe2 <= PelecHeatMax) #&& (P1h<=PWaterCompressorMax) && (P12h<=PWaterCompressorMax) #&& (P1l <= PheatPumpLowMax) && (P12l <= PheatPumpLowMax)
				#println("放热2 Pe1=$Pe1 PelecHeatMax=$PelecHeatMax,Pe2=$Pe2,PelecHeatMax=$PelecHeatMax")
			    return (P1+P2+Pe)*dt,true,P1,P2,Pe
            else
				#println("放热3")
                return 9999.0,false,9999.0,9999.0,9999.0
            end
        end

        #蓄热温度的末态高于用热温度,放热全程热泵补热;或放热到蓄热的最低温度,先放热后热泵
		if (Tsaim-dT_EvaporationStandard>=Tuse)
			#P1 = min((heatLoad-P2-Psout) / COP1,P1Only)
            #=OOTest001=#
            P1=0
			Pe = heatLoad-P2-Psout-P1 * COP1
			#P1h=P1*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P1l=P1-P1h
            if (Pe <= PelecHeatMax) #&& (P1h<=PWaterCompressorMax) #&& (P1l <= PheatPumpLowMax)
			    return (P1+P2+Pe)*dt,true,P1,P2,Pe
            else
                return 9999.0,false,9999.0,9999.0,9999.0
            end
		elseif Tsstart-dT_EvaporationStandard<=Tuse#蓄热温度的初态低于用热温度，放热全程电热补热
			P1 = 0
			Pe = heatLoad-P2-Psout
            if Pe <= PelecHeatMax
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
		# 主要看水蒸气压缩机有没有余量
		P3h = min(PWaterCompressorMax-P1hOnly,cpm_h * (Tsaim - Tsstart) / dt / COP2_design / (1+cp_cw/latentHeat*(Tuse-0.5*(Tsstart+Tsaim))))

		P3 = P3h*COP2_design/COP1
		Pe = heatLoad+cpm_h / dt * (Tsaim - Tsstart) - (P3 + P1Only)*COP1
		
		#P3l = P3 - P3h
		if (Pe <= PelecHeatMax) #&& (P3h<=PWaterCompressorMax) #&& (P3l <= PheatPumpLowMax)
			#println("蓄热1 Pe=$Pe,PelecHeatMax=$PelecHeatMax")
			return (P1Only + P3 + Pe) * dt, true, P1Only, P3, Pe
		else
			#println("蓄热2 Pe=$Pe,PelecHeatMax=$PelecHeatMax")
			return 9999.0,false,9999.0,9999.0,9999.0
		end
	end

	"""计算蓄热温度相对较高时升温的电度情况(考虑蓄热温度超过热泵温度上限)"""
	function powerCalculate_highTs(Tsaim, Tsstart, dt)
		COP3value = COP3(TWaste, (Tsaim + Tsstart) / 2 + dT_EvaporationStandard)
		COPWatervalue = COP2(TCompressorIn, (Tsaim + Tsstart) / 2 + dT_EvaporationStandard)
		function mode1(dt1)# 模式1：只有电热加热给蓄热储热
			Pe = cpm_h / dt1 * (Tsaim - Tsstart)+heatLoad-P1Only*COP1
			if Pe <= PelecHeatMax
				return true, P1Only, 0.0, Pe
			else
				return false, 9999.0, 9999.0, 9999.0
			end
			
		end

		# 模式2：热泵也给蓄热储热（此时热泵COP按照蓄热温度计算）
		function mode2(Tsaim2, Tsstart2, dt2)
			# 如果水蒸气压缩机容量没有富余的话，应该优先用热泵为工厂供热，而不是为蓄热罐储热
			Psin=cpm_h / dt2 * (Tsaim2 - Tsstart2)#蓄热储热功率
			P3hin=min(Psin/COPWatervalue,PWaterCompressorMax)
            #=OOTest100=#
			#k=cp_cw/latentHeat * (Tsaim2 - Tsstart2)#闪蒸回收热比例
            k=0
			u_Pm=(heatLoad/COPWatervalue+PWaterCompressorMax)/(1+k)#线性规划求解的一个分类指标,用来描述用热负荷的大小
			# 先看闪蒸回收热是否够系统需求,如果够的话多余的热直接放掉
			if P3hin*k*COPWatervalue>=heatLoad
				return true, 0.0, P3hin, Psin-P3hin*COPWatervalue
			end
			
			if P3hin>=PWaterCompressorMax && u_Pm>=PWaterCompressorMax
				#蓄热负荷和用热负荷都大,任意进行热泵配置
				P1h=PWaterCompressorMax
				P3h=0.0
			elseif P3hin<PWaterCompressorMax && u_Pm>PWaterCompressorMax
				#蓄热负荷小，用热负荷大,蓄热占满，剩余的用热
				P3h=P3hin
				P1h=min(PWaterCompressorMax,(P3hin+heatLoad/COPWatervalue)/(1+k))-P3h
			elseif P3hin>PWaterCompressorMax && u_Pm<PWaterCompressorMax
				#蓄热负荷大，用热负荷小,蓄热按照供热程度占满
				P1h=0.0
				P3h=heatLoad/COPWatervalue/k
			else
				#蓄热负荷和用热负荷都小
				P3h=min(heatLoad/COPWatervalue/k,P3hin)
				P1h=min(heatLoad/COPWatervalue/k,(P3hin+heatLoad/COPWatervalue)/(1+k))-P3h
			end

			Pe = cpm_h / dt2 * (Tsaim2 - Tsstart2) + heatLoad - (P1h+P3h) * COPWatervalue
			if Pe <= PelecHeatMax
				# flag P1 P3 Pe
				return true, P1h*COPWatervalue/COP1, P3h*COPWatervalue/COP1, Pe
			else
				return false, 9999.0, 9999.0, 9999.0
			end
		end
		#COP3value = COP3(TWaste, (min(Tsaim,TcChangeToElec - dT_EvaporationStandard) + Tsstart) / 2 + dT_EvaporationStandard)
		
		Ptotal = 9999.0
		P1res = 0.0
		P3res = 0.0
		Peres = 0.0
		flag = false
		# 模式1：电热加热
		flag, P1res, P3res, Peres= mode1(dt)
		Ptotal=P1res+P3res+Peres

		# 模式2：热泵加电热,又分两种模式：末态温度不大于电加热温度界限与末态温度大于电加热温度界线
		if (Tsaim <= TcChangeToElec - dT_EvaporationStandard)	# 该时间层的末态温度小于电加热温度界限
			flag_2, P1_2, P3_2, Pe_2= mode2(Tsaim,Tsstart,dt)
			Ptotal_2=P1_2+P3_2+Pe_2
		elseif (Tsaim > TcChangeToElec - dT_EvaporationStandard>Tsstart) # 该时间层的温度跨越电加热温度界限
			tmid=dt*(TcChangeToElec - dT_EvaporationStandard-Tsstart)/(Tsaim - Tsstart)
			flag_21, P1_21, P3_21, Pe_21= mode2(TcChangeToElec - dT_EvaporationStandard,Tsstart,tmid)
			flag_22, P1_22, P3_22, Pe_22= mode2(Tsaim,TcChangeToElec - dT_EvaporationStandard,dt-tmid)
			flag_2=flag_21 && flag_22
			P1_2=P1_21*tmid/dt+P1_22*(dt-tmid)/dt
			P3_2=P3_21*tmid/dt+P3_22*(dt-tmid)/dt
			Pe_2=Pe_21*tmid/dt+Pe_22*(dt-tmid)/dt
			Ptotal_2=P1_2+P3_2+Pe_2
		else# 该时间层初态温度高于电加热温度界线，切换模式1
			flag_2, P1_2, P3_2, Pe_2=mode1(dt)
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

	for (j,Tsaim) in enumerate(TsListEnd)
		for (i,Tsstart) in enumerate(TsListStart)
			# 先计算温度下降的功率,存在一边电加热一边开2号模式的情况
			if Tsstart<Tsmin || Tsaim>Tsmax
				C[i, j] = 9999.0
				P1Matrix[i, j] = 9999.0
				P3Matrix[i, j] = 9999.0
				PeMatrix[i, j] = 9999.0
			elseif Tsaim < Tsstart
				costValue,flag,P1,P2,Pe=powerCalculate_disCharge(Tsaim, Tsstart, dt)
				if flag
					C[i, j]=costValue
					P1Matrix[i, j] = P1
					P2Matrix[i, j] = P2
					PeMatrix[i, j] = Pe
				end
			elseif Tsaim == Tsstart
				temp = heatLoad-P1Only*COP1
				if temp <= PelecHeatMax
					P1Matrix[i, j] = P1Only
					PeMatrix[i, j] = temp
					C[i, j]=(P1Matrix[i, j]+PeMatrix[i, j])*dt
				end
				#=
				P1Matrix[i, j] = P1Only
				PeMatrix[i, j] = heatLoad-P1Only*COP1
				C[i, j]=(P1Matrix[i, j]+PeMatrix[i, j])*dt
				=#
			elseif Tsaim > Tsstart
				if Tsaim <= Tuse - dT_EvaporationStandard
					C1, flag, P1, P3, Pe = powerCalculate_lowTs(Tsaim, Tsstart, dt)
					if flag
						C[i, j] = C1
						P1Matrix[i, j] = P1
						P3Matrix[i, j] = P3
						PeMatrix[i, j] = Pe
					end
					#C[i, j] = flag ? C1 : 99
				elseif Tsstart < Tuse - dT_EvaporationStandard < Tsaim
					dt1 = (Tuse - dT_EvaporationStandard - Tsstart) / (Tsaim - Tsstart) * dt
					dt2 = dt - dt1
					C1, flag1, P11, P31, Pe1 = powerCalculate_lowTs(Tuse - dT_EvaporationStandard, Tsstart, dt1)
					C2, flag2, P12, P32, Pe2 = powerCalculate_highTs(Tsaim, Tuse - dT_EvaporationStandard, dt2)
					flag = flag1 && flag2
					if flag
						C[i, j] = C1 + C2
						P1Matrix[i, j] = (P11 * dt1 + P12 * dt2) / dt
						P3Matrix[i, j] = (P31 * dt1 + P32 * dt2) / dt
						PeMatrix[i, j] = (Pe1 * dt1 + Pe2 * dt2) / dt
					end
				elseif Tsstart >= Tuse - dT_EvaporationStandard
					C1, flag, P1, P3, Pe = powerCalculate_highTs(Tsaim, Tsstart, dt)
					if flag
						C[i, j] = C1
						P1Matrix[i, j] = P1
						P3Matrix[i, j] = P3
						PeMatrix[i, j] = Pe
					end
				end
			end
		end
	end

	return C, P1Matrix, P2Matrix, P3Matrix, PeMatrix
end

"""001"""
function getStateTransitionCost_SingleStep(
	::OO001;
	#C::Matrix, P1Matrix::Matrix, P2Matrix::Matrix, P3Matrix::Matrix, PeMatrix::Matrix;
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
	TCompressorIn::Real,# 压缩机吸气温度

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容

	# 设备运行约束
	PheatPumpMax::Real,# 热泵最大功率
    PelecHeatMax::Real,# 电锅炉最大功率
    PWaterCompressorMax::Real,# 水蒸气压缩机的最大功率
	#PheatPumpLowMax::Real,# 低温热泵的最大功率,一般来说，水蒸气压缩机和低温热泵同时运行时，由于低温热泵按照最劣工况设计，所以水蒸气压缩机会先触碰运行边界，不需要考虑低温热泵的边界问题
    Tsmin::Real,# 蓄热的最小温度
    
	# 求解参数
	TsListStart::Vector,# 状态参数高温蓄热温度起始值列表
	TsListEnd::Vector,# 状态参数高温蓄热温度结束值列表
	dt::Real = 1.0,# 时间步长
	Tsmax::Real = 220.0
)
	COP1 = COPOverlap(TWaste, Tuse)
	COP2 = COPWater
	COP2_design = COP2(TCompressorIn,Tuse)
	COP3 = COPOverlap
	heatLoadPumpMax=PWaterCompressorMax*COP2_design	# 热泵设计的供热功率

	nT = length(TsListStart)# 温度步数

	#C[i,j]表示温度从TsList[i]到TsList[j]时的最低功率；在负载、环境不变的情况下，C[i,j]是不变的
	C = fill(9999.0, nT, nT)# 状态转移矩阵
	P1Matrix = zeros(nT, nT)# 状态转移功率1参数
	P2Matrix = zeros(nT, nT)# 状态转移功率2参数
	P3Matrix = zeros(nT, nT)# 状态转移功率3参数
	PeMatrix = zeros(nT, nT)# 状态转移功率电加热参数

	
	# 只用热泵供热时的功率
	P1Only = min(heatLoad, heatLoadPumpMax)/ COP1# 只用热泵供热时的功率
	P1hOnly = P1Only*COP1/COP2_design	#只用热泵供热时的压缩机功率

	"""计算不同蓄热温度下蓄热温度降低的电度"""
	function powerCalculate_disCharge(Tsaim, Tsstart, dt)
		COP2value=COP2((Tsaim + Tsstart) / 2 - dT_EvaporationStandard, Tuse)
		Psout=cpm_h*(Tsstart-Tsaim)/dt
		P2=Psout / (COP2value-1)
		if (P2+Psout>heatLoad) || (P2 > PWaterCompressorMax) # 降不到这个温度
			return 9999.0,false,9999.0,9999.0,9999.0
		end
        #先看是不是温度到下限了
        if Tsaim == Tsmin
            P1 = min((heatLoad-P2-Psout) / COP1,PheatPumpMax)
            Pe = heatLoad-P2-Psout-P1 * COP1
            P21=min(PWaterCompressorMax,heatLoad/COP2value)  # 前半段蓄热供热最大P2功率
            Pe1=heatLoad-P21*COP2value
            P12 = min(heatLoad / COP1,P1Only)   # 后半段热泵供热的功率,P1Only保证了压缩机功率不超过设计功率，所以后面不用比较压缩机功率了
			Pe2 = heatLoad-P12 * COP1
			#P1h=P1*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P1l=P1-P1h
			#P12h=P12*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P12l=P12-P12h
            if (Pe1 <= PelecHeatMax) && (Pe2 <= PelecHeatMax) #&& (P1h<=PWaterCompressorMax) && (P12h<=PWaterCompressorMax) #&& (P1l <= PheatPumpLowMax) && (P12l <= PheatPumpLowMax)
				#println("放热2 Pe1=$Pe1 PelecHeatMax=$PelecHeatMax,Pe2=$Pe2,PelecHeatMax=$PelecHeatMax")
			    return (P1+P2+Pe)*dt,true,P1,P2,Pe
            else
				#println("放热3")
                return 9999.0,false,9999.0,9999.0,9999.0
            end
        end

        #蓄热温度的末态高于用热温度,放热全程热泵补热;或放热到蓄热的最低温度,先放热后热泵
		if (Tsaim-dT_EvaporationStandard>=Tuse)
			P1 = min((heatLoad-P2-Psout) / COP1,P1Only)
            #=OOTest001=#
            #P1=0
			Pe = heatLoad-P2-Psout-P1 * COP1
			#P1h=P1*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P1l=P1-P1h
            if (Pe <= PelecHeatMax) #&& (P1h<=PWaterCompressorMax) #&& (P1l <= PheatPumpLowMax)
			    return (P1+P2+Pe)*dt,true,P1,P2,Pe
            else
                return 9999.0,false,9999.0,9999.0,9999.0
            end
		elseif Tsstart-dT_EvaporationStandard<=Tuse#蓄热温度的初态低于用热温度，放热全程电热补热
			P1 = 0
			Pe = heatLoad-P2-Psout
            if Pe <= PelecHeatMax
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
		# 主要看水蒸气压缩机有没有余量
		P3h = min(PWaterCompressorMax-P1hOnly,cpm_h * (Tsaim - Tsstart) / dt / COP2_design / (1+cp_cw/latentHeat*(Tuse-0.5*(Tsstart+Tsaim))))

		P3 = P3h*COP2_design/COP1
		Pe = heatLoad+cpm_h / dt * (Tsaim - Tsstart) - (P3 + P1Only)*COP1
		
		#P3l = P3 - P3h
		if (Pe <= PelecHeatMax) #&& (P3h<=PWaterCompressorMax) #&& (P3l <= PheatPumpLowMax)
			#println("蓄热1 Pe=$Pe,PelecHeatMax=$PelecHeatMax")
			return (P1Only + P3 + Pe) * dt, true, P1Only, P3, Pe
		else
			#println("蓄热2 Pe=$Pe,PelecHeatMax=$PelecHeatMax")
			return 9999.0,false,9999.0,9999.0,9999.0
		end
	end

	"""计算蓄热温度相对较高时升温的电度情况(考虑蓄热温度超过热泵温度上限)"""
	function powerCalculate_highTs(Tsaim, Tsstart, dt)
		COP3value = COP3(TWaste, (Tsaim + Tsstart) / 2 + dT_EvaporationStandard)
		COPWatervalue = COP2(TCompressorIn, (Tsaim + Tsstart) / 2 + dT_EvaporationStandard)
		function mode1(dt1)# 模式1：只有电热加热给蓄热储热
			Pe = cpm_h / dt1 * (Tsaim - Tsstart)+heatLoad-P1Only*COP1
			if Pe <= PelecHeatMax
				return true, P1Only, 0.0, Pe
			else
				return false, 9999.0, 9999.0, 9999.0
			end
			
		end

		# 模式2：热泵也给蓄热储热（此时热泵COP按照蓄热温度计算）
		function mode2(Tsaim2, Tsstart2, dt2)
			# 如果水蒸气压缩机容量没有富余的话，应该优先用热泵为工厂供热，而不是为蓄热罐储热
			Psin=cpm_h / dt2 * (Tsaim2 - Tsstart2)#蓄热储热功率
			P3hin=min(Psin/COPWatervalue,PWaterCompressorMax)
            #=OOTest100=#
			#k=cp_cw/latentHeat * (Tsaim2 - Tsstart2)#闪蒸回收热比例
            k=0
			u_Pm=(heatLoad/COPWatervalue+PWaterCompressorMax)/(1+k)#线性规划求解的一个分类指标,用来描述用热负荷的大小
			# 先看闪蒸回收热是否够系统需求,如果够的话多余的热直接放掉
			if P3hin*k*COPWatervalue>=heatLoad
				return true, 0.0, P3hin, Psin-P3hin*COPWatervalue
			end
			
			if P3hin>=PWaterCompressorMax && u_Pm>=PWaterCompressorMax
				#蓄热负荷和用热负荷都大,任意进行热泵配置
				P1h=PWaterCompressorMax
				P3h=0.0
			elseif P3hin<PWaterCompressorMax && u_Pm>PWaterCompressorMax
				#蓄热负荷小，用热负荷大,蓄热占满，剩余的用热
				P3h=P3hin
				P1h=min(PWaterCompressorMax,(P3hin+heatLoad/COPWatervalue)/(1+k))-P3h
			elseif P3hin>PWaterCompressorMax && u_Pm<PWaterCompressorMax
				#蓄热负荷大，用热负荷小,蓄热按照供热程度占满
				P1h=0.0
				P3h=heatLoad/COPWatervalue/k
			else
				#蓄热负荷和用热负荷都小
				P3h=min(heatLoad/COPWatervalue/k,P3hin)
				P1h=min(heatLoad/COPWatervalue/k,(P3hin+heatLoad/COPWatervalue)/(1+k))-P3h
			end

			Pe = cpm_h / dt2 * (Tsaim2 - Tsstart2) + heatLoad - (P1h+P3h) * COPWatervalue
			if Pe <= PelecHeatMax
				# flag P1 P3 Pe
				return true, P1h*COPWatervalue/COP1, P3h*COPWatervalue/COP1, Pe
			else
				return false, 9999.0, 9999.0, 9999.0
			end
		end
		#COP3value = COP3(TWaste, (min(Tsaim,TcChangeToElec - dT_EvaporationStandard) + Tsstart) / 2 + dT_EvaporationStandard)
		
		Ptotal = 9999.0
		P1res = 0.0
		P3res = 0.0
		Peres = 0.0
		flag = false
		# 模式1：电热加热
		flag, P1res, P3res, Peres= mode1(dt)
		Ptotal=P1res+P3res+Peres

		# 模式2：热泵加电热,又分两种模式：末态温度不大于电加热温度界限与末态温度大于电加热温度界线
		if (Tsaim <= TcChangeToElec - dT_EvaporationStandard)	# 该时间层的末态温度小于电加热温度界限
			flag_2, P1_2, P3_2, Pe_2= mode2(Tsaim,Tsstart,dt)
			Ptotal_2=P1_2+P3_2+Pe_2
		elseif (Tsaim > TcChangeToElec - dT_EvaporationStandard>Tsstart) # 该时间层的温度跨越电加热温度界限
			tmid=dt*(TcChangeToElec - dT_EvaporationStandard-Tsstart)/(Tsaim - Tsstart)
			flag_21, P1_21, P3_21, Pe_21= mode2(TcChangeToElec - dT_EvaporationStandard,Tsstart,tmid)
			flag_22, P1_22, P3_22, Pe_22= mode2(Tsaim,TcChangeToElec - dT_EvaporationStandard,dt-tmid)
			flag_2=flag_21 && flag_22
			P1_2=P1_21*tmid/dt+P1_22*(dt-tmid)/dt
			P3_2=P3_21*tmid/dt+P3_22*(dt-tmid)/dt
			Pe_2=Pe_21*tmid/dt+Pe_22*(dt-tmid)/dt
			Ptotal_2=P1_2+P3_2+Pe_2
		else# 该时间层初态温度高于电加热温度界线，切换模式1
			flag_2, P1_2, P3_2, Pe_2=mode1(dt)
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

	for (j,Tsaim) in enumerate(TsListEnd)
		for (i,Tsstart) in enumerate(TsListStart)
			# 先计算温度下降的功率,存在一边电加热一边开2号模式的情况
			if Tsstart<Tsmin || Tsaim>Tsmax
				C[i, j] = 9999.0
				P1Matrix[i, j] = 9999.0
				P3Matrix[i, j] = 9999.0
				PeMatrix[i, j] = 9999.0
			elseif Tsaim < Tsstart
				costValue,flag,P1,P2,Pe=powerCalculate_disCharge(Tsaim, Tsstart, dt)
				if flag
					C[i, j]=costValue
					P1Matrix[i, j] = P1
					P2Matrix[i, j] = P2
					PeMatrix[i, j] = Pe
				end
			elseif Tsaim == Tsstart
				temp = heatLoad-P1Only*COP1
				if temp <= PelecHeatMax
					P1Matrix[i, j] = P1Only
					PeMatrix[i, j] = temp
					C[i, j]=(P1Matrix[i, j]+PeMatrix[i, j])*dt
				end
				#=
				P1Matrix[i, j] = P1Only
				PeMatrix[i, j] = heatLoad-P1Only*COP1
				C[i, j]=(P1Matrix[i, j]+PeMatrix[i, j])*dt
				=#
			elseif Tsaim > Tsstart
				if Tsaim <= Tuse - dT_EvaporationStandard
					C1, flag, P1, P3, Pe = powerCalculate_lowTs(Tsaim, Tsstart, dt)
					if flag
						C[i, j] = C1
						P1Matrix[i, j] = P1
						P3Matrix[i, j] = P3
						PeMatrix[i, j] = Pe
					end
					#C[i, j] = flag ? C1 : 99
				elseif Tsstart < Tuse - dT_EvaporationStandard < Tsaim
					dt1 = (Tuse - dT_EvaporationStandard - Tsstart) / (Tsaim - Tsstart) * dt
					dt2 = dt - dt1
					C1, flag1, P11, P31, Pe1 = powerCalculate_lowTs(Tuse - dT_EvaporationStandard, Tsstart, dt1)
					C2, flag2, P12, P32, Pe2 = powerCalculate_highTs(Tsaim, Tuse - dT_EvaporationStandard, dt2)
					flag = flag1 && flag2
					if flag
						C[i, j] = C1 + C2
						P1Matrix[i, j] = (P11 * dt1 + P12 * dt2) / dt
						P3Matrix[i, j] = (P31 * dt1 + P32 * dt2) / dt
						PeMatrix[i, j] = (Pe1 * dt1 + Pe2 * dt2) / dt
					end
				elseif Tsstart >= Tuse - dT_EvaporationStandard
					C1, flag, P1, P3, Pe = powerCalculate_highTs(Tsaim, Tsstart, dt)
					if flag
						C[i, j] = C1
						P1Matrix[i, j] = P1
						P3Matrix[i, j] = P3
						PeMatrix[i, j] = Pe
					end
				end
			end
		end
	end

	return C, P1Matrix, P2Matrix, P3Matrix, PeMatrix
end

"""010"""
function getStateTransitionCost_SingleStep(
	::OO010;
	#C::Matrix, P1Matrix::Matrix, P2Matrix::Matrix, P3Matrix::Matrix, PeMatrix::Matrix;
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
	TCompressorIn::Real,# 压缩机吸气温度

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容

	# 设备运行约束
	PheatPumpMax::Real,# 热泵最大功率
    PelecHeatMax::Real,# 电锅炉最大功率
    PWaterCompressorMax::Real,# 水蒸气压缩机的最大功率
	#PheatPumpLowMax::Real,# 低温热泵的最大功率,一般来说，水蒸气压缩机和低温热泵同时运行时，由于低温热泵按照最劣工况设计，所以水蒸气压缩机会先触碰运行边界，不需要考虑低温热泵的边界问题
    Tsmin::Real,# 蓄热的最小温度
    
	# 求解参数
	TsListStart::Vector,# 状态参数高温蓄热温度起始值列表
	TsListEnd::Vector,# 状态参数高温蓄热温度结束值列表
	dt::Real = 1.0,# 时间步长
	Tsmax::Real = 220.0
)
	COP1 = COPOverlap(TWaste, Tuse)
	COP2 = COPWater
	COP2_design = COP2(TCompressorIn,Tuse)
	COP3 = COPOverlap
	heatLoadPumpMax=PWaterCompressorMax*COP2_design	# 热泵设计的供热功率

	nT = length(TsListStart)# 温度步数

	#C[i,j]表示温度从TsList[i]到TsList[j]时的最低功率；在负载、环境不变的情况下，C[i,j]是不变的
	C = fill(9999.0, nT, nT)# 状态转移矩阵
	P1Matrix = zeros(nT, nT)# 状态转移功率1参数
	P2Matrix = zeros(nT, nT)# 状态转移功率2参数
	P3Matrix = zeros(nT, nT)# 状态转移功率3参数
	PeMatrix = zeros(nT, nT)# 状态转移功率电加热参数
	
	# 只用热泵供热时的功率
	P1Only = min(heatLoad, heatLoadPumpMax)/ COP1# 只用热泵供热时的功率
	P1hOnly = P1Only*COP1/COP2_design	#只用热泵供热时的压缩机功率

	"""计算不同蓄热温度下蓄热温度降低的电度"""
	function powerCalculate_disCharge(Tsaim, Tsstart, dt)
		COP2value=COP2((Tsaim + Tsstart) / 2 - dT_EvaporationStandard, Tuse)        
        # 蓄热温度下降时，表现的功率
        PsDrop=cpm_h*(Tsstart-Tsaim)/dt
		Psout=cpm_h*(Tsstart-Tsaim)/dt
		P2=Psout / (COP2value-1)
		if (P2+Psout>heatLoad) || (P2 > PWaterCompressorMax) # 降不到这个温度
			return 9999.0,false,9999.0,9999.0,9999.0
		end
        #先看是不是温度到下限了
        if Tsaim == Tsmin
            P1 = min((heatLoad-P2-Psout) / COP1,PheatPumpMax)
            Pe = heatLoad-P2-Psout-P1 * COP1
            P21=min(PWaterCompressorMax,heatLoad/COP2value)  # 前半段蓄热供热最大P2功率
            Pe1=heatLoad-P21*COP2value
            P12 = min(heatLoad / COP1,P1Only)   # 后半段热泵供热的功率,P1Only保证了压缩机功率不超过设计功率，所以后面不用比较压缩机功率了
			Pe2 = heatLoad-P12 * COP1
			#P1h=P1*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P1l=P1-P1h
			#P12h=P12*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P12l=P12-P12h
            if (Pe1 <= PelecHeatMax) && (Pe2 <= PelecHeatMax) #&& (P1h<=PWaterCompressorMax) && (P12h<=PWaterCompressorMax) #&& (P1l <= PheatPumpLowMax) && (P12l <= PheatPumpLowMax)
				#println("放热2 Pe1=$Pe1 PelecHeatMax=$PelecHeatMax,Pe2=$Pe2,PelecHeatMax=$PelecHeatMax")
			    return (P1+P2+Pe)*dt,true,P1,P2,Pe
            else
				#println("放热3")
                return 9999.0,false,9999.0,9999.0,9999.0
            end
        end

    #蓄热温度的末态高于用热温度,放热全程热泵补热;或放热到蓄热的最低温度,先放热后热泵
		if (Tsaim-dT_EvaporationStandard>=Tuse)
			#P1 = min((heatLoad-P2-Psout) / COP1,P1Only)
            #=OOTest001=#
            P1=0
			Pe = heatLoad-P2-Psout-P1 * COP1
			#P1h=P1*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P1l=P1-P1h
            if (Pe <= PelecHeatMax) #&& (P1h<=PWaterCompressorMax) #&& (P1l <= PheatPumpLowMax)
			    return (P1+P2+Pe)*dt,true,P1,P2,Pe
            else
                return 9999.0,false,9999.0,9999.0,9999.0
            end
		elseif Tsstart-dT_EvaporationStandard<=Tuse#蓄热温度的初态低于用热温度，放热全程电热补热
            #热回收比例
            k=cp_cw/latentHeat*(Tuse-(Tsstart+Tsaim)/2)
			P1 = 0
			Pe = heatLoad*(1-k)-P2-Psout
            if Pe <= PelecHeatMax
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
		# 主要看水蒸气压缩机有没有余量
		P3h = min(PWaterCompressorMax-P1hOnly,cpm_h * (Tsaim - Tsstart) / dt / COP2_design / (1+cp_cw/latentHeat*(Tuse-0.5*(Tsstart+Tsaim))))

		P3 = P3h*COP2_design/COP1
		Pe = heatLoad+cpm_h / dt * (Tsaim - Tsstart) - (P3 + P1Only)*COP1
		
		#P3l = P3 - P3h
		if (Pe <= PelecHeatMax) #&& (P3h<=PWaterCompressorMax) #&& (P3l <= PheatPumpLowMax)
			#println("蓄热1 Pe=$Pe,PelecHeatMax=$PelecHeatMax")
			return (P1Only + P3 + Pe) * dt, true, P1Only, P3, Pe
		else
			#println("蓄热2 Pe=$Pe,PelecHeatMax=$PelecHeatMax")
			return 9999.0,false,9999.0,9999.0,9999.0
		end
	end

	"""计算蓄热温度相对较高时升温的电度情况(考虑蓄热温度超过热泵温度上限)"""
	function powerCalculate_highTs(Tsaim, Tsstart, dt)
		COP3value = COP3(TWaste, (Tsaim + Tsstart) / 2 + dT_EvaporationStandard)
		COPWatervalue = COP2(TCompressorIn, (Tsaim + Tsstart) / 2 + dT_EvaporationStandard)
		function mode1(dt1)# 模式1：只有电热加热给蓄热储热
			Pe = cpm_h / dt1 * (Tsaim - Tsstart)+heatLoad-P1Only*COP1
			if Pe <= PelecHeatMax
				return true, P1Only, 0.0, Pe
			else
				return false, 9999.0, 9999.0, 9999.0
			end
			
		end

		# 模式2：热泵也给蓄热储热（此时热泵COP按照蓄热温度计算）
		function mode2(Tsaim2, Tsstart2, dt2)
			# 如果水蒸气压缩机容量没有富余的话，应该优先用热泵为工厂供热，而不是为蓄热罐储热
			Psin=cpm_h / dt2 * (Tsaim2 - Tsstart2)#蓄热储热功率
			P3hin=min(Psin/COPWatervalue,PWaterCompressorMax)
            #=OOTest100=#
			#k=cp_cw/latentHeat * (Tsaim2 - Tsstart2)#闪蒸回收热比例
            k=0
			u_Pm=(heatLoad/COPWatervalue+PWaterCompressorMax)/(1+k)#线性规划求解的一个分类指标,用来描述用热负荷的大小
			# 先看闪蒸回收热是否够系统需求,如果够的话多余的热直接放掉
			if P3hin*k*COPWatervalue>=heatLoad
				return true, 0.0, P3hin, Psin-P3hin*COPWatervalue
			end
			
			if P3hin>=PWaterCompressorMax && u_Pm>=PWaterCompressorMax
				#蓄热负荷和用热负荷都大,任意进行热泵配置
				P1h=PWaterCompressorMax
				P3h=0.0
			elseif P3hin<PWaterCompressorMax && u_Pm>PWaterCompressorMax
				#蓄热负荷小，用热负荷大,蓄热占满，剩余的用热
				P3h=P3hin
				P1h=min(PWaterCompressorMax,(P3hin+heatLoad/COPWatervalue)/(1+k))-P3h
			elseif P3hin>PWaterCompressorMax && u_Pm<PWaterCompressorMax
				#蓄热负荷大，用热负荷小,蓄热按照供热程度占满
				P1h=0.0
				P3h=heatLoad/COPWatervalue/k
			else
				#蓄热负荷和用热负荷都小
				P3h=min(heatLoad/COPWatervalue/k,P3hin)
				P1h=min(heatLoad/COPWatervalue/k,(P3hin+heatLoad/COPWatervalue)/(1+k))-P3h
			end

			Pe = cpm_h / dt2 * (Tsaim2 - Tsstart2) + heatLoad - (P1h+P3h) * COPWatervalue
			if Pe <= PelecHeatMax
				# flag P1 P3 Pe
				return true, P1h*COPWatervalue/COP1, P3h*COPWatervalue/COP1, Pe
			else
				return false, 9999.0, 9999.0, 9999.0
			end
		end
		#COP3value = COP3(TWaste, (min(Tsaim,TcChangeToElec - dT_EvaporationStandard) + Tsstart) / 2 + dT_EvaporationStandard)
		
		Ptotal = 9999.0
		P1res = 0.0
		P3res = 0.0
		Peres = 0.0
		flag = false
		# 模式1：电热加热
		flag, P1res, P3res, Peres= mode1(dt)
		Ptotal=P1res+P3res+Peres

		# 模式2：热泵加电热,又分两种模式：末态温度不大于电加热温度界限与末态温度大于电加热温度界线
		if (Tsaim <= TcChangeToElec - dT_EvaporationStandard)	# 该时间层的末态温度小于电加热温度界限
			flag_2, P1_2, P3_2, Pe_2= mode2(Tsaim,Tsstart,dt)
			Ptotal_2=P1_2+P3_2+Pe_2
		elseif (Tsaim > TcChangeToElec - dT_EvaporationStandard>Tsstart) # 该时间层的温度跨越电加热温度界限
			tmid=dt*(TcChangeToElec - dT_EvaporationStandard-Tsstart)/(Tsaim - Tsstart)
			flag_21, P1_21, P3_21, Pe_21= mode2(TcChangeToElec - dT_EvaporationStandard,Tsstart,tmid)
			flag_22, P1_22, P3_22, Pe_22= mode2(Tsaim,TcChangeToElec - dT_EvaporationStandard,dt-tmid)
			flag_2=flag_21 && flag_22
			P1_2=P1_21*tmid/dt+P1_22*(dt-tmid)/dt
			P3_2=P3_21*tmid/dt+P3_22*(dt-tmid)/dt
			Pe_2=Pe_21*tmid/dt+Pe_22*(dt-tmid)/dt
			Ptotal_2=P1_2+P3_2+Pe_2
		else# 该时间层初态温度高于电加热温度界线，切换模式1
			flag_2, P1_2, P3_2, Pe_2=mode1(dt)
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

	for (j,Tsaim) in enumerate(TsListEnd)
		for (i,Tsstart) in enumerate(TsListStart)
			# 先计算温度下降的功率,存在一边电加热一边开2号模式的情况
			if Tsstart<Tsmin || Tsaim>Tsmax
				C[i, j] = 9999.0
				P1Matrix[i, j] = 9999.0
				P3Matrix[i, j] = 9999.0
				PeMatrix[i, j] = 9999.0
			elseif Tsaim < Tsstart
				costValue,flag,P1,P2,Pe=powerCalculate_disCharge(Tsaim, Tsstart, dt)
				if flag
					C[i, j]=costValue
					P1Matrix[i, j] = P1
					P2Matrix[i, j] = P2
					PeMatrix[i, j] = Pe
				end
			elseif Tsaim == Tsstart
				temp = heatLoad-P1Only*COP1
				if temp <= PelecHeatMax
					P1Matrix[i, j] = P1Only
					PeMatrix[i, j] = temp
					C[i, j]=(P1Matrix[i, j]+PeMatrix[i, j])*dt
				end
				#=
				P1Matrix[i, j] = P1Only
				PeMatrix[i, j] = heatLoad-P1Only*COP1
				C[i, j]=(P1Matrix[i, j]+PeMatrix[i, j])*dt
				=#
			elseif Tsaim > Tsstart
				if Tsaim <= Tuse - dT_EvaporationStandard
					C1, flag, P1, P3, Pe = powerCalculate_lowTs(Tsaim, Tsstart, dt)
					if flag
						C[i, j] = C1
						P1Matrix[i, j] = P1
						P3Matrix[i, j] = P3
						PeMatrix[i, j] = Pe
					end
					#C[i, j] = flag ? C1 : 99
				elseif Tsstart < Tuse - dT_EvaporationStandard < Tsaim
					dt1 = (Tuse - dT_EvaporationStandard - Tsstart) / (Tsaim - Tsstart) * dt
					dt2 = dt - dt1
					C1, flag1, P11, P31, Pe1 = powerCalculate_lowTs(Tuse - dT_EvaporationStandard, Tsstart, dt1)
					C2, flag2, P12, P32, Pe2 = powerCalculate_highTs(Tsaim, Tuse - dT_EvaporationStandard, dt2)
					flag = flag1 && flag2
					if flag
						C[i, j] = C1 + C2
						P1Matrix[i, j] = (P11 * dt1 + P12 * dt2) / dt
						P3Matrix[i, j] = (P31 * dt1 + P32 * dt2) / dt
						PeMatrix[i, j] = (Pe1 * dt1 + Pe2 * dt2) / dt
					end
				elseif Tsstart >= Tuse - dT_EvaporationStandard
					C1, flag, P1, P3, Pe = powerCalculate_highTs(Tsaim, Tsstart, dt)
					if flag
						C[i, j] = C1
						P1Matrix[i, j] = P1
						P3Matrix[i, j] = P3
						PeMatrix[i, j] = Pe
					end
				end
			end
		end
	end

	return C, P1Matrix, P2Matrix, P3Matrix, PeMatrix
end

"""011"""
function getStateTransitionCost_SingleStep(
	::OO011;
	#C::Matrix, P1Matrix::Matrix, P2Matrix::Matrix, P3Matrix::Matrix, PeMatrix::Matrix;
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
	TCompressorIn::Real,# 压缩机吸气温度

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容

	# 设备运行约束
	PheatPumpMax::Real,# 热泵最大功率
    PelecHeatMax::Real,# 电锅炉最大功率
    PWaterCompressorMax::Real,# 水蒸气压缩机的最大功率
	#PheatPumpLowMax::Real,# 低温热泵的最大功率,一般来说，水蒸气压缩机和低温热泵同时运行时，由于低温热泵按照最劣工况设计，所以水蒸气压缩机会先触碰运行边界，不需要考虑低温热泵的边界问题
    Tsmin::Real,# 蓄热的最小温度
    
	# 求解参数
	TsListStart::Vector,# 状态参数高温蓄热温度起始值列表
	TsListEnd::Vector,# 状态参数高温蓄热温度结束值列表
	dt::Real = 1.0,# 时间步长
	Tsmax::Real = 220.0
)
	COP1 = COPOverlap(TWaste, Tuse)
	COP2 = COPWater
	COP2_design = COP2(TCompressorIn,Tuse)
	COP3 = COPOverlap
	heatLoadPumpMax=PWaterCompressorMax*COP2_design	# 热泵设计的供热功率

	nT = length(TsListStart)# 温度步数

	#C[i,j]表示温度从TsList[i]到TsList[j]时的最低功率；在负载、环境不变的情况下，C[i,j]是不变的
	C = fill(9999.0, nT, nT)# 状态转移矩阵
	P1Matrix = zeros(nT, nT)# 状态转移功率1参数
	P2Matrix = zeros(nT, nT)# 状态转移功率2参数
	P3Matrix = zeros(nT, nT)# 状态转移功率3参数
	PeMatrix = zeros(nT, nT)# 状态转移功率电加热参数
	
	# 只用热泵供热时的功率
	P1Only = min(heatLoad, heatLoadPumpMax)/ COP1# 只用热泵供热时的功率
	P1hOnly = P1Only*COP1/COP2_design	#只用热泵供热时的压缩机功率

	"""计算不同蓄热温度下蓄热温度降低的电度"""
	function powerCalculate_disCharge(Tsaim, Tsstart, dt)
		COP2value=COP2((Tsaim + Tsstart) / 2 - dT_EvaporationStandard, Tuse)        
        # 蓄热温度下降时，表现的功率
        #=以下是系统为OOTest010的修改=#
        PsDrop=cpm_h*(Tsstart-Tsaim)/dt
		Psout=cpm_h*(Tsstart-Tsaim)/dt
		P2=Psout / (COP2value-1)
		if (P2+Psout>heatLoad) || (P2 > PWaterCompressorMax) # 降不到这个温度
			return 9999.0,false,9999.0,9999.0,9999.0
		end
        #先看是不是温度到下限了
        if Tsaim == Tsmin
            P1 = min((heatLoad-P2-Psout) / COP1,PheatPumpMax)
            Pe = heatLoad-P2-Psout-P1 * COP1
            P21=min(PWaterCompressorMax,heatLoad/COP2value)  # 前半段蓄热供热最大P2功率
            Pe1=heatLoad-P21*COP2value
            P12 = min(heatLoad / COP1,P1Only)   # 后半段热泵供热的功率,P1Only保证了压缩机功率不超过设计功率，所以后面不用比较压缩机功率了
			Pe2 = heatLoad-P12 * COP1
			#P1h=P1*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P1l=P1-P1h
			#P12h=P12*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P12l=P12-P12h
            if (Pe1 <= PelecHeatMax) && (Pe2 <= PelecHeatMax) #&& (P1h<=PWaterCompressorMax) && (P12h<=PWaterCompressorMax) #&& (P1l <= PheatPumpLowMax) && (P12l <= PheatPumpLowMax)
				#println("放热2 Pe1=$Pe1 PelecHeatMax=$PelecHeatMax,Pe2=$Pe2,PelecHeatMax=$PelecHeatMax")
			    return (P1+P2+Pe)*dt,true,P1,P2,Pe
            else
				#println("放热3")
                return 9999.0,false,9999.0,9999.0,9999.0
            end
        end

    #蓄热温度的末态高于用热温度,放热全程热泵补热;或放热到蓄热的最低温度,先放热后热泵
		if (Tsaim-dT_EvaporationStandard>=Tuse)
			P1 = min((heatLoad-P2-Psout) / COP1,P1Only)
            #=OOTest001=#
            #P1=0
			Pe = heatLoad-P2-Psout-P1 * COP1
			#P1h=P1*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P1l=P1-P1h
            if (Pe <= PelecHeatMax) #&& (P1h<=PWaterCompressorMax) #&& (P1l <= PheatPumpLowMax)
			    return (P1+P2+Pe)*dt,true,P1,P2,Pe
            else
                return 9999.0,false,9999.0,9999.0,9999.0
            end
		elseif Tsstart-dT_EvaporationStandard<=Tuse#蓄热温度的初态低于用热温度，放热全程电热补热
            #热回收比例
            k=cp_cw/latentHeat*(Tuse-(Tsstart+Tsaim)/2)
			P1 = 0
			Pe = heatLoad*(1-k)-P2-Psout
            if Pe <= PelecHeatMax
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
		# 主要看水蒸气压缩机有没有余量
		P3h = min(PWaterCompressorMax-P1hOnly,cpm_h * (Tsaim - Tsstart) / dt / COP2_design / (1+cp_cw/latentHeat*(Tuse-0.5*(Tsstart+Tsaim))))

		P3 = P3h*COP2_design/COP1
		Pe = heatLoad+cpm_h / dt * (Tsaim - Tsstart) - (P3 + P1Only)*COP1
		
		#P3l = P3 - P3h
		if (Pe <= PelecHeatMax) #&& (P3h<=PWaterCompressorMax) #&& (P3l <= PheatPumpLowMax)
			#println("蓄热1 Pe=$Pe,PelecHeatMax=$PelecHeatMax")
			return (P1Only + P3 + Pe) * dt, true, P1Only, P3, Pe
		else
			#println("蓄热2 Pe=$Pe,PelecHeatMax=$PelecHeatMax")
			return 9999.0,false,9999.0,9999.0,9999.0
		end
	end

	"""计算蓄热温度相对较高时升温的电度情况(考虑蓄热温度超过热泵温度上限)"""
	function powerCalculate_highTs(Tsaim, Tsstart, dt)
		COP3value = COP3(TWaste, (Tsaim + Tsstart) / 2 + dT_EvaporationStandard)
		COPWatervalue = COP2(TCompressorIn, (Tsaim + Tsstart) / 2 + dT_EvaporationStandard)
		function mode1(dt1)# 模式1：只有电热加热给蓄热储热
			Pe = cpm_h / dt1 * (Tsaim - Tsstart)+heatLoad-P1Only*COP1
			if Pe <= PelecHeatMax
				return true, P1Only, 0.0, Pe
			else
				return false, 9999.0, 9999.0, 9999.0
			end
			
		end

		# 模式2：热泵也给蓄热储热（此时热泵COP按照蓄热温度计算）
		function mode2(Tsaim2, Tsstart2, dt2)
			# 如果水蒸气压缩机容量没有富余的话，应该优先用热泵为工厂供热，而不是为蓄热罐储热
			Psin=cpm_h / dt2 * (Tsaim2 - Tsstart2)#蓄热储热功率
			P3hin=min(Psin/COPWatervalue,PWaterCompressorMax)
            #=OOTest100=#
			#k=cp_cw/latentHeat * (Tsaim2 - Tsstart2)#闪蒸回收热比例
            k=0
			u_Pm=(heatLoad/COPWatervalue+PWaterCompressorMax)/(1+k)#线性规划求解的一个分类指标,用来描述用热负荷的大小
			# 先看闪蒸回收热是否够系统需求,如果够的话多余的热直接放掉
			if P3hin*k*COPWatervalue>=heatLoad
				return true, 0.0, P3hin, Psin-P3hin*COPWatervalue
			end
			
			if P3hin>=PWaterCompressorMax && u_Pm>=PWaterCompressorMax
				#蓄热负荷和用热负荷都大,任意进行热泵配置
				P1h=PWaterCompressorMax
				P3h=0.0
			elseif P3hin<PWaterCompressorMax && u_Pm>PWaterCompressorMax
				#蓄热负荷小，用热负荷大,蓄热占满，剩余的用热
				P3h=P3hin
				P1h=min(PWaterCompressorMax,(P3hin+heatLoad/COPWatervalue)/(1+k))-P3h
			elseif P3hin>PWaterCompressorMax && u_Pm<PWaterCompressorMax
				#蓄热负荷大，用热负荷小,蓄热按照供热程度占满
				P1h=0.0
				P3h=heatLoad/COPWatervalue/k
			else
				#蓄热负荷和用热负荷都小
				P3h=min(heatLoad/COPWatervalue/k,P3hin)
				P1h=min(heatLoad/COPWatervalue/k,(P3hin+heatLoad/COPWatervalue)/(1+k))-P3h
			end

			Pe = cpm_h / dt2 * (Tsaim2 - Tsstart2) + heatLoad - (P1h+P3h) * COPWatervalue
			if Pe <= PelecHeatMax
				# flag P1 P3 Pe
				return true, P1h*COPWatervalue/COP1, P3h*COPWatervalue/COP1, Pe
			else
				return false, 9999.0, 9999.0, 9999.0
			end
		end
		#COP3value = COP3(TWaste, (min(Tsaim,TcChangeToElec - dT_EvaporationStandard) + Tsstart) / 2 + dT_EvaporationStandard)
		
		Ptotal = 9999.0
		P1res = 0.0
		P3res = 0.0
		Peres = 0.0
		flag = false
		# 模式1：电热加热
		flag, P1res, P3res, Peres= mode1(dt)
		Ptotal=P1res+P3res+Peres

		# 模式2：热泵加电热,又分两种模式：末态温度不大于电加热温度界限与末态温度大于电加热温度界线
		if (Tsaim <= TcChangeToElec - dT_EvaporationStandard)	# 该时间层的末态温度小于电加热温度界限
			flag_2, P1_2, P3_2, Pe_2= mode2(Tsaim,Tsstart,dt)
			Ptotal_2=P1_2+P3_2+Pe_2
		elseif (Tsaim > TcChangeToElec - dT_EvaporationStandard>Tsstart) # 该时间层的温度跨越电加热温度界限
			tmid=dt*(TcChangeToElec - dT_EvaporationStandard-Tsstart)/(Tsaim - Tsstart)
			flag_21, P1_21, P3_21, Pe_21= mode2(TcChangeToElec - dT_EvaporationStandard,Tsstart,tmid)
			flag_22, P1_22, P3_22, Pe_22= mode2(Tsaim,TcChangeToElec - dT_EvaporationStandard,dt-tmid)
			flag_2=flag_21 && flag_22
			P1_2=P1_21*tmid/dt+P1_22*(dt-tmid)/dt
			P3_2=P3_21*tmid/dt+P3_22*(dt-tmid)/dt
			Pe_2=Pe_21*tmid/dt+Pe_22*(dt-tmid)/dt
			Ptotal_2=P1_2+P3_2+Pe_2
		else# 该时间层初态温度高于电加热温度界线，切换模式1
			flag_2, P1_2, P3_2, Pe_2=mode1(dt)
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

	for (j,Tsaim) in enumerate(TsListEnd)
		for (i,Tsstart) in enumerate(TsListStart)
			# 先计算温度下降的功率,存在一边电加热一边开2号模式的情况
			if Tsstart<Tsmin || Tsaim>Tsmax
				C[i, j] = 9999.0
				P1Matrix[i, j] = 9999.0
				P3Matrix[i, j] = 9999.0
				PeMatrix[i, j] = 9999.0
			elseif Tsaim < Tsstart
				costValue,flag,P1,P2,Pe=powerCalculate_disCharge(Tsaim, Tsstart, dt)
				if flag
					C[i, j]=costValue
					P1Matrix[i, j] = P1
					P2Matrix[i, j] = P2
					PeMatrix[i, j] = Pe
				end
			elseif Tsaim == Tsstart
				temp = heatLoad-P1Only*COP1
				if temp <= PelecHeatMax
					P1Matrix[i, j] = P1Only
					PeMatrix[i, j] = temp
					C[i, j]=(P1Matrix[i, j]+PeMatrix[i, j])*dt
				end
				#=
				P1Matrix[i, j] = P1Only
				PeMatrix[i, j] = heatLoad-P1Only*COP1
				C[i, j]=(P1Matrix[i, j]+PeMatrix[i, j])*dt
				=#
			elseif Tsaim > Tsstart
				if Tsaim <= Tuse - dT_EvaporationStandard
					C1, flag, P1, P3, Pe = powerCalculate_lowTs(Tsaim, Tsstart, dt)
					if flag
						C[i, j] = C1
						P1Matrix[i, j] = P1
						P3Matrix[i, j] = P3
						PeMatrix[i, j] = Pe
					end
					#C[i, j] = flag ? C1 : 99
				elseif Tsstart < Tuse - dT_EvaporationStandard < Tsaim
					dt1 = (Tuse - dT_EvaporationStandard - Tsstart) / (Tsaim - Tsstart) * dt
					dt2 = dt - dt1
					C1, flag1, P11, P31, Pe1 = powerCalculate_lowTs(Tuse - dT_EvaporationStandard, Tsstart, dt1)
					C2, flag2, P12, P32, Pe2 = powerCalculate_highTs(Tsaim, Tuse - dT_EvaporationStandard, dt2)
					flag = flag1 && flag2
					if flag
						C[i, j] = C1 + C2
						P1Matrix[i, j] = (P11 * dt1 + P12 * dt2) / dt
						P3Matrix[i, j] = (P31 * dt1 + P32 * dt2) / dt
						PeMatrix[i, j] = (Pe1 * dt1 + Pe2 * dt2) / dt
					end
				elseif Tsstart >= Tuse - dT_EvaporationStandard
					C1, flag, P1, P3, Pe = powerCalculate_highTs(Tsaim, Tsstart, dt)
					if flag
						C[i, j] = C1
						P1Matrix[i, j] = P1
						P3Matrix[i, j] = P3
						PeMatrix[i, j] = Pe
					end
				end
			end
		end
	end

	return C, P1Matrix, P2Matrix, P3Matrix, PeMatrix
end

"""100"""
function getStateTransitionCost_SingleStep(
	::OO100;
	#C::Matrix, P1Matrix::Matrix, P2Matrix::Matrix, P3Matrix::Matrix, PeMatrix::Matrix;
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
	TCompressorIn::Real,# 压缩机吸气温度

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容

	# 设备运行约束
	PheatPumpMax::Real,# 热泵最大功率
    PelecHeatMax::Real,# 电锅炉最大功率
    PWaterCompressorMax::Real,# 水蒸气压缩机的最大功率
	#PheatPumpLowMax::Real,# 低温热泵的最大功率,一般来说，水蒸气压缩机和低温热泵同时运行时，由于低温热泵按照最劣工况设计，所以水蒸气压缩机会先触碰运行边界，不需要考虑低温热泵的边界问题
    Tsmin::Real,# 蓄热的最小温度
    
	# 求解参数
	TsListStart::Vector,# 状态参数高温蓄热温度起始值列表
	TsListEnd::Vector,# 状态参数高温蓄热温度结束值列表
	dt::Real = 1.0,# 时间步长
	Tsmax::Real = 220.0
)
	COP1 = COPOverlap(TWaste, Tuse)
	COP2 = COPWater
	COP2_design = COP2(TCompressorIn,Tuse)
	COP3 = COPOverlap
	heatLoadPumpMax=PWaterCompressorMax*COP2_design	# 热泵设计的供热功率

	nT = length(TsListStart)# 温度步数

	#C[i,j]表示温度从TsList[i]到TsList[j]时的最低功率；在负载、环境不变的情况下，C[i,j]是不变的
	C = fill(9999.0, nT, nT)# 状态转移矩阵
	P1Matrix = zeros(nT, nT)# 状态转移功率1参数
	P2Matrix = zeros(nT, nT)# 状态转移功率2参数
	P3Matrix = zeros(nT, nT)# 状态转移功率3参数
	PeMatrix = zeros(nT, nT)# 状态转移功率电加热参数

	
	# 只用热泵供热时的功率
	P1Only = min(heatLoad, heatLoadPumpMax)/ COP1# 只用热泵供热时的功率
	P1hOnly = P1Only*COP1/COP2_design	#只用热泵供热时的压缩机功率

	"""计算不同蓄热温度下蓄热温度降低的电度"""
	function powerCalculate_disCharge(Tsaim, Tsstart, dt)
		COP2value=COP2((Tsaim + Tsstart) / 2 - dT_EvaporationStandard, Tuse)
		Psout=cpm_h*(Tsstart-Tsaim)/dt
		P2=Psout / (COP2value-1)
		if (P2+Psout>heatLoad) || (P2 > PWaterCompressorMax) # 降不到这个温度
			return 9999.0,false,9999.0,9999.0,9999.0
		end
        #先看是不是温度到下限了
        if Tsaim == Tsmin
            P1 = min((heatLoad-P2-Psout) / COP1,PheatPumpMax)
            Pe = heatLoad-P2-Psout-P1 * COP1
            P21=min(PWaterCompressorMax,heatLoad/COP2value)  # 前半段蓄热供热最大P2功率
            Pe1=heatLoad-P21*COP2value
            P12 = min(heatLoad / COP1,P1Only)   # 后半段热泵供热的功率,P1Only保证了压缩机功率不超过设计功率，所以后面不用比较压缩机功率了
			Pe2 = heatLoad-P12 * COP1
			#P1h=P1*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P1l=P1-P1h
			#P12h=P12*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P12l=P12-P12h
            if (Pe1 <= PelecHeatMax) && (Pe2 <= PelecHeatMax) #&& (P1h<=PWaterCompressorMax) && (P12h<=PWaterCompressorMax) #&& (P1l <= PheatPumpLowMax) && (P12l <= PheatPumpLowMax)
				#println("放热2 Pe1=$Pe1 PelecHeatMax=$PelecHeatMax,Pe2=$Pe2,PelecHeatMax=$PelecHeatMax")
			    return (P1+P2+Pe)*dt,true,P1,P2,Pe
            else
				#println("放热3")
                return 9999.0,false,9999.0,9999.0,9999.0
            end
        end

        #蓄热温度的末态高于用热温度,放热全程热泵补热;或放热到蓄热的最低温度,先放热后热泵
		if (Tsaim-dT_EvaporationStandard>=Tuse)
			#P1 = min((heatLoad-P2-Psout) / COP1,P1Only)
            #=OOTest001=#
            P1=0
			Pe = heatLoad-P2-Psout-P1 * COP1
			#P1h=P1*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P1l=P1-P1h
            if (Pe <= PelecHeatMax) #&& (P1h<=PWaterCompressorMax) #&& (P1l <= PheatPumpLowMax)
			    return (P1+P2+Pe)*dt,true,P1,P2,Pe
            else
                return 9999.0,false,9999.0,9999.0,9999.0
            end
		elseif Tsstart-dT_EvaporationStandard<=Tuse#蓄热温度的初态低于用热温度，放热全程电热补热
			P1 = 0
			Pe = heatLoad-P2-Psout
            if Pe <= PelecHeatMax
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
		# 主要看水蒸气压缩机有没有余量
		P3h = min(PWaterCompressorMax-P1hOnly,cpm_h * (Tsaim - Tsstart) / dt / COP2_design / (1+cp_cw/latentHeat*(Tuse-0.5*(Tsstart+Tsaim))))

		P3 = P3h*COP2_design/COP1
		Pe = heatLoad+cpm_h / dt * (Tsaim - Tsstart) - (P3 + P1Only)*COP1
		
		#P3l = P3 - P3h
		if (Pe <= PelecHeatMax) #&& (P3h<=PWaterCompressorMax) #&& (P3l <= PheatPumpLowMax)
			#println("蓄热1 Pe=$Pe,PelecHeatMax=$PelecHeatMax")
			return (P1Only + P3 + Pe) * dt, true, P1Only, P3, Pe
		else
			#println("蓄热2 Pe=$Pe,PelecHeatMax=$PelecHeatMax")
			return 9999.0,false,9999.0,9999.0,9999.0
		end
	end

	"""计算蓄热温度相对较高时升温的电度情况(考虑蓄热温度超过热泵温度上限)"""
	function powerCalculate_highTs(Tsaim, Tsstart, dt)
		COP3value = COP3(TWaste, (Tsaim + Tsstart) / 2 + dT_EvaporationStandard)
		COPWatervalue = COP2(TCompressorIn, (Tsaim + Tsstart) / 2 + dT_EvaporationStandard)
		function mode1(dt1)# 模式1：只有电热加热给蓄热储热
			Pe = cpm_h / dt1 * (Tsaim - Tsstart)+heatLoad-P1Only*COP1
			if Pe <= PelecHeatMax
				return true, P1Only, 0.0, Pe
			else
				return false, 9999.0, 9999.0, 9999.0
			end
			
		end

		# 模式2：热泵也给蓄热储热（此时热泵COP按照蓄热温度计算）
		function mode2(Tsaim2, Tsstart2, dt2)
			# 如果水蒸气压缩机容量没有富余的话，应该优先用热泵为工厂供热，而不是为蓄热罐储热
			Psin=cpm_h / dt2 * (Tsaim2 - Tsstart2)#蓄热储热功率
			P3hin=min(Psin/COPWatervalue,PWaterCompressorMax)
            #=OOTest100=#
			k=cp_cw/latentHeat * (Tsaim2 - Tsstart2)#闪蒸回收热比例
            #k=0
			u_Pm=(heatLoad/COPWatervalue+PWaterCompressorMax)/(1+k)#线性规划求解的一个分类指标,用来描述用热负荷的大小
			# 先看闪蒸回收热是否够系统需求,如果够的话多余的热直接放掉
			if P3hin*k*COPWatervalue>=heatLoad
				return true, 0.0, P3hin, Psin-P3hin*COPWatervalue
			end
			
			if P3hin>=PWaterCompressorMax && u_Pm>=PWaterCompressorMax
				#蓄热负荷和用热负荷都大,任意进行热泵配置
				P1h=PWaterCompressorMax
				P3h=0.0
			elseif P3hin<PWaterCompressorMax && u_Pm>PWaterCompressorMax
				#蓄热负荷小，用热负荷大,蓄热占满，剩余的用热
				P3h=P3hin
				P1h=min(PWaterCompressorMax,(P3hin+heatLoad/COPWatervalue)/(1+k))-P3h
			elseif P3hin>PWaterCompressorMax && u_Pm<PWaterCompressorMax
				#蓄热负荷大，用热负荷小,蓄热按照供热程度占满
				P1h=0.0
				P3h=heatLoad/COPWatervalue/k
			else
				#蓄热负荷和用热负荷都小
				P3h=min(heatLoad/COPWatervalue/k,P3hin)
				P1h=min(heatLoad/COPWatervalue/k,(P3hin+heatLoad/COPWatervalue)/(1+k))-P3h
			end

			Pe = cpm_h / dt2 * (Tsaim2 - Tsstart2) + heatLoad - (P1h+P3h) * COPWatervalue
			if Pe <= PelecHeatMax
				# flag P1 P3 Pe
				return true, P1h*COPWatervalue/COP1, P3h*COPWatervalue/COP1, Pe
			else
				return false, 9999.0, 9999.0, 9999.0
			end
		end
		#COP3value = COP3(TWaste, (min(Tsaim,TcChangeToElec - dT_EvaporationStandard) + Tsstart) / 2 + dT_EvaporationStandard)
		
		Ptotal = 9999.0
		P1res = 0.0
		P3res = 0.0
		Peres = 0.0
		flag = false
		# 模式1：电热加热
		flag, P1res, P3res, Peres= mode1(dt)
		Ptotal=P1res+P3res+Peres

		# 模式2：热泵加电热,又分两种模式：末态温度不大于电加热温度界限与末态温度大于电加热温度界线
		if (Tsaim <= TcChangeToElec - dT_EvaporationStandard)	# 该时间层的末态温度小于电加热温度界限
			flag_2, P1_2, P3_2, Pe_2= mode2(Tsaim,Tsstart,dt)
			Ptotal_2=P1_2+P3_2+Pe_2
		elseif (Tsaim > TcChangeToElec - dT_EvaporationStandard>Tsstart) # 该时间层的温度跨越电加热温度界限
			tmid=dt*(TcChangeToElec - dT_EvaporationStandard-Tsstart)/(Tsaim - Tsstart)
			flag_21, P1_21, P3_21, Pe_21= mode2(TcChangeToElec - dT_EvaporationStandard,Tsstart,tmid)
			flag_22, P1_22, P3_22, Pe_22= mode2(Tsaim,TcChangeToElec - dT_EvaporationStandard,dt-tmid)
			flag_2=flag_21 && flag_22
			P1_2=P1_21*tmid/dt+P1_22*(dt-tmid)/dt
			P3_2=P3_21*tmid/dt+P3_22*(dt-tmid)/dt
			Pe_2=Pe_21*tmid/dt+Pe_22*(dt-tmid)/dt
			Ptotal_2=P1_2+P3_2+Pe_2
		else# 该时间层初态温度高于电加热温度界线，切换模式1
			flag_2, P1_2, P3_2, Pe_2=mode1(dt)
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

	for (j,Tsaim) in enumerate(TsListEnd)
		for (i,Tsstart) in enumerate(TsListStart)
			# 先计算温度下降的功率,存在一边电加热一边开2号模式的情况
			if Tsstart<Tsmin || Tsaim>Tsmax
				C[i, j] = 9999.0
				P1Matrix[i, j] = 9999.0
				P3Matrix[i, j] = 9999.0
				PeMatrix[i, j] = 9999.0
			elseif Tsaim < Tsstart
				costValue,flag,P1,P2,Pe=powerCalculate_disCharge(Tsaim, Tsstart, dt)
				if flag
					C[i, j]=costValue
					P1Matrix[i, j] = P1
					P2Matrix[i, j] = P2
					PeMatrix[i, j] = Pe
				end
			elseif Tsaim == Tsstart
				temp = heatLoad-P1Only*COP1
				if temp <= PelecHeatMax
					P1Matrix[i, j] = P1Only
					PeMatrix[i, j] = temp
					C[i, j]=(P1Matrix[i, j]+PeMatrix[i, j])*dt
				end
				#=
				P1Matrix[i, j] = P1Only
				PeMatrix[i, j] = heatLoad-P1Only*COP1
				C[i, j]=(P1Matrix[i, j]+PeMatrix[i, j])*dt
				=#
			elseif Tsaim > Tsstart
				if Tsaim <= Tuse - dT_EvaporationStandard
					C1, flag, P1, P3, Pe = powerCalculate_lowTs(Tsaim, Tsstart, dt)
					if flag
						C[i, j] = C1
						P1Matrix[i, j] = P1
						P3Matrix[i, j] = P3
						PeMatrix[i, j] = Pe
					end
					#C[i, j] = flag ? C1 : 99
				elseif Tsstart < Tuse - dT_EvaporationStandard < Tsaim
					dt1 = (Tuse - dT_EvaporationStandard - Tsstart) / (Tsaim - Tsstart) * dt
					dt2 = dt - dt1
					C1, flag1, P11, P31, Pe1 = powerCalculate_lowTs(Tuse - dT_EvaporationStandard, Tsstart, dt1)
					C2, flag2, P12, P32, Pe2 = powerCalculate_highTs(Tsaim, Tuse - dT_EvaporationStandard, dt2)
					flag = flag1 && flag2
					if flag
						C[i, j] = C1 + C2
						P1Matrix[i, j] = (P11 * dt1 + P12 * dt2) / dt
						P3Matrix[i, j] = (P31 * dt1 + P32 * dt2) / dt
						PeMatrix[i, j] = (Pe1 * dt1 + Pe2 * dt2) / dt
					end
				elseif Tsstart >= Tuse - dT_EvaporationStandard
					C1, flag, P1, P3, Pe = powerCalculate_highTs(Tsaim, Tsstart, dt)
					if flag
						C[i, j] = C1
						P1Matrix[i, j] = P1
						P3Matrix[i, j] = P3
						PeMatrix[i, j] = Pe
					end
				end
			end
		end
	end

	return C, P1Matrix, P2Matrix, P3Matrix, PeMatrix
end

"""101"""
function getStateTransitionCost_SingleStep(
	::OO101;
	#C::Matrix, P1Matrix::Matrix, P2Matrix::Matrix, P3Matrix::Matrix, PeMatrix::Matrix;
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
	TCompressorIn::Real,# 压缩机吸气温度

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容

	# 设备运行约束
	PheatPumpMax::Real,# 热泵最大功率
    PelecHeatMax::Real,# 电锅炉最大功率
    PWaterCompressorMax::Real,# 水蒸气压缩机的最大功率
	#PheatPumpLowMax::Real,# 低温热泵的最大功率,一般来说，水蒸气压缩机和低温热泵同时运行时，由于低温热泵按照最劣工况设计，所以水蒸气压缩机会先触碰运行边界，不需要考虑低温热泵的边界问题
    Tsmin::Real,# 蓄热的最小温度
    
	# 求解参数
	TsListStart::Vector,# 状态参数高温蓄热温度起始值列表
	TsListEnd::Vector,# 状态参数高温蓄热温度结束值列表
	dt::Real = 1.0,# 时间步长
	Tsmax::Real = 220.0
)
	COP1 = COPOverlap(TWaste, Tuse)
	COP2 = COPWater
	COP2_design = COP2(TCompressorIn,Tuse)
	COP3 = COPOverlap
	heatLoadPumpMax=PWaterCompressorMax*COP2_design	# 热泵设计的供热功率

	nT = length(TsListStart)# 温度步数

	#C[i,j]表示温度从TsList[i]到TsList[j]时的最低功率；在负载、环境不变的情况下，C[i,j]是不变的
	C = fill(9999.0, nT, nT)# 状态转移矩阵
	P1Matrix = zeros(nT, nT)# 状态转移功率1参数
	P2Matrix = zeros(nT, nT)# 状态转移功率2参数
	P3Matrix = zeros(nT, nT)# 状态转移功率3参数
	PeMatrix = zeros(nT, nT)# 状态转移功率电加热参数

	
	# 只用热泵供热时的功率
	P1Only = min(heatLoad, heatLoadPumpMax)/ COP1# 只用热泵供热时的功率
	P1hOnly = P1Only*COP1/COP2_design	#只用热泵供热时的压缩机功率

	"""计算不同蓄热温度下蓄热温度降低的电度"""
	function powerCalculate_disCharge(Tsaim, Tsstart, dt)
		COP2value=COP2((Tsaim + Tsstart) / 2 - dT_EvaporationStandard, Tuse)
		Psout=cpm_h*(Tsstart-Tsaim)/dt
		P2=Psout / (COP2value-1)
		if (P2+Psout>heatLoad) || (P2 > PWaterCompressorMax) # 降不到这个温度
			return 9999.0,false,9999.0,9999.0,9999.0
		end
        #先看是不是温度到下限了
        if Tsaim == Tsmin
            P1 = min((heatLoad-P2-Psout) / COP1,PheatPumpMax)
            Pe = heatLoad-P2-Psout-P1 * COP1
            P21=min(PWaterCompressorMax,heatLoad/COP2value)  # 前半段蓄热供热最大P2功率
            Pe1=heatLoad-P21*COP2value
            P12 = min(heatLoad / COP1,P1Only)   # 后半段热泵供热的功率,P1Only保证了压缩机功率不超过设计功率，所以后面不用比较压缩机功率了
			Pe2 = heatLoad-P12 * COP1
			#P1h=P1*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P1l=P1-P1h
			#P12h=P12*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P12l=P12-P12h
            if (Pe1 <= PelecHeatMax) && (Pe2 <= PelecHeatMax) #&& (P1h<=PWaterCompressorMax) && (P12h<=PWaterCompressorMax) #&& (P1l <= PheatPumpLowMax) && (P12l <= PheatPumpLowMax)
				#println("放热2 Pe1=$Pe1 PelecHeatMax=$PelecHeatMax,Pe2=$Pe2,PelecHeatMax=$PelecHeatMax")
			    return (P1+P2+Pe)*dt,true,P1,P2,Pe
            else
				#println("放热3")
                return 9999.0,false,9999.0,9999.0,9999.0
            end
        end

        #蓄热温度的末态高于用热温度,放热全程热泵补热;或放热到蓄热的最低温度,先放热后热泵
		if (Tsaim-dT_EvaporationStandard>=Tuse)
			P1 = min((heatLoad-P2-Psout) / COP1,P1Only)
            #=OOTest001=#
            #P1=0
			Pe = heatLoad-P2-Psout-P1 * COP1
			#P1h=P1*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P1l=P1-P1h
            if (Pe <= PelecHeatMax) #&& (P1h<=PWaterCompressorMax) #&& (P1l <= PheatPumpLowMax)
			    return (P1+P2+Pe)*dt,true,P1,P2,Pe
            else
                return 9999.0,false,9999.0,9999.0,9999.0
            end
		elseif Tsstart-dT_EvaporationStandard<=Tuse#蓄热温度的初态低于用热温度，放热全程电热补热
			P1 = 0
			Pe = heatLoad-P2-Psout
            if Pe <= PelecHeatMax
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
		# 主要看水蒸气压缩机有没有余量
		P3h = min(PWaterCompressorMax-P1hOnly,cpm_h * (Tsaim - Tsstart) / dt / COP2_design / (1+cp_cw/latentHeat*(Tuse-0.5*(Tsstart+Tsaim))))

		P3 = P3h*COP2_design/COP1
		Pe = heatLoad+cpm_h / dt * (Tsaim - Tsstart) - (P3 + P1Only)*COP1
		
		#P3l = P3 - P3h
		if (Pe <= PelecHeatMax) #&& (P3h<=PWaterCompressorMax) #&& (P3l <= PheatPumpLowMax)
			#println("蓄热1 Pe=$Pe,PelecHeatMax=$PelecHeatMax")
			return (P1Only + P3 + Pe) * dt, true, P1Only, P3, Pe
		else
			#println("蓄热2 Pe=$Pe,PelecHeatMax=$PelecHeatMax")
			return 9999.0,false,9999.0,9999.0,9999.0
		end
	end

	"""计算蓄热温度相对较高时升温的电度情况(考虑蓄热温度超过热泵温度上限)"""
	function powerCalculate_highTs(Tsaim, Tsstart, dt)
		COP3value = COP3(TWaste, (Tsaim + Tsstart) / 2 + dT_EvaporationStandard)
		COPWatervalue = COP2(TCompressorIn, (Tsaim + Tsstart) / 2 + dT_EvaporationStandard)
		function mode1(dt1)# 模式1：只有电热加热给蓄热储热
			Pe = cpm_h / dt1 * (Tsaim - Tsstart)+heatLoad-P1Only*COP1
			if Pe <= PelecHeatMax
				return true, P1Only, 0.0, Pe
			else
				return false, 9999.0, 9999.0, 9999.0
			end
			
		end

		# 模式2：热泵也给蓄热储热（此时热泵COP按照蓄热温度计算）
		function mode2(Tsaim2, Tsstart2, dt2)
			# 如果水蒸气压缩机容量没有富余的话，应该优先用热泵为工厂供热，而不是为蓄热罐储热
			Psin=cpm_h / dt2 * (Tsaim2 - Tsstart2)#蓄热储热功率
			P3hin=min(Psin/COPWatervalue,PWaterCompressorMax)
            #=OOTest100=#
			k=cp_cw/latentHeat * (Tsaim2 - Tsstart2)#闪蒸回收热比例
            #k=0
			u_Pm=(heatLoad/COPWatervalue+PWaterCompressorMax)/(1+k)#线性规划求解的一个分类指标,用来描述用热负荷的大小
			# 先看闪蒸回收热是否够系统需求,如果够的话多余的热直接放掉
			if P3hin*k*COPWatervalue>=heatLoad
				return true, 0.0, P3hin, Psin-P3hin*COPWatervalue
			end
			
			if P3hin>=PWaterCompressorMax && u_Pm>=PWaterCompressorMax
				#蓄热负荷和用热负荷都大,任意进行热泵配置
				P1h=PWaterCompressorMax
				P3h=0.0
			elseif P3hin<PWaterCompressorMax && u_Pm>PWaterCompressorMax
				#蓄热负荷小，用热负荷大,蓄热占满，剩余的用热
				P3h=P3hin
				P1h=min(PWaterCompressorMax,(P3hin+heatLoad/COPWatervalue)/(1+k))-P3h
			elseif P3hin>PWaterCompressorMax && u_Pm<PWaterCompressorMax
				#蓄热负荷大，用热负荷小,蓄热按照供热程度占满
				P1h=0.0
				P3h=heatLoad/COPWatervalue/k
			else
				#蓄热负荷和用热负荷都小
				P3h=min(heatLoad/COPWatervalue/k,P3hin)
				P1h=min(heatLoad/COPWatervalue/k,(P3hin+heatLoad/COPWatervalue)/(1+k))-P3h
			end

			Pe = cpm_h / dt2 * (Tsaim2 - Tsstart2) + heatLoad - (P1h+P3h) * COPWatervalue
			if Pe <= PelecHeatMax
				# flag P1 P3 Pe
				return true, P1h*COPWatervalue/COP1, P3h*COPWatervalue/COP1, Pe
			else
				return false, 9999.0, 9999.0, 9999.0
			end
		end
		#COP3value = COP3(TWaste, (min(Tsaim,TcChangeToElec - dT_EvaporationStandard) + Tsstart) / 2 + dT_EvaporationStandard)
		
		Ptotal = 9999.0
		P1res = 0.0
		P3res = 0.0
		Peres = 0.0
		flag = false
		# 模式1：电热加热
		flag, P1res, P3res, Peres= mode1(dt)
		Ptotal=P1res+P3res+Peres

		# 模式2：热泵加电热,又分两种模式：末态温度不大于电加热温度界限与末态温度大于电加热温度界线
		if (Tsaim <= TcChangeToElec - dT_EvaporationStandard)	# 该时间层的末态温度小于电加热温度界限
			flag_2, P1_2, P3_2, Pe_2= mode2(Tsaim,Tsstart,dt)
			Ptotal_2=P1_2+P3_2+Pe_2
		elseif (Tsaim > TcChangeToElec - dT_EvaporationStandard>Tsstart) # 该时间层的温度跨越电加热温度界限
			tmid=dt*(TcChangeToElec - dT_EvaporationStandard-Tsstart)/(Tsaim - Tsstart)
			flag_21, P1_21, P3_21, Pe_21= mode2(TcChangeToElec - dT_EvaporationStandard,Tsstart,tmid)
			flag_22, P1_22, P3_22, Pe_22= mode2(Tsaim,TcChangeToElec - dT_EvaporationStandard,dt-tmid)
			flag_2=flag_21 && flag_22
			P1_2=P1_21*tmid/dt+P1_22*(dt-tmid)/dt
			P3_2=P3_21*tmid/dt+P3_22*(dt-tmid)/dt
			Pe_2=Pe_21*tmid/dt+Pe_22*(dt-tmid)/dt
			Ptotal_2=P1_2+P3_2+Pe_2
		else# 该时间层初态温度高于电加热温度界线，切换模式1
			flag_2, P1_2, P3_2, Pe_2=mode1(dt)
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

	for (j,Tsaim) in enumerate(TsListEnd)
		for (i,Tsstart) in enumerate(TsListStart)
			# 先计算温度下降的功率,存在一边电加热一边开2号模式的情况
			if Tsstart<Tsmin || Tsaim>Tsmax
				C[i, j] = 9999.0
				P1Matrix[i, j] = 9999.0
				P3Matrix[i, j] = 9999.0
				PeMatrix[i, j] = 9999.0
			elseif Tsaim < Tsstart
				costValue,flag,P1,P2,Pe=powerCalculate_disCharge(Tsaim, Tsstart, dt)
				if flag
					C[i, j]=costValue
					P1Matrix[i, j] = P1
					P2Matrix[i, j] = P2
					PeMatrix[i, j] = Pe
				end
			elseif Tsaim == Tsstart
				temp = heatLoad-P1Only*COP1
				if temp <= PelecHeatMax
					P1Matrix[i, j] = P1Only
					PeMatrix[i, j] = temp
					C[i, j]=(P1Matrix[i, j]+PeMatrix[i, j])*dt
				end
				#=
				P1Matrix[i, j] = P1Only
				PeMatrix[i, j] = heatLoad-P1Only*COP1
				C[i, j]=(P1Matrix[i, j]+PeMatrix[i, j])*dt
				=#
			elseif Tsaim > Tsstart
				if Tsaim <= Tuse - dT_EvaporationStandard
					C1, flag, P1, P3, Pe = powerCalculate_lowTs(Tsaim, Tsstart, dt)
					if flag
						C[i, j] = C1
						P1Matrix[i, j] = P1
						P3Matrix[i, j] = P3
						PeMatrix[i, j] = Pe
					end
					#C[i, j] = flag ? C1 : 99
				elseif Tsstart < Tuse - dT_EvaporationStandard < Tsaim
					dt1 = (Tuse - dT_EvaporationStandard - Tsstart) / (Tsaim - Tsstart) * dt
					dt2 = dt - dt1
					C1, flag1, P11, P31, Pe1 = powerCalculate_lowTs(Tuse - dT_EvaporationStandard, Tsstart, dt1)
					C2, flag2, P12, P32, Pe2 = powerCalculate_highTs(Tsaim, Tuse - dT_EvaporationStandard, dt2)
					flag = flag1 && flag2
					if flag
						C[i, j] = C1 + C2
						P1Matrix[i, j] = (P11 * dt1 + P12 * dt2) / dt
						P3Matrix[i, j] = (P31 * dt1 + P32 * dt2) / dt
						PeMatrix[i, j] = (Pe1 * dt1 + Pe2 * dt2) / dt
					end
				elseif Tsstart >= Tuse - dT_EvaporationStandard
					C1, flag, P1, P3, Pe = powerCalculate_highTs(Tsaim, Tsstart, dt)
					if flag
						C[i, j] = C1
						P1Matrix[i, j] = P1
						P3Matrix[i, j] = P3
						PeMatrix[i, j] = Pe
					end
				end
			end
		end
	end

	return C, P1Matrix, P2Matrix, P3Matrix, PeMatrix
end

"""110"""
function getStateTransitionCost_SingleStep(
	::OO110;
	#C::Matrix, P1Matrix::Matrix, P2Matrix::Matrix, P3Matrix::Matrix, PeMatrix::Matrix;
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
	TCompressorIn::Real,# 压缩机吸气温度

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容

	# 设备运行约束
	PheatPumpMax::Real,# 热泵最大功率
    PelecHeatMax::Real,# 电锅炉最大功率
    PWaterCompressorMax::Real,# 水蒸气压缩机的最大功率
	#PheatPumpLowMax::Real,# 低温热泵的最大功率,一般来说，水蒸气压缩机和低温热泵同时运行时，由于低温热泵按照最劣工况设计，所以水蒸气压缩机会先触碰运行边界，不需要考虑低温热泵的边界问题
    Tsmin::Real,# 蓄热的最小温度
    
	# 求解参数
	TsListStart::Vector,# 状态参数高温蓄热温度起始值列表
	TsListEnd::Vector,# 状态参数高温蓄热温度结束值列表
	dt::Real = 1.0,# 时间步长
	Tsmax::Real = 220.0
)
	COP1 = COPOverlap(TWaste, Tuse)
	COP2 = COPWater
	COP2_design = COP2(TCompressorIn,Tuse)
	COP3 = COPOverlap
	heatLoadPumpMax=PWaterCompressorMax*COP2_design	# 热泵设计的供热功率

	nT = length(TsListStart)# 温度步数

	#C[i,j]表示温度从TsList[i]到TsList[j]时的最低功率；在负载、环境不变的情况下，C[i,j]是不变的
	C = fill(9999.0, nT, nT)# 状态转移矩阵
	P1Matrix = zeros(nT, nT)# 状态转移功率1参数
	P2Matrix = zeros(nT, nT)# 状态转移功率2参数
	P3Matrix = zeros(nT, nT)# 状态转移功率3参数
	PeMatrix = zeros(nT, nT)# 状态转移功率电加热参数

	
	# 只用热泵供热时的功率
	P1Only = min(heatLoad, heatLoadPumpMax)/ COP1# 只用热泵供热时的功率
	P1hOnly = P1Only*COP1/COP2_design	#只用热泵供热时的压缩机功率

	"""计算不同蓄热温度下蓄热温度降低的电度"""
	function powerCalculate_disCharge(Tsaim, Tsstart, dt)
		COP2value=COP2((Tsaim + Tsstart) / 2 - dT_EvaporationStandard, Tuse)
		Psout=cpm_h*(Tsstart-Tsaim)/dt
		P2=Psout / (COP2value-1)
		if (P2+Psout>heatLoad) || (P2 > PWaterCompressorMax) # 降不到这个温度
			return 9999.0,false,9999.0,9999.0,9999.0
		end
        #先看是不是温度到下限了
        if Tsaim == Tsmin
            P1 = min((heatLoad-P2-Psout) / COP1,PheatPumpMax)
            Pe = heatLoad-P2-Psout-P1 * COP1
            P21=min(PWaterCompressorMax,heatLoad/COP2value)  # 前半段蓄热供热最大P2功率
            Pe1=heatLoad-P21*COP2value
            P12 = min(heatLoad / COP1,P1Only)   # 后半段热泵供热的功率,P1Only保证了压缩机功率不超过设计功率，所以后面不用比较压缩机功率了
			Pe2 = heatLoad-P12 * COP1
			#P1h=P1*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P1l=P1-P1h
			#P12h=P12*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P12l=P12-P12h
            if (Pe1 <= PelecHeatMax) && (Pe2 <= PelecHeatMax) #&& (P1h<=PWaterCompressorMax) && (P12h<=PWaterCompressorMax) #&& (P1l <= PheatPumpLowMax) && (P12l <= PheatPumpLowMax)
				#println("放热2 Pe1=$Pe1 PelecHeatMax=$PelecHeatMax,Pe2=$Pe2,PelecHeatMax=$PelecHeatMax")
			    return (P1+P2+Pe)*dt,true,P1,P2,Pe
            else
				#println("放热3")
                return 9999.0,false,9999.0,9999.0,9999.0
            end
        end

        #蓄热温度的末态高于用热温度,放热全程热泵补热;或放热到蓄热的最低温度,先放热后热泵
		if (Tsaim-dT_EvaporationStandard>=Tuse)
			#P1 = min((heatLoad-P2-Psout) / COP1,P1Only)
            #=OOTest001=#
            P1=0
			Pe = heatLoad-P2-Psout-P1 * COP1
			#P1h=P1*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P1l=P1-P1h
            if (Pe <= PelecHeatMax) #&& (P1h<=PWaterCompressorMax) #&& (P1l <= PheatPumpLowMax)
			    return (P1+P2+Pe)*dt,true,P1,P2,Pe
            else
                return 9999.0,false,9999.0,9999.0,9999.0
            end
		elseif Tsstart-dT_EvaporationStandard<=Tuse#蓄热温度的初态低于用热温度，放热全程电热补热
			P1 = 0
			Pe = heatLoad-P2-Psout
            if Pe <= PelecHeatMax
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
		# 主要看水蒸气压缩机有没有余量
		P3h = min(PWaterCompressorMax-P1hOnly,cpm_h * (Tsaim - Tsstart) / dt / COP2_design / (1+cp_cw/latentHeat*(Tuse-0.5*(Tsstart+Tsaim))))

		P3 = P3h*COP2_design/COP1
		Pe = heatLoad+cpm_h / dt * (Tsaim - Tsstart) - (P3 + P1Only)*COP1
		
		#P3l = P3 - P3h
		if (Pe <= PelecHeatMax) #&& (P3h<=PWaterCompressorMax) #&& (P3l <= PheatPumpLowMax)
			#println("蓄热1 Pe=$Pe,PelecHeatMax=$PelecHeatMax")
			return (P1Only + P3 + Pe) * dt, true, P1Only, P3, Pe
		else
			#println("蓄热2 Pe=$Pe,PelecHeatMax=$PelecHeatMax")
			return 9999.0,false,9999.0,9999.0,9999.0
		end
	end

	"""计算蓄热温度相对较高时升温的电度情况(考虑蓄热温度超过热泵温度上限)"""
	function powerCalculate_highTs(Tsaim, Tsstart, dt)
		COP3value = COP3(TWaste, (Tsaim + Tsstart) / 2 + dT_EvaporationStandard)
		COPWatervalue = COP2(TCompressorIn, (Tsaim + Tsstart) / 2 + dT_EvaporationStandard)
		function mode1(dt1)# 模式1：只有电热加热给蓄热储热
			Pe = cpm_h / dt1 * (Tsaim - Tsstart)+heatLoad-P1Only*COP1
			if Pe <= PelecHeatMax
				return true, P1Only, 0.0, Pe
			else
				return false, 9999.0, 9999.0, 9999.0
			end
			
		end

		# 模式2：热泵也给蓄热储热（此时热泵COP按照蓄热温度计算）
		function mode2(Tsaim2, Tsstart2, dt2)
			# 如果水蒸气压缩机容量没有富余的话，应该优先用热泵为工厂供热，而不是为蓄热罐储热
			Psin=cpm_h / dt2 * (Tsaim2 - Tsstart2)#蓄热储热功率
			P3hin=min(Psin/COPWatervalue,PWaterCompressorMax)
            #=OOTest100=#
			k=cp_cw/latentHeat * (Tsaim2 - Tsstart2)#闪蒸回收热比例
            #k=0
			u_Pm=(heatLoad/COPWatervalue+PWaterCompressorMax)/(1+k)#线性规划求解的一个分类指标,用来描述用热负荷的大小
			# 先看闪蒸回收热是否够系统需求,如果够的话多余的热直接放掉
			if P3hin*k*COPWatervalue>=heatLoad
				return true, 0.0, P3hin, Psin-P3hin*COPWatervalue
			end
			
			if P3hin>=PWaterCompressorMax && u_Pm>=PWaterCompressorMax
				#蓄热负荷和用热负荷都大,任意进行热泵配置
				P1h=PWaterCompressorMax
				P3h=0.0
			elseif P3hin<PWaterCompressorMax && u_Pm>PWaterCompressorMax
				#蓄热负荷小，用热负荷大,蓄热占满，剩余的用热
				P3h=P3hin
				P1h=min(PWaterCompressorMax,(P3hin+heatLoad/COPWatervalue)/(1+k))-P3h
			elseif P3hin>PWaterCompressorMax && u_Pm<PWaterCompressorMax
				#蓄热负荷大，用热负荷小,蓄热按照供热程度占满
				P1h=0.0
				P3h=heatLoad/COPWatervalue/k
			else
				#蓄热负荷和用热负荷都小
				P3h=min(heatLoad/COPWatervalue/k,P3hin)
				P1h=min(heatLoad/COPWatervalue/k,(P3hin+heatLoad/COPWatervalue)/(1+k))-P3h
			end

			Pe = cpm_h / dt2 * (Tsaim2 - Tsstart2) + heatLoad - (P1h+P3h) * COPWatervalue
			if Pe <= PelecHeatMax
				# flag P1 P3 Pe
				return true, P1h*COPWatervalue/COP1, P3h*COPWatervalue/COP1, Pe
			else
				return false, 9999.0, 9999.0, 9999.0
			end
		end
		#COP3value = COP3(TWaste, (min(Tsaim,TcChangeToElec - dT_EvaporationStandard) + Tsstart) / 2 + dT_EvaporationStandard)
		
		Ptotal = 9999.0
		P1res = 0.0
		P3res = 0.0
		Peres = 0.0
		flag = false
		# 模式1：电热加热
		flag, P1res, P3res, Peres= mode1(dt)
		Ptotal=P1res+P3res+Peres

		# 模式2：热泵加电热,又分两种模式：末态温度不大于电加热温度界限与末态温度大于电加热温度界线
		if (Tsaim <= TcChangeToElec - dT_EvaporationStandard)	# 该时间层的末态温度小于电加热温度界限
			flag_2, P1_2, P3_2, Pe_2= mode2(Tsaim,Tsstart,dt)
			Ptotal_2=P1_2+P3_2+Pe_2
		elseif (Tsaim > TcChangeToElec - dT_EvaporationStandard>Tsstart) # 该时间层的温度跨越电加热温度界限
			tmid=dt*(TcChangeToElec - dT_EvaporationStandard-Tsstart)/(Tsaim - Tsstart)
			flag_21, P1_21, P3_21, Pe_21= mode2(TcChangeToElec - dT_EvaporationStandard,Tsstart,tmid)
			flag_22, P1_22, P3_22, Pe_22= mode2(Tsaim,TcChangeToElec - dT_EvaporationStandard,dt-tmid)
			flag_2=flag_21 && flag_22
			P1_2=P1_21*tmid/dt+P1_22*(dt-tmid)/dt
			P3_2=P3_21*tmid/dt+P3_22*(dt-tmid)/dt
			Pe_2=Pe_21*tmid/dt+Pe_22*(dt-tmid)/dt
			Ptotal_2=P1_2+P3_2+Pe_2
		else# 该时间层初态温度高于电加热温度界线，切换模式1
			flag_2, P1_2, P3_2, Pe_2=mode1(dt)
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

	for (j,Tsaim) in enumerate(TsListEnd)
		for (i,Tsstart) in enumerate(TsListStart)
			# 先计算温度下降的功率,存在一边电加热一边开2号模式的情况
			if Tsstart<Tsmin || Tsaim>Tsmax
				C[i, j] = 9999.0
				P1Matrix[i, j] = 9999.0
				P3Matrix[i, j] = 9999.0
				PeMatrix[i, j] = 9999.0
			elseif Tsaim < Tsstart
				costValue,flag,P1,P2,Pe=powerCalculate_disCharge(Tsaim, Tsstart, dt)
				if flag
					C[i, j]=costValue
					P1Matrix[i, j] = P1
					P2Matrix[i, j] = P2
					PeMatrix[i, j] = Pe
				end
			elseif Tsaim == Tsstart
				temp = heatLoad-P1Only*COP1
				if temp <= PelecHeatMax
					P1Matrix[i, j] = P1Only
					PeMatrix[i, j] = temp
					C[i, j]=(P1Matrix[i, j]+PeMatrix[i, j])*dt
				end
				#=
				P1Matrix[i, j] = P1Only
				PeMatrix[i, j] = heatLoad-P1Only*COP1
				C[i, j]=(P1Matrix[i, j]+PeMatrix[i, j])*dt
				=#
			elseif Tsaim > Tsstart
				if Tsaim <= Tuse - dT_EvaporationStandard
					C1, flag, P1, P3, Pe = powerCalculate_lowTs(Tsaim, Tsstart, dt)
					if flag
						C[i, j] = C1
						P1Matrix[i, j] = P1
						P3Matrix[i, j] = P3
						PeMatrix[i, j] = Pe
					end
					#C[i, j] = flag ? C1 : 99
				elseif Tsstart < Tuse - dT_EvaporationStandard < Tsaim
					dt1 = (Tuse - dT_EvaporationStandard - Tsstart) / (Tsaim - Tsstart) * dt
					dt2 = dt - dt1
					C1, flag1, P11, P31, Pe1 = powerCalculate_lowTs(Tuse - dT_EvaporationStandard, Tsstart, dt1)
					C2, flag2, P12, P32, Pe2 = powerCalculate_highTs(Tsaim, Tuse - dT_EvaporationStandard, dt2)
					flag = flag1 && flag2
					if flag
						C[i, j] = C1 + C2
						P1Matrix[i, j] = (P11 * dt1 + P12 * dt2) / dt
						P3Matrix[i, j] = (P31 * dt1 + P32 * dt2) / dt
						PeMatrix[i, j] = (Pe1 * dt1 + Pe2 * dt2) / dt
					end
				elseif Tsstart >= Tuse - dT_EvaporationStandard
					C1, flag, P1, P3, Pe = powerCalculate_highTs(Tsaim, Tsstart, dt)
					if flag
						C[i, j] = C1
						P1Matrix[i, j] = P1
						P3Matrix[i, j] = P3
						PeMatrix[i, j] = Pe
					end
				end
			end
		end
	end

	return C, P1Matrix, P2Matrix, P3Matrix, PeMatrix
end

"""111"""
function getStateTransitionCost_SingleStep(
	::OO111;
	#C::Matrix, P1Matrix::Matrix, P2Matrix::Matrix, P3Matrix::Matrix, PeMatrix::Matrix;
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
	TCompressorIn::Real,# 压缩机吸气温度

	# 高温蓄热参数
	cpm_h::Real,# 高温蓄热热容

	# 设备运行约束
	PheatPumpMax::Real,# 热泵最大功率
    PelecHeatMax::Real,# 电锅炉最大功率
    PWaterCompressorMax::Real,# 水蒸气压缩机的最大功率
	#PheatPumpLowMax::Real,# 低温热泵的最大功率,一般来说，水蒸气压缩机和低温热泵同时运行时，由于低温热泵按照最劣工况设计，所以水蒸气压缩机会先触碰运行边界，不需要考虑低温热泵的边界问题
    Tsmin::Real,# 蓄热的最小温度
    
	# 求解参数
	TsListStart::Vector,# 状态参数高温蓄热温度起始值列表
	TsListEnd::Vector,# 状态参数高温蓄热温度结束值列表
	dt::Real = 1.0,# 时间步长
	Tsmax::Real = 220.0
)
	COP1 = COPOverlap(TWaste, Tuse)
	COP2 = COPWater
	COP2_design = COP2(TCompressorIn,Tuse)
	COP3 = COPOverlap
	heatLoadPumpMax=PWaterCompressorMax*COP2_design	# 热泵设计的供热功率

	nT = length(TsListStart)# 温度步数

	#C[i,j]表示温度从TsList[i]到TsList[j]时的最低功率；在负载、环境不变的情况下，C[i,j]是不变的
	C = fill(9999.0, nT, nT)# 状态转移矩阵
	P1Matrix = zeros(nT, nT)# 状态转移功率1参数
	P2Matrix = zeros(nT, nT)# 状态转移功率2参数
	P3Matrix = zeros(nT, nT)# 状态转移功率3参数
	PeMatrix = zeros(nT, nT)# 状态转移功率电加热参数

	
	# 只用热泵供热时的功率
	P1Only = min(heatLoad, heatLoadPumpMax)/ COP1# 只用热泵供热时的功率
	P1hOnly = P1Only*COP1/COP2_design	#只用热泵供热时的压缩机功率

	"""计算不同蓄热温度下蓄热温度降低的电度"""
	function powerCalculate_disCharge(Tsaim, Tsstart, dt)
		COP2value=COP2((Tsaim + Tsstart) / 2 - dT_EvaporationStandard, Tuse)
		Psout=cpm_h*(Tsstart-Tsaim)/dt
		P2=Psout / (COP2value-1)
		if (P2+Psout>heatLoad) || (P2 > PWaterCompressorMax) # 降不到这个温度
			return 9999.0,false,9999.0,9999.0,9999.0
		end
        #先看是不是温度到下限了
        if Tsaim == Tsmin
            P1 = min((heatLoad-P2-Psout) / COP1,PheatPumpMax)
            Pe = heatLoad-P2-Psout-P1 * COP1
            P21=min(PWaterCompressorMax,heatLoad/COP2value)  # 前半段蓄热供热最大P2功率
            Pe1=heatLoad-P21*COP2value
            P12 = min(heatLoad / COP1,P1Only)   # 后半段热泵供热的功率,P1Only保证了压缩机功率不超过设计功率，所以后面不用比较压缩机功率了
			Pe2 = heatLoad-P12 * COP1
			#P1h=P1*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P1l=P1-P1h
			#P12h=P12*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P12l=P12-P12h
            if (Pe1 <= PelecHeatMax) && (Pe2 <= PelecHeatMax) #&& (P1h<=PWaterCompressorMax) && (P12h<=PWaterCompressorMax) #&& (P1l <= PheatPumpLowMax) && (P12l <= PheatPumpLowMax)
				#println("放热2 Pe1=$Pe1 PelecHeatMax=$PelecHeatMax,Pe2=$Pe2,PelecHeatMax=$PelecHeatMax")
			    return (P1+P2+Pe)*dt,true,P1,P2,Pe
            else
				#println("放热3")
                return 9999.0,false,9999.0,9999.0,9999.0
            end
        end

        #蓄热温度的末态高于用热温度,放热全程热泵补热;或放热到蓄热的最低温度,先放热后热泵
		if (Tsaim-dT_EvaporationStandard>=Tuse)
			P1 = min((heatLoad-P2-Psout) / COP1,P1Only)
            #=OOTest001=#
            #P1=0
			Pe = heatLoad-P2-Psout-P1 * COP1
			#P1h=P1*COP1/COP2((Tsaim+Tsstart)/2-dT_EvaporationStandard,Tuse)
			#P1l=P1-P1h
            if (Pe <= PelecHeatMax) #&& (P1h<=PWaterCompressorMax) #&& (P1l <= PheatPumpLowMax)
			    return (P1+P2+Pe)*dt,true,P1,P2,Pe
            else
                return 9999.0,false,9999.0,9999.0,9999.0
            end
		elseif Tsstart-dT_EvaporationStandard<=Tuse#蓄热温度的初态低于用热温度，放热全程电热补热
			P1 = 0
			Pe = heatLoad-P2-Psout
            if Pe <= PelecHeatMax
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
		# 主要看水蒸气压缩机有没有余量
		P3h = min(PWaterCompressorMax-P1hOnly,cpm_h * (Tsaim - Tsstart) / dt / COP2_design / (1+cp_cw/latentHeat*(Tuse-0.5*(Tsstart+Tsaim))))

		P3 = P3h*COP2_design/COP1
		Pe = heatLoad+cpm_h / dt * (Tsaim - Tsstart) - (P3 + P1Only)*COP1
		
		#P3l = P3 - P3h
		if (Pe <= PelecHeatMax) #&& (P3h<=PWaterCompressorMax) #&& (P3l <= PheatPumpLowMax)
			#println("蓄热1 Pe=$Pe,PelecHeatMax=$PelecHeatMax")
			return (P1Only + P3 + Pe) * dt, true, P1Only, P3, Pe
		else
			#println("蓄热2 Pe=$Pe,PelecHeatMax=$PelecHeatMax")
			return 9999.0,false,9999.0,9999.0,9999.0
		end
	end

	"""计算蓄热温度相对较高时升温的电度情况(考虑蓄热温度超过热泵温度上限)"""
	function powerCalculate_highTs(Tsaim, Tsstart, dt)
		COP3value = COP3(TWaste, (Tsaim + Tsstart) / 2 + dT_EvaporationStandard)
		COPWatervalue = COP2(TCompressorIn, (Tsaim + Tsstart) / 2 + dT_EvaporationStandard)
		function mode1(dt1)# 模式1：只有电热加热给蓄热储热
			Pe = cpm_h / dt1 * (Tsaim - Tsstart)+heatLoad-P1Only*COP1
			if Pe <= PelecHeatMax
				return true, P1Only, 0.0, Pe
			else
				return false, 9999.0, 9999.0, 9999.0
			end
			
		end

		# 模式2：热泵也给蓄热储热（此时热泵COP按照蓄热温度计算）
		function mode2(Tsaim2, Tsstart2, dt2)
			# 如果水蒸气压缩机容量没有富余的话，应该优先用热泵为工厂供热，而不是为蓄热罐储热
			Psin=cpm_h / dt2 * (Tsaim2 - Tsstart2)#蓄热储热功率
			P3hin=min(Psin/COPWatervalue,PWaterCompressorMax)
            #=OOTest100=#
			k=cp_cw/latentHeat * (Tsaim2 - Tsstart2)#闪蒸回收热比例
            #k=0
			u_Pm=(heatLoad/COPWatervalue+PWaterCompressorMax)/(1+k)#线性规划求解的一个分类指标,用来描述用热负荷的大小
			# 先看闪蒸回收热是否够系统需求,如果够的话多余的热直接放掉
			if P3hin*k*COPWatervalue>=heatLoad
				return true, 0.0, P3hin, Psin-P3hin*COPWatervalue
			end
			
			if P3hin>=PWaterCompressorMax && u_Pm>=PWaterCompressorMax
				#蓄热负荷和用热负荷都大,任意进行热泵配置
				P1h=PWaterCompressorMax
				P3h=0.0
			elseif P3hin<PWaterCompressorMax && u_Pm>PWaterCompressorMax
				#蓄热负荷小，用热负荷大,蓄热占满，剩余的用热
				P3h=P3hin
				P1h=min(PWaterCompressorMax,(P3hin+heatLoad/COPWatervalue)/(1+k))-P3h
			elseif P3hin>PWaterCompressorMax && u_Pm<PWaterCompressorMax
				#蓄热负荷大，用热负荷小,蓄热按照供热程度占满
				P1h=0.0
				P3h=heatLoad/COPWatervalue/k
			else
				#蓄热负荷和用热负荷都小
				P3h=min(heatLoad/COPWatervalue/k,P3hin)
				P1h=min(heatLoad/COPWatervalue/k,(P3hin+heatLoad/COPWatervalue)/(1+k))-P3h
			end

			Pe = cpm_h / dt2 * (Tsaim2 - Tsstart2) + heatLoad - (P1h+P3h) * COPWatervalue
			if Pe <= PelecHeatMax
				# flag P1 P3 Pe
				return true, P1h*COPWatervalue/COP1, P3h*COPWatervalue/COP1, Pe
			else
				return false, 9999.0, 9999.0, 9999.0
			end
		end
		#COP3value = COP3(TWaste, (min(Tsaim,TcChangeToElec - dT_EvaporationStandard) + Tsstart) / 2 + dT_EvaporationStandard)
		
		Ptotal = 9999.0
		P1res = 0.0
		P3res = 0.0
		Peres = 0.0
		flag = false
		# 模式1：电热加热
		flag, P1res, P3res, Peres= mode1(dt)
		Ptotal=P1res+P3res+Peres

		# 模式2：热泵加电热,又分两种模式：末态温度不大于电加热温度界限与末态温度大于电加热温度界线
		if (Tsaim <= TcChangeToElec - dT_EvaporationStandard)	# 该时间层的末态温度小于电加热温度界限
			flag_2, P1_2, P3_2, Pe_2= mode2(Tsaim,Tsstart,dt)
			Ptotal_2=P1_2+P3_2+Pe_2
		elseif (Tsaim > TcChangeToElec - dT_EvaporationStandard>Tsstart) # 该时间层的温度跨越电加热温度界限
			tmid=dt*(TcChangeToElec - dT_EvaporationStandard-Tsstart)/(Tsaim - Tsstart)
			flag_21, P1_21, P3_21, Pe_21= mode2(TcChangeToElec - dT_EvaporationStandard,Tsstart,tmid)
			flag_22, P1_22, P3_22, Pe_22= mode2(Tsaim,TcChangeToElec - dT_EvaporationStandard,dt-tmid)
			flag_2=flag_21 && flag_22
			P1_2=P1_21*tmid/dt+P1_22*(dt-tmid)/dt
			P3_2=P3_21*tmid/dt+P3_22*(dt-tmid)/dt
			Pe_2=Pe_21*tmid/dt+Pe_22*(dt-tmid)/dt
			Ptotal_2=P1_2+P3_2+Pe_2
		else# 该时间层初态温度高于电加热温度界线，切换模式1
			flag_2, P1_2, P3_2, Pe_2=mode1(dt)
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

	for (j,Tsaim) in enumerate(TsListEnd)
		for (i,Tsstart) in enumerate(TsListStart)
			# 先计算温度下降的功率,存在一边电加热一边开2号模式的情况
			if Tsstart<Tsmin || Tsaim>Tsmax
				C[i, j] = 9999.0
				P1Matrix[i, j] = 9999.0
				P3Matrix[i, j] = 9999.0
				PeMatrix[i, j] = 9999.0
			elseif Tsaim < Tsstart
				costValue,flag,P1,P2,Pe=powerCalculate_disCharge(Tsaim, Tsstart, dt)
				if flag
					C[i, j]=costValue
					P1Matrix[i, j] = P1
					P2Matrix[i, j] = P2
					PeMatrix[i, j] = Pe
				end
			elseif Tsaim == Tsstart
				temp = heatLoad-P1Only*COP1
				if temp <= PelecHeatMax
					P1Matrix[i, j] = P1Only
					PeMatrix[i, j] = temp
					C[i, j]=(P1Matrix[i, j]+PeMatrix[i, j])*dt
				end
				#=
				P1Matrix[i, j] = P1Only
				PeMatrix[i, j] = heatLoad-P1Only*COP1
				C[i, j]=(P1Matrix[i, j]+PeMatrix[i, j])*dt
				=#
			elseif Tsaim > Tsstart
				if Tsaim <= Tuse - dT_EvaporationStandard
					C1, flag, P1, P3, Pe = powerCalculate_lowTs(Tsaim, Tsstart, dt)
					if flag
						C[i, j] = C1
						P1Matrix[i, j] = P1
						P3Matrix[i, j] = P3
						PeMatrix[i, j] = Pe
					end
					#C[i, j] = flag ? C1 : 99
				elseif Tsstart < Tuse - dT_EvaporationStandard < Tsaim
					dt1 = (Tuse - dT_EvaporationStandard - Tsstart) / (Tsaim - Tsstart) * dt
					dt2 = dt - dt1
					C1, flag1, P11, P31, Pe1 = powerCalculate_lowTs(Tuse - dT_EvaporationStandard, Tsstart, dt1)
					C2, flag2, P12, P32, Pe2 = powerCalculate_highTs(Tsaim, Tuse - dT_EvaporationStandard, dt2)
					flag = flag1 && flag2
					if flag
						C[i, j] = C1 + C2
						P1Matrix[i, j] = (P11 * dt1 + P12 * dt2) / dt
						P3Matrix[i, j] = (P31 * dt1 + P32 * dt2) / dt
						PeMatrix[i, j] = (Pe1 * dt1 + Pe2 * dt2) / dt
					end
				elseif Tsstart >= Tuse - dT_EvaporationStandard
					C1, flag, P1, P3, Pe = powerCalculate_highTs(Tsaim, Tsstart, dt)
					if flag
						C[i, j] = C1
						P1Matrix[i, j] = P1
						P3Matrix[i, j] = P3
						PeMatrix[i, j] = Pe
					end
				end
			end
		end
	end

	return C, P1Matrix, P2Matrix, P3Matrix, PeMatrix
end
