"""
基于物性的COP计算
"""
function COPTe_Tc(Te, Tc, eta_s, refrigerant)
	h1 = CoolProp.PropsSI("H", "T", Te + 273.15, "Q", 1, refrigerant)
	s1 = CoolProp.PropsSI("S", "T", Te + 273.15, "Q", 1, refrigerant)
	p2 = CoolProp.PropsSI("P", "T", Tc + 273.15, "Q", 1, refrigerant)
	h2 = CoolProp.PropsSI("H", "S", s1, "P", p2, refrigerant)
	wt = (h2 - h1) / eta_s  # 理论压缩功等于绝热压缩功除以绝热效率
	h3 = CoolProp.PropsSI("H", "T", Tc + 273.15, "Q", 0, refrigerant)
	return (h1 - h3) / wt + 1 # 制热循环效率
end


"""
生成COP的函数、COP的梯度和海森矩阵
"""
function getCOP_g_h(
	minTe::Real,# 蒸发温度下限
	maxTe::Real,# 蒸发温度上限
	minTc::Real,# 冷凝温度下限
	maxTc::Real,# 冷凝温度上限
	refrigerant::String,# 工质
	maxCOP::Real,# 最大COP
	eta_s::Real,# 绝热效率
	dT::Real,# 插值步长
)
	TcChangeToElec = maxTc# 冷凝温度转换到电热的阈值
	# 生成基于物性的COP计算函数
	# 然后进行样条插值,并计算梯度和海森矩阵
	minTe = floor(minTe, digits = 1)
	maxTe = ceil(maxTe, digits = 1)
	minTc = floor(minTc, digits = 1)
	maxTc = ceil(maxTc, digits = 1)
	dT = floor(dT, digits = 1)
	TeList = minTe:dT:maxTe
	TcList = minTc:dT:maxTc

	step = round(Int, dT * 10)

	#Threads.@threads for (j,Tc) in enumerate(TcList)
	if refrigerant in ["R134a"]
        fp=joinpath(pwd(), "src", "refrigerantPropertys", "R134a_10_80_20_100_0.1.csv")
		if isfile(fp)
            dfTemp = CSV.read(fp, DataFrame)
        else
            throw(error("文件$(fp)不存在"))
        end
		iStart = round(Int, (minTe - 10) * 10 + 1)
		iEnd = round(Int, (maxTe - 10) * 10 + 1)
		jStart = round(Int, (minTc - 20) * 10 + 1)
		jEnd = round(Int, (maxTc - 20) * 10 + 1)
	elseif refrigerant in ["water", "Water"]
        fp=joinpath(pwd(), "src", "refrigerantPropertys", "water_70_190_70_190_0.1.csv")
		if isfile(fp)
            dfTemp = CSV.read(fp, DataFrame)
        else
            throw(error("文件$(fp)不存在"))
        end
		iStart = round(Int, (minTe - 70) * 10 + 1)
		iEnd = round(Int, (maxTe - 70) * 10 + 1)
		jStart = round(Int, (minTc - 70) * 10 + 1)
		jEnd = round(Int, (maxTc - 70) * 10 + 1)
		#@info any(isnan.(dfTemp))
	end
	# println(minTe," ", maxTe," ", minTc," ", maxTc)
	# println(size(dfTemp))
	# println(iStart," ", iEnd," ", jStart," ", jEnd," ", step)
	COPMatrix = Matrix(dfTemp[iStart:step:iEnd, jStart:step:jEnd]) * eta_s .+ 1.0

	m, n = size(COPMatrix)
	count = 0
	for j ∈ 1:n
		for i ∈ 1:m
			if COPMatrix[i, j] >= maxCOP || COPMatrix[i, j] <= 1
				COPMatrix[i, j] = maxCOP
			end
			count += 1
		end
	end

	itpCOP = interpolate(COPMatrix, BSpline(Cubic(Line(OnGrid()))))
	sitpCOP = scale(itpCOP, TeList, TcList)
	#println("minTe: $minTe, maxTe: $maxTe, minTc: $minTc, maxTc: $maxTc ")
	function COPfunction(x::T...)::T where {T <: Real}
		Te = min(maxTe, max(minTe, x[1]))
		Tc = min(maxTc, max(minTc, x[2]))
		COP = sitpCOP(Te, Tc)
		if x[2] > TcChangeToElec
			return 1.0
		end
		if COP > maxCOP || COP <= 0
			return maxCOP
		end
		return COP
	end

	function COPfunction_g(g::AbstractVector{T}, x::T...)::Nothing where {T <: Real}
		Te = min(maxTe, max(minTe, x[1]))
		Tc = min(maxTc, max(minTc, x[2]))
		COP = sitpCOP(Te, Tc)
		if (x[2] > TcChangeToElec) || (COP > maxCOP) || (COP <= 0)
			g[1], g[2] = zeros(2)
			return nothing
		end
		g[1], g[2] = Interpolations.gradient(sitpCOP, Te, Tc) |> Vector
		return nothing
	end

	function COPfunction_h(H::AbstractMatrix{T}, x::T...)::Nothing where {T <: Real}
		Te = min(maxTe, max(minTe, x[1]))
		Tc = min(maxTc, max(minTc, x[2]))
		COP = sitpCOP(Te, Tc)
		if (x[2] > TcChangeToElec) || (COP > maxCOP) || (COP <= 0)
			H[1, 1], H[2, 1], _, H[2, 2] = zeros(4)
			return nothing
		end
		H[1, 1], H[2, 1], _, H[2, 2] = Interpolations.hessian(sitpCOP, Te, Tc) |> Matrix

		H = LowerTriangular(H)

		return nothing
	end

	return COPfunction, COPfunction_g, COPfunction_h
end


"""
生成系统的6个COP函数
"""
function getCOPFunction(
	Tair::Vector,
	dTair::Real,
	dT_l::Real,
	TWaste::Real,
	maxTcLow::Real,
	dTlc_he::Real,
	maxTeh::Real,
	maxTcHigh::Real,
	refrigerantHigh::String,
	refrigerantLow::String,
	maxCOP::Real,
	eta_s::Real,
	COPInterpolateGap::Real,
)
	TeAirSource = Tair .- dTair                # 空气源蒸发器温度
	minTel = minimum(TeAirSource)
	maxTel = max(TWaste - dT_l, TeAirSource...)
	minTcl = minimum(Tair)
	minTeh = maxTcLow - dTlc_he
	minTch = minTeh + 0.1

	COPh, COPh_g, COPh_h = getCOP_g_h(minTeh, maxTeh, minTch, maxTcHigh, refrigerantHigh, maxCOP, eta_s, COPInterpolateGap)

	COPl, COPl_g, COPl_h = getCOP_g_h(minTel, maxTel, minTcl, maxTcLow, refrigerantLow, maxCOP, eta_s, COPInterpolateGap)
	return COPh, COPh_g, COPh_h,COPl, COPl_g, COPl_h
end

"""
计算复叠系统的COP
"""
function getOverlapCOP_calculation(
	Tair::Vector,
	dTair::Real,
	dT_l::Real,
	TWaste::Real,
	maxTcLow::Real,
	dTlc_he::Real,
	maxTeh::Real,
	maxTcHigh::Real,
	refrigerantHigh::String,
	refrigerantLow::String,
	maxCOP::Real,
	eta_s::Real,
	COPInterpolateGap::Real,
)
	COPh, COPh_g, COPh_h, COPl, COPl_g, COPl_h = getCOPFunction(
        Tair, dTair, dT_l, TWaste, maxTcLow, dTlc_he, maxTeh, maxTcHigh,
		refrigerantHigh,
		refrigerantLow,
		maxCOP,
		eta_s,
		COPInterpolateGap,
	)

    """
    计算总COP的变化步长
    """
    function getdC(Te,Tc,Tm,dT_l)
        # 计算初始值
        cl=COPl(Te,Tm+dT_l/2)
        ch=COPh(Tm-dT_l/2,Tc)
        c=cl*ch/(cl+ch-1)
        # 计算导数
        dcl=COPl_g(Te,Tm+dT_l/2)[2]
        dch=COPh_g(Tm-dT_l/2,Tc)[1]
        # 计算二阶导数
        ddcl=COPl_h(Te,Tm+dT_l/2)[2,2]
        ddch=COPh_h(Tm-dT_l/2,Tc)[1,1]
        ddc=((ddcl*(ch-c)+ddch*(cl-c))+2*(dcl*dch-dc*(dcl+dch)))/(cl+ch-1)
        
        # 返回步长
        return -dc/ddc
    end
    
    maxTm = maxTcLow - dT_l/2
    minTeh = maxTcLow - dTlc_he
    minTm = minTeh + dT_l/2

    
end

function getOverlapCOP_fixMidTemperature(
	Tair::Vector,
	dTair::Real,
	dT_l::Real,
	TWaste::Real,
	maxTcLow::Real,
	dTlc_he::Real,
	maxTeh::Real,
	maxTcHigh::Real,
	refrigerantHigh::String,
	refrigerantLow::String,
	maxCOP::Real,
	eta_s::Real,
	COPInterpolateGap::Real,
)
	COPh, COPh_g, COPh_h, COPl, COPl_g, COPl_h = getCOPFunction(
        Tair, dTair, dT_l, TWaste, maxTcLow, dTlc_he, maxTeh, maxTcHigh,
		refrigerantHigh,
		refrigerantLow,
		maxCOP,
		eta_s,
		COPInterpolateGap,
	)

    """
    计算总COP的变化步长
    """
    function getdC(Te,Tc,Tm,dT_l)
        # 计算初始值
        cl=COPl(Te,Tm+dT_l/2)
        ch=COPh(Tm-dT_l/2,Tc)
        c=cl*ch/(cl+ch-1)
        # 计算导数
        dcl=COPl_g(Te,Tm+dT_l/2)[2]
        dch=COPh_g(Tm-dT_l/2,Tc)[1]
        # 计算二阶导数
        ddcl=COPl_h(Te,Tm+dT_l/2)[2,2]
        ddch=COPh_h(Tm-dT_l/2,Tc)[1,1]
        ddc=((ddcl*(ch-c)+ddch*(cl-c))+2*(dcl*dch-dc*(dcl+dch)))/(cl+ch-1)
        
        # 返回步长
        return -dc/ddc
    end
    
    maxTm = maxTcLow - dT_l/2
    minTeh = maxTcLow - dTlc_he
    minTm = minTeh + dT_l/2

    
end