"""
记录制冷剂的文件信息
"""
struct Refrigerant
	refrigerant::String
	minTe::Real
	maxTe::Real
	minTc::Real
	maxTc::Real
	step::Real
end

function Base.getproperty(r::Refrigerant, s::Symbol)
	if s == :filePath
		return joinpath(pwd(),"refrigerantPropertys",r.fileName)
	elseif s== :fileName
		fn=fieldnames(Refrigerant)
		temp=map(x->r.x,fn) .|> string
		return join(temp,"_")*".csv"
	else
		return r.s
	end
end

refR134a=Refrigerant("R134a", -10.0, 80.0, 20.0,100.0, 0.1)
refWater=Refrigerant("water", 70.0, 190.0, 70.0, 190.0, 0.1)
refNH3=Refrigerant("Ammonia", -10.0, 125.0, 20.0,125.0, 0.1)
refR1233zdE=Refrigerant("R1233zdE", -10.0, 155.0, 20.0,155.0, 0.1)

"""读取COP文件,如果没找到文件可以自动生成"""
function readCOPFile(r::Refrigerant;generateFile::Bool=true)
	if isfile(r.filePath)
		return CSV.read(r.filePath, DataFrame)
	elseif generateFile
		return generateCOPFile(r)
	else
		throw(error("文件$(r.filePath)不存在,使用generateCOPFile生成文件"))
		return DataFrame()
	end
end

"""生成记录COP的文件路径"""
function generateCOPFilePath()
	fp=joinpath(pwd(), "refrigerantPropertys")
	if !isdir(fp)
		@info "生成目录"*fp
		mkdir(fp)
	end
	@info "目录"*fp*"已存在"
end

"""生成记录COP的文件路径"""
function generateCOPFile(r::Refrigerant;renew=true)
	generateCOPFilePath()
	if !renew && isfile(r.filePath)
		@info "文件"*r.filePath*"已存在"
		return
	end
	CoolProp.PropsSI("H", "T", 100 + 273.15, "Q", 1, "water")

	# 计算并生成COP文件
	@info "文件"*r.filePath*"不存在，开始计算"
	TeList=collect(r.minTe:r.step:r.maxTe)
	TcList=collect(r.minTc:r.step:r.maxTc)
	nTe=length(TeList)
	nTc=length(TcList)
	COPMatrix=fill(9999.0,nTe, nTc)
	Threads.@threads for (i,Tc) in enumerate(TcList)
		for (j,Te) in enumerate(TeList[TeList .< Tc])
			COPMatrix[i,j]=COPTe_Tc(Te, Tc, 1.0, r.refrigerant)
		end
	end
	df=DataFrame(COPMatrix, TcList)
	insertcols!(df,1,:Te=>TeList)
	CSV.write(r.filePath, df)
	@info "文件"*r.filePath*"已生成"
	return df
end

"""
基于物性的COP计算
"""
function COPTe_Tc(Te::Real, Tc::Real, eta_s::Real, refrigerant::String)
	h1 = CoolProp.PropsSI("H", "T", Te + 273.15, "Q", 1, refrigerant)
	s1 = CoolProp.PropsSI("S", "T", Te + 273.15, "Q", 1, refrigerant)
	p2 = CoolProp.PropsSI("P", "T", Tc + 273.15, "Q", 1, refrigerant)
	h2 = CoolProp.PropsSI("H", "S", s1, "P", p2, refrigerant)
	wt = (h2 - h1) / eta_s  # 理论压缩功等于绝热压缩功除以绝热效率
	h3 = CoolProp.PropsSI("H", "T", Tc + 273.15, "Q", 0, refrigerant)
	return (h1 - h3) / wt + 1 # 制热循环效率
end

"""在原表的基础上取出目标温度范围和步长的表格"""
function getCOPTable(
	minTe::Real,
	maxTe::Real,
	minTc::Real,
	maxTc::Real,
	refrigerant::Refrigerant,
	dT::Real,
)
	dT=0.1
	minTe = floor(minTe, digits = 1)
	maxTe = ceil(maxTe, digits = 1)
	minTc = floor(minTc, digits = 1)
	maxTc = ceil(maxTc, digits = 1)
	dT = floor(dT, digits = 1)
	TeList = minTe:dT:maxTe
	TcList = minTc:dT:maxTc

	step = round(Int, dT * 10)

	#Threads.@threads for (j,Tc) in enumerate(TcList)
	dfTemp=readCOPFile(refrigerant)
	iStart = round(Int, (minTe - refrigerant.minTe) * 10 + 1)
	iEnd = round(Int, (maxTe - refrigerant.minTe) * 10 + 1)
	jStart = round(Int, (minTc - refrigerant.minTc) * 10 + 1)
	jEnd = round(Int, (maxTc - refrigerant.minTc) * 10 + 1)
	COPMatrix = Matrix(dfTemp[iStart:step:iEnd, jStart:step:jEnd])
	return COPMatrix,TeList,TcList
end

"""生成COP的线性插值函数"""
function getCOP(
	minTe::Real,# 蒸发温度下限,这里是实际设计中的蒸发冷凝温度界限
	maxTe::Real,# 蒸发温度上限
	minTc::Real,# 冷凝温度下限
	maxTc::Real,# 冷凝温度上限
	refrigerant::Refrigerant,# 工质
	maxCOP::Real,# 最大COP
	eta_s::Real,# 绝热效率
	dT::Real,# 插值步长
)
	TcChangeToElec = maxTc# 冷凝温度转换到电热的阈值
	# 生成基于物性的COP计算函数
	# 然后进行样条插值,并计算梯度和海森矩阵
	COPMatrix,TeList,TcList=getCOPTable(minTe,maxTe,minTc,maxTc,refrigerant,dT)
	COPMatrix = COPMatrix * eta_s .+ (1.0-eta_s)
	COPMatrix[COPMatrix .>= maxCOP] .= maxCOP

	sitpCOP = linear_interpolation((TeList,TcList),COPMatrix)
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

	return COPfunction
end

"""
生成COP的函数、COP的梯度和海森矩阵
"""
function getCOP_g_h(
	minTe::Real,# 蒸发温度下限,这里是实际设计中的蒸发冷凝温度界限
	maxTe::Real,# 蒸发温度上限
	minTc::Real,# 冷凝温度下限
	maxTc::Real,# 冷凝温度上限
	refrigerant::Refrigerant,# 工质
	maxCOP::Real,# 最大COP
	eta_s::Real,# 绝热效率
	dT::Real,# 插值步长
)
	TcChangeToElec = maxTc# 冷凝温度转换到电热的阈值
	# 生成基于物性的COP计算函数
	# 然后进行样条插值,并计算梯度和海森矩阵
	COPMatrix,TeList,TcList=getCOPTable(minTe,maxTe,minTc,maxTc,refrigerant,dT)
	COPMatrix = COPMatrix * eta_s .+ (1.0-eta_s)
	COPMatrix[COPMatrix .>= maxCOP] .= maxCOP
	COPMatrix[COPMatrix <= 0.0] .= maxCOP

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
		# if COP > maxCOP || COP <= 0
		# 	return maxCOP
		# end
		return COP
	end

	function COPfunction_g(g::AbstractVector{T}, x::T...)::Nothing where {T <: Real}
		Te = min(maxTe, max(minTe, x[1]))
		Tc = min(maxTc, max(minTc, x[2]))
		COP = sitpCOP(Te, Tc)
		if (x[2] > TcChangeToElec) #|| (COP > maxCOP) || (COP <= 0)
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
		if (x[2] > TcChangeToElec) #|| (COP > maxCOP) || (COP <= 0)
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
记录两个循环复叠的文件信息
"""
struct OverlapRefrigerant
	refrigerantLow::Refrigerant	# 低温循环工质
	refrigerantHigh::Refrigerant	# 高温循环工质
	minTeHigh::Real
	maxTcLow::Real
	midTDifference::Real	#复叠温差
	step::Real	# 插值步长
	OverlapRefrigerant(refrigerantLow::Refrigerant,refrigerantHigh::Refrigerant,minTeHigh::Real,maxTcLow::Real,midTDifference::Real,step::Real)= maxTcLow-minTeHigh>midTDifference ? new(refrigerantLow,refrigerantHigh,minTeHigh,maxTcLow,midTDifference,step) : throw(error("蒸发器温度与冷凝温度差值必须大于复叠温差"))
end

function Base.getproperty(r::OverlapRefrigerant, s::Symbol)
	if s == :fileName
		temp=string.([
			r.refrigerantLow.refrigerant,
			r.refrigerantHigh.refrigerant,
			r.minTeHigh,
			r.maxTcLow,
			r.midTDifference,
		])
		return join(temp,"_")*".csv"
	elseif s == :fileNameCOP
		return "OverlapCOP"*r.fileName
	elseif s == :filePathCOP
		return joinpath(pwd(),"refrigerantPropertys",r.fileNameCOP)
	elseif s == :fileNameMidT
		return "OverlapMidT"*r.fileName
	elseif s == :filePathMidT
		return joinpath(pwd(),"refrigerantPropertys",r.fileNameMidT)
	else
		return r.s
	end
end

function readCOPFile(or::OverlapRefrigerant;generateFile::Bool=true)
	if isfile(or.filePathCOP)
		return CSV.read(or.filePathCOP, DataFrame)
	elseif generateFile
		return generateCOPFile(or)
	else
		throw(error("文件$(or.filePathCOP)不存在,使用generateCOPFile生成文件"))
		return DataFrame()
	end
end

"""生成COP文件"""
function generateCOPFile(or::OverlapRefrigerant)
	
end

"""计算给定中间温度范围内的最优COP与中间温度"""
function generateCOP(
	or::OverlapRefrigerant,
	TeList::Vector,
	TcList::Vector;
	maxCOP::Real=21.0,
	eta_s::Real=0.7,
	dT::Real=0.1,
	abserr::Real=1e-6,
	maxIter::Int=100,
	saveOverlapCOP::Bool=true
)
	if maximum(TeList) > or.minTeHigh
		throw(error("蒸发温度范围错误"))
	end
	if minimum(TcList) < or.maxTcLow
		throw(error("冷凝温度范围错误"))
	end
	# 读取两个循环的COP函数
	COPl,COPl_g,COPl_h=getCOP_g_h(
		or.refrigerantLow.minTe,
		or.minTeHigh,
		or.minTeHigh,
		or.maxTcLow,
		or.refrigerantLow,
		maxCOP,
		eta_s,# 绝热效率
		dT,# 插值步长
	)
	COPh,COPh_g,COPh_h=getCOP_g_h(
		or.minTeHigh,
		or.maxTcLow,
		or.maxTcLow,
		or.refrigerantHigh.maxTc,
		or.refrigerantHigh,
		maxCOP,
		eta_s,# 绝热效率
		dT,# 插值步长
	)
	dT_l = or.midTDifference/2
	function getCOPValues(Te,Tc,Tm)
		cl=COPl(Te,Tm+dT_l)
        ch=COPh(Tm-dT_l,Tc)
        c=cl*ch/(cl+ch-1)
		return c
	end
	function getStep(Te,Tc,Tm)
        # 计算初始值
        cl=COPl(Te,Tm+dT_l)
        ch=COPh(Tm-dT_l,Tc)
        c=cl*ch/(cl+ch-1)
        # 计算导数
        dcl=COPl_g(Te,Tm+dT_l)[2]
        dch=COPh_g(Tm-dT_l,Tc)[1]
        # 计算二阶导数
        ddcl=COPl_h(Te,Tm+dT_l)[2,2]
        ddch=COPh_h(Tm-dT_l,Tc)[1,1]
        ddc=((ddcl*(ch-c)+ddch*(cl-c))+2*(dcl*dch-dc*(dcl+dch)))/(cl+ch-1)
        
        # 返回步长
        return -dc/ddc,c
    end
	nTe=length(TeList)
	nTc=length(TcList)
	COPMatrix=zeros(nTe,nTc)
	MidTMatrix=zeros(nTe,nTc)
	@info "正在计算"*r.fileName*"..."
	Threads.@threads for (j,Tc) in enumerate(TcList)
		for (i,Te) in enumerate(TeList)
			Tm=0.5*(or.minTeHigh+or.maxTcLow)
			dMidT=1.0
			count=0
			while abs(dMidT)>abserr && count < maxIter
				dMidT=getStep(Te,Tc,Tm)
				if dMidT >0 && Tm > or.maxTcLow - dT_l - abserr
					Tm=or.maxTcLow
					break
				end
				if dMidT < 0 && Tm < or.minTeHigh + dT_l + abserr
					Tm=or.minTeHigh
					break
				end
				Tm=max(or.minTeHigh + dT_l,min(Tm+dMidT,or.maxTcLow - dT_l))
				count+=1
			end
			MidTMatrix[i,j]=Tm
			COPMatrix[i,j]=getCOPValues(Te,Tc,Tm)
		end
	end

	if saveOverlapCOP
		dfCOP=DataFrame(COPMatrix, TcList)
		insertcols!(dfCOP,1,:Te=>TeList)
		CSV.write(r.filePathCOP, dfCOP)

		dfMidT=DataFrame(MidTMatrix, TcList)
		insertcols!(dfMidT,1,:Te=>TeList)
		CSV.write(r.filePathMidT, dfMidT)
		@info "文件"*r.fileName*"已生成"
	end
	
	return dfCOP,dfMidT
end

"""
计算双循环复叠系统的最优COP
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

"""
计算双循环复叠系统固定中间温度的COP
"""
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