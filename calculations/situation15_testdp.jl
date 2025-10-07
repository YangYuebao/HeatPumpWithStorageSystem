function forwardSolve(VForward::Vector, C::Matrix, nT::Int)
    VForwardNew = fill(99999.0, nT)
	lastTsIndex = zeros(nT)
	for j in 1:nT
		for i in 1:nT
            temp = VForward[i] + C[i, j]
			if temp < VForwardNew[j]
				VForwardNew[j] = temp
				lastTsIndex[j] = i
			end
		end
	end
	return VForwardNew, lastTsIndex
end

test = [
    2  22  13  20  22
   24  23  20   1  24
    9  10   7   6  13
    5  17  15  22  15
    1  21  20   6  22
]
VForward = zeros(5)

forwardSolve(VForward,test,5)
