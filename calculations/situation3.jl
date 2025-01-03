using Revise
using Plots, CSV, DataFrames
using HeatPumpWithStorageSystem

#=
可用工质：
R134a
R245fa
R1234ze(E)
R1233zd(E)
Ammonia 氨气
Water
=#
begin
    hourlyTariff = zeros(24)
    hourlyTariff[1:6] .= 0.3340
    hourlyTariff[7:10] .= 0.7393
    hourlyTariff[11:13] .= 1.2360
    hourlyTariff[14:17] .= 0.7393
    hourlyTariff[18:22] .= 1.2360
    hourlyTariff[23:24] .= 0.3340

    Twastein = 40
    Twasteout = 30
    Tuse = 80
    dTuse = 10

    heatStorageOutEfficiency = 0.98
    heatStorageInEfficiency = 1.0

    heatConsumptionPower = 1.0     # 每小时用热功率kW
    workingHours = 16
    workingStartHour = 8

    refrigerantList = ["R134a", "R245fa", "R1234ze(E)", "R1233zd(E)", "Ammonia", "Water"]
    maxHeatStorageHour=18    #蓄热时长最大值
    maxheatStorageInputHour = 4
end

dfSummary = DataFrame("蓄热时长"=>0:maxHeatStorageHour)
for refrigerant in refrigerantList
    dailyCost = zeros(maxHeatStorageHour+1)
    filePath=joinpath(pwd(), "calculations", "situation3","rf$(refrigerant)")
    filePathPower = joinpath(filePath,"power")
    filePathHeatStorage = joinpath(filePath,"heatStorage")
    if !isdir(filePath)
        mkdir(filePath)
    end
    if !isdir(filePathPower)
        mkdir(filePathPower)
    end
    if !isdir(filePathHeatStorage)
        mkdir(filePathHeatStorage)
    end
    for heatStorageHour in 0:maxHeatStorageHour
        heatStorageCapacity = heatConsumptionPower * heatStorageHour / heatStorageOutEfficiency      # 蓄热量kWh(相变蓄热)
        COP, heatConsumptionPowerList, heatStorageCapacityConstraint, heatpumpPowerConstraint, hourlyTariffList, _,_,heatStorageVelocity = generateSystemCoff(
            HeatPumpStoragePhaseChange();
            refrigerant=refrigerant,
            Twastein=Twastein,
            Twasteout=Twasteout,
            Tuse=Tuse,
            heatStorageCapacity=heatStorageCapacity,
            heatConsumptionPower=heatConsumptionPower,
            heatStorageOutEfficiency=heatStorageOutEfficiency,
            heatStorageInEfficiency=heatStorageInEfficiency,
            hourlyTariff=hourlyTariff,
            workingHours=workingHours,
            workingStartHour=workingStartHour,
            maxheatStorageInputHour=maxheatStorageInputHour,
            dTuse=dTuse
        )

        costList, PList, heatStorageList, heatStorageOutput, heatpumpOutput, COPList = generateAndSolve(HeatPumpStoragePhaseChange(), MinimizeCost();
            simulationDays=7,
            COP=COP,
            heatConsumptionPowerList=heatConsumptionPowerList,
            heatStorageCapacityConstraint=heatStorageCapacityConstraint,
            heatpumpPowerConstraint=heatpumpPowerConstraint,
            hourlyTariffList=hourlyTariffList,
            heatStorageOutEfficiency=heatStorageOutEfficiency,
            heatStorageInEfficiency=heatStorageInEfficiency,
            heatStorageVelocity=heatStorageVelocity
        )

        dailyCost[heatStorageHour+1] = sum(costList)

        plt3 = plot(1:25, vcat(PList, PList[1]), title="$(heatStorageHour)h Power rf=$(refrigerant)", xlabel="Hour", ylabel="kW", legend=:none)
        plt4 = plot(1:25, vcat(heatStorageList, heatStorageList[1]), title="$(heatStorageHour)h Heat Storage  rf=$(refrigerant)", xlabel="Hour", ylabel="kWh", legend=:none)
        #plt5 = plot(1:25, vcat(hourlyTariff, hourlyTariff[1]), title="Hourly Tariff", xlabel="Hour", ylabel="kW", legend=:none)

        df = DataFrame(
            "当前小时"=>0:23,
            "当前小时电价"=>hourlyTariff,
            "当小时电费" => costList,
            "当小时电度" => PList,
            "时刻蓄热容量" => heatStorageList,
            #"蓄热供热量" => heatStorageOutput,
            "热泵制热量" => heatpumpOutput,
            "热泵COP" => COPList,
            "热量需求" => heatConsumptionPowerList
        )
     
        CSV.write(joinpath(filePath, "每小每千瓦时热需求_工况3_工质$(refrigerant)_蓄热$(heatStorageHour)h.csv"), df)
        savefig(plt3, joinpath(filePathPower, "工况3_功率变化_工质$(refrigerant)_蓄热$(heatStorageHour)h.png"))
        savefig(plt4, joinpath(filePathHeatStorage, "工况3_蓄热变化_工质$(refrigerant)_蓄热$(heatStorageHour)h.png"))
    end
    dfSummary[!,refrigerant]=dailyCost
end

plt5 = plot(0:maxHeatStorageHour,Matrix(dfSummary[:,2:end]),label=reshape(refrigerantList,:,length(refrigerantList)),xlabel="Heat Storage Hour",ylabel = "Daily Cost")
CSV.write(joinpath(pwd(), "calculations", "situation3", "summary.csv"), dfSummary)
savefig(plt5, joinpath(pwd(), "calculations", "situation3", "summary.png"))

