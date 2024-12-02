
using Pkg
Pkg.add("Plots") 
Pkg.instantiate()

#using Revise
using Plots
using HeatPumpWithStorageSystem

hourlyTariff=ones(24)*0.77
hourlyTariff[4:6].*=1.6
hourlyTariff[11:15].*=1.6
hourlyTariff[16:23].*=0.4

COP,heatConsumptionPowerList,heatStorageCapacityConstraint,heatpumpPowerConstraint,hourlyTariffList=generateSystemCoff(HeatPumpStoragePhaseChange();refrigerant="water",hourlyTariff=hourlyTariff)

cost, Plist,heatStorageList=generateAndSolve(HeatPumpStoragePhaseChange(),MinimizeCost();
    simulationDays=7,
    COP=COP,
    heatConsumptionPowerList=heatConsumptionPowerList,
    heatStorageCapacityConstraint=heatStorageCapacityConstraint,
    heatpumpPowerConstraint=heatpumpPowerConstraint,
    hourlyTariffList=hourlyTariffList
)


plt1=plot(1:7*24,Plist,title="Power of Heatpump",xlabel="Hour",ylabel="kW",legend=:none)
plt2=plot(1:7*24,heatStorageList,title="Heat Storage",xlabel="Hour",ylabel="kWh",legend=:none)


plt3=plot(8:32,Plist[3*24+1:4*24+1],title="Power of Heatpump",xlabel="Hour",ylabel="kW",legend=:none)
plt4=plot(8:32,heatStorageList[3*24+1:4*24+1],title="Heat Storage",xlabel="Hour",ylabel="kWh",legend=:none)
plt5=plot(1:24,hourlyTariff,title="Hourly Tariff",xlabel="Hour",ylabel="kW",legend=:none)


begin
Te=70
Tc=180
refrigerant="water"
h1 = CoolProp.PropsSI("H","T",Te+273.15,"Q",1,refrigerant)
s1 = CoolProp.PropsSI("S","T",Te+273.15,"Q",1,refrigerant)
p2 = CoolProp.PropsSI("P","T",Tc+273.15,"Q",1,refrigerant)
h2 = CoolProp.PropsSI("H","S",s1,"P",p2,refrigerant)
wt = (h2-h1)/0.7
h3 = CoolProp.PropsSI("H","T",Tc+273.15,"Q",0,refrigerant)
COP = (h1-h3)/wt + 1
end

a=2
function haha()
    a=3
    function testf1(x)
        return x^2+1+a
    end
    return testf1
end
