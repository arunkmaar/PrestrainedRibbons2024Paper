print("Loading")
enr = load("enr")
# print("Step-1: decreasing separation")
lvals = list((x * (-0.25) for x in range(150,350)))


enrForward = load("enr")
enrInitial = enrForward

runModEnergy = run(enr)

# # ###------
#runModEnergyForward = runModEnergy;

#runModEnergy = run(enr,ICP =['ea','mExt','l2Tau'], UZSTOP={'ea': -58.9},NPR=0)

runModEnergy = runModEnergy + run(runModEnergy('UZ'),ISP=0,ICP =['ea','mExt','l2Tau'], UZSTOP={'ea': -45.9})
runModEnergy = runModEnergy + run(runModEnergy('UZ'),ISP=0,ICP =['mExt','ea','l2Tau'], UZSTOP={'mExt': 0.00},DS='-',NPR=1000)
print("step 0")
runModEnergy = run(runModEnergy('UZ'),ICP =['ea','mExt','l2Tau'], ISP=0,UZSTOP={'ea': -50.111}, DS='-',UZR={'ea':lvals},e="enrMZero")
runModEnergyForward = runModEnergy
lvals = list((x * (-0.25) for x in range(80,200)))
print("step 1")
runModEnergy += run(runModEnergy('UZ'),ICP =['ea','mExt','l2Tau'], ISP=0,UZSTOP={'ea': -32.8},e="enrMZero", DS=1.E-4, UZR={'ea':lvals},NMX=12000,NPR=100,NTST=400,NCOL=6)

print("step 2")
runModEnergy2 =  run(runModEnergy('UZ'),ICP =['ea','mExt','l2Tau'], ISP=0,UZSTOP={'ea': {-31.18}},e="enrMZero", DS=1.E-5,DSMAX=1.E-1,DSMIN=1.E-14, UZR={'ea':lvals},NMX=120000,NPR=200,STOP={'LP1','BP50'})

print("near bifurcation")
runModEnergy2 +=  run(runModEnergy2('UZ'),ICP =['ea','mExt','l2Tau'], ISP=2,UZSTOP={'ea': {-29.0}},e="enrMZero", DS=1.E-9,DSMAX=1.E-1,DSMIN=1.E-14, UZR={'ea':lvals},NMX=120000,NPR=1,STOP={'LP1'},IAD=1)
runModEnergy += runModEnergy2


#runModEnergy = runModEnergy + run(runModEnergy('UZ'),ICP =['ea','mExt','l2Tau'], ISP=2,UZSTOP={'ea': -52.1930},e="enrMZero", DS=1.E-4, UZR={'ea':lvals},NMX=1200)

#runModEnergy1 =  run(runModEnergy('UZ'),ICP =['ea','mExt','l2Tau'], ISP=2,UZSTOP={'ea': -50.},e="enrMZero",NPR=100, NMX=60000,DS = 1.E-10, DSMIN = 1.E-15, DSMAX = 2.E-7,NTST=200,NCOL=4,STOP={'LP1','BP1'},EPSL = 1e-10, EPSU = 1e-10, EPSS = 1e-7)

#runModEnergy += runModEnergy1


# --------

# lvals = list((x * (-0.5) for x in range(1,250)))
# runModEnergy = runModEnergy + run(runModEnergy('UZ'),ICP =['ea','mExt','l2Tau'], ISW=-1,ISP=2,UZSTOP={'ea': 0}, DS=1.E-7, UZR={'ea':lvals})
# runModEnergy =  run(runModEnergy('UZ'),ICP =['ea','mExt','l2Tau'], UZSTOP={'ea': -52.9}, UZR={'ea':lvals},DS='-')
# runModEnergy = relabel(runModEnergy)

# runModEnergy =  run(runModEnergy('BP1'),ISW=-1, NPR=0, ISP=2, UZSTOP={'ea': -90.001},UZR={'ea':lvals}, IAD = 1)
# runModEnergy = relabel(runModEnergy)

# lvalsInitial = list((x * (5) for x in range(1,250)))
# enrInitial= load("enr")
# runModEnergyInitial = run(enrInitial,UZR={'ea':lvals}, UZSTOP={'ea': -51.01})

# runModEnergy += runModEnergyInitial




# branchpoints = runModEnergy2("BP")
# ii = 0
# for solution in branchpoints:
#     branch1Curve=  run(solution, ISW=-1, NPR=300, ISP=0,NMX=4000, UZSTOP={'ea': -0.}, UZR={-1:[-29.2,-30.]}, IAD = 1,DS=1.E-8,DSMIN=1.E-13,DSMAX=1.E-3)
#     ii+=1
#     print("bp point:", ii)
#     runModEnergy2 += branch1Curve
#     if(ii==500):
#         break

#runModEnergy += run(runModEnergy2('BP19'), ISW=-1, NPR=300, ISP=0,NMX=4000, UZSTOP={'ea': -0.}, UZR={-1:[-29.2,-30.]}, IAD = 1,DS=1.E-8,DSMIN=-1.E-13,DSMAX=1.E-4)

# # #runModEnergy = runModEnergy+ run(runModEnergy('BP19'), ISW=-1, NPR=0, ISP=0, UZSTOP={'ea': -60.555})
# runModEnergy = relabel(runModEnergy)


###-------------


lvals = list((x * (-1.) for x in range(50,160)))
runModEnergyForward =  run(runModEnergy('UZ'),ICP =['ea','mExt','l2Tau'], ISP=0,UZSTOP={'ea': -153.08},e="enrMZero", DS=-1.E-2, UZR={'ea':lvals},NMX=40000,NPR=0,STOP={})
runModEnergy += runModEnergyForward

lvals3 = list(((-x) for x in range(1,100)))
runModInitial =  run(enrInitial,ICP =['ea','mExt','l2Tau'], ISP=2,UZSTOP={'ea': -31.148},e="enrMZero", DS=-1.E-3, UZR={'ea':lvals3},NMX=12000,NPR=0,DSMAX=1E-1)

runModEnergy += runModInitial
# ####---------------



# ######-------------------
#lvals2 = list((x * (-5) for x in range(1,100)))
#runModEnergyForward = run(enrForward,ICP =['ea','mExt','l2Tau'], ISP=2,UZSTOP={'ea': -0.219},UZR={'ea':lvals2})

#runModEnergy += runModEnergyForward

# # #runModEnergy = runModEnergy + run(runModEnergy('UZ'), NPR=0, ICP = ['thetaExt','l2Tau','ea'], UZSTOP={'thetaExt': 0.0}, e="enr3", DS="-",EPSL = 1e-10, EPSU = 1e-10, EPSS = 1e-10)
# # runModEnergy = relabel(runModEnergy)
# # # runModEnergy = run(runModEnergy('UZ'),ISP=2, NPR=500, ICP = ['eb','l2Kappa','thetaExt'], UZSTOP={'eb': -31.0}, e="enr4",EPSL = 1e-10, EPSU = 1e-10, EPSS = 1e-10, UZR={'eb':lvals})
# # # runModEnergy = relabel(runModEnergy)

# # # #runModEnergy = run(runModEnergy('UZ'),ISP=0, ICP = ['eb','l2Kappa','thetaExt'], UZSTOP={'eb': -120.2}, e="enr4", DS="-", UZR={'eb':lvals})


# # # runModEnergy = runModEnergy+ run(runModEnergy('UZ'),ISP=2, ICP = ['eb','l2Kappa','thetaExt'],NPR=0, UZSTOP={'eb': -0.}, e="enr4", DS=1.E-6, DSMAX=1.,EPSL = 1e-6, EPSU = 1e-6, EPSS = 1e-6, UZR={})
# # # runModEnergy = relabel(runModEnergy)

# # # #run(runModEnergy('BP1'),ISW=-1, ISP=2, ICP = ['eb','l2Kappa','thetaExt'],NPR=0, UZSTOP={'eb': -0}, e="enr4", DS=1.E-6, DSMAX=1.)
# ###----

runModEnergy = relabel(runModEnergy)
save(runModEnergy,'runModEnergy')
plot(runModEnergy,use_labels=False,use_symbols=False)
# # #plot(runModEnergy)
#plot(runModEnergy)
wait()

