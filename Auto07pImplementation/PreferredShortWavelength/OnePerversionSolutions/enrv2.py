print("Loading")
enr = load("enr")
# print("Step-1: decreasing separation")
lvals = list((x * (-0.5) for x in range(100,250)))


enrForward = load("enr")

runModEnergy = run(enr)

# # ###------


#runModEnergy = run(enr,ICP =['ea','mExt','l2Tau'], UZSTOP={'ea': -58.9},NPR=0)

runModEnergy = runModEnergy + run(runModEnergy('UZ'),ICP =['ea','mExt','l2Tau'], UZSTOP={'ea': -57.9})
runModEnergy = runModEnergy + run(runModEnergy('UZ'),ICP =['mExt','ea','l2Tau'], UZSTOP={'mExt': 0.00},DS='-',NPR=1000)
#runModEnergy = run(runModEnergy('UZ'),ICP =['ea','mExt','l2Tau'], ISP=0,UZSTOP={'ea': -90.9}, DS='-',UZR={'ea':lvals})
runModEnergy = run(runModEnergy('UZ'),ICP =['ea','mExt','l2Tau'], UZSTOP={'ea': -55.30},e="enrMZero", DS=1.E-4, UZR={-1:-80.11},NMX=1200,NPR=10,ISP=2)
runModEnergyForward = runModEnergy


#runModEnergy = runModEnergy + run(runModEnergy('UZ'),ICP =['ea','mExt','l2Tau'], ISP=2,UZSTOP={'ea': -52.1930},e="enrMZero", DS=1.E-4, UZR={'ea':lvals},NMX=1200)

runModEnergy1 =  run(runModEnergy('UZ'),ICP =['ea','mExt','l2Tau'], ISP=2,UZSTOP={'ea': -55.2552},e="enrMZero",NPR=1, NMX=60000,DS = 1.E-10, DSMIN = 1.E-15, DSMAX = 2.E-3,NTST=200,NCOL=4,STOP={'LP1'})

#runModEnergy1 +=  run(runModEnergy1('UZ'),ICP =['ea','mExt','l2Tau'], ISP=2,UZSTOP={'ea': -55.0},e="enrMZero",NPR=1, NMX=60000,DS = 1.E-10, DSMIN = 1.E-15, DSMAX = 2.E-5,NTST=200,NCOL=4,STOP={'LP1'})
runModEnergy += runModEnergy1


#--------

# lvals = list((x * (-0.5) for x in range(1,250)))
# runModEnergy = runModEnergy + run(runModEnergy('UZ'),ICP =['ea','mExt','l2Tau'], ISW=-1,ISP=2,UZSTOP={'ea': 0}, DS=1.E-7, UZR={'ea':lvals})
#runModEnergy =  run(runModEnergy('UZ'),ICP =['ea','mExt','l2Tau'], UZSTOP={'ea': -52.9}, UZR={'ea':lvals},DS='-')
#runModEnergy = relabel(runModEnergy)

# runModEnergy =  run(runModEnergy('BP1'),ISW=-1, NPR=0, ISP=2, UZSTOP={'ea': -90.001},UZR={'ea':lvals}, IAD = 1)
# runModEnergy = relabel(runModEnergy)

# lvalsInitial = list((x * (5) for x in range(1,250)))
# enrInitial= load("enr")
# runModEnergyInitial = run(enrInitial,UZR={'ea':lvals}, UZSTOP={'ea': -51.01})

# runModEnergy += runModEnergyInitial




# branchpoints = runModEnergy("BP")
# ii = 0
# for solution in branchpoints:
#     branch1Curve=  run(solution, ISW=-1, NPR=0, ISP=0,NMX=4000, UZSTOP={'ea': -62.}, UZR={-1:[1,2,-50,-0,-85,-75,-90]}, IAD = 1)
#     ii+=1
#     print("bp point:", ii)
#     runModEnergy += branch1Curve
#     if(ii==50):
#          break


# # #runModEnergy = runModEnergy+ run(runModEnergy('BP19'), ISW=-1, NPR=0, ISP=0, UZSTOP={'ea': -60.555})
# runModEnergy = relabel(runModEnergy)
runModEnergyForward = run(runModEnergyForward('UZ'), UZSTOP={'ea': -85.11}, DS='-', UZR={'ea':lvals})
#lvals3 = list(((-0.5*x -120) for x in range(1,150)))
#runModEnergyForward = runModEnergyForward + run(runModEnergyForward('UZ'), UZSTOP={'ea': -60.1}, UZR={'ea':lvals3})


runModEnergy += runModEnergyForward
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
#plot(runModEnergy,use_labels=False,use_symbols=False)
plot(runModEnergy)
#plot(runModEnergy)
wait()

