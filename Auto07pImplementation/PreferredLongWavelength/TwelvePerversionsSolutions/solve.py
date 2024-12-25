# Start of the guess step
print("Starting guess step")

# Load the initial guess step
guessStep = load("guessStep")

# Save the initial state with external moment mExt = 0 for later use
planarRibState = guessStep

print("Run guess step with increasing sinusoidal external moment (amplitude specified in c.guessStep)")
runGuessStep = run(guessStep)

print("At fixed external moment, decreasing separation between ribbon ends to epsilon = -45.9")
runGuessStep = run(
    runGuessStep('UZ'),
    ICP=['epsilon', 'mExt', 'l2Tau'],
    UZSTOP={'epsilon': -45.9},
    ISP=0
)

print("Set external moment to zero")
runGuessStep = run(
    runGuessStep('UZ'),
    ICP=['mExt', 'epsilon', 'l2Tau'],
    UZSTOP={'mExt': 0.00},
    DS='-',
    NPR=1000
)

print("Trace the one perversion solution branch by changing separation between ribbon ends")

# Define epsilon values to print ribbon states
epsilonValues = [x * (-0.25) for x in range(170, 210)]

print("Use changeSeparation.f90 that sets mExt = 0, and change separation to epsilon=-50.11 afterwards")
changeSeparationToEMinus50 = run(
    runGuessStep('UZ'),
    ICP=['epsilon', 'mExt', 'l2Tau'],
    ISP=0,
    UZSTOP={'epsilon': -50.111},
    DS='-',
    UZR={'epsilon': epsilonValues},
    e="changeSeparation"
)

epsilonValues = [x * (-0.25) for x in range(80, 200)]
changeSeparation = run(
    changeSeparationToEMinus50('UZ'),
    ICP=['epsilon', 'mExt', 'l2Tau'],
    ISP=0,
    UZSTOP={'epsilon': -32.8},
    e="changeSeparation",
    DS=1.E-4,
    UZR={'epsilon': epsilonValues},
    NMX=12000,
    NPR=0,
    NTST=400,
    NCOL=6
)

changeSeparation += run(
    changeSeparation('UZ'),
    UZSTOP={'epsilon': -31.18},
    DS=1.E-5,
    DSMAX=1.E-1,
)

changeSeparation += run(
    changeSeparation('UZ'),
    ISP=2,
    UZSTOP={'epsilon': -29.0},
    DS=1.E-9,
    DSMAX=1.E-2,
    STOP='LP1'
)

# Merge results in changeSeparation
changeSeparation = merge(changeSeparation)

print("decreasing separation between ribbon ends")

# Define epsilon values for decreasing separation between epsilon=-50 and epsilon=-154
epsilonValues = [x * (-3.) for x in range(17, 52)]

changeSeparation += run(
    changeSeparationToEMinus50('UZ'),
    UZSTOP={'epsilon': -153.08},
    DS=-1.E-2,
    UZR={'epsilon': epsilonValues},
    NMX=40000,
    NPR=0
)

# Define epsilon values for planar ribbon states
epsilonValues = [(-x) for x in range(1, 32)]

changeSeparation += run(
    planarRibState,
    ICP=['epsilon', 'mExt', 'l2Tau'],
    ISP=2,
    UZSTOP={'epsilon': -31.148},
    e="changeSeparation",
    DS=-1.E-3,
    UZR={'epsilon': epsilonValues},
    NMX=12000,
    DSMAX=1.E-1
)

print("Gather solutions in one object and relabel points")
changeSeparation = merge(changeSeparation)
changeSeparation = relabel(changeSeparation)

print("Save and plot")
save(changeSeparation, 'changeSeparation')
bifPlot = plot(changeSeparation, use_labels=False, use_symbols=False)
bifPlot.config(bifurcation_y=['l2Tau'])
bifPlot.config(bifurcation_x=['epsilon'])

# Clean up and wait for the next step
clean()
wait()
