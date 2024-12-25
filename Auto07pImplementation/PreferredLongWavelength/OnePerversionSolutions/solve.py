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
    UZSTOP={'epsilon': -45.9}
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


print("Use changeSeparation.f90 that sets mExt = 0, and change separation to epsilon=-50.11 afterwards")
changeSeparationToEMinus50 = run(
    runGuessStep('UZ'),
    ICP=['epsilon', 'mExt', 'l2Tau'],
    ISP=0,
    UZSTOP={'epsilon': -50.111},
    DS='-',
    e="changeSeparation"
)

# Define epsilon values to print ribbon states
epsilonValues = [x * (-0.25) for x in range(108, 200)]

changeSeparation = run(
    changeSeparationToEMinus50('UZ'),
    ICP=['epsilon', 'mExt', 'l2Tau'],
    ISP=0,
    UZSTOP={'epsilon': -29.15},
    e="changeSeparation",
    DS=1.E-4,
    UZR={'epsilon': epsilonValues},
    NMX=12000,
    NPR=100
)

changeSeparation += run(
    changeSeparation('UZ'),
    ICP=['epsilon', 'mExt', 'l2Tau'],
    ISP=2,
    UZSTOP={'epsilon': {-20.05, -29.0437}},
    e="changeSeparation",
    DS=1.E-4,
    DSMAX=0.5,
    DSMIN=1.E-14,
    UZR={'epsilon': epsilonValues},
    NMX=120000,
    NPR=400,
    STOP={'LP1'}
)

changeSeparation += run(
    changeSeparation('UZ'),
    ICP=['epsilon', 'mExt', 'l2Tau'],
    ISP=2,
    UZSTOP={'epsilon': {-20.0}},
    e="changeSeparation",
    DS=1.E-9,
    DSMAX=1.E-3,
    DSMIN=1.E-14,
    UZR={'epsilon': epsilonValues},
    NMX=120000,
    NPR=0,
    STOP={'LP1'},
    IAD=1
)

# Merge results in changeSeparation
changeSeparation = merge(changeSeparation)

print("decreasing separation between ribbon ends")

# Define epsilon values for decreasing separation between epsilon=-50 and epsilon=-154
epsilonValues = [x * (-1.) for x in range(50, 154)]

changeSeparation += run(
    changeSeparationToEMinus50('UZ'),
    ICP=['epsilon', 'mExt', 'l2Tau'],
    ISP=0,
    UZSTOP={'epsilon': -153.08},
    e="changeSeparation",
    DS=-1.E-2,
    UZR={'epsilon': epsilonValues},
    NMX=12000,
    NPR=0
)

# Define epsilon values for planar ribbon states
epsilonValues = [(-2. * x) for x in range(1, 15)]

changeSeparation += run(
    planarRibState,
    ICP=['epsilon', 'mExt', 'l2Tau'],
    ISP=2,
    UZSTOP={'epsilon': -29.04},
    e="changeSeparation",
    DS=-1.E-3,
    UZR={'epsilon': epsilonValues},
    NMX=12000,
    NPR=0,
    DSMAX=.1E-1
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
