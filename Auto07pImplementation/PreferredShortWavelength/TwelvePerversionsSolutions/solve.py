# Start of the guess step
print("Starting guess step")

# Load the initial guess step
guessStep = load("guessStep")

# Save the initial state with external moment mExt = 0 for later use
planarRibState = guessStep

print("Run guess step with increasing sinusoidal external moment (amplitude specified in c.guessStep)")
runGuessStep = run(guessStep)

print("At fixed external moment, decreasing separation between ribbon ends to epsilon = -58.9")
runGuessStep = run(
    runGuessStep('UZ'),
    ICP=['epsilon', 'mExt', 'l2Tau'],
    UZSTOP={'epsilon': -58.9}
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

# Define epsilon values at which we want to print the ribbon states
epsilonValues = [x * (-0.25) for x in range(200, 450)]

print("Use changeSeparation.f90 that sets mExt = 0, and change separation to epsilon=-57.8 afterwards")
changeSeparationToEMinus57d8 = run(
    runGuessStep('UZ'),
    ICP=['epsilon', 'mExt', 'l2Tau'],
    ISP=0,
    UZSTOP={'epsilon': -57.8},
    DS=1.E-4,
    UZR={'epsilon': epsilonValues},
    NMX=1200,
    e="changeSeparation"
)

changeSeparation = run(
    changeSeparationToEMinus57d8('UZ'),
    ISP=2,
    UZSTOP={'epsilon': -52.69},
    NPR=1,
    STOP='BP28'
)


# Merge results in changeSeparation
changeSeparation = merge(changeSeparation)

print("decreasing separation between ribbon ends")

# Define epsilon values for decreasing separation
epsilonValues = [-x for x in range(50, 100)]

changeSeparation += run(
    changeSeparationToEMinus57d8('UZ'),
    UZSTOP={'epsilon': -97.053},
    DS='-'
)

# Define epsilon values for planar ribbon states
epsilonValues = [(-5. * x) for x in range(1, 15)]
tempEpsilonValues = [(0.75- x) for x in range(50, 100)]
epsilonValues.extend(tempEpsilonValues)

changeSeparation += run(
    planarRibState,
    ICP=['epsilon', 'mExt', 'l2Tau'],
    ISP=2,
    UZSTOP={'epsilon': -54.9830},
    e="changeSeparation",
    UZR={'epsilon': epsilonValues},
    DS=-1.E-4,
    DSMAX=1E-2,
    STOP='BP3'
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
