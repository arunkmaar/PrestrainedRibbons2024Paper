# Start of the guess step
print("Starting guess step")

# Load the initial guess step
guessStep = load("guessStep")

# Save the initial state with external moment mExt = 0 for later use
planarRibState = guessStep

print("Run guess step with increasing sinusoidal external moment (amplitude specified in c.guessStep)")
runGuessStep = run(guessStep)

print("At fixed external moment, decreasing separation between ribbon ends to epsilon = -57.9")
runGuessStep = run(
    runGuessStep('UZ'),
    ICP=['epsilon', 'mExt', 'l2Tau'],
    UZSTOP={'epsilon': -57.9}
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

print("Use changeSeparation.f90 that sets mExt = 0, and change separation to epsilon=-55.30 afterwards")
pointToRestart = run(
    runGuessStep('UZ'),
    ICP=['epsilon', 'mExt', 'l2Tau'],
    UZSTOP={'epsilon': -55.30},
    e="changeSeparation",
    DS=1.E-4,
    NMX=1200,
    NPR=10,
    ISP=2
)

changeSeparation = run(
    pointToRestart('UZ'),
    ISP=2,
    UZSTOP={'epsilon': -55.2552},
    NPR=1,
    DS=1.E-10,
    DSMIN=1.E-15,
    DSMAX=2.E-3,
    STOP='LP1'
)


# Merge results in changeSeparation
changeSeparation = merge(changeSeparation)

print("changing separation between ribbon ends")
# Define epsilon values at which we want to print the ribbon states
epsilonValues = [x * (-0.5) for x in range(100, 250)]

changeSeparation += run(
    pointToRestart('UZ'),
    UZSTOP={'epsilon': -85.11},
    DS='-',
    UZR={'epsilon': epsilonValues}
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
