import os
import sys

def check_residual(goldResidual, testResidual, maxDiff, relDiff=1e-7, precisionThreshold=1e-7):
    res=False
    if(abs(goldResidual) < precisionThreshold):
        print("Using absolute error")
        # Check residual using absolute error
        # as values are to small to check reliably
        # with relative error due to machine precision
        if(abs(goldResidual - testResidual) < maxDiff):
            res=True
    else:
        print("Using relative error")
        if(abs(goldResidual - testResidual) / abs(goldResidual) < relDiff):
            res=True

    return res

numIterPassed = False
goldConvergenceHistoryPassed = True

goldFileName =  sys.argv[1]
testLogFileName =  sys.argv[2]
maxDiff = float(sys.argv[3])

## Read gold file and extract pure convergence history
goldFile = open(goldFileName, 'rt')
goldFileContent = goldFile.readlines()
goldConvergenceHistory = []
for line in goldFileContent:
    if not line.startswith('#'):
        goldConvergenceHistory.append(line)
goldFile.close()

## Read current log file and extract pure convergence history
testLogFile = open(testLogFileName, 'rt')
testLogFileContent = testLogFile.readlines()
testLogConvergenceHistory = []
for line in testLogFileContent:
    if not line.startswith('#'):
        testLogConvergenceHistory.append(line)
testLogFile.close()

# Check for same number of iterations
if len(goldConvergenceHistory) == len(testLogConvergenceHistory):
    print('Number of iterations: OK')
    numIterPassed = True
else:
    print('Number of iterations: WRONG -- deviation from value in gold file.')

# Check residual in each iteration
if numIterPassed:
    for iteration in range(len(goldConvergenceHistory)):
        goldResidual = float((goldConvergenceHistory[iteration].split())[1])
        testResidual = float((testLogConvergenceHistory[iteration].split())[1])

        if check_residual(goldResidual, testResidual, maxDiff):
            print(('Residual in iteration {}: OK').format(str(iteration)))
        else:
            print(('Residual in iteration {}: WRONG -- deviation from value in gold file.').format(str(iteration)))
            goldConvergenceHistoryPassed = False
            break

# Print final result to be processed by Tribits
if numIterPassed and goldConvergenceHistoryPassed:
    print('End Result: TEST PASSED')
else:
    print('End Result: TEST FAILED')
