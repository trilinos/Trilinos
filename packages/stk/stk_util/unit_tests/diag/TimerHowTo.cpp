#include <gtest/gtest.h>
#include <stk_util/diag/PrintTimer.hpp>
#include <stk_util/diag/Timer.hpp>

namespace
{

const double tolerance = 5e-2;

void doWork()
{
    ::usleep(1e5);
}

bool isNear(double a, double b, double tolerance)
{
    bool isNear = false;
    double diff = a - b;
    if(diff > -tolerance && diff < tolerance)
    {
        isNear = true;
    }
    return isNear;
}

bool approximatelyEqualAsNumbers(const std::string &expectedWord, const std::string &actualWord, double tol)
{
    std::istringstream expectedStream(expectedWord);
    std::istringstream actualStream(actualWord);
    double expectedDouble, actualDouble;
    expectedStream >> expectedDouble;
    actualStream >> actualDouble;
    bool isEqual = false;
    if(!expectedStream.fail() && !actualStream.fail() && isNear(expectedDouble, actualDouble, tol))
    {
        isEqual = true;
    }
    return isEqual;
}

bool areStringsEqualWithNumberTolerance(const std::string &expectedString, const std::string &actualString, double tol)
{
    std::istringstream expectedStream(expectedString);
    std::istringstream actualStream(actualString);
    bool isEqual = true;
    while(!expectedStream.fail())
    {
        std::string expectedWord;
        expectedStream >> expectedWord;
        std::string actualWord;
        actualStream >> actualWord;
        if(!expectedStream.fail() && expectedWord != "SKIP")
        {
            if(expectedWord != actualWord && !approximatelyEqualAsNumbers(expectedWord, actualWord, tol))
            {
                std::cerr << "Expected: '" << expectedWord << "' but got '" << actualWord << "'" << std::endl;
                isEqual = false;
            }
        }
    }
    return isEqual;
}

TEST(StkDiagTimerHowTo, useTheSimplestTimer)
{
    unsigned enabledTimerMask = 1;
    stk::diag::Timer rootTimer = createRootTimer("totalTestRuntime", stk::diag::TimerSet(enabledTimerMask));

    stk::diag::TimeBlock totalTestRuntime(rootTimer);
    doWork();

    std::ostringstream outputStream;
    bool printTimingsOnlySinceLastPrint = false;
    stk::diag::printTimersTable(outputStream, rootTimer, stk::diag::METRICS_ALL, printTimingsOnlySinceLastPrint);

    std::string expectedOutput = "                                                                      \
                Timer                                Count       CPU Time              Wall Time        \
            ---------------------------------------- ----- --------------------- ---------------------  \
            totalTestRuntime                              1        0.001 SKIP        0.100 (100.0%)     \
            ";
    EXPECT_TRUE(areStringsEqualWithNumberTolerance(expectedOutput, outputStream.str(), tolerance));

    stk::diag::deleteRootTimer(rootTimer);
}

TEST(StkDiagTimerHowTo, printTwice)
{
    unsigned enabledTimerMask = 1;
    stk::diag::Timer rootTimer = createRootTimer("totalTestRuntime", stk::diag::TimerSet(enabledTimerMask));

    stk::diag::TimeBlock totalTestRuntime(rootTimer);
    doWork();

    std::ostringstream outputStream;
    bool printTimingsOnlySinceLastPrint = false;
    stk::diag::printTimersTable(outputStream, rootTimer, stk::diag::METRICS_ALL, printTimingsOnlySinceLastPrint);

    doWork();
    stk::diag::printTimersTable(outputStream, rootTimer, stk::diag::METRICS_ALL, printTimingsOnlySinceLastPrint);

    std::string expectedOutput = "                                                                      \
                 Timer                   Count       CPU Time              Wall Time                    \
            ---------------------------------------- ----- --------------------- ---------------------  \
            totalTestRuntime                             1        0.000  SKIP        0.100 (100.0%)     \
                                                                                                        \
                             Timer                   Count       CPU Time              Wall Time        \
            ---------------------------------------- ----- --------------------- ---------------------  \
            totalTestRuntime                             1        0.001 SKIP        0.200 (100.0%)      \
            ";
    EXPECT_TRUE(areStringsEqualWithNumberTolerance(expectedOutput, outputStream.str(), tolerance));

    stk::diag::deleteRootTimer(rootTimer);
}

TEST(StkDiagTimerHowTo, useChildTimers)
{
    unsigned enabledTimerMask = 3;
    stk::diag::Timer rootTimer = createRootTimer("totalTestRuntime", stk::diag::TimerSet(enabledTimerMask));
    stk::diag::TimeBlock totalTestRuntime(rootTimer);

    unsigned child1TimerMask = 1;
    const std::string childName1 = "childTimer1";
    stk::diag::Timer childTimer1(childName1, child1TimerMask, rootTimer);
    unsigned child2TimerMask = 2;
    const std::string childName2 = "childTimer2";
    stk::diag::Timer childTimer2(childName2, child2TimerMask, rootTimer);

    {
        stk::diag::TimeBlock timeStuffInThisScope(childTimer1);
        stk::diag::TimeBlock timeStuffInThisScopeAgain(childTimer2);
        doWork();
    }

    std::ostringstream outputStream;
    bool printTimingsOnlySinceLastPrint = false;
    stk::diag::printTimersTable(outputStream, rootTimer, stk::diag::METRICS_ALL, printTimingsOnlySinceLastPrint);

    {
        stk::diag::TimeBlock timeStuffInThisScope(childTimer1);
        doWork();
    }

    stk::diag::printTimersTable(outputStream, rootTimer, stk::diag::METRICS_ALL, printTimingsOnlySinceLastPrint);

    std::string expectedOutput = "                                                                      \
                 Timer                   Count       CPU Time              Wall Time                    \
            ---------------------------------------- ----- --------------------- ---------------------  \
            totalTestRuntime                             1        0.000  SKIP        0.100 SKIP         \
              childTimer1                                1        0.000  SKIP        0.100 SKIP         \
              childTimer2                                1        0.000  SKIP        0.100 SKIP         \
                                                                                                        \
                             Timer                   Count       CPU Time              Wall Time        \
            ---------------------------------------- ----- --------------------- ---------------------  \
            totalTestRuntime                             1        0.000  SKIP        0.200 SKIP         \
              childTimer1                                2        0.000  SKIP        0.200 SKIP         \
              childTimer2                                1        0.000  SKIP        0.100 SKIP         \
            ";
    EXPECT_TRUE(areStringsEqualWithNumberTolerance(expectedOutput, outputStream.str(), tolerance));

    stk::diag::deleteRootTimer(rootTimer);
}

TEST(StkDiagTimerHowTo, disableChildTimers)
{
    unsigned enabledTimerMask = 1;
    stk::diag::Timer rootTimer = createRootTimer("totalTestRuntime", stk::diag::TimerSet(enabledTimerMask));
    stk::diag::TimeBlock totalTestRuntime(rootTimer);

    unsigned enabledCHhildTimerMask = enabledTimerMask;
    unsigned disabledCHhildTimerMask = 2;
    const std::string disabledTimerString = "disabledTimer";
    stk::diag::Timer disabledTimer(disabledTimerString, disabledCHhildTimerMask, rootTimer);
    const std::string enabledTimerString = "enabledTimer";
    stk::diag::Timer enabledTimer(enabledTimerString, enabledCHhildTimerMask, rootTimer);

    {
        stk::diag::TimeBlock timeStuffInThisScope(disabledTimer);
        stk::diag::TimeBlock timeStuffInThisScopeAgain(enabledTimer);
        doWork();
    }

    std::ostringstream outputStream;
    bool printTimingsOnlySinceLastPrint = false;
    stk::diag::printTimersTable(outputStream, rootTimer, stk::diag::METRICS_ALL, printTimingsOnlySinceLastPrint);

    {
        stk::diag::TimeBlock timeStuffInThisScope(disabledTimer);
        doWork();
    }

    stk::diag::printTimersTable(outputStream, rootTimer, stk::diag::METRICS_ALL, printTimingsOnlySinceLastPrint);

    std::string expectedOutput = "                                                                      \
                 Timer                   Count       CPU Time              Wall Time                    \
            ---------------------------------------- ----- --------------------- ---------------------  \
            totalTestRuntime                             1        0.000  SKIP        0.100 SKIP         \
              enabledTimer                               1        0.000  SKIP        0.100 SKIP         \
                                                                                                        \
                             Timer                   Count       CPU Time              Wall Time        \
            ---------------------------------------- ----- --------------------- ---------------------  \
            totalTestRuntime                             1        0.000  SKIP        0.200 SKIP         \
              enabledTimer                               1        0.000  SKIP        0.100 SKIP         \
            ";
    EXPECT_TRUE(areStringsEqualWithNumberTolerance(expectedOutput, outputStream.str(), tolerance));

    stk::diag::deleteRootTimer(rootTimer);
}
}
