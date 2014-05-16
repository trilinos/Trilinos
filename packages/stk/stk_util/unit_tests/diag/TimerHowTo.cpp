#include <gtest/gtest.h>
#include <bitset>
#include <stk_util/diag/PrintTimer.hpp>
#include <stk_util/diag/Timer.hpp>
#include <comparison/stringAndNumberComparisons.h>

namespace
{

const double tolerance = 5e-2;

void doWork()
{
    ::usleep(1e5);
}

TEST(StkDiagTimerHowTo, useTheRootTimer)
{
    stk::diag::TimerSet enabledTimerSet(0);
    stk::diag::Timer rootTimer = createRootTimer("totalTestRuntime", enabledTimerSet);

    stk::diag::TimeBlock totalTestRuntime(rootTimer);
    doWork();

    std::ostringstream outputStream;
    bool printTimingsOnlySinceLastPrint = false;
    stk::diag::printTimersTable(outputStream, rootTimer, stk::diag::METRICS_ALL, printTimingsOnlySinceLastPrint);

    std::string expectedOutput = "                                                                      \
                Timer                                Count       CPU Time              Wall Time        \
            ---------------------------------------- ----- --------------------- ---------------------  \
            totalTestRuntime                              1        0.001 SKIP        0.100 SKIP     \
            ";
    EXPECT_TRUE(unitTestUtils::areStringsEqualWithToleranceForNumbers(expectedOutput, outputStream.str(), tolerance));

    stk::diag::deleteRootTimer(rootTimer);
}

TEST(StkDiagTimerHowTo, useChildTimers)
{
    enum {CHILDMASK1 = 1, CHILDMASK2 = 2};
    stk::diag::TimerSet enabledTimerSet(CHILDMASK1 | CHILDMASK2);
    stk::diag::Timer rootTimer = createRootTimer("totalTestRuntime", enabledTimerSet);
    stk::diag::TimeBlock totalTestRuntime(rootTimer);

    const std::string childName1 = "childTimer1";
    stk::diag::Timer childTimer1(childName1, CHILDMASK1, rootTimer);
    const std::string childName2 = "childTimer2";
    stk::diag::Timer childTimer2(childName2, CHILDMASK2, rootTimer);

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
    EXPECT_TRUE(unitTestUtils::areStringsEqualWithToleranceForNumbers(expectedOutput, outputStream.str(), tolerance));

    stk::diag::deleteRootTimer(rootTimer);
}

TEST(StkDiagTimerHowTo, disableChildTimers)
{
    enum {CHILDMASK1 = 1, CHILDMASK2 = 2};
    stk::diag::TimerSet enabledTimerSet(CHILDMASK2);
    stk::diag::Timer rootTimer = createRootTimer("totalTestRuntime", enabledTimerSet);
    stk::diag::TimeBlock totalTestRuntime(rootTimer);

    const std::string disabledTimerString = "disabledTimer";
    stk::diag::Timer disabledTimer(disabledTimerString, CHILDMASK1, rootTimer);
    const std::string enabledTimerString = "enabledTimer";
    stk::diag::Timer enabledTimer(enabledTimerString, CHILDMASK2, rootTimer);

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
    EXPECT_TRUE(unitTestUtils::areStringsEqualWithToleranceForNumbers(expectedOutput, outputStream.str(), tolerance));

    stk::diag::deleteRootTimer(rootTimer);
}

}
