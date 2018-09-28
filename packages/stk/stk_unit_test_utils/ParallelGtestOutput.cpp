// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "ParallelGtestOutput.hpp"
#include <gtest/gtest.h>            // for TestInfo, UnitTest, etc
#include <gtest/gtest-message.h>
#include <stdarg.h>                 // for va_end, va_list, va_start
#include <stdio.h>                  // for printf, vprintf, fflush, NULL, etc
#include <string>                   // for string
#include "gtest/gtest-test-part.h"  // for TestPartResult
#include "mpi.h"                    // for MPI_Comm, ompi_communicator_t, etc
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################


namespace stk
{
namespace unit_test_util
{
// Macro for referencing flags.
#ifdef GTEST_FLAG
# undef GTEST_FLAG
#endif
#define GTEST_FLAG(name) ::testing::FLAGS_gtest_##name

enum GTestColor {
  COLOR_DEFAULT,
  COLOR_RED,
  COLOR_GREEN,
  COLOR_YELLOW
};

bool should_print_time() {
  return GTEST_FLAG(print_time);
}

const char* GetAnsiColorCode(GTestColor color) {
  switch (color) {
    case COLOR_RED:     return "1";
    case COLOR_GREEN:   return "2";
    case COLOR_YELLOW:  return "3";
    default:            return NULL;
  };
}

void ColoredPrintf(GTestColor color, const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);

    static const bool in_color_mode = true;
    const bool use_color = in_color_mode && (color != COLOR_DEFAULT);

    if (!use_color)
    {
        vprintf(fmt, args);
        va_end(args);
        return;
    }

    printf("\033[0;3%sm", GetAnsiColorCode(color));
    vprintf(fmt, args);
    printf("\033[m"); // Resets the terminal to default.
        va_end(args);
}

class MinimalistPrinter : public ::testing::EmptyTestEventListener
{
public:
    MinimalistPrinter(int procId, MPI_Comm comm) :
            mProcId(procId), mNumFails(0), m_comm(comm)
    {
    }

private:

    int mProcId;
    int mNumFails;
    MPI_Comm m_comm;

    std::string get_filename_for_print(const char* filepath)
    {
        std::string filename = filepath;
        return filename.substr(filename.find_last_of("/") + 1);
    }

    virtual void OnTestStart(const ::testing::TestInfo& test_info)
    {
        if(mProcId == 0)
        {
            printf("*** Starting test %s.%s from %s:%d\n",
                   test_info.test_case_name(),
                   test_info.name(),
                   get_filename_for_print(test_info.file()).c_str(),
                   test_info.line());
        }
    }

    // Called after a failed assertion or a SUCCEED() invocation.
    virtual void OnTestPartResult(
            const ::testing::TestPartResult& test_part_result)
    {
//        if(mProcId == 0)
        {
            printf("%s on proc %d in %s:%d\n%s\n",
                    test_part_result.failed() ? "*** Failure" : "*** Success",
                    mProcId,
                    test_part_result.file_name(),
                    test_part_result.line_number(),
                    test_part_result.summary());
        }
    }

    // Called after a test ends.
    virtual void OnTestEnd(const ::testing::TestInfo& test_info)
    {
        int numFailuresThisProc = 0;
        if(!test_info.result()->Passed())
        {
            numFailuresThisProc = 1;
        }
        int numTotalFailures = -1;
        int root = 0;
        MPI_Reduce(&numFailuresThisProc, &numTotalFailures, 1, MPI_INT, MPI_SUM, root, m_comm);
        if(mProcId == 0)
        {
            if(numTotalFailures == 0)
            {
                ColoredPrintf(COLOR_GREEN, "[       OK ] ");
            }
            else
            {
                int numProcs = -1;
                MPI_Comm_size(m_comm, &numProcs);
                ColoredPrintf(COLOR_RED, "[  FAILED  ] ");
                printf("on %d of %d procs ", numTotalFailures, numProcs);
                mNumFails++;
            }
            printf("%s.%s", test_info.test_case_name(), test_info.name());
            if ( should_print_time() )
            {
                printf(" (%s ms)\n", ::testing::internal::StreamableToString(
                     test_info.result()->elapsed_time()).c_str());
            }
            else
            {
                printf("\n");
            }
            fflush(stdout);
        }
    }

    void OnTestIterationEnd(const ::testing::UnitTest& unit_test, int iteration)
    {
        if(mProcId == 0)
        {
            ColoredPrintf(COLOR_GREEN, "[  PASSED  ] ");
            printf("%d tests.\n", unit_test.test_to_run_count() - mNumFails);

            if(mNumFails > 0)
            {
                ColoredPrintf(COLOR_RED, "[  FAILED  ] ");
                printf("%d tests:\n", mNumFails);
                print_failed_tests(unit_test);
            }
            if ( should_print_time() )
            {
                printf("*** Total elapsed time: %s ms.",
                       ::testing::internal::StreamableToString(unit_test.elapsed_time()).c_str());
            }
            printf("\n");
        }
    }

    void print_failed_tests(const ::testing::UnitTest& unit_test)
    {
        for(int i = 0; i < unit_test.total_test_case_count(); i++)
        {
            const ::testing::TestCase* testCase = unit_test.GetTestCase(i);
            if(testCase->Failed())
                print_failed_tests_in_case(testCase);
        }
    }

    void print_failed_tests_in_case(const ::testing::TestCase* testCase)
    {
        for(int j = 0; j < testCase->total_test_count(); j++)
        {
            const ::testing::TestInfo* testInfo = testCase->GetTestInfo(j);
            if(testInfo->result()->Failed())
                print_failed_test_name(testCase, testInfo);
        }
    }

    void print_failed_test_name(const ::testing::TestCase* testCase, const ::testing::TestInfo* testInfo)
    {
        ColoredPrintf(COLOR_RED, "[  FAILED  ] ");
        printf("--gtest_filter=%s.%s\n", testCase->name(), testInfo->name());
    }
};

void create_parallel_output_with_comm(int procId, MPI_Comm comm)
{
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
    listeners.Append(new MinimalistPrinter(procId, comm));
    delete listeners.Release(listeners.default_result_printer());
}

void create_parallel_output(int procId)
{
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
    listeners.Append(new MinimalistPrinter(procId, MPI_COMM_WORLD));
    delete listeners.Release(listeners.default_result_printer());
}

}
}
