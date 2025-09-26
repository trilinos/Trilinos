// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "ParallelGtestOutput.hpp"
#include <gtest/gtest.h>            // for TestInfo, UnitTest, etc
#include <gtest/gtest-message.h>
#include <stdarg.h>                 // for va_end, va_list, va_start
#include <stdio.h>                  // for printf, vprintf, fflush, NULL, etc
#include <unistd.h>
#include <stk_util/environment/EnvData.hpp>
#include <stk_util/stk_config.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/util/string_case_compare.hpp>
#include <string>                   // for string
#include "gtest/gtest-test-part.h"  // for TestPartResult

#ifdef STK_HAS_MPI

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

bool should_print_time() {
  return GTEST_FLAG(print_time);
}

class MinimalistPrinter : public ::testing::EmptyTestEventListener
{
public:
    MinimalistPrinter(int procId, MPI_Comm comm) :
            mProcId(procId), mNumFails(0), m_comm(comm)
    {
    }

protected:

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
#ifdef STK_BUILT_FOR_SIERRA
            printf("*** Starting test %s.%s from %s:%d\n",
                   test_info.test_case_name(),
                   test_info.name(),
                   get_filename_for_print(test_info.file()).c_str(),
                   test_info.line());
#else
//older versions of gtest don't have TestInfo::file() nor TestInfo::line()
            printf("*** Starting test %s.%s\n",
                   test_info.test_case_name(),
                   test_info.name());
#endif
        }
    }

    std::string test_result_string(const ::testing::TestPartResult & test_part_result) const
    {
        if (test_part_result.skipped())
        {
            return "*** Skipped";
        }
        else if(test_part_result.failed())
        {
          return  "*** Failure";
        }
        else return "*** Success";
    }

    // Called after a failed assertion or a SUCCEED() invocation.
    virtual void OnTestPartResult(
            const ::testing::TestPartResult& test_part_result)
    {
        printf("%s on proc %d in %s:%d\n%s\n",
                test_result_string(test_part_result).c_str(),
                mProcId,
                test_part_result.file_name(),
                test_part_result.line_number(),
                test_part_result.summary());
    }

    bool test_has_failed(const ::testing::TestInfo & test_info)
    {
        if (test_info.result() == nullptr) 
        {
            return true;
        }

        if (test_info.result()->Skipped() || test_info.result()->Passed())
        {
            return false;
        }

        return true;
    }

    int get_num_failures_across_procs(const ::testing::TestInfo& test_info)
    {
        int numFailuresThisProc = 0;
        if(test_has_failed(test_info))
        {
            numFailuresThisProc = 1;
        }
        int numTotalFailures = -1;
        int root = 0;
        MPI_Reduce(&numFailuresThisProc, &numTotalFailures, 1, MPI_INT, MPI_SUM, root, m_comm);
        return numTotalFailures;
    }

    // Called after a test ends.
    virtual void OnTestEnd(const ::testing::TestInfo& test_info)
    {
        const int numTotalFailures = get_num_failures_across_procs(test_info);
        if(mProcId == 0)
        {
            if(numTotalFailures == 0)
            {
                print_in_color_for_terminal("[       OK ] ", colorGreen);
            }
            else
            {
                int numProcs = stk::parallel_machine_size(m_comm);
                print_failed("on " + std::to_string(numTotalFailures) + " of " + std::to_string(numProcs) + " procs ");
                mNumFails++;
            }
            printf("%s.%s", test_info.test_case_name(), test_info.name());
            if ( should_print_time() )
            {
#ifdef STK_BUILT_FOR_SIERRA
                size_t millis = test_info.result() != nullptr ? test_info.result()->elapsed_time() : 0;
#else
                size_t millis = 0;
#endif
                printf(" (%s ms)\n", ::testing::internal::StreamableToString(millis).c_str());
            }
            else
            {
                printf("\n");
            }
            fflush(stdout);
        }
    }

    void OnTestIterationEnd(const ::testing::UnitTest& unit_test, int /*iteration*/)
    {
        std::vector<std::string> failedTestNames = get_failed_test_names(unit_test);
        collect_failed_test_names_from_all_procs(failedTestNames);

        if(mProcId == 0)
        {
            print_passed(std::to_string(unit_test.test_to_run_count() - mNumFails) + " tests.");
            if(mNumFails > 0)
            {
                print_failed(std::to_string(mNumFails) + " tests:");
                for(const std::string & testName : failedTestNames)
                    print_failed("--gtest_filter=" + testName);
            }
            if ( should_print_time() )
            {
                printf("*** Total elapsed time: %s ms.",
                       ::testing::internal::StreamableToString(unit_test.elapsed_time()).c_str());
            }
            printf("\n");
        }
    }

    std::vector<std::string> get_failed_test_names(const ::testing::UnitTest& unit_test)
    {
        std::vector<std::string> failedTestNames;
        for(int i = 0; i < unit_test.total_test_case_count(); i++)
        {
            const ::testing::TestCase* testCase = unit_test.GetTestCase(i);
            if(testCase->Failed())
                append_failed_test_names_for_case(testCase, failedTestNames);
        }
        return failedTestNames;
    }

    void append_failed_test_names_for_case(const ::testing::TestCase* testCase, std::vector<std::string>& failedTestNames)
    {
        for(int j = 0; j < testCase->total_test_count(); j++)
        {
            const ::testing::TestInfo* testInfo = testCase->GetTestInfo(j);
            if(testInfo->result()->Failed())
                failedTestNames.push_back(std::string {testCase->name()} + "." + testInfo->name());
        }
    }

    void collect_failed_test_names_from_all_procs(std::vector<std::string>& failedTestNames)
    {
        if(stk::parallel_machine_size(m_comm) > 1)
        {
            std::vector<std::string> localNames {failedTestNames};
            stk::parallel_vector_concat(m_comm, localNames, failedTestNames);
            stk::util::sort_and_unique(failedTestNames);
        }
    }

    void print_failed(const std::string &message)
    {
        print_in_color_for_terminal("[  FAILED  ] ", colorRed);
        printf("%s\n", message.c_str());
    }

    void print_passed(const std::string &message)
    {
        print_in_color_for_terminal("[  PASSED  ] ", colorGreen);
        printf("%s\n", message.c_str());
    }

    bool user_requests_color()
    {
      std::string c = GTEST_FLAG(color);
      // follows the convention from gtest itself
      return stk::equal_case(c, "yes") || stk::equal_case(c, "true") || stk::equal_case(c, "t") ||
             stk::equal_case(c, "1");
    }

    bool is_auto_color()
    {
      std::string c = GTEST_FLAG(color);
      // follows the convention from gtest itself
      return stk::equal_case(c, "auto");
    }

    bool is_terminal_output()
    {
      return isatty(fileno(stdout));
    }

    void print_in_color_for_terminal(const std::string &message, const char *colorString)
    {
      if (user_requests_color() || (is_auto_color() && is_terminal_output())) {
        printf("%s%s%s", colorString, message.c_str(), colorDefault);
      } else {
        printf("%s", message.c_str());
      }
    }

    static constexpr const char *colorDefault = "\033[0m";
    static constexpr const char *colorRed = "\033[31m";
    static constexpr const char *colorGreen = "\033[32m";
};

class OutputCapturer
{
public:
    OutputCapturer(stk::ParallelMachine comm)
    {
        originalCout = std::cout.rdbuf();
        std::cout.rdbuf(capturedOutput.rdbuf());

        originalOutputP0 = stk::EnvData::instance().m_outputP0;
        if (stk::parallel_machine_rank(comm) == 0)
            stk::EnvData::instance().m_outputP0 = &capturedOutput;
        else
            stk::EnvData::instance().m_outputP0 = &stk::EnvData::instance().m_outputNull;

        stk::parallel_machine_barrier(comm);
    }
    ~OutputCapturer()
    {
        std::cout.rdbuf(originalCout);
        stk::EnvData::instance().m_outputP0 = originalOutputP0;
    }

    std::string get_and_reset_captured_output() const
    {
        const std::string s = capturedOutput.str();
        capturedOutput.rdbuf()->str("");
        return s;
    }

private:
    std::ostream *originalOutputP0;
    std::streambuf *originalCout;
    std::ostringstream capturedOutput;
};

class MinimalistPrinterOnlyOutputOnFailure : public MinimalistPrinter
{
public:
    using MinimalistPrinter::MinimalistPrinter;

    virtual void OnTestStart(const ::testing::TestInfo& test_info) override
    {
        mOutputCapturer.reset(new OutputCapturer(m_comm));
        MinimalistPrinter::OnTestStart(test_info);
    }

    virtual void OnTestPartResult(const ::testing::TestPartResult& test_part_result) override
    {
        printf("%s", mOutputCapturer->get_and_reset_captured_output().c_str());
        MinimalistPrinter::OnTestPartResult(test_part_result);
    }

    virtual void OnTestEnd(const ::testing::TestInfo& test_info) override
    {
        const int numTotalFailures = get_num_failures_across_procs(test_info);
        if(numTotalFailures > 0)
            printf("%s", mOutputCapturer->get_and_reset_captured_output().c_str());
        MinimalistPrinter::OnTestEnd(test_info);
        mOutputCapturer.reset();
    }

protected:
    std::unique_ptr<OutputCapturer> mOutputCapturer;
};

template <typename PrinterType>
void create_parallel_output_using_printer(int procId, MPI_Comm comm)
{
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
    listeners.Append(new PrinterType(procId, comm));
    delete listeners.Release(listeners.default_result_printer());
}

void create_parallel_output_only_on_failure(int procId, MPI_Comm comm)
{
    create_parallel_output_using_printer<MinimalistPrinterOnlyOutputOnFailure>(procId, comm);
}

void create_parallel_output_with_comm(int procId, MPI_Comm comm)
{
    create_parallel_output_using_printer<MinimalistPrinter>(procId, comm);
}

void create_parallel_output(int procId)
{
    create_parallel_output_with_comm(procId, MPI_COMM_WORLD);
}

}
}

#endif

