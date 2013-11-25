#include <gtest/gtest.h>
#include <mpi.h>

int gl_argc=0;
char** gl_argv=0;
int proc_id=-1;

class MinimalistPrinter : public ::testing::EmptyTestEventListener
{
    // Called before a test starts.
    virtual void OnTestStart(const ::testing::TestInfo& test_info)
    {
        if (proc_id == 0 )
        {
            printf("*** Test %s.%s starting.\n",
             test_info.test_case_name(), test_info.name());
        }
    }

    // Called after a failed assertion or a SUCCEED() invocation.
    virtual void OnTestPartResult(const ::testing::TestPartResult& test_part_result)
    {
        if (proc_id == 0)
        {
            printf("%s in %s:%d\n%s\n",
             test_part_result.failed() ? "*** Failure" : "Success",
             test_part_result.file_name(),
             test_part_result.line_number(),
             test_part_result.summary());
        }
    }

    // Called after a test ends.
    virtual void OnTestEnd(const ::testing::TestInfo& test_info)
    {
        if (proc_id == 0)
        {
            printf("*** Test %s.%s ending.\n",
             test_info.test_case_name(), test_info.name());
        }
    }
};

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    testing::InitGoogleTest(&argc, argv);

    gl_argc = argc;
    gl_argv = argv;

    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
    delete listeners.Release(listeners.default_result_printer());
    listeners.Append(new MinimalistPrinter);

    int returnVal = RUN_ALL_TESTS();
    
    MPI_Finalize();
    return returnVal;
}
