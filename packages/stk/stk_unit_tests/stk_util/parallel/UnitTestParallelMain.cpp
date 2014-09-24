// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <gtest/gtest.h>                // for TestInfo, etc
#include <mpi.h>                        // for MPI_Comm_rank, MPI_Finalize, etc
#include <stdarg.h>                     // for va_end, va_list, va_start
#include <stdio.h>                      // for printf, vprintf, fflush, etc
#include "gtest/gtest-test-part.h"      // for TestPartResult

int gl_argc = 0;
char** gl_argv = 0;

enum GTestColor {
  COLOR_DEFAULT,
  COLOR_RED,
  COLOR_GREEN,
  COLOR_YELLOW
};

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
    MinimalistPrinter(int procId) :
            mProcId(procId), mNumFails(0)
    {
    }

private:

    int mProcId;
    int mNumFails;

    virtual void OnTestStart(const ::testing::TestInfo& test_info)
    {
        if(mProcId == 0)
        {
            printf("*** Starting test %s.%s.\n",
                    test_info.test_case_name(), test_info.name());
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
        MPI_Reduce(&numFailuresThisProc, &numTotalFailures, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
        if(mProcId == 0)
        {
            if(numTotalFailures == 0)
            {
                ColoredPrintf(COLOR_GREEN, "[       OK ] ");
            }
            else
            {
                int numProcs = -1;
                MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
                ColoredPrintf(COLOR_RED, "[  FAILED  ] ");
                printf("on %d of %d procs ", numTotalFailures, numProcs);
                mNumFails++;
            }
            printf("%s.%s\n", test_info.test_case_name(), test_info.name());
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
                printf("%d tests.\n", mNumFails);
            }
        }
    }
};

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int procId = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &procId);

    testing::InitGoogleTest(&argc, argv);

    gl_argc = argc;
    gl_argv = argv;

    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
    listeners.Append(new MinimalistPrinter(procId));
    delete listeners.Release(listeners.default_result_printer());

    int returnVal = RUN_ALL_TESTS();

    MPI_Finalize();

    return returnVal;
}
