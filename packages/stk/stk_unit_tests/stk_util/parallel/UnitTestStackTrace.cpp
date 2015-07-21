#include <gtest/gtest.h>
#include "stk_util/parallel/DebugTool.hpp"


std::string knownFunction()
{
    return getStackTrace();
}

TEST(getStackTrace, testStackTrace)
{
#if defined(__GNUC__) && !defined(__ICC) && !defined(NDEBUG)
    // for gcc debug testing only
    std::string trace = knownFunction();
    std::string nameOfFunction = "knownFunction";
    std::size_t found = trace.find(nameOfFunction);
    EXPECT_TRUE(found != std::string::npos);
#endif
}

