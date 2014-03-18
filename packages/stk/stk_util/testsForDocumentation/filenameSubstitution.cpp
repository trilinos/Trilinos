#include <gtest/gtest.h>
#include <string>
#include <mpi.h>
#include <stk_util/diag/FileUtils.hpp>

namespace
{
  TEST(StkUtilHowTo, filenameSubstitution)
  {
    std::string file_name_to_change_in_place = "%B-%P.txt";
    std::string expected_result = "stdin-1.txt";
    
    stk::util::filename_substitution(file_name_to_change_in_place);
    EXPECT_STREQ(file_name_to_change_in_place.c_str(), expected_result.c_str());
  }
}
