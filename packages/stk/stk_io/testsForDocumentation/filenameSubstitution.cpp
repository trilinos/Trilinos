#include <gtest/gtest.h>
#include <string>
#include <mpi.h>
#include <stk_io/Utils.hpp>

namespace
{
  TEST(StkMeshIoBrokerHowTo, filenameSubstitution)
  {
    std::string file_name = "%B-%P.e";
    std::string substituted_name = "stdin-1.e";
    
    stk::io::Utils::filename_substitution(file_name);
    EXPECT_STREQ(file_name.c_str(), substituted_name.c_str());
  }
}
