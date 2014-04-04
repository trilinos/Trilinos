#include <gtest/gtest.h>
#include <string>
#include <mpi.h>
#include <stk_util/environment/ProgramOptions.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <stk_util/environment/FileUtils.hpp>

namespace
{
  TEST(StkUtilHowTo, useFilenameSubstitutionWithNoCommandLineOptions)
  {
    const std::string default_base_filename = "stdin";
    const std::string numProcsString = "1";
    const std::string expected_filename = default_base_filename + "-" + numProcsString + ".e";

    std::string file_name = "%B-%P.e";
    stk::util::filename_substitution(file_name);
    EXPECT_EQ(expected_filename, file_name);
  }

  void setFilenameInCommandLineOptions(const std::string &filename)
  {
    boost::program_options::variables_map &command_line_options = stk::get_variables_map();
    command_line_options.insert(std::make_pair("input-deck", boost::program_options::variable_value(filename, false)));
    stk::EnvData::instance().m_inputFile = filename;
  }
  TEST(StkUtilHowTo, useFilenameSubstitutionWithFileComingFromCommandLineOptions)
  {
    const std::string base_filename = "myfile";
    const std::string full_filename = "/path/to/" + base_filename + ".g";
    setFilenameInCommandLineOptions(full_filename);

    const std::string numProcsString = "1";
    const std::string expected_filename = base_filename + "-" + numProcsString + ".e";

    std::string file_name = "%B-%P.e";
    stk::util::filename_substitution(file_name);

    EXPECT_EQ(expected_filename, file_name);
  }
}
