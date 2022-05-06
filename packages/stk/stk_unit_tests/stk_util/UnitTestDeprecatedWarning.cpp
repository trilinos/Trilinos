/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "gtest/gtest.h"
#include "stk_util/environment/Env.hpp"
#include "stk_util/environment/RuntimeDoomed.hpp"
#include "stk_util/registry/DeprecationWarning.hpp"
#include "stk_util/registry/VersionNumber.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include <ostream>
#include <sstream>
#include <string>

namespace tftk
{
namespace util
{

namespace
{
std::ostringstream & test_stream()
{
  static std::ostringstream s_testStringStream;

  return s_testStringStream;
}

void test_report_handler(const char * message, int type)
{
  test_stream() << "Message type " << type << ": " << message << std::endl;
}

class DeprecatedWarning : public ::testing::Test
{
public:
  DeprecatedWarning()
  {
    stk::reset_doomed_count();
    test_stream().str({});
    orig_report_handler = stk::set_report_handler(&test_report_handler);

    orig_cerr = std::cerr.rdbuf(redirected_cerr.rdbuf());
  }

  ~DeprecatedWarning()
  {
    stk::set_report_handler(orig_report_handler);
    std::cerr.rdbuf(orig_cerr);
  }

  void set_current_version(int major, int minor)
  {
    stk::util::VersionNumber::set_current_version(stk::util::VersionNumber{major, minor});
  }

  std::stringstream redirected_cerr;

private:
  stk::REH orig_report_handler;
  std::streambuf * orig_cerr;
};
}

TEST_F(DeprecatedWarning, BeforeErrorVersion)
{
  set_current_version(4, 48);
  stk::reset_doomed_count();

  stk::util::VersionNumber removal_version(4, 51);
  stk::util::DeprecationWarning(removal_version) << "    my feature";

  EXPECT_FALSE(stk::is_doomed());
  std::string expected_message(
      "Deprecated feature removed in Version 4.51 detected.\n    my feature\n");
  auto actual_message = test_stream().str();
  EXPECT_TRUE(actual_message.find(expected_message) != std::string::npos);
  if( sierra::Env::parallel_rank() == 0 ) {
    EXPECT_EQ("WARNING: " + expected_message, redirected_cerr.str());
  }
  else {
    EXPECT_EQ("", redirected_cerr.str());
  }
}

TEST_F(DeprecatedWarning, ErrorVersion)
{
  set_current_version(4, 48);

  stk::util::VersionNumber removal_version(4, 47);
  stk::util::DeprecationWarning(removal_version) << "    my feature";

  EXPECT_TRUE(stk::is_doomed());
  std::string expected_message(
      "Deprecated feature removed in Version 4.47 detected.\n    my feature\n");
  EXPECT_EQ("Message type 1: " + expected_message, test_stream().str());
}

TEST_F(DeprecatedWarning, SetInReleaseVersion)
{
  set_current_version(4, 48);

  stk::util::VersionNumber removal_version(4, 48);
  EXPECT_ANY_THROW(stk::util::DeprecationWarning(removal_version) << "    my feature");
}
}
}
