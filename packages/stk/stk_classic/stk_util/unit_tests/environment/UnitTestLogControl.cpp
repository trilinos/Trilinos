/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <sstream>

#include <iostream>
#include <stdexcept>
#include <stk_util/environment/LogControl.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

namespace {

std::ostringstream os;

} // namespace <empty>

STKUNIT_UNIT_TEST(UnitTestLogControl, UnitTest)
{
  stk_classic::RuleMap rule_map;

  rule_map.addLogControlRule("system", stk_classic::LogControlRuleInterval(2));
//  rule_map.addLogControlRule("transient", stk_classic::LogControlRuleInterval(3));
  rule_map.addLogControlRule("nonlinear", stk_classic::LogControlRuleInterval(5));

  int work_count = 0;

  stk_classic::LogControl log_control_i(os, *rule_map.getLogControlRule("system"));
  for (int i = 0; i < 10; ++i) {
    log_control_i.next();

    os << "Running system " << i << std::endl;

    stk_classic::LogControl log_control_j(os, stk_classic::LogControlRuleInterval(3));
    for (int j = 0; j < 10; ++j) {
      log_control_j.next();
        
      os << "  Running transient " << j << std::endl;

      stk_classic::LogControl log_control_k(os, *rule_map.getLogControlRule("nonlinear"));
      for (int k = 0; k < 10; ++k) {
        log_control_k.next();
        
        os << "    Running nonlinear " << k << std::endl;

        os << "      Work count " <<  work_count << std::endl;

        ++work_count;
      }
    }
  }

  std::string result =
    "Running system 0\n"
    "  Running transient 0\n"
    "    Running nonlinear 0\n"
    "      Work count 0\n"
    "    Running nonlinear 5\n"
    "      Work count 5\n"
    "    Running nonlinear 9\n"
    "      Work count 9\n"
    "  Running transient 3\n"
    "    Running nonlinear 0\n"
    "      Work count 30\n"
    "    Running nonlinear 5\n"
    "      Work count 35\n"
    "    Running nonlinear 9\n"
    "      Work count 39\n"
    "  Running transient 6\n"
    "    Running nonlinear 0\n"
    "      Work count 60\n"
    "    Running nonlinear 5\n"
    "      Work count 65\n"
    "    Running nonlinear 9\n"
    "      Work count 69\n"
    "  Running transient 9\n"
    "    Running nonlinear 0\n"
    "      Work count 90\n"
    "    Running nonlinear 5\n"
    "      Work count 95\n"
    "    Running nonlinear 9\n"
    "      Work count 99\n"
    "Running system 2\n"
    "  Running transient 0\n"
    "    Running nonlinear 0\n"
    "      Work count 200\n"
    "    Running nonlinear 5\n"
    "      Work count 205\n"
    "    Running nonlinear 9\n"
    "      Work count 209\n"
    "  Running transient 3\n"
    "    Running nonlinear 0\n"
    "      Work count 230\n"
    "    Running nonlinear 5\n"
    "      Work count 235\n"
    "    Running nonlinear 9\n"
    "      Work count 239\n"
    "  Running transient 6\n"
    "    Running nonlinear 0\n"
    "      Work count 260\n"
    "    Running nonlinear 5\n"
    "      Work count 265\n"
    "    Running nonlinear 9\n"
    "      Work count 269\n"
    "  Running transient 9\n"
    "    Running nonlinear 0\n"
    "      Work count 290\n"
    "    Running nonlinear 5\n"
    "      Work count 295\n"
    "    Running nonlinear 9\n"
    "      Work count 299\n"
    "Running system 4\n"
    "  Running transient 0\n"
    "    Running nonlinear 0\n"
    "      Work count 400\n"
    "    Running nonlinear 5\n"
    "      Work count 405\n"
    "    Running nonlinear 9\n"
    "      Work count 409\n"
    "  Running transient 3\n"
    "    Running nonlinear 0\n"
    "      Work count 430\n"
    "    Running nonlinear 5\n"
    "      Work count 435\n"
    "    Running nonlinear 9\n"
    "      Work count 439\n"
    "  Running transient 6\n"
    "    Running nonlinear 0\n"
    "      Work count 460\n"
    "    Running nonlinear 5\n"
    "      Work count 465\n"
    "    Running nonlinear 9\n"
    "      Work count 469\n"
    "  Running transient 9\n"
    "    Running nonlinear 0\n"
    "      Work count 490\n"
    "    Running nonlinear 5\n"
    "      Work count 495\n"
    "    Running nonlinear 9\n"
    "      Work count 499\n"
    "Running system 6\n"
    "  Running transient 0\n"
    "    Running nonlinear 0\n"
    "      Work count 600\n"
    "    Running nonlinear 5\n"
    "      Work count 605\n"
    "    Running nonlinear 9\n"
    "      Work count 609\n"
    "  Running transient 3\n"
    "    Running nonlinear 0\n"
    "      Work count 630\n"
    "    Running nonlinear 5\n"
    "      Work count 635\n"
    "    Running nonlinear 9\n"
    "      Work count 639\n"
    "  Running transient 6\n"
    "    Running nonlinear 0\n"
    "      Work count 660\n"
    "    Running nonlinear 5\n"
    "      Work count 665\n"
    "    Running nonlinear 9\n"
    "      Work count 669\n"
    "  Running transient 9\n"
    "    Running nonlinear 0\n"
    "      Work count 690\n"
    "    Running nonlinear 5\n"
    "      Work count 695\n"
    "    Running nonlinear 9\n"
    "      Work count 699\n"
    "Running system 8\n"
    "  Running transient 0\n"
    "    Running nonlinear 0\n"
    "      Work count 800\n"
    "    Running nonlinear 5\n"
    "      Work count 805\n"
    "    Running nonlinear 9\n"
    "      Work count 809\n"
    "  Running transient 3\n"
    "    Running nonlinear 0\n"
    "      Work count 830\n"
    "    Running nonlinear 5\n"
    "      Work count 835\n"
    "    Running nonlinear 9\n"
    "      Work count 839\n"
    "  Running transient 6\n"
    "    Running nonlinear 0\n"
    "      Work count 860\n"
    "    Running nonlinear 5\n"
    "      Work count 865\n"
    "    Running nonlinear 9\n"
    "      Work count 869\n"
    "  Running transient 9\n"
    "    Running nonlinear 0\n"
    "      Work count 890\n"
    "    Running nonlinear 5\n"
    "      Work count 895\n"
    "    Running nonlinear 9\n"
    "      Work count 899\n";
  
  STKUNIT_ASSERT_EQUAL((result == os.str()), true);
}

