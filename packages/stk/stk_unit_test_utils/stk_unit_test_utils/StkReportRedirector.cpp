/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/


#include <iostream>
#include <stk_unit_test_utils/StkReportRedirector.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace stk
{
namespace unit_test_util
{
namespace {
std::ostringstream &test_stream() {
  static std::ostringstream s_testStringStream;
  return s_testStringStream;
}

void test_report_handler(const char *message, int type) {
  test_stream() << "Message type " << type << ": " << message << std::endl;
}
}  // namespace

StkReportRedirector::StkReportRedirector() {
    test_stream().str({});
    orig_report_handler = stk::set_report_handler(&test_report_handler);
}

std::string StkReportRedirector::get_string() const
{
  return test_stream().str();
}

StkReportRedirector::~StkReportRedirector() { stk::set_report_handler(orig_report_handler); }

}
}
