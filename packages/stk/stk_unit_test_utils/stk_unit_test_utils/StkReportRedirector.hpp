/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef STK_STK_UNIT_TEST_UTILS_STKREPORTREDIRECTOR_HPP_
#define STK_STK_UNIT_TEST_UTILS_STKREPORTREDIRECTOR_HPP_
#include <stk_util/util/ReportHandler.hpp>
#include "stk_util/stk_config.h"

namespace stk
{
namespace unit_test_util
{

class StkReportRedirector {
 public:
  StkReportRedirector();
  ~StkReportRedirector();
  std::string get_string() const;

 private:
  stk::REH orig_report_handler;
};

namespace simple_fields {

class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
StkReportRedirector : public stk::unit_test_util::StkReportRedirector {};

} // namespace simple_fields

}
}


#endif /* STK_STK_UNIT_TEST_UTILS_STKREPORTREDIRECTOR_HPP_ */
