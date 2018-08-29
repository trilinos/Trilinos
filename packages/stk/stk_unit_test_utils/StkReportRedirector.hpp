/*--------------------------------------------------------------------*/
/*    Copyright 2003, 2017 Sandia Corporation.                        */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef STK_STK_UNIT_TEST_UTILS_STKREPORTREDIRECTOR_HPP_
#define STK_STK_UNIT_TEST_UTILS_STKREPORTREDIRECTOR_HPP_
#include <stk_util/util/ReportHandler.hpp>


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

}
}


#endif /* STK_STK_UNIT_TEST_UTILS_STKREPORTREDIRECTOR_HPP_ */
