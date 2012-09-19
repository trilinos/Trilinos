/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_mms/stk_mms.h>
#include <stk_util/diag/UserPlugin.hpp>

namespace stk {
namespace mms {
namespace regression_tests {

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(exactScalarSolution, euler3D)
{
  const double xb=1.234, yb=2.345, zb=2.3, tb = 23.4;

  FAD_Type x(4, 0, xb), y(4, 1, yb), z(4, 2, zb), t(4, 3, tb);

  const std::string file_name = "./mms_user.so";
  const std::string register_function_name = "register_euler_density";
  
  sierra::Plugin::Registry::rootInstance().
    registerDL(file_name.c_str(), register_function_name.c_str());

  const std::string function_name = "euler_density";
  scalar_FAD_function * densitySub 
    = sierra::Plugin::UserSubroutine<scalar_FAD_function>::getFunction(function_name.c_str());

  FAD_Type rho = (*densitySub) (x,y,z,t);

  // check values
  STKUNIT_EXPECT_DOUBLE_EQ(rho.val(), 1.0+(xb-tb)*yb*zb);

  // check derivatives
  STKUNIT_EXPECT_DOUBLE_EQ(rho.dx(_X), yb*zb);  
  STKUNIT_EXPECT_DOUBLE_EQ(rho.dx(_Y), (xb-tb)*zb);  
  STKUNIT_EXPECT_DOUBLE_EQ(rho.dx(_Z), (xb-tb)*yb);  
  STKUNIT_EXPECT_DOUBLE_EQ(rho.dx(_T), -yb*zb);  
}
  
}// namespace regression_tests
}// namespace mms
}// namespace stk

