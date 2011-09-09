/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/PrintTable.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/ExceptionWatch.hpp>
#include <stk_percept/fixtures/Fixture.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <math.h>

namespace stk {
namespace percept {
namespace unit_tests {

//=============================================================================
//=============================================================================
//=============================================================================

#ifndef REDS
STKUNIT_UNIT_TEST(function, stringFunction_xy_basic)
{
  EXCEPTWATCH;

  // start_demo_stringFunction_xy_basic
  double x=1.234, y=2.345, z=0.0;

  StringFunction sfxy(" x - y ");        // an analytic function defined by the expression x - y
  MDArray xyz(3);                        // a rank-one array with the input point - MDArray is a typedef for Intrepid::FieldContainer<double>
  xyz(0) = x;
  xyz(1) = y;
  xyz(2) = z;
  double time = 0.0;
  MDArray result(1);                     // a rank-one array to hold the result
  sfxy(xyz, result, time);               // invoke the string function
  std::cout << result(0) << std::endl;   // this will print the result of -1.111

  // now do it more easily with a helper function
  evalPrint(x,y,z,time, sfxy);           // helper function that will evaluate sfxy at (1.234,2.345,0,0) and print result of -1.111
  // end_demo
}

STKUNIT_UNIT_TEST(function, stringFunction_xy_basic_1)
{
  EXCEPTWATCH;

  StringFunction sfx("x");
  StringFunction sfy("y");
  StringFunction sfxy("x-y");
  double x = 1.234;
  double y = 5.678;
  double z = 0;
  double t = 0;
  double xy = x-y;

  //StringFunction sfxy1== sfx-sfy;
  evalPrint(1,2,3,0, sfxy);
  double vx = eval(x, y, z, t, sfx);
  std::cout << "x = " << x << " vx= " << vx << std::endl;
  STKUNIT_EXPECT_DOUBLE_EQ(vx, x);
  double vy = eval(x, y, z, t, sfy);
  STKUNIT_EXPECT_DOUBLE_EQ(vy, y);
  double vxy = eval(x, y, z, t, sfxy);
  STKUNIT_EXPECT_DOUBLE_EQ(vxy, xy);

}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(function, stringFunction_xy_basic_2)
{
  EXCEPTWATCH;
  //StringFunction sfx("x","sfx", Dimensions(), Dimensions(), 1, 4, 1);
  //         StringFunction sfx("x","sfx", Dimensions<4>(), Dimensions<1>());
  //         StringFunction sfx("x","sfx", Dimensions<4>(), Dimensions<2,2>());

  //StringFunction sfx("x","sfx", Dimensions(3), Dimensions(2,2));

#define TEST_SF_NAMED_ARG 0
#if TEST_SF_NAMED_ARG
  // these should not compile and are tests of the "named argument" proxy
  StringFunction sftestNA1("x","sftestNA1", Dimensions(3), Dimensions(2,3));
  StringFunction sftestNA2("x",std::string("sftestNA2"), Dimensions(3), Dimensions(2,3));
  sftestNA1.setDomainDimensions(Dimensions(3));
  sftestNA2.setDomainDimensions(Dimensions(3));
#else
  StringFunction sftestNA("x", Name("sftestNA"), Dimensions(3), Dimensions(2,3));
  sftestNA.setDomainDimensions(Dimensions(3));
#endif
  StringFunction sftest("x",Name("sftest"), Dimensions(3), Dimensions(2,3));
  MDArray sftest_domain = sftest.getNewDomain();
  MDArray sftest_codomain = sftest.getNewCodomain();
  STKUNIT_EXPECT_EQ(sftest_domain.rank(), 1);
  STKUNIT_EXPECT_EQ(sftest_domain.dimension(0), 3);
  STKUNIT_EXPECT_EQ(sftest_codomain.rank(), 2);
  STKUNIT_EXPECT_EQ(sftest_codomain.dimension(0), 2);
  STKUNIT_EXPECT_EQ(sftest_codomain.dimension(1), 3);

  StringFunction sfx("x",Name("sfx"));
  STKUNIT_EXPECT_EQ(sfx.getDomainDimensions().size(), 1u);
  STKUNIT_EXPECT_EQ(sfx.getDomainDimensions()[0], 3);
  STKUNIT_EXPECT_EQ(sfx.getCodomainDimensions().size(), 1u);
  STKUNIT_EXPECT_EQ(sfx.getCodomainDimensions()[0], 1);

  StringFunction sfy("y");
  StringFunction sfxy("x-y");
  double x = 1.234;
  double y = 5.678;
  double z = 0;
  double t = 0;
  double xy = x-y;
  //StringFunction sfxy1== sfx-sfy;
  evalPrint(1,2,3,0, sfxy);
  double vx = eval(x, y, z, t, sfx);
  std::cout << "x = " << x << " vx= " << vx << std::endl;
  STKUNIT_EXPECT_DOUBLE_EQ(vx, x);
  double vy = eval(x, y, z, t, sfy);
  STKUNIT_EXPECT_DOUBLE_EQ(vy, y);
  double vxy = eval(x, y, z, t, sfxy);
  STKUNIT_EXPECT_DOUBLE_EQ(vxy, xy);

}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(function, stringFunction_test_alias)
{
  EXCEPTWATCH;
  StringFunction sfx("x", Name("sfx"), Dimensions(3), Dimensions(1));
  StringFunction sfy("y", Name("sfy"), Dimensions(3), Dimensions(1));
  StringFunction sfxy("x-y", Name("sfxy"), Dimensions(3), Dimensions(1));
  StringFunction sfembedded("sfxy", Name("sfembedded"), Dimensions(3), Dimensions(1));
  double x = 1.234;
  double y = 5.678;
  double z = 0;
  double t = 0;
  double xy = x-y;
  //StringFunction sfxy1== sfx-sfy;
  evalPrint(1,2,3,0, sfxy);
  double vx = eval(x, y, z, t, sfx);
  std::cout << "x = " << x << " vx= " << vx << std::endl;
  STKUNIT_EXPECT_DOUBLE_EQ(vx, x);
  double vy = eval(x, y, z, t, sfy);
  STKUNIT_EXPECT_DOUBLE_EQ(vy, y);
  double vxy = eval(x, y, z, t, sfxy);
  STKUNIT_EXPECT_DOUBLE_EQ(vxy, xy);

  // embedded function
  std::cout << "sfembedded = ..." << sfembedded << std::endl;
  evalPrint(1,2,3,0, sfembedded);
  std::cout << "sfembedded = " << eval(x, y, z, t, sfembedded) << std::endl;

  double vxy1 = eval(x, y, z, t, sfembedded);
  STKUNIT_EXPECT_DOUBLE_EQ(vxy1, xy);

  sfembedded.addAlias("sfalias");
  StringFunction sftestalias("sfalias", Name("sftestalias"));
  double vxy2 = eval(x, y, z, t, sftestalias);
  std::cout << "sftestalias = " << eval(x, y, z, t, sftestalias) << std::endl;
  STKUNIT_EXPECT_DOUBLE_EQ(vxy2, xy);
}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(function, stringFunction_vector_valued)
{
  EXCEPTWATCH;
  double x = 1.234;
  double y = 5.678;
  double z = 3.456;
  double t = 0;

  bool didCatch = false;
  try {
    StringFunction sfv0("v[0]=x; v[1]=y; v[2]=z; x", Name("sfv"), Dimensions(1,4), Dimensions(1,3) );
    evalVec3Print(1,2,3,0, sfv0);
  }
  catch (...)
  {
    didCatch = true;
    std::cout << "TEST::function::stringFunctionVector: expected to catch this since dom/codomain dimensions should be rank-1" << std::endl;
  }
  STKUNIT_EXPECT_TRUE(didCatch);

  StringFunction sfv("v[0]=x*y*z; v[1]=y; v[2]=z; x", Name("sfv"), Dimensions(3), Dimensions(3) );
  evalVec3Print(1.234, 2.345e-3, 3.456e+5, 0., sfv);
  MDArray vec = evalVec3(x, y, z, t, sfv);
  std::cout << " x = " << x
            << " y = " << y
            << " z = " << z
            << " val = " << (vec[0]*vec[1]*vec[2]) << std::endl;
  STKUNIT_EXPECT_DOUBLE_EQ(vec[0], x*y*z);
  STKUNIT_EXPECT_DOUBLE_EQ(vec[1], y);
  STKUNIT_EXPECT_DOUBLE_EQ(vec[2], z);

}

enum {NPTS = 4};
static double testpoints[NPTS][4] = {
  {0.1234,  -0.5678, 0.9, 0.812},
  {0.1234e-3,  -0.5678e-5, 0.9e+8, 0.812e-4},
  {.101, 102., 10201., 0.0122},
  {0.003, -100001.1, 44.1, 3.}
};

static double testpoints_fd[NPTS][4] = {
  {0.1234,  -0.5678, 0.9, 0.812},
  {0.1234e-3,  -0.5678e-5, 0.9e-3, 0.812e-4},
  {101, 102., 10.2, 0.0122},
  {0.003, .002, -0.0011,0}
};


// the extra name argument is only to get around an intel compiler warning
#define DO_OP_STKUNIT_UNIT_TEST(OP,name)                                             \
{                                                                       \
  double xy ## name            = x OP              y;                   \
  StringFunction sfxy ## name ("x"  QUOTE(OP)  "y");                    \
  StringFunction sfxy1 ## name = sfx OP sfy;                            \
  double vxy ## name           = eval(x, y, z, t, sfxy ## name );       \
  STKUNIT_EXPECT_DOUBLE_EQ(vxy ## name , xy ## name);                           \
  double vxy1 ## name          = eval(x, y, z, t, sfxy1 ## name);       \
  STKUNIT_EXPECT_DOUBLE_EQ(vxy1 ## name, xy ## name);                           \
}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(function, stringFunction_arithmetic_ops)
{
  EXCEPTWATCH;
  for (unsigned ipts = 0; ipts < NPTS; ipts++)
  {
    double x = testpoints[ipts][0];
    double y = testpoints[ipts][1];
    double z = testpoints[ipts][2];
    double t = testpoints[ipts][3];

    // start_demo_stringFunction_arithmetic_ops
    StringFunction sfx("x");                               //python:  sfx = percept.StringFunction("x")
    StringFunction sfy("y");
    StringFunction sfxy(" x - y ");

    // build a new StringFunction from two existing ones
    StringFunction sfxy2 = sfx - sfy;
    double         xy    = x - y;
    double         vxy   = eval(x, y, z, t, sfxy);
    double         vxy2  = eval(x, y, z, t, sfxy2);

    // the two different functions should give the same result
    STKUNIT_EXPECT_DOUBLE_EQ(vxy, vxy2);

    // and they should give the same result as C++
    STKUNIT_EXPECT_DOUBLE_EQ(vxy, xy);

    // end_demo

    StringFunction sfx1("x");
    StringFunction sfy1("y");
    double vx = eval(x, y, z, t, sfx1);
    STKUNIT_EXPECT_DOUBLE_EQ(vx, x);
    double vy = eval(x, y, z, t, sfy1);
    STKUNIT_EXPECT_DOUBLE_EQ(vy, y);

    DO_OP_STKUNIT_UNIT_TEST(-, _minus);
    DO_OP_STKUNIT_UNIT_TEST(+, _plus);
    DO_OP_STKUNIT_UNIT_TEST(*, _mult);
    DO_OP_STKUNIT_UNIT_TEST(/, _div);
  }
}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(function, stringFunction_derivative)
{
  EXCEPTWATCH;
  for (unsigned ipts = 0; ipts < NPTS; ipts++)
  {
    double x = testpoints[ipts][0];
    double y = testpoints[ipts][1];
    double z = testpoints[ipts][2];
    double t = testpoints[ipts][3];

    // start_demo_stringFunction_derivative
    StringFunction sfxy(" x - y ");
    StringFunction dsfxy_y("-1");
    MDArrayString dy(1,1);
    dy(0,0)="y";
    std::string dy1[1][1] = {{"y"}};
    std::cout << "dy1= " << dy1[0][0] << std::endl;
    //Teuchos::RCP<Function> dsfxy_y_1 = sfxy.derivative(MDArrayString_from(dy1));
    //Teuchos::RCP<Function> dsfxy_y_1 = sfxy.derivative(dy);
    Teuchos::RCP<Function> dsfxy_y_1 = sfxy.derivative_test(dy);


    double         dvxy   = eval(x, y, z, t, *dsfxy_y_1);
    double         dvxy1  = eval(x, y, z, t, dsfxy_y);

    // the two different functions should give the same result
    STKUNIT_EXPECT_DOUBLE_EQ(dvxy, dvxy1);

    // and they should give the same result as C++
    STKUNIT_EXPECT_DOUBLE_EQ(dvxy, -1.0);

    // end_demo
  }
}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(function, stringFunction_derivative_1)
{
  EXCEPTWATCH;
  for (unsigned ipts = 0; ipts < NPTS; ipts++)
  {
    double x = testpoints[ipts][0];
    double y = testpoints[ipts][1];
    double z = testpoints[ipts][2];
    double t = testpoints[ipts][3];

    double eps=1.e-6;
    double eps_loc = eps*(std::fabs(x) + std::fabs(y) + std::fabs(z) + std::fabs(t))/4.0;

    // start_demo_stringFunction_derivative_1
    StringFunction sfxy(" x - y ");
    StringFunction dsfxy_grad("v[0]=1; v[1]= -1; v[2]=0", Name("test"), Dimensions(3), Dimensions(3) );
    MDArrayString dxyz(3,1);
    dxyz(0,0)="x"; dxyz(1,0)="y"; dxyz(2,0)="z";
    std::string grad[] = {"1","-1","0"};
    sfxy.setGradientStrings(grad, 3);
    Teuchos::RCP<Function> dsfxy_grad_1  = sfxy.derivative_test(dxyz);
    Teuchos::RCP<Function> dsfxy_grad_fd = sfxy.derivative_test_fd(dxyz, eps_loc);
    Teuchos::RCP<Function> dsfxy_grad_2  = sfxy.derivative(dxyz);

    MDArray dvxy1   = evalVec3(x, y, z, t, *dsfxy_grad_1);
    MDArray dvxy_fd = evalVec3(x, y, z, t, *dsfxy_grad_fd);
    MDArray dvxy2   = evalVec3(x, y, z, t, *dsfxy_grad_2);
    MDArray dvxy    = evalVec3(x, y, z, t, dsfxy_grad);

    // the two different functions should give the same result
    for (int ii = 0; ii < 3; ii++)
    {
      STKUNIT_EXPECT_DOUBLE_EQ(dvxy(ii), dvxy1(ii));
      STKUNIT_EXPECT_DOUBLE_EQ(dvxy(ii), dvxy2(ii));
      STKUNIT_EXPECT_DOUBLE_EQ_APPROX(dvxy(ii), dvxy_fd(ii));
    }

    // and they should give the same result as C++
    STKUNIT_EXPECT_DOUBLE_EQ(dvxy(0), 1.0);
    STKUNIT_EXPECT_DOUBLE_EQ(dvxy(1), -1.0);

    // end_demo
  }
}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(function, stringFunction_derivative_2)
{
  EXCEPTWATCH;
  for (unsigned ipts = 0; ipts < NPTS; ipts++)
  {
    double x = testpoints_fd[ipts][0];
    double y = testpoints_fd[ipts][1];
    double z = testpoints_fd[ipts][2];
    double t = testpoints_fd[ipts][3];

    double eps=1.e-10;
    double eps_loc = eps*(std::fabs(x) + std::fabs(y) + std::fabs(z) + std::fabs(t))/4.0;

    // start_demo_stringFunction_derivative_2
    StringFunction sf(" sin(x*y*z*z) " );
    std::string grad[] = {"y*z*z*cos(x*y*z*z)",  "x*z*z*cos(x*y*z*z)", "2*x*y*z*cos(x*y*z*z)"};
    std::string gradv = "v[0]="+grad[0]+"; v[1]="+grad[1]+" ; v[2]="+grad[2]+";";
    StringFunction dsf_grad(gradv.c_str(), Name("test"), Dimensions(3), Dimensions(3) );
    MDArrayString dxyz(3,1);
    dxyz(0,0)="x"; dxyz(1,0)="y"; dxyz(2,0)="z";
    sf.setGradientStrings(grad, 3);
    Teuchos::RCP<Function> dsf_grad_fd = sf.derivative_test_fd(dxyz, eps_loc);
    Teuchos::RCP<Function> dsf_grad_2  = sf.derivative(dxyz);

    MDArray dv_fd = evalVec3(x, y, z, t, *dsf_grad_fd);
    MDArray dv2   = evalVec3(x, y, z, t, *dsf_grad_2);
    MDArray dv    = evalVec3(x, y, z, t, dsf_grad);

    // the two different functions should give the same result
    for (int ii = 0; ii < 3; ii++)
    {
      //std::cout << "\n ii= " << ii << "\n"<< std::endl;
      STKUNIT_EXPECT_DOUBLE_EQ(dv(ii), dv2(ii));
      if (std::fabs(dv(ii)-dv_fd(ii)) > 0.5*(std::fabs(dv_fd(ii))+std::fabs(dv(ii)))*1.e-6)
      {
        std::cout << "\nii = " << ii << " x= "<<x<<" y= "<<y<<" z= " << z << " expected= " << dv(ii) << " actual= " << dv_fd(ii) << std::endl;
      }

      STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(dv(ii), dv_fd(ii), 1.e-4);
    }


    // end_demo
  }
}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(function, stringFunction_multiplePoints)
{
  EXCEPTWATCH;
  MDArray points(NPTS, 3);
  MDArray output(NPTS, 1);
  MDArray output_expect(NPTS, 1);

  StringFunction sf1("x+y*z");
  for (unsigned ipts = 0; ipts < NPTS; ipts++)
  {
    double x = testpoints[ipts][0];
    double y = testpoints[ipts][1];
    double z = testpoints[ipts][2];
    double t = testpoints[ipts][3];
    points(ipts, 0) = x;
    points(ipts, 1) = y;
    points(ipts, 2) = z;
    //points(ipts, 3) = t;

    //std::cout << "stringFunction_op: ipts= " << ipts << std::endl;

    double vx = eval(x, y, z, t, sf1);
    STKUNIT_EXPECT_DOUBLE_EQ(vx, x+y*z);
    output_expect(ipts, 0) = vx;
  }
  //StringFunction sf2(sf1.getFunctionString().c_str(), Name("sf2"), Dimensions(NPTS, 4), Dimensions(NPTS, 1));
  StringFunction sf2(sf1.getFunctionString().c_str(), Name("sf2"), Dimensions( 3), Dimensions( 1));
  sf2(points, output, 0.0);
  for (unsigned ipts = 0; ipts < NPTS; ipts++)
  {
    STKUNIT_EXPECT_DOUBLE_EQ(output(ipts, 0), output_expect(ipts, 0));
  }
}

//=============================================================================
//=============================================================================
//=============================================================================
STKUNIT_UNIT_TEST(function, stringFunction_expressions)
{
  EXCEPTWATCH;
  /*

    PI, E

    (*this)["exp"] = new CFunction1(std::exp);
    (*this)["ln"] = new CFunction1(std::log);
    (*this)["log"] = new CFunction1(std::log);
    (*this)["log10"] = new CFunction1(std::log10);
    (*this)["pow"] = new CFunction2(std::pow);
    (*this)["sqrt"] = new CFunction1(std::sqrt);
    (*this)["erfc"] = new CFunction1(erfc);
    (*this)["erf"] = new CFunction1(erf);

    (*this)["acos"] = new CFunction1(std::acos);
    (*this)["asin"] = new CFunction1(std::asin);
    (*this)["atan"] = new CFunction1(std::atan);
    (*this)["atan2"] = new CFunction2(std::atan2);
    (*this)["ceil"] = new CFunction1(std::ceil);
    (*this)["cos"] = new CFunction1(std::cos);
    (*this)["cosh"] = new CFunction1(std::cosh);
    (*this)["floor"] = new CFunction1(std::floor);
    (*this)["sin"] = new CFunction1(std::sin);
    (*this)["sinh"] = new CFunction1(std::sinh);
    (*this)["tan"] = new CFunction1(std::tan);
    (*this)["tanh"] = new CFunction1(std::tanh);

    (*this)["abs"] = new CFunction1(std::fabs);
    (*this)["fabs"] = new CFunction1(std::fabs);
    (*this)["deg"] = new CFunction1(deg);
    (*this)["mod"] = new CFunction2(std::fmod);
    (*this)["fmod"] = new CFunction2(std::fmod);
    (*this)["ipart"] = new CFunction1(ipart);
    (*this)["fpart"] = new CFunction1(fpart);
    (*this)["max"] = new CFunction2(max);
    (*this)["min"] = new CFunction2(min);
    (*this)["poltorectx"] = new CFunction2(poltorectx);
    (*this)["poltorecty"] = new CFunction2(poltorecty);
    (*this)["rad"] = new CFunction1(rad);
    (*this)["recttopola"] = new CFunction2(recttopola);
    (*this)["recttopolr"] = new CFunction2(recttopolr);

    OPCODE_UNDEFINED,
    OPCODE_CONSTANT,
    OPCODE_RVALUE,
    OPCODE_STATEMENT,
    OPCODE_ARGUMENT,

    OPCODE_TIERNARY,
    OPCODE_MULTIPLY,
    OPCODE_DIVIDE,
    OPCODE_MODULUS,
    OPCODE_ADD,
    OPCODE_SUBTRACT,
    OPCODE_UNARY_MINUS,
    OPCODE_FUNCTION,

    OPCODE_EQUAL,
    OPCODE_NOT_EQUAL,
    OPCODE_LESS,
    OPCODE_GREATER,
    OPCODE_LESS_EQUAL,
    OPCODE_GREATER_EQUAL,

    OPCODE_UNARY_NOT,
    OPCODE_LOGICAL_AND,
    OPCODE_LOGICAL_OR,

    OPCODE_ASSIGN
  */



#define EXPR_TO_TEST1 (exp(x)+log(x)+log10(x)+pow(x,y)+sqrt(x)+erfc(x)+erf(x)+acos(x)+ \
                       asin(x)+atan(x)+atan2(x,z)+cos(x)+cosh(x)+sin(x)+sinh(x)+tan(x)+tanh(x)+abs(y)+fabs(y))
#define EXPR_TO_TEST2 (x/y*z-t+(4*x)-(1.23e-3/z))
#define EXPR_TO_TEST3 (4 % 2)
#define EXPR_TO_TEST4 (-z)
#define EXPR_TO_TEST5 (exp(E))
#define EXPR_TO_TEST6 (PI)
#define EXPR_TO_TEST7 (atan2(x,z))
#define EXPR_TO_TEST8 (sin(x+y))

#define EXPR_TO_TEST1A (exp(x)+log(x)+log10(x)+pow(x,y)+sqrt(x)+erfc(x)+erf(x)+acos(x)+asin(x))
#define EXPR_TO_TEST1B (atan(x)+atan2(x,z)+cos(x)+cosh(x)+sin(x)+sinh(x)+tan(x)+tanh(x)+abs(y)+fabs(y))


#define DO_SF_STKUNIT_UNIT_TEST(expr)                                                \
  { using namespace std;                                                \
    const char* str= QUOTE(expr);                                       \
    StringFunction sf(str);                                             \
    double v_loc = eval(x, y, z, t, sf);                                \
    double ve_loc = expr;                                               \
    /* std::cout << "expr= " << str << " x = " << x << " ve= " << ve_loc << " v= " << v_loc << std::endl; */ \
    STKUNIT_EXPECT_DOUBLE_EQ(v_loc, ve_loc);                                    \
                                                                        \
  }


#define DO_SF_TEST_IPT(expr,ipt)                        \
  if(0) std::cout << "ipt= " << ipt << std::endl;       \
  DO_SF_STKUNIT_UNIT_TEST(expr)

  {
    double x = 0.1234;
    double y = -0.5678;
    double z = 0.9;
    double t = 0.812;
    double PI = M_PI;
    double E = M_E;

    StringFunction sf1("x+y");
    double ve=x+y;
    double v = eval(x,y,z,t,sf1);
    std::cout << "x= " << x << " y= " << y << " v= " << v << " ve= " << ve << std::endl;

    DO_SF_STKUNIT_UNIT_TEST(EXPR_TO_TEST1);
    DO_SF_STKUNIT_UNIT_TEST(EXPR_TO_TEST2);
    DO_SF_STKUNIT_UNIT_TEST(EXPR_TO_TEST3);
    DO_SF_STKUNIT_UNIT_TEST(EXPR_TO_TEST4);
    DO_SF_STKUNIT_UNIT_TEST(EXPR_TO_TEST5);
    DO_SF_STKUNIT_UNIT_TEST(EXPR_TO_TEST6);
  }


  for (unsigned ipts = 0; ipts < NPTS; ipts++)
  {
    double x = testpoints[ipts][0];
    double y = testpoints[ipts][1];
    double z = testpoints[ipts][2];
    double t = testpoints[ipts][3];
    double PI = M_PI;
    double E = M_E;

    DO_SF_TEST_IPT(EXPR_TO_TEST1,ipts);
    DO_SF_TEST_IPT(EXPR_TO_TEST2,ipts);
    DO_SF_TEST_IPT(EXPR_TO_TEST3,ipts);
    DO_SF_TEST_IPT(EXPR_TO_TEST4,ipts);
    DO_SF_TEST_IPT(EXPR_TO_TEST5,ipts);
    DO_SF_TEST_IPT(EXPR_TO_TEST6,ipts);
  }

  //exit(1);
}

#define DO_SF_TIMING_TEST_STRING(expr,numIt)                    \
{ using namespace std;                                          \
  const char* str= QUOTE(expr);                                 \
  StringFunction sf(str);                                       \
  double val=0.0;                                               \
  for (unsigned it=0; it < numIt; it++)                         \
  { val += eval(x, y, z, t, sf); }                              \
  valGlobal += val; /* force compiler to do something */        \
  if (0) std::cout << "string val= " << val << std::endl;       \
}

#define DO_SF_TIMING_TEST_CPP(expr,numIt)                       \
{ using namespace std;                                          \
  double val=0.0;                                               \
  for (unsigned it=0; it < numIt; it++)                         \
  { val += expr; }                                              \
  if (0) std::cout << "cpp val= " << val << std::endl;          \
  valGlobal += val; /* force compiler to do something */        \
}


#define COL_SEP "|"
//#define COL_SEP " "
#define EXPR_CELL_WIDTH (80)
#define TIME_IT1(expr,numIt)                                            \
{                                                                       \
  double cpp_time=0;                                                    \
  double string_time=0;                                                 \
  TIME_IT(DO_SF_TIMING_TEST_CPP(expr,numIt),cpp_time);                  \
  TIME_IT(DO_SF_TIMING_TEST_STRING(expr,numIt),string_time);            \
  if (0) std::cout << "for expression= " << QUOTE(expr) << " timings= " << std::endl; \
  if (0) std::cout << "cpp_time= " << cpp_time << " string_time= " << string_time << " ratio= " << string_time/cpp_time << std::endl; \
  table << COL_SEP <<  stk::cell_width(EXPR_CELL_WIDTH) << QUOTE(expr) << COL_SEP << cpp_time << COL_SEP << string_time << COL_SEP << (string_time/cpp_time)<< COL_SEP << stk::end_row; \
}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(function, stringFunction_timing)
{
  EXCEPTWATCH;
  //unsigned numIt = 1024*1024;
  //unsigned numIt = 100*1024;
  unsigned numIt = 1024;
  double valGlobal=0;
  for (unsigned ipts = 0; ipts < NPTS; ipts++)
  {
    double x = testpoints[ipts][0];
    double y = testpoints[ipts][1];
    double z = testpoints[ipts][2];
    double t = testpoints[ipts][3];
    //double PI = M_PI;
    //double E = M_E;

    stk::PrintTable table;
    std::ostringstream msg; msg << "Timings for " << numIt << " iterations for  point # " << ipts << " x,y,z,t= " << x << " "<< y << " "<< z << " "<< t << " ";
    table.setTitle(msg.str());

    table << COL_SEP << stk::cell_width(EXPR_CELL_WIDTH) << "expression" << COL_SEP
          << "cpp time" << COL_SEP
          << "string time" << COL_SEP
          << "ratio" << COL_SEP
          << stk::end_header;

    TIME_IT1(EXPR_TO_TEST1A, numIt);
    TIME_IT1(EXPR_TO_TEST1B, numIt);
    TIME_IT1(EXPR_TO_TEST2, (numIt*100));
    TIME_IT1(EXPR_TO_TEST8, numIt);

    std::cout << table;
  }
  std::cout << "valGlobal= " << valGlobal << std::endl;
}
#endif

} // namespace unit_tests
} // namespace percept
} // namespace stk

