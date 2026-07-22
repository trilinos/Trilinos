// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include "Sacado.hpp"

typedef Sacado::Fad::DFad<double> FadType;
typedef Sacado::Fad::SFad<FadType,1> HessianType;

/** Seed to compute a second directional derivative.This is used
  * to compute (given a function f)
  *
  *    H(x) * v = d/dt(\nabla f(x+t*v))|_{t=0}
  *       
  * where x is the value to differentiate about, v is the direction. This function
  * seeds for one value of the independent variable x[i] in the direction v[i]. What
  * is returned from this function is the independent variable to pass the derivative
  * to.
  *
  * \param[in] num_vars Size of the vector x (and v)
  * \param[in] index Index of variable that this is seeding
  * \param[in] xi Value of the x vector to seed with
  * \param[in] vi Value of the v vector to seed with
  */
inline HessianType seed_second_deriv(int num_vars,int index,double xi,double vi)
{
  typedef HessianType SecondFadType;

  SecondFadType x = SecondFadType(1,FadType(num_vars,index,xi));
  x.fastAccessDx(0) = vi;

  return x;
}

//**********************************************************************
double func(double x, double y)
{
  // Grad = y * std::cos(x*y)
  //      = x * std::cos(x*y)-0.25*std::sin(y)
  
  // Hess*v = dy * std::cos(x*y) - y * ( dx * sin(x*y) + dy * sin(x*y))
  //        = dx * std::cos(x*y) - x * ( dx * sin(x*y) + dy * sin(x*y)) - 0.25 * dy * std::cos(y)
  //
  return std::sin(x*y)+0.25*std::cos(y);
}

std::vector<double> hess_func(double x, double y,double dx, double dy)
{
  std::vector<double> v(2);
  v[0] = dy * std::cos(x*y) - y * ( y*dx * sin(x*y) + x*dy * sin(x*y));
  v[1] = dx * std::cos(x*y) - x * ( y*dx * sin(x*y) + x*dy * sin(x*y)) - 0.25 * dy * std::cos(y);

  return v;
}

//**********************************************************************
TEUCHOS_UNIT_TEST(hessian_test_k,correctness)
{
  typedef HessianType ScalarT;
 
  using Teuchos::RCP;
  using Teuchos::rcp;


  double x_val = 0.25;
  double y_val = 0.5;
  double dx_val = 2.0;
  double dy_val = 3.0;

  PHX::View<ScalarT*> x("x",5);
  PHX::View<ScalarT*> y("y",5);
  PHX::View<ScalarT*> dx("dx",5);
  PHX::View<ScalarT*> dy("dy",5);
  PHX::View<ScalarT*> result("result",5);

  for(int i=0;i<5;++i) {
    dx(i) = ScalarT(dx_val);
    dy(i) = ScalarT(dy_val);
    x(i) = seed_second_deriv(2,0,x_val,dx_val);
    y(i) = seed_second_deriv(2,1,y_val,dy_val);
  }

  for(int i=0;i<5;++i)
    result(i) = std::sin(x(i)*y(i))+0.25*std::cos(y(i));

  for(int i=0;i<5;i++) {
    double x_val  = Sacado::scalarValue(x(i));
    double y_val  = Sacado::scalarValue(y(i));
    double dx_val = Sacado::scalarValue(dx(i));
    double dy_val = Sacado::scalarValue(dy(i));
    double f = func(x_val,y_val);
    std::vector<double> hess = hess_func(x_val,y_val,dx_val,dy_val);

    ScalarT r = result(i);

    TEST_EQUALITY(Sacado::scalarValue(r),f);
    TEST_EQUALITY(r.fastAccessDx(0).fastAccessDx(0),hess[0]);
    TEST_EQUALITY(r.fastAccessDx(0).fastAccessDx(1),hess[1]);

    out << "RESULT = " << r << std::endl;
  }
}
