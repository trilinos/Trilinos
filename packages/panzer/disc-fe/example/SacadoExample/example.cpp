// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <vector>
#include <iostream>
#include <unistd.h>

#include "Sacado.hpp"

typedef Sacado::Fad::DFad<double> FadType;
typedef Sacado::Fad::DFad<Sacado::Fad::SFad<double,1> > SecondFadType;

template <typename Scalar>
Scalar func(const Scalar & x,const Scalar & y)
{
  // Grad = y * std::cos(x*y)
  //      = x * std::cos(x*y)-0.25*std::sin(y)

  // Hess*v = dy * std::cos(x*y) - y * ( dx * sin(x*y) + dy * sin(x*y))
  //        = dx * std::cos(x*y) - x * ( dx * sin(x*y) + dy * sin(x*y)) - 0.25 * dy * std::cos(y)
  return std::sin(x*y)+0.25*std::cos(y);
}

void grad_func(const double & x,const double & y,double g[2])
{
  g[0] = y * std::cos(x*y);
  g[1] = x * std::cos(x*y)-0.25*std::sin(y);
}

template <typename Scalar>
void func_dx(const Scalar & x,const Scalar & y,Scalar final[2])
{
  typedef Sacado::Fad::DFad<Scalar> Fad;
 
  Fad x_fad; x_fad.resizeAndZero(2);
  Fad y_fad; y_fad.resizeAndZero(2);
 
  x_fad.val() = x;
  x_fad.fastAccessDx(0) = 1.0;

  y_fad.val() = y;
  y_fad.fastAccessDx(1) = 1.0;

  Fad result = func(x_fad,y_fad);

  final[0] = result.fastAccessDx(0);
  final[1] = result.fastAccessDx(1);
}

std::vector<double> hess_func(double x, double y,double dx, double dy)
{
  std::vector<double> v(2);
  v[0] = dy * std::cos(x*y) - y * ( y*dx * sin(x*y) + x*dy * sin(x*y));
  v[1] = dx * std::cos(x*y) - x * ( y*dx * sin(x*y) + x*dy * sin(x*y)) - 0.25 * dy * std::cos(y);
  
  return v;
}

/** Seed to compute a second directional derivative.This is used
  * to compute (given a function f)
  * 
  * H(x) * v = d/dt(\nabla f(x+t*v))|_{t=0}
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
inline SecondFadType seed_second_deriv(int num_vars,int index,double xi,double vi)
{
  SecondFadType x = SecondFadType(num_vars,index,xi);
  x.val().fastAccessDx(0) = vi;
  
  return x;
}

int main(int /* argc */, char*[] /* argv[] */)
{
  std::cout << "Outputting" << std::endl;

  // directions
  double x_val = 0.25;
  double y_val = 0.5;
  double dx = 2.0;
  double dy = 3.0;

  // value
  {
    double x = x_val;
    double y = y_val;
    double f = func(x,y);

    std::cout << "Value:  " << f << std::endl;
  }

  std::cout << std::endl;

  // first derivative
  {
    FadType x = FadType(2,0,x_val);
    FadType y = FadType(2,1,y_val);
    FadType f = func(x,y);

    std::cout << "First:  " << f << std::endl;

    double gf[2];
    grad_func(x_val,y_val,gf);
    std::cout << "First Exact:  " << gf[0] << " " << gf[1] << std::endl;
  }

  std::cout << std::endl;

  // second derivative
  {
    SecondFadType x = seed_second_deriv(2,0,x_val,dx);
    SecondFadType y = seed_second_deriv(2,1,y_val,dy);
    SecondFadType f = func(x,y);

    std::cout << "Second: " << f << std::endl;

    std::vector<double> hess =  hess_func(x_val,y_val,dx,dy);
    std::cout << "Second Exact: " << hess[0]  << " " << hess[1] << std::endl;
  }

  return 0;
}
