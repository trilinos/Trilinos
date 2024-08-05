// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_11.cpp
    \brief Verify that the implementation of the Coleman-Li Trust-Region
           model passes derivative checks 
*/

#include "ROL_ColemanLiModel.hpp"
#include "ROL_HS2.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_RandomVector.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
//#include <fenv.h>

typedef double RealT;

int main(int argc, char *argv[]) {
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  typedef ROL::Vector<RealT>          V;
  typedef ROL::Objective<RealT>       OBJ;
  typedef ROL::BoundConstraint<RealT> CON; 

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  // Save the format state of the original std::cout.
  ROL::nullstream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  RealT zero(0);

  ROL::Ptr<V>   x0;
  ROL::Ptr<V>   x;
  ROL::Ptr<V>   g;
  ROL::Ptr<OBJ> obj;
  ROL::Ptr<CON> con;
  ROL::Ptr<OBJ> model;  

  ROL::ZOO::getHS2<RealT> HS2;
  obj = HS2.getObjective();
  con = HS2.getBoundConstraint();
  x0  = HS2.getInitialGuess();
  x   = HS2.getSolution();

  g = x->dual().clone();

  // Need to evaluate the gradient to construct the model
  obj->gradient(*g,*x,zero);

  model = ROL::makePtr<ROL::ColemanLiModel<RealT>>(*obj,*con,*x,*g);

  ROL::Ptr<V> s = x->clone();
  ROL::Ptr<V> v = x->clone();
  ROL::Ptr<V> u = x->clone();

  ROL::RandomizeVector(*s,-1.0,1.0);
  ROL::RandomizeVector(*u,-1.0,1.0);
  ROL::RandomizeVector(*v,-1.0,1.0);

  model->checkGradient(*s,*v);
  model->checkHessVec(*s,*v);
  model->checkHessSym(*s,*u,*v);

  return 0;
}


