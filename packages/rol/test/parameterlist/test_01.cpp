// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ROL_ParameterList.hpp"
#include <iostream>


int main( int argc, char* argv[] ) {

  auto os_ptr = ROL::makeStreamPtr(std::cout, argc);
  auto& os = *os_ptr;
  int errorFlag = 0;

  using ROL::detail::operator<<;

  try {

    auto parlist = ROL::ParameterList();

 
    auto& genlist = parlist.sublist("General");
    genlist.set("Inexact Objective Function",false);
    genlist.set("Inexact Gradient",false);
    genlist.set("Inexact Hessian-Times-A-Vector",false);

    {
      std::vector<int> ivalues = { 1,2,3,4 };
      std::string key = "Example Integer Values";
      genlist.set(key,ivalues);
//      auto value = genlist.get<std::vector<int>>(key);
      auto value = ROL::getArrayFromStringParameter<int>(genlist,key);
      os << key << " = " << value << std::endl;
    }

    auto& secant = genlist.sublist("Secant");
    secant.set("Type","Limited-Memory BFGS");
    secant.set("Use as Preconditioner",true);
    secant.set("Use as Hessian",false);
    secant.set("Maximum Storage",10);
    
    auto& krylov = genlist.sublist("Krylov");
    {
      std::string key("Iteration Limit");
      auto maxit = krylov.get(key,10);
      os << key << " = " << maxit << std::endl;
    }

  krylov.set("Type","Conjugate Gradients");
    krylov.set("Absolute Tolerance",1.e-4);
    krylov.set("Relative Tolerance",1.e-2);
    krylov.set("Iteration Limit",100);

   
    auto& step = parlist.sublist("Step");
    auto& ls = step.sublist("Line Search");
    ls.set("Function Evaluation Limit", 20);
    ls.set("Sufficient Decrease Tolerance", 1.e-4);

    auto& dm = ls.sublist("Descent Method");
    dm.set("Type","Newton-Krylov");

    auto& cc = ls.sublist("Curvature Condition");
    cc.set("Type","Strong Wolfe Conditions");
    cc.set("General Parameter","0.9");
    cc.set("Generalized Wolfe Parameter","0.6");

    auto& lsm = ls.sublist("Line-Search Method");
    lsm.set("Type","Cubic Interpolation");
    lsm.set("Backtracking Rate","0.5");

    {
      std::string key("Iteration Limit");
      auto maxit = krylov.get<int>(key,10);
      os << key << " = " << maxit << std::endl;
    }
//    {
//      std::string key("Type");
//      auto type = krylov.get<std::string>(key);
//      os << key << " = " << type << std::endl;
//    }

   
    os << parlist; 
  } catch( std::exception& e ) {
    errorFlag = -1000;
    os << e.what() << std::endl;
  };
  
  if( errorFlag ) std::cout << "End Result: TEST FAILED\n";
  else            std::cout << "End Result: TEST PASSED\n"; 

  return 0;
}

