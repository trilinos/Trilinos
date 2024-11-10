// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_GlobalMPISession.hpp"
#include "compressed_sensing.hpp"
#include "ROL_OptimizationSolver.hpp"
#include <iostream>

using OrdinalT = int;
using RealT = double;
using ROL::Ptr;
using ROL::makePtr;

using VectorT        = ROL::Vector<RealT>;
using ObjectiveT     = ROL::Objective<RealT>;
using ConstraintT    = ROL::Constraint<RealT>;

using TeuchosVectorT = ROL::TeuchosVector<OrdinalT,RealT>;

int main( int argc, char* argv[] ) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  ROL::Ptr<std::ostream> os_ptr;
  ROL::nullstream bhs;
  os_ptr = ROL::makePtrFromRef(argc > 1 ? std::cout : bhs); 
  auto& os = *os_ptr;

  int errorFlag{0};

  try { 
  
    std::string paramfile = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile(paramfile);
    
    int rows{4}, cols{12};
  
    Ptr<ObjectiveT> obj = make_CompressedSensingObjective<OrdinalT,RealT>(rows,cols);
    Ptr<VectorT>    x   = makePtr<TeuchosVectorT>(cols);
  
    auto v = x->clone();
    x->randomize(-1.0,1.0);
    v->randomize(-1.0,1.0);
  
    RealT gamma = 1e-2;
  
    auto opt = ROL::make_l1_penalty_problem( obj, x, gamma );
    ROL::OptimizationSolver<RealT> solver(*opt, *parlist);
    solver.solve(os);
  } catch( std::exception& e ) {
    os << e.what() << std::endl;
    errorFlag = -1000;
  }; // end try

  if(errorFlag==0)
    std::cout << "End Result: TEST PASSED" << std::endl;
  else
    std::cout << "End Result: TEST FAILED" << std::endl;

  return 0;
}
