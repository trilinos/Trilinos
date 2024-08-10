// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_02.cpp
    \brief Example demonstrating optimzation over a complex field using
           the Rayleigh quotient for a Hermitian matrix as an objective
           and the constraint that the solution be orthogonal to 
           the eigenvector corresponding to the minimum eigenvalue
*/


#include <cmath>

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_LinearConstraint.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_SingletonVector.hpp"

#include "ROL_ComplexStdVector.hpp"
#include "HermitianMatrix.hpp"
#include "RayleighQuotient.hpp"
#include "OrthogonalProjector.hpp"


int main(int argc, char *argv[]) {

  using namespace  ROL;
  using RealT     = double;              

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Ptr<std::ostream> os_ptr;
  nullstream bhs; // outputs nothing
  if (iprint > 0)
    os_ptr = ROL::makePtrFromRef(std::cout);
  else
    os_ptr = ROL::makePtrFromRef(bhs);

  auto& os = *os_ptr;

  int errorFlag  = 0;

  // *** Example body.
 
  try {

    auto A   = HermitianMatrix<RealT>::example_Matrix();
    auto d   = HermitianMatrix<RealT>::example_eigenvalues();
    auto Z   = HermitianMatrix<RealT>::example_eigenvectors();
    auto N   = A.size();
    auto P   = makePtr<OrthogonalProjector<RealT>>(Z[0]);
    auto v0  = makePtr<ComplexStdVector<RealT>>(N,std::complex<RealT>(0,0));
    auto z   = makePtr<ComplexStdVector<RealT>>(N);
    auto obj = makePtr<RayleighQuotient<RealT>>(A);
    auto con = makePtr<LinearConstraint<RealT>>(P,v0);
    auto l   = makePtr<ComplexStdVector<RealT>>(1);

    os << "Test Matrix A:" << std::endl;
    os << A;

    os << "The eigendecomposition AZ=ZD has the eigenvalues " << std::endl;
    os << "D = diag([" << d[0] << "," << d[1] << "," 
                       << d[2] << "," << d[3] << "])" << std::endl;
 
    ParameterList parlist;
//    auto& steplist = parlist.sublist("Step");
//    steplist.set("Type", "");
//
//    auto& trlist = steplist.sublist("Trust Region");
//    trlist.set("Subproblem Solver", "Truncated CG");
        
    z->randomize();

    auto problem = OptimizationProblem<RealT>(obj,z,con,l);
    auto solver  = OptimizationSolver<RealT>(problem,parlist);
    solver.solve(os); 

    os << "Optimal z:" << std::endl;
    
    z->set(*(problem.getSolutionVector()));
    z->print(os);
    
    auto tol = std::sqrt(ROL_EPSILON<RealT>());

    auto value_opt   = obj->value(*z,tol);
    auto value_error = std::abs(d[0]-value_opt);

    errorFlag += value_error > tol;

    os << "rho(z_opt): " << std::setprecision(12) << value_opt << std::endl;
//    os << "Minimum eigenvalue error: " << std::setprecision(12) << value_error << std::endl;
        

  }
  catch (std::logic_error& err) {
    os << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}

