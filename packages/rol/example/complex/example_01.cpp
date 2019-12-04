// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Example demonstrating optimzation over a complex field using
           the Rayleigh quotient for a Hermitian matrix as an objective

           Random (complex-valued) initial guesses for \f$z\$f are
           expected to converge to the eigenvector corresponding to 
           the minimum eigenvalue.
*/


#include <cmath>

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_OptimizationSolver.hpp"
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

    auto A = HermitianMatrix<RealT>::example_Matrix();
    auto d = HermitianMatrix<RealT>::example_eigenvalues();
    auto Z = HermitianMatrix<RealT>::example_eigenvectors();
    auto N = A.size();

    os << "Test Matrix A:" << std::endl;
    os << A;

    os << "The eigendecomposition AZ=ZD has the eigenvalues " << std::endl;
    os << "D = diag([" << d[0] << "," << d[1] << "," 
                       << d[2] << "," << d[3] << "])" << std::endl;
 
    ParameterList parlist;
    auto& steplist = parlist.sublist("Step");
    steplist.set("Type", "Trust Region");

    auto& trlist = steplist.sublist("Trust Region");
    trlist.set("Subproblem Solver", "Truncated CG");
        
    auto z   = makePtr<ComplexStdVector<RealT>>(N);
    auto obj = makePtr<RayleighQuotient<RealT>>(A);

    z->randomize();

    auto problem = OptimizationProblem<RealT>(obj,z);
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
    os << "Minimum eigenvalue error: " << std::setprecision(12) << value_error << std::endl;
        
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

