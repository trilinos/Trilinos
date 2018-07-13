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
    \brief Shows how to solve the inverse Poisson problem using trust-region
           methods with dense Hessian diagnostics.
*/

#define USE_HESSVEC 1

#include "ROL_Types.hpp"
#include "ROL_PoissonInversion.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>
#include <algorithm>

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  // *** Example body.

  try {

    int dim = 128; // Set problem dimension.
    ROL::ZOO::Objective_PoissonInversion<RealT> obj(dim, 1e-6);

    // Define algorithm.
    ROL::ParameterList parlist;
    std::string stepname = "Trust Region";
    parlist.sublist("Step").sublist(stepname).set("Subproblem Solver", "Truncated CG");
    parlist.sublist("General").sublist("Krylov").set("Iteration Limit",50);
    parlist.sublist("General").sublist("Krylov").set("Relative Tolerance",1e-2);
    parlist.sublist("General").sublist("Krylov").set("Absolute Tolerance",1e-4);
    parlist.sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist.sublist("Status Test").set("Step Tolerance",1.e-14);
    parlist.sublist("Status Test").set("Iteration Limit",100);
    ROL::Algorithm<RealT> algo(stepname,parlist);

    // Iteration vector.
    ROL::Ptr<std::vector<RealT> > x_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    // Set initial guess.
    for (int i=0; i<dim; i++) {
      (*x_ptr)[i] = 0.1;
    }
    ROL::StdVector<RealT> x(x_ptr);

    // Run algorithm.
    algo.run(x, obj, true, *outStream);

    // Compute dense Hessian matrix. 
    Teuchos::SerialDenseMatrix<int, RealT> H(x.dimension(), x.dimension());
    H = ROL::computeDenseHessian<RealT>(obj, x);
    //H.print(*outStream);

    // Compute and print eigenvalues.
    std::vector<std::vector<RealT> > eigenvals = ROL::computeEigenvalues<RealT>(H);

    *outStream << "\nEigenvalues:\n";
    for (unsigned i=0; i<(eigenvals[0]).size(); i++) {
      if (i==0) {
        *outStream << std::right
                   << std::setw(28) << "Real"
                   << std::setw(28) << "Imag"
                   << "\n";
      }
      *outStream << std::scientific << std::setprecision(16) << std::right
                 << std::setw(28) << (eigenvals[0])[i]
                 << std::setw(28) << (eigenvals[1])[i]
                 << "\n";
    }

    // Compute and print generalized eigenvalues.
    Teuchos::SerialDenseMatrix<int, RealT> M = computeDotMatrix(x);
    //M.print(*outStream);
    std::vector<std::vector<RealT> > genEigenvals = ROL::computeGenEigenvalues<RealT>(H, M);

    *outStream << "\nGeneralized eigenvalues:\n";
    for (unsigned i=0; i<(genEigenvals[0]).size(); i++) {
      if (i==0) {
        *outStream << std::right
                   << std::setw(28) << "Real"
                   << std::setw(28) << "Imag"
                   << "\n";
      }
      *outStream << std::scientific << std::setprecision(16) << std::right
                 << std::setw(28) << (genEigenvals[0])[i]
                 << std::setw(28) << (genEigenvals[1])[i]
                 << "\n";
    }

    // Sort and compare eigenvalues and generalized eigenvalues - should be close.
    std::sort((eigenvals[0]).begin(), (eigenvals[0]).end());
    std::sort((eigenvals[1]).begin(), (eigenvals[1]).end());
    std::sort((genEigenvals[0]).begin(), (genEigenvals[0]).end());
    std::sort((genEigenvals[1]).begin(), (genEigenvals[1]).end());

    RealT errtol = std::sqrt(ROL::ROL_EPSILON<RealT>());
    for (unsigned i=0; i<(eigenvals[0]).size(); i++) {
      if ( std::abs( (genEigenvals[0])[i] - (eigenvals[0])[i] ) > errtol*((eigenvals[0])[i]+ROL::ROL_THRESHOLD<RealT>()) ) {
        errorFlag++;
        *outStream << std::scientific << std::setprecision(20) << "Real genEigenvals - eigenvals (" << i << ") = " << std::abs( (genEigenvals[0])[i] - (eigenvals[0])[i] ) << " > " << errtol*((eigenvals[0])[i]+1e4*ROL::ROL_THRESHOLD<RealT>()) << "\n";
      }
      if ( std::abs( (genEigenvals[1])[i] - (eigenvals[1])[i] ) > errtol*((eigenvals[1])[i]+ROL::ROL_THRESHOLD<RealT>()) ) {
        errorFlag++;
        *outStream << std::scientific << std::setprecision(20) << "Imag genEigenvals - eigenvals (" << i << ") = " << std::abs( (genEigenvals[1])[i] - (eigenvals[1])[i] ) << " > " << errtol*((eigenvals[1])[i]+ROL::ROL_THRESHOLD<RealT>()) << "\n";
      }
    }

    // Compute inverse of Hessian.
    Teuchos::SerialDenseMatrix<int, RealT> invH = ROL::computeInverse<RealT>(H);
    Teuchos::SerialDenseMatrix<int, RealT> HinvH(H);

    // Multiply with Hessian and verify that it gives the identity (l2 dot matrix M from above).
    HinvH.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, H, invH, 0.0);
    //*outStream << std::scientific << std::setprecision(6); HinvH.print(*outStream);
    HinvH -= M;
    if (HinvH.normOne() > errtol) {
        errorFlag++;
        *outStream << std::scientific << std::setprecision(20) << "1-norm of H*inv(H) - I = " << HinvH.normOne() << " > " << errtol << "\n";      
    }

    // Use Newton algorithm with line search.
    stepname = "Line Search";
    parlist.sublist("Step").sublist(stepname).sublist("Descent Method").set("Type", "Newton-Krylov");
    ROL::Algorithm<RealT> newton_algo(stepname,parlist);

    // Reset initial guess.
    for (int i=0; i<dim; i++) {
      (*x_ptr)[i] = 0.1;
    }

    // Run Newton algorithm.
    newton_algo.run(x, obj, true, *outStream);

    ROL::Ptr<const ROL::AlgorithmState<RealT> > new_state = newton_algo.getState();
    ROL::Ptr<const ROL::AlgorithmState<RealT> > old_state = algo.getState();
    *outStream << "old_optimal_value = " << old_state->value << std::endl;
    *outStream << "new_optimal_value = " << new_state->value << std::endl;
    if ( std::abs(new_state->value - old_state->value) / std::abs(old_state->value) > errtol ) {
        errorFlag++;
        *outStream << std::scientific << std::setprecision(20) << "\nabs(new_optimal_value - old_optimal_value) / abs(old_optimal_value)  = " << std::abs(new_state->value - old_state->value) / std::abs(old_state->value) << " > " << errtol << "\n";      
    }

  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}

