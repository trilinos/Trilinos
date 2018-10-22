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

/*! \file  test_01.cpp
    \brief Test derivative checks.
*/

//#define USE_HESSVEC 0

#include "ROL_StdVector.hpp"
#include "ROL_GetTestProblems.hpp"
#include "ROL_HelperFunctions.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

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

  // Save the format state of the original std::cout.
  ROL::nullstream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  int errorFlag  = 0;

  // Specify interval on which to generate uniform random numbers.
  RealT left = -1.0, right = 1.0;

  // *** Test body.

  try {

    int dim = 128;
    ROL::Ptr<std::vector<RealT> > x_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT> > y_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::Ptr<std::vector<RealT> > z_ptr = ROL::makePtr<std::vector<RealT>>(dim, 0.0);
    ROL::StdVector<RealT> x(x_ptr);
    ROL::StdVector<RealT> y(y_ptr);
    ROL::StdVector<RealT> z(z_ptr);

    // set x,y
    for (int i=0; i<dim; i++) {
      (*x_ptr)[i] = 10.0* (1.0 + (RealT)rand() / (RealT)RAND_MAX);
      (*y_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*z_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
    }

    //ROL::ZOO::Objective_Rosenbrock<RealT> obj;
    ROL::ZOO::Objective_PoissonInversion<RealT> obj(dim,1.e-6);
    //ROL::ZOO::Objective_SumOfSquares<RealT> obj;
    //ROL::ZOO::Objective_LeastSquares<RealT> obj;

    std::vector<std::vector<RealT> > gCheck = obj.checkGradient(x, y);

    for (unsigned i=0; i<gCheck.size(); i++) {
      if (i==0) {
       	*outStream << std::right
                   << std::setw(20) << "Step size"
                   << std::setw(20) << "grad'*dir"
                   << std::setw(20) << "FD approx"
                   << std::setw(20) << "abs error"
       	       	   << "\n";
      }
      *outStream << std::scientific << std::setprecision(8) << std::right
                 << std::setw(20) << gCheck[i][0]
                 << std::setw(20) << gCheck[i][1]
                 << std::setw(20) << gCheck[i][2]
                 << std::setw(20) << gCheck[i][3]
                 << "\n";
    }

    *outStream << "\n";
    std::vector<std::vector<RealT> > hvCheck = obj.checkHessVec(x, y);

    for (unsigned i=0; i<hvCheck.size(); i++) {
      if (i==0) {
        *outStream << std::right
                   << std::setw(20) << "Step size"
                   << std::setw(20) << "norm(Hess*vec)"
                   << std::setw(20) << "norm(FD approx)"
                   << std::setw(20) << "norm(abs error)"
                   << "\n";
      }
      *outStream << std::scientific << std::setprecision(8) << std::right
                 << std::setw(20) << hvCheck[i][0]
                 << std::setw(20) << hvCheck[i][1]
                 << std::setw(20) << hvCheck[i][2]
                 << std::setw(20) << hvCheck[i][3]
                 << "\n";
    }

    *outStream << "\n";
    std::vector<RealT> hsymCheck = obj.checkHessSym(x, y, z);

    *outStream << std::right
               << std::setw(20) << "<w, H(x)v>"
               << std::setw(20) << "<v, H(x)w>"
               << std::setw(20) << "abs error"
               << "\n";
    *outStream << std::scientific << std::setprecision(8) << std::right
               << std::setw(20) << hsymCheck[0]
               << std::setw(20) << hsymCheck[1]
               << std::setw(20) << hsymCheck[2]
               << "\n";

    Teuchos::SerialDenseMatrix<int, RealT> H(x.dimension(), x.dimension());
    H = ROL::computeDenseHessian(obj, x);
    //H.print(*outStream);

    std::vector<std::vector<RealT> > eigenvals = ROL::computeEigenvalues(H);

    *outStream << "\n";
    for (unsigned i=0; i<(eigenvals[0]).size(); i++) {
      if (i==0) {
        *outStream << std::right
                   << std::setw(20) << "Real"
                   << std::setw(20) << "Imag"
                   << "\n";
      }
      *outStream << std::scientific << std::setprecision(8) << std::right
                 << std::setw(20) << (eigenvals[0])[i]
                 << std::setw(20) << (eigenvals[1])[i]
                 << "\n";
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

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return 0;

}

