// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "BrusselatorProblemInterface.H"
#include "LOCA_Parameter_Vector.H"
#include "NOX_LAPACK_Vector.H"
#include "NOX_LAPACK_Matrix.H"
#include "LOCA_GlobalData.H"
#include "NOX_Utils.H"

BrusselatorProblemInterface::BrusselatorProblemInterface(
		    const Teuchos::RCP<LOCA::GlobalData>& global_data,
		    int N, double a, 
		    double b, 
		    double d1,
		    double d2,
		    std::ofstream& file)  : 
  globalData(global_data),
  initialGuess(2*N),
  alpha(a),
  beta(b),
  D1(d1),
  D2(d2),
  n(N),
  outputFile(file)
{
  for (int i=0; i<n; i++) {
    initialGuess(i) = alpha;
    initialGuess(n+i) = beta/alpha;
  }
}

const NOX::LAPACK::Vector&
BrusselatorProblemInterface::getInitialGuess()
{
  return initialGuess;
}

bool
BrusselatorProblemInterface::computeF(NOX::LAPACK::Vector& f, 
				      const NOX::LAPACK::Vector &x)
{
  double h = 1.0 / double(n-1);
  double hh = h*h;

  f(0) = x(0) - alpha;
  for (int i=1; i<n-1; i++)
    f(i) = D1*(x(i-1) - 2.0*x(i) + x(i+1)) / hh
      + alpha - (beta + 1.0)*x(i) + x(i)*x(i)*x(i+n);
  f(n-1) = x(n-1) - alpha;

  f(n) = x(n) - beta/alpha;
  for (int i=n+1; i<2*n-1; i++)
    f(i) = D2*(x(i-1) - 2.0*x(i) + x(i+1)) / hh
      + beta*x(i-n) - x(i-n)*x(i-n)*x(i);
  f(2*n-1) = x(2*n-1) - beta/alpha;
  
  return true;
}

bool
BrusselatorProblemInterface::computeJacobian(NOX::LAPACK::Matrix<double>& J, 
					     const NOX::LAPACK::Vector & x)
{
  double h = 1.0 / double(n-1);
  double hh = h*h;

  J(0,0) = 1.0;
  for (int i=1; i<n-1; i++) {
    J(i,i-1) = D1/hh;
    J(i,i) = -2.0*D1/hh - (beta + 1.0) + 2.0*x(i)*x(i+n);
    J(i,i+1) = D1/hh;
    J(i,i+n) = x(i)*x(i);
  }
  J(n-1,n-1) = 1.0;

  J(n,n) = 1.0;
  for (int i=n+1; i<2*n-1; i++) {
    J(i,i-n) = beta - 2.0*x(i-n)*x(i);
    J(i,i-1) = D2/hh;
    J(i,i) = -2.0*D2/hh - x(i-n)*x(i-n);
    J(i,i+1) = D2/hh;
  }
  J(2*n-1,2*n-1) = 1.0;

  return true;
}

bool
BrusselatorProblemInterface::computeShiftedMatrix(
				    double alpha, double beta,
				    const NOX::LAPACK::Vector& x,
				    NOX::LAPACK::Matrix<double>& A)
{
  bool res = true;
  if (alpha != 0.0) {
    res = computeJacobian(A, x);
    A.scale(alpha);
  }
  else
    A.scale(0.0);
  if (beta != 0.0) {
    // Note:  the mass matrix terms for the BCs are zero
    for (int i=1; i<n-1; i++) {
      A(i,i) += beta;
      A(i+n,i+n) += beta;
    }
  }
  return res;
}

void
BrusselatorProblemInterface::setParams(const LOCA::ParameterVector& p) {
  alpha = p.getValue("alpha");
  beta = p.getValue("beta");
  D1 = p.getValue("D1");
  D2 = p.getValue("D2");
}

void
BrusselatorProblemInterface::printSolution(const NOX::LAPACK::Vector &x,
					   const double conParam)
{
   if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails)) {
     globalData->locaUtils->out() 
       << "At parameter value: " 
       << globalData->locaUtils->sciformat(conParam)
       << "   the solution vector is\n";
     
     for (int i=0; i<6; i++) 
       globalData->locaUtils->out() 
	 << " " << globalData->locaUtils->sciformat(x(i));
     globalData->locaUtils->out() << " ...";
     for (int i=n-2; i<n; i++)  
       globalData->locaUtils->out() 
	 << " " << globalData->locaUtils->sciformat(x(i));
     globalData->locaUtils->out() << std::endl;
     
     for (int i=n; i<n+6; i++) 
       globalData->locaUtils->out() 
	 << " " << globalData->locaUtils->sciformat(x(i));
     globalData->locaUtils->out() << " ...";
     for (int i=2*n-2; i<2*n; i++)  
       globalData->locaUtils->out() 
	 << " " << globalData->locaUtils->sciformat(x(i));
     globalData->locaUtils->out() << std::endl; 
   }

   outputFile << conParam << " ";
   for (int i=0; i<2*n; i++)
     outputFile << x(i) << " ";
   outputFile << std::endl << std::endl;
}

