// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos.hpp"
#include "Stokhos_UnitTestHelpers.hpp"
#include <iomanip>
#ifdef HAVE_STOKHOS_FORUQTK
#include "Stokhos_gaussq.h"
#endif

using std::cout;
using std::setw;
using std::endl;



namespace Stokhos {
using namespace Teuchos;

class JacobiTester
{
public:
  /** */
  JacobiTester(int quadOrder)
    : quad_()
    {
      /* We'll set up a Gauss-Legendre quadrature rule for doing the
       * integrations. We don't want to use Gauss-Jacobi quadrature here,
       * because errors in the Jacobi basis will cause errors in the
       * G-J quadrature as well. */
      Array< RCP<const OneDOrthogPolyBasis<int,double> > > qBases(1);
      qBases[0] = rcp(new LegendreBasis<int,double>(quadOrder));

      RCP<const CompletePolynomialBasis<int,double> > qBasis =
        rcp(new CompletePolynomialBasis<int,double>(qBases));

      quad_ = rcp(new TensorProductQuadrature<int, double>(qBasis, quadOrder));
    }

  /** Compute the inner product 
      \f[
      \int_{-1}^1 (1-x)^\alpha (1+x)^\beta
      P_n^{(\alpha,\beta)}(x) P_m^{(\alpha,\beta)}(x) \, dx
      \f]
      using Gauss-Legendre quadrature and compare to the exact solution. 

      Note: the gamma function appearing in the exact solution is
      computed by the function tgamma(x), not the more logical gamma(x).
  */
  bool testInnerProduct(double alpha, double beta, int nMax) const 
    {
      JacobiBasis<int, double> basis(nMax, alpha, beta, false);

      Array<Array<double> > qp = quad_->getQuadPoints();
      Array<double> qw = quad_->getQuadWeights();
      bool pass = true;

      for (int n=0; n<=nMax; n++)
      {
        double nFact = tgamma(n+1.0);
        cout << "n=" << n << endl;
        for (double x=-1.0; x<=1.0; x+=0.25)
        {
          cout << setw(20) << x << setw(20) << basis.evaluate(x, n)
               << endl;
        }
        for (int m=0; m<=nMax; m++)
        {
          double sum = 0.0;
          for (int q=0; q<qw.size(); q++)
          {
            double x = qp[q][0];
            double w = qw[q] * pow(1-x,alpha)*pow(1+x,beta);
            double Pn = basis.evaluate(x, n);
            double Pm = basis.evaluate(x, m);
            sum += 2.0*w*Pn*Pm;
          }
          double exact = 0.0;
          if (n==m)
            exact = pow(2.0, alpha+beta+1.0)/(2.0*n+alpha+beta+1.0)
              * tgamma(n+alpha+1.0)*tgamma(n+beta+1.0)
              /tgamma(n+alpha+beta+1.0)/nFact;
          double err = fabs(exact - sum);
          cout << setw(4) << n << setw(4) << m 
               << setw(20) << exact << setw(20) << sum << setw(20) << err
               << endl;
          /* Use a fairly loose tolerance because the Gauss-Legendre 
           * quadrature won't be exact when either alpha or beta is not
           * an integer */
          if (err > 1.0e-6) 
          {
            pass = false;
            cout << "***** FAIL ******" << endl;
          }
        }
      }
      return pass;
    }
private:
  RCP<Quadrature<int, double> > quad_;
};



}

int main( int argc, char* argv[] ) 
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  Stokhos::JacobiTester tester(400);
  
  bool allOK = true;
  for (double alpha = 0.75; alpha <= 2.0; alpha += 0.25)
  {
    for (double beta = 0.75; beta <= alpha; beta += 0.25)
    {
      cout << "alpha=" << setw(20) << alpha 
           << " beta=" << setw(20) << beta << endl;
      bool ok = tester.testInnerProduct(alpha, beta, 8);
      allOK = allOK && ok;
    }
  }
  
  if (allOK==true)
  {
    cout << "Jacobi tests PASSED!" << endl;
    cout << "End Result: TEST PASSED" << endl;
    return 0;
  }
  else
  {
    cout << "Jacobi tests FAILED ***" << endl;
    return 1;
  }
}
