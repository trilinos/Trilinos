// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_BurkardtRulesDef.hpp
    \brief  Definition file for integration rules provided by John Burkardt.
	    <A HREF="http://people.sc.fsu.edu/~jburkardt/cpp_src/sandia_rules/sandia_rules.html">
	    <\A>
    \author Created by D. Kouri and D. Ridzal.
 */

#ifdef _MSC_VER
#include "winmath.h"
#endif

//# include <cstdlib>
//# include <iomanip>
//# include <iostream>
//# include <cmath>
//# include <ctime>
//# include <Teuchos_ScalarTraitsDecl.hpp>

namespace Intrepid {

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::chebyshev1_compute ( int n, Scalar x[], Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    CHEBYSHEV1_COMPUTE computes a Chebyshev type 1 quadrature rule.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X) / sqrt ( 1 - x^2 ) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the order.
//    1 <= N.
//
//    Output, Scalar X[N], the abscissas.
//
//    Output, Scalar W[N], the weights.
//
{
  if (n<1) {
    std::cerr << "\n";
    std::cerr << "CHEBYSHEV1_COMPUTE - Fatal error!\n";
    std::cerr << "  Illegal value of N = " << n << "\n";
    std::exit (1);
  }

  for (int i=0;i<n;i++) {
    w[i] = M_PI/(Scalar)(n);
  }
  for (int i=0;i<n;i++) {
    x[i] = std::cos(M_PI*(Scalar)(2*n-1-2*i)/(Scalar)(2*n));
  }
  if ((n%2)==1) {
    x[(n-1)/2] = 0.0;
  }

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::chebyshev1_compute_points ( int n, Scalar x[] )
//****************************************************************************
//
//  Purpose:
//
//    CHEBYSHEV1_COMPUTE_POINTS computes Chebyshev type 1 quadrature points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the order.
//    1 <= N.
//
//    Output, Scalar X[N], the abscissas.
//
{
  if (n<1) {
    std::cerr << "\n";
    std::cerr << "CHEBYSHEV1_COMPUTE_POINTS - Fatal error!\n";
    std::cerr << "  Illegal value of N = " << n << "\n";
    std::exit(1);
  }

  for (int i=0;i<n;i++) {
    x[i] = std::cos(M_PI*(Scalar)(2*n-1-2*i)/(Scalar)(2*n));
  }
  if ((n%2)==1) {
    x[(n-1)/2] = 0.0;
  }

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::chebyshev1_compute_weights ( int n, Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    CHEBYSHEV1_COMPUTE_WEIGHTS computes Chebyshev type 1 quadrature weights.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the order.
//    1 <= N.
//
//    Output, Scalar W[N], the weights.
//
{
  if (n<1) {
    std::cerr << "\n";
    std::cerr << "CHEBYSHEV1_COMPUTE_WEIGHTS - Fatal error!\n";
    std::cerr << "  Illegal value of N = " << n << "\n";
    std::exit(1);
  }

  for (int i=0;i<n;i++) {
    w[i] = M_PI/(Scalar)n;
  }

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::chebyshev2_compute ( int n, Scalar x[], Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    CHEBYSHEV2_COMPUTE computes a Chebyshev type 2 quadrature rule.
//
//  Discussion:
//
//    The integral:
//
//      integral ( -1 <= x <= 1 ) f(x)  sqrt ( 1 - x^2 )  dx
//
//    The quadrature rule:
//
//      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the order.
//    1 <= N.
//
//    Output, Scalar X[N], the abscissas.
//
//    Output, Scalar W[N], the weights.
//
{
  Scalar angle;

  if (n<1) {
    std::cerr << "\n";
    std::cerr << "CHEBYSHEV2_COMPUTE - Fatal error!\n";
    std::cerr << "  Illegal value of N = " << n << "\n";
    std::exit(1);
  }

  for (int i=0;i<n;i++) {
    angle = M_PI*(Scalar)(n-i)/(Scalar)(n+1);
    w[i]  = M_PI/(Scalar)(n+1)*std::pow(std::sin(angle),2);
    x[i]  = std::cos(angle);
  }

  if ((n%2)==1) {
    x[(n-1)/2] = 0.0;
  }

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::chebyshev2_compute_points ( int n, Scalar x[] )
//****************************************************************************
//
//  Purpose:
//
//    CHEBYSHEV2_COMPUTE_POINTS computes Chebyshev type 2 quadrature points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the order.
//    1 <= N.
//
//    Output, Scalar X[N], the abscissas.
//
{
  Scalar angle;

  if (n<1) {
    std::cerr << "\n";
    std::cerr << "CHEBYSHEV2_COMPUTE_POINTS - Fatal error!\n";
    std::cerr << "  Illegal value of N = " << n << "\n";
    std::exit(1);
  }

  for (int i=0;i<n;i++) {
    angle = M_PI*(Scalar)(n-i)/(Scalar)(n+1);
    x[i]  = std::cos(angle);
  }

  if ((n%2)==1) {
    x[(n-1)/2] = 0.0;
  }

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::chebyshev2_compute_weights ( int n, Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    CHEBYSHEV2_COMPUTE_WEIGHTS computes Chebyshev type 2 quadrature weights.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the order.
//    1 <= N.
//
//    Output, Scalar W[N], the weights.
//
{
  Scalar angle;

  if (n<1) {
    std::cerr << "\n";
    std::cerr << "CHEBYSHEV2_COMPUTE_WEIGHTS - Fatal error!\n";
    std::cerr << "  Illegal value of N = " << n << "\n";
    std::exit(1);
  }

  for (int i=0;i<n;i++) {
    angle = M_PI*(Scalar)(n-i)/(Scalar)(n+1);
    w[i]  = M_PI/(Scalar)(n+1)*std::pow(std::sin(angle),2);
  }

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::clenshaw_curtis_compute ( int n, Scalar x[], Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis quadrature rule.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//    1 <= N.
//
//    Output, Scalar X[N], the abscissas.
//
//    Output, Scalar W[N], the weights.
//
{
  Scalar b, theta;
  int i, j;

  if (n<1) {
    std::cerr << "\n";
    std::cerr << "CLENSHAW_CURTIS_COMPUTE - Fatal error!\n";
    std::cerr << "  Illegal value of N = " << n << "\n";
    std::exit(1);
  }
  else if (n==1) {
    x[0] = 0.0;
    w[0] = 2.0;
  }
  else {
    for (i=0;i<n;i++) {
      x[i] = std::cos((Scalar)(n-1-i)*M_PI/(Scalar)(n-1));
    }
    x[0] = -1.0;
    if ((n%2)==1) {
      x[(n-1)/2] = 0.0;
    }
    x[n-1] = +1.0;

    for (i=0;i<n;i++) {
      theta = (Scalar)i*M_PI/(Scalar)(n-1);

      w[i] = 1.0;

      for (j=1;j<=(n-1)/2;j++) {
        if (2*j==(n-1)) {
          b = 1.0;
        }
        else {
          b = 2.0;
        }

        w[i] = w[i]-b*std::cos(2.0*(Scalar)(j)*theta)/(Scalar)(4*j*j-1);
      }
    }

    w[0] = w[0]/(Scalar)(n-1);
    for (i=1;i<n-1;i++) {
      w[i] = 2.0*w[i]/(Scalar)(n-1);
    }
    w[n-1] = w[n-1]/(Scalar)(n-1);
  }

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::clenshaw_curtis_compute_points ( int n, Scalar x[] )
//****************************************************************************
//
//  Purpose:
//
//    CLENSHAW_CURTIS_COMPUTE_POINTS computes Clenshaw Curtis quadrature points.
//
//  Discussion:
//
//    Our convention is that the abscissas are numbered from left to right.
//
//    This rule is defined on [-1,1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, Scalar X[N], the abscissas.
//
{
  int index;

  if (n<1) {
    std::cerr << "\n";
    std::cerr << "CLENSHAW_CURTIS_COMPUTE_POINTS - Fatal error!\n";
    std::cerr << "  N < 1.\n";
    std::exit(1);
  }
  else if (n==1) {
    x[0] = 0.0;
  }
  else {
    for (index=1;index<=n;index++) {
      x[index-1] = std::cos((Scalar)(n-index)*M_PI/(Scalar)(n-1));
    }
    x[0] = -1.0;
    if ((n%2)==1) {
      x[(n-1)/2] = 0.0;
    }
    x[n-1] = +1.0;
  }
  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::clenshaw_curtis_compute_weights ( int n, Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    CLENSHAW_CURTIS_COMPUTE_WEIGHTS computes Clenshaw Curtis quadrature weights.
//
//  Discussion:
//
//    The user must preallocate space for the output array W.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Charles Clenshaw, Alan Curtis,
//    A Method for Numerical Integration on an Automatic Computer,
//    Numerische Mathematik,
//    Volume 2, Number 1, December 1960, pages 197-205.
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, Scalar W[N], the weights.
//
{
  Scalar b, theta;
  int i, j;

  if (n<1) {
    std::cerr << "\n";
    std::cerr << "CLENSHAW_CURTIS_COMPUTE_WEIGHTS - Fatal error!\n";
    std::cerr << "  N < 1.\n";
    std::exit(1);
  }
  else if (n==1) {
    w[0] = 2.0;
    return;
  }

  for (i=1;i<=n;i++) {
    theta = (Scalar)(i-1)*M_PI/(Scalar)(n-1);

    w[i-1] = 1.0;

    for (j=1;j<=(n-1)/2;j++) {
      if (2*j==(n-1)) {
        b = 1.0;
      }
      else {
        b = 2.0;
      }

      w[i-1] = w[i-1]-b*std::cos(2.0*(Scalar)j*theta)/(Scalar)(4*j*j-1);
    }
  }

  w[0] = w[0]/(Scalar)(n-1);
  for (i=1;i<n-1;i++) {
    w[i] = 2.0*w[i]/(Scalar)(n-1);
  }
  w[n-1] = w[n-1]/(Scalar)(n-1);

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::fejer2_compute ( int n, Scalar x[], Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    FEJER2_COMPUTE computes a Fejer type 2 rule.
//
//  Discussion:
//
//    Our convention is that the abscissas are numbered from left to right.
//
//    The rule is defined on [-1,1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//    1 <= N.
//
//    Output, Scalar X[N], the abscissas.
//
//    Output, Scalar W[N], the weights.
//
{
  Scalar p, theta;

  if (n<1) {
    std::cerr << "\n";
    std::cerr << "FEJER2_COMPUTE - Fatal error!\n";
    std::cerr << "  Illegal value of N = " << n << "\n";
    std::exit(1);
  }
  else if (n==1) {
    x[0] = 0.0;
    w[0] = 2.0;
    return;
  }

  for (int i=0;i<n;i++) {
    x[i] =  std::cos((Scalar)(n-i)*M_PI/(Scalar)(n+1));
  }
  if ((n%2)==1) {
    x[(n-1)/2] = 0.0;
  }

  if (n==2) {
    w[0] = 1.0;
    w[1] = 1.0;
  }
  else {
    for (int i=0;i<n;i++) {
      theta = (Scalar)(n-i)*M_PI/(Scalar)(n+1);

      w[i] = 1.0;

      for (int j=1;j<=((n-1)/2);j++) {
        w[i] = w[i]-2.0*std::cos(2.0*(Scalar)j*theta)/(Scalar)(4*j*j-1);
      }
      p = 2.0*(Scalar)(((n+1)/2))-1.0;
      w[i] = w[i]-std::cos((p+1.0)*theta)/p;
    }
    for (int i=0;i<n;i++) {
      w[i] = 2.0*w[i]/(Scalar)(n+1);
    }
  }
  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::fejer2_compute_points ( int n, Scalar x[] )
//****************************************************************************
//
//  Purpose:
//
//    FEJER2_COMPUTE_POINTS computes Fejer type 2 quadrature points.
//
//  Discussion:
//
//    Our convention is that the abscissas are numbered from left to right.
//
//    The rule is defined on [-1,1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//    1 <= N.
//
//    Output, Scalar X[N], the abscissas.
//
{
  int i;

  if (n<1) {
    std::cerr << "\n";
    std::cerr << "FEJER2_COMPUTE_POINTS - Fatal error!\n";
    std::cerr << "  N < 1.\n";
    std::exit(1);
  }
  else if (n==1) {
    x[0] = 0.0;
  }
  else {
    for (i=1;i<=n;i++) {
      x[i-1] = std::cos((Scalar)(n+1-i)*M_PI/(Scalar)(n+1));
    }
    if ((n%2)==1) {
      x[(n-1)/2] = 0.0;
    }
  }
  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::fejer2_compute_weights ( int n, Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    FEJER2_COMPUTE_WEIGHTS computes Fejer type 2 quadrature weights.
//
//  Discussion:
//
//    The user must preallocate space for the output array W.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//    Walter Gautschi,
//    Numerical Quadrature in the Presence of a Singularity,
//    SIAM Journal on Numerical Analysis,
//    Volume 4, Number 3, 1967, pages 357-362.
//
//    Joerg Waldvogel,
//    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
//    BIT Numerical Mathematics,
//    Volume 43, Number 1, 2003, pages 1-18.
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, Scalar W[N], the weights.
//
{
  int i, j;
  Scalar p, theta;

  if (n<1) {
    std::cerr << "\n";
    std::cerr << "FEJER2_COMPUTE_WEIGHTS - Fatal error!\n";
    std::cerr << "  N < 1.\n";
    std::exit(1);
  }
  else if (n==1) {
    w[0] = 2.0;
  }
  else if (n==2) {
    w[0] = 1.0;
    w[1] = 1.0;
  }
  else {
    for (i=1;i<=n;i++) {
      theta = (Scalar)(n+1-i)*M_PI/(Scalar)(n+1);

      w[i-1] = 1.0;

      for (j=1;j<=((n-1)/2);j++) {
        w[i-1] = w[i-1]-2.0*std::cos(2.0*(Scalar)j*theta)/(Scalar)(4*j*j-1);
      }
      p = 2.0*(Scalar)(((n+1)/2))-1.0;
      w[i-1] = w[i-1]-std::cos((p+1.0)*theta)/p;
    }
    for (i=0;i<n;i++) {
      w[i] = 2.0*w[i]/(Scalar)(n+1);
    }
  }
  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::hermite_compute ( int n, Scalar x[], Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    HERMITE_COMPUTE computes a Gauss-Hermite quadrature rule.
//
//  Discussion:
//
//    The code uses an algorithm by Elhay and Kautsky.
//
//    The abscissas are the zeros of the N-th order Hermite polynomial.
//
//    The integral:
//
//      integral ( -oo < x < +oo ) exp ( - x * x ) * f(x) dx
//
//    The quadrature rule:
//
//      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2011
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int N, the number of abscissas.
//
//    Output, Scalar X[N], the abscissas.
//
//    Output, Scalar W[N], the weights.
//
{
  Scalar* bj;
//
//  Define the zero-th moment.
//
  Scalar zemu = std::sqrt(M_PI); 
//
//  Define the Jacobi matrix.
//
  bj = new Scalar[n];

  for (int i=0;i<n;i++) {
    bj[i] = std::sqrt((Scalar)(i+1)/2.0);
  }

  for (int i=0;i<n;i++) {
    x[i] = 0.0;
  }

  w[0] = std::sqrt (zemu);
  for (int i=1;i<n;i++) {
    w[i] = 0.0;
  }
//
//  Diagonalize the Jacobi matrix.
//
  IntrepidBurkardtRules::imtqlx ( n, x, bj, w );

  for (int i=0;i<n;i++) {
    w[i] = w[i]*w[i];
  }
  
  // Ensure that zero is actually zero.
  if (n%2) {
    int ind = (int)((Scalar)n/2.0);
    x[ind]  = 0.0;
  }

  delete [] bj;

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::hermite_compute_points ( int order, Scalar x[] )
//****************************************************************************
//
//  Purpose:
//
//    HERMITE_COMPUTE_POINTS computes Hermite quadrature points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order.
//
//    Output, Scalar X[ORDER], the abscissas.
//
{
  Scalar *w;  w = new Scalar[order];
  IntrepidBurkardtRules::hermite_compute ( order, x, w );
  delete [] w;

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::hermite_compute_weights ( int order, Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    HERMITE_COMPUTE_WEIGHTS computes Hermite quadrature weights.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order.
//
//    Output, Scalar W[ORDER], the weights.
//
{
  Scalar *x; x = new Scalar[order];
  IntrepidBurkardtRules::hermite_compute ( order, x, w );
  delete [] x;

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::hermite_genz_keister_lookup ( int n, Scalar x[], Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    HERMITE_GENZ_KEISTER_LOOKUP looks up a Genz-Keister Hermite rule.
//
//  Discussion:
//
//    The integral:
//
//      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
//
//    The quadrature rule:
//
//      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
//
//    A nested family of rules for the Hermite integration problem
//    was produced by Genz and Keister.  The structure of the nested
//    family was denoted by 1+2+6+10+16, that is, it comprised rules
//    of successive orders O = 1, 3, 9, 19, and 35.
//
//    The precisions of these rules are P = 1, 5, 15, 29, and 51.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz, Bradley Keister,
//    Fully symmetric interpolatory rules for multiple integrals
//    over infinite regions with Gaussian weight,
//    Journal of Computational and Applied Mathematics,
//    Volume 71, 1996, pages 299-309
//
//    Florian Heiss, Viktor Winschel,
//    Likelihood approximation by numerical integration on sparse grids,
//    Journal of Econometrics,
//    Volume 144, 2008, pages 62-80.
//
//    Thomas Patterson,
//    The Optimal Addition of Points to Quadrature Formulae,
//    Mathematics of Computation,
//    Volume 22, Number 104, October 1968, pages 847-856.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be 1, 3, 9, 19, or 35.
//
//    Output, Scalar X[N], the abscissas.
//
//    Output, Scalar W[N], the weights.
//
{
  IntrepidBurkardtRules::hermite_genz_keister_lookup_points ( n, x );
  IntrepidBurkardtRules::hermite_genz_keister_lookup_weights ( n, w );

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::hermite_genz_keister_lookup_points ( int n, Scalar x[] )
//****************************************************************************
//
//  Purpose:
//
//    HERMITE_GENZ_KEISTER_LOOKUP_POINTS looks up Genz-Keister Hermite abscissas.
//
//  Discussion:
//
//    The integral:
//
//      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
//
//    The quadrature rule:
//
//      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
//
//    A nested family of rules for the Hermite integration problem
//    was produced by Genz and Keister.  The structure of the nested
//    family was denoted by 1+2+6+10+16, that is, it comprised rules
//    of successive orders O = 1, 3, 9, 19, and 35.
//
//    The precisions of these rules are P = 1, 5, 15, 29, and 51.
//
//    Three related families begin the same way, but end with a different final
//    rule.  As a convenience, this function includes these final rules as well:
//
//    Designation  Orders       Precisions
//
//    1+2+6+10+16, 1,3,9,19,35  1,5,15,29,51
//    1+2+6+10+18  1,3,9,19,37  1,5,15,29,55
//    1+2+6+10+22  1,3,9,19,41  1,5,15,29,63
//    1+2+6+10+24  1,3,9,19,43  1,5,15,29,67
//
//    Some of the data in this function was kindly supplied directly by
//    Alan Genz on 24 April 2011.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz, Bradley Keister,
//    Fully symmetric interpolatory rules for multiple integrals
//    over infinite regions with Gaussian weight,
//    Journal of Computational and Applied Mathematics,
//    Volume 71, 1996, pages 299-309
//
//    Florian Heiss, Viktor Winschel,
//    Likelihood approximation by numerical integration on sparse grids,
//    Journal of Econometrics,
//    Volume 144, 2008, pages 62-80.
//
//    Thomas Patterson,
//    The Optimal Addition of Points to Quadrature Formulae,
//    Mathematics of Computation,
//    Volume 22, Number 104, October 1968, pages 847-856.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be 1, 3, 9, 19, 35, 27, 41, or 43.
//
//    Output, Scalar X[N], the abscissas.
//
{
  if (n==1) {
    x[ 0] =   0.0000000000000000E+00;
  }
  else if (n==3) {
    x[ 0] =  -1.2247448713915889E+00;
    x[ 1] =   0.0000000000000000E+00;
    x[ 2] =   1.2247448713915889E+00;
  }
  else if (n==9) {
    x[ 0] =  -2.9592107790638380E+00;
    x[ 1] =  -2.0232301911005157E+00;
    x[ 2] =  -1.2247448713915889E+00;
    x[ 3] =  -5.2403354748695763E-01;
    x[ 4] =   0.0000000000000000E+00;
    x[ 5] =   5.2403354748695763E-01;
    x[ 6] =   1.2247448713915889E+00;
    x[ 7] =   2.0232301911005157E+00;
    x[ 8] =   2.9592107790638380E+00;
  }
  else if (n==19) {
    x[ 0] =  -4.4995993983103881E+00;
    x[ 1] =  -3.6677742159463378E+00;
    x[ 2] =  -2.9592107790638380E+00;
    x[ 3] =  -2.2665132620567876E+00;
    x[ 4] =  -2.0232301911005157E+00;
    x[ 5] =  -1.8357079751751868E+00;
    x[ 6] =  -1.2247448713915889E+00;
    x[ 7] =  -8.7004089535290285E-01;
    x[ 8] =  -5.2403354748695763E-01;
    x[ 9] =   0.0000000000000000E+00;
    x[10] =   5.2403354748695763E-01;
    x[11] =   8.7004089535290285E-01;
    x[12] =   1.2247448713915889E+00;
    x[13] =   1.8357079751751868E+00;
    x[14] =   2.0232301911005157E+00;
    x[15] =   2.2665132620567876E+00;
    x[16] =   2.9592107790638380E+00;
    x[17] =   3.6677742159463378E+00;
    x[18] =   4.4995993983103881E+00;
  }
  else if (n==35) {
    x[ 0] =  -6.3759392709822356E+00;
    x[ 1] =  -5.6432578578857449E+00;
    x[ 2] =  -5.0360899444730940E+00;
    x[ 3] =  -4.4995993983103881E+00;
    x[ 4] =  -4.0292201405043713E+00;
    x[ 5] =  -3.6677742159463378E+00;
    x[ 6] =  -3.3491639537131945E+00;
    x[ 7] =  -2.9592107790638380E+00;
    x[ 8] =  -2.5705583765842968E+00;
    x[ 9] =  -2.2665132620567876E+00;
    x[10] =  -2.0232301911005157E+00;
    x[11] =  -1.8357079751751868E+00;
    x[12] =  -1.5794121348467671E+00;
    x[13] =  -1.2247448713915889E+00;
    x[14] =  -8.7004089535290285E-01;
    x[15] =  -5.2403354748695763E-01;
    x[16] =  -1.7606414208200893E-01;
    x[17] =   0.0000000000000000E+00;
    x[18] =   1.7606414208200893E-01;
    x[19] =   5.2403354748695763E-01;
    x[20] =   8.7004089535290285E-01;
    x[21] =   1.2247448713915889E+00;
    x[22] =   1.5794121348467671E+00;
    x[23] =   1.8357079751751868E+00;
    x[24] =   2.0232301911005157E+00;
    x[25] =   2.2665132620567876E+00;
    x[26] =   2.5705583765842968E+00;
    x[27] =   2.9592107790638380E+00;
    x[28] =   3.3491639537131945E+00;
    x[29] =   3.6677742159463378E+00;
    x[30] =   4.0292201405043713E+00;
    x[31] =   4.4995993983103881E+00;
    x[32] =   5.0360899444730940E+00;
    x[33] =   5.6432578578857449E+00;
    x[34] =   6.3759392709822356E+00;
  }
  else if (n==37) {
    x[ 0] =  -6.853200069757519;
    x[ 1] =  -6.124527854622158;
    x[ 2] =  -5.521865209868350;
    x[ 3] =  -4.986551454150765;
    x[ 4] =  -4.499599398310388;
    x[ 5] =  -4.057956316089741;
    x[ 6] =  -3.667774215946338;
    x[ 7] =  -3.315584617593290;
    x[ 8] =  -2.959210779063838;
    x[ 9] =  -2.597288631188366;
    x[10] =  -2.266513262056788;
    x[11] =  -2.023230191100516;
    x[12] =  -1.835707975175187;
    x[13] =  -1.561553427651873;
    x[14] =  -1.224744871391589;
    x[15] =  -0.870040895352903;
    x[16] =  -0.524033547486958;
    x[17] =  -0.214618180588171;
    x[18] =   0.000000000000000;
    x[19] =   0.214618180588171;
    x[20] =   0.524033547486958;
    x[21] =   0.870040895352903;
    x[22] =   1.224744871391589;
    x[23] =   1.561553427651873;
    x[24] =   1.835707975175187;
    x[25] =   2.023230191100516;
    x[26] =   2.266513262056788;
    x[27] =   2.597288631188366;
    x[28] =   2.959210779063838;
    x[29] =   3.315584617593290;
    x[30] =   3.667774215946338;
    x[31] =   4.057956316089741;
    x[32] =   4.499599398310388;
    x[33] =   4.986551454150765;
    x[34] =   5.521865209868350;
    x[35] =   6.124527854622158;
    x[36] =   6.853200069757519;
  }
  else if (n==41) {
    x[ 0] =  -7.251792998192644;
    x[ 1] =  -6.547083258397540;
    x[ 2] =  -5.961461043404500;
    x[ 3] =  -5.437443360177798;
    x[ 4] =  -4.953574342912980;
    x[ 5] =  -4.4995993983103881;
    x[ 6] =  -4.070919267883068;
    x[ 7] =  -3.6677742159463378;
    x[ 8] =  -3.296114596212218;
    x[ 9] =  -2.9592107790638380;
    x[10] =  -2.630415236459871;
    x[11] =  -2.2665132620567876;
    x[12] =  -2.043834754429505;
    x[13] =  -2.0232301911005157;
    x[14] =  -1.8357079751751868;
    x[15] =  -1.585873011819188;
    x[16] =  -1.2247448713915889;
    x[17] =  -0.87004089535290285;
    x[18] =  -0.52403354748695763;
    x[19] =  -0.195324784415805;
    x[20] =   0.0000000000000000;
    x[21] =   0.195324784415805;
    x[22] =   0.52403354748695763;
    x[23] =   0.87004089535290285;
    x[24] =   1.2247448713915889;
    x[25] =   1.585873011819188;
    x[26] =   1.8357079751751868;
    x[27] =   2.0232301911005157;
    x[28] =   2.043834754429505;
    x[29] =   2.2665132620567876;
    x[30] =   2.630415236459871;
    x[31] =   2.9592107790638380;
    x[32] =   3.296114596212218;
    x[33] =   3.6677742159463378;
    x[34] =   4.070919267883068;
    x[35] =   4.4995993983103881;
    x[36] =   4.953574342912980;
    x[37] =   5.437443360177798;
    x[38] =   5.961461043404500;
    x[39] =   6.547083258397540;
    x[40] =   7.251792998192644;
  }
  else if (n==43) {
    x[ 0] = -10.167574994881873;
    x[ 1] =  -7.231746029072501;
    x[ 2] =  -6.535398426382995;
    x[ 3] =  -5.954781975039809;
    x[ 4] =  -5.434053000365068;
    x[ 5] =  -4.952329763008589;
    x[ 6] =  -4.4995993983103881;
    x[ 7] =  -4.071335874253583;
    x[ 8] =  -3.6677742159463378;
    x[ 9] =  -3.295265921534226;
    x[10] =  -2.9592107790638380;
    x[11] =  -2.633356763661946;
    x[12] =  -2.2665132620567876;
    x[13] =  -2.089340389294661;
    x[14] =  -2.0232301911005157;
    x[15] =  -1.8357079751751868;
    x[16] =  -1.583643465293944;
    x[17] =  -1.2247448713915889;
    x[18] =  -0.87004089535290285;
    x[19] =  -0.52403354748695763;
    x[20] =  -0.196029453662011;
    x[21] =   0.0000000000000000;
    x[22] =   0.196029453662011;
    x[23] =   0.52403354748695763;
    x[24] =   0.87004089535290285;
    x[25] =   1.2247448713915889;
    x[26] =   1.583643465293944;
    x[27] =   1.8357079751751868;
    x[28] =   2.0232301911005157;
    x[29] =   2.089340389294661;
    x[30] =   2.2665132620567876;
    x[31] =   2.633356763661946;
    x[32] =   2.9592107790638380;
    x[33] =   3.295265921534226;
    x[34] =   3.6677742159463378;
    x[35] =   4.071335874253583;
    x[36] =   4.4995993983103881;
    x[37] =   4.952329763008589;
    x[38] =   5.434053000365068;
    x[39] =   5.954781975039809;
    x[40] =   6.535398426382995;
    x[41] =   7.231746029072501;
    x[42] =  10.167574994881873;
  }
  else {
    std::cerr << "\n";
    std::cerr << "HERMITE_GENZ_KEISTER_LOOKUP_POINTS - Fatal error!\n";
    std::cerr << "  Illegal input value of N.\n";
    std::cerr << "  N must be 1, 3, 9, 19, 35, 37, 41 or 43.\n";
    std::exit(1);
  }
  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::hermite_genz_keister_lookup_weights ( int n, Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    HERMITE_GENZ_KEISTER_LOOKUP_WEIGHTS looks up Genz-Keister Hermite weights.
//
//  Discussion:
//
//    The integral:
//
//      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
//
//    The quadrature rule:
//
//      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
//
//    A nested family of rules for the Hermite integration problem
//    was produced by Genz and Keister.  The structure of the nested
//    family was denoted by 1+2+6+10+16, that is, it comprised rules
//    of successive orders O = 1, 3, 9, 19, and 35.
//
//    The precisions of these rules are P = 1, 5, 15, 29, and 51.
//
//    Three related families begin the same way, but end with a different final
//    rule.  As a convenience, this function includes these final rules as well:
//
//    Designation  Orders       Precisions
//
//    1+2+6+10+16, 1,3,9,19,35  1,5,15,29,51
//    1+2+6+10+18  1,3,9,19,37  1,5,15,29,55
//    1+2+6+10+22  1,3,9,19,41  1,5,15,29,63
//    1+2+6+10+24  1,3,9,19,43  1,5,15,29,67
//
//    Some of the data in this function was kindly supplied directly by
//    Alan Genz on 24 April 2011.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz, Bradley Keister,
//    Fully symmetric interpolatory rules for multiple integrals
//    over infinite regions with Gaussian weight,
//    Journal of Computational and Applied Mathematics,
//    Volume 71, 1996, pages 299-309
//
//    Florian Heiss, Viktor Winschel,
//    Likelihood approximation by numerical integration on sparse grids,
//    Journal of Econometrics,
//    Volume 144, 2008, pages 62-80.
//
//    Thomas Patterson,
//    The Optimal Addition of Points to Quadrature Formulae,
//    Mathematics of Computation,
//    Volume 22, Number 104, October 1968, pages 847-856.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be 1, 3, 9, 19, 35, 37, 41, or 43.
//
//    Output, Scalar W[N], the weights.
//
{
  if (n==1) {
    w[ 0] =   1.7724538509055159E+00;
  }
  else if (n==3) {
    w[ 0] =   2.9540897515091930E-01;
    w[ 1] =   1.1816359006036772E+00;
    w[ 2] =   2.9540897515091930E-01;
  }
  else if (n==9) {
    w[ 0] =   1.6708826306882348E-04;
    w[ 1] =   1.4173117873979098E-02;
    w[ 2] =   1.6811892894767771E-01;
    w[ 3] =   4.7869428549114124E-01;
    w[ 4] =   4.5014700975378197E-01;
    w[ 5] =   4.7869428549114124E-01;
    w[ 6] =   1.6811892894767771E-01;
    w[ 7] =   1.4173117873979098E-02;
    w[ 8] =   1.6708826306882348E-04;
  }
  else if (n==19) {
    w[ 0] =   1.5295717705322357E-09;
    w[ 1] =   1.0802767206624762E-06;
    w[ 2] =   1.0656589772852267E-04;
    w[ 3] =   5.1133174390883855E-03;
    w[ 4] =  -1.1232438489069229E-02;
    w[ 5] =   3.2055243099445879E-02;
    w[ 6] =   1.1360729895748269E-01;
    w[ 7] =   1.0838861955003017E-01;
    w[ 8] =   3.6924643368920851E-01;
    w[ 9] =   5.3788160700510168E-01;
    w[10] =   3.6924643368920851E-01;
    w[11] =   1.0838861955003017E-01;
    w[12] =   1.1360729895748269E-01;
    w[13] =   3.2055243099445879E-02;
    w[14] =  -1.1232438489069229E-02;
    w[15] =   5.1133174390883855E-03;
    w[16] =   1.0656589772852267E-04;
    w[17] =   1.0802767206624762E-06;
    w[18] =   1.5295717705322357E-09;
  }
  else if (n==35) {
    w[ 0] =   1.8684014894510604E-18;
    w[ 1] =   9.6599466278563243E-15;
    w[ 2] =   5.4896836948499462E-12;
    w[ 3] =   8.1553721816916897E-10;
    w[ 4] =   3.7920222392319532E-08;
    w[ 5] =   4.3737818040926989E-07;
    w[ 6] =   4.8462799737020461E-06;
    w[ 7] =   6.3328620805617891E-05;
    w[ 8] =   4.8785399304443770E-04;
    w[ 9] =   1.4515580425155904E-03;
    w[10] =   4.0967527720344047E-03;
    w[11] =   5.5928828911469180E-03;
    w[12] =   2.7780508908535097E-02;
    w[13] =   8.0245518147390893E-02;
    w[14] =   1.6371221555735804E-01;
    w[15] =   2.6244871488784277E-01;
    w[16] =   3.3988595585585218E-01;
    w[17] =   9.1262675363737921E-04;
    w[18] =   3.3988595585585218E-01;
    w[19] =   2.6244871488784277E-01;
    w[20] =   1.6371221555735804E-01;
    w[21] =   8.0245518147390893E-02;
    w[22] =   2.7780508908535097E-02;
    w[23] =   5.5928828911469180E-03;
    w[24] =   4.0967527720344047E-03;
    w[25] =   1.4515580425155904E-03;
    w[26] =   4.8785399304443770E-04;
    w[27] =   6.3328620805617891E-05;
    w[28] =   4.8462799737020461E-06;
    w[29] =   4.3737818040926989E-07;
    w[30] =   3.7920222392319532E-08;
    w[31] =   8.1553721816916897E-10;
    w[32] =   5.4896836948499462E-12;
    w[33] =   9.6599466278563243E-15;
    w[34] =   1.8684014894510604E-18;
  }
  else if (n==37) {
    w[ 0] = 0.19030350940130498E-20;
    w[ 1] = 0.187781893143728947E-16;
    w[ 2] = 0.182242751549129356E-13;
    w[ 3] = 0.45661763676186859E-11;
    w[ 4] = 0.422525843963111041E-09;
    w[ 5] = 0.16595448809389819E-07;
    w[ 6] = 0.295907520230744049E-06;
    w[ 7] = 0.330975870979203419E-05;
    w[ 8] = 0.32265185983739747E-04;
    w[ 9] = 0.234940366465975222E-03;
    w[10] = 0.985827582996483824E-03;
    w[11] = 0.176802225818295443E-02;
    w[12] = 0.43334988122723492E-02;
    w[13] = 0.15513109874859354E-01;
    w[14] = 0.442116442189845444E-01;
    w[15] = 0.937208280655245902E-01;
    w[16] = 0.143099302896833389E+00;
    w[17] = 0.147655710402686249E+00;
    w[18] = 0.968824552928425499E-01;
    w[19] = 0.147655710402686249E+00;
    w[20] = 0.143099302896833389E+00;
    w[21] = 0.937208280655245902E-01;
    w[22] = 0.442116442189845444E-01;
    w[23] = 0.15513109874859354E-01;
    w[24] = 0.43334988122723492E-02;
    w[25] = 0.176802225818295443E-02;
    w[26] = 0.985827582996483824E-03;
    w[27] = 0.234940366465975222E-03;
    w[28] = 0.32265185983739747E-04;
    w[29] = 0.330975870979203419E-05;
    w[30] = 0.295907520230744049E-06;
    w[31] = 0.16595448809389819E-07;
    w[32] = 0.422525843963111041E-09;
    w[33] = 0.45661763676186859E-11;
    w[34] = 0.182242751549129356E-13;
    w[35] = 0.187781893143728947E-16;
    w[36] = 0.19030350940130498E-20;
  }
  else if (n==41) {
    w[ 0] =   0.664195893812757801E-23;
    w[ 1] =   0.860427172512207236E-19;
    w[ 2] =   0.1140700785308509E-15;
    w[ 3] =   0.408820161202505983E-13;
    w[ 4] =   0.581803393170320419E-11;
    w[ 5] =   0.400784141604834759E-09;
    w[ 6] =   0.149158210417831408E-07;
    w[ 7] =   0.315372265852264871E-06;
    w[ 8] =   0.381182791749177506E-05;
    w[ 9] =   0.288976780274478689E-04;
    w[10] =   0.189010909805097887E-03;
    w[11] =   0.140697424065246825E-02;
    w[12] = - 0.144528422206988237E-01;
    w[13] =   0.178852543033699732E-01;
    w[14] =   0.705471110122962612E-03;
    w[15] =   0.165445526705860772E-01;
    w[16] =   0.45109010335859128E-01;
    w[17] =   0.928338228510111845E-01;
    w[18] =   0.145966293895926429E+00;
    w[19] =   0.165639740400529554E+00;
    w[20] =   0.562793426043218877E-01;
    w[21] =   0.165639740400529554E+00;
    w[22] =   0.145966293895926429E+00;
    w[23] =   0.928338228510111845E-01;
    w[24] =   0.45109010335859128E-01;
    w[25] =   0.165445526705860772E-01;
    w[26] =   0.705471110122962612E-03;
    w[27] =   0.178852543033699732E-01;
    w[28] = - 0.144528422206988237E-01;
    w[29] =   0.140697424065246825E-02;
    w[30] =   0.189010909805097887E-03;
    w[31] =   0.288976780274478689E-04;
    w[32] =   0.381182791749177506E-05;
    w[33] =   0.315372265852264871E-06;
    w[34] =   0.149158210417831408E-07;
    w[35] =   0.400784141604834759E-09;
    w[36] =   0.581803393170320419E-11;
    w[37] =   0.408820161202505983E-13;
    w[38] =   0.1140700785308509E-15;
    w[39] =   0.860427172512207236E-19;
    w[40] =   0.664195893812757801E-23;
  }
  else if (n==43) {
    w[ 0] =   0.546191947478318097E-37;
    w[ 1] =   0.87544909871323873E-23;
    w[ 2] =   0.992619971560149097E-19;
    w[ 3] =   0.122619614947864357E-15;
    w[ 4] =   0.421921851448196032E-13;
    w[ 5] =   0.586915885251734856E-11;
    w[ 6] =   0.400030575425776948E-09;
    w[ 7] =   0.148653643571796457E-07;
    w[ 8] =   0.316018363221289247E-06;
    w[ 9] =   0.383880761947398577E-05;
    w[10] =   0.286802318064777813E-04;
    w[11] =   0.184789465688357423E-03;
    w[12] =   0.150909333211638847E-02;
    w[13] = - 0.38799558623877157E-02;
    w[14] =   0.67354758901013295E-02;
    w[15] =   0.139966252291568061E-02;
    w[16] =   0.163616873493832402E-01;
    w[17] =   0.450612329041864976E-01;
    w[18] =   0.928711584442575456E-01;
    w[19] =   0.145863292632147353E+00;
    w[20] =   0.164880913687436689E+00;
    w[21] =   0.579595986101181095E-01;
    w[22] =   0.164880913687436689E+00;
    w[23] =   0.145863292632147353E+00;
    w[24] =   0.928711584442575456E-01;
    w[25] =   0.450612329041864976E-01;
    w[26] =   0.163616873493832402E-01;
    w[27] =   0.139966252291568061E-02;
    w[28] =   0.67354758901013295E-02;
    w[29] = - 0.38799558623877157E-02;
    w[30] =   0.150909333211638847E-02;
    w[31] =   0.184789465688357423E-03;
    w[32] =   0.286802318064777813E-04;
    w[33] =   0.383880761947398577E-05;
    w[34] =   0.316018363221289247E-06;
    w[35] =   0.148653643571796457E-07;
    w[36] =   0.400030575425776948E-09;
    w[37] =   0.586915885251734856E-11;
    w[38] =   0.421921851448196032E-13;
    w[39] =   0.122619614947864357E-15;
    w[40] =   0.992619971560149097E-19;
    w[41] =   0.87544909871323873E-23;
    w[42] =   0.546191947478318097E-37;
  }
  else {
    std::cerr << "\n";
    std::cerr << "HERMITE_GENZ_KEISTER_LOOKUP_WEIGHTS - Fatal error!\n";
    std::cerr << "  Illegal input value of N.\n";
    std::cerr << "  N must be 1, 3, 9, 19, 35, 37, 41 or 43.\n";
    std::exit(1);
  }
  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::hermite_lookup ( int n, Scalar x[], Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    HERMITE_LOOKUP looks up abscissas and weights for Gauss-Hermite quadrature.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Vladimir Krylov,
//    Approximate Calculation of Integrals,
//    Dover, 2006,
//    ISBN: 0486445798.
//    LC: QA311.K713.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3,
//    LC: QA47.M315.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be between 1 and 20.
//
//    Output, Scalar X[N], the abscissas.
//
//    Output, Scalar W[N], the weights.
//
{
  IntrepidBurkardtRules::hermite_lookup_points ( n, x );
  IntrepidBurkardtRules::hermite_lookup_weights ( n, w );

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::hermite_lookup_points ( int n, Scalar x[] )
//****************************************************************************
//
//  Purpose:
//
//    HERMITE_LOOKUP_POINTS looks up abscissas for Hermite quadrature.
//
//  Discussion:
//
//    The integral:
//
//      integral ( -oo < x < +oo ) exp ( - x * x ) * f(x) dx
//
//    The quadrature rule:
//
//      sum ( 1 <= i <= n ) w(i) * f ( x(i) ).
//
//    Mathematica can numerically estimate the abscissas
//    of order N to P digits by the command:
//
//      NSolve [ HermiteH [ n, x ] == 0, x, p ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Vladimir Krylov,
//    Approximate Calculation of Integrals,
//    Dover, 2006,
//    ISBN: 0486445798,
//    LC: QA311.K713.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3,
//    LC: QA47.M315.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be between 1 and 20.
//
//    Output, Scalar X[N], the abscissas.
//
{
  if (n==1) {
    x[ 0] = 0.0;
  }
  else if(n==2) {
    x[ 0] = - 0.707106781186547524400844362105E+00;
    x[ 1] =   0.707106781186547524400844362105E+00;
  }
  else if (n==3) {
    x[ 0] = - 0.122474487139158904909864203735E+01;
    x[ 1] =   0.0E+00;
    x[ 2] =   0.122474487139158904909864203735E+01;
  }
  else if (n==4) {
    x[ 0] = - 0.165068012388578455588334111112E+01;
    x[ 1] = - 0.524647623275290317884060253835E+00;
    x[ 2] =   0.524647623275290317884060253835E+00;
    x[ 3] =   0.165068012388578455588334111112E+01;
  }
  else if (n==5) {
    x[ 0] = - 0.202018287045608563292872408814E+01;
    x[ 1] = - 0.958572464613818507112770593893E+00;
    x[ 2] =   0.0E+00;
    x[ 3] =   0.958572464613818507112770593893E+00;
    x[ 4] =   0.202018287045608563292872408814E+01;
  }
  else if (n==6) {
    x[ 0] = - 0.235060497367449222283392198706E+01;
    x[ 1] = - 0.133584907401369694971489528297E+01;
    x[ 2] = - 0.436077411927616508679215948251E+00;
    x[ 3] =   0.436077411927616508679215948251E+00;
    x[ 4] =   0.133584907401369694971489528297E+01;
    x[ 5] =   0.235060497367449222283392198706E+01;
  }
  else if (n==7) {
    x[ 0] = - 0.265196135683523349244708200652E+01;
    x[ 1] = - 0.167355162876747144503180139830E+01;
    x[ 2] = - 0.816287882858964663038710959027E+00;
    x[ 3] =   0.0E+00;
    x[ 4] =   0.816287882858964663038710959027E+00;
    x[ 5] =   0.167355162876747144503180139830E+01;
    x[ 6] =   0.265196135683523349244708200652E+01;
  }
  else if (n==8) {
    x[ 0] = - 0.293063742025724401922350270524E+01;
    x[ 1] = - 0.198165675669584292585463063977E+01;
    x[ 2] = - 0.115719371244678019472076577906E+01;
    x[ 3] = - 0.381186990207322116854718885584E+00;
    x[ 4] =   0.381186990207322116854718885584E+00;
    x[ 5] =   0.115719371244678019472076577906E+01;
    x[ 6] =   0.198165675669584292585463063977E+01;
    x[ 7] =   0.293063742025724401922350270524E+01;
  }
  else if (n==9) {
    x[ 0] = - 0.319099320178152760723004779538E+01;
    x[ 1] = - 0.226658058453184311180209693284E+01;
    x[ 2] = - 0.146855328921666793166701573925E+01;
    x[ 3] = - 0.723551018752837573322639864579E+00;
    x[ 4] =   0.0E+00;
    x[ 5] =   0.723551018752837573322639864579E+00;
    x[ 6] =   0.146855328921666793166701573925E+01;
    x[ 7] =   0.226658058453184311180209693284E+01;
    x[ 8] =   0.319099320178152760723004779538E+01;
  }
  else if (n==10) {
    x[ 0] =  - 0.343615911883773760332672549432E+01;
    x[ 1] =  - 0.253273167423278979640896079775E+01;
    x[ 2] =  - 0.175668364929988177345140122011E+01;
    x[ 3] =  - 0.103661082978951365417749191676E+01;
    x[ 4] =  - 0.342901327223704608789165025557E+00;
    x[ 5] =    0.342901327223704608789165025557E+00;
    x[ 6] =    0.103661082978951365417749191676E+01;
    x[ 7] =    0.175668364929988177345140122011E+01;
    x[ 8] =    0.253273167423278979640896079775E+01;
    x[ 9] =    0.343615911883773760332672549432E+01;
  }
  else if (n==11) {
    x[ 0] =  - 0.366847084655958251845837146485E+01;
    x[ 1] =  - 0.278329009978165177083671870152E+01;
    x[ 2] =  - 0.202594801582575533516591283121E+01;
    x[ 3] =  - 0.132655708449493285594973473558E+01;
    x[ 4] =  - 0.656809566882099765024611575383E+00;
    x[ 5] =    0.0E+00;
    x[ 6] =    0.656809566882099765024611575383E+00;
    x[ 7] =    0.132655708449493285594973473558E+01;
    x[ 8] =    0.202594801582575533516591283121E+01;
    x[ 9] =    0.278329009978165177083671870152E+01;
    x[10] =    0.366847084655958251845837146485E+01;
  }
  else if (n==12) {
    x[ 0] =  - 0.388972489786978191927164274724E+01;
    x[ 1] =  - 0.302063702512088977171067937518E+01;
    x[ 2] =  - 0.227950708050105990018772856942E+01;
    x[ 3] =  - 0.159768263515260479670966277090E+01;
    x[ 4] =  - 0.947788391240163743704578131060E+00;
    x[ 5] =  - 0.314240376254359111276611634095E+00;
    x[ 6] =    0.314240376254359111276611634095E+00;
    x[ 7] =    0.947788391240163743704578131060E+00;
    x[ 8] =    0.159768263515260479670966277090E+01;
    x[ 9] =    0.227950708050105990018772856942E+01;
    x[10] =    0.302063702512088977171067937518E+01;
    x[11] =    0.388972489786978191927164274724E+01;
  }
  else if (n==13) {
    x[ 0] =  - 0.410133759617863964117891508007E+01;
    x[ 1] =  - 0.324660897837240998812205115236E+01;
    x[ 2] =  - 0.251973568567823788343040913628E+01;
    x[ 3] =  - 0.185310765160151214200350644316E+01;
    x[ 4] =  - 0.122005503659074842622205526637E+01;
    x[ 5] =  - 0.605763879171060113080537108602E+00;
    x[ 6] =    0.0E+00;
    x[ 7] =    0.605763879171060113080537108602E+00;
    x[ 8] =    0.122005503659074842622205526637E+01;
    x[ 9] =    0.185310765160151214200350644316E+01;
    x[10] =    0.251973568567823788343040913628E+01;
    x[11] =    0.324660897837240998812205115236E+01;
    x[12] =    0.410133759617863964117891508007E+01;
  }
  else if (n==14) {
    x[ 0] =  - 0.430444857047363181262129810037E+01;
    x[ 1] =  - 0.346265693360227055020891736115E+01;
    x[ 2] =  - 0.274847072498540256862499852415E+01;
    x[ 3] =  - 0.209518325850771681573497272630E+01;
    x[ 4] =  - 0.147668273114114087058350654421E+01;
    x[ 5] =  - 0.878713787329399416114679311861E+00;
    x[ 6] =  - 0.291745510672562078446113075799E+00;
    x[ 7] =    0.291745510672562078446113075799E+00;
    x[ 8] =    0.878713787329399416114679311861E+00;
    x[ 9] =    0.147668273114114087058350654421E+01;
    x[10] =    0.209518325850771681573497272630E+01;
    x[11] =    0.274847072498540256862499852415E+01;
    x[12] =    0.346265693360227055020891736115E+01;
    x[13] =    0.430444857047363181262129810037E+01;
  }
  else if (n==15) {
    x[ 0] =  - 0.449999070730939155366438053053E+01;
    x[ 1] =  - 0.366995037340445253472922383312E+01;
    x[ 2] =  - 0.296716692790560324848896036355E+01;
    x[ 3] =  - 0.232573248617385774545404479449E+01;
    x[ 4] =  - 0.171999257518648893241583152515E+01;
    x[ 5] =  - 0.113611558521092066631913490556E+01;
    x[ 6] =  - 0.565069583255575748526020337198E+00;
    x[ 7] =    0.0E+00;
    x[ 8] =    0.565069583255575748526020337198E+00;
    x[ 9] =    0.113611558521092066631913490556E+01;
    x[10] =    0.171999257518648893241583152515E+01;
    x[11] =    0.232573248617385774545404479449E+01;
    x[12] =    0.296716692790560324848896036355E+01;
    x[13] =    0.366995037340445253472922383312E+01;
    x[14] =    0.449999070730939155366438053053E+01;
  }
  else if (n==16) {
    x[ 0] =  - 0.468873893930581836468849864875E+01;
    x[ 1] =  - 0.386944790486012269871942409801E+01;
    x[ 2] =  - 0.317699916197995602681399455926E+01;
    x[ 3] =  - 0.254620215784748136215932870545E+01;
    x[ 4] =  - 0.195178799091625397743465541496E+01;
    x[ 5] =  - 0.138025853919888079637208966969E+01;
    x[ 6] =  - 0.822951449144655892582454496734E+00;
    x[ 7] =  - 0.273481046138152452158280401965E+00;
    x[ 8] =    0.273481046138152452158280401965E+00;
    x[ 9] =    0.822951449144655892582454496734E+00;
    x[10] =    0.138025853919888079637208966969E+01;
    x[11] =    0.195178799091625397743465541496E+01;
    x[12] =    0.254620215784748136215932870545E+01;
    x[13] =    0.317699916197995602681399455926E+01;
    x[14] =    0.386944790486012269871942409801E+01;
    x[15] =    0.468873893930581836468849864875E+01;
  }
  else if (n==17) {
    x[ 0] =  - 0.487134519367440308834927655662E+01;
    x[ 1] =  - 0.406194667587547430689245559698E+01;
    x[ 2] =  - 0.337893209114149408338327069289E+01;
    x[ 3] =  - 0.275776291570388873092640349574E+01;
    x[ 4] =  - 0.217350282666662081927537907149E+01;
    x[ 5] =  - 0.161292431422123133311288254454E+01;
    x[ 6] =  - 0.106764872574345055363045773799E+01;
    x[ 7] =  - 0.531633001342654731349086553718E+00;
    x[ 8] =    0.0E+00;
    x[ 9] =    0.531633001342654731349086553718E+00;
    x[10] =    0.106764872574345055363045773799E+01;
    x[11] =    0.161292431422123133311288254454E+01;
    x[12] =    0.217350282666662081927537907149E+01;
    x[13] =    0.275776291570388873092640349574E+01;
    x[14] =    0.337893209114149408338327069289E+01;
    x[15] =    0.406194667587547430689245559698E+01;
    x[16] =    0.487134519367440308834927655662E+01;
  }
  else if (n==18) {
    x[ 0] =  - 0.504836400887446676837203757885E+01;
    x[ 1] =  - 0.424811787356812646302342016090E+01;
    x[ 2] =  - 0.357376906848626607950067599377E+01;
    x[ 3] =  - 0.296137750553160684477863254906E+01;
    x[ 4] =  - 0.238629908916668600026459301424E+01;
    x[ 5] =  - 0.183553160426162889225383944409E+01;
    x[ 6] =  - 0.130092085838961736566626555439E+01;
    x[ 7] =  - 0.776682919267411661316659462284E+00;
    x[ 8] =  - 0.258267750519096759258116098711E+00;
    x[ 9] =    0.258267750519096759258116098711E+00;
    x[10] =    0.776682919267411661316659462284E+00;
    x[11] =    0.130092085838961736566626555439E+01;
    x[12] =    0.183553160426162889225383944409E+01;
    x[13] =    0.238629908916668600026459301424E+01;
    x[14] =    0.296137750553160684477863254906E+01;
    x[15] =    0.357376906848626607950067599377E+01;
    x[16] =    0.424811787356812646302342016090E+01;
    x[17] =    0.504836400887446676837203757885E+01;
  }
  else if (n==19) {
    x[ 0] =  - 0.522027169053748216460967142500E+01;
    x[ 1] =  - 0.442853280660377943723498532226E+01;
    x[ 2] =  - 0.376218735196402009751489394104E+01;
    x[ 3] =  - 0.315784881834760228184318034120E+01;
    x[ 4] =  - 0.259113378979454256492128084112E+01;
    x[ 5] =  - 0.204923170985061937575050838669E+01;
    x[ 6] =  - 0.152417061939353303183354859367E+01;
    x[ 7] =  - 0.101036838713431135136859873726E+01;
    x[ 8] =  - 0.503520163423888209373811765050E+00;
    x[ 9] =    0.0E+00;
    x[10] =    0.503520163423888209373811765050E+00;
    x[11] =    0.101036838713431135136859873726E+01;
    x[12] =    0.152417061939353303183354859367E+01;
    x[13] =    0.204923170985061937575050838669E+01;
    x[14] =    0.259113378979454256492128084112E+01;
    x[15] =    0.315784881834760228184318034120E+01;
    x[16] =    0.376218735196402009751489394104E+01;
    x[17] =    0.442853280660377943723498532226E+01;
    x[18] =    0.522027169053748216460967142500E+01;
  }
  else if (n==20) {
    x[ 0] =  - 0.538748089001123286201690041068E+01;
    x[ 1] =  - 0.460368244955074427307767524898E+01;
    x[ 2] =  - 0.394476404011562521037562880052E+01;
    x[ 3] =  - 0.334785456738321632691492452300E+01;
    x[ 4] =  - 0.278880605842813048052503375640E+01;
    x[ 5] =  - 0.225497400208927552308233334473E+01;
    x[ 6] =  - 0.173853771211658620678086566214E+01;
    x[ 7] =  - 0.123407621539532300788581834696E+01;
    x[ 8] =  - 0.737473728545394358705605144252E+00;
    x[ 9] =  - 0.245340708300901249903836530634E+00;
    x[10] =    0.245340708300901249903836530634E+00;
    x[11] =    0.737473728545394358705605144252E+00;
    x[12] =    0.123407621539532300788581834696E+01;
    x[13] =    0.173853771211658620678086566214E+01;
    x[14] =    0.225497400208927552308233334473E+01;
    x[15] =    0.278880605842813048052503375640E+01;
    x[16] =    0.334785456738321632691492452300E+01;
    x[17] =    0.394476404011562521037562880052E+01;
    x[18] =    0.460368244955074427307767524898E+01;
    x[19] =    0.538748089001123286201690041068E+01;
  }
  else {
    std::cerr << "\n";
    std::cerr << "HERMITE_LOOKUP_POINTS - Fatal error!\n";
    std::cerr << "  Illegal value of N = " << n << "\n";
    std::cerr << "  Legal values are 1 through 20.\n";
    std::exit(1);
  }

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::hermite_lookup_weights ( int n, Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    HERMITE_LOOKUP_WEIGHTS looks up weights for Hermite quadrature.
//
//  Discussion:
//
//    The integral:
//
//      integral ( -oo < x < +oo ) exp ( - x * x ) * f(x) dx
//
//    The quadrature rule:
//
//      sum ( 1 <= i <= n ) w(i) * f ( x(i) ).
//
//    Mathematica can numerically estimate the abscissas
//    of order N to P digits by the command:
//
//      NSolve [ HermiteH [ n, x ] == 0, x, p ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Vladimir Krylov,
//    Approximate Calculation of Integrals,
//    Dover, 2006,
//    ISBN: 0486445798,
//    LC: QA311.K713.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3,
//    LC: QA47.M315.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be between 1 and 20.
//
//    Output, Scalar W[N], the weights.
//
{
  if (n==1) {
    w[ 0] = 1.77245385090551602729816748334;
  }
  else if (n==2) {
    w[ 0] = 0.886226925452758013649083741671E+00;
    w[ 1] = 0.886226925452758013649083741671E+00;
  }
  else if (n==3) {
    w[ 0] = 0.295408975150919337883027913890E+00;
    w[ 1] = 0.118163590060367735153211165556E+01;
    w[ 2] = 0.295408975150919337883027913890E+00;
  }
  else if (n==4) {
    w[ 0] = 0.813128354472451771430345571899E-01;
    w[ 1] = 0.804914090005512836506049184481E+00;
    w[ 2] = 0.804914090005512836506049184481E+00;
    w[ 3] = 0.813128354472451771430345571899E-01;
  }
  else if (n==5) {
    w[ 0] = 0.199532420590459132077434585942E-01;
    w[ 1] = 0.393619323152241159828495620852E+00;
    w[ 2] = 0.945308720482941881225689324449E+00;
    w[ 3] = 0.393619323152241159828495620852E+00;
    w[ 4] = 0.199532420590459132077434585942E-01;
  }
  else if (n==6) {
    w[ 0] = 0.453000990550884564085747256463E-02;
    w[ 1] = 0.157067320322856643916311563508E+00;
    w[ 2] = 0.724629595224392524091914705598E+00;
    w[ 3] = 0.724629595224392524091914705598E+00;
    w[ 4] = 0.157067320322856643916311563508E+00;
    w[ 5] = 0.453000990550884564085747256463E-02;
  }
  else if (n==7) {
    w[ 0] = 0.971781245099519154149424255939E-03;
    w[ 1] = 0.545155828191270305921785688417E-01;
    w[ 2] = 0.425607252610127800520317466666E+00;
    w[ 3] = 0.810264617556807326764876563813E+00;
    w[ 4] = 0.425607252610127800520317466666E+00;
    w[ 5] = 0.545155828191270305921785688417E-01;
    w[ 6] = 0.971781245099519154149424255939E-03;
  }
  else if (n==8) {
    w[ 0] = 0.199604072211367619206090452544E-03;
    w[ 1] = 0.170779830074134754562030564364E-01;
    w[ 2] = 0.207802325814891879543258620286E+00;
    w[ 3] = 0.661147012558241291030415974496E+00;
    w[ 4] = 0.661147012558241291030415974496E+00;
    w[ 5] = 0.207802325814891879543258620286E+00;
    w[ 6] = 0.170779830074134754562030564364E-01;
    w[ 7] = 0.199604072211367619206090452544E-03;
  }
  else if (n==9) {
    w[ 0] = 0.396069772632643819045862946425E-04;
    w[ 1] = 0.494362427553694721722456597763E-02;
    w[ 2] = 0.884745273943765732879751147476E-01;
    w[ 3] = 0.432651559002555750199812112956E+00;
    w[ 4] = 0.720235215606050957124334723389E+00;
    w[ 5] = 0.432651559002555750199812112956E+00;
    w[ 6] = 0.884745273943765732879751147476E-01;
    w[ 7] = 0.494362427553694721722456597763E-02;
    w[ 8] = 0.396069772632643819045862946425E-04;
  }
  else if (n==10) {
    w[ 0] =  0.764043285523262062915936785960E-05;
    w[ 1] =  0.134364574678123269220156558585E-02;
    w[ 2] =  0.338743944554810631361647312776E-01;
    w[ 3] =  0.240138611082314686416523295006E+00;
    w[ 4] =  0.610862633735325798783564990433E+00;
    w[ 5] =  0.610862633735325798783564990433E+00;
    w[ 6] =  0.240138611082314686416523295006E+00;
    w[ 7] =  0.338743944554810631361647312776E-01;
    w[ 8] =  0.134364574678123269220156558585E-02;
    w[ 9] =  0.764043285523262062915936785960E-05;
  }
  else if (n==11) {
    w[ 0] =  0.143956039371425822033088366032E-05;
    w[ 1] =  0.346819466323345510643413772940E-03;
    w[ 2] =  0.119113954449115324503874202916E-01;
    w[ 3] =  0.117227875167708503381788649308E+00;
    w[ 4] =  0.429359752356125028446073598601E+00;
    w[ 5] =  0.654759286914591779203940657627E+00;
    w[ 6] =  0.429359752356125028446073598601E+00;
    w[ 7] =  0.117227875167708503381788649308E+00;
    w[ 8] =  0.119113954449115324503874202916E-01;
    w[ 9] =  0.346819466323345510643413772940E-03;
    w[10] =  0.143956039371425822033088366032E-05;
  }
  else if (n==12) {
    w[ 0] =  0.265855168435630160602311400877E-06;
    w[ 1] =  0.857368704358785865456906323153E-04;
    w[ 2] =  0.390539058462906185999438432620E-02;
    w[ 3] =  0.516079856158839299918734423606E-01;
    w[ 4] =  0.260492310264161129233396139765E+00;
    w[ 5] =  0.570135236262479578347113482275E+00;
    w[ 6] =  0.570135236262479578347113482275E+00;
    w[ 7] =  0.260492310264161129233396139765E+00;
    w[ 8] =  0.516079856158839299918734423606E-01;
    w[ 9] =  0.390539058462906185999438432620E-02;
    w[10] =  0.857368704358785865456906323153E-04;
    w[11] =  0.265855168435630160602311400877E-06;
  }
  else if (n==13) {
    w[ 0] =  0.482573185007313108834997332342E-07;
    w[ 1] =  0.204303604027070731248669432937E-04;
    w[ 2] =  0.120745999271938594730924899224E-02;
    w[ 3] =  0.208627752961699392166033805050E-01;
    w[ 4] =  0.140323320687023437762792268873E+00;
    w[ 5] =  0.421616296898543221746893558568E+00;
    w[ 6] =  0.604393187921161642342099068579E+00;
    w[ 7] =  0.421616296898543221746893558568E+00;
    w[ 8] =  0.140323320687023437762792268873E+00;
    w[ 9] =  0.208627752961699392166033805050E-01;
    w[10] =  0.120745999271938594730924899224E-02;
    w[11] =  0.204303604027070731248669432937E-04;
    w[12] =  0.482573185007313108834997332342E-07;
  }
  else if (n==14) {
    w[ 0] =  0.862859116812515794532041783429E-08;
    w[ 1] =  0.471648435501891674887688950105E-05;
    w[ 2] =  0.355092613551923610483661076691E-03;
    w[ 3] =  0.785005472645794431048644334608E-02;
    w[ 4] =  0.685055342234652055387163312367E-01;
    w[ 5] =  0.273105609064246603352569187026E+00;
    w[ 6] =  0.536405909712090149794921296776E+00;
    w[ 7] =  0.536405909712090149794921296776E+00;
    w[ 8] =  0.273105609064246603352569187026E+00;
    w[ 9] =  0.685055342234652055387163312367E-01;
    w[10] =  0.785005472645794431048644334608E-02;
    w[11] =  0.355092613551923610483661076691E-03;
    w[12] =  0.471648435501891674887688950105E-05;
    w[13] =  0.862859116812515794532041783429E-08;
  }
  else if (n==15) {
    w[ 0] =  0.152247580425351702016062666965E-08;
    w[ 1] =  0.105911554771106663577520791055E-05;
    w[ 2] =  0.100004441232499868127296736177E-03;
    w[ 3] =  0.277806884291277589607887049229E-02;
    w[ 4] =  0.307800338725460822286814158758E-01;
    w[ 5] =  0.158488915795935746883839384960E+00;
    w[ 6] =  0.412028687498898627025891079568E+00;
    w[ 7] =  0.564100308726417532852625797340E+00;
    w[ 8] =  0.412028687498898627025891079568E+00;
    w[ 9] =  0.158488915795935746883839384960E+00;
    w[10] =  0.307800338725460822286814158758E-01;
    w[11] =  0.277806884291277589607887049229E-02;
    w[12] =  0.100004441232499868127296736177E-03;
    w[13] =  0.105911554771106663577520791055E-05;
    w[14] =  0.152247580425351702016062666965E-08;
  }
  else if (n==16) {
    w[ 0] =  0.265480747401118224470926366050E-09;
    w[ 1] =  0.232098084486521065338749423185E-06;
    w[ 2] =  0.271186009253788151201891432244E-04;
    w[ 3] =  0.932284008624180529914277305537E-03;
    w[ 4] =  0.128803115355099736834642999312E-01;
    w[ 5] =  0.838100413989858294154207349001E-01;
    w[ 6] =  0.280647458528533675369463335380E+00;
    w[ 7] =  0.507929479016613741913517341791E+00;
    w[ 8] =  0.507929479016613741913517341791E+00;
    w[ 9] =  0.280647458528533675369463335380E+00;
    w[10] =  0.838100413989858294154207349001E-01;
    w[11] =  0.128803115355099736834642999312E-01;
    w[12] =  0.932284008624180529914277305537E-03;
    w[13] =  0.271186009253788151201891432244E-04;
    w[14] =  0.232098084486521065338749423185E-06;
    w[15] =  0.265480747401118224470926366050E-09;
  }
  else if (n==17) {
    w[ 0] =  0.458057893079863330580889281222E-10;
    w[ 1] =  0.497707898163079405227863353715E-07;
    w[ 2] =  0.711228914002130958353327376218E-05;
    w[ 3] =  0.298643286697753041151336643059E-03;
    w[ 4] =  0.506734995762753791170069495879E-02;
    w[ 5] =  0.409200341495762798094994877854E-01;
    w[ 6] =  0.172648297670097079217645196219E+00;
    w[ 7] =  0.401826469470411956577635085257E+00;
    w[ 8] =  0.530917937624863560331883103379E+00;
    w[ 9] =  0.401826469470411956577635085257E+00;
    w[10] =  0.172648297670097079217645196219E+00;
    w[11] =  0.409200341495762798094994877854E-01;
    w[12] =  0.506734995762753791170069495879E-02;
    w[13] =  0.298643286697753041151336643059E-03;
    w[14] =  0.711228914002130958353327376218E-05;
    w[15] =  0.497707898163079405227863353715E-07;
    w[16] =  0.458057893079863330580889281222E-10;
  }
  else if (n==18) {
    w[ 0] =  0.782819977211589102925147471012E-11;
    w[ 1] =  0.104672057957920824443559608435E-07;
    w[ 2] =  0.181065448109343040959702385911E-05;
    w[ 3] =  0.918112686792940352914675407371E-04;
    w[ 4] =  0.188852263026841789438175325426E-02;
    w[ 5] =  0.186400423875446519219315221973E-01;
    w[ 6] =  0.973017476413154293308537234155E-01;
    w[ 7] =  0.284807285669979578595606820713E+00;
    w[ 8] =  0.483495694725455552876410522141E+00;
    w[ 9] =  0.483495694725455552876410522141E+00;
    w[10] =  0.284807285669979578595606820713E+00;
    w[11] =  0.973017476413154293308537234155E-01;
    w[12] =  0.186400423875446519219315221973E-01;
    w[13] =  0.188852263026841789438175325426E-02;
    w[14] =  0.918112686792940352914675407371E-04;
    w[15] =  0.181065448109343040959702385911E-05;
    w[16] =  0.104672057957920824443559608435E-07;
    w[17] =  0.782819977211589102925147471012E-11;
  }
  else if (n==19) {
    w[ 0] =  0.132629709449851575185289154385E-11;
    w[ 1] =  0.216305100986355475019693077221E-08;
    w[ 2] =  0.448824314722312295179447915594E-06;
    w[ 3] =  0.272091977631616257711941025214E-04;
    w[ 4] =  0.670877521407181106194696282100E-03;
    w[ 5] =  0.798886677772299020922211491861E-02;
    w[ 6] =  0.508103869090520673569908110358E-01;
    w[ 7] =  0.183632701306997074156148485766E+00;
    w[ 8] =  0.391608988613030244504042313621E+00;
    w[ 9] =  0.502974888276186530840731361096E+00;
    w[10] =  0.391608988613030244504042313621E+00;
    w[11] =  0.183632701306997074156148485766E+00;
    w[12] =  0.508103869090520673569908110358E-01;
    w[13] =  0.798886677772299020922211491861E-02;
    w[14] =  0.670877521407181106194696282100E-03;
    w[15] =  0.272091977631616257711941025214E-04;
    w[16] =  0.448824314722312295179447915594E-06;
    w[17] =  0.216305100986355475019693077221E-08;
    w[18] =  0.132629709449851575185289154385E-11;
  }
  else if (n==20) {
    w[ 0] =  0.222939364553415129252250061603E-12;
    w[ 1] =  0.439934099227318055362885145547E-09;
    w[ 2] =  0.108606937076928169399952456345E-06;
    w[ 3] =  0.780255647853206369414599199965E-05;
    w[ 4] =  0.228338636016353967257145917963E-03;
    w[ 5] =  0.324377334223786183218324713235E-02;
    w[ 6] =  0.248105208874636108821649525589E-01;
    w[ 7] =  0.109017206020023320013755033535E+00;
    w[ 8] =  0.286675505362834129719659706228E+00;
    w[ 9] =  0.462243669600610089650328639861E+00;
    w[10] =  0.462243669600610089650328639861E+00;
    w[11] =  0.286675505362834129719659706228E+00;
    w[12] =  0.109017206020023320013755033535E+00;
    w[13] =  0.248105208874636108821649525589E-01;
    w[14] =  0.324377334223786183218324713235E-02;
    w[15] =  0.228338636016353967257145917963E-03;
    w[16] =  0.780255647853206369414599199965E-05;
    w[17] =  0.108606937076928169399952456345E-06;
    w[18] =  0.439934099227318055362885145547E-09;
    w[19] =  0.222939364553415129252250061603E-12;
  }
  else {
    std::cerr << "\n";
    std::cerr << "HERMITE_LOOKUP_WEIGHTS - Fatal error!\n";
    std::cerr << "  Illegal value of N = " << n << "\n";
    std::cerr << "  Legal values are 1 through 20.\n";
    std::exit(1);
  }

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::imtqlx ( int n, Scalar d[], Scalar e[], Scalar z[] )
//****************************************************************************
//
//  Purpose:
//
//    IMTQLX diagonalizes a symmetric tridiagonal matrix.
//
//  Discussion:
//
//    This routine is a slightly modified version of the EISPACK routine to 
//    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
//
//    The authors thank the authors of EISPACK for permission to use this
//    routine. 
//
//    It has been modified to produce the product Q' * Z, where Z is an input 
//    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
//    The changes consist (essentially) of applying the orthogonal transformations
//    directly to Z as they are generated.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//    Roger Martin, James Wilkinson,
//    The Implicit QL Algorithm,
//    Numerische Mathematik,
//    Volume 12, Number 5, December 1968, pages 377-383.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, Scalar D(N), the diagonal entries of the matrix.
//    On output, the information in D has been overwritten.
//
//    Input/output, Scalar E(N), the subdiagonal entries of the 
//    matrix, in entries E(1) through E(N-1).  On output, the information in
//    E has been overwritten.
//
//    Input/output, Scalar Z(N).  On input, a vector.  On output,
//    the value of Q' * Z, where Q is the matrix that diagonalizes the
//    input symmetric tridiagonal matrix.
//
{
  Scalar b = 0, c = 0, f = 0, g = 0, p = 0, r = 0, s = 0;
  int i = 0, ii = 0, j = 0, k = 0, l = 0, mml = 0, m = 0, itn = 30;
  Scalar prec = IntrepidBurkardtRules::r8_epsilon(1.0); //2.22E-16;

  if (n==1) {
    return;
  }

  e[n-1] = 0.0;

  for (l=1;l<=n;l++) {
    j = 0;
    for ( ; ; ) {
      for (m=l;m<=n;m++) {
        if (m==n) {
          break;
        }
        if (std::abs(e[m-1])<=prec*(std::abs(d[m-1])+std::abs(d[m]))) {
          break;
        }
      }
      p = d[l-1];
      if (m==l) {
        break;
      }
      if (itn<=j) {
        std::cerr << "\n";
        std::cerr << "IMTQLX - Fatal error!\n";
        std::cerr << "  Iteration limit exceeded\n";
        std::exit(1);
      }
      j   = j+1;
      g   = (d[l]-p)/(2.0*e[l-1]);
      r   = std::sqrt(g*g+1.0);
      g   = d[m-1]-p+e[l-1]/(g+std::abs(r)*IntrepidBurkardtRules::r8_sign(g));
      s   = 1.0;
      c   = 1.0;
      p   = 0.0;
      mml = m-l;

      for (ii=1;ii<=mml;ii++) {
        i = m-ii;
        f = s*e[i-1];
        b = c*e[i-1];

        if (std::abs(g)<=std::abs(f)) {
          c    = g/f;
          r    = std::sqrt(c*c+1.0);
          e[i] = f*r;
          s    = 1.0/r;
          c    = c*s;
        }
        else {
          s    = f/g;
          r    = std::sqrt(s*s+1.0);
          e[i] = g*r;
          c    = 1.0/r;
          s    = s*c;
        }
        g      = d[i]-p;
        r      = (d[i-1]-g)*s+2.0*c*b;
        p      = s*r;
        d[i]   = g+p;
        g      = c*r-b;
        f      = z[i];
        z[i]   = s*z[i-1]+c*f;
        z[i-1] = c*z[i-1]-s*f;
      }
      d[l-1] = d[l-1]-p;
      e[l-1] = g;
      e[m-1] = 0.0;
    }
  }
//
//  Sorting.
//
  for (ii=2;ii<=m;ii++) {
    i = ii-1;
    k = i;
    p = d[i-1];

    for (j=ii;j<=n;j++) {
      if (d[j-1]<p) {
         k = j;
         p = d[j-1];
      }
    }

    if (k!=i) {
      d[k-1] = d[i-1];
      d[i-1] = p;
      p      = z[i-1];
      z[i-1] = z[k-1];
      z[k-1] = p;
    }
  }
  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::laguerre_compute ( int n, Scalar x[], Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    LAGUERRE_COMPUTE: Laguerre quadrature rule by the Elhay-Kautsky method.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 April 2011
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, Scalar X[N], the abscissas.
//
//    Output, Scalar W[N], the weights.
//
{
  Scalar *bj;
  int i;
  Scalar zemu;
//
//  Define the zero-th moment.
//
  zemu = 1.0;
//
//  Define the Jacobi matrix.
//
  bj = new Scalar[n];

  for (i=0;i<n;i++) {
    bj[i] = (Scalar)(i+1);
  }

  for (i=0;i<n;i++) {
    x[i] = (Scalar)(2*i+1);
  }

  w[0] = std::sqrt(zemu);

  for (i=1;i<n;i++) {
    w[i] = 0.0;
  }
//
//  Diagonalize the Jacobi matrix.
//
  IntrepidBurkardtRules::imtqlx(n,x,bj,w);

  for (i=0;i<n;i++) {
    w[i] = w[i]*w[i];
  }

  delete [] bj;

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::laguerre_compute_points ( int order, Scalar x[] )
//****************************************************************************
//
//  Purpose:
//
//    LAGUERRE_COMPUTE_POINTS computes Laguerre quadrature points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order.
//
//    Output, Scalar X[ORDER], the abscissas.
//
{
  Scalar *w; w = new Scalar[order];
  IntrepidBurkardtRules::laguerre_compute ( order, x, w );
  delete [] w;

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::laguerre_compute_weights ( int order, Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    LAGUERRE_COMPUTE_WEIGHTS computes Laguerre quadrature weights.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order.
//
//    Output, Scalar W[ORDER], the weights.
//
{
  Scalar *x; x = new Scalar[order];
  IntrepidBurkardtRules::laguerre_compute ( order, x, w );
  delete [] x;

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::laguerre_lookup ( int n, Scalar x[], Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    LAGUERRE_LOOKUP looks up abscissas and weights for Laguerre quadrature.
//
//  Discussion:
//
//    The abscissas are the zeroes of the Laguerre polynomial L(N)(X).
//
//    The integral:
//
//      Integral ( 0 <= X < +oo ) exp ( -X ) * F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= N ) W(I) * f ( X(I) )
//
//    The integral:
//
//      Integral ( 0 <= X < +oo ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= N ) W(I) * exp ( X(I) ) * f ( X(I) )
//
//    Mathematica can numerically estimate the abscissas for the
//    n-th order polynomial to p digits of precision by the command:
//
//      NSolve [ LaguerreL[n,x] == 0, x, p ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Vladimir Krylov,
//    Approximate Calculation of Integrals,
//    Dover, 2006,
//    ISBN: 0486445798,
//    LC: QA311.K713.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be between 1 and 20.
//
//    Output, Scalar X[N], the abscissas.
//
//    Output, Scalar W[N], the weights.
//
{
  IntrepidBurkardtRules::laguerre_lookup_points ( n, x );
  IntrepidBurkardtRules::laguerre_lookup_weights ( n, w );

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::laguerre_lookup_points ( int n, Scalar x[] )
//****************************************************************************
//
//  Purpose:
//
//    LAGUERRE_LOOKUP_POINTS looks up abscissas for Laguerre quadrature.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Vladimir Krylov,
//    Approximate Calculation of Integrals,
//    Dover, 2006,
//    ISBN: 0486445798,
//    LC: QA311.K713.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be between 1 and 20.
//
//    Output, Scalar X[N], the abscissas.
//
{
  if (n==1) {
    x[ 0] =  1.00000000000000000000000000000E+00;
  }
  else if (n==2) {
    x[ 0] = 0.585786437626904951198311275790E+00;
    x[ 1] = 3.41421356237309504880168872421E+00;
  }
  else if (n==3) {
    x[ 0] = 0.415774556783479083311533873128E+00;
    x[ 1] = 2.29428036027904171982205036136E+00;
    x[ 2] = 6.28994508293747919686641576551E+00;
  }
  else if (n==4) {
    x[ 0] = 0.322547689619392311800361459104E+00;
    x[ 1] = 1.74576110115834657568681671252E+00;
    x[ 2] = 4.53662029692112798327928538496E+00;
    x[ 3] = 9.39507091230113312923353644342E+00;
  }
  else if (n==5) {
    x[ 0] = 0.263560319718140910203061943361E+00;
    x[ 1] = 1.41340305910651679221840798019E+00;
    x[ 2] = 3.59642577104072208122318658878E+00;
    x[ 3] = 7.08581000585883755692212418111E+00;
    x[ 4] = 12.6408008442757826594332193066E+00;
  }
  else if (n==6) {
    x[ 0] = 0.222846604179260689464354826787E+00;
    x[ 1] = 1.18893210167262303074315092194E+00;
    x[ 2] = 2.99273632605931407769132528451E+00;
    x[ 3] = 5.77514356910451050183983036943E+00;
    x[ 4] = 9.83746741838258991771554702994E+00;
    x[ 5] = 15.9828739806017017825457915674E+00;
  }
  else if (n==7) {
    x[ 0] = 0.193043676560362413838247885004E+00;
    x[ 1] = 1.02666489533919195034519944317E+00;
    x[ 2] = 2.56787674495074620690778622666E+00;
    x[ 3] = 4.90035308452648456810171437810E+00;
    x[ 4] = 8.18215344456286079108182755123E+00;
    x[ 5] = 12.7341802917978137580126424582E+00;
    x[ 6] = 19.3957278622625403117125820576E+00;
  }
  else if (n==8) {
    x[ 0] = 0.170279632305100999788861856608E+00;
    x[ 1] = 0.903701776799379912186020223555E+00;
    x[ 2] = 2.25108662986613068930711836697E+00;
    x[ 3] = 4.26670017028765879364942182690E+00;
    x[ 4] = 7.04590540239346569727932548212E+00;
    x[ 5] = 10.7585160101809952240599567880E+00;
    x[ 6] = 15.7406786412780045780287611584E+00;
    x[ 7] = 22.8631317368892641057005342974E+00;
  }
  else if (n==9) {
    x[ 0] = 0.152322227731808247428107073127E+00;
    x[ 1] = 0.807220022742255847741419210952E+00;
    x[ 2] = 2.00513515561934712298303324701E+00;
    x[ 3] = 3.78347397333123299167540609364E+00;
    x[ 4] = 6.20495677787661260697353521006E+00;
    x[ 5] = 9.37298525168757620180971073215E+00;
    x[ 6] = 13.4662369110920935710978818397E+00;
    x[ 7] = 18.8335977889916966141498992996E+00;
    x[ 8] = 26.3740718909273767961410072937E+00;
  }
  else if (n==10) {
    x[ 0] = 0.137793470540492430830772505653E+00;
    x[ 1] = 0.729454549503170498160373121676E+00;
    x[ 2] = 1.80834290174031604823292007575E+00;
    x[ 3] = 3.40143369785489951448253222141E+00;
    x[ 4] = 5.55249614006380363241755848687E+00;
    x[ 5] = 8.33015274676449670023876719727E+00;
    x[ 6] = 11.8437858379000655649185389191E+00;
    x[ 7] = 16.2792578313781020995326539358E+00;
    x[ 8] = 21.9965858119807619512770901956E+00;
    x[ 9] = 29.9206970122738915599087933408E+00;
  }
  else if (n==11) {
    x[ 0] = 0.125796442187967522675794577516E+00;
    x[ 1] = 0.665418255839227841678127839420E+00;
    x[ 2] = 1.64715054587216930958700321365E+00;
    x[ 3] = 3.09113814303525495330195934259E+00;
    x[ 4] = 5.02928440157983321236999508366E+00;
    x[ 5] = 7.50988786380661681941099714450E+00;
    x[ 6] = 10.6059509995469677805559216457E+00;
    x[ 7] = 14.4316137580641855353200450349E+00;
    x[ 8] = 19.1788574032146786478174853989E+00;
    x[ 9] = 25.2177093396775611040909447797E+00;
    x[10] = 33.4971928471755372731917259395E+00;
  }
  else if (n==12) {
    x[ 0] = 0.115722117358020675267196428240E+00;
    x[ 1] = 0.611757484515130665391630053042E+00;
    x[ 2] = 1.51261026977641878678173792687E+00;
    x[ 3] = 2.83375133774350722862747177657E+00;
    x[ 4] = 4.59922763941834848460572922485E+00;
    x[ 5] = 6.84452545311517734775433041849E+00;
    x[ 6] = 9.62131684245686704391238234923E+00;
    x[ 7] = 13.0060549933063477203460524294E+00;
    x[ 8] = 17.1168551874622557281840528008E+00;
    x[ 9] = 22.1510903793970056699218950837E+00;
    x[10] = 28.4879672509840003125686072325E+00;
    x[11] = 37.0991210444669203366389142764E+00;
  }
  else if (n==13) {
    x[ 0] = 0.107142388472252310648493376977E+00;
    x[ 1] = 0.566131899040401853406036347177E+00;
    x[ 2] = 1.39856433645101971792750259921E+00;
    x[ 3] = 2.61659710840641129808364008472E+00;
    x[ 4] = 4.23884592901703327937303389926E+00;
    x[ 5] = 6.29225627114007378039376523025E+00;
    x[ 6] = 8.81500194118697804733348868036E+00;
    x[ 7] = 11.8614035888112425762212021880E+00;
    x[ 8] = 15.5107620377037527818478532958E+00;
    x[ 9] = 19.8846356638802283332036594634E+00;
    x[10] = 25.1852638646777580842970297823E+00;
    x[11] = 31.8003863019472683713663283526E+00;
    x[12] = 40.7230086692655795658979667001E+00;
  }
  else if (n==14) {
    x[ 0] = 0.0997475070325975745736829452514E+00;
    x[ 1] = 0.526857648851902896404583451502E+00;
    x[ 2] = 1.30062912125149648170842022116E+00;
    x[ 3] = 2.43080107873084463616999751038E+00;
    x[ 4] = 3.93210282229321888213134366778E+00;
    x[ 5] = 5.82553621830170841933899983898E+00;
    x[ 6] = 8.14024014156514503005978046052E+00;
    x[ 7] = 10.9164995073660188408130510904E+00;
    x[ 8] = 14.2108050111612886831059780825E+00;
    x[ 9] = 18.1048922202180984125546272083E+00;
    x[10] = 22.7233816282696248232280886985E+00;
    x[11] = 28.2729817232482056954158923218E+00;
    x[12] = 35.1494436605924265828643121364E+00;
    x[13] = 44.3660817111174230416312423666E+00;
  }
  else if (n==15) {
    x[ 0] = 0.0933078120172818047629030383672E+00;
    x[ 1] = 0.492691740301883908960101791412E+00;
    x[ 2] = 1.21559541207094946372992716488E+00;
    x[ 3] = 2.26994952620374320247421741375E+00;
    x[ 4] = 3.66762272175143727724905959436E+00;
    x[ 5] = 5.42533662741355316534358132596E+00;
    x[ 6] = 7.56591622661306786049739555812E+00;
    x[ 7] = 10.1202285680191127347927394568E+00;
    x[ 8] = 13.1302824821757235640991204176E+00;
    x[ 9] = 16.6544077083299578225202408430E+00;
    x[10] = 20.7764788994487667729157175676E+00;
    x[11] = 25.6238942267287801445868285977E+00;
    x[12] = 31.4075191697539385152432196202E+00;
    x[13] = 38.5306833064860094162515167595E+00;
    x[14] = 48.0260855726857943465734308508E+00;
  }
  else if (n==16) {
    x[ 0] = 0.0876494104789278403601980973401E+00;
    x[ 1] = 0.462696328915080831880838260664E+00;
    x[ 2] = 1.14105777483122685687794501811E+00;
    x[ 3] = 2.12928364509838061632615907066E+00;
    x[ 4] = 3.43708663389320664523510701675E+00;
    x[ 5] = 5.07801861454976791292305830814E+00;
    x[ 6] = 7.07033853504823413039598947080E+00;
    x[ 7] = 9.43831433639193878394724672911E+00;
    x[ 8] = 12.2142233688661587369391246088E+00;
    x[ 9] = 15.4415273687816170767647741622E+00;
    x[10] = 19.1801568567531348546631409497E+00;
    x[11] = 23.5159056939919085318231872752E+00;
    x[12] = 28.5787297428821403675206137099E+00;
    x[13] = 34.5833987022866258145276871778E+00;
    x[14] = 41.9404526476883326354722330252E+00;
    x[15] = 51.7011603395433183643426971197E+00;
  }
  else if (n==17) {
    x[ 0] = 0.0826382147089476690543986151980E+00;
    x[ 1] = 0.436150323558710436375959029847E+00;
    x[ 2] = 1.07517657751142857732980316755E+00;
    x[ 3] = 2.00519353164923224070293371933E+00;
    x[ 4] = 3.23425612404744376157380120696E+00;
    x[ 5] = 4.77351351370019726480932076262E+00;
    x[ 6] = 6.63782920536495266541643929703E+00;
    x[ 7] = 8.84668551116980005369470571184E+00;
    x[ 8] = 11.4255293193733525869726151469E+00;
    x[ 9] = 14.4078230374813180021982874959E+00;
    x[10] = 17.8382847307011409290658752412E+00;
    x[11] = 21.7782682577222653261749080522E+00;
    x[12] = 26.3153178112487997766149598369E+00;
    x[13] = 31.5817716804567331343908517497E+00;
    x[14] = 37.7960938374771007286092846663E+00;
    x[15] = 45.3757165339889661829258363215E+00;
    x[16] = 55.3897517898396106640900199790E+00;
  }
  else if (n==18) {
    x[ 0] = 0.0781691666697054712986747615334E+00;
    x[ 1] = 0.412490085259129291039101536536E+00;
    x[ 2] = 1.01652017962353968919093686187E+00;
    x[ 3] = 1.89488850996976091426727831954E+00;
    x[ 4] = 3.05435311320265975115241130719E+00;
    x[ 5] = 4.50420553888989282633795571455E+00;
    x[ 6] = 6.25672507394911145274209116326E+00;
    x[ 7] = 8.32782515660563002170470261564E+00;
    x[ 8] = 10.7379900477576093352179033397E+00;
    x[ 9] = 13.5136562075550898190863812108E+00;
    x[10] = 16.6893062819301059378183984163E+00;
    x[11] = 20.3107676262677428561313764553E+00;
    x[12] = 24.4406813592837027656442257980E+00;
    x[13] = 29.1682086625796161312980677805E+00;
    x[14] = 34.6279270656601721454012429438E+00;
    x[15] = 41.0418167728087581392948614284E+00;
    x[16] = 48.8339227160865227486586093290E+00;
    x[17] = 59.0905464359012507037157810181E+00;
  }
  else if (n==19) {
    x[ 0] = 0.0741587837572050877131369916024E+00;
    x[ 1] = 0.391268613319994607337648350299E+00;
    x[ 2] = 0.963957343997958058624878377130E+00;
    x[ 3] = 1.79617558206832812557725825252E+00;
    x[ 4] = 2.89365138187378399116494713237E+00;
    x[ 5] = 4.26421553962776647436040018167E+00;
    x[ 6] = 5.91814156164404855815360191408E+00;
    x[ 7] = 7.86861891533473373105668358176E+00;
    x[ 8] = 10.1324237168152659251627415800E+00;
    x[ 9] = 12.7308814638423980045092979656E+00;
    x[10] = 15.6912783398358885454136069861E+00;
    x[11] = 19.0489932098235501532136429732E+00;
    x[12] = 22.8508497608294829323930586693E+00;
    x[13] = 27.1606693274114488789963947149E+00;
    x[14] = 32.0691222518622423224362865906E+00;
    x[15] = 37.7129058012196494770647508283E+00;
    x[16] = 44.3173627958314961196067736013E+00;
    x[17] = 52.3129024574043831658644222420E+00;
    x[18] = 62.8024231535003758413504690673E+00;
  }
  else if (n==20) {
    x[ 0] = 0.0705398896919887533666890045842E+00;
    x[ 1] = 0.372126818001611443794241388761E+00;
    x[ 2] = 0.916582102483273564667716277074E+00;
    x[ 3] = 1.70730653102834388068768966741E+00;
    x[ 4] = 2.74919925530943212964503046049E+00;
    x[ 5] = 4.04892531385088692237495336913E+00;
    x[ 6] = 5.61517497086161651410453988565E+00;
    x[ 7] = 7.45901745367106330976886021837E+00;
    x[ 8] = 9.59439286958109677247367273428E+00;
    x[ 9] = 12.0388025469643163096234092989E+00;
    x[10] = 14.8142934426307399785126797100E+00;
    x[11] = 17.9488955205193760173657909926E+00;
    x[12] = 21.4787882402850109757351703696E+00;
    x[13] = 25.4517027931869055035186774846E+00;
    x[14] = 29.9325546317006120067136561352E+00;
    x[15] = 35.0134342404790000062849359067E+00;
    x[16] = 40.8330570567285710620295677078E+00;
    x[17] = 47.6199940473465021399416271529E+00;
    x[18] = 55.8107957500638988907507734445E+00;
    x[19] = 66.5244165256157538186403187915E+00;
  }
  else {
    std::cerr << "\n";
    std::cerr << "LAGUERRE_LOOKUP_POINTS - Fatal error!\n";
    std::cerr << "  Illegal value of N = " << n << "\n";
    std::cerr << "  Legal values are 1 through 20.\n";
    std::exit(1);
  }

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::laguerre_lookup_weights ( int n, Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    LAGUERRE_LOOKUP_WEIGHTS looks up weights for Laguerre quadrature.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Vladimir Krylov,
//    Approximate Calculation of Integrals,
//    Dover, 2006,
//    ISBN: 0486445798,
//    LC: QA311.K713.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be between 1 and 20.
//
//    Output, Scalar W[N], the weights.
//
{
  if (n==1) {
    w[ 0] =  1.00000000000000000000000000000E+00;
  }
  else if (n==2) {
    w[ 0] = 0.85355339059327376220042218105E+00;
    w[ 1] = 0.146446609406726237799577818948E+00;
  }
  else if (n==3) {
    w[ 0] = 0.71109300992917301544959019114E+00;
    w[ 1] = 0.27851773356924084880144488846E+00;
    w[ 2] = 0.010389256501586135748964920401E+00;
  }
  else if (n==4) {
    w[ 0] = 0.60315410434163360163596602382E+00;
    w[ 1] = 0.35741869243779968664149201746E+00;
    w[ 2] = 0.03888790851500538427243816816E+00;
    w[ 3] = 0.0005392947055613274501037905676E+00;
  }
  else if (n==5) {
    w[ 0] = 0.52175561058280865247586092879E+00;
    w[ 1] = 0.3986668110831759274541333481E+00;
    w[ 2] = 0.0759424496817075953876533114E+00;
    w[ 3] = 0.00361175867992204845446126257E+00;
    w[ 4] = 0.00002336997238577622789114908455E+00;
  }
  else if (n==6) {
    w[ 0] = 0.45896467394996359356828487771E+00;
    w[ 1] = 0.4170008307721209941133775662E+00;
    w[ 2] = 0.1133733820740449757387061851E+00;
    w[ 3] = 0.01039919745314907489891330285E+00;
    w[ 4] = 0.000261017202814932059479242860E+00;
    w[ 5] = 8.98547906429621238825292053E-07;
  }
  else if (n==7) {
    w[ 0] = 0.40931895170127390213043288002E+00;
    w[ 1] = 0.4218312778617197799292810054E+00;
    w[ 2] = 0.1471263486575052783953741846E+00;
    w[ 3] = 0.0206335144687169398657056150E+00;
    w[ 4] = 0.00107401014328074552213195963E+00;
    w[ 5] = 0.0000158654643485642012687326223E+00;
    w[ 6] = 3.17031547899558056227132215E-08;
  }
  else if (n==8) {
    w[ 0] = 0.36918858934163752992058283938E+00;
    w[ 1] = 0.4187867808143429560769785813E+00;
    w[ 2] = 0.175794986637171805699659867E+00;
    w[ 3] = 0.033343492261215651522132535E+00;
    w[ 4] = 0.0027945362352256725249389241E+00;
    w[ 5] = 0.00009076508773358213104238501E+00;
    w[ 6] = 8.4857467162725315448680183E-07;
    w[ 7] = 1.04800117487151038161508854E-09;
  }
  else if (n==9) {
    w[ 0] = 0.336126421797962519673467717606E+00;
    w[ 1] = 0.411213980423984387309146942793E+00;
    w[ 2] = 0.199287525370885580860575607212E+00;
    w[ 3] = 0.0474605627656515992621163600479E+00;
    w[ 4] = 0.00559962661079458317700419900556E+00;
    w[ 5] = 0.000305249767093210566305412824291E+00;
    w[ 6] = 6.59212302607535239225572284875E-06;
    w[ 7] = 4.1107693303495484429024104033E-08;
    w[ 8] = 3.29087403035070757646681380323E-11;
  }
  else if (n==10) {
    w[ 0] = 0.30844111576502014154747083468E+00;
    w[ 1] = 0.4011199291552735515157803099E+00;
    w[ 2] = 0.218068287611809421588648523E+00;
    w[ 3] = 0.062087456098677747392902129E+00;
    w[ 4] = 0.009501516975181100553839072E+00;
    w[ 5] = 0.0007530083885875387754559644E+00;
    w[ 6] = 0.00002825923349599565567422564E+00;
    w[ 7] = 4.249313984962686372586577E-07;
    w[ 8] = 1.839564823979630780921535E-09;
    w[ 9] = 9.911827219609008558377547E-13;
  }
  else if (n==11) {
    w[ 0] = 0.28493321289420060505605102472E+00;
    w[ 1] = 0.3897208895278493779375535080E+00;
    w[ 2] = 0.232781831848991333940223796E+00;
    w[ 3] = 0.076564453546196686400854179E+00;
    w[ 4] = 0.014393282767350695091863919E+00;
    w[ 5] = 0.001518880846484873069847776E+00;
    w[ 6] = 0.0000851312243547192259720424E+00;
    w[ 7] = 2.29240387957450407857683E-06;
    w[ 8] = 2.48635370276779587373391E-08;
    w[ 9] = 7.71262693369132047028153E-11;
    w[10] = 2.883775868323623861597778E-14;
  }
  else if (n==12) {
    w[ 0] = 0.26473137105544319034973889206E+00;
    w[ 1] = 0.3777592758731379820244905567E+00;
    w[ 2] = 0.244082011319877564254870818E+00;
    w[ 3] = 0.09044922221168093072750549E+00;
    w[ 4] = 0.02010238115463409652266129E+00;
    w[ 5] = 0.002663973541865315881054158E+00;
    w[ 6] = 0.000203231592662999392121433E+00;
    w[ 7] = 8.3650558568197987453363E-06;
    w[ 8] = 1.66849387654091026116990E-07;
    w[ 9] = 1.34239103051500414552392E-09;
    w[10] = 3.06160163503502078142408E-12;
    w[11] = 8.148077467426241682473119E-16;
  }
  else if (n==13) {
    w[ 0] = 0.24718870842996262134624918596E+00;
    w[ 1] = 0.3656888229005219453067175309E+00;
    w[ 2] = 0.252562420057658502356824289E+00;
    w[ 3] = 0.10347075802418370511421863E+00;
    w[ 4] = 0.02643275441556161577815877E+00;
    w[ 5] = 0.00422039604025475276555209E+00;
    w[ 6] = 0.000411881770472734774892473E+00;
    w[ 7] = 0.0000235154739815532386882897E+00;
    w[ 8] = 7.3173116202490991040105E-07;
    w[ 9] = 1.10884162570398067979151E-08;
    w[10] = 6.7708266922058988406462E-11;
    w[11] = 1.15997995990507606094507E-13;
    w[12] = 2.245093203892758415991872E-17;
  }
  else if (n==14) {
    w[ 0] = 0.23181557714486497784077486110E+00;
    w[ 1] = 0.3537846915975431518023313013E+00;
    w[ 2] = 0.258734610245428085987320561E+00;
    w[ 3] = 0.11548289355692321008730499E+00;
    w[ 4] = 0.03319209215933736003874996E+00;
    w[ 5] = 0.00619286943700661021678786E+00;
    w[ 6] = 0.00073989037786738594242589E+00;
    w[ 7] = 0.000054907194668416983785733E+00;
    w[ 8] = 2.4095857640853774967578E-06;
    w[ 9] = 5.801543981676495180886E-08;
    w[10] = 6.819314692484974119616E-10;
    w[11] = 3.2212077518948479398089E-12;
    w[12] = 4.2213524405165873515980E-15;
    w[13] = 6.05237502228918880839871E-19;
  }
  else if (n==15) {
    w[ 0] = 0.21823488594008688985641323645E+00;
    w[ 1] = 0.3422101779228833296389489568E+00;
    w[ 2] = 0.263027577941680097414812275E+00;
    w[ 3] = 0.12642581810593053584303055E+00;
    w[ 4] = 0.04020686492100091484158548E+00;
    w[ 5] = 0.00856387780361183836391576E+00;
    w[ 6] = 0.00121243614721425207621921E+00;
    w[ 7] = 0.00011167439234425194199258E+00;
    w[ 8] = 6.459926762022900924653E-06;
    w[ 9] = 2.226316907096272630332E-07;
    w[10] = 4.227430384979365007351E-09;
    w[11] = 3.921897267041089290385E-11;
    w[12] = 1.4565152640731264063327E-13;
    w[13] = 1.4830270511133013354616E-16;
    w[14] = 1.60059490621113323104998E-20;
  }
  else if (n==16) {
    w[ 0] = 0.20615171495780099433427363674E+00;
    w[ 1] = 0.3310578549508841659929830987E+00;
    w[ 2] = 0.265795777644214152599502021E+00;
    w[ 3] = 0.13629693429637753997554751E+00;
    w[ 4] = 0.0473289286941252189780623E+00;
    w[ 5] = 0.0112999000803394532312490E+00;
    w[ 6] = 0.0018490709435263108642918E+00;
    w[ 7] = 0.00020427191530827846012602E+00;
    w[ 8] = 0.00001484458687398129877135E+00;
    w[ 9] = 6.828319330871199564396E-07;
    w[10] = 1.881024841079673213882E-08;
    w[11] = 2.862350242973881619631E-10;
    w[12] = 2.127079033224102967390E-12;
    w[13] = 6.297967002517867787174E-15;
    w[14] = 5.050473700035512820402E-18;
    w[15] = 4.1614623703728551904265E-22;
  }
  else if (n==17) {
    w[ 0] = 0.19533220525177083214592729770E+00;
    w[ 1] = 0.3203753572745402813366256320E+00;
    w[ 2] = 0.267329726357171097238809604E+00;
    w[ 3] = 0.14512985435875862540742645E+00;
    w[ 4] = 0.0544369432453384577793806E+00;
    w[ 5] = 0.0143572977660618672917767E+00;
    w[ 6] = 0.0026628247355727725684324E+00;
    w[ 7] = 0.0003436797271562999206118E+00;
    w[ 8] = 0.00003027551783782870109437E+00;
    w[ 9] = 1.768515053231676895381E-06;
    w[10] = 6.57627288681043332199E-08;
    w[11] = 1.469730932159546790344E-09;
    w[12] = 1.81691036255544979555E-11;
    w[13] = 1.095401388928687402976E-13;
    w[14] = 2.617373882223370421551E-16;
    w[15] = 1.6729356931461546908502E-19;
    w[16] = 1.06562631627404278815253E-23;
  }
  else if (n==18) {
    w[ 0] = 0.18558860314691880562333775228E+00;
    w[ 1] = 0.3101817663702252936495975957E+00;
    w[ 2] = 0.267866567148536354820854395E+00;
    w[ 3] = 0.15297974746807490655384308E+00;
    w[ 4] = 0.0614349178609616527076780E+00;
    w[ 5] = 0.0176872130807729312772600E+00;
    w[ 6] = 0.0036601797677599177980266E+00;
    w[ 7] = 0.0005406227870077353231284E+00;
    w[ 8] = 0.0000561696505121423113818E+00;
    w[ 9] = 4.01530788370115755859E-06;
    w[10] = 1.91466985667567497969E-07;
    w[11] = 5.8360952686315941292E-09;
    w[12] = 1.07171126695539012773E-10;
    w[13] = 1.08909871388883385562E-12;
    w[14] = 5.38666474837830887608E-15;
    w[15] = 1.049865978035703408779E-17;
    w[16] = 5.405398451631053643566E-21;
    w[17] = 2.6916532692010286270838E-25;
  }
  else if (n==19) {
    w[ 0] = 0.17676847491591250225103547981E+00;
    w[ 1] = 0.3004781436072543794821568077E+00;
    w[ 2] = 0.267599547038175030772695441E+00;
    w[ 3] = 0.15991337213558021678551215E+00;
    w[ 4] = 0.0682493799761491134552355E+00;
    w[ 5] = 0.0212393076065443249244062E+00;
    w[ 6] = 0.0048416273511483959672501E+00;
    w[ 7] = 0.0008049127473813667665946E+00;
    w[ 8] = 0.0000965247209315350170843E+00;
    w[ 9] = 8.20730525805103054409E-06;
    w[10] = 4.8305667247307725394E-07;
    w[11] = 1.90499136112328569994E-08;
    w[12] = 4.8166846309280615577E-10;
    w[13] = 7.3482588395511443768E-12;
    w[14] = 6.2022753875726163989E-14;
    w[15] = 2.54143084301542272372E-16;
    w[16] = 4.07886129682571235007E-19;
    w[17] = 1.707750187593837061004E-22;
    w[18] = 6.715064649908189959990E-27;
  }
  else if (n==20) {
    w[ 0] = 0.168746801851113862149223899689E+00;
    w[ 1] = 0.291254362006068281716795323812E+00;
    w[ 2] = 0.266686102867001288549520868998E+00;
    w[ 3] = 0.166002453269506840031469127816E+00;
    w[ 4] = 0.0748260646687923705400624639615E+00;
    w[ 5] = 0.0249644173092832210728227383234E+00;
    w[ 6] = 0.00620255084457223684744754785395E+00;
    w[ 7] = 0.00114496238647690824203955356969E+00;
    w[ 8] = 0.000155741773027811974779809513214E+00;
    w[ 9] = 0.0000154014408652249156893806714048E+00;
    w[10] = 1.08648636651798235147970004439E-06;
    w[11] = 5.33012090955671475092780244305E-08;
    w[12] = 1.7579811790505820035778763784E-09;
    w[13] = 3.72550240251232087262924585338E-11;
    w[14] = 4.76752925157819052449488071613E-13;
    w[15] = 3.37284424336243841236506064991E-15;
    w[16] = 1.15501433950039883096396247181E-17;
    w[17] = 1.53952214058234355346383319667E-20;
    w[18] = 5.28644272556915782880273587683E-24;
    w[19] = 1.65645661249902329590781908529E-28;
  }
  else {
    std::cerr << "\n";
    std::cerr << "LAGUERRE_LOOKUP_WEIGHTS - Fatal error!\n";
    std::cerr << "  Illegal value of N = " << n << "\n";
    std::cerr << "  Legal values are 1 through 20.\n";
    std::exit(1);
  }

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::legendre_compute ( int n, Scalar x[], Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    LEGENDRE_COMPUTE: Legendre quadrature rule by the Elhay-Kautsky method.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2011
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, Scalar X[N], the abscissas.
//
//    Output, Scalar W[N], the weights.
//
{
  Scalar *bj;
  int i;
  Scalar zemu;
//
//  Define the zero-th moment.
//
  zemu = 2.0;
//
//  Define the Jacobi matrix.
//
  bj = new Scalar[n];

  for (i=0;i<n;i++) {
    bj[i] = (Scalar)((i+1)*(i+1))/(Scalar)(4*(i+1)*(i+1)-1);
    bj[i] = std::sqrt(bj[i]);
  }

  for (i=0;i<n;i++) {
    x[i] = 0.0;
  }

  w[0] = std::sqrt(zemu);

  for (i=1;i<n;i++) {
    w[i] = 0.0;
  }
//
//  Diagonalize the Jacobi matrix.
//
  IntrepidBurkardtRules::imtqlx(n,x,bj,w);

  for (i=0;i<n;i++) {
    w[i] = w[i]*w[i];
  }

  // Ensure that zero is actually zero.
  if (n%2) {
    int ind = (int)((Scalar)n/2.0);
    x[ind]  = 0.0;
  }

  delete [] bj;

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::legendre_compute_points ( int n, Scalar x[] )
//****************************************************************************
//
//  Purpose:
//
//    LEGENDRE_COMPUTE_POINTS computes Legendre quadrature points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, Scalar X[N], the abscissas.
//
{
  Scalar *w; w= new Scalar[n];
  IntrepidBurkardtRules::legendre_compute ( n, x, w );
  delete [] w;

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::legendre_compute_weights ( int n, Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    LEGENDRE_COMPUTE_WEIGHTS computes Legendre quadrature weights.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, Scalar W[N], the weights.
//
{
  Scalar *x; x = new Scalar[n];
  IntrepidBurkardtRules::legendre_compute ( n, x, w );
  delete [] x;

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::legendre_lookup ( int n, Scalar x[], Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    LEGENDRE_LOOKUP looks up abscissas and weights for Gauss-Legendre quadrature.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Vladimir Krylov,
//    Approximate Calculation of Integrals,
//    Dover, 2006,
//    ISBN: 0486445798.
//    LC: QA311.K713.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3,
//    LC: QA47.M315.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be between 1 and 33.
//
//    Output, Scalar X[N], the abscissas.
//
//    Output, Scalar W[N], the abscissas.
//
{
  IntrepidBurkardtRules::legendre_lookup_points ( n, x );
  IntrepidBurkardtRules::legendre_lookup_weights ( n, w );

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::legendre_lookup_points ( int n, Scalar x[] )
//****************************************************************************
//
//  Purpose:
//
//    LEGENDRE_LOOKUP_POINTS looks up abscissas for Gauss-Legendre quadrature.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Vladimir Krylov,
//    Approximate Calculation of Integrals,
//    Dover, 2006,
//    ISBN: 0486445798.
//    LC: QA311.K713.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3,
//    LC: QA47.M315.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be between 1 and 33.
//
//    Output, Scalar X[N], the abscissas.
//
{
  if (n==1) {
    x[ 0] = 0.000000000000000000000000000000;
  }
  else if (n==2) {
    x[ 0] = -0.577350269189625764509148780502;
    x[ 1] =  0.577350269189625764509148780502;
  }
  else if (n==3) {
    x[ 0] = -0.774596669241483377035853079956;
    x[ 1] =  0.000000000000000000000000000000;
    x[ 2] =  0.774596669241483377035853079956;
  }
  else if (n==4) {
    x[ 0] = -0.861136311594052575223946488893;
    x[ 1] = -0.339981043584856264802665759103;
    x[ 2] =  0.339981043584856264802665759103;
    x[ 3] =  0.861136311594052575223946488893;
  }
  else if (n==5) {
    x[ 0] = -0.906179845938663992797626878299;
    x[ 1] = -0.538469310105683091036314420700;
    x[ 2] =  0.000000000000000000000000000000;
    x[ 3] =  0.538469310105683091036314420700;
    x[ 4] =  0.906179845938663992797626878299;
  }
  else if (n==6) {
    x[ 0] = -0.932469514203152027812301554494;
    x[ 1] = -0.661209386466264513661399595020;
    x[ 2] = -0.238619186083196908630501721681;
    x[ 3] =  0.238619186083196908630501721681;
    x[ 4] =  0.661209386466264513661399595020;
    x[ 5] =  0.932469514203152027812301554494;
  }
  else if (n==7) {
    x[ 0] = -0.949107912342758524526189684048;
    x[ 1] = -0.741531185599394439863864773281;
    x[ 2] = -0.405845151377397166906606412077;
    x[ 3] =  0.000000000000000000000000000000;
    x[ 4] =  0.405845151377397166906606412077;
    x[ 5] =  0.741531185599394439863864773281;
    x[ 6] =  0.949107912342758524526189684048;
  }
  else if (n==8) {
    x[ 0] = -0.960289856497536231683560868569;
    x[ 1] = -0.796666477413626739591553936476;
    x[ 2] = -0.525532409916328985817739049189;
    x[ 3] = -0.183434642495649804939476142360;
    x[ 4] =  0.183434642495649804939476142360;
    x[ 5] =  0.525532409916328985817739049189;
    x[ 6] =  0.796666477413626739591553936476;
    x[ 7] =  0.960289856497536231683560868569;
  }
  else if (n==9) {
    x[ 0] = -0.968160239507626089835576203;
    x[ 1] = -0.836031107326635794299429788;
    x[ 2] = -0.613371432700590397308702039;
    x[ 3] = -0.324253423403808929038538015;
    x[ 4] =  0.000000000000000000000000000;
    x[ 5] =  0.324253423403808929038538015;
    x[ 6] =  0.613371432700590397308702039;
    x[ 7] =  0.836031107326635794299429788;
    x[ 8] =  0.968160239507626089835576203;
  }
  else if (n==10) {
    x[ 0] = -0.973906528517171720077964012;
    x[ 1] = -0.865063366688984510732096688;
    x[ 2] = -0.679409568299024406234327365;
    x[ 3] = -0.433395394129247190799265943;
    x[ 4] = -0.148874338981631210884826001;
    x[ 5] =  0.148874338981631210884826001;
    x[ 6] =  0.433395394129247190799265943;
    x[ 7] =  0.679409568299024406234327365;
    x[ 8] =  0.865063366688984510732096688;
    x[ 9] =  0.973906528517171720077964012;
  }
  else if (n==11) {
    x[ 0] = -0.978228658146056992803938001;
    x[ 1] = -0.887062599768095299075157769;
    x[ 2] = -0.730152005574049324093416252;
    x[ 3] = -0.519096129206811815925725669;
    x[ 4] = -0.269543155952344972331531985;
    x[ 5] =  0.000000000000000000000000000;
    x[ 6] =  0.269543155952344972331531985;
    x[ 7] =  0.519096129206811815925725669;
    x[ 8] =  0.730152005574049324093416252;
    x[ 9] =  0.887062599768095299075157769;
    x[10] =  0.978228658146056992803938001;
  }
  else if (n==12) {
    x[ 0] = -0.981560634246719250690549090;
    x[ 1] = -0.904117256370474856678465866;
    x[ 2] = -0.769902674194304687036893833;
    x[ 3] = -0.587317954286617447296702419;
    x[ 4] = -0.367831498998180193752691537;
    x[ 5] = -0.125233408511468915472441369;
    x[ 6] =  0.125233408511468915472441369;
    x[ 7] =  0.367831498998180193752691537;
    x[ 8] =  0.587317954286617447296702419;
    x[ 9] =  0.769902674194304687036893833;
    x[10] =  0.904117256370474856678465866;
    x[11] =  0.981560634246719250690549090;
  }
  else if (n==13) {
    x[ 0] = -0.984183054718588149472829449;
    x[ 1] = -0.917598399222977965206547837;
    x[ 2] = -0.801578090733309912794206490;
    x[ 3] = -0.642349339440340220643984607;
    x[ 4] = -0.448492751036446852877912852;
    x[ 5] = -0.230458315955134794065528121;
    x[ 6] =  0.000000000000000000000000000;
    x[ 7] =  0.230458315955134794065528121;
    x[ 8] =  0.448492751036446852877912852;
    x[ 9] =  0.642349339440340220643984607;
    x[10] =  0.80157809073330991279420649;
    x[11] =  0.91759839922297796520654784;
    x[12] =  0.98418305471858814947282945;
  }
  else if (n==14) {
    x[ 0] = -0.986283808696812338841597267;
    x[ 1] = -0.928434883663573517336391139;
    x[ 2] = -0.827201315069764993189794743;
    x[ 3] = -0.687292904811685470148019803;
    x[ 4] = -0.515248636358154091965290719;
    x[ 5] = -0.319112368927889760435671824;
    x[ 6] = -0.108054948707343662066244650;
    x[ 7] =  0.108054948707343662066244650;
    x[ 8] =  0.31911236892788976043567182;
    x[ 9] =  0.51524863635815409196529072;
    x[10] =  0.68729290481168547014801980;
    x[11] =  0.82720131506976499318979474;
    x[12] =  0.92843488366357351733639114;
    x[13] =  0.98628380869681233884159727;
  }
  else if (n==15) {
    x[ 0] = -0.987992518020485428489565719;
    x[ 1] = -0.937273392400705904307758948;
    x[ 2] = -0.848206583410427216200648321;
    x[ 3] = -0.724417731360170047416186055;
    x[ 4] = -0.570972172608538847537226737;
    x[ 5] = -0.394151347077563369897207371;
    x[ 6] = -0.201194093997434522300628303;
    x[ 7] =  0.00000000000000000000000000;
    x[ 8] =  0.20119409399743452230062830;
    x[ 9] =  0.39415134707756336989720737;
    x[10] =  0.57097217260853884753722674;
    x[11] =  0.72441773136017004741618605;
    x[12] =  0.84820658341042721620064832;
    x[13] =  0.93727339240070590430775895;
    x[14] =  0.98799251802048542848956572;
  }
  else if (n==16) {
    x[ 0] = -0.989400934991649932596154173;
    x[ 1] = -0.944575023073232576077988416;
    x[ 2] = -0.865631202387831743880467898;
    x[ 3] = -0.755404408355003033895101195;
    x[ 4] = -0.617876244402643748446671764;
    x[ 5] = -0.458016777657227386342419443;
    x[ 6] = -0.281603550779258913230460501;
    x[ 7] = -0.09501250983763744018531934;
    x[ 8] =  0.09501250983763744018531934;
    x[ 9] =  0.28160355077925891323046050;
    x[10] =  0.45801677765722738634241944;
    x[11] =  0.61787624440264374844667176;
    x[12] =  0.75540440835500303389510119;
    x[13] =  0.86563120238783174388046790;
    x[14] =  0.94457502307323257607798842;
    x[15] =  0.98940093499164993259615417;
  }
  else if (n==17) {
    x[ 0] = -0.990575475314417335675434020;
    x[ 1] = -0.950675521768767761222716958;
    x[ 2] = -0.880239153726985902122955694;
    x[ 3] = -0.781514003896801406925230056;
    x[ 4] = -0.657671159216690765850302217;
    x[ 5] = -0.512690537086476967886246569;
    x[ 6] = -0.35123176345387631529718552;
    x[ 7] = -0.17848418149584785585067749;
    x[ 8] =  0.00000000000000000000000000;
    x[ 9] =  0.17848418149584785585067749;
    x[10] =  0.35123176345387631529718552;
    x[11] =  0.51269053708647696788624657;
    x[12] =  0.65767115921669076585030222;
    x[13] =  0.78151400389680140692523006;
    x[14] =  0.88023915372698590212295569;
    x[15] =  0.95067552176876776122271696;
    x[16] =  0.99057547531441733567543402;
  }
  else if (n==18) {
    x[ 0] = -0.991565168420930946730016005;
    x[ 1] = -0.955823949571397755181195893;
    x[ 2] = -0.892602466497555739206060591;
    x[ 3] = -0.803704958972523115682417455;
    x[ 4] = -0.691687043060353207874891081;
    x[ 5] = -0.55977083107394753460787155;
    x[ 6] = -0.41175116146284264603593179;
    x[ 7] = -0.25188622569150550958897285;
    x[ 8] = -0.08477501304173530124226185;
    x[ 9] =  0.08477501304173530124226185;
    x[10] =  0.25188622569150550958897285;
    x[11] =  0.41175116146284264603593179;
    x[12] =  0.55977083107394753460787155;
    x[13] =  0.69168704306035320787489108;
    x[14] =  0.80370495897252311568241746;
    x[15] =  0.89260246649755573920606059;
    x[16] =  0.95582394957139775518119589;
    x[17] =  0.99156516842093094673001600;
  }
  else if (n==19) {
    x[ 0] = -0.992406843843584403189017670;
    x[ 1] = -0.960208152134830030852778841;
    x[ 2] = -0.903155903614817901642660929;
    x[ 3] = -0.822714656537142824978922487;
    x[ 4] = -0.72096617733522937861709586;
    x[ 5] = -0.60054530466168102346963816;
    x[ 6] = -0.46457074137596094571726715;
    x[ 7] = -0.31656409996362983199011733;
    x[ 8] = -0.16035864564022537586809612;
    x[ 9] =  0.00000000000000000000000000;
    x[10] =  0.16035864564022537586809612;
    x[11] =  0.31656409996362983199011733;
    x[12] =  0.46457074137596094571726715;
    x[13] =  0.60054530466168102346963816;
    x[14] =  0.72096617733522937861709586;
    x[15] =  0.82271465653714282497892249;
    x[16] =  0.90315590361481790164266093;
    x[17] =  0.96020815213483003085277884;
    x[18] =  0.99240684384358440318901767;
  }
  else if (n==20) {
    x[ 0] = -0.993128599185094924786122388;
    x[ 1] = -0.963971927277913791267666131;
    x[ 2] = -0.912234428251325905867752441;
    x[ 3] = -0.83911697182221882339452906;
    x[ 4] = -0.74633190646015079261430507;
    x[ 5] = -0.63605368072651502545283670;
    x[ 6] = -0.51086700195082709800436405;
    x[ 7] = -0.37370608871541956067254818;
    x[ 8] = -0.22778585114164507808049620;
    x[ 9] = -0.07652652113349733375464041;
    x[10] =  0.07652652113349733375464041;
    x[11] =  0.22778585114164507808049620;
    x[12] =  0.37370608871541956067254818;
    x[13] =  0.51086700195082709800436405;
    x[14] =  0.63605368072651502545283670;
    x[15] =  0.74633190646015079261430507;
    x[16] =  0.83911697182221882339452906;
    x[17] =  0.91223442825132590586775244;
    x[18] =  0.96397192727791379126766613;
    x[19] =  0.99312859918509492478612239;
  }
  else if (n==21) {
    x[ 0] =  -0.99375217062038950026024204;
    x[ 1] =  -0.96722683856630629431662221;
    x[ 2] =  -0.92009933415040082879018713;
    x[ 3] =  -0.85336336458331728364725064;
    x[ 4] =  -0.76843996347567790861587785;
    x[ 5] =  -0.66713880419741231930596667;
    x[ 6] =  -0.55161883588721980705901880;
    x[ 7] =  -0.42434212020743878357366889;
    x[ 8] =  -0.28802131680240109660079252;
    x[ 9] =  -0.14556185416089509093703098;
    x[10] =   0.00000000000000000000000000;
    x[11] =  +0.14556185416089509093703098;
    x[12] =  +0.28802131680240109660079252;
    x[13] =  +0.42434212020743878357366889;
    x[14] =  +0.55161883588721980705901880;
    x[15] =  +0.66713880419741231930596667;
    x[16] =  +0.76843996347567790861587785;
    x[17] =  +0.85336336458331728364725064;
    x[18] =  +0.92009933415040082879018713;
    x[19] =  +0.96722683856630629431662221;
    x[20] =  +0.99375217062038950026024204;
  }
  else if (n==22) {
    x[ 0] = -0.99429458548239929207303142;
    x[ 1] = -0.97006049783542872712395099;
    x[ 2] = -0.92695677218717400052069294;
    x[ 3] = -0.86581257772030013653642564;
    x[ 4] = -0.78781680597920816200427796;
    x[ 5] = -0.69448726318668278005068984;
    x[ 6] = -0.58764040350691159295887693;
    x[ 7] = -0.46935583798675702640633071;
    x[ 8] = -0.34193582089208422515814742;
    x[ 9] = -0.20786042668822128547884653;
    x[10] = -0.06973927331972222121384180;
    x[11] =  0.06973927331972222121384180;
    x[12] =  0.20786042668822128547884653;
    x[13] =  0.34193582089208422515814742;
    x[14] =  0.46935583798675702640633071;
    x[15] =  0.58764040350691159295887693;
    x[16] =  0.69448726318668278005068984;
    x[17] =  0.78781680597920816200427796;
    x[18] =  0.86581257772030013653642564;
    x[19] =  0.92695677218717400052069294;
    x[20] =  0.97006049783542872712395099;
    x[21] =  0.99429458548239929207303142;
  }
  else if (n==23) {
    x[ 0] = -0.99476933499755212352392572;
    x[ 1] = -0.97254247121811523195602408;
    x[ 2] = -0.93297108682601610234919699;
    x[ 3] = -0.87675235827044166737815689;
    x[ 4] = -0.80488840161883989215111841;
    x[ 5] = -0.71866136313195019446162448;
    x[ 6] = -0.61960987576364615638509731;
    x[ 7] = -0.50950147784600754968979305;
    x[ 8] = -0.39030103803029083142148887;
    x[ 9] = -0.26413568097034493053386954;
    x[10] = -0.13325682429846611093174268;
    x[11] =  0.00000000000000000000000000;
    x[12] =  0.13325682429846611093174268;
    x[13] =  0.26413568097034493053386954;
    x[14] =  0.39030103803029083142148887;
    x[15] =  0.50950147784600754968979305;
    x[16] =  0.61960987576364615638509731;
    x[17] =  0.71866136313195019446162448;
    x[18] =  0.80488840161883989215111841;
    x[19] =  0.87675235827044166737815689;
    x[20] =  0.93297108682601610234919699;
    x[21] =  0.97254247121811523195602408;
    x[22] =  0.99476933499755212352392572;
  }
  else if (n==24) {
    x[ 0] = -0.99518721999702136017999741;
    x[ 1] = -0.97472855597130949819839199;
    x[ 2] = -0.93827455200273275852364900;
    x[ 3] = -0.88641552700440103421315434;
    x[ 4] = -0.82000198597390292195394987;
    x[ 5] = -0.74012419157855436424382810;
    x[ 6] = -0.64809365193697556925249579;
    x[ 7] = -0.54542147138883953565837562;
    x[ 8] = -0.43379350762604513848708423;
    x[ 9] = -0.31504267969616337438679329;
    x[10] = -0.19111886747361630915863982;
    x[11] = -0.06405689286260562608504308;
    x[12] =  0.06405689286260562608504308;
    x[13] =  0.19111886747361630915863982;
    x[14] =  0.31504267969616337438679329;
    x[15] =  0.43379350762604513848708423;
    x[16] =  0.54542147138883953565837562;
    x[17] =  0.64809365193697556925249579;
    x[18] =  0.74012419157855436424382810;
    x[19] =  0.82000198597390292195394987;
    x[20] =  0.88641552700440103421315434;
    x[21] =  0.93827455200273275852364900;
    x[22] =  0.97472855597130949819839199;
    x[23] =  0.99518721999702136017999741;
  }
  else if (n==25) {
    x[ 0] = -0.99555696979049809790878495;
    x[ 1] = -0.97666392145951751149831539;
    x[ 2] = -0.94297457122897433941401117;
    x[ 3] = -0.89499199787827536885104201;
    x[ 4] = -0.83344262876083400142102111;
    x[ 5] = -0.75925926303735763057728287;
    x[ 6] = -0.67356636847346836448512063;
    x[ 7] = -0.57766293024122296772368984;
    x[ 8] = -0.47300273144571496052218212;
    x[ 9] = -0.36117230580938783773582173;
    x[10] = -0.24386688372098843204519036;
    x[11] = -0.12286469261071039638735982;
    x[12] =  0.00000000000000000000000000;
    x[13] =  0.12286469261071039638735982;
    x[14] =  0.24386688372098843204519036;
    x[15] =  0.36117230580938783773582173;
    x[16] =  0.47300273144571496052218212;
    x[17] =  0.57766293024122296772368984;
    x[18] =  0.67356636847346836448512063;
    x[19] =  0.75925926303735763057728287;
    x[20] =  0.83344262876083400142102111;
    x[21] =  0.89499199787827536885104201;
    x[22] =  0.94297457122897433941401117;
    x[23] =  0.97666392145951751149831539;
    x[24] =  0.99555696979049809790878495;
  }
  else if (n==26) {
    x[ 0] = -0.99588570114561692900321696;
    x[ 1] = -0.97838544595647099110058035;
    x[ 2] = -0.94715906666171425013591528;
    x[ 3] = -0.90263786198430707421766560;
    x[ 4] = -0.84544594278849801879750706;
    x[ 5] = -0.77638594882067885619296725;
    x[ 6] = -0.69642726041995726486381391;
    x[ 7] = -0.60669229301761806323197875;
    x[ 8] = -0.50844071482450571769570306;
    x[ 9] = -0.40305175512348630648107738;
    x[10] = -0.29200483948595689514283538;
    x[11] = -0.17685882035689018396905775;
    x[12] = -0.05923009342931320709371858;
    x[13] =  0.05923009342931320709371858;
    x[14] =  0.17685882035689018396905775;
    x[15] =  0.29200483948595689514283538;
    x[16] =  0.40305175512348630648107738;
    x[17] =  0.50844071482450571769570306;
    x[18] =  0.60669229301761806323197875;
    x[19] =  0.69642726041995726486381391;
    x[20] =  0.77638594882067885619296725;
    x[21] =  0.84544594278849801879750706;
    x[22] =  0.90263786198430707421766560;
    x[23] =  0.94715906666171425013591528;
    x[24] =  0.97838544595647099110058035;
    x[25] =  0.99588570114561692900321696;
  }
  else if (n==27) {
    x[ 0] = -0.99617926288898856693888721;
    x[ 1] = -0.97992347596150122285587336;
    x[ 2] = -0.95090055781470500685190803;
    x[ 3] = -0.90948232067749110430064502;
    x[ 4] = -0.85620790801829449030273722;
    x[ 5] = -0.79177163907050822714439734;
    x[ 6] = -0.71701347373942369929481621;
    x[ 7] = -0.63290797194649514092773464;
    x[ 8] = -0.54055156457945689490030094;
    x[ 9] = -0.44114825175002688058597416;
    x[10] = -0.33599390363850889973031903;
    x[11] = -0.22645936543953685885723911;
    x[12] = -0.11397258560952996693289498;
    x[13] =  0.00000000000000000000000000;
    x[14] =  0.11397258560952996693289498;
    x[15] =  0.22645936543953685885723911;
    x[16] =  0.33599390363850889973031903;
    x[17] =  0.44114825175002688058597416;
    x[18] =  0.54055156457945689490030094;
    x[19] =  0.63290797194649514092773464;
    x[20] =  0.71701347373942369929481621;
    x[21] =  0.79177163907050822714439734;
    x[22] =  0.85620790801829449030273722;
    x[23] =  0.90948232067749110430064502;
    x[24] =  0.95090055781470500685190803;
    x[25] =  0.97992347596150122285587336;
    x[26] =  0.99617926288898856693888721;
  }
  else if (n==28) {
    x[ 0] = -0.99644249757395444995043639;
    x[ 1] = -0.98130316537087275369455995;
    x[ 2] = -0.95425928062893819725410184;
    x[ 3] = -0.91563302639213207386968942;
    x[ 4] = -0.86589252257439504894225457;
    x[ 5] = -0.80564137091717917144788596;
    x[ 6] = -0.73561087801363177202814451;
    x[ 7] = -0.65665109403886496121989818;
    x[ 8] = -0.56972047181140171930800328;
    x[ 9] = -0.47587422495511826103441185;
    x[10] = -0.37625151608907871022135721;
    x[11] = -0.27206162763517807767682636;
    x[12] = -0.16456928213338077128147178;
    x[13] = -0.05507928988403427042651653;
    x[14] =  0.05507928988403427042651653;
    x[15] =  0.16456928213338077128147178;
    x[16] =  0.27206162763517807767682636;
    x[17] =  0.37625151608907871022135721;
    x[18] =  0.47587422495511826103441185;
    x[19] =  0.56972047181140171930800328;
    x[20] =  0.65665109403886496121989818;
    x[21] =  0.73561087801363177202814451;
    x[22] =  0.80564137091717917144788596;
    x[23] =  0.86589252257439504894225457;
    x[24] =  0.91563302639213207386968942;
    x[25] =  0.95425928062893819725410184;
    x[26] =  0.98130316537087275369455995;
    x[27] =  0.99644249757395444995043639;
  }
  else if (n==29) {
    x[ 0] = -0.99667944226059658616319153;
    x[ 1] = -0.98254550526141317487092602;
    x[ 2] = -0.95728559577808772579820804;
    x[ 3] = -0.92118023295305878509375344;
    x[ 4] = -0.87463780492010279041779342;
    x[ 5] = -0.81818548761525244498957221;
    x[ 6] = -0.75246285173447713391261008;
    x[ 7] = -0.67821453760268651515618501;
    x[ 8] = -0.59628179713822782037958621;
    x[ 9] = -0.50759295512422764210262792;
    x[10] = -0.41315288817400866389070659;
    x[11] = -0.31403163786763993494819592;
    x[12] = -0.21135228616600107450637573;
    x[13] = -0.10627823013267923017098239;
    x[14] =  0.00000000000000000000000000;
    x[15] =  0.10627823013267923017098239;
    x[16] =  0.21135228616600107450637573;
    x[17] =  0.31403163786763993494819592;
    x[18] =  0.41315288817400866389070659;
    x[19] =  0.50759295512422764210262792;
    x[20] =  0.59628179713822782037958621;
    x[21] =  0.67821453760268651515618501;
    x[22] =  0.75246285173447713391261008;
    x[23] =  0.81818548761525244498957221;
    x[24] =  0.87463780492010279041779342;
    x[25] =  0.92118023295305878509375344;
    x[26] =  0.95728559577808772579820804;
    x[27] =  0.98254550526141317487092602;
    x[28] =  0.99667944226059658616319153;
  }
  else if (n==30) {
    x[ 0] = -0.99689348407464954027163005;
    x[ 1] = -0.98366812327974720997003258;
    x[ 2] = -0.96002186496830751221687103;
    x[ 3] = -0.92620004742927432587932428;
    x[ 4] = -0.88256053579205268154311646;
    x[ 5] = -0.82956576238276839744289812;
    x[ 6] = -0.76777743210482619491797734;
    x[ 7] = -0.69785049479331579693229239;
    x[ 8] = -0.62052618298924286114047756;
    x[ 9] = -0.53662414814201989926416979;
    x[10] = -0.44703376953808917678060990;
    x[11] = -0.35270472553087811347103721;
    x[12] = -0.25463692616788984643980513;
    x[13] = -0.15386991360858354696379467;
    x[14] = -0.05147184255531769583302521;
    x[15] =  0.05147184255531769583302521;
    x[16] =  0.15386991360858354696379467;
    x[17] =  0.25463692616788984643980513;
    x[18] =  0.35270472553087811347103721;
    x[19] =  0.44703376953808917678060990;
    x[20] =  0.53662414814201989926416979;
    x[21] =  0.62052618298924286114047756;
    x[22] =  0.69785049479331579693229239;
    x[23] =  0.76777743210482619491797734;
    x[24] =  0.82956576238276839744289812;
    x[25] =  0.88256053579205268154311646;
    x[26] =  0.92620004742927432587932428;
    x[27] =  0.96002186496830751221687103;
    x[28] =  0.98366812327974720997003258;
    x[29] =  0.99689348407464954027163005;
  }
  else if (n==31) {
    x[ 0] = -0.99708748181947707405562655;
    x[ 1] = -0.98468590966515248400246517;
    x[ 2] = -0.96250392509294966178905240;
    x[ 3] = -0.93075699789664816495694576;
    x[ 4] = -0.88976002994827104337419201;
    x[ 5] = -0.83992032014626734008690454;
    x[ 6] = -0.78173314841662494040636002;
    x[ 7] = -0.71577678458685328390597087;
    x[ 8] = -0.64270672292426034618441820;
    x[ 9] = -0.56324916140714926272094492;
    x[10] = -0.47819378204490248044059404;
    x[11] = -0.38838590160823294306135146;
    x[12] = -0.29471806998170161661790390;
    x[13] = -0.19812119933557062877241300;
    x[14] = -0.09955531215234152032517479;
    x[15] =  0.00000000000000000000000000;
    x[16] =  0.09955531215234152032517479;
    x[17] =  0.19812119933557062877241300;
    x[18] =  0.29471806998170161661790390;
    x[19] =  0.38838590160823294306135146;
    x[20] =  0.47819378204490248044059404;
    x[21] =  0.56324916140714926272094492;
    x[22] =  0.64270672292426034618441820;
    x[23] =  0.71577678458685328390597087;
    x[24] =  0.78173314841662494040636002;
    x[25] =  0.83992032014626734008690454;
    x[26] =  0.88976002994827104337419201;
    x[27] =  0.93075699789664816495694576;
    x[28] =  0.96250392509294966178905240;
    x[29] =  0.98468590966515248400246517;
    x[30] =  0.99708748181947707405562655;
  }
  else if (n==32) {
    x[ 0] = -0.99726386184948156354498113;
    x[ 1] = -0.98561151154526833540017504;
    x[ 2] = -0.96476225558750643077381193;
    x[ 3] = -0.93490607593773968917091913;
    x[ 4] = -0.89632115576605212396530724;
    x[ 5] = -0.84936761373256997013369300;
    x[ 6] = -0.79448379596794240696309730;
    x[ 7] = -0.73218211874028968038742667;
    x[ 8] = -0.66304426693021520097511517;
    x[ 9] = -0.58771575724076232904074548;
    x[10] = -0.50689990893222939002374747;
    x[11] = -0.42135127613063534536411944;
    x[12] = -0.33186860228212764977991681;
    x[13] = -0.23928736225213707454460321;
    x[14] = -0.14447196158279649348518637;
    x[15] = -0.04830766568773831623481257;
    x[16] =  0.04830766568773831623481257;
    x[17] =  0.14447196158279649348518637;
    x[18] =  0.23928736225213707454460321;
    x[19] =  0.33186860228212764977991681;
    x[20] =  0.42135127613063534536411944;
    x[21] =  0.50689990893222939002374747;
    x[22] =  0.58771575724076232904074548;
    x[23] =  0.66304426693021520097511517;
    x[24] =  0.73218211874028968038742667;
    x[25] =  0.79448379596794240696309730;
    x[26] =  0.84936761373256997013369300;
    x[27] =  0.89632115576605212396530724;
    x[28] =  0.93490607593773968917091913;
    x[29] =  0.96476225558750643077381193;
    x[30] =  0.98561151154526833540017504;
    x[31] =  0.99726386184948156354498113;
  }
  else if (n==33) {
    x[ 0] = -0.99742469424645521726616802;
    x[ 1] = -0.98645572623064248811037570;
    x[ 2] = -0.96682290968999276892837771;
    x[ 3] = -0.93869437261116835035583512;
    x[ 4] = -0.90231676774343358304053133;
    x[ 5] = -0.85800965267650406464306148;
    x[ 6] = -0.80616235627416658979620087;
    x[ 7] = -0.74723049644956215785905512;
    x[ 8] = -0.68173195996974278626821595;
    x[ 9] = -0.61024234583637902730728751;
    x[10] = -0.53338990478634764354889426;
    x[11] = -0.45185001727245069572599328;
    x[12] = -0.36633925774807334107022062;
    x[13] = -0.27760909715249702940324807;
    x[14] = -0.18643929882799157233579876;
    x[15] = -0.09363106585473338567074292;
    x[16] =  0.00000000000000000000000000;
    x[17] =  0.09363106585473338567074292;
    x[18] =  0.18643929882799157233579876;
    x[19] =  0.27760909715249702940324807;
    x[20] =  0.36633925774807334107022062;
    x[21] =  0.45185001727245069572599328;
    x[22] =  0.53338990478634764354889426;
    x[23] =  0.61024234583637902730728751;
    x[24] =  0.68173195996974278626821595;
    x[25] =  0.74723049644956215785905512;
    x[26] =  0.80616235627416658979620087;
    x[27] =  0.85800965267650406464306148;
    x[28] =  0.90231676774343358304053133;
    x[29] =  0.93869437261116835035583512;
    x[30] =  0.96682290968999276892837771;
    x[31] =  0.98645572623064248811037570;
    x[32] =  0.99742469424645521726616802;
  }
  else {
    std::cerr << "\n";
    std::cerr << "LEGENDRE_LOOKUP_POINTS - Fatal error!\n";
    std::cerr << "  Illegal value of N = " << n << "\n";
    std::cerr << "  Legal values are 1 through 33.\n";
    std::exit(1);
  }
  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::legendre_lookup_weights ( int n, Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    LEGENDRE_LOOKUP_WEIGHTS looks up weights for Gauss-Legendre quadrature.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Vladimir Krylov,
//    Approximate Calculation of Integrals,
//    Dover, 2006,
//    ISBN: 0486445798.
//    LC: QA311.K713.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3,
//    LC: QA47.M315.
//
//  Parameters:
//
//    Input, int N, the order.
//    N must be between 1 and 33.
//
//    Output, Scalar W[N], the weights.
//
{
  if (n==1) {
    w[ 0] = 2.000000000000000000000000000000;
  }
  else if (n==2) {
    w[ 0] = 1.000000000000000000000000000000;
    w[ 1] = 1.000000000000000000000000000000;
  }
  else if (n==3) {
    w[ 0] = 0.555555555555555555555555555556;
    w[ 1] = 0.888888888888888888888888888889;
    w[ 2] = 0.555555555555555555555555555556;
  }
  else if (n==4) {
    w[ 0] = 0.347854845137453857373063949222;
    w[ 1] = 0.652145154862546142626936050778;
    w[ 2] = 0.652145154862546142626936050778;
    w[ 3] = 0.347854845137453857373063949222;
  }
  else if (n==5) {
    w[ 0] = 0.236926885056189087514264040720;
    w[ 1] = 0.478628670499366468041291514836;
    w[ 2] = 0.568888888888888888888888888889;
    w[ 3] = 0.478628670499366468041291514836;
    w[ 4] = 0.236926885056189087514264040720;
  }
  else if (n==6) {
    w[ 0] = 0.171324492379170345040296142173;
    w[ 1] = 0.360761573048138607569833513838;
    w[ 2] = 0.467913934572691047389870343990;
    w[ 3] = 0.467913934572691047389870343990;
    w[ 4] = 0.360761573048138607569833513838;
    w[ 5] = 0.171324492379170345040296142173;
  }
  else if (n==7) {
    w[ 0] = 0.129484966168869693270611432679;
    w[ 1] = 0.279705391489276667901467771424;
    w[ 2] = 0.381830050505118944950369775489;
    w[ 3] = 0.417959183673469387755102040816;
    w[ 4] = 0.381830050505118944950369775489;
    w[ 5] = 0.279705391489276667901467771424;
    w[ 6] = 0.129484966168869693270611432679;
  }
  else if (n==8) {
    w[ 0] = 0.101228536290376259152531354310;
    w[ 1] = 0.222381034453374470544355994426;
    w[ 2] = 0.313706645877887287337962201987;
    w[ 3] = 0.362683783378361982965150449277;
    w[ 4] = 0.362683783378361982965150449277;
    w[ 5] = 0.313706645877887287337962201987;
    w[ 6] = 0.222381034453374470544355994426;
    w[ 7] = 0.101228536290376259152531354310;
  }
  else if (n==9) {
    w[ 0] = 0.081274388361574411971892158111;
    w[ 1] = 0.18064816069485740405847203124;
    w[ 2] = 0.26061069640293546231874286942;
    w[ 3] = 0.31234707704000284006863040658;
    w[ 4] = 0.33023935500125976316452506929;
    w[ 5] = 0.31234707704000284006863040658;
    w[ 6] = 0.26061069640293546231874286942;
    w[ 7] = 0.18064816069485740405847203124;
    w[ 8] = 0.081274388361574411971892158111;
  }
  else if (n==10) {
    w[ 0] = 0.066671344308688137593568809893;
    w[ 1] = 0.14945134915058059314577633966;
    w[ 2] = 0.21908636251598204399553493423;
    w[ 3] = 0.26926671930999635509122692157;
    w[ 4] = 0.29552422471475287017389299465;
    w[ 5] = 0.29552422471475287017389299465;
    w[ 6] = 0.26926671930999635509122692157;
    w[ 7] = 0.21908636251598204399553493423;
    w[ 8] = 0.14945134915058059314577633966;
    w[ 9] = 0.066671344308688137593568809893;
  }
  else if (n==11) {
    w[ 0] = 0.055668567116173666482753720443;
    w[ 1] = 0.12558036946490462463469429922;
    w[ 2] = 0.18629021092773425142609764143;
    w[ 3] = 0.23319376459199047991852370484;
    w[ 4] = 0.26280454451024666218068886989;
    w[ 5] = 0.27292508677790063071448352834;
    w[ 6] = 0.26280454451024666218068886989;
    w[ 7] = 0.23319376459199047991852370484;
    w[ 8] = 0.18629021092773425142609764143;
    w[ 9] = 0.12558036946490462463469429922;
    w[10] = 0.055668567116173666482753720443;
  }
  else if (n==12) {
    w[ 0] = 0.047175336386511827194615961485;
    w[ 1] = 0.10693932599531843096025471819;
    w[ 2] = 0.16007832854334622633465252954;
    w[ 3] = 0.20316742672306592174906445581;
    w[ 4] = 0.23349253653835480876084989892;
    w[ 5] = 0.24914704581340278500056243604;
    w[ 6] = 0.24914704581340278500056243604;
    w[ 7] = 0.23349253653835480876084989892;
    w[ 8] = 0.20316742672306592174906445581;
    w[ 9] = 0.16007832854334622633465252954;
    w[10] = 0.10693932599531843096025471819;
    w[11] = 0.047175336386511827194615961485;
  }
  else if (n==13) {
    w[ 0] = 0.040484004765315879520021592201;
    w[ 1] = 0.092121499837728447914421775954;
    w[ 2] = 0.13887351021978723846360177687;
    w[ 3] = 0.17814598076194573828004669200;
    w[ 4] = 0.20781604753688850231252321931;
    w[ 5] = 0.22628318026289723841209018604;
    w[ 6] = 0.23255155323087391019458951527;
    w[ 7] = 0.22628318026289723841209018604;
    w[ 8] = 0.20781604753688850231252321931;
    w[ 9] = 0.17814598076194573828004669200;
    w[10] = 0.13887351021978723846360177687;
    w[11] = 0.092121499837728447914421775954;
    w[12] = 0.040484004765315879520021592201;
  }
  else if (n==14) {
    w[ 0] = 0.035119460331751863031832876138;
    w[ 1] = 0.08015808715976020980563327706;
    w[ 2] = 0.12151857068790318468941480907;
    w[ 3] = 0.15720316715819353456960193862;
    w[ 4] = 0.18553839747793781374171659013;
    w[ 5] = 0.20519846372129560396592406566;
    w[ 6] = 0.21526385346315779019587644332;
    w[ 7] = 0.21526385346315779019587644332;
    w[ 8] = 0.20519846372129560396592406566;
    w[ 9] = 0.18553839747793781374171659013;
    w[10] = 0.15720316715819353456960193862;
    w[11] = 0.12151857068790318468941480907;
    w[12] = 0.08015808715976020980563327706;
    w[13] = 0.035119460331751863031832876138;
  }
  else if (n==15) {
    w[ 0] = 0.030753241996117268354628393577;
    w[ 1] = 0.070366047488108124709267416451;
    w[ 2] = 0.107159220467171935011869546686;
    w[ 3] = 0.13957067792615431444780479451;
    w[ 4] = 0.16626920581699393355320086048;
    w[ 5] = 0.18616100001556221102680056187;
    w[ 6] = 0.19843148532711157645611832644;
    w[ 7] = 0.20257824192556127288062019997;
    w[ 8] = 0.19843148532711157645611832644;
    w[ 9] = 0.18616100001556221102680056187;
    w[10] = 0.16626920581699393355320086048;
    w[11] = 0.13957067792615431444780479451;
    w[12] = 0.107159220467171935011869546686;
    w[13] = 0.070366047488108124709267416451;
    w[14] = 0.030753241996117268354628393577;
  }
  else if (n==16) {
    w[ 0] = 0.027152459411754094851780572456;
    w[ 1] = 0.062253523938647892862843836994;
    w[ 2] = 0.09515851168249278480992510760;
    w[ 3] = 0.12462897125553387205247628219;
    w[ 4] = 0.14959598881657673208150173055;
    w[ 5] = 0.16915651939500253818931207903;
    w[ 6] = 0.18260341504492358886676366797;
    w[ 7] = 0.18945061045506849628539672321;
    w[ 8] = 0.18945061045506849628539672321;
    w[ 9] = 0.18260341504492358886676366797;
    w[10] = 0.16915651939500253818931207903;
    w[11] = 0.14959598881657673208150173055;
    w[12] = 0.12462897125553387205247628219;
    w[13] = 0.09515851168249278480992510760;
    w[14] = 0.062253523938647892862843836994;
    w[15] = 0.027152459411754094851780572456;
  }
  else if (n==17) {
    w[ 0] = 0.024148302868547931960110026288;
    w[ 1] = 0.055459529373987201129440165359;
    w[ 2] = 0.085036148317179180883535370191;
    w[ 3] = 0.111883847193403971094788385626;
    w[ 4] = 0.13513636846852547328631998170;
    w[ 5] = 0.15404576107681028808143159480;
    w[ 6] = 0.16800410215645004450997066379;
    w[ 7] = 0.17656270536699264632527099011;
    w[ 8] = 0.17944647035620652545826564426;
    w[ 9] = 0.17656270536699264632527099011;
    w[10] = 0.16800410215645004450997066379;
    w[11] = 0.15404576107681028808143159480;
    w[12] = 0.13513636846852547328631998170;
    w[13] = 0.111883847193403971094788385626;
    w[14] = 0.085036148317179180883535370191;
    w[15] = 0.055459529373987201129440165359;
    w[16] = 0.024148302868547931960110026288;
  }
  else if (n==18) {
    w[ 0] = 0.021616013526483310313342710266;
    w[ 1] = 0.049714548894969796453334946203;
    w[ 2] = 0.07642573025488905652912967762;
    w[ 3] = 0.10094204410628716556281398492;
    w[ 4] = 0.12255520671147846018451912680;
    w[ 5] = 0.14064291467065065120473130375;
    w[ 6] = 0.15468467512626524492541800384;
    w[ 7] = 0.16427648374583272298605377647;
    w[ 8] = 0.16914238296314359184065647013;
    w[ 9] = 0.16914238296314359184065647013;
    w[10] = 0.16427648374583272298605377647;
    w[11] = 0.15468467512626524492541800384;
    w[12] = 0.14064291467065065120473130375;
    w[13] = 0.12255520671147846018451912680;
    w[14] = 0.10094204410628716556281398492;
    w[15] = 0.07642573025488905652912967762;
    w[16] = 0.049714548894969796453334946203;
    w[17] = 0.021616013526483310313342710266;
  }
  else if (n==19) {
    w[ 0] = 0.019461788229726477036312041464;
    w[ 1] = 0.044814226765699600332838157402;
    w[ 2] = 0.069044542737641226580708258006;
    w[ 3] = 0.091490021622449999464462094124;
    w[ 4] = 0.111566645547333994716023901682;
    w[ 5] = 0.12875396253933622767551578486;
    w[ 6] = 0.14260670217360661177574610944;
    w[ 7] = 0.15276604206585966677885540090;
    w[ 8] = 0.15896884339395434764995643946;
    w[ 9] = 0.16105444984878369597916362532;
    w[10] = 0.15896884339395434764995643946;
    w[11] = 0.15276604206585966677885540090;
    w[12] = 0.14260670217360661177574610944;
    w[13] = 0.12875396253933622767551578486;
    w[14] = 0.111566645547333994716023901682;
    w[15] = 0.091490021622449999464462094124;
    w[16] = 0.069044542737641226580708258006;
    w[17] = 0.044814226765699600332838157402;
    w[18] = 0.019461788229726477036312041464;
  }
  else if (n==20) {
    w[ 0] = 0.017614007139152118311861962352;
    w[ 1] = 0.040601429800386941331039952275;
    w[ 2] = 0.062672048334109063569506535187;
    w[ 3] = 0.08327674157670474872475814322;
    w[ 4] = 0.10193011981724043503675013548;
    w[ 5] = 0.11819453196151841731237737771;
    w[ 6] = 0.13168863844917662689849449975;
    w[ 7] = 0.14209610931838205132929832507;
    w[ 8] = 0.14917298647260374678782873700;
    w[ 9] = 0.15275338713072585069808433195;
    w[10] = 0.15275338713072585069808433195;
    w[11] = 0.14917298647260374678782873700;
    w[12] = 0.14209610931838205132929832507;
    w[13] = 0.13168863844917662689849449975;
    w[14] = 0.11819453196151841731237737771;
    w[15] = 0.10193011981724043503675013548;
    w[16] = 0.08327674157670474872475814322;
    w[17] = 0.062672048334109063569506535187;
    w[18] = 0.040601429800386941331039952275;
    w[19] = 0.017614007139152118311861962352;
  }
  else if (n==21) {
    w[ 0] =   0.016017228257774333324224616858;
    w[ 1] =   0.036953789770852493799950668299;
    w[ 2] =   0.057134425426857208283635826472;
    w[ 3] =   0.076100113628379302017051653300;
    w[ 4] =   0.093444423456033861553289741114;
    w[ 5] =   0.108797299167148377663474578070;
    w[ 6] =   0.12183141605372853419536717713;
    w[ 7] =   0.13226893863333746178105257450;
    w[ 8] =   0.13988739479107315472213342387;
    w[ 9] =   0.14452440398997005906382716655;
    w[10] =   0.14608113364969042719198514768;
    w[11] =   0.14452440398997005906382716655;
    w[12] =   0.13988739479107315472213342387;
    w[13] =   0.13226893863333746178105257450;
    w[14] =   0.12183141605372853419536717713;
    w[15] =   0.108797299167148377663474578070;
    w[16] =   0.093444423456033861553289741114;
    w[17] =   0.076100113628379302017051653300;
    w[18] =   0.057134425426857208283635826472;
    w[19] =   0.036953789770852493799950668299;
    w[20] =   0.016017228257774333324224616858;
  }
  else if (n==22) {
    w[ 0] = 0.014627995298272200684991098047;
    w[ 1] = 0.033774901584814154793302246866;
    w[ 2] = 0.052293335152683285940312051273;
    w[ 3] = 0.06979646842452048809496141893;
    w[ 4] = 0.08594160621706772741444368137;
    w[ 5] = 0.10041414444288096493207883783;
    w[ 6] = 0.11293229608053921839340060742;
    w[ 7] = 0.12325237681051242428556098615;
    w[ 8] = 0.13117350478706237073296499253;
    w[ 9] = 0.13654149834601517135257383123;
    w[10] = 0.13925187285563199337541024834;
    w[11] = 0.13925187285563199337541024834;
    w[12] = 0.13654149834601517135257383123;
    w[13] = 0.13117350478706237073296499253;
    w[14] = 0.12325237681051242428556098615;
    w[15] = 0.11293229608053921839340060742;
    w[16] = 0.10041414444288096493207883783;
    w[17] = 0.08594160621706772741444368137;
    w[18] = 0.06979646842452048809496141893;
    w[19] = 0.052293335152683285940312051273;
    w[20] = 0.033774901584814154793302246866;
    w[21] = 0.014627995298272200684991098047;
  }
  else if (n==23) {
    w[ 0] = 0.013411859487141772081309493459;
    w[ 1] = 0.030988005856979444310694219642;
    w[ 2] = 0.048037671731084668571641071632;
    w[ 3] = 0.064232421408525852127169615159;
    w[ 4] = 0.079281411776718954922892524742;
    w[ 5] = 0.092915766060035147477018617370;
    w[ 6] = 0.104892091464541410074086185015;
    w[ 7] = 0.11499664022241136494164351293;
    w[ 8] = 0.12304908430672953046757840067;
    w[ 9] = 0.12890572218808214997859533940;
    w[10] = 0.13246203940469661737164246470;
    w[11] = 0.13365457218610617535145711055;
    w[12] = 0.13246203940469661737164246470;
    w[13] = 0.12890572218808214997859533940;
    w[14] = 0.12304908430672953046757840067;
    w[15] = 0.11499664022241136494164351293;
    w[16] = 0.104892091464541410074086185015;
    w[17] = 0.092915766060035147477018617370;
    w[18] = 0.079281411776718954922892524742;
    w[19] = 0.064232421408525852127169615159;
    w[20] = 0.048037671731084668571641071632;
    w[21] = 0.030988005856979444310694219642;
    w[22] = 0.013411859487141772081309493459;
  }
  else if (n==24) {
    w[ 0] = 0.012341229799987199546805667070;
    w[ 1] = 0.028531388628933663181307815952;
    w[ 2] = 0.044277438817419806168602748211;
    w[ 3] = 0.059298584915436780746367758500;
    w[ 4] = 0.07334648141108030573403361525;
    w[ 5] = 0.08619016153195327591718520298;
    w[ 6] = 0.09761865210411388826988066446;
    w[ 7] = 0.10744427011596563478257734245;
    w[ 8] = 0.11550566805372560135334448391;
    w[ 9] = 0.12167047292780339120446315348;
    w[10] = 0.12583745634682829612137538251;
    w[11] = 0.12793819534675215697405616522;
    w[12] = 0.12793819534675215697405616522;
    w[13] = 0.12583745634682829612137538251;
    w[14] = 0.12167047292780339120446315348;
    w[15] = 0.11550566805372560135334448391;
    w[16] = 0.10744427011596563478257734245;
    w[17] = 0.09761865210411388826988066446;
    w[18] = 0.08619016153195327591718520298;
    w[19] = 0.07334648141108030573403361525;
    w[20] = 0.059298584915436780746367758500;
    w[21] = 0.044277438817419806168602748211;
    w[22] = 0.028531388628933663181307815952;
    w[23] = 0.012341229799987199546805667070;
  }
  else if (n==25) {
    w[ 0] = 0.0113937985010262879479029641132;
    w[ 1] = 0.026354986615032137261901815295;
    w[ 2] = 0.040939156701306312655623487712;
    w[ 3] = 0.054904695975835191925936891541;
    w[ 4] = 0.068038333812356917207187185657;
    w[ 5] = 0.080140700335001018013234959669;
    w[ 6] = 0.091028261982963649811497220703;
    w[ 7] = 0.100535949067050644202206890393;
    w[ 8] = 0.108519624474263653116093957050;
    w[ 9] = 0.11485825914571164833932554587;
    w[10] = 0.11945576353578477222817812651;
    w[11] = 0.12224244299031004168895951895;
    w[12] = 0.12317605372671545120390287308;
    w[13] = 0.12224244299031004168895951895;
    w[14] = 0.11945576353578477222817812651;
    w[15] = 0.11485825914571164833932554587;
    w[16] = 0.108519624474263653116093957050;
    w[17] = 0.100535949067050644202206890393;
    w[18] = 0.091028261982963649811497220703;
    w[19] = 0.080140700335001018013234959669;
    w[20] = 0.068038333812356917207187185657;
    w[21] = 0.054904695975835191925936891541;
    w[22] = 0.040939156701306312655623487712;
    w[23] = 0.026354986615032137261901815295;
    w[24] = 0.0113937985010262879479029641132;
  }
  else if (n==26) {
    w[ 0] = 0.010551372617343007155651187685;
    w[ 1] = 0.024417851092631908789615827520;
    w[ 2] = 0.037962383294362763950303141249;
    w[ 3] = 0.050975825297147811998319900724;
    w[ 4] = 0.063274046329574835539453689907;
    w[ 5] = 0.07468414976565974588707579610;
    w[ 6] = 0.08504589431348523921044776508;
    w[ 7] = 0.09421380035591414846366488307;
    w[ 8] = 0.10205916109442542323841407025;
    w[ 9] = 0.10847184052857659065657942673;
    w[10] = 0.11336181654631966654944071844;
    w[11] = 0.11666044348529658204466250754;
    w[12] = 0.11832141527926227651637108570;
    w[13] = 0.11832141527926227651637108570;
    w[14] = 0.11666044348529658204466250754;
    w[15] = 0.11336181654631966654944071844;
    w[16] = 0.10847184052857659065657942673;
    w[17] = 0.10205916109442542323841407025;
    w[18] = 0.09421380035591414846366488307;
    w[19] = 0.08504589431348523921044776508;
    w[20] = 0.07468414976565974588707579610;
    w[21] = 0.063274046329574835539453689907;
    w[22] = 0.050975825297147811998319900724;
    w[23] = 0.037962383294362763950303141249;
    w[24] = 0.024417851092631908789615827520;
    w[25] = 0.010551372617343007155651187685;
  }
  else if (n==27) {
    w[ 0] = 0.0097989960512943602611500550912;
    w[ 1] = 0.022686231596180623196034206447;
    w[ 2] = 0.035297053757419711022578289305;
    w[ 3] = 0.047449412520615062704096710114;
    w[ 4] = 0.058983536859833599110300833720;
    w[ 5] = 0.069748823766245592984322888357;
    w[ 6] = 0.079604867773057771263074959010;
    w[ 7] = 0.088423158543756950194322802854;
    w[ 8] = 0.096088727370028507565652646558;
    w[ 9] = 0.102501637817745798671247711533;
    w[10] = 0.107578285788533187212162984427;
    w[11] = 0.111252488356845192672163096043;
    w[12] = 0.113476346108965148620369948092;
    w[13] = 0.11422086737895698904504573690;
    w[14] = 0.113476346108965148620369948092;
    w[15] = 0.111252488356845192672163096043;
    w[16] = 0.107578285788533187212162984427;
    w[17] = 0.102501637817745798671247711533;
    w[18] = 0.096088727370028507565652646558;
    w[19] = 0.088423158543756950194322802854;
    w[20] = 0.079604867773057771263074959010;
    w[21] = 0.069748823766245592984322888357;
    w[22] = 0.058983536859833599110300833720;
    w[23] = 0.047449412520615062704096710114;
    w[24] = 0.035297053757419711022578289305;
    w[25] = 0.022686231596180623196034206447;
    w[26] = 0.0097989960512943602611500550912;
  }
  else if (n==28) {
    w[ 0] = 0.009124282593094517738816153923;
    w[ 1] = 0.021132112592771259751500380993;
    w[ 2] = 0.032901427782304379977630819171;
    w[ 3] = 0.044272934759004227839587877653;
    w[ 4] = 0.055107345675716745431482918227;
    w[ 5] = 0.06527292396699959579339756678;
    w[ 6] = 0.07464621423456877902393188717;
    w[ 7] = 0.08311341722890121839039649824;
    w[ 8] = 0.09057174439303284094218603134;
    w[ 9] = 0.09693065799792991585048900610;
    w[10] = 0.10211296757806076981421663851;
    w[11] = 0.10605576592284641791041643700;
    w[12] = 0.10871119225829413525357151930;
    w[13] = 0.11004701301647519628237626560;
    w[14] = 0.11004701301647519628237626560;
    w[15] = 0.10871119225829413525357151930;
    w[16] = 0.10605576592284641791041643700;
    w[17] = 0.10211296757806076981421663851;
    w[18] = 0.09693065799792991585048900610;
    w[19] = 0.09057174439303284094218603134;
    w[20] = 0.08311341722890121839039649824;
    w[21] = 0.07464621423456877902393188717;
    w[22] = 0.06527292396699959579339756678;
    w[23] = 0.055107345675716745431482918227;
    w[24] = 0.044272934759004227839587877653;
    w[25] = 0.032901427782304379977630819171;
    w[26] = 0.021132112592771259751500380993;
    w[27] = 0.009124282593094517738816153923;
  }
  else if (n==29) {
    w[ 0] = 0.0085169038787464096542638133022;
    w[ 1] = 0.019732085056122705983859801640;
    w[ 2] = 0.030740492202093622644408525375;
    w[ 3] = 0.041402062518682836104830010114;
    w[ 4] = 0.051594826902497923912594381180;
    w[ 5] = 0.061203090657079138542109848024;
    w[ 6] = 0.070117933255051278569581486949;
    w[ 7] = 0.078238327135763783828144888660;
    w[ 8] = 0.085472257366172527545344849297;
    w[ 9] = 0.091737757139258763347966411077;
    w[10] = 0.096963834094408606301900074883;
    w[11] = 0.101091273759914966121820546907;
    w[12] = 0.104073310077729373913328471285;
    w[13] = 0.105876155097320941406591327852;
    w[14] = 0.10647938171831424424651112691;
    w[15] = 0.105876155097320941406591327852;
    w[16] = 0.104073310077729373913328471285;
    w[17] = 0.101091273759914966121820546907;
    w[18] = 0.096963834094408606301900074883;
    w[19] = 0.091737757139258763347966411077;
    w[20] = 0.085472257366172527545344849297;
    w[21] = 0.078238327135763783828144888660;
    w[22] = 0.070117933255051278569581486949;
    w[23] = 0.061203090657079138542109848024;
    w[24] = 0.051594826902497923912594381180;
    w[25] = 0.041402062518682836104830010114;
    w[26] = 0.030740492202093622644408525375;
    w[27] = 0.019732085056122705983859801640;
    w[28] = 0.0085169038787464096542638133022;
  }
  else if (n==30) {
    w[ 0] = 0.007968192496166605615465883475;
    w[ 1] = 0.018466468311090959142302131912;
    w[ 2] = 0.028784707883323369349719179611;
    w[ 3] = 0.038799192569627049596801936446;
    w[ 4] = 0.048402672830594052902938140423;
    w[ 5] = 0.057493156217619066481721689402;
    w[ 6] = 0.06597422988218049512812851512;
    w[ 7] = 0.07375597473770520626824385002;
    w[ 8] = 0.08075589522942021535469493846;
    w[ 9] = 0.08689978720108297980238753072;
    w[10] = 0.09212252223778612871763270709;
    w[11] = 0.09636873717464425963946862635;
    w[12] = 0.09959342058679526706278028210;
    w[13] = 0.10176238974840550459642895217;
    w[14] = 0.10285265289355884034128563671;
    w[15] = 0.10285265289355884034128563671;
    w[16] = 0.10176238974840550459642895217;
    w[17] = 0.09959342058679526706278028210;
    w[18] = 0.09636873717464425963946862635;
    w[19] = 0.09212252223778612871763270709;
    w[20] = 0.08689978720108297980238753072;
    w[21] = 0.08075589522942021535469493846;
    w[22] = 0.07375597473770520626824385002;
    w[23] = 0.06597422988218049512812851512;
    w[24] = 0.057493156217619066481721689402;
    w[25] = 0.048402672830594052902938140423;
    w[26] = 0.038799192569627049596801936446;
    w[27] = 0.028784707883323369349719179611;
    w[28] = 0.018466468311090959142302131912;
    w[29] = 0.007968192496166605615465883475;
  }
  else if (n==31) {
    w[ 0] = 0.0074708315792487758586968750322;
    w[ 1] = 0.017318620790310582463157996087;
    w[ 2] = 0.027009019184979421800608708092;
    w[ 3] = 0.036432273912385464024392010468;
    w[ 4] = 0.045493707527201102902315857895;
    w[ 5] = 0.054103082424916853711666259087;
    w[ 6] = 0.062174786561028426910343543687;
    w[ 7] = 0.069628583235410366167756126255;
    w[ 8] = 0.076390386598776616426357674901;
    w[ 9] = 0.082392991761589263903823367432;
    w[10] = 0.087576740608477876126198069695;
    w[11] = 0.091890113893641478215362871607;
    w[12] = 0.095290242912319512807204197488;
    w[13] = 0.097743335386328725093474010979;
    w[14] = 0.099225011226672307874875514429;
    w[15] = 0.09972054479342645142753383373;
    w[16] = 0.099225011226672307874875514429;
    w[17] = 0.097743335386328725093474010979;
    w[18] = 0.095290242912319512807204197488;
    w[19] = 0.091890113893641478215362871607;
    w[20] = 0.087576740608477876126198069695;
    w[21] = 0.082392991761589263903823367432;
    w[22] = 0.076390386598776616426357674901;
    w[23] = 0.069628583235410366167756126255;
    w[24] = 0.062174786561028426910343543687;
    w[25] = 0.054103082424916853711666259087;
    w[26] = 0.045493707527201102902315857895;
    w[27] = 0.036432273912385464024392010468;
    w[28] = 0.027009019184979421800608708092;
    w[29] = 0.017318620790310582463157996087;
    w[30] = 0.0074708315792487758586968750322;
  }
  else if (n==32) {
    w[ 0] = 0.007018610009470096600407063739;
    w[ 1] = 0.016274394730905670605170562206;
    w[ 2] = 0.025392065309262059455752589789;
    w[ 3] = 0.034273862913021433102687732252;
    w[ 4] = 0.042835898022226680656878646606;
    w[ 5] = 0.050998059262376176196163244690;
    w[ 6] = 0.058684093478535547145283637300;
    w[ 7] = 0.06582222277636184683765006371;
    w[ 8] = 0.07234579410884850622539935648;
    w[ 9] = 0.07819389578707030647174091883;
    w[10] = 0.08331192422694675522219907460;
    w[11] = 0.08765209300440381114277146275;
    w[12] = 0.09117387869576388471286857711;
    w[13] = 0.09384439908080456563918023767;
    w[14] = 0.09563872007927485941908200220;
    w[15] = 0.09654008851472780056676483006;
    w[16] = 0.09654008851472780056676483006;
    w[17] = 0.09563872007927485941908200220;
    w[18] = 0.09384439908080456563918023767;
    w[19] = 0.09117387869576388471286857711;
    w[20] = 0.08765209300440381114277146275;
    w[21] = 0.08331192422694675522219907460;
    w[22] = 0.07819389578707030647174091883;
    w[23] = 0.07234579410884850622539935648;
    w[24] = 0.06582222277636184683765006371;
    w[25] = 0.058684093478535547145283637300;
    w[26] = 0.050998059262376176196163244690;
    w[27] = 0.042835898022226680656878646606;
    w[28] = 0.034273862913021433102687732252;
    w[29] = 0.025392065309262059455752589789;
    w[30] = 0.016274394730905670605170562206;
    w[31] = 0.007018610009470096600407063739;
  }
  else if (n==33) {
    w[ 0] = 0.0066062278475873780586492352085;
    w[ 1] = 0.015321701512934676127945768534;
    w[ 2] = 0.023915548101749480350533257529;
    w[ 3] = 0.032300358632328953281561447250;
    w[ 4] = 0.040401541331669591563409790527;
    w[ 5] = 0.048147742818711695670146880138;
    w[ 6] = 0.055470846631663561284944495439;
    w[ 7] = 0.062306482530317480031627725771;
    w[ 8] = 0.068594572818656712805955073015;
    w[ 9] = 0.074279854843954149342472175919;
    w[10] = 0.079312364794886738363908384942;
    w[11] = 0.083647876067038707613928014518;
    w[12] = 0.087248287618844337607281670945;
    w[13] = 0.090081958660638577239743705500;
    w[14] = 0.092123986643316846213240977717;
    w[15] = 0.093356426065596116160999126274;
    w[16] = 0.09376844616020999656730454155;
    w[17] = 0.093356426065596116160999126274;
    w[18] = 0.092123986643316846213240977717;
    w[19] = 0.090081958660638577239743705500;
    w[20] = 0.087248287618844337607281670945;
    w[21] = 0.083647876067038707613928014518;
    w[22] = 0.079312364794886738363908384942;
    w[23] = 0.074279854843954149342472175919;
    w[24] = 0.068594572818656712805955073015;
    w[25] = 0.062306482530317480031627725771;
    w[26] = 0.055470846631663561284944495439;
    w[27] = 0.048147742818711695670146880138;
    w[28] = 0.040401541331669591563409790527;
    w[29] = 0.032300358632328953281561447250;
    w[30] = 0.023915548101749480350533257529;
    w[31] = 0.015321701512934676127945768534;
    w[32] = 0.0066062278475873780586492352085;
  }
  else {
    std::cerr << "\n";
    std::cerr << "LEGENDRE_LOOKUP_WEIGHTS - Fatal error!\n";
    std::cerr << "  Illegal value of N = " << n << "\n";
    std::cerr << "  Legal values are 1 through 33.\n";
    std::exit(1);
  }
  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::patterson_lookup ( int n, Scalar x[], Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    PATTERSON_LOOKUP looks up Patterson quadrature points and weights.
//
//  Discussion:
//
//    Our convention is that the abscissas are numbered from left to right.
//
//    The rule is defined on [-1,1],
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Prem Kythe, Michael Schaeferkotter,
//    Handbook of Computational Methods for Integration,
//    Chapman and Hall, 2004,
//    ISBN: 1-58488-428-2,
//    LC: QA299.3.K98.
//
//    Thomas Patterson,
//    The Optimal Addition of Points to Quadrature Formulae,
//    Mathematics of Computation,
//    Volume 22, Number 104, October 1968, pages 847-856.
//
//  Parameters:
//
//    Input, int N, the order.
//    Legal values are 1, 3, 7, 15, 31, 63, 127 and 255.
//
//    Output, Scalar X[N], the abscissas.
//
//    Output, Scalar W[N], the weights.
//
{
  IntrepidBurkardtRules::patterson_lookup_points ( n, x );
  IntrepidBurkardtRules::patterson_lookup_weights ( n, w );

  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::patterson_lookup_points ( int n, Scalar x[] )
//****************************************************************************
//
//  Purpose:
//
//    PATTERSON_LOOKUP_POINTS looks up Patterson quadrature points.
//
//  Discussion:
//
//    Our convention is that the abscissas are numbered from left to right.
//
//    The rule is defined on [-1,1],
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 December 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Prem Kythe, Michael Schaeferkotter,
//    Handbook of Computational Methods for Integration,
//    Chapman and Hall, 2004,
//    ISBN: 1-58488-428-2,
//    LC: QA299.3.K98.
//
//    Thomas Patterson,
//    The Optimal Addition of Points to Quadrature Formulae,
//    Mathematics of Computation,
//    Volume 22, Number 104, October 1968, pages 847-856.
//
//  Parameters:
//
//    Input, int N, the order.
//    Legal values are 1, 3, 7, 15, 31, 63, 127 and 255.
//
//    Output, Scalar X[N], the abscissas.
//
{
  if (n==1) {
    x[ 0] =  0.0;
  }
  else if (n==3) {
    x[ 0] = -0.77459666924148337704;
    x[ 1] =  0.0;
    x[ 2] =  0.77459666924148337704;
  }
  else if (n==7) {
    x[ 0] = -0.96049126870802028342;
    x[ 1] = -0.77459666924148337704;
    x[ 2] = -0.43424374934680255800;
    x[ 3] =  0.0;
    x[ 4] =  0.43424374934680255800;
    x[ 5] =  0.77459666924148337704;
    x[ 6] =  0.96049126870802028342;
  }
  else if (n==15) {
    x[ 0] = -0.99383196321275502221;
    x[ 1] = -0.96049126870802028342;
    x[ 2] = -0.88845923287225699889;
    x[ 3] = -0.77459666924148337704;
    x[ 4] = -0.62110294673722640294;
    x[ 5] = -0.43424374934680255800;
    x[ 6] = -0.22338668642896688163;
    x[ 7] =  0.0;
    x[ 8] =  0.22338668642896688163;
    x[ 9] =  0.43424374934680255800;
    x[10] =  0.62110294673722640294;
    x[11] =  0.77459666924148337704;
    x[12] =  0.88845923287225699889;
    x[13] =  0.96049126870802028342;
    x[14] =  0.99383196321275502221;
  }
  else if (n==31) {
    x[ 0] = -0.99909812496766759766;
    x[ 1] = -0.99383196321275502221;
    x[ 2] = -0.98153114955374010687;
    x[ 3] = -0.96049126870802028342;
    x[ 4] = -0.92965485742974005667;
    x[ 5] = -0.88845923287225699889;
    x[ 6] = -0.83672593816886873550;
    x[ 7] = -0.77459666924148337704;
    x[ 8] = -0.70249620649152707861;
    x[ 9] = -0.62110294673722640294;
    x[10] = -0.53131974364437562397;
    x[11] = -0.43424374934680255800;
    x[12] = -0.33113539325797683309;
    x[13] = -0.22338668642896688163;
    x[14] = -0.11248894313318662575;
    x[15] =  0.0;
    x[16] =  0.11248894313318662575;
    x[17] =  0.22338668642896688163;
    x[18] =  0.33113539325797683309;
    x[19] =  0.43424374934680255800;
    x[20] =  0.53131974364437562397;
    x[21] =  0.62110294673722640294;
    x[22] =  0.70249620649152707861;
    x[23] =  0.77459666924148337704;
    x[24] =  0.83672593816886873550;
    x[25] =  0.88845923287225699889;
    x[26] =  0.92965485742974005667;
    x[27] =  0.96049126870802028342;
    x[28] =  0.98153114955374010687;
    x[29] =  0.99383196321275502221;
    x[30] =  0.99909812496766759766;
  }
  else if (n==63) {
    x[ 0] = -0.99987288812035761194;
    x[ 1] = -0.99909812496766759766;
    x[ 2] = -0.99720625937222195908;
    x[ 3] = -0.99383196321275502221;
    x[ 4] = -0.98868475754742947994;
    x[ 5] = -0.98153114955374010687;
    x[ 6] = -0.97218287474858179658;
    x[ 7] = -0.96049126870802028342;
    x[ 8] = -0.94634285837340290515;
    x[ 9] = -0.92965485742974005667;
    x[10] = -0.91037115695700429250;
    x[11] = -0.88845923287225699889;
    x[12] = -0.86390793819369047715;
    x[13] = -0.83672593816886873550;
    x[14] = -0.80694053195021761186;
    x[15] = -0.77459666924148337704;
    x[16] = -0.73975604435269475868;
    x[17] = -0.70249620649152707861;
    x[18] = -0.66290966002478059546;
    x[19] = -0.62110294673722640294;
    x[20] = -0.57719571005204581484;
    x[21] = -0.53131974364437562397;
    x[22] = -0.48361802694584102756;
    x[23] = -0.43424374934680255800;
    x[24] = -0.38335932419873034692;
    x[25] = -0.33113539325797683309;
    x[26] = -0.27774982202182431507;
    x[27] = -0.22338668642896688163;
    x[28] = -0.16823525155220746498;
    x[29] = -0.11248894313318662575;
    x[30] = -0.056344313046592789972;
    x[31] =  0.0;
    x[32] =  0.056344313046592789972;
    x[33] =  0.11248894313318662575;
    x[34] =  0.16823525155220746498;
    x[35] =  0.22338668642896688163;
    x[36] =  0.27774982202182431507;
    x[37] =  0.33113539325797683309;
    x[38] =  0.38335932419873034692;
    x[39] =  0.43424374934680255800;
    x[40] =  0.48361802694584102756;
    x[41] =  0.53131974364437562397;
    x[42] =  0.57719571005204581484;
    x[43] =  0.62110294673722640294;
    x[44] =  0.66290966002478059546;
    x[45] =  0.70249620649152707861;
    x[46] =  0.73975604435269475868;
    x[47] =  0.77459666924148337704;
    x[48] =  0.80694053195021761186;
    x[49] =  0.83672593816886873550;
    x[50] =  0.86390793819369047715;
    x[51] =  0.88845923287225699889;
    x[52] =  0.91037115695700429250;
    x[53] =  0.92965485742974005667;
    x[54] =  0.94634285837340290515;
    x[55] =  0.96049126870802028342;
    x[56] =  0.97218287474858179658;
    x[57] =  0.98153114955374010687;
    x[58] =  0.98868475754742947994;
    x[59] =  0.99383196321275502221;
    x[60] =  0.99720625937222195908;
    x[61] =  0.99909812496766759766;
    x[62] =  0.99987288812035761194;
  }
  else if (n==127) {
    x[  0] = -0.99998243035489159858;
    x[  1] = -0.99987288812035761194;
    x[  2] = -0.99959879967191068325;
    x[  3] = -0.99909812496766759766;
    x[  4] = -0.99831663531840739253;
    x[  5] = -0.99720625937222195908;
    x[  6] = -0.99572410469840718851;
    x[  7] = -0.99383196321275502221;
    x[  8] = -0.99149572117810613240;
    x[  9] = -0.98868475754742947994;
    x[ 10] = -0.98537149959852037111;
    x[ 11] = -0.98153114955374010687;
    x[ 12] = -0.97714151463970571416;
    x[ 13] = -0.97218287474858179658;
    x[ 14] = -0.96663785155841656709;
    x[ 15] = -0.96049126870802028342;
    x[ 16] = -0.95373000642576113641;
    x[ 17] = -0.94634285837340290515;
    x[ 18] = -0.93832039777959288365;
    x[ 19] = -0.92965485742974005667;
    x[ 20] = -0.92034002547001242073;
    x[ 21] = -0.91037115695700429250;
    x[ 22] = -0.89974489977694003664;
    x[ 23] = -0.88845923287225699889;
    x[ 24] = -0.87651341448470526974;
    x[ 25] = -0.86390793819369047715;
    x[ 26] = -0.85064449476835027976;
    x[ 27] = -0.83672593816886873550;
    x[ 28] = -0.82215625436498040737;
    x[ 29] = -0.80694053195021761186;
    x[ 30] = -0.79108493379984836143;
    x[ 31] = -0.77459666924148337704;
    x[ 32] = -0.75748396638051363793;
    x[ 33] = -0.73975604435269475868;
    x[ 34] = -0.72142308537009891548;
    x[ 35] = -0.70249620649152707861;
    x[ 36] = -0.68298743109107922809;
    x[ 37] = -0.66290966002478059546;
    x[ 38] = -0.64227664250975951377;
    x[ 39] = -0.62110294673722640294;
    x[ 40] = -0.59940393024224289297;
    x[ 41] = -0.57719571005204581484;
    x[ 42] = -0.55449513263193254887;
    x[ 43] = -0.53131974364437562397;
    x[ 44] = -0.50768775753371660215;
    x[ 45] = -0.48361802694584102756;
    x[ 46] = -0.45913001198983233287;
    x[ 47] = -0.43424374934680255800;
    x[ 48] = -0.40897982122988867241;
    x[ 49] = -0.38335932419873034692;
    x[ 50] = -0.35740383783153215238;
    x[ 51] = -0.33113539325797683309;
    x[ 52] = -0.30457644155671404334;
    x[ 53] = -0.27774982202182431507;
    x[ 54] = -0.25067873030348317661;
    x[ 55] = -0.22338668642896688163;
    x[ 56] = -0.19589750271110015392;
    x[ 57] = -0.16823525155220746498;
    x[ 58] = -0.14042423315256017459;
    x[ 59] = -0.11248894313318662575;
    x[ 60] = -0.084454040083710883710;
    x[ 61] = -0.056344313046592789972;
    x[ 62] = -0.028184648949745694339;
    x[ 63] =  0.0;
    x[ 64] =  0.028184648949745694339;
    x[ 65] =  0.056344313046592789972;
    x[ 66] =  0.084454040083710883710;
    x[ 67] =  0.11248894313318662575;
    x[ 68] =  0.14042423315256017459;
    x[ 69] =  0.16823525155220746498;
    x[ 70] =  0.19589750271110015392;
    x[ 71] =  0.22338668642896688163;
    x[ 72] =  0.25067873030348317661;
    x[ 73] =  0.27774982202182431507;
    x[ 74] =  0.30457644155671404334;
    x[ 75] =  0.33113539325797683309;
    x[ 76] =  0.35740383783153215238;
    x[ 77] =  0.38335932419873034692;
    x[ 78] =  0.40897982122988867241;
    x[ 79] =  0.43424374934680255800;
    x[ 80] =  0.45913001198983233287;
    x[ 81] =  0.48361802694584102756;
    x[ 82] =  0.50768775753371660215;
    x[ 83] =  0.53131974364437562397;
    x[ 84] =  0.55449513263193254887;
    x[ 85] =  0.57719571005204581484;
    x[ 86] =  0.59940393024224289297;
    x[ 87] =  0.62110294673722640294;
    x[ 88] =  0.64227664250975951377;
    x[ 89] =  0.66290966002478059546;
    x[ 90] =  0.68298743109107922809;
    x[ 91] =  0.70249620649152707861;
    x[ 92] =  0.72142308537009891548;
    x[ 93] =  0.73975604435269475868;
    x[ 94] =  0.75748396638051363793;
    x[ 95] =  0.77459666924148337704;
    x[ 96] =  0.79108493379984836143;
    x[ 97] =  0.80694053195021761186;
    x[ 98] =  0.82215625436498040737;
    x[ 99] =  0.83672593816886873550;
    x[100] =  0.85064449476835027976;
    x[101] =  0.86390793819369047715;
    x[102] =  0.87651341448470526974;
    x[103] =  0.88845923287225699889;
    x[104] =  0.89974489977694003664;
    x[105] =  0.91037115695700429250;
    x[106] =  0.92034002547001242073;
    x[107] =  0.92965485742974005667;
    x[108] =  0.93832039777959288365;
    x[109] =  0.94634285837340290515;
    x[110] =  0.95373000642576113641;
    x[111] =  0.96049126870802028342;
    x[112] =  0.96663785155841656709;
    x[113] =  0.97218287474858179658;
    x[114] =  0.97714151463970571416;
    x[115] =  0.98153114955374010687;
    x[116] =  0.98537149959852037111;
    x[117] =  0.98868475754742947994;
    x[118] =  0.99149572117810613240;
    x[119] =  0.99383196321275502221;
    x[120] =  0.99572410469840718851;
    x[121] =  0.99720625937222195908;
    x[122] =  0.99831663531840739253;
    x[123] =  0.99909812496766759766;
    x[124] =  0.99959879967191068325;
    x[125] =  0.99987288812035761194;
    x[126] =  0.99998243035489159858;
  }
  else if (n==255) {
    x[  0] = -0.99999759637974846462;
    x[  1] = -0.99998243035489159858;
    x[  2] = -0.99994399620705437576;
    x[  3] = -0.99987288812035761194;
    x[  4] = -0.99976049092443204733;
    x[  5] = -0.99959879967191068325;
    x[  6] = -0.99938033802502358193;
    x[  7] = -0.99909812496766759766;
    x[  8] = -0.99874561446809511470;
    x[  9] = -0.99831663531840739253;
    x[ 10] = -0.99780535449595727456;
    x[ 11] = -0.99720625937222195908;
    x[ 12] = -0.99651414591489027385;
    x[ 13] = -0.99572410469840718851;
    x[ 14] = -0.99483150280062100052;
    x[ 15] = -0.99383196321275502221;
    x[ 16] = -0.99272134428278861533;
    x[ 17] = -0.99149572117810613240;
    x[ 18] = -0.99015137040077015918;
    x[ 19] = -0.98868475754742947994;
    x[ 20] = -0.98709252795403406719;
    x[ 21] = -0.98537149959852037111;
    x[ 22] = -0.98351865757863272876;
    x[ 23] = -0.98153114955374010687;
    x[ 24] = -0.97940628167086268381;
    x[ 25] = -0.97714151463970571416;
    x[ 26] = -0.97473445975240266776;
    x[ 27] = -0.97218287474858179658;
    x[ 28] = -0.96948465950245923177;
    x[ 29] = -0.96663785155841656709;
    x[ 30] = -0.96364062156981213252;
    x[ 31] = -0.96049126870802028342;
    x[ 32] = -0.95718821610986096274;
    x[ 33] = -0.95373000642576113641;
    x[ 34] = -0.95011529752129487656;
    x[ 35] = -0.94634285837340290515;
    x[ 36] = -0.94241156519108305981;
    x[ 37] = -0.93832039777959288365;
    x[ 38] = -0.93406843615772578800;
    x[ 39] = -0.92965485742974005667;
    x[ 40] = -0.92507893290707565236;
    x[ 41] = -0.92034002547001242073;
    x[ 42] = -0.91543758715576504064;
    x[ 43] = -0.91037115695700429250;
    x[ 44] = -0.90514035881326159519;
    x[ 45] = -0.89974489977694003664;
    x[ 46] = -0.89418456833555902286;
    x[ 47] = -0.88845923287225699889;
    x[ 48] = -0.88256884024734190684;
    x[ 49] = -0.87651341448470526974;
    x[ 50] = -0.87029305554811390585;
    x[ 51] = -0.86390793819369047715;
    x[ 52] = -0.85735831088623215653;
    x[ 53] = -0.85064449476835027976;
    x[ 54] = -0.84376688267270860104;
    x[ 55] = -0.83672593816886873550;
    x[ 56] = -0.82952219463740140018;
    x[ 57] = -0.82215625436498040737;
    x[ 58] = -0.81462878765513741344;
    x[ 59] = -0.80694053195021761186;
    x[ 60] = -0.79909229096084140180;
    x[ 61] = -0.79108493379984836143;
    x[ 62] = -0.78291939411828301639;
    x[ 63] = -0.77459666924148337704;
    x[ 64] = -0.76611781930376009072;
    x[ 65] = -0.75748396638051363793;
    x[ 66] = -0.74869629361693660282;
    x[ 67] = -0.73975604435269475868;
    x[ 68] = -0.73066452124218126133;
    x[ 69] = -0.72142308537009891548;
    x[ 70] = -0.71203315536225203459;
    x[ 71] = -0.70249620649152707861;
    x[ 72] = -0.69281376977911470289;
    x[ 73] = -0.68298743109107922809;
    x[ 74] = -0.67301883023041847920;
    x[ 75] = -0.66290966002478059546;
    x[ 76] = -0.65266166541001749610;
    x[ 77] = -0.64227664250975951377;
    x[ 78] = -0.63175643771119423041;
    x[ 79] = -0.62110294673722640294;
    x[ 80] = -0.61031811371518640016;
    x[ 81] = -0.59940393024224289297;
    x[ 82] = -0.58836243444766254143;
    x[ 83] = -0.57719571005204581484;
    x[ 84] = -0.56590588542365442262;
    x[ 85] = -0.55449513263193254887;
    x[ 86] = -0.54296566649831149049;
    x[ 87] = -0.53131974364437562397;
    x[ 88] = -0.51955966153745702199;
    x[ 89] = -0.50768775753371660215;
    x[ 90] = -0.49570640791876146017;
    x[ 91] = -0.48361802694584102756;
    x[ 92] = -0.47142506587165887693;
    x[ 93] = -0.45913001198983233287;
    x[ 94] = -0.44673538766202847374;
    x[ 95] = -0.43424374934680255800;
    x[ 96] = -0.42165768662616330006;
    x[ 97] = -0.40897982122988867241;
    x[ 98] = -0.39621280605761593918;
    x[ 99] = -0.38335932419873034692;
    x[100] = -0.37042208795007823014;
    x[101] = -0.35740383783153215238;
    x[102] = -0.34430734159943802278;
    x[103] = -0.33113539325797683309;
    x[104] = -0.31789081206847668318;
    x[105] = -0.30457644155671404334;
    x[106] = -0.29119514851824668196;
    x[107] = -0.27774982202182431507;
    x[108] = -0.26424337241092676194;
    x[109] = -0.25067873030348317661;
    x[110] = -0.23705884558982972721;
    x[111] = -0.22338668642896688163;
    x[112] = -0.20966523824318119477;
    x[113] = -0.19589750271110015392;
    x[114] = -0.18208649675925219825;
    x[115] = -0.16823525155220746498;
    x[116] = -0.15434681148137810869;
    x[117] = -0.14042423315256017459;
    x[118] = -0.12647058437230196685;
    x[119] = -0.11248894313318662575;
    x[120] = -0.098482396598119202090;
    x[121] = -0.084454040083710883710;
    x[122] = -0.070406976042855179063;
    x[123] = -0.056344313046592789972;
    x[124] = -0.042269164765363603212;
    x[125] = -0.028184648949745694339;
    x[126] = -0.014093886410782462614;
    x[127] =  0.0;
    x[128] =  0.014093886410782462614;
    x[129] =  0.028184648949745694339;
    x[130] =  0.042269164765363603212;
    x[131] =  0.056344313046592789972;
    x[132] =  0.070406976042855179063;
    x[133] =  0.084454040083710883710;
    x[134] =  0.098482396598119202090;
    x[135] =  0.11248894313318662575;
    x[136] =  0.12647058437230196685;
    x[137] =  0.14042423315256017459;
    x[138] =  0.15434681148137810869;
    x[139] =  0.16823525155220746498;
    x[140] =  0.18208649675925219825;
    x[141] =  0.19589750271110015392;
    x[142] =  0.20966523824318119477;
    x[143] =  0.22338668642896688163;
    x[144] =  0.23705884558982972721;
    x[145] =  0.25067873030348317661;
    x[146] =  0.26424337241092676194;
    x[147] =  0.27774982202182431507;
    x[148] =  0.29119514851824668196;
    x[149] =  0.30457644155671404334;
    x[150] =  0.31789081206847668318;
    x[151] =  0.33113539325797683309;
    x[152] =  0.34430734159943802278;
    x[153] =  0.35740383783153215238;
    x[154] =  0.37042208795007823014;
    x[155] =  0.38335932419873034692;
    x[156] =  0.39621280605761593918;
    x[157] =  0.40897982122988867241;
    x[158] =  0.42165768662616330006;
    x[159] =  0.43424374934680255800;
    x[160] =  0.44673538766202847374;
    x[161] =  0.45913001198983233287;
    x[162] =  0.47142506587165887693;
    x[163] =  0.48361802694584102756;
    x[164] =  0.49570640791876146017;
    x[165] =  0.50768775753371660215;
    x[166] =  0.51955966153745702199;
    x[167] =  0.53131974364437562397;
    x[168] =  0.54296566649831149049;
    x[169] =  0.55449513263193254887;
    x[170] =  0.56590588542365442262;
    x[171] =  0.57719571005204581484;
    x[172] =  0.58836243444766254143;
    x[173] =  0.59940393024224289297;
    x[174] =  0.61031811371518640016;
    x[175] =  0.62110294673722640294;
    x[176] =  0.63175643771119423041;
    x[177] =  0.64227664250975951377;
    x[178] =  0.65266166541001749610;
    x[179] =  0.66290966002478059546;
    x[180] =  0.67301883023041847920;
    x[181] =  0.68298743109107922809;
    x[182] =  0.69281376977911470289;
    x[183] =  0.70249620649152707861;
    x[184] =  0.71203315536225203459;
    x[185] =  0.72142308537009891548;
    x[186] =  0.73066452124218126133;
    x[187] =  0.73975604435269475868;
    x[188] =  0.74869629361693660282;
    x[189] =  0.75748396638051363793;
    x[190] =  0.76611781930376009072;
    x[191] =  0.77459666924148337704;
    x[192] =  0.78291939411828301639;
    x[193] =  0.79108493379984836143;
    x[194] =  0.79909229096084140180;
    x[195] =  0.80694053195021761186;
    x[196] =  0.81462878765513741344;
    x[197] =  0.82215625436498040737;
    x[198] =  0.82952219463740140018;
    x[199] =  0.83672593816886873550;
    x[200] =  0.84376688267270860104;
    x[201] =  0.85064449476835027976;
    x[202] =  0.85735831088623215653;
    x[203] =  0.86390793819369047715;
    x[204] =  0.87029305554811390585;
    x[205] =  0.87651341448470526974;
    x[206] =  0.88256884024734190684;
    x[207] =  0.88845923287225699889;
    x[208] =  0.89418456833555902286;
    x[209] =  0.89974489977694003664;
    x[210] =  0.90514035881326159519;
    x[211] =  0.91037115695700429250;
    x[212] =  0.91543758715576504064;
    x[213] =  0.92034002547001242073;
    x[214] =  0.92507893290707565236;
    x[215] =  0.92965485742974005667;
    x[216] =  0.93406843615772578800;
    x[217] =  0.93832039777959288365;
    x[218] =  0.94241156519108305981;
    x[219] =  0.94634285837340290515;
    x[220] =  0.95011529752129487656;
    x[221] =  0.95373000642576113641;
    x[222] =  0.95718821610986096274;
    x[223] =  0.96049126870802028342;
    x[224] =  0.96364062156981213252;
    x[225] =  0.96663785155841656709;
    x[226] =  0.96948465950245923177;
    x[227] =  0.97218287474858179658;
    x[228] =  0.97473445975240266776;
    x[229] =  0.97714151463970571416;
    x[230] =  0.97940628167086268381;
    x[231] =  0.98153114955374010687;
    x[232] =  0.98351865757863272876;
    x[233] =  0.98537149959852037111;
    x[234] =  0.98709252795403406719;
    x[235] =  0.98868475754742947994;
    x[236] =  0.99015137040077015918;
    x[237] =  0.99149572117810613240;
    x[238] =  0.99272134428278861533;
    x[239] =  0.99383196321275502221;
    x[240] =  0.99483150280062100052;
    x[241] =  0.99572410469840718851;
    x[242] =  0.99651414591489027385;
    x[243] =  0.99720625937222195908;
    x[244] =  0.99780535449595727456;
    x[245] =  0.99831663531840739253;
    x[246] =  0.99874561446809511470;
    x[247] =  0.99909812496766759766;
    x[248] =  0.99938033802502358193;
    x[249] =  0.99959879967191068325;
    x[250] =  0.99976049092443204733;
    x[251] =  0.99987288812035761194;
    x[252] =  0.99994399620705437576;
    x[253] =  0.99998243035489159858;
    x[254] =  0.99999759637974846462;
  }
  else {
    std::cerr << "\n";
    std::cerr << "PATTERSON_LOOKUP_POINTS - Fatal error!\n";
    std::cerr << "  Unexpected value of N = " << n << "\n";
    std::exit(1);
  }
  return;
}

//****************************************************************************
template<class Scalar> 
void IntrepidBurkardtRules::patterson_lookup_weights ( int n, Scalar w[] )
//****************************************************************************
//
//  Purpose:
//
//    PATTERSON_LOOKUP_WEIGHTS looks up Patterson quadrature weights.
//
//  Discussion:
//
//    The allowed orders are 1, 3, 7, 15, 31, 63, 127 and 255.
//
//    The weights are positive, symmetric and should sum to 2.
//
//    The user must preallocate space for the output array W.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 December 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int N, the order.
//    Legal values are 1, 3, 7, 15, 31, 63, 127 or 255.
//
//    Output, Scalar W[N], the weights.
//
{
  if (n==1) {
    w[  0] = 2.0;
  }
  else if (n==3) {
    w[  0] = 0.555555555555555555556;
    w[  1] = 0.888888888888888888889;
    w[  2] = 0.555555555555555555556;
  }
  else if (n==7) {
    w[  0] = 0.104656226026467265194;
    w[  1] = 0.268488089868333440729;
    w[  2] = 0.401397414775962222905;
    w[  3] = 0.450916538658474142345;
    w[  4] = 0.401397414775962222905;
    w[  5] = 0.268488089868333440729;
    w[  6] = 0.104656226026467265194;
  }
  else if (n==15) {
    w[  0] = 0.0170017196299402603390;
    w[  1] = 0.0516032829970797396969;
    w[  2] = 0.0929271953151245376859;
    w[  3] = 0.134415255243784220360;
    w[  4] = 0.171511909136391380787;
    w[  5] = 0.200628529376989021034;
    w[  6] = 0.219156858401587496404;
    w[  7] = 0.225510499798206687386;
    w[  8] = 0.219156858401587496404;
    w[  9] = 0.200628529376989021034;
    w[ 10] = 0.171511909136391380787;
    w[ 11] = 0.134415255243784220360;
    w[ 12] = 0.0929271953151245376859;
    w[ 13] = 0.0516032829970797396969;
    w[ 14] = 0.0170017196299402603390;
  }
  else if (n==31) {
    w[  0] = 0.00254478079156187441540;
    w[  1] = 0.00843456573932110624631;
    w[  2] = 0.0164460498543878109338;
    w[  3] = 0.0258075980961766535646;
    w[  4] = 0.0359571033071293220968;
    w[  5] = 0.0464628932617579865414;
    w[  6] = 0.0569795094941233574122;
    w[  7] = 0.0672077542959907035404;
    w[  8] = 0.0768796204990035310427;
    w[  9] = 0.0857559200499903511542;
    w[ 10] = 0.0936271099812644736167;
    w[ 11] = 0.100314278611795578771;
    w[ 12] = 0.105669893580234809744;
    w[ 13] = 0.109578421055924638237;
    w[ 14] = 0.111956873020953456880;
    w[ 15] = 0.112755256720768691607;
    w[ 16] = 0.111956873020953456880;
    w[ 17] = 0.109578421055924638237;
    w[ 18] = 0.105669893580234809744;
    w[ 19] = 0.100314278611795578771;
    w[ 20] = 0.0936271099812644736167;
    w[ 21] = 0.0857559200499903511542;
    w[ 22] = 0.0768796204990035310427;
    w[ 23] = 0.0672077542959907035404;
    w[ 24] = 0.0569795094941233574122;
    w[ 25] = 0.0464628932617579865414;
    w[ 26] = 0.0359571033071293220968;
    w[ 27] = 0.0258075980961766535646;
    w[ 28] = 0.0164460498543878109338;
    w[ 29] = 0.00843456573932110624631;
    w[ 30] = 0.00254478079156187441540;
  }
  else if (n==63) {
    w[  0] = 0.000363221481845530659694;
    w[  1] = 0.00126515655623006801137;
    w[  2] = 0.00257904979468568827243;
    w[  3] = 0.00421763044155885483908;
    w[  4] = 0.00611550682211724633968;
    w[  5] = 0.00822300795723592966926;
    w[  6] = 0.0104982469096213218983;
    w[  7] = 0.0129038001003512656260;
    w[  8] = 0.0154067504665594978021;
    w[  9] = 0.0179785515681282703329;
    w[ 10] = 0.0205942339159127111492;
    w[ 11] = 0.0232314466399102694433;
    w[ 12] = 0.0258696793272147469108;
    w[ 13] = 0.0284897547458335486125;
    w[ 14] = 0.0310735511116879648799;
    w[ 15] = 0.0336038771482077305417;
    w[ 16] = 0.0360644327807825726401;
    w[ 17] = 0.0384398102494555320386;
    w[ 18] = 0.0407155101169443189339;
    w[ 19] = 0.0428779600250077344929;
    w[ 20] = 0.0449145316536321974143;
    w[ 21] = 0.0468135549906280124026;
    w[ 22] = 0.0485643304066731987159;
    w[ 23] = 0.0501571393058995374137;
    w[ 24] = 0.0515832539520484587768;
    w[ 25] = 0.0528349467901165198621;
    w[ 26] = 0.0539054993352660639269;
    w[ 27] = 0.0547892105279628650322;
    w[ 28] = 0.0554814043565593639878;
    w[ 29] = 0.0559784365104763194076;
    w[ 30] = 0.0562776998312543012726;
    w[ 31] = 0.0563776283603847173877;
    w[ 32] = 0.0562776998312543012726;
    w[ 33] = 0.0559784365104763194076;
    w[ 34] = 0.0554814043565593639878;
    w[ 35] = 0.0547892105279628650322;
    w[ 36] = 0.0539054993352660639269;
    w[ 37] = 0.0528349467901165198621;
    w[ 38] = 0.0515832539520484587768;
    w[ 39] = 0.0501571393058995374137;
    w[ 40] = 0.0485643304066731987159;
    w[ 41] = 0.0468135549906280124026;
    w[ 42] = 0.0449145316536321974143;
    w[ 43] = 0.0428779600250077344929;
    w[ 44] = 0.0407155101169443189339;
    w[ 45] = 0.0384398102494555320386;
    w[ 46] = 0.0360644327807825726401;
    w[ 47] = 0.0336038771482077305417;
    w[ 48] = 0.0310735511116879648799;
    w[ 49] = 0.0284897547458335486125;
    w[ 50] = 0.0258696793272147469108;
    w[ 51] = 0.0232314466399102694433;
    w[ 52] = 0.0205942339159127111492;
    w[ 53] = 0.0179785515681282703329;
    w[ 54] = 0.0154067504665594978021;
    w[ 55] = 0.0129038001003512656260;
    w[ 56] = 0.0104982469096213218983;
    w[ 57] = 0.00822300795723592966926;
    w[ 58] = 0.00611550682211724633968;
    w[ 59] = 0.00421763044155885483908;
    w[ 60] = 0.00257904979468568827243;
    w[ 61] = 0.00126515655623006801137;
    w[ 62] = 0.000363221481845530659694;
  }
  else if (n==127) {
    w[  0] = 0.0000505360952078625176247;
    w[  1] = 0.000180739564445388357820;
    w[  2] = 0.000377746646326984660274;
    w[  3] = 0.000632607319362633544219;
    w[  4] = 0.000938369848542381500794;
    w[  5] = 0.00128952408261041739210;
    w[  6] = 0.00168114286542146990631;
    w[  7] = 0.00210881524572663287933;
    w[  8] = 0.00256876494379402037313;
    w[  9] = 0.00305775341017553113613;
    w[ 10] = 0.00357289278351729964938;
    w[ 11] = 0.00411150397865469304717;
    w[ 12] = 0.00467105037211432174741;
    w[ 13] = 0.00524912345480885912513;
    w[ 14] = 0.00584344987583563950756;
    w[ 15] = 0.00645190005017573692280;
    w[ 16] = 0.00707248999543355546805;
    w[ 17] = 0.00770337523327974184817;
    w[ 18] = 0.00834283875396815770558;
    w[ 19] = 0.00898927578406413572328;
    w[ 20] = 0.00964117772970253669530;
    w[ 21] = 0.0102971169579563555237;
    w[ 22] = 0.0109557333878379016480;
    w[ 23] = 0.0116157233199551347270;
    w[ 24] = 0.0122758305600827700870;
    w[ 25] = 0.0129348396636073734547;
    w[ 26] = 0.0135915710097655467896;
    w[ 27] = 0.0142448773729167743063;
    w[ 28] = 0.0148936416648151820348;
    w[ 29] = 0.0155367755558439824399;
    w[ 30] = 0.0161732187295777199419;
    w[ 31] = 0.0168019385741038652709;
    w[ 32] = 0.0174219301594641737472;
    w[ 33] = 0.0180322163903912863201;
    w[ 34] = 0.0186318482561387901863;
    w[ 35] = 0.0192199051247277660193;
    w[ 36] = 0.0197954950480974994880;
    w[ 37] = 0.0203577550584721594669;
    w[ 38] = 0.0209058514458120238522;
    w[ 39] = 0.0214389800125038672465;
    w[ 40] = 0.0219563663053178249393;
    w[ 41] = 0.0224572658268160987071;
    w[ 42] = 0.0229409642293877487608;
    w[ 43] = 0.0234067774953140062013;
    w[ 44] = 0.0238540521060385400804;
    w[ 45] = 0.0242821652033365993580;
    w[ 46] = 0.0246905247444876769091;
    w[ 47] = 0.0250785696529497687068;
    w[ 48] = 0.0254457699654647658126;
    w[ 49] = 0.0257916269760242293884;
    w[ 50] = 0.0261156733767060976805;
    w[ 51] = 0.0264174733950582599310;
    w[ 52] = 0.0266966229274503599062;
    w[ 53] = 0.0269527496676330319634;
    w[ 54] = 0.0271855132296247918192;
    w[ 55] = 0.0273946052639814325161;
    w[ 56] = 0.0275797495664818730349;
    w[ 57] = 0.0277407021782796819939;
    w[ 58] = 0.0278772514766137016085;
    w[ 59] = 0.0279892182552381597038;
    w[ 60] = 0.0280764557938172466068;
    w[ 61] = 0.0281388499156271506363;
    w[ 62] = 0.0281763190330166021307;
    w[ 63] = 0.0281888141801923586938;
    w[ 64] = 0.0281763190330166021307;
    w[ 65] = 0.0281388499156271506363;
    w[ 66] = 0.0280764557938172466068;
    w[ 67] = 0.0279892182552381597038;
    w[ 68] = 0.0278772514766137016085;
    w[ 69] = 0.0277407021782796819939;
    w[ 70] = 0.0275797495664818730349;
    w[ 71] = 0.0273946052639814325161;
    w[ 72] = 0.0271855132296247918192;
    w[ 73] = 0.0269527496676330319634;
    w[ 74] = 0.0266966229274503599062;
    w[ 75] = 0.0264174733950582599310;
    w[ 76] = 0.0261156733767060976805;
    w[ 77] = 0.0257916269760242293884;
    w[ 78] = 0.0254457699654647658126;
    w[ 79] = 0.0250785696529497687068;
    w[ 80] = 0.0246905247444876769091;
    w[ 81] = 0.0242821652033365993580;
    w[ 82] = 0.0238540521060385400804;
    w[ 83] = 0.0234067774953140062013;
    w[ 84] = 0.0229409642293877487608;
    w[ 85] = 0.0224572658268160987071;
    w[ 86] = 0.0219563663053178249393;
    w[ 87] = 0.0214389800125038672465;
    w[ 88] = 0.0209058514458120238522;
    w[ 89] = 0.0203577550584721594669;
    w[ 90] = 0.0197954950480974994880;
    w[ 91] = 0.0192199051247277660193;
    w[ 92] = 0.0186318482561387901863;
    w[ 93] = 0.0180322163903912863201;
    w[ 94] = 0.0174219301594641737472;
    w[ 95] = 0.0168019385741038652709;
    w[ 96] = 0.0161732187295777199419;
    w[ 97] = 0.0155367755558439824399;
    w[ 98] = 0.0148936416648151820348;
    w[ 99] = 0.0142448773729167743063;
    w[100] = 0.0135915710097655467896;
    w[101] = 0.0129348396636073734547;
    w[102] = 0.0122758305600827700870;
    w[103] = 0.0116157233199551347270;
    w[104] = 0.0109557333878379016480;
    w[105] = 0.0102971169579563555237;
    w[106] = 0.00964117772970253669530;
    w[107] = 0.00898927578406413572328;
    w[108] = 0.00834283875396815770558;
    w[109] = 0.00770337523327974184817;
    w[110] = 0.00707248999543355546805;
    w[111] = 0.00645190005017573692280;
    w[112] = 0.00584344987583563950756;
    w[113] = 0.00524912345480885912513;
    w[114] = 0.00467105037211432174741;
    w[115] = 0.00411150397865469304717;
    w[116] = 0.00357289278351729964938;
    w[117] = 0.00305775341017553113613;
    w[118] = 0.00256876494379402037313;
    w[119] = 0.00210881524572663287933;
    w[120] = 0.00168114286542146990631;
    w[121] = 0.00128952408261041739210;
    w[122] = 0.000938369848542381500794;
    w[123] = 0.000632607319362633544219;
    w[124] = 0.000377746646326984660274;
    w[125] = 0.000180739564445388357820;
    w[126] = 0.0000505360952078625176247;
  }
  else if (n==255) {
    w[  0] = 0.69379364324108267170E-05;
    w[  1] = 0.25157870384280661489E-04;
    w[  2] = 0.53275293669780613125E-04;
    w[  3] = 0.90372734658751149261E-04;
    w[  4] = 0.13575491094922871973E-03;
    w[  5] = 0.18887326450650491366E-03;
    w[  6] = 0.24921240048299729402E-03;
    w[  7] = 0.31630366082226447689E-03;
    w[  8] = 0.38974528447328229322E-03;
    w[  9] = 0.46918492424785040975E-03;
    w[ 10] = 0.55429531493037471492E-03;
    w[ 11] = 0.64476204130572477933E-03;
    w[ 12] = 0.74028280424450333046E-03;
    w[ 13] = 0.84057143271072246365E-03;
    w[ 14] = 0.94536151685852538246E-03;
    w[ 15] = 0.10544076228633167722E-02;
    w[ 16] = 0.11674841174299594077E-02;
    w[ 17] = 0.12843824718970101768E-02;
    w[ 18] = 0.14049079956551446427E-02;
    w[ 19] = 0.15288767050877655684E-02;
    w[ 20] = 0.16561127281544526052E-02;
    w[ 21] = 0.17864463917586498247E-02;
    w[ 22] = 0.19197129710138724125E-02;
    w[ 23] = 0.20557519893273465236E-02;
    w[ 24] = 0.21944069253638388388E-02;
    w[ 25] = 0.23355251860571608737E-02;
    w[ 26] = 0.24789582266575679307E-02;
    w[ 27] = 0.26245617274044295626E-02;
    w[ 28] = 0.27721957645934509940E-02;
    w[ 29] = 0.29217249379178197538E-02;
    w[ 30] = 0.30730184347025783234E-02;
    w[ 31] = 0.32259500250878684614E-02;
    w[ 32] = 0.33803979910869203823E-02;
    w[ 33] = 0.35362449977167777340E-02;
    w[ 34] = 0.36933779170256508183E-02;
    w[ 35] = 0.38516876166398709241E-02;
    w[ 36] = 0.40110687240750233989E-02;
    w[ 37] = 0.41714193769840788528E-02;
    w[ 38] = 0.43326409680929828545E-02;
    w[ 39] = 0.44946378920320678616E-02;
    w[ 40] = 0.46573172997568547773E-02;
    w[ 41] = 0.48205888648512683476E-02;
    w[ 42] = 0.49843645647655386012E-02;
    w[ 43] = 0.51485584789781777618E-02;
    w[ 44] = 0.53130866051870565663E-02;
    w[ 45] = 0.54778666939189508240E-02;
    w[ 46] = 0.56428181013844441585E-02;
    w[ 47] = 0.58078616599775673635E-02;
    w[ 48] = 0.59729195655081658049E-02;
    w[ 49] = 0.61379152800413850435E-02;
    w[ 50] = 0.63027734490857587172E-02;
    w[ 51] = 0.64674198318036867274E-02;
    w[ 52] = 0.66317812429018878941E-02;
    w[ 53] = 0.67957855048827733948E-02;
    w[ 54] = 0.69593614093904229394E-02;
    w[ 55] = 0.71224386864583871532E-02;
    w[ 56] = 0.72849479805538070639E-02;
    w[ 57] = 0.74468208324075910174E-02;
    w[ 58] = 0.76079896657190565832E-02;
    w[ 59] = 0.77683877779219912200E-02;
    w[ 60] = 0.79279493342948491103E-02;
    w[ 61] = 0.80866093647888599710E-02;
    w[ 62] = 0.82443037630328680306E-02;
    w[ 63] = 0.84009692870519326354E-02;
    w[ 64] = 0.85565435613076896192E-02;
    w[ 65] = 0.87109650797320868736E-02;
    w[ 66] = 0.88641732094824942641E-02;
    w[ 67] = 0.90161081951956431600E-02;
    w[ 68] = 0.91667111635607884067E-02;
    w[ 69] = 0.93159241280693950932E-02;
    w[ 70] = 0.94636899938300652943E-02;
    w[ 71] = 0.96099525623638830097E-02;
    w[ 72] = 0.97546565363174114611E-02;
    w[ 73] = 0.98977475240487497440E-02;
    w[ 74] = 0.10039172044056840798E-01;
    w[ 75] = 0.10178877529236079733E-01;
    w[ 76] = 0.10316812330947621682E-01;
    w[ 77] = 0.10452925722906011926E-01;
    w[ 78] = 0.10587167904885197931E-01;
    w[ 79] = 0.10719490006251933623E-01;
    w[ 80] = 0.10849844089337314099E-01;
    w[ 81] = 0.10978183152658912470E-01;
    w[ 82] = 0.11104461134006926537E-01;
    w[ 83] = 0.11228632913408049354E-01;
    w[ 84] = 0.11350654315980596602E-01;
    w[ 85] = 0.11470482114693874380E-01;
    w[ 86] = 0.11588074033043952568E-01;
    w[ 87] = 0.11703388747657003101E-01;
    w[ 88] = 0.11816385890830235763E-01;
    w[ 89] = 0.11927026053019270040E-01;
    w[ 90] = 0.12035270785279562630E-01;
    w[ 91] = 0.12141082601668299679E-01;
    w[ 92] = 0.12244424981611985899E-01;
    w[ 93] = 0.12345262372243838455E-01;
    w[ 94] = 0.12443560190714035263E-01;
    w[ 95] = 0.12539284826474884353E-01;
    w[ 96] = 0.12632403643542078765E-01;
    w[ 97] = 0.12722884982732382906E-01;
    w[ 98] = 0.12810698163877361967E-01;
    w[ 99] = 0.12895813488012114694E-01;
    w[100] = 0.12978202239537399286E-01;
    w[101] = 0.13057836688353048840E-01;
    w[102] = 0.13134690091960152836E-01;
    w[103] = 0.13208736697529129966E-01;
    w[104] = 0.13279951743930530650E-01;
    w[105] = 0.13348311463725179953E-01;
    w[106] = 0.13413793085110098513E-01;
    w[107] = 0.13476374833816515982E-01;
    w[108] = 0.13536035934956213614E-01;
    w[109] = 0.13592756614812395910E-01;
    w[110] = 0.13646518102571291428E-01;
    w[111] = 0.13697302631990716258E-01;
    w[112] = 0.13745093443001896632E-01;
    w[113] = 0.13789874783240936517E-01;
    w[114] = 0.13831631909506428676E-01;
    w[115] = 0.13870351089139840997E-01;
    w[116] = 0.13906019601325461264E-01;
    w[117] = 0.13938625738306850804E-01;
    w[118] = 0.13968158806516938516E-01;
    w[119] = 0.13994609127619079852E-01;
    w[120] = 0.14017968039456608810E-01;
    w[121] = 0.14038227896908623303E-01;
    w[122] = 0.14055382072649964277E-01;
    w[123] = 0.14069424957813575318E-01;
    w[124] = 0.14080351962553661325E-01;
    w[125] = 0.14088159516508301065E-01;
    w[126] = 0.14092845069160408355E-01;
    w[127] = 0.14094407090096179347E-01;
    w[128] = 0.14092845069160408355E-01;
    w[129] = 0.14088159516508301065E-01;
    w[130] = 0.14080351962553661325E-01;
    w[131] = 0.14069424957813575318E-01;
    w[132] = 0.14055382072649964277E-01;
    w[133] = 0.14038227896908623303E-01;
    w[134] = 0.14017968039456608810E-01;
    w[135] = 0.13994609127619079852E-01;
    w[136] = 0.13968158806516938516E-01;
    w[137] = 0.13938625738306850804E-01;
    w[138] = 0.13906019601325461264E-01;
    w[139] = 0.13870351089139840997E-01;
    w[140] = 0.13831631909506428676E-01;
    w[141] = 0.13789874783240936517E-01;
    w[142] = 0.13745093443001896632E-01;
    w[143] = 0.13697302631990716258E-01;
    w[144] = 0.13646518102571291428E-01;
    w[145] = 0.13592756614812395910E-01;
    w[146] = 0.13536035934956213614E-01;
    w[147] = 0.13476374833816515982E-01;
    w[148] = 0.13413793085110098513E-01;
    w[149] = 0.13348311463725179953E-01;
    w[150] = 0.13279951743930530650E-01;
    w[151] = 0.13208736697529129966E-01;
    w[152] = 0.13134690091960152836E-01;
    w[153] = 0.13057836688353048840E-01;
    w[154] = 0.12978202239537399286E-01;
    w[155] = 0.12895813488012114694E-01;
    w[156] = 0.12810698163877361967E-01;
    w[157] = 0.12722884982732382906E-01;
    w[158] = 0.12632403643542078765E-01;
    w[159] = 0.12539284826474884353E-01;
    w[160] = 0.12443560190714035263E-01;
    w[161] = 0.12345262372243838455E-01;
    w[162] = 0.12244424981611985899E-01;
    w[163] = 0.12141082601668299679E-01;
    w[164] = 0.12035270785279562630E-01;
    w[165] = 0.11927026053019270040E-01;
    w[166] = 0.11816385890830235763E-01;
    w[167] = 0.11703388747657003101E-01;
    w[168] = 0.11588074033043952568E-01;
    w[169] = 0.11470482114693874380E-01;
    w[170] = 0.11350654315980596602E-01;
    w[171] = 0.11228632913408049354E-01;
    w[172] = 0.11104461134006926537E-01;
    w[173] = 0.10978183152658912470E-01;
    w[174] = 0.10849844089337314099E-01;
    w[175] = 0.10719490006251933623E-01;
    w[176] = 0.10587167904885197931E-01;
    w[177] = 0.10452925722906011926E-01;
    w[178] = 0.10316812330947621682E-01;
    w[179] = 0.10178877529236079733E-01;
    w[180] = 0.10039172044056840798E-01;
    w[181] = 0.98977475240487497440E-02;
    w[182] = 0.97546565363174114611E-02;
    w[183] = 0.96099525623638830097E-02;
    w[184] = 0.94636899938300652943E-02;
    w[185] = 0.93159241280693950932E-02;
    w[186] = 0.91667111635607884067E-02;
    w[187] = 0.90161081951956431600E-02;
    w[188] = 0.88641732094824942641E-02;
    w[189] = 0.87109650797320868736E-02;
    w[190] = 0.85565435613076896192E-02;
    w[191] = 0.84009692870519326354E-02;
    w[192] = 0.82443037630328680306E-02;
    w[193] = 0.80866093647888599710E-02;
    w[194] = 0.79279493342948491103E-02;
    w[195] = 0.77683877779219912200E-02;
    w[196] = 0.76079896657190565832E-02;
    w[197] = 0.74468208324075910174E-02;
    w[198] = 0.72849479805538070639E-02;
    w[199] = 0.71224386864583871532E-02;
    w[200] = 0.69593614093904229394E-02;
    w[201] = 0.67957855048827733948E-02;
    w[202] = 0.66317812429018878941E-02;
    w[203] = 0.64674198318036867274E-02;
    w[204] = 0.63027734490857587172E-02;
    w[205] = 0.61379152800413850435E-02;
    w[206] = 0.59729195655081658049E-02;
    w[207] = 0.58078616599775673635E-02;
    w[208] = 0.56428181013844441585E-02;
    w[209] = 0.54778666939189508240E-02;
    w[210] = 0.53130866051870565663E-02;
    w[211] = 0.51485584789781777618E-02;
    w[212] = 0.49843645647655386012E-02;
    w[213] = 0.48205888648512683476E-02;
    w[214] = 0.46573172997568547773E-02;
    w[215] = 0.44946378920320678616E-02;
    w[216] = 0.43326409680929828545E-02;
    w[217] = 0.41714193769840788528E-02;
    w[218] = 0.40110687240750233989E-02;
    w[219] = 0.38516876166398709241E-02;
    w[220] = 0.36933779170256508183E-02;
    w[221] = 0.35362449977167777340E-02;
    w[222] = 0.33803979910869203823E-02;
    w[223] = 0.32259500250878684614E-02;
    w[224] = 0.30730184347025783234E-02;
    w[225] = 0.29217249379178197538E-02;
    w[226] = 0.27721957645934509940E-02;
    w[227] = 0.26245617274044295626E-02;
    w[228] = 0.24789582266575679307E-02;
    w[229] = 0.23355251860571608737E-02;
    w[230] = 0.21944069253638388388E-02;
    w[231] = 0.20557519893273465236E-02;
    w[232] = 0.19197129710138724125E-02;
    w[233] = 0.17864463917586498247E-02;
    w[234] = 0.16561127281544526052E-02;
    w[235] = 0.15288767050877655684E-02;
    w[236] = 0.14049079956551446427E-02;
    w[237] = 0.12843824718970101768E-02;
    w[238] = 0.11674841174299594077E-02;
    w[239] = 0.10544076228633167722E-02;
    w[240] = 0.94536151685852538246E-03;
    w[241] = 0.84057143271072246365E-03;
    w[242] = 0.74028280424450333046E-03;
    w[243] = 0.64476204130572477933E-03;
    w[244] = 0.55429531493037471492E-03;
    w[245] = 0.46918492424785040975E-03;
    w[246] = 0.38974528447328229322E-03;
    w[247] = 0.31630366082226447689E-03;
    w[248] = 0.24921240048299729402E-03;
    w[249] = 0.18887326450650491366E-03;
    w[250] = 0.13575491094922871973E-03;
    w[251] = 0.90372734658751149261E-04;
    w[252] = 0.53275293669780613125E-04;
    w[253] = 0.25157870384280661489E-04;
    w[254] = 0.69379364324108267170E-05;
  }
  else {
    std::cerr << "\n";
    std::cerr << "PATTERSON_LOOKUP_WEIGHTS - Fatal error!\n";
    std::cerr << "  Unexpected value of N = " << n << ".\n";
    std::exit(1);
  }
  return;
}

//***************************************************************************
template<class Scalar>
void IntrepidBurkardtRules::trapezoidal_compute ( int n, Scalar x[], Scalar w[] ) 
//***************************************************************************
{
  if (n==1) {
    x[0] = 0.0;
    w[0] = 2.0;
  }
  else {
    Scalar h = 1.0/((Scalar)n-1.0);
    for (int i=0; i<n; i++) {
      x[i] = -1.0 + (Scalar)i*h*2.0;
      if (i==0||i==n-1) {
	w[i] = h;
      }
      else {
	w[i] = 2.0*h;
      }
    }
  }
  return;
}

//***************************************************************************
template<class Scalar>
void IntrepidBurkardtRules::trapezoidal_compute_points ( int n, Scalar x[] ) 
//***************************************************************************
{
  if (n==1) {
    x[0] = 0.0;
  }
  else {
    Scalar h = 1.0/((Scalar)n-1.0);
    for (int i=0; i<n; i++) {
      x[i] = -1.0 + (Scalar)i*h*2.0;
    }
  }
  return;
}

//***************************************************************************
template<class Scalar>
void IntrepidBurkardtRules::trapezoidal_compute_weights ( int n, Scalar w[] ) 
//***************************************************************************
{
  if (n==1) {
    w[0] = 2.0;
  }
  else {
    Scalar h = 1.0/((Scalar)n-1.0);
    for (int i=0; i<n; i++) {
      if (i==0||i==n-1) {
	w[i] = h;
      }
      else {
	w[i] = 2.0*h;
      }
    }
  }
  return;
}

//****************************************************************************
template<class Scalar> 
Scalar IntrepidBurkardtRules::r8_epsilon(Scalar one)
//****************************************************************************
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, Scalar R8_EPSILON, the R8 round-off unit.
//
{
  Scalar value; value = 1.0;

  while (1.0<(Scalar)(1.0+value)) {
    value = value / 2.0;
  }

  value = 2.0 * value;

  return value;
}

//****************************************************************************
template<class Scalar> 
Scalar IntrepidBurkardtRules::r8_sign ( Scalar x )
//****************************************************************************
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, Scalar X, the number whose sign is desired.
//
//    Output, Scalar R8_SIGN, the sign of X.
//
{
  Scalar value;

  if (x<0.0) {
    value = -1.0;
  }
  else {
    value = 1.0;
  }
  return value;
}

} // end of namespace Intrepid

