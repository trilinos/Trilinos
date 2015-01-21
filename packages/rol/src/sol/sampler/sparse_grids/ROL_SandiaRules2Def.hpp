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


namespace ROL {

//****************************************************************************80
void SandiaRules2::ccn_points ( int n, int dim, double x[] )
//****************************************************************************80
//
//  Purpose:
//
//    CCN_POINTS computes nested Clenshaw Curtis quadrature points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double X[N], the abscissas.
//
{
  SandiaRules::ccn_compute_points ( n, x );
  return;
}

//****************************************************************************80
void SandiaRules2::ccn_weights ( int n, int dim, double w[] )
//****************************************************************************80
//
//  Purpose:
//
//    CCN_WEIGHTS computes nested Clenshaw Curtis quadrature weights.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double W[N], the weights.
//
{
  SandiaRules::ccn_compute_weights ( n, w );
  return;
}

//****************************************************************************80
void SandiaRules2::clenshaw_curtis_points ( int n, int dim, double x[] )
//****************************************************************************80
//
//  Purpose:
//
//    CLENSHAW_CURTIS_POINTS computes Clenshaw Curtis quadrature points.
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
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double X[N], the abscissas.
//
{
  int index;
  double pi = 3.141592653589793;

  if ( n < 1 )
  {
    std::cerr << "\n";
    std::cerr << "CLENSHAW_CURTIS_POINTS - Fatal error!\n";
    std::cerr << "  Order is less than 1.\n";
    std::exit ( 1 );
  }
  else if ( n == 1 )
  {
    x[0] = 0.0;
  }
  else
  {
    for ( index = 1; index <= n; index++ )
    {
      x[index-1] =  std::cos ( ( double ) ( n - index ) * pi
                             / ( double ) ( n - 1     ) );
    }
    x[0] = -1.0;
    if ( ( n % 2 ) == 1 )
    {
      x[(n-1)/2] = 0.0;
    }
    x[n-1] = +1.0;
  }
  return;
}

//****************************************************************************80
void SandiaRules2::clenshaw_curtis_weights ( int n, int dim, double w[] )
//****************************************************************************80
//
//  Purpose:
//
//    CLENSHAW_CURTIS_WEIGHTS computes Clenshaw Curtis quadrature weights.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double W[N], the weights.
//
{
  double b;
  int i;
  int j;
  double pi = 3.141592653589793;
  double theta;

  if ( n < 1 )
  {
    std::cerr << "\n";
    std::cerr << "CLENSHAW_CURTIS_WEIGHTS - Fatal error!\n";
    std::cerr << "  Order is less than 1.\n";
    std::exit ( 1 );
  }
  else if ( n == 1 )
  {
    w[0] = 2.0;
    return;
  }

  for ( i = 1; i <= n; i++ )
  {
    theta = ( double ) ( i - 1 ) * pi / ( double ) ( n - 1 );

    w[i-1] = 1.0;

    for ( j = 1; j <= ( n - 1 ) / 2; j++ )
    {
      if ( 2 * j == ( n - 1 ) )
      {
        b = 1.0;
      }
      else
      {
        b = 2.0;
      }

      w[i-1] = w[i-1] - b *  std::cos ( 2.0 * ( double ) ( j ) * theta )
           / ( double ) ( 4 * j * j - 1 );
    }
  }

  w[0] = w[0] / ( double ) ( n - 1 );
  for ( i = 1; i < n - 1; i++ )
  {
    w[i] = 2.0 * w[i] / ( double ) ( n - 1 );
  }
  w[n-1] = w[n-1] / ( double ) ( n - 1 );

  return;
}

//****************************************************************************80
void SandiaRules2::fejer2_points ( int n, int dim, double x[] )
//****************************************************************************80
//
//  Purpose:
//
//    FEJER2_POINTS computes Fejer type 2 quadrature points.
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
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double X[N], the abscissas.
//
{
  int index;
  double pi = 3.141592653589793;

  if ( n < 1 )
  {
    std::cerr << "\n";
    std::cerr << "FEJER2_POINTS - Fatal error!\n";
    std::cerr << "  Order is less than 1.\n";
    std::exit ( 1 );
  }
  else if ( n == 1 )
  {
    x[0] = 0.0;
  }
  else
  {
    for ( index = 1; index <= n; index++ )
    {
      x[index-1] =  std::cos ( ( double ) ( n + 1 - index ) * pi
                             / ( double ) ( n + 1 ) );
    }
    if ( ( n % 2 ) == 1 )
    {
      x[(n-1)/2] = 0.0;
    }
  }
  return;
}

//****************************************************************************80
void SandiaRules2::fejer2_weights ( int n, int dim, double w[] )
//****************************************************************************80
//
//  Purpose:
//
//    FEJER2_WEIGHTS computes Fejer type 2 quadrature weights.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double W[N], the weights.
//
{
  int i;
  int j;
  double p;
  double pi = 3.141592653589793;
  double theta;

  if ( n < 1 )
  {
    std::cerr << "\n";
    std::cerr << "FEJER2_WEIGHTS - Fatal error!\n";
    std::cerr << "  N < 1.\n";
    std::exit ( 1 );
  }
  else if ( n == 1 )
  {
    w[0] = 2.0;
  }
  else if ( n == 2 )
  {
    w[0] = 1.0;
    w[1] = 1.0;
  }
  else
  {
    for ( i = 1; i <= n; i++ )
    {
      theta = ( double ) ( n + 1 - i ) * pi
            / ( double ) ( n + 1 );

      w[i-1] = 1.0;

      for ( j = 1; j <= ( ( n - 1 ) / 2 ); j++ )
      {
        w[i-1] = w[i-1] - 2.0 *  std::cos ( 2.0 * ( double ) ( j ) * theta )
          / ( double ) ( 4 * j * j - 1 );
      }
      p = 2.0 * ( double ) ( ( ( n + 1 ) / 2 ) ) - 1.0;
      w[i-1] = w[i-1] -  std::cos ( ( p + 1.0 ) * theta ) / p;
    }
    for ( i = 0; i < n; i++ )
    {
      w[i] = 2.0 * w[i] / ( double ) ( n + 1 );
    }
  }
  return;
}

//****************************************************************************80
void SandiaRules2::gen_hermite_points ( int n, int dim, double x[] )
//****************************************************************************80
//
//  Purpose:
//
//    GEN_HERMITE_POINTS: Generalized Hermite quadrature points.
//
//  Discussion:
//
//    This function assumes the existence of a function:
//      double parameter ( int dim, int offset )
//    which can supply the value of the ALPHA parameter.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double X[N], the abscissas.
//
{
  double alpha;

  alpha = parameter ( dim, 0 );

  SandiaRules::gen_hermite_compute_points ( n, alpha, x );

  return;
}

//****************************************************************************80
void SandiaRules2::gen_hermite_weights ( int n, int dim, double w[] )
//****************************************************************************80
//
//  Purpose:
//
//    GEN_HERMITE_WEIGHTS: Generalized Hermite quadrature weights.
//
//  Discussion:
//
//    This function assumes the existence of a function:
//      double parameter ( int dim, int offset )
//    which can supply the value of the ALPHA parameter.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double W[N], the weights.
//
{
  double alpha;
  alpha = parameter ( dim, 0 );
  SandiaRules::gen_hermite_compute_weights ( n, alpha, w );
  return;
}

//****************************************************************************80
void SandiaRules2::gen_laguerre_points ( int n, int dim, double x[] )
//****************************************************************************80
//
//  Purpose:
//
//    GEN_LAGUERRE_POINTS: Generalized Laguerre quadrature points.
//
//  Discussion:
//
//    This function assumes the existence of a function:
//      double parameter ( int dim, int offset )
//    which can supply the value of the ALPHA parameter.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double X[N], the abscissas.
//
{
  double alpha;
  alpha = parameter ( dim, 0 );
  SandiaRules::gen_laguerre_compute_points ( n, alpha, x );
  return;
}

//****************************************************************************80
void SandiaRules2::gen_laguerre_weights ( int n, int dim, double w[] )
//****************************************************************************80
//
//  Purpose:
//
//    GEN_LAGUERRE_WEIGHTS: Generalized Laguerre quadrature weights.
//
//  Discussion:
//
//    This function assumes the existence of a function:
//      double parameter ( int dim, int offset )
//    which can supply the value of the ALPHA parameter.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double W[N], the weights.
//
{
  double alpha;
  alpha = parameter ( dim, 0 );
  SandiaRules::gen_laguerre_compute_weights ( n, alpha, w );
  return;
}

//****************************************************************************80
void SandiaRules2::hcc_points ( int n, int dim, double x[] )
//****************************************************************************80
//
//  Purpose:
//
//    HCC_POINTS computes Hermite-Cubic-Chebyshev-Spacing quadrature points.
//
//  Discussion:
//
//    For the HCC rule, we assume that an interval has been divided by
//    M nodes X into Chebyshev-spaced subintervals, and that at each
//    abscissa both function and derivative information is available.
//    The piecewise cubic Hermite interpolant is constructed for this data.
//    The quadrature rule uses N = 2 * M abscissas, where each node is
//    listed twice, and the weights occur in pairs, with the first multiplying
//    the function value and the second the derivative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double X[N], the abscissas.
//
{
  SandiaRules::hcc_compute_points ( n, x );
  return;
}

//****************************************************************************80
void SandiaRules2::hcc_weights ( int n, int dim, double w[] )
//****************************************************************************80
//
//  Purpose:
//
//    HCC_WEIGHTS computes Hermite-Cubic-Chebyshev-Spacing quadrature weights.
//
//  Discussion:
//
//    For the HCC rule, we assume that an interval has been divided by
//    M nodes X into Chebyshev-spaced subintervals, and that at each
//    abscissa both function and derivative information is available.
//    The piecewise cubic Hermite interpolant is constructed for this data.
//    The quadrature rule uses N = 2 * M abscissas, where each node is
//    listed twice, and the weights occur in pairs, with the first multiplying
//    the function value and the second the derivative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double W[N], the weights.
//
{
  SandiaRules::hcc_compute_weights ( n, w );
  return;
}

//****************************************************************************80
void SandiaRules2::hce_points ( int n, int dim, double x[] )
//****************************************************************************80
//
//  Purpose:
//
//    HCE_POINTS computes Hermite-Cubic-Equal-Spacing quadrature points.
//
//  Discussion:
//
//    For the HCE rule, we assume that an interval has been divided by
//    M nodes X into equally spaced subintervals, and that at each
//    abscissa both function and derivative information is available.
//    The piecewise cubic Hermite interpolant is constructed for this data.
//    The quadrature rule uses N = 2 * M abscissas, where each node is
//    listed twice, and the weights occur in pairs, with the first multiplying
//    the function value and the second the derivative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double X[N], the abscissas.
//
{
  SandiaRules::hce_compute_points ( n, x );
  return;
}

//****************************************************************************80
void SandiaRules2::hce_weights ( int n, int dim, double w[] )
//****************************************************************************80
//
//  Purpose:
//
//    HCE_WEIGHTS computes Hermite-Cubic-Equal-Spacing quadrature weights.
//
//  Discussion:
//
//    For the HCE rule, we assume that an interval has been divided by
//    M nodes X into equally spaced subintervals, and that at each
//    abscissa both function and derivative information is available.
//    The piecewise cubic Hermite interpolant is constructed for this data.
//    The quadrature rule uses N = 2 * M abscissas, where each node is
//    listed twice, and the weights occur in pairs, with the first multiplying
//    the function value and the second the derivative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double W[N], the weights.
//
{
  SandiaRules::hce_compute_weights ( n, w );
  return;
}

//****************************************************************************80
void SandiaRules2::hermite_genz_keister_points ( int n, int dim, double x[] )
//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_GENZ_KEISTER_POINTS looks up Genz-Keister Hermite abscissas.
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
//  Parameters:
//
//    Input, int N, the order.
//    The allowed orders are 1, 3, 9, 19 and 35.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double X[N], the abscissas.
//
{
  SandiaRules::hermite_genz_keister_lookup_points ( n, x );
  return;
}

//****************************************************************************80
void SandiaRules2::hermite_genz_keister_weights ( int n, int dim, double w[] )
//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_GENZ_KEISTER_WEIGHTS looks up Genz-Keister Hermite weights.
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
//  Parameters:
//
//    Input, int N, the order.
//    The allowed orders are 1, 3, 9, 19 and 35.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double W[N], the weights.
//
{
  SandiaRules::hermite_genz_keister_lookup_weights ( n, w );
  return;
}

//****************************************************************************80
void SandiaRules2::hermite_points ( int n, int dim, double x[] )
//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POINTS computes Hermite quadrature points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double X[N], the abscissas.
//
{
  SandiaRules::hermite_compute_points ( n, x );
  return;
}

//****************************************************************************80
void SandiaRules2::hermite_weights ( int n, int dim, double w[] )
//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_WEIGHTS computes Hermite quadrature weights.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double W[N], the weights.
//
{
  SandiaRules::hermite_compute_weights ( n, w );
  return;
}

//****************************************************************************80
void SandiaRules2::jacobi_points ( int n, int dim, double x[] )
//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_POINTS computes Jacobi quadrature points.
//
//  Discussion:
//
//    This function assumes the existence of a function:
//      double parameter ( int dim, int offset )
//    which can supply the values of the ALPHA and BETA parameters.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double X[N], the abscissas.
//
{
  double alpha;
  double beta;

  alpha = parameter ( dim, 0 );
  beta  = parameter ( dim, 1 );

  SandiaRules::jacobi_compute_points ( n, alpha, beta, x );

  return;
}

//****************************************************************************80
void SandiaRules2::jacobi_weights ( int n, int dim, double w[] )
//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_WEIGHTS computes Jacobi quadrature weights.
//
//  Discussion:
//
//    This function assumes the existence of a function:
//      double parameter ( int dim, int offset )
//    which can supply the values of the ALPHA and BETA parameters.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double W[N], the weights.
//
{
  double alpha;
  double beta;

  alpha = parameter ( dim, 0 );
  beta  = parameter ( dim, 1 );

  SandiaRules::jacobi_compute_weights ( n, alpha, beta, w );

  return;
}

//****************************************************************************80
void SandiaRules2::laguerre_points ( int n, int dim, double x[] )
//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_POINTS computes Laguerre quadrature points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double X[N], the abscissas.
//
{
  laguerre_compute_points ( n, x );
  return;
}

//****************************************************************************80
void SandiaRules2::laguerre_weights ( int n, int dim, double w[] )
//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_WEIGHTS computes Laguerre quadrature weights.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double W[N], the weights.
//
{
  SandiaRules::laguerre_compute_weights ( n, w );
  return;
}

//****************************************************************************80
void SandiaRules2::legendre_points ( int n, int dim, double x[] )
//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POINTS computes Legendre quadrature points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double X[N], the abscissas.
//
{
  SandiaRules::legendre_compute_points ( n, x );
  return;
}

//****************************************************************************80
void SandiaRules2::legendre_weights ( int n, int dim, double w[] )
//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_WEIGHTS computes Legendre quadrature weights.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double W[N], the weights.
//
{
  SandiaRules::legendre_compute_weights ( n, w );
  return;
}

//****************************************************************************80
void SandiaRules2::ncc_points ( int n, int dim, double x[] )
//****************************************************************************80
//
//  Purpose:
//
//    NCC_POINTS computes Newton Cotes Closed quadrature points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double X[N], the abscissas.
//
{
  SandiaRules::ncc_compute_points ( n, x );
  return;
}

//****************************************************************************80
void SandiaRules2::ncc_weights ( int n, int dim, double w[] )
//****************************************************************************80
//
//  Purpose:
//
//    NCC_WEIGHTS computes Newton Cotes Closed quadrature weights.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double W[N], the weights.
//
{
  SandiaRules::ncc_compute_weights ( n, w );
}

//****************************************************************************80
void SandiaRules2::nco_points ( int n, int dim, double x[] )
//****************************************************************************80
//
//  Purpose:
//
//    NCO_POINTS computes Newton Cotes Open quadrature points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double X[N], the abscissas.
//
{
  SandiaRules::nco_compute_points ( n, x );
  return;
}

//****************************************************************************80
void SandiaRules2::nco_weights ( int n, int dim, double w[] )
//****************************************************************************80
//
//  Purpose:
//
//    NCO_WEIGHTS computes Newton Cotes Open quadrature weights.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double W[N], the weights.
//
{
  SandiaRules::nco_compute_weights ( n, w );
  return;
}

//****************************************************************************80
void SandiaRules2::patterson_points ( int n, int dim, double x[] )
//****************************************************************************80
//
//  Purpose:
//
//    PATTERSON_POINTS looks up Patterson quadrature points.
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
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double X[N], the abscissas.
//
{
  SandiaRules::patterson_lookup_points ( n, x );
  return;
}

//****************************************************************************80
void SandiaRules2::patterson_weights ( int n, int dim, double w[] )
//****************************************************************************80
//
//  Purpose:
//
//    PATTERSON_WEIGHTS looks up Patterson quadrature weights.
//
//  Discussion:
//
//    The allowed orders are 1, 3, 7, 15, 31, 63 and 127.
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
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int DIM, the spatial dimension represented by this rule,
//    in cases where a multidimensional product rule is being formed.
//
//    Output, double W[N], the weights.
//
{
  SandiaRules::patterson_lookup_weights ( n, w );
  return;
}

//****************************************************************************80
double SandiaRules2::parameter ( int dim, int offset )
//****************************************************************************80
//
//  Purpose:
//
//    PARAMETER is a user-supplied routine to retrieve parameters.
//
//  Discussion:
//
//    The declaration for this function is in SANDIA_RULES.H
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM, the spatial dimension.
//
//    Input, int OFFSET, the offset of the parameter within the 
//    spatial dimension.
//
//    Output, double PARAMETER, the value of the OFFSET-th parameter
//    associated with the DIM-th dimension.
//
{
  return 0.0;
}

} // namespace ROL
