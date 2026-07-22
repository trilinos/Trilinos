// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// high_order_example
//
// usage:
//   high_order_example
//
// output:
//   computes a high-order derivative of a simple function along a small set
//   of directions/modes using DFad, and then computes the same thing through
//   the multivariate chain rule.

#include <iostream>

#include "Sacado.hpp"

typedef Sacado::Fad::DFad<double> FadType;

// Function to compute derivatives of
template <typename Scalar>
Scalar func(const int m, const Scalar x[]) {
  using std::exp;
  Scalar z = 1.0;
  for (int i=0; i<m; ++i)
    z *= x[i];
  return exp(z);
}

// Build nested Fad type for computing derivatives up to order N
template <int N>
struct MakeFad {
  // Replace double in FadType with MakeFad<N-1>::type
  typedef typename MakeFad<N-1>::type nested_type;
  typedef typename Sacado::mpl::apply<FadType,nested_type>::type type;

  // Initialize Fad for computation of full derivative
  static type apply(const int n, const int i, const double x) {
    return type(n,i,MakeFad<N-1>::apply(n,i,x));
  }

  // Initialize Fad object for derivative-matrix product
  static type apply(const int n, const double x, const double v[]) {
    type x_fad(n,MakeFad<N-1>::apply(n,x,v));
    for (int i=0; i<n; ++i)
      x_fad.fastAccessDx(i) = v[i];
    return x_fad;
  }
};
template <>
struct MakeFad<1> {
  typedef FadType type;

  // Initialize Fad for computation of full derivative
  static type apply(const int n, const int i, const double x) {
    return type(n,i,x);
  }

  // Initialize Fad object for derivative-matrix product
  static type apply(const int n, const double x, const double v[]) {
    type x_fad(n,x);
    for (int i=0; i<n; ++i)
      x_fad.fastAccessDx(i) = v[i];
    return x_fad;
  }
};

// Extract d^r/(dx_k_1 ... dx_k_r) where the sequence k_1 ... k_r is indicated
// by the passed iterator (e.g., iterator into a std::vector)
template <typename TermIterator>
double
extract_derivative(const double& x, TermIterator term, TermIterator term_end)
{
  return x;
}
template <typename T, typename TermIterator>
double
extract_derivative(const T& x, TermIterator term, TermIterator term_end)
{
  // If there are no terms, return value
  if (term == term_end)
    return Sacado::ScalarValue<T>::eval(x);

  // Get the first term
  auto k = *term;

  // Extract remaining terms
  return extract_derivative(x.fastAccessDx(k), ++term, term_end);
}

int main(int argc, char **argv)
{
  const int deg = 6; // Order of derivative to compute
  const int m = 3;   // Number of coordinates
  const int n = 2;   // Number of modes
  const double x0[m] = { 1.0, 1.0, 1.0 };    // Expansion point
  const double V[m][n] = { {0.1,  0.2},
                           {0.3,  0.4},
                           {0.5,  0.6} };    // Mode coordinates

  // Derivative term we wish to extract.  This corresponds to the derivative
  // d^6/(ds_0 ds_0 ds_0 ds_1 ds_1 ds_1) = d^6/(ds_0^3 ds_1^3)
  std::vector<int> term = { 0, 0, 0, 1, 1, 1 };

  // For y = f(x(s)), x(s) = x0 + V*s, compute
  // d^r y / (ds_i_1 ... ds_i_r) where (i_1,...,i_r) is indicated by term
  typedef typename MakeFad<deg>::type NestedFadType;
  NestedFadType x_fad[m];
  for (int i=0; i<m; ++i)
    x_fad[i] = MakeFad<deg>::apply(n,x0[i],V[i]);
  const NestedFadType f_fad = func(m,x_fad);
  const double z1 = extract_derivative(f_fad, term.begin(), term.end());

  // Now compute the same thing by first computing
  // d^k y / (dx_j_1 .. dx_j_k) for 0 <= k <= deg and j_1,...,j_k = 0,...,m-1
  NestedFadType x_fad2[m];
  for (int i=0; i<m; ++i)
    x_fad2[i] = MakeFad<deg>::apply(m,i,x0[i]);
  const NestedFadType f_fad2 = func(m,x_fad2);

  // Now compute multivariate chain-rule:
  // d^r y / (ds_i_1 ... ds_i_r) =
  //  \sum_{j_1,...,j_r=0}^{m-1} d^r y/(dx_j_1 .. dx_j_r) *
  //     V[j_1][i_1]...V[j_r][i_r].
  // This requires iterating (j_1,...,j_r) over the m^r-dimensional tensor
  // product space [0,m-1] x ... x [0,m-1]
  double z2 = 0.0;
  const int r = term.size();
  std::vector<int> t(r, 0);
  bool finished = false;
  while (!finished) {
    const double z = extract_derivative(f_fad2, t.begin(), t.end());
    double c = 1.0;
    for (int i=0; i<r; ++i) {
      c *= V[t[i]][term[i]];
    }
    z2 += z*c;

    // Increment t to next term in tensor product space
    ++t[0];
    int j=0;
    while (j<r-1 && t[j] >= m) {
      t[j] = 0;
      ++j;
      ++t[j];
    }
    if (t[r-1] == m)
      finished = true;
  }

  const double error = std::abs(z1-z2)/std::abs(z2);
  std::cout << "z (reduced) = " << z1 << " z (full) = " << z2
            << " error = " << error << std::endl;

  const double tol = 1.0e-14;
  if (error < tol) {
    std::cout << "\nExample passed!" << std::endl;
    return 0;
  }

  std::cout <<"\nSomething is wrong, example failed!" << std::endl;
  return 1;
}
