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
// Questions: Alejandro Mota (amota@sandia.gov)
//
// ************************************************************************
// @HEADER

#include "gtest/gtest.h"
#include "Intrepid_MiniTensor_Solvers.h"

int
main(int ac, char * av[])
{
  Kokkos::initialize();

  ::testing::InitGoogleTest(&ac, av);

  auto const
  retval = RUN_ALL_TESTS();

  Kokkos::finalize();

  return retval;
}

namespace Intrepid {

namespace {
//
// Define some nonlinear systems (NLS) to test nonlinear solution methods.
//

//
//
//
template<typename S>
class SquareRootNLS : public Intrepid::Function_Base<SquareRootNLS<S>, S>
{
public:
  SquareRootNLS(S const c) : c_(c)
  {
  }

  static constexpr
  Intrepid::Index
  DIMENSION{1};

  static constexpr
  char const * const
  NAME{"Square Root"};

  // Default value.
  template<typename T, Intrepid::Index N>
  T
  value(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<SquareRootNLS<S>, S>::value(*this, x);
  }

  // Explicit gradient.
  template<typename T, Intrepid::Index N>
  Intrepid::Vector<T, N>
  gradient(Intrepid::Vector<T, N> const & x) const
  {
    Intrepid::Index const
    dimension = x.get_dimension();

    assert(dimension == DIMENSION);

    Intrepid::Vector<T, N>
    r(dimension);

    r(0) = x(0) * x(0) - c_;

    return r;
  }

  // Default AD hessian.
  template<typename T, Intrepid::Index N>
  Intrepid::Tensor<T, N>
  hessian(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<SquareRootNLS<S>, S>::hessian(*this, x);
  }

private:
  S const
  c_{0.0};
};

//
//
//
template<typename S>
class QuadraticNLS : public Intrepid::Function_Base<QuadraticNLS<S>, S>
{
public:
  QuadraticNLS(S const a, S const b, S const c) :  a_(a), b_(b), c_(c)
  {
  }

  static constexpr
  Intrepid::Index
  DIMENSION{2};

  static constexpr
  char const * const
  NAME{"Quadratic"};

  // Default value.
  template<typename T, Intrepid::Index N>
  T
  value(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<QuadraticNLS<S>, S>::value(*this, x);
  }

  // Explicit gradient.
  template<typename T, Intrepid::Index N>
  Intrepid::Vector<T, N>
  gradient(Intrepid::Vector<T, N> const & x) const
  {
    Intrepid::Index const
    dimension = x.get_dimension();

    assert(dimension == DIMENSION);

    Intrepid::Vector<T, N>
    r(dimension);

    r(0) = 2.0 * c_ * (x(0) - a_);
    r(1) = 2.0 * c_ * (x(1) - b_);

    return r;
  }

  // Default AD hessian.
  template<typename T, Intrepid::Index N>
  Intrepid::Tensor<T, N>
  hessian(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<QuadraticNLS<S>, S>::hessian(*this, x);
  }

private:
  S const
  a_{0.0};

  S const
  b_{0.0};

  S const
  c_{0.0};
};

//
//
//
template<typename S>
class GaussianNLS : public Intrepid::Function_Base<GaussianNLS<S>, S>
{
public:
  GaussianNLS(S const a, S const b, S const c) : a_(a), b_(b), c_(c)
  {
  }

  static constexpr
  Intrepid::Index
  DIMENSION{2};

  static constexpr
  char const * const
  NAME{"Inverted Gaussian"};

  // Default value.
  template<typename T, Intrepid::Index N>
  T
  value(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<GaussianNLS<S>, S>::value(*this, x);
  }

  // Explicit gradient.
  template<typename T, Intrepid::Index N>
  Intrepid::Vector<T, N>
  gradient(Intrepid::Vector<T, N> const & x) const
  {
    Intrepid::Index const
    dimension = x.get_dimension();

    assert(dimension == DIMENSION);

    Intrepid::Vector<T, N>
    r(dimension);

    T const
    xa = (x(0) - a_) * c_;

    T const
    xb = (x(1) - b_) * c_;

    T const
    e = std::exp(- xa * xa - xb * xb);

    r(0) = 2.0 * xa * e * c_ * c_;
    r(1) = 2.0 * xb * e * c_ * c_;

    return r;
  }

  // Default AD hessian.
  template<typename T, Intrepid::Index N>
  Intrepid::Tensor<T, N>
  hessian(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<GaussianNLS<S>, S>::hessian(*this, x);
  }

private:
  S const
  a_{0.0};

  S const
  b_{0.0};

  S const
  c_{0.0};
};

//
//
//
template<typename S>
class BananaNLS : public Intrepid::Function_Base<BananaNLS<S>, S>
{
public:

  BananaNLS()
  {
  }

  static constexpr
  Intrepid::Index
  DIMENSION{2};

  static constexpr
  char const * const
  NAME{"Rosenbrock's Banana"};

  // Default value.
  template<typename T, Intrepid::Index N>
  T
  value(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<BananaNLS<S>, S>::value(*this, x);
  }

  // Explicit gradient.
  template<typename T, Intrepid::Index N>
  Intrepid::Vector<T, N>
  gradient(Intrepid::Vector<T, N> const & x) const
  {
    Intrepid::Index const
    dimension = x.get_dimension();

    assert(dimension == DIMENSION);

    Intrepid::Vector<T, N>
    r(dimension);

    r(0) = 2.0 * (x(0) - 1.0) + 400.0 * x(0) * (x(0) * x(0) - x(1));
    r(1) = 200.0 * (x(1) - x(0) * x(0));

    return r;
  }

  // Default AD hessian.
  template<typename T, Intrepid::Index N>
  Intrepid::Tensor<T, N>
  hessian(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<BananaNLS<S>, S>::hessian(*this, x);
  }

};

//
//
//
template<typename S>
class MatyasNLS : public Intrepid::Function_Base<MatyasNLS<S>, S>
{
public:

  MatyasNLS() {}

  static constexpr
  Intrepid::Index
  DIMENSION{2};

  static constexpr
  char const * const
  NAME{"Matyas"};

  // Default value.
  template<typename T, Intrepid::Index N>
  T
  value(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<MatyasNLS<S>, S>::value(*this, x);
  }

  // Explicit gradient.
  template<typename T, Intrepid::Index N>
  Intrepid::Vector<T, N>
  gradient(Intrepid::Vector<T, N> const & x) const
  {
    Intrepid::Index const
    dimension = x.get_dimension();

    assert(dimension == DIMENSION);

    Intrepid::Vector<T, N>
    r(dimension);

    r(0) = (13.0 * x(0) - 12.0 * x(1)) / 25.0;
    r(1) = (13.0 * x(1) - 12.0 * x(0)) / 25.0;

    return r;
  }

  // Default AD hessian.
  template<typename T, Intrepid::Index N>
  Intrepid::Tensor<T, N>
  hessian(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<MatyasNLS<S>, S>::hessian(*this, x);
  }

};

//
//
//
template<typename S>
class McCormickNLS : public Intrepid::Function_Base<McCormickNLS<S>, S>
{
public:

  McCormickNLS() {}

  static constexpr
  Intrepid::Index
  DIMENSION{2};

  static constexpr
  char const * const
  NAME{"McCormick"};

  // Default value.
  template<typename T, Intrepid::Index N>
  T
  value(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<McCormickNLS<S>, S>::value(*this, x);
  }

  // Explicit gradient.
  template<typename T, Intrepid::Index N>
  Intrepid::Vector<T, N>
  gradient(Intrepid::Vector<T, N> const & x) const
  {
    Intrepid::Index const
    dimension = x.get_dimension();

    assert(dimension == DIMENSION);

    Intrepid::Vector<T, N>
    r(dimension);

    r(0) = std::cos(x(0) + x(1)) + 2.0 * x(0) - 2.0 * x(1) - 1.5;
    r(1) = std::cos(x(0) + x(1)) - 2.0 * x(0) + 2.0 * x(1) + 2.5;

    return r;
  }

  // Default AD hessian.
  template<typename T, Intrepid::Index N>
  Intrepid::Tensor<T, N>
  hessian(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<McCormickNLS<S>, S>::hessian(*this, x);
  }

};

//
//
//
template<typename S>
class StyblinskiTangNLS : public Intrepid::Function_Base<StyblinskiTangNLS<S>, S>
{
public:

  StyblinskiTangNLS() {}

  static constexpr
  Intrepid::Index
  DIMENSION{2};

  static constexpr
  char const * const
  NAME{"Styblinski-Tang"};

  // Default value.
  template<typename T, Intrepid::Index N>
  T
  value(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<StyblinskiTangNLS<S>, S>::value(*this, x);
  }

  // Explicit gradient.
  template<typename T, Intrepid::Index N>
  Intrepid::Vector<T, N>
  gradient(Intrepid::Vector<T, N> const & x) const
  {
    Intrepid::Index const
    dimension = x.get_dimension();

    assert(dimension == DIMENSION);

    Intrepid::Vector<T, N>
    r(dimension);

    r(0) = 2.0 * x(0) * x(0) * x(0) - 16.0 * x(0) + 2.5;
    r(1) = 2.0 * x(1) * x(1) * x(1) - 16.0 * x(1) + 2.5;

    return r;
  }

  // Default AD hessian.
  template<typename T, Intrepid::Index N>
  Intrepid::Tensor<T, N>
  hessian(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<StyblinskiTangNLS<S>, S>::hessian(*this, x);
  }

};

//
// Define some nonlinear functions (NLF) to test nonlinear optimization methods.
//

//
// Paraboloid of revolution
//
template<typename S>
class Paraboloid : public Intrepid::Function_Base<Paraboloid<S>, S>
{
public:

  Paraboloid() {}

  static constexpr
  Intrepid::Index
  DIMENSION{2};

  static constexpr
  char const * const
  NAME{"Paraboloid"};

  // Explicit value.
  template<typename T, Intrepid::Index N>
  T
  value(Intrepid::Vector<T, N> const & x)
  {
    Intrepid::Index const
    dimension = x.get_dimension();

    assert(dimension == DIMENSION);

    T const
    f = (x(0) * x(0) + x(1) * x(1));

    return f;
  }

  // Default AD gradient.
  template<typename T, Intrepid::Index N>
  Intrepid::Vector<T, N>
  gradient(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<Paraboloid<S>, S>::gradient(*this, x);
  }

  // Default AD hessian.
  template<typename T, Intrepid::Index N>
  Intrepid::Tensor<T, N>
  hessian(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<Paraboloid<S>, S>::hessian(*this, x);
  }

};

//
// Beale's function
//
template<typename S>
class Beale : public Intrepid::Function_Base<Beale<S>, S>
{
public:

  Beale() {}

  static constexpr
  Intrepid::Index
  DIMENSION{2};

  static constexpr
  char const * const
  NAME{"Beale"};

  // Explicit value.
  template<typename T, Intrepid::Index N>
  T
  value(Intrepid::Vector<T, N> const & X)
  {
    Intrepid::Index const
    dimension = X.get_dimension();

    assert(dimension == DIMENSION);

    T const &
    x = X(0);

    T const &
    y = X(1);

    T const
    a = 1.5 - x + x * y;

    T const
    b = 2.25 - x + x * y * y;

    T const
    c = 2.625 - x + x * y * y * y;

    T const
    f = a * a + b * b + c * c;

    return f;
  }

  // Default AD gradient.
  template<typename T, Intrepid::Index N>
  Intrepid::Vector<T, N>
  gradient(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<Beale<S>, S>::gradient(*this, x);
  }

  // Default AD hessian.
  template<typename T, Intrepid::Index N>
  Intrepid::Tensor<T, N>
  hessian(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<Beale<S>, S>::hessian(*this, x);
  }

};

//
// Booth's function
//
template<typename S>
class Booth : public Intrepid::Function_Base<Booth<S>, S>
{
public:

  Booth() {}

  static constexpr
  Intrepid::Index
  DIMENSION{2};

  static constexpr
  char const * const
  NAME{"Booth"};

  // Explicit value.
  template<typename T, Intrepid::Index N>
  T
  value(Intrepid::Vector<T, N> const & X)
  {
    Intrepid::Index const
    dimension = X.get_dimension();

    assert(dimension == DIMENSION);

    T const &
    x = X(0);

    T const &
    y = X(1);

    T const
    a = x + 2 * y - 7;

    T const
    b = 2 * x + y - 5;

    T const
    f = a * a + b * b;

    return f;
  }

  // Default AD gradient.
  template<typename T, Intrepid::Index N>
  Intrepid::Vector<T, N>
  gradient(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<Booth<S>, S>::gradient(*this, x);
  }

  // Default AD hessian.
  template<typename T, Intrepid::Index N>
  Intrepid::Tensor<T, N>
  hessian(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<Booth<S>, S>::hessian(*this, x);
  }

};

//
// Goldstein-Price function
//
template<typename S>
class GoldsteinPrice : public Intrepid::Function_Base<GoldsteinPrice<S>, S>
{
public:

  GoldsteinPrice() {}

  static constexpr
  Intrepid::Index
  DIMENSION{2};

  static constexpr
  char const * const
  NAME{"Goldstein-Price"};

  // Explicit value.
  template<typename T, Intrepid::Index N>
  T
  value(Intrepid::Vector<T, N> const & X)
  {
    Intrepid::Index const
    dimension = X.get_dimension();

    assert(dimension == DIMENSION);

    T const &
    x = X(0);

    T const &
    y = X(1);

    T const
    a = x + y + 1;

    T const
    b = 19 - 14 * x + 3 * x * x - 14 * y + 6 * x * y + 3 * y * y;

    T const
    c = 2 * x - 3 * y;

    T const
    d = 18 - 32 * x + 12 * x * x + 48 * y - 36 * x * y + 27 * y * y;

    T const
    e = 1 + a * a * b;

    T const
    f = 30 + c * c * d;

    T const
    fn = e * f;

    return fn;
  }

  // Default AD gradient.
  template<typename T, Intrepid::Index N>
  Intrepid::Vector<T, N>
  gradient(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<GoldsteinPrice<S>, S>::gradient(*this, x);
  }

  // Default AD hessian.
  template<typename T, Intrepid::Index N>
  Intrepid::Tensor<T, N>
  hessian(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<GoldsteinPrice<S>, S>::hessian(*this, x);
  }

};
  
//
// Failure function to test failed mechanism
//
template<typename S>
class Failure : public Intrepid::Function_Base<Failure<S>, S>
{
public:

  Failure() {}

  static constexpr
  Intrepid::Index
  DIMENSION{1};

  static constexpr
  char const * const
  NAME{"Failure"};

  // Explicit value.
  template<typename T, Intrepid::Index N>
  T
  value(Intrepid::Vector<T, N> const & X)
  {
    // Get a convenience reference to the failed flag in case it is used more
    // than once.
    bool &
    failed = Intrepid::Function_Base<Failure<S>, S>::failed;

    // Set the flag to signal that an unrecoverable error happened.

    failed = true;

    T const
    fn = 0.0;

    return fn;
  }

  // Default AD gradient.
  template<typename T, Intrepid::Index N>
  Intrepid::Vector<T, N>
  gradient(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<Failure<S>, S>::gradient(*this, x);
  }

  // Default AD hessian.
  template<typename T, Intrepid::Index N>
  Intrepid::Tensor<T, N>
  hessian(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<Failure<S>, S>::hessian(*this, x);
  }

};

//
// Non-monotonic function to test monotonicity enforcement.
//
template<typename S>
class Mesa : public Intrepid::Function_Base<Mesa<S>, S>
{
public:

  Mesa() {}

  static constexpr
  Intrepid::Index
  DIMENSION{1};

  static constexpr
  char const * const
  NAME{"Mesa"};

  // Explicit value.
  template<typename T, Intrepid::Index N>
  T
  value(Intrepid::Vector<T, N> const & X)
  {
    T const &
    x = X(0);

    T
    y = x * x;

    if (-1.0 <= x && x <= 1.0) {
      y = y + 100.0;
    }

    return y;
  }

  // Default AD gradient.
  template<typename T, Intrepid::Index N>
  Intrepid::Vector<T, N>
  gradient(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<Mesa<S>, S>::gradient(*this, x);
  }

  // Default AD hessian.
  template<typename T, Intrepid::Index N>
  Intrepid::Tensor<T, N>
  hessian(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<Mesa<S>, S>::hessian(*this, x);
  }

};

//
// Function to test boundedness or residual enforcement.
//
template<typename S>
class Sigmoid : public Intrepid::Function_Base<Sigmoid<S>, S>
{
public:

  Sigmoid() {}

  static constexpr
  Intrepid::Index
  DIMENSION{1};

  static constexpr
  char const * const
  NAME{"Sigmoid"};

  // Explicit value.
  template<typename T, Intrepid::Index N>
  T
  value(Intrepid::Vector<T, N> const & X)
  {
    T const &
    x = X(0);

    T const
    x2 = x * x;

    T const
    x4 = x2 * x2;

    T const
    x8 = x4 * x4;

    T const
    x16 = x8 * x8;

    T const
    x32 = x16 * x16;

    T
    y = x * x32;

    return y;
  }

  // Default AD gradient.
  template<typename T, Intrepid::Index N>
  Intrepid::Vector<T, N>
  gradient(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<Sigmoid<S>, S>::gradient(*this, x);
  }

  // Default AD hessian.
  template<typename T, Intrepid::Index N>
  Intrepid::Tensor<T, N>
  hessian(Intrepid::Vector<T, N> const & x)
  {
    return Intrepid::Function_Base<Sigmoid<S>, S>::hessian(*this, x);
  }

};

//
// Test the solution methods by themselves.
//

// Test one function with one method.
template <typename STEP, typename FN, typename T, Intrepid::Index N>
bool
solveFNwithSTEP(STEP & step_method, FN & function, Intrepid::Vector<T, N> & x)
{
  Intrepid::Minimizer<T, N>
  minimizer;

  minimizer.solve(step_method, function, x);

  minimizer.printReport(std::cout);

  return minimizer.converged;
}

// Test one system with various methods.
template <typename FN, typename T, Intrepid::Index N>
bool
solveFN(FN & function, Intrepid::Vector<T, N> const & x)
{
  bool
  all_ok = true;

  Intrepid::Vector<T, N>
  y;

  Intrepid::NewtonStep<FN, T, N>
  newton_step;

  y = x;

  bool const
  newton_ok = solveFNwithSTEP(newton_step, function, y);

  all_ok = all_ok && newton_ok;

  Intrepid::TrustRegionStep<FN, T, N>
  trust_region_step;

  y = x;

  bool const
  trust_region_ok = solveFNwithSTEP(trust_region_step, function, y);

  all_ok = all_ok && trust_region_ok;

  Intrepid::ConjugateGradientStep<FN, T, N>
  pcg_step;

  y = x;

  bool const
  pcg_ok = solveFNwithSTEP(pcg_step, function, y);

  all_ok = all_ok && pcg_ok;

  Intrepid::LineSearchRegularizedStep<FN, T, N>
  line_search_step;

  y = x;

  bool const
  line_search_ok = solveFNwithSTEP(line_search_step, function, y);

  all_ok = all_ok && line_search_ok;

  return all_ok;
}

// Test various systems with various methods.
bool testSystemsAndMethods()
{
  constexpr Intrepid::Index
  max_dimension{2};

  bool
  all_ok = true;

  Intrepid::Vector<Real, max_dimension>
  x;

  SquareRootNLS<Real>
  square_root(2.0);

  x.set_dimension(SquareRootNLS<Real>::DIMENSION);

  x(0) = 10.0;

  bool const
  square_root_ok = solveFN(square_root, x);

  all_ok = all_ok && square_root_ok;

  QuadraticNLS<Real>
  quadratic(10.0, 15.0, 1.0);

  x.set_dimension(QuadraticNLS<Real>::DIMENSION);

  x(0) = -15.0;
  x(1) = -10.0;

  bool const
  quadratic_ok = solveFN(quadratic, x);

  all_ok = all_ok && quadratic_ok;

  GaussianNLS<Real>
  gaussian(1.0, 2.0, 0.125);

  x.set_dimension(GaussianNLS<Real>::DIMENSION);

  x(0) = 0.0;
  x(1) = 0.0;

  bool const
  gaussian_ok = solveFN(gaussian, x);

  all_ok = all_ok && gaussian_ok;

  BananaNLS<Real>
  banana;

  x.set_dimension(BananaNLS<Real>::DIMENSION);

  x(0) = 0.0;
  x(1) = 3.0;

  bool const
  banana_ok = solveFN(banana, x);

  all_ok = all_ok && banana_ok;

  MatyasNLS<Real>
  matyas;

  x.set_dimension(MatyasNLS<Real>::DIMENSION);

  x(0) = 10.0;
  x(1) =  0.0;

  bool const
  matyas_ok = solveFN(matyas, x);

  all_ok = all_ok && matyas_ok;

  McCormickNLS<Real>
  mccormick;

  x.set_dimension(McCormickNLS<Real>::DIMENSION);

  x(0) = -0.5;
  x(1) = -1.5;

  bool const
  mccormick_ok = solveFN(mccormick, x);

  all_ok = all_ok && mccormick_ok;

  StyblinskiTangNLS<Real>
  styblinski_tang;

  x.set_dimension(StyblinskiTangNLS<Real>::DIMENSION);

  x(0) = -4.0;
  x(1) = -4.0;

  bool const
  styblinski_tang_ok = solveFN(styblinski_tang, x);

  all_ok = all_ok && styblinski_tang_ok;

  Paraboloid<Real>
  paraboloid;

  x.set_dimension(Paraboloid<Real>::DIMENSION);

  x(0) = 128.0;
  x(1) = 256.0;;

  bool const
  paraboloid_ok = solveFN(paraboloid, x);

  all_ok = all_ok && paraboloid_ok;

  Beale<Real>
  beale;

  x.set_dimension(Beale<Real>::DIMENSION);

  x(0) = -4.5;
  x(1) = -4.5;

  bool const
  beale_ok = solveFN(beale, x);

  all_ok = all_ok && beale_ok;

  Booth<Real>
  booth;

  x.set_dimension(Booth<Real>::DIMENSION);

  x(0) = -10.0;
  x(1) = -10.0;

  bool const
  booth_ok = solveFN(booth, x);

  all_ok = all_ok && booth_ok;

  GoldsteinPrice<Real>
  goldstein_price;

  x.set_dimension(GoldsteinPrice<Real>::DIMENSION);

  x(0) = 2.0;
  x(1) = 2.0;

  bool const
  goldstein_price_ok = solveFN(goldstein_price, x);

  all_ok = all_ok && goldstein_price_ok;

  return all_ok;
}

} // anonymous namespace

TEST(NonlinearSystems, NonlinearMethods)
{
  bool const
  passed = testSystemsAndMethods();

  ASSERT_EQ(passed, true);
}

TEST(Testing, OptimizationMethods)
{
  constexpr Intrepid::Index
  dimension{2};

  using MIN = Intrepid::Minimizer<Real, dimension>;
  using FN = BananaNLS<Real>;
  using STEP = Intrepid::NewtonStep<FN, Real, dimension>;

  MIN
  minimizer;

  FN
  banana;

  STEP
  step;

  Intrepid::Vector<Real, dimension>
  x;

  x(0) = 0.0;
  x(1) = 3.0;

  minimizer.solve(step, banana, x);

  minimizer.printReport(std::cout);

  ASSERT_EQ(minimizer.converged, true);
}

TEST(Testing, ValueGradientHessian)
{
  constexpr Intrepid::Index
  dimension{2};

  Paraboloid<Real>
  p;

  Intrepid::Vector<Real, dimension> const
  x(0.0, 0.0);

  Real const
  f = p.value(x);

  Intrepid::Vector<Real, dimension> const
  df = p.gradient(x);

  Intrepid::Tensor<Real, dimension> const
  ddf = p.hessian(x);

  std::cout << "Point   : " << x << '\n';
  std::cout << "Value   : " << f << '\n';
  std::cout << "Gradient: " << df << '\n';
  std::cout << "Hessian : " << ddf << '\n';

  Intrepid::Tensor<Real, dimension> const
  I = Intrepid::identity<Real, dimension>(dimension);

  Real const
  error = std::sqrt(f) + Intrepid::norm(df) + Intrepid::norm(ddf - 2.0 * I);

  ASSERT_LE(error, Intrepid::machine_epsilon<Real>());
}

TEST(Testing, MixedStorage)
{
  Intrepid::Index const
  dimension{2};

  std::cout << '\n';

  Intrepid::Vector<Real, 3>
  v(1.0, 2.0, 3.0);

  v.set_dimension(dimension);

  std::cout << "Vector   : " << v << '\n';

  Intrepid::Tensor<Real, 3>
  A(1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 7.0, 8.0, 9.0);

  A.set_dimension(dimension);

  std::cout << "Tensor   : " << A << '\n';

  Intrepid::Matrix<Real, 3, 4>
  B(Intrepid::ONES);

  B.set_dimensions(4, 2);

  std::cout << "Matrix   : " << B << '\n';

  bool const
  passed = v.get_dimension() == dimension && A.get_dimension() == dimension &&
    B.get_num_rows() == 4 && B.get_num_cols() == 2;

  ASSERT_EQ(passed, true);
}

TEST(Testing, FailedFlag)
{
  constexpr Intrepid::Index
  dimension{1};

  using MIN = Intrepid::Minimizer<Real, dimension>;
  using FN = Failure<Real>;
  using STEP = Intrepid::NewtonStep<FN, Real, dimension>;

  MIN
  minimizer;

  FN
  fn;

  STEP
  step;

  Intrepid::Vector<Real, dimension>
  x;

  x(0) = 0.0;

  minimizer.solve(step, fn, x);

  ASSERT_EQ(minimizer.failed, true);
}

TEST(Testing, Monotonicity)
{
  constexpr Intrepid::Index
  dimension{1};

  using MIN = Intrepid::Minimizer<Real, dimension>;
  using FN = Mesa<Real>;
  using STEP = Intrepid::NewtonStep<FN, Real, dimension>;

  MIN
  minimizer;

  minimizer.enforce_monotonicity = true;

  FN
  fn;

  STEP
  step;

  Intrepid::Vector<Real, dimension>
  x;

  x(0) = 2.0;

  minimizer.solve(step, fn, x);

  ASSERT_EQ(minimizer.monotonic, false);
  ASSERT_EQ(minimizer.failed, true);
}

TEST(Testing, Boundedness)
{
  constexpr Intrepid::Index
  dimension{1};

  using MIN = Intrepid::Minimizer<Real, dimension>;
  using FN = Sigmoid<Real>;
  using STEP = Intrepid::NewtonStep<FN, Real, dimension>;

  MIN
  minimizer;

  minimizer.enforce_boundedness = true;

  FN
  fn;

  STEP
  step;

  Intrepid::Vector<Real, dimension>
  x;

  x(0) = 0.5;

  minimizer.solve(step, fn, x);

  minimizer.printReport(std::cout);

  ASSERT_EQ(minimizer.bounded, true);
  ASSERT_EQ(minimizer.failed, false);
}

} // namespace Intrepid
