// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_TestFunctions_h)
#define MiniTensor_TestFunctions_h

#include "MiniTensor_Solvers.h"

namespace minitensor {

/// \addtogroup minitensor_solvers
/// @{

//
// Define some nonlinear systems (NLS) to test nonlinear solution methods.
//

//
//
//
///
/// Square root NLS: the residual is \f$ r(x) = x^2 - c \f$,
/// whose root is \f$ x = \sqrt{c} \f$. One-dimensional test
/// for nonlinear system solvers.
///
template<typename S, Index M = 1>
class SquareRoot : public Function_Base<SquareRoot<S, M>, S, M>
{
public:

  ///
  /// Constructor.
  /// \param c Value whose square root is sought.
  ///
  SquareRoot(S const c) : c_(c)
  {
  }

  ///
  /// Function name.
  ///
  static constexpr
  char const * const
  NAME{"Square Root"};

  ///
  /// Base class type.
  ///
  using Base = Function_Base<SquareRoot<S, M>, S, M>;

  // Default value.
  ///
  /// Objective value.
  ///
  template<typename T, Index N>
  T
  value(Vector<T, N> const & x)
  {
    return Base::value(*this, x);
  }

  // Explicit gradient.
  ///
  /// Gradient of the objective.
  ///
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x) const
  {
    Index const
    dimension = x.get_dimension();

    assert(dimension == Base::DIMENSION);

    Vector<T, N>
    r(dimension);

    r(0) = x(0) * x(0) - c_;

    return r;
  }

  // Default AD hessian.
  ///
  /// Hessian of the objective.
  ///
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

private:
  ///
  /// Value whose square root is sought.
  ///
  S const
  c_{0.0};
};

//
//
//
///
/// Quadratic NLS: the residual is the gradient of
/// \f$ c((x_1-a)^2 + (x_2-b)^2) \f$. Two-dimensional, convex,
/// single minimum at \f$ (a, b) \f$.
///
template<typename S, Index M = 2>
class Quadratic : public Function_Base<Quadratic<S, M>, S, M>
{
public:
  ///
  /// Constructor.
  /// \param a First coordinate of the minimizer.
  /// \param b Second coordinate of the minimizer.
  /// \param c Curvature scale factor.
  ///
  Quadratic(S const a, S const b, S const c) :  a_(a), b_(b), c_(c)
  {
  }

  ///
  /// Function name.
  ///
  static constexpr
  char const * const
  NAME{"Quadratic"};

  ///
  /// Base class type.
  ///
  using Base = Function_Base<Quadratic<S, M>, S, M>;

  // Default value.
  ///
  /// Objective value.
  ///
  template<typename T, Index N>
  T
  value(Vector<T, N> const & x)
  {
    return Base::value(*this, x);
  }

  // Explicit gradient.
  ///
  /// Gradient of the objective.
  ///
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x) const
  {
    Index const
    dimension = x.get_dimension();

    assert(dimension == Base::DIMENSION);

    Vector<T, N>
    r(dimension);

    r(0) = 2.0 * c_ * (x(0) - a_);
    r(1) = 2.0 * c_ * (x(1) - b_);

    return r;
  }

  // Default AD hessian.
  ///
  /// Hessian of the objective.
  ///
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

private:
  ///
  /// First coordinate of the minimizer.
  ///
  S const
  a_{0.0};

  ///
  /// Second coordinate of the minimizer.
  ///
  S const
  b_{0.0};

  ///
  /// Curvature scale factor.
  ///
  S const
  c_{0.0};
};

//
//
//
///
/// Inverted Gaussian NLS: the residual is the gradient of
/// \f$ -e^{-c^2((x_1-a)^2+(x_2-b)^2)} \f$ scaled by
/// \f$ c \f$. Two-dimensional; minimum at \f$ (a, b) \f$,
/// nearly flat away from it, stressing globalization.
///
template<typename S, Index M = 2>
class Gaussian : public Function_Base<Gaussian<S, M>, S, M>
{
public:
  ///
  /// Constructor.
  /// \param a First coordinate of the center (minimizer).
  /// \param b Second coordinate of the center (minimizer).
  /// \param c Inverse width of the Gaussian.
  ///
  Gaussian(S const a, S const b, S const c) : a_(a), b_(b), c_(c)
  {
  }

  ///
  /// Function name.
  ///
  static constexpr
  char const * const
  NAME{"Inverted Gaussian"};

  ///
  /// Base class type.
  ///
  using Base = Function_Base<Gaussian<S, M>, S, M>;

  // Default value.
  ///
  /// Objective value.
  ///
  template<typename T, Index N>
  T
  value(Vector<T, N> const & x)
  {
    return Base::value(*this, x);
  }

  // Explicit gradient.
  ///
  /// Gradient of the objective.
  ///
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x) const
  {
    Index const
    dimension = x.get_dimension();

    assert(dimension == Base::DIMENSION);

    Vector<T, N>
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
  ///
  /// Hessian of the objective.
  ///
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

private:
  ///
  /// First coordinate of the center.
  ///
  S const
  a_{0.0};

  ///
  /// Second coordinate of the center.
  ///
  S const
  b_{0.0};

  ///
  /// Inverse width of the Gaussian.
  ///
  S const
  c_{0.0};
};

//
//
//
///
/// Rosenbrock's banana function as an NLS: the residual is
/// the gradient of \f$ (1-x_1)^2 + 100(x_2-x_1^2)^2 \f$.
/// Two-dimensional; its curved, flat-bottomed valley stresses
/// step control. Minimum at \f$ (1, 1) \f$.
///
template<typename S, Index M = 2>
class Banana : public Function_Base<Banana<S, M>, S, M>
{
public:

  ///
  /// Constructor.
  ///
  Banana()
  {
  }

  ///
  /// Function name.
  ///
  static constexpr
  char const * const
  NAME{"Rosenbrock's Banana"};

  ///
  /// Base class type.
  ///
  using Base = Function_Base<Banana<S, M>, S, M>;

  // Default value.
  ///
  /// Objective value.
  ///
  template<typename T, Index N>
  T
  value(Vector<T, N> const & x)
  {
    return Base::value(*this, x);
  }

  // Explicit gradient.
  ///
  /// Gradient of the objective.
  ///
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x) const
  {
    Index const
    dimension = x.get_dimension();

    assert(dimension == Base::DIMENSION);

    Vector<T, N>
    r(dimension);

    r(0) = 2.0 * (x(0) - 1.0) + 400.0 * x(0) * (x(0) * x(0) - x(1));
    r(1) = 200.0 * (x(1) - x(0) * x(0));

    return r;
  }

  // Default AD hessian.
  ///
  /// Hessian of the objective.
  ///
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
//
//
///
/// Matyas function NLS: the residual is the gradient of
/// \f$ 0.26(x_1^2+x_2^2) - 0.48 x_1 x_2 \f$.
/// Two-dimensional, convex but ill-conditioned; minimum at
/// the origin.
///
template<typename S, Index M = 2>
class Matyas : public Function_Base<Matyas<S, M>, S, M>
{
public:

  ///
  /// Constructor.
  ///
  Matyas() {}

  ///
  /// Function name.
  ///
  static constexpr
  char const * const
  NAME{"Matyas"};

  ///
  /// Base class type.
  ///
  using Base = Function_Base<Matyas<S, M>, S, M>;

  // Default value.
  ///
  /// Objective value.
  ///
  template<typename T, Index N>
  T
  value(Vector<T, N> const & x)
  {
    return Base::value(*this, x);
  }

  // Explicit gradient.
  ///
  /// Gradient of the objective.
  ///
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x) const
  {
    Index const
    dimension = x.get_dimension();

    assert(dimension == Base::DIMENSION);

    Vector<T, N>
    r(dimension);

    r(0) = (13.0 * x(0) - 12.0 * x(1)) / 25.0;
    r(1) = (13.0 * x(1) - 12.0 * x(0)) / 25.0;

    return r;
  }

  // Default AD hessian.
  ///
  /// Hessian of the objective.
  ///
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
//
//
///
/// McCormick function NLS: the residual is the gradient of
/// \f$ \sin(x_1+x_2) + (x_1-x_2)^2 - 1.5 x_1
/// + 2.5 x_2 + 1 \f$. Two-dimensional with multiple minima;
/// global minimum near \f$ (-0.547, -1.547) \f$.
///
template<typename S, Index M = 2>
class McCormick : public Function_Base<McCormick<S, M>, S, M>
{
public:

  ///
  /// Constructor.
  ///
  McCormick() {}

  ///
  /// Function name.
  ///
  static constexpr
  char const * const
  NAME{"McCormick"};

  ///
  /// Base class type.
  ///
  using Base = Function_Base<McCormick<S, M>, S, M>;

  // Default value.
  ///
  /// Objective value.
  ///
  template<typename T, Index N>
  T
  value(Vector<T, N> const & x)
  {
    return Base::value(*this, x);
  }

  // Explicit gradient.
  ///
  /// Gradient of the objective.
  ///
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x) const
  {
    Index const
    dimension = x.get_dimension();

    assert(dimension == Base::DIMENSION);

    Vector<T, N>
    r(dimension);

    r(0) = std::cos(x(0) + x(1)) + 2.0 * x(0) - 2.0 * x(1) - 1.5;
    r(1) = std::cos(x(0) + x(1)) - 2.0 * x(0) + 2.0 * x(1) + 2.5;

    return r;
  }

  // Default AD hessian.
  ///
  /// Hessian of the objective.
  ///
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
//
//
///
/// Styblinski-Tang function NLS: the residual is the
/// gradient of
/// \f$ \frac{1}{2}\sum_{i=1}^{2}
/// (x_i^4 - 16 x_i^2 + 5 x_i) \f$.
/// Two-dimensional with multiple local minima; global
/// minimum at \f$ x_i \approx -2.9035 \f$.
///
template<typename S, Index M = 2>
class StyblinskiTang : public Function_Base<StyblinskiTang<S, M>, S, M>
{
public:

  ///
  /// Constructor.
  ///
  StyblinskiTang() {}

  ///
  /// Function name.
  ///
  static constexpr
  char const * const
  NAME{"Styblinski-Tang"};

  ///
  /// Base class type.
  ///
  using Base = Function_Base<StyblinskiTang<S, M>, S, M>;

  // Default value.
  ///
  /// Objective value.
  ///
  template<typename T, Index N>
  T
  value(Vector<T, N> const & x)
  {
    return Base::value(*this, x);
  }

  // Explicit gradient.
  ///
  /// Gradient of the objective.
  ///
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x) const
  {
    Index const
    dimension = x.get_dimension();

    assert(dimension == Base::DIMENSION);

    Vector<T, N>
    r(dimension);

    r(0) = 2.0 * x(0) * x(0) * x(0) - 16.0 * x(0) + 2.5;
    r(1) = 2.0 * x(1) * x(1) * x(1) - 16.0 * x(1) + 2.5;

    return r;
  }

  // Default AD hessian.
  ///
  /// Hessian of the objective.
  ///
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
// Define some nonlinear functions (NLF) to test nonlinear optimization methods.
//

//
// Paraboloid of revolution
//
///
/// Paraboloid of revolution
/// \f$ (x_1-x_c)^2 + (x_2-y_c)^2 \f$. Two-dimensional,
/// convex, single minimum at \f$ (x_c, y_c) \f$.
///
template<typename S, Index M = 2>
class Paraboloid : public Function_Base<Paraboloid<S, M>, S, M>
{
public:

  ///
  /// Constructor.
  /// \param xc First coordinate of the minimizer.
  /// \param yc Second coordinate of the minimizer.
  ///
  Paraboloid(S xc = 0.0, S yc = 0.0) : xc_(xc), yc_(yc)
  {
  }

  ///
  /// Function name.
  ///
  static constexpr
  char const * const
  NAME{"Paraboloid"};

  ///
  /// Base class type.
  ///
  using Base = Function_Base<Paraboloid<S, M>, S, M>;

  // Explicit value.
  ///
  /// Objective value.
  ///
  template<typename T, Index N>
  T
  value(Vector<T, N> const & x)
  {
    assert(x.get_dimension() == Base::DIMENSION);

    T const
    a = x(0) - xc_;

    T const
    b = x(1) - yc_;

    T const
    f = a * a + b * b;

    return f;
  }

  // Default AD gradient.
  ///
  /// Gradient of the objective.
  ///
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

  // Default AD hessian.
  ///
  /// Hessian of the objective.
  ///
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

private:
  ///
  /// First coordinate of the minimizer.
  ///
  S
  xc_{0.0};

  ///
  /// Second coordinate of the minimizer.
  ///
  S
  yc_{0.0};
};

//
//
//
///
/// Rosenbrock's function \f$ (a-x_1)^2 + b(x_2-x_1^2)^2 \f$
/// with defaults \f$ a=1, b=100 \f$. Two-dimensional; its
/// curved banana valley stresses line search and step
/// control. Minimum at \f$ (a, a^2) \f$.
///
template<typename S, Index M = 2>
class Rosenbrock : public Function_Base<Rosenbrock<S, M>, S, M>
{
public:

  ///
  /// Constructor.
  /// \param a Location parameter; the minimizer is
  /// \f$ (a, a^2) \f$.
  /// \param b Weight of the valley (coupling) term.
  ///
  Rosenbrock(S a = 1.0, S b = 100.0) : a_(a), b_(b)
  {
  }

  ///
  /// Function name.
  ///
  static constexpr
  char const * const
  NAME{"Rosenbrock's Function 2D"};

  ///
  /// Base class type.
  ///
  using Base = Function_Base<Rosenbrock<S, M>, S, M>;

  // Explicit value.
  ///
  /// Objective value.
  ///
  template<typename T, Index N>
  T
  value(Vector<T, N> const & x)
  {
    T const
    a = (a_ - x(0));

    T const
    b = (x(1) - x(0) * x(0));

    return a * a + b_ * b * b;
  }

  // Default AD gradient.
  ///
  /// Gradient of the objective.
  ///
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

  // Default AD hessian.
  ///
  /// Hessian of the objective.
  ///
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

private:
  ///
  /// Location parameter.
  ///
  S
  a_{1.0};

  ///
  /// Weight of the valley (coupling) term.
  ///
  S
  b_{100.0};
};

//
// Beale's function
//
///
/// Beale's function \f$ (1.5-x+xy)^2 + (2.25-x+xy^2)^2
/// + (2.625-x+xy^3)^2 \f$. Two-dimensional with sharp
/// curved valleys; minimum at \f$ (3, 0.5) \f$.
///
template<typename S, Index M = 2>
class Beale : public Function_Base<Beale<S, M>, S, M>
{
public:

  ///
  /// Constructor.
  ///
  Beale() {}

  ///
  /// Function name.
  ///
  static constexpr
  char const * const
  NAME{"Beale"};

  ///
  /// Base class type.
  ///
  using Base = Function_Base<Beale<S, M>, S, M>;

  // Explicit value.
  ///
  /// Objective value.
  ///
  template<typename T, Index N>
  T
  value(Vector<T, N> const & X)
  {
    assert(X.get_dimension() == Base::DIMENSION);

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
  ///
  /// Gradient of the objective.
  ///
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

  // Default AD hessian.
  ///
  /// Hessian of the objective.
  ///
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
// Booth's function
//
///
/// Booth's function \f$ (x+2y-7)^2 + (2x+y-5)^2 \f$.
/// Two-dimensional, convex quadratic; minimum at
/// \f$ (1, 3) \f$.
///
template<typename S, Index M = 2>
class Booth : public Function_Base<Booth<S, M>, S, M>
{
public:

  ///
  /// Constructor.
  ///
  Booth() {}

  ///
  /// Function name.
  ///
  static constexpr
  char const * const
  NAME{"Booth"};

  ///
  /// Base class type.
  ///
  using Base = Function_Base<Booth<S, M>, S, M>;

  // Explicit value.
  ///
  /// Objective value.
  ///
  template<typename T, Index N>
  T
  value(Vector<T, N> const & X)
  {
    assert(X.get_dimension() == Base::DIMENSION);

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
  ///
  /// Gradient of the objective.
  ///
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

  // Default AD hessian.
  ///
  /// Hessian of the objective.
  ///
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
// Goldstein-Price function
//
///
/// Goldstein-Price function
/// \f$ [1+(x+y+1)^2(19-14x+3x^2-14y+6xy+3y^2)] \f$
/// \f$ [30+(2x-3y)^2(18-32x+12x^2+48y-36xy+27y^2)] \f$.
/// Two-dimensional with several local minima; global
/// minimum of 3 at \f$ (0, -1) \f$.
///
template<typename S, Index M = 2>
class GoldsteinPrice : public Function_Base<GoldsteinPrice<S, M>, S, M>
{
public:

  ///
  /// Constructor.
  ///
  GoldsteinPrice() {}

  ///
  /// Function name.
  ///
  static constexpr
  char const * const
  NAME{"Goldstein-Price"};

  ///
  /// Base class type.
  ///
  using Base = Function_Base<GoldsteinPrice<S, M>, S, M>;

  // Explicit value.
  ///
  /// Objective value.
  ///
  template<typename T, Index N>
  T
  value(Vector<T, N> const & X)
  {
    assert(X.get_dimension() == Base::DIMENSION);

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
  ///
  /// Gradient of the objective.
  ///
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

  // Default AD hessian.
  ///
  /// Hessian of the objective.
  ///
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
// Failure function to test failed mechanism
//
///
/// Function that always reports failure: value() sets the
/// failed flag and returns zero. One-dimensional; exercises
/// the solvers' failure detection mechanism.
///
template<typename S, Index M = 1>
class Failure : public Function_Base<Failure<S, M>, S, M>
{
public:

  ///
  /// Constructor.
  ///
  Failure() {}

  ///
  /// Function name.
  ///
  static constexpr
  char const * const
  NAME{"Failure"};

  ///
  /// Base class type.
  ///
  using Base = Function_Base<Failure<S, M>, S, M>;

  // Explicit value.
  ///
  /// Objective value; flags an unrecoverable failure and
  /// returns zero.
  ///
  template<typename T, Index N>
  T
  value(Vector<T, N> const & X)
  {
    // Set the flag to signal that an unrecoverable error happened.
    this->set_failed("Testing failure mechanism");

    T const
    fn = 0.0;

    return fn;
  }

  // Default AD gradient.
  ///
  /// Gradient of the objective.
  ///
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

  // Default AD hessian.
  ///
  /// Hessian of the objective.
  ///
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
// Non-monotonic function to test monotonicity enforcement.
//
///
/// Mesa function: \f$ x^2 \f$ plus a plateau of height 100
/// on \f$ [-1, 1] \f$. One-dimensional, discontinuous and
/// non-monotonic near the origin; tests monotonicity
/// enforcement.
///
template<typename S, Index M = 1>
class Mesa : public Function_Base<Mesa<S, M>, S, M>
{
public:

  ///
  /// Constructor.
  ///
  Mesa() {}

  ///
  /// Function name.
  ///
  static constexpr
  char const * const
  NAME{"Mesa"};

  ///
  /// Base class type.
  ///
  using Base = Function_Base<Mesa<S, M>, S, M>;

  // Explicit value.
  ///
  /// Objective value.
  ///
  template<typename T, Index N>
  T
  value(Vector<T, N> const & X)
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
  ///
  /// Gradient of the objective.
  ///
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

  // Default AD hessian.
  ///
  /// Hessian of the objective.
  ///
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
// Function to test boundedness or residual enforcement.
//
///
/// Steep sigmoid-like monomial \f$ x^{33} \f$.
/// One-dimensional, extremely flat near the origin and steep
/// beyond \f$ |x| = 1 \f$; tests boundedness and residual
/// enforcement.
///
template<typename S, Index M = 1>
class Sigmoid : public Function_Base<Sigmoid<S, M>, S, M>
{
public:

  ///
  /// Constructor.
  ///
  Sigmoid() {}

  ///
  /// Function name.
  ///
  static constexpr
  char const * const
  NAME{"Sigmoid"};

  ///
  /// Base class type.
  ///
  using Base = Function_Base<Sigmoid<S, M>, S, M>;

  // Explicit value.
  ///
  /// Objective value.
  ///
  template<typename T, Index N>
  T
  value(Vector<T, N> const & X)
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
  ///
  /// Gradient of the objective.
  ///
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

  // Default AD hessian.
  ///
  /// Hessian of the objective.
  ///
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
// Functions to test constraint interface.
//

//
// Identity
//
///
/// Identity map equality constraint \f$ c(x) = x \f$ with
/// NC constraints in NV variables.
///
template<typename S, Index NC, Index NV>
class Identity : public Equality_Constraint<Identity<S, NC, NV>, S, NC, NV>
{
public:

  ///
  /// Constructor.
  ///
  Identity() {}

  ///
  /// Function name.
  ///
  static constexpr
  char const * const
  NAME{"Identity Map"};

  ///
  /// Base class type.
  ///
  using Base = Equality_Constraint<Identity<S, NC, NV>, S, NC, NV>;

  // Explicit value.
  ///
  /// Constraint value.
  ///
  template<typename T, Index N>
  Vector<T, NC>
  value(Vector<T, N> const & x)
  {
    assert(x.get_dimension() == NV);
    return x;
  }

  // Default AD gradient.
  ///
  /// Gradient (Jacobian) of the constraint.
  ///
  template<typename T, Index N>
  Matrix<T, NC, NV>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }
};

//
// A nonlinear function
//
///
/// Nonlinear equality constraints, 3 equations in 5
/// variables: \f$ x \cdot x - 10 \f$,
/// \f$ x_2 x_3 - 5 x_4 x_5 \f$,
/// \f$ x_1^3 + x_2^3 + 1 \f$.
///
template<typename S, Index NC = 3, Index NV = 5>
class Nonlinear01 : public Equality_Constraint<Nonlinear01<S, NC, NV>, S, NC, NV>
{
public:

  ///
  /// Constructor.
  ///
  Nonlinear01() {}

  ///
  /// Function name.
  ///
  static constexpr
  char const * const
  NAME{"Nonlinear 01"};

  ///
  /// Base class type.
  ///
  using Base = Equality_Constraint<Nonlinear01<S, NC, NV>, S, NC, NV>;

  // Explicit value.
  ///
  /// Constraint value.
  ///
  template<typename T, Index N = 5>
  Vector<T, NC>
  value(Vector<T, N> const & x)
  {
    assert(x.get_dimension() == NV);

    Vector<T, NC>
    c(Filler::ZEROS);

    c(0) = dot(x, x) - 10.0;

    c(1) = x(1) * x(2) - 5.0 * x(3) * x(4);

    c(2) = x(0) * x(0) * x(0) + x(1) * x(1) * x(1) + 1.0;

    return c;
  }

  // Default AD gradient.
  ///
  /// Gradient (Jacobian) of the constraint.
  ///
  template<typename T, Index N = 5>
  Matrix<T, NC, NV>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }
};

//
// Circumference feasible region
//
///
/// Circumference equality constraint
/// \f$ r^2 - \|x - c\|^2 = 0 \f$: the feasible set is
/// the circle of radius \f$ r \f$ centered at
/// \f$ c = (x_c, y_c) \f$.
///
template<typename S, Index NC = 1, Index NV = 2>
class Circumference : public Equality_Constraint<Circumference<S, NC, NV>, S, NC, NV>
{
public:

  ///
  /// Constructor.
  /// \param r Radius of the circumference.
  /// \param xc First coordinate of the center.
  /// \param yc Second coordinate of the center.
  ///
  Circumference(S const r, S const xc = S(0.0), S const yc = S(0.0)) : r_(r)
  {
    c_(0) = xc;
    c_(1) = yc;
  }

  ///
  /// Function name.
  ///
  static constexpr
  char const * const
  NAME{"Circumference"};

  ///
  /// Base class type.
  ///
  using Base = Equality_Constraint<Circumference<S, NC, NV>, S, NC, NV>;

  // Explicit value.
  ///
  /// Constraint value.
  ///
  template<typename T, Index N = 2>
  Vector<T, NC>
  value(Vector<T, N> const & x)
  {
    assert(x.get_dimension() == NV);

    Vector<T, NC>
    f(Filler::ZEROS);

    f(0) = r_ * r_ - norm_square(x - c_);

    return f;
  }

  // Default AD gradient.
  ///
  /// Gradient (Jacobian) of the constraint.
  ///
  template<typename T, Index N = 2>
  Matrix<T, NC, NV>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

private:
  ///
  /// Radius of the circumference.
  ///
  S
  r_{0.0};

  ///
  /// Center of the circumference.
  ///
  Vector<S, NV>
  c_;
};

//
// Circle feasible region
//
///
/// Circle (disk) inequality constraint
/// \f$ r^2 - \|x - c\|^2 \geq 0 \f$: the feasible set
/// is the closed disk of radius \f$ r \f$ centered at
/// \f$ c = (x_c, y_c) \f$.
///
template<typename S, Index NC = 1, Index NV = 2>
class Circle : public Inequality_Constraint<Circle<S, NC, NV>, S, NC, NV>
{
public:

  ///
  /// Constructor.
  /// \param r Radius of the disk.
  /// \param xc First coordinate of the center.
  /// \param yc Second coordinate of the center.
  ///
  Circle(S const r, S const xc = S(0.0), S const yc = S(0.0)) : r_(r)
  {
    c_(0) = xc;
    c_(1) = yc;
  }

  ///
  /// Function name.
  ///
  static constexpr
  char const * const
  NAME{"Circle constraint"};

  ///
  /// Base class type.
  ///
  using Base = Inequality_Constraint<Circle<S, NC, NV>, S, NC, NV>;

  // Explicit value.
  ///
  /// Constraint value.
  ///
  template<typename T, Index N = 2>
  Vector<T, NC>
  value(Vector<T, N> const & x)
  {
    assert(x.get_dimension() == NV);

    Vector<T, NC>
    f(Filler::ZEROS);

    f(0) = r_ * r_ - norm_square(x - c_);

    return f;
  }

  // Default AD gradient.
  ///
  /// Gradient (Jacobian) of the constraint.
  ///
  template<typename T, Index N = 2>
  Matrix<T, NC, NV>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

private:
  ///
  /// Radius of the disk.
  ///
  S
  r_{0.0};

  ///
  /// Center of the disk.
  ///
  Vector<S, NV>
  c_;
};

/// @}
} // namespace minitensor

#endif // MiniTensor_TestFunctions_h
