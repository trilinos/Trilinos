// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Quaternion_h)
#define MiniTensor_Quaternion_h

// Quaternions and conversions among rotation representations:
// rotation matrix (rt), unit quaternion (q) and principal rotation
// vector (rv). Ported from Norma.jl, which in turn adopted them from
// the ONDAP FEM code. See: Object-oriented finite-element dynamic
// simulation of geometrically nonlinear space structures, Victor
// Balopoulos, Ph.D. dissertation, Cornell University, 1997.
#include "MiniTensor_Rotations.h"

namespace minitensor {

namespace impl {

//
// Compute sin(x)/x with asymptotic expansions near zero to avoid
// division by small values. This function is frequently encountered
// in rotational algebra when evaluating exp/log maps of
// skew-symmetric matrices.
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
sin_x_over_x(T const & x)
{
  using S = typename Sacado::ScalarType<T>::type;

  T const
  y = std::abs(x);

  S const
  y_val = Sacado::ScalarValue<T>::eval(y);

  S const
  epsilon2 = std::sqrt(machine_epsilon<T>());

  S const
  epsilon4 = std::sqrt(epsilon2);

  if (y_val > epsilon4) {
    return std::sin(y) / y;
  } else if (y_val > epsilon2) {
    return 1.0 - y * y / 6.0;
  }

  return T(1.0);
}

} // namespace impl

/// \addtogroup minitensor_rotations
/// @{

///
/// Quaternion class. Scalar-first convention: the quaternion
/// \f$ q = (q_s; q_v) \f$ has scalar part \f$ q_s \f$ and vector part
/// \f$ q_v \in R^3 \f$. A unit quaternion encodes the rotation by
/// angle \f$ \theta \f$ about the unit axis \f$ n \f$ as
/// \f$ q_s = \cos(\theta/2) \f$, \f$ q_v = \sin(\theta/2)\, n \f$.
///
template<typename T>
class Quaternion
{
public:

  ///
  /// Default construction: identity quaternion \f$ (1; 0, 0, 0) \f$.
  ///
  KOKKOS_INLINE_FUNCTION
  Quaternion() : s_(1.0), v_(Filler::ZEROS)
  {
  }

  ///
  /// Construction from scalar and vector parts.
  /// \param s scalar part
  /// \param v vector part
  ///
  KOKKOS_INLINE_FUNCTION
  Quaternion(T const & s, Vector<T, 3> const & v) : s_(s), v_(v)
  {
  }

  ///
  /// Construction from four components, scalar first.
  /// \param s scalar part
  /// \param v0 first component of the vector part
  /// \param v1 second component of the vector part
  /// \param v2 third component of the vector part
  ///
  KOKKOS_INLINE_FUNCTION
  Quaternion(T const & s, T const & v0, T const & v1, T const & v2) :
      s_(s), v_(v0, v1, v2)
  {
  }

  ///
  /// Construction from a 4-vector in scalar-first order.
  /// \param q 4-vector \f$ (q_s, q_{v0}, q_{v1}, q_{v2}) \f$
  ///
  explicit
  KOKKOS_INLINE_FUNCTION
  Quaternion(Vector<T, 4> const & q) : s_(q(0)), v_(q(1), q(2), q(3))
  {
  }

  ///
  /// Scalar part.
  ///
  KOKKOS_INLINE_FUNCTION
  T const &
  scalar() const
  {
    return s_;
  }

  ///
  /// Scalar part, mutable.
  ///
  KOKKOS_INLINE_FUNCTION
  T &
  scalar()
  {
    return s_;
  }

  ///
  /// Vector part.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, 3> const &
  vector() const
  {
    return v_;
  }

  ///
  /// Vector part, mutable.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, 3> &
  vector()
  {
    return v_;
  }

  ///
  /// Indexed access, scalar first.
  /// \param i component index: 0 is the scalar part, 1..3 the vector part
  ///
  KOKKOS_INLINE_FUNCTION
  T const &
  operator()(Index const i) const
  {
    assert(i < 4);
    return i == 0 ? s_ : v_(i - 1);
  }

  ///
  /// Indexed access, scalar first, mutable.
  /// \param i component index: 0 is the scalar part, 1..3 the vector part
  ///
  KOKKOS_INLINE_FUNCTION
  T &
  operator()(Index const i)
  {
    assert(i < 4);
    return i == 0 ? s_ : v_(i - 1);
  }

  ///
  /// Components as a 4-vector in scalar-first order.
  ///
  KOKKOS_INLINE_FUNCTION
  Vector<T, 4>
  to_vector() const
  {
    Vector<T, 4>
    q(4);

    q(0) = s_;
    q(1) = v_(0);
    q(2) = v_(1);
    q(3) = v_(2);

    return q;
  }

private:

  ///
  /// Scalar part.
  ///
  T
  s_{};

  ///
  /// Vector part.
  ///
  Vector<T, 3>
  v_;
};

///
/// Hamilton product of two quaternions
/// \param p quaternion
/// \param q quaternion
/// \return \f$ p q = (p_s q_s - p_v \cdot q_v;
/// \; p_s q_v + q_s p_v + p_v \times q_v) \f$
///
template<typename T>
KOKKOS_INLINE_FUNCTION
Quaternion<T>
operator*(Quaternion<T> const & p, Quaternion<T> const & q);

///
/// Quaternion equality, tested component-wise.
/// \param p quaternion
/// \param q quaternion
///
template<typename T>
KOKKOS_INLINE_FUNCTION
bool
operator==(Quaternion<T> const & p, Quaternion<T> const & q);

///
/// Quaternion inequality, tested component-wise.
/// \param p quaternion
/// \param q quaternion
///
template<typename T>
KOKKOS_INLINE_FUNCTION
bool
operator!=(Quaternion<T> const & p, Quaternion<T> const & q);

///
/// Quaternion conjugate
/// \param q quaternion
/// \return \f$ \bar{q} = (q_s; -q_v) \f$
///
template<typename T>
KOKKOS_INLINE_FUNCTION
Quaternion<T>
conjugate(Quaternion<T> const & q);

///
/// Quaternion norm
/// \param q quaternion
/// \return \f$ \sqrt{q_s^2 + q_v \cdot q_v} \f$
///
template<typename T>
KOKKOS_INLINE_FUNCTION
T
norm(Quaternion<T> const & q);

///
/// Quaternion inverse
/// \param q quaternion
/// \return \f$ q^{-1} = \bar{q} / |q|^2 \f$
///
template<typename T>
KOKKOS_INLINE_FUNCTION
Quaternion<T>
inverse(Quaternion<T> const & q);

///
/// Quaternion normalized to unit norm
/// \param q quaternion
/// \return \f$ q / |q| \f$
///
template<typename T>
KOKKOS_INLINE_FUNCTION
Quaternion<T>
unit(Quaternion<T> const & q);

///
/// Quaternion output, scalar-first components separated by commas.
/// \param os output stream
/// \param q quaternion
/// \return os
///
template<typename T>
std::ostream &
operator<<(std::ostream & os, Quaternion<T> const & q);

///
/// Quaternion of a rotation matrix, using Spurrier's singularity-free
/// algorithm. Selects the largest of the trace and the diagonal
/// entries of the matrix and computes the quaternion components
/// accordingly, which ensures a stable evaluation for all rotation
/// angles including \f$ \theta = \pi \f$.
/// \param R rotation matrix with \f$ R \in SO(3) \f$
/// \return unit quaternion \f$ q \f$ such that rt_of_q(q) recovers R
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Quaternion<T>
q_of_rt(Tensor<T, N> const & R);

///
/// Principal rotation vector of a quaternion. The sign of the
/// quaternion is first normalized so that \f$ q_s \geq 0 \f$, which
/// selects the principal rotation \f$ |rv| \leq \pi \f$. Small
/// rotations use an alternate evaluation to avoid numerical
/// instability.
/// \param q unit quaternion
/// \return rotation vector \f$ rv \f$ with \f$ |rv| \leq \pi \f$
///
template<typename T>
KOKKOS_INLINE_FUNCTION
Vector<T, 3>
rv_of_q(Quaternion<T> const & q);

///
/// Quaternion of a rotation vector, computed as
/// \f$ q_s = \cos(|rv|/2) \f$,
/// \f$ q_v = \frac{1}{2}\,\psi(|rv|/2)\, rv \f$ with
/// \f$ \psi(x) = \sin(x)/x \f$, using asymptotic expansions of
/// \f$ \psi \f$ for small angles.
/// \param rv rotation vector
/// \return unit quaternion
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Quaternion<T>
q_of_rv(Vector<T, N> const & rv);

///
/// Rotation matrix of a quaternion, based on the identity
/// \f$ R = 2 q_v \otimes q_v + 2 q_s \mathrm{skew}(q_v) + (2 q_s^2 - 1) I \f$.
/// \param q unit quaternion
/// \return rotation matrix \f$ R \in SO(3) \f$
///
template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor<T, 3>
rt_of_q(Quaternion<T> const & q);

///
/// Principal rotation vector of a rotation matrix, computed through
/// the quaternion representation as rv_of_q(q_of_rt(R)). Unlike
/// log_rotation, this path is singularity-free for all rotation
/// angles including \f$ \theta = \pi \f$.
/// \param R rotation matrix with \f$ R \in SO(3) \f$
/// \return rotation vector \f$ rv \f$ with \f$ |rv| \leq \pi \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, 3>
rv_of_rt(Tensor<T, N> const & R);

///
/// Rotation matrix of a rotation vector, computed through the
/// quaternion representation as rt_of_q(q_of_rv(rv)).
/// \param rv rotation vector
/// \return rotation matrix \f$ R \in SO(3) \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, 3>
rt_of_rv(Vector<T, N> const & rv);

///
/// Rotation vector equivalent to old (representing the same rotation)
/// but as close as possible to prev. Resolves the \f$ 2 \pi k \f$
/// ambiguity of the rotation-vector representation, which enforces
/// continuity in incremental rotation updates such as time stepping.
/// \param old rotation vector to continue
/// \param prev previous rotation vector
/// \return rotation vector equivalent to old and closest to prev
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
rv_continue(Vector<T, N> const & old, Vector<T, N> const & prev);

//
// Hamilton product of two quaternions.
//
template<typename T>
KOKKOS_INLINE_FUNCTION
Quaternion<T>
operator*(Quaternion<T> const & p, Quaternion<T> const & q)
{
  T const
  s = p.scalar() * q.scalar() - dot(p.vector(), q.vector());

  Vector<T, 3> const
  v = p.scalar() * q.vector() + q.scalar() * p.vector() +
      cross(p.vector(), q.vector());

  return Quaternion<T>(s, v);
}

//
// Quaternion equality.
//
template<typename T>
KOKKOS_INLINE_FUNCTION
bool
operator==(Quaternion<T> const & p, Quaternion<T> const & q)
{
  return p.scalar() == q.scalar() && p.vector() == q.vector();
}

//
// Quaternion inequality.
//
template<typename T>
KOKKOS_INLINE_FUNCTION
bool
operator!=(Quaternion<T> const & p, Quaternion<T> const & q)
{
  return !(p == q);
}

//
// Quaternion conjugate.
//
template<typename T>
KOKKOS_INLINE_FUNCTION
Quaternion<T>
conjugate(Quaternion<T> const & q)
{
  return Quaternion<T>(q.scalar(), -q.vector());
}

//
// Quaternion norm.
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
norm(Quaternion<T> const & q)
{
  return std::sqrt(
      q.scalar() * q.scalar() + dot(q.vector(), q.vector()));
}

//
// Quaternion inverse.
//
template<typename T>
KOKKOS_INLINE_FUNCTION
Quaternion<T>
inverse(Quaternion<T> const & q)
{
  T const
  norm_square = q.scalar() * q.scalar() + dot(q.vector(), q.vector());

  return Quaternion<T>(q.scalar() / norm_square, -q.vector() / norm_square);
}

//
// Quaternion normalized to unit norm.
//
template<typename T>
KOKKOS_INLINE_FUNCTION
Quaternion<T>
unit(Quaternion<T> const & q)
{
  T const
  norm_q = norm(q);

  return Quaternion<T>(q.scalar() / norm_q, q.vector() / norm_q);
}

//
// Quaternion output.
//
template<typename T>
std::ostream &
operator<<(std::ostream & os, Quaternion<T> const & q)
{
  os << std::scientific << std::setprecision(17);

  os << std::setw(24) << q.scalar();

  for (Index i = 0; i < 3; ++i) {
    os << "," << std::setw(24) << q.vector()(i);
  }

  return os;
}

//
// Quaternion of a rotation matrix, Spurrier's algorithm.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Quaternion<T>
q_of_rt(Tensor<T, N> const & R)
{
  assert(R.get_dimension() == 3);

  //firewalls, make sure R \in SO(3)
  assert(norm(dot_t(R, R) - eye<T, N>(3)) <
         std::max(1.0e-12 * norm(R), 1.0e-12));
  assert(std::abs(det(R) - 1.0) <
         std::max(1.0e-12 * norm(R), 1.0e-12));

  using S = typename Sacado::ScalarType<T>::type;

  T const
  trace_R = trace(R);

  // Select the largest of the trace and the diagonal entries. The
  // sentinel index 3 denotes the trace.
  S
  max_value = Sacado::ScalarValue<T>::eval(trace_R);

  Index
  max_index = 3;

  for (Index i = 0; i < 3; ++i) {

    S const
    diagonal_value = Sacado::ScalarValue<T>::eval(R(i, i));

    if (diagonal_value > max_value) {
      max_value = diagonal_value;
      max_index = i;
    }
  }

  Quaternion<T>
  q;

  switch (max_index) {

  case 3:
  {
    T const
    root = std::sqrt(trace_R + 1.0);

    T const
    factor = 0.5 / root;

    q = Quaternion<T>(
        0.5 * root,
        factor * (R(2, 1) - R(1, 2)),
        factor * (R(0, 2) - R(2, 0)),
        factor * (R(1, 0) - R(0, 1)));
  }
  break;

  case 2:
  {
    T const
    root = std::sqrt(2.0 * R(2, 2) + 1.0 - trace_R);

    T const
    factor = 0.5 / root;

    q = Quaternion<T>(
        factor * (R(1, 0) - R(0, 1)),
        factor * (R(0, 2) + R(2, 0)),
        factor * (R(1, 2) + R(2, 1)),
        0.5 * root);
  }
  break;

  case 1:
  {
    T const
    root = std::sqrt(2.0 * R(1, 1) + 1.0 - trace_R);

    T const
    factor = 0.5 / root;

    q = Quaternion<T>(
        factor * (R(0, 2) - R(2, 0)),
        factor * (R(0, 1) + R(1, 0)),
        0.5 * root,
        factor * (R(1, 2) + R(2, 1)));
  }
  break;

  case 0:
  {
    T const
    root = std::sqrt(2.0 * R(0, 0) + 1.0 - trace_R);

    T const
    factor = 0.5 / root;

    q = Quaternion<T>(
        factor * (R(2, 1) - R(1, 2)),
        0.5 * root,
        factor * (R(0, 1) + R(1, 0)),
        factor * (R(0, 2) + R(2, 0)));
  }
  break;

  }

  return q;
}

//
// Principal rotation vector of a quaternion.
//
template<typename T>
KOKKOS_INLINE_FUNCTION
Vector<T, 3>
rv_of_q(Quaternion<T> const & q)
{
  // Normalize the sign so that the scalar part is non-negative,
  // which selects the principal rotation |rv| <= pi.
  bool const
  negative = Sacado::ScalarValue<T>::eval(q.scalar()) < 0.0;

  T const
  qs = negative == true ? T(-q.scalar()) : q.scalar();

  Vector<T, 3>
  qv(q.vector()(0), q.vector()(1), q.vector()(2));

  if (negative == true) {
    for (Index i = 0; i < 3; ++i) {
      qv(i) = -qv(i);
    }
  }

  T const
  qv_norm = norm(qv);

  // Small rotation: rv = 2 qv to leading order.
  if (Sacado::ScalarValue<T>::eval(qv_norm) <
      std::sqrt(machine_epsilon<T>())) {
    return 2.0 * qv;
  }

  T
  rv_norm;

  // asin is accurate for small vector parts, acos otherwise.
  if (Sacado::ScalarValue<T>::eval(qv_norm) < std::sqrt(0.5)) {
    rv_norm = 2.0 * std::asin(qv_norm);
  } else {
    rv_norm = 2.0 * std::acos(qs);
  }

  return rv_norm / qv_norm * qv;
}

//
// Quaternion of a rotation vector.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Quaternion<T>
q_of_rv(Vector<T, N> const & rv)
{
  assert(rv.get_dimension() == 3);

  T const
  half_norm = 0.5 * norm(rv);

  T const
  factor = 0.5 * impl::sin_x_over_x(half_norm);

  return Quaternion<T>(
      std::cos(half_norm),
      factor * rv(0),
      factor * rv(1),
      factor * rv(2));
}

//
// Rotation matrix of a quaternion.
//
template<typename T>
KOKKOS_INLINE_FUNCTION
Tensor<T, 3>
rt_of_q(Quaternion<T> const & q)
{
  T const &
  qs = q.scalar();

  Vector<T, 3> const &
  qv = q.vector();

  Tensor<T, 3> const
  R = 2.0 * dyad(qv, qv) + 2.0 * qs * skew(qv) +
      (2.0 * qs * qs - 1.0) * identity<T, 3>(3);

  return R;
}

//
// Principal rotation vector of a rotation matrix.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, 3>
rv_of_rt(Tensor<T, N> const & R)
{
  return rv_of_q(q_of_rt(R));
}

//
// Rotation matrix of a rotation vector.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, 3>
rt_of_rv(Vector<T, N> const & rv)
{
  return rt_of_q(q_of_rv(rv));
}

//
// Rotation vector equivalent to old but closest to prev.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
rv_continue(Vector<T, N> const & old, Vector<T, N> const & prev)
{
  assert(old.get_dimension() == 3);
  assert(prev.get_dimension() == 3);

  using S = typename Sacado::ScalarType<T>::type;

  T const
  norm_old = norm(old);

  Vector<T, N>
  direction(old.get_dimension());

  T
  projection;

  if (Sacado::ScalarValue<T>::eval(norm_old) > 0.0) {
    direction = old / norm_old;
    projection = dot(direction, prev);
  } else {
    projection = norm(prev);
    if (Sacado::ScalarValue<T>::eval(projection) == 0.0) {
      return old;
    }
    direction = prev / projection;
  }

  if (Sacado::ScalarValue<T>::eval(projection) == 0.0) {
    return old;
  }

  S const
  pi = std::acos(S(-1.0));

  // The equivalent rotation vectors are (|old| + 2 pi k) direction;
  // their signed coordinate along direction closest to that of prev
  // determines the number of turns k.
  S const
  number_turns = std::round(
      0.5 *
      (Sacado::ScalarValue<T>::eval(projection) -
       Sacado::ScalarValue<T>::eval(norm_old)) / pi);

  if (number_turns == 0.0) {
    return old;
  }

  return (2.0 * number_turns * pi + norm_old) * direction;
}

/// @}
} // namespace minitensor

#endif // MiniTensor_Quaternion_h
