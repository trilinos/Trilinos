// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Rotations_h)
#define MiniTensor_Rotations_h

// Rotation algebra: SO(N) log and exp maps.
#include "MiniTensor_MatrixFunctions.h"
#include "MiniTensor_Norms.h"

namespace minitensor {

/// \addtogroup minitensor_rotations
/// @{

///
/// Logarithmic map of a rotation
/// \param R with \f$ R \in SO(3) \f$
/// \return \f$ r = \log R \f$ with \f$ r \in so(3) \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
log_rotation(Tensor<T, N> const & R);

///
/// Logarithmic map of a 180 degree rotation
/// \param R with \f$ R \in SO(3) \f$
/// \return \f$ r = \log R \f$ with \f$ r \in so(3) \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
log_rotation_pi(Tensor<T, N> const & R);

///
/// Exponential map of a skew-symmetric tensor
/// \param r \f$ r \in so(3) \f$
/// \return \f$ R = \exp R \f$ with \f$ R \in SO(3) \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
exp_skew_symmetric(Tensor<T, N> const & r);

///
/// Axial vector of a skew-symmetric tensor, the inverse of
/// skew(Vector). Defined for 3D only.
/// \param W skew-symmetric tensor \f$ W \in so(3) \f$
/// \return \f$ w \f$ such that \f$ W u = w \times u \f$ for all \f$ u \f$
///
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
vee(Tensor<T, N> const & W);

//
// R^N logarithmic map of a rotation.
// \param R with \f$ R \in SO(N) \f$
// \return \f$ r = \log R \f$ with \f$ r \in so(N) \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
log_rotation(Tensor<T, N> const & R)
{
  Index const
  dimension = R.get_dimension();

  //firewalls, make sure R \in SO(N)
  assert(norm(dot_t(R,R) - eye<T, N>(dimension)) <
         std::max(1.0e-12 * norm(R), 1.0e-12));
  assert(std::abs(det(R) - 1.0) <
         std::max(1.0e-12 * norm(R), 1.0e-12));
  // acos requires input between -1 and +1
  T
  cosine = 0.5 * (trace(R) - 1.0);

  if (cosine < -1.0) {
    cosine = -1.0;
  } else if(cosine > 1.0) {
    cosine = 1.0;
  }
  T
  theta = std::acos(cosine);

  Tensor<T, N>
  r(dimension);

  switch (dimension) {

    default:
      MT_ERROR_EXIT("Logarithm of SO(N) N != 2,3 not implemented.");
      break;

    case 3:
      if (theta == 0.0) {

        r = zero<T, N>(3);

      } else if (std::abs(cosine + 1.0) < 10.0 * machine_epsilon<T>())  {

        r = log_rotation_pi(R);

      } else {

        r = theta / std::sin(theta) * skew(R);

      }
      break;

    case 2:
      r(0,0) = 0.0;
      r(0,1) = -theta;
      r(1,0) = theta;
      r(1,1) = 0.0;
      break;

    case 1:
      r(0,0) = 0.0;
      break;

  }

  return r;
}

// R^N Logarithmic map of a 180-degree rotation.
// \param R with \f$ R \in SO(N) \f$
// \return \f$ r = \log R \f$ with \f$ r \in so(N) \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
log_rotation_pi(Tensor<T, N> const & R)
{
  Index const
  dimension = R.get_dimension();

  // set firewall to make sure the rotation is indeed 180 degrees
  assert(std::abs(trace(R) + 1.0) < 10.0 * machine_epsilon<T>());

  Tensor<T, N>
  r(dimension);

  switch (dimension) {

    default:
      MT_ERROR_EXIT("Logarithm of SO(N) N != 2,3 not implemented.");
      break;

    case 3:
    {
      Vector<T, N>
      normal(3);

      Tensor<T, N> const
      B = R - identity<T, N>(3);

      Vector<T, N> const
      u = row(B, 0);

      Vector<T, N> const
      v = row(B, 1);

      normal = cross(u, v);

      if (norm(normal) < machine_epsilon<T>()) {

        Vector<T, N> const
        w = row(B, 2);

        normal = cross(u, w);

        if (norm(normal) < machine_epsilon<T>()) {
          MT_ERROR_EXIT("Cannot determine rotation vector of rotation.");
        }

      }

      normal = unit(normal);

      r.fill(Filler::ZEROS);
      r(0,1) = -normal(2);
      r(0,2) =  normal(1);
      r(1,0) =  normal(2);
      r(1,2) = -normal(0);
      r(2,0) = -normal(1);
      r(2,1) =  normal(0);

      T const
      pi = std::acos(-1.0);

      r = pi * r;
    }
    break;

    case 2:
    {
      T theta = std::acos(-1.0);
      if (R(0,0) > 0.0) {
        theta = -theta;
      }

      r(0,0) = 0.0;
      r(0,1) = -theta;
      r(1,0) = theta;
      r(1,1) = 0.0;
    }
    break;

  }

  return r;
}

//
// R^N exponential map of a skew-symmetric tensor.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
exp_skew_symmetric(Tensor<T, N> const & r)
{
  // Check whether skew-symmetry holds
  assert(norm(sym(r)) < std::max(1.0e-12 * norm(r), 1.0e-12));

  Index const
  dimension = r.get_dimension();

  Tensor<T, N>
  R = identity<T, N>(dimension);

  T
  theta = 0.0;

  switch (dimension) {

    default:
      R = exp(r);
      break;

    case 3:
      theta = std::sqrt(r(2,1)*r(2,1)+r(0,2)*r(0,2)+r(1,0)*r(1,0));

      //Check whether norm == 0. If so, return identity.
      if (theta >= machine_epsilon<T>()) {
        R += sin(theta) / theta * r +
            (1.0 - cos(theta)) / (theta * theta) * r * r;
      }
      break;

    case 2:
      theta = r(1,0);

      {
        T const
        c = std::cos(theta);

        T const
        s = std::sin(theta);

        R(0,0) = c;
        R(0,1) = -s;
        R(1,0) = s;
        R(1,1) = c;
      }

      break;

    case 1:
      R(0,0) = 1.0;
      break;

 }

  return R;
}

//
// R^N axial vector of a skew-symmetric tensor, undefined for N != 3.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
vee(Tensor<T, N> const & W)
{
  // Check whether skew-symmetry holds
  assert(norm(sym(W)) < std::max(1.0e-12 * norm(W), 1.0e-12));

  Index const
  dimension = W.get_dimension();

  Vector<T, N>
  w(dimension);

  switch (dimension) {

  case 3:
    w(0) = W(2, 1);
    w(1) = W(0, 2);
    w(2) = W(1, 0);
    break;

  default:
    MT_ERROR_EXIT("Axial vector from tensor defined for 3D only");
    break;

  }

  return w;
}

/// @}
} // namespace minitensor

#endif // MiniTensor_Rotations_h
