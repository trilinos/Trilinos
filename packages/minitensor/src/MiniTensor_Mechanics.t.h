// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Mechanics_t_h)
#define MiniTensor_Mechanics_t_h

namespace minitensor {

//
// Push forward covariant vector
// \param \f$ F, u \f$
// \return \f$ F^{-T} u \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
push_forward_covariant(Tensor<T, N> const & F, Vector<T, N> const & u)
{
  Index const
  dimension = F.get_dimension();

  Vector<T, N>
  v(dimension);

  T const
  J = det(F);

  assert(J > 0.0);

  switch (dimension) {

    default:
      MT_ERROR_EXIT("Supports only 2D and 3D.");
      break;

    case 3:
      v(0) = (
          (-F(1,2)*F(2,1) + F(1,1)*F(2,2)) * u(0) +
          ( F(1,2)*F(2,0) - F(1,0)*F(2,2)) * u(1) +
          (-F(1,1)*F(2,0) + F(1,0)*F(2,1)) * u(2)) / J;

      v(1) = (
          ( F(0,2)*F(2,1) - F(0,1)*F(2,2)) * u(0) +
          (-F(0,2)*F(2,0) + F(0,0)*F(2,2)) * u(1) +
          ( F(0,1)*F(2,0) - F(0,0)*F(2,1)) * u(2)) / J;

      v(2) = (
          (-F(0,2)*F(1,1) + F(0,1)*F(1,2)) * u(0) +
          ( F(0,2)*F(1,0) - F(0,0)*F(1,2)) * u(1) +
          (-F(0,1)*F(1,0) + F(0,0)*F(1,1)) * u(2)) / J;

      break;

    case 2:
      v(0) = ( F(1,1) * u(0) - F(1,0) * u(1)) / J;
      v(1) = (-F(0,1) * u(0) + F(0,0) * u(1)) / J;
      break;

  }

  return v;
}

//
// Pull back covariant vector
// \param \f$ F, v \f$
// \return \f$ F^T v \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
pull_back_covariant(Tensor<T, N> const & F, Vector<T, N> const & u)
{
  Index const
  dimension = F.get_dimension();

  Vector<T, N>
  v(dimension);

  switch (dimension) {

    default:
      MT_ERROR_EXIT("Supports only 2D and 3D.");
      break;

    case 3:
      v(0) = F(0,0) * u(0) + F(1,0) * u(1) + F(2,0) * u(2);
      v(1) = F(0,1) * u(0) + F(1,1) * u(1) + F(2,1) * u(2);
      v(2) = F(0,2) * u(0) + F(1,2) * u(1) + F(2,2) * u(2);

      break;

    case 2:
      v(0) = F(0,0) * u(0) + F(1,0) * u(1);
      v(1) = F(0,1) * u(0) + F(1,1) * u(1);

      break;

  }

  return v;
}

//
// Push forward contravariant vector
// \param \f$ F, u \f$
// \return \f$ F u \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
push_forward_contravariant(Tensor<T, N> const & F, Vector<T, N> const & u)
{
  Index const
  dimension = F.get_dimension();

  Vector<T, N>
  v(dimension);

  switch (dimension) {

    default:
      MT_ERROR_EXIT("Supports only 2D and 3D.");
      break;

    case 3:
      v(0) = F(0,0) * u(0) + F(0,1) * u(1) + F(0,2) * u(2);
      v(1) = F(1,0) * u(0) + F(1,1) * u(1) + F(1,2) * u(2);
      v(2) = F(2,0) * u(0) + F(2,1) * u(1) + F(2,2) * u(2);

      break;

    case 2:
      v(0) = F(0,0) * u(0) + F(0,1) * u(1);
      v(1) = F(1,0) * u(0) + F(1,1) * u(1);

      break;

  }

  return v;
}

//
// Pull back contravariant vector
// \param \f$ F, u \f$
// \return \f$ F^{-1} u \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
pull_back_contravariant(Tensor<T, N> const & F, Vector<T, N> const & u)
{
  Index const
  dimension = F.get_dimension();

  Vector<T, N>
  v(dimension);

  T const
  J = det(F);

  assert(J > 0.0);

  switch (dimension) {

    default:
      MT_ERROR_EXIT("Supports only 2D and 3D.");
      break;

    case 3:
      v(0) = (
          (-F(1,2)*F(2,1) + F(1,1)*F(2,2)) * u(0) +
          ( F(0,2)*F(2,1) - F(0,1)*F(2,2)) * u(1) +
          (-F(0,2)*F(1,1) + F(0,1)*F(1,2)) * u(2)) / J;

      v(1) = (
          ( F(1,2)*F(2,0) - F(1,0)*F(2,2)) * u(0) +
          (-F(0,2)*F(2,0) + F(0,0)*F(2,2)) * u(1) +
          ( F(0,2)*F(1,0) - F(0,0)*F(1,2)) * u(2)) / J;

      v(2) = (
          (-F(1,1)*F(2,0) + F(1,0)*F(2,1)) * u(0) +
          ( F(0,1)*F(2,0) - F(0,0)*F(2,1)) * u(1) +
          (-F(0,1)*F(1,0) + F(0,0)*F(1,1)) * u(2)) / J;

      break;

    case 2:
      v(0) = ( F(1,1) * u(0) - F(0,1) * u(1)) / J;
      v(1) = (-F(1,0) * u(0) + F(0,0) * u(1)) / J;
      break;

  }

  return v;
}

//
// Push forward covariant tensor
// \param \f$ F, A \f$
// \return \f$ F^{-T} A F^{-1} \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
push_forward_covariant(Tensor<T, N> const & F, Tensor<T, N> const & A)
{
  Index const
  dimension = F.get_dimension();

  Tensor<T, N>
  G(dimension);

  T const
  J = det(F);

  assert(J > 0.0);

  switch (dimension) {

    default:
      MT_ERROR_EXIT("Supports only 2D and 3D.");
      break;

    case 3:
      G(0,0) = (-F(1,2)*F(2,1) + F(1,1)*F(2,2)) / J;
      G(0,1) = ( F(0,2)*F(2,1) - F(0,1)*F(2,2)) / J;
      G(0,2) = (-F(0,2)*F(1,1) + F(0,1)*F(1,2)) / J;

      G(1,0) = ( F(1,2)*F(2,0) - F(1,0)*F(2,2)) / J;
      G(1,1) = (-F(0,2)*F(2,0) + F(0,0)*F(2,2)) / J;
      G(1,2) = ( F(0,2)*F(1,0) - F(0,0)*F(1,2)) / J;

      G(2,0) = (-F(1,1)*F(2,0) + F(1,0)*F(2,1)) / J;
      G(2,1) = ( F(0,1)*F(2,0) - F(0,0)*F(2,1)) / J;
      G(2,2) = (-F(0,1)*F(1,0) + F(0,0)*F(1,1)) / J;
      break;

    case 2:
      G(0,0) =  F(1,1) / J;
      G(0,1) = -F(0,1) / J;

      G(1,0) = -F(1,0) / J;
      G(1,1) =  F(0,0) / J;
      break;

  }

  return t_dot(G, dot(A, G));
}

//
// Pull back covariant tensor
// \param \f$ F, A \f$
// \return \f$ F^T A F\f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
pull_back_covariant(Tensor<T, N> const & F, Tensor<T, N> const & A)
{
  return t_dot(F, dot(A, F));
}

//
// Push forward contravariant tensor
// \param \f$ F, A \f$
// \return \f$ F A F^T \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
push_forward_contravariant(Tensor<T, N> const & F, Tensor<T, N> const & A)
{
  return dot_t(dot(F, A), F);
}

//
// Pull back contravariant tensor
// \param \f$ F, A \f$
// \return \f$ F^{-1} A F^{-T} \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
pull_back_contravariant(Tensor<T, N> const & F, Tensor<T, N> const & A)
{
  Index const
  dimension = F.get_dimension();

  Tensor<T, N>
  G(dimension);

  T const
  J = det(F);

  assert(J > 0.0);

  switch (dimension) {

    default:
      MT_ERROR_EXIT("Supports only 2D and 3D.");
      break;

    case 3:
      G(0,0) = (-F(1,2)*F(2,1) + F(1,1)*F(2,2)) / J;
      G(0,1) = ( F(0,2)*F(2,1) - F(0,1)*F(2,2)) / J;
      G(0,2) = (-F(0,2)*F(1,1) + F(0,1)*F(1,2)) / J;

      G(1,0) = ( F(1,2)*F(2,0) - F(1,0)*F(2,2)) / J;
      G(1,1) = (-F(0,2)*F(2,0) + F(0,0)*F(2,2)) / J;
      G(1,2) = ( F(0,2)*F(1,0) - F(0,0)*F(1,2)) / J;

      G(2,0) = (-F(1,1)*F(2,0) + F(1,0)*F(2,1)) / J;
      G(2,1) = ( F(0,1)*F(2,0) - F(0,0)*F(2,1)) / J;
      G(2,2) = (-F(0,1)*F(1,0) + F(0,0)*F(1,1)) / J;
      break;

    case 2:
      G(0,0) =  F(1,1) / J;
      G(0,1) = -F(0,1) / J;

      G(1,0) = -F(1,0) / J;
      G(1,1) =  F(0,0) / J;
      break;

  }

  return dot_t(dot(G, A), G);
}

//
// Piola transformation for vector
// \param \f$ F, u \f$
// \return \f$ \det F F^{-1} u \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
piola(Tensor<T, N> const & F, Vector<T, N> const & u)
{
  Index const
  dimension = F.get_dimension();

  Vector<T, N>
  v(dimension);

  switch (dimension) {

    default:
      MT_ERROR_EXIT("Supports only 2D and 3D.");
      break;

    case 3:
      v(0) = (
          (-F(1,2)*F(2,1) + F(1,1)*F(2,2)) * u(0) +
          ( F(0,2)*F(2,1) - F(0,1)*F(2,2)) * u(1) +
          (-F(0,2)*F(1,1) + F(0,1)*F(1,2)) * u(2));

      v(1) = (
          ( F(1,2)*F(2,0) - F(1,0)*F(2,2)) * u(0) +
          (-F(0,2)*F(2,0) + F(0,0)*F(2,2)) * u(1) +
          ( F(0,2)*F(1,0) - F(0,0)*F(1,2)) * u(2));

      v(2) = (
          (-F(1,1)*F(2,0) + F(1,0)*F(2,1)) * u(0) +
          ( F(0,1)*F(2,0) - F(0,0)*F(2,1)) * u(1) +
          (-F(0,1)*F(1,0) + F(0,0)*F(1,1)) * u(2));

      break;

    case 2:
      v(0) = ( F(1,1) * u(0) - F(0,1) * u(1));
      v(1) = (-F(1,0) * u(0) + F(0,0) * u(1));
      break;

  }

  return v;
}

//
// Inverse Piola transformation for vector
// \param \f$ F, u \f$
// \return \f$ (\det F)^{-1} F u \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<T, N>
piola_inverse(Tensor<T, N> const & F, Vector<T, N> const & u)
{
  Index const
  dimension = F.get_dimension();

  Vector<T, N>
  v(dimension);

  T const
  J = det(F);

  assert(J > 0.0);

  switch (dimension) {

    default:
      MT_ERROR_EXIT("Supports only 2D and 3D.");
      break;

    case 3:
      v(0) = (F(0,0) * u(0) + F(0,1) * u(1) + F(0,2) * u(2)) / J;
      v(1) = (F(1,0) * u(0) + F(1,1) * u(1) + F(1,2) * u(2)) / J;
      v(2) = (F(2,0) * u(0) + F(2,1) * u(1) + F(2,2) * u(2)) / J;

      break;

    case 2:
      v(0) = (F(0,0) * u(0) + F(0,1) * u(1)) / J;
      v(1) = (F(1,0) * u(0) + F(1,1) * u(1)) / J;

      break;

  }

  return v;
}

//
// Piola transformation for tensor, applied on second
// index. Useful for transforming Cauchy stress to 1PK stress.
// \param \f$ F, \sigma \f$
// \return \f$ \det F \sigma F^{-T} \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
piola(Tensor<T, N> const & F, Tensor<T, N> const & sigma)
{
  Index const
  dimension = F.get_dimension();

  Tensor<T, N>
  G(dimension);

  switch (dimension) {

    default:
      MT_ERROR_EXIT("Supports only 2D and 3D.");
      break;

    case 3:
      G(0,0) = (-F(1,2)*F(2,1) + F(1,1)*F(2,2));
      G(0,1) = ( F(0,2)*F(2,1) - F(0,1)*F(2,2));
      G(0,2) = (-F(0,2)*F(1,1) + F(0,1)*F(1,2));

      G(1,0) = ( F(1,2)*F(2,0) - F(1,0)*F(2,2));
      G(1,1) = (-F(0,2)*F(2,0) + F(0,0)*F(2,2));
      G(1,2) = ( F(0,2)*F(1,0) - F(0,0)*F(1,2));

      G(2,0) = (-F(1,1)*F(2,0) + F(1,0)*F(2,1));
      G(2,1) = ( F(0,1)*F(2,0) - F(0,0)*F(2,1));
      G(2,2) = (-F(0,1)*F(1,0) + F(0,0)*F(1,1));
      break;

    case 2:
      G(0,0) =  F(1,1);
      G(0,1) = -F(0,1);

      G(1,0) = -F(1,0);
      G(1,1) =  F(0,0);
      break;

  }

  return dot_t(sigma, G);
}

//
// Inverse Piola transformation for tensor, applied on second
// index. Useful for transforming 1PK stress to Cauchy stress.
// \param \f$ F, P \f$
// \return \f$ (\det F)^{-1} P F^T \f$
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<T, N>
piola_inverse(Tensor<T, N> const & F, Tensor<T, N> const & P)
{
  T const
  J = det(F);

  assert(J > 0.0);

  return dot_t(P, F) / J;
}

//
// Smallest eigenvalue by inverse iteration.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
T
smallest_eigenvalue(Tensor<T, N> const & A)
{
  Tensor<T, N>
  B = inverse(A);

  T const
  tolerance = machine_epsilon<T>();

  Index const
  dimension = A.get_dimension();

  Vector<T, N>
  v(dimension, Filler::ONES);

  Index const
  maximum_iterations = 128;

  T
  relative_error = 1.0;

  Index
  k = 0;

  while (relative_error > tolerance && k < maximum_iterations) {

    Vector<T, N> const
    w = v;

    v = unit(B * w);

    relative_error = norm(v - w) / norm(w);

    ++k;
  }

  return v * A * v;
}

//
// Check strict ellipticity condition for 4th-order tensor.
// Assume A has major symmetries.
//
template<typename T, Index N>
KOKKOS_INLINE_FUNCTION
bool
check_strict_ellipticity(Tensor4<T, N> const & A)
{
  // Convert to 2nd-order tensor
  Tensor<T, dimension_square<N>::value> const
  B(A);

  // Check bounds for eigenvalues
  T const
  lower_bound = bounds_eigenvalues(B).first;

  if (lower_bound > 0.0) {
    return true;
  }

  // Get eigenvalue closest to zero only
  T const
  smallest_eigenvalue = smallest_eigenvavlue(B);

  if (smallest_eigenvalue > 0.0) {
    return true;
  }

  return false;
}

//
// Check strong ellipticity condition for 4th-order tensor.
// Assume A has major and minor symmetries.
//
template<typename T, Index N>
std::pair<bool, Vector<T, N>>
check_strong_ellipticity(Tensor4<T, N> const & A)
{
  bool
  is_elliptic = true;

  Index const
  dimension = A.get_dimension();

  Vector<T, N>
  eigenvector(dimension, Filler::ONES);

  eigenvector /= dimension;

  Index const
  maximum_iterarions = 128;

  T const
  tolerance = machine_epsilon<T>();

  T
  error = 1.0;

#if defined(KOKKOS_ENABLE_CUDA)
  T
  prev_eigenvalue = DBL_MAX;
#else
  using S = typename Sacado::ScalarType<T>::type;

  T
  prev_eigenvalue = std::numeric_limits<S>::max();
#endif

  T
  curr_eigenvalue = prev_eigenvalue;

  Index
  iteration = 0;

  while (error > tolerance && iteration < maximum_iterarions) {

    Tensor<T, N>
    Q = dot2(eigenvector, dot(A, eigenvector));

    Tensor<T, N>
    V;

    Tensor<T, N>
    D;

    std::tie(V, D) = eig_sym(Q);

    curr_eigenvalue = D(dimension - 1, dimension - 1);

    eigenvector = col(V, dimension - 1);

    error = std::abs(prev_eigenvalue) / std::abs(curr_eigenvalue) - 1.0;

    prev_eigenvalue = curr_eigenvalue;

    ++iteration;
  }

  if (curr_eigenvalue <= 0.0) {
    is_elliptic = false;
  }

  return std::make_pair(is_elliptic, eigenvector);
}

} // namespace minitensor

#endif // MiniTensor_Mechanics_t_h
