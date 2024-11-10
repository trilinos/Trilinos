// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Tensor3_t_h)
#define MiniTensor_Tensor3_t_h

namespace minitensor {

//
// \return \f$ B = A : u := B_i = A_{ijk} u_{jk} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Vector<typename Promote<S, T>::type, N>
dotdot(Tensor3<T, N> const & A, Tensor<S, N> const & u)
{
  Index const
  dimension = A.get_dimension();

  assert(u.get_dimension() == dimension);

  Vector<typename Promote<S, T>::type, N>
  B(dimension);

  for (Index i = 0; i < dimension; ++i) {

    typename Promote<S, T>::type
    s = 0.0;

    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        s += A(i,j,k) * u(j,k);
      }
    }
    B(i) = s;
  }

  return B;
}

//
// \return \f$ B = A \cdot u := B_{ij} = A_{ijp} u_p \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dot(Tensor3<T, N> const & A, Vector<S, N> const & u)
{
  Index const
  dimension = A.get_dimension();

  assert(u.get_dimension() == dimension);

  Tensor<typename Promote<S, T>::type, N>
  B(N);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {

      typename Promote<S, T>::type
      s = 0.0;

      for (Index p = 0; p < dimension; ++p) {
        s += A(i,j,p) * u(p);
      }
      B(i,j) = s;
    }
  }

  return B;
}

//
// \return \f$ B = u \cdot A := B_{ij} = u_p A{pij} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dot(Vector<S, N> const & u, Tensor3<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  assert(u.get_dimension() == dimension);

  Tensor<typename Promote<S, T>::type, N>
  B(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {

      typename Promote<S, T>::type
      s = 0.0;

      for (Index p = 0; p < dimension; ++p) {
        s += u(p) * A(p,i,j);
      }
      B(i,j) = s;
    }
  }

  return B;
}


//
// \return \f$ B = A \cdot u := B_{ij} = A_{ipj} u_p \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dot2(Tensor3<T, N> const & A, Vector<S, N> const & u)
{
  Index const
  dimension = A.get_dimension();

  assert(u.get_dimension() == dimension);

  Tensor<typename Promote<S, T>::type, N>
  B(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {

      typename Promote<S, T>::type
      s = 0.0;

      for (Index p = 0; p < dimension; ++p) {
        s += A(i,p,j) * u(p);
      }
      B(i,j) = s;
    }
  }

  return B;
}

//
// \return \f$ B = u \cdot A := B_{ij} = u_p A_{ipj} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor<typename Promote<S, T>::type, N>
dot2(Vector<S, N> const & u, Tensor3<T, N> const & A)
{
  return dot2(A, u);
}

//
// \return \f$ C = A \cdot B := C_{ijk} = A_{ijp} B_{pk} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
dot(Tensor3<T, N> const & A, Tensor<S, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor3<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index k = 0; k < dimension; ++k) {
      for (Index j = 0; j < dimension; ++j) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < dimension; ++p) {
          s += A(i,j,p) * B(p,k);
        }
        C(i,j,k) = s;
      }
    }
  }

  return C;
}

//
// \return \f$ C = A \cdot B := C_{ijk} = A_{ip} B_{pjk} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
dot(Tensor<S, N> const & A, Tensor3<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor3<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index k = 0; k < dimension; ++k) {
      for (Index j = 0; j < dimension; ++j) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < dimension; ++p) {
          s += A(i,p) * B(p,j,k);
        }
        C(i,j,k) = s;
      }
    }
  }

  return C;
}

//
// \return \f$ C = A \cdot B := C_{ijk} = A_{ipj} B_{pk} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
dot2(Tensor3<T, N> const & A, Tensor<S, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor3<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index k = 0; k < dimension; ++k) {
      for (Index j = 0; j < dimension; ++j) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < dimension; ++p) {
          s += A(i,p,j) * B(p,k);
        }
        C(i,j,k) = s;
      }
    }
  }

  return C;
}


//
// \return \f$ C = A \cdot B := C_{ijk} = A_{ip} B_{jpk} \f$
//
template<typename S, typename T, Index N>
KOKKOS_INLINE_FUNCTION
Tensor3<typename Promote<S, T>::type, N>
dot2(Tensor<S, N> const & A, Tensor3<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor3<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index k = 0; k < dimension; ++k) {
      for (Index j = 0; j < dimension; ++j) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < dimension; ++p) {
          s += A(i,p) * B(j,p,k);
        }
        C(i,j,k) = s;
      }
    }
  }

  return C;
}


//
// 3rd-order tensor input
// \param A 3rd-order tensor
// \param is input stream
// \return is input stream
//
template<typename T, Index N>
std::istream &
operator>>(std::istream & is, Tensor3<T, N> & A)
{
  Index const
  dimension = A.get_dimension();

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        is >> A(i,j,k);
      }
    }
  }

  return is;
}

//
// 3rd-order tensor output
// \param A 3rd-order tensor
// \param os output stream
// \return os output stream
//
template<typename T, Index N>
std::ostream &
operator<<(std::ostream & os, Tensor3<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  if (dimension == 0) {
    return os;
  }

  os << std::scientific << std::setprecision(17);

  for (Index i = 0; i < dimension; ++i) {

    for (Index j = 0; j < dimension; ++j) {

      os << std::setw(24) << A(i,j,0);

      for (Index k = 1; k < dimension; ++k) {
        os << "," << std::setw(24) << A(i,j,k);
      }

      os << std::endl;

    }

    os << std::endl;
    os << std::endl;

  }

  return os;
}

} // namespace minitensor

#endif // MiniTensor_Tensor3_t_h
