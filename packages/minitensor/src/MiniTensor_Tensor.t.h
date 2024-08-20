// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Tensor_t_h)
#define MiniTensor_Tensor_t_h

namespace minitensor {

//
// tensor input
//
template<typename T, Index N>
std::istream &
operator>>(std::istream & is, Tensor<T, N> & A)
{

  Index const
  dimension = A.get_dimension();

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      is >> A(i,j);
    }
  }

  return is;
}

//
// tensor output
//
template<typename T, Index N>
std::ostream &
operator<<(std::ostream & os, Tensor<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  if (dimension == 0) {
    return os;
  }

  os << std::scientific << std::setprecision(17);

  for (Index i = 0; i < dimension; ++i) {

    os << std::setw(24) << A(i,0);

    for (Index j = 1; j < dimension; ++j) {
      os << "," << std::setw(24) << A(i,j);
    }

    os << std::endl;
  }

  return os;
}

} // namespace minitensor

#endif // MiniTensor_Tensor_t_h
