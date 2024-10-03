// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Matrix_t_h)
#define MiniTensor_Matrix_t_h

namespace minitensor {

//
// Matrix input
//
template<typename T, Index M, Index N>
std::istream &
operator>>(std::istream & is, Matrix<T, M, N> & A)
{
  Index const
  num_rows = A.get_num_rows();

  Index const
  num_cols = A.get_num_cols();

  for (Index i = 0; i < num_rows; ++i) {
    for (Index j = 0; j < num_cols; ++j) {
      is >> A(i,j);
    }
  }

  return is;
}

//
// Matrix output
//
template<typename T, Index M, Index N>
std::ostream &
operator<<(std::ostream & os, Matrix<T, M, N> const & A)
{
  Index const
  num_rows = A.get_num_rows();

  Index const
  num_cols = A.get_num_cols();

  Index const
  dimension = num_rows * num_cols;

  if (dimension == 0) {
    return os;
  }

  os << std::scientific << std::setprecision(17);

  for (Index i = 0; i < num_rows; ++i) {

    os << std::setw(24) << A(i,0);

    for (Index j = 1; j < num_cols; ++j) {
      os << "," << std::setw(24) << A(i,j);
    }

    os << std::endl;
  }

  return os;
}

} // namespace minitensor

#endif // MiniTensor_Matrix_t_h
