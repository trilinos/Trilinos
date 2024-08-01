// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Vector_t_h)
#define MiniTensor_Vector_t_h

namespace minitensor {

//
// R^N vector input
// \param u vector
// \param is input stream
// \return is input stream
//
template<typename T, Index N>
std::istream &
operator>>(std::istream & is, Vector<T, N> & u)
{
  Index const
  dimension = u.get_dimension();

  for (Index i = 0; i < dimension; ++i) {
    is >> u(i);
  }

  return is;
}

//
// R^N vector output
// \param u vector
// \param os output stream
// \return os output stream
//
template<typename T, Index N>
std::ostream &
operator<<(std::ostream & os, Vector<T, N> const & u)
{
  Index const
  dimension = u.get_dimension();

  if (dimension == 0) {
    return os;
  }

  os << std::scientific << std::setprecision(17);

  os << std::setw(24) << u(0);

  for (Index i = 1; i < dimension; ++i) {
    os << "," << std::setw(24) << u(i);
  }

  return os;
}

} // namespace minitensor

#endif // MiniTensor_Vector_t_h
