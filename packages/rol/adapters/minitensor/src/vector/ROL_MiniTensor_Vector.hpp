// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(ROL_MiniTensor_Vector_hpp)
#define ROL_MiniTensor_Vector_hpp

#include "MiniTensor.h"
#include "ROL_Vector.hpp"

namespace ROL {

///
/// ROL container for MiniTensor Vector
///
template<typename T, minitensor::Index N>
class MiniTensorVector : public Vector<T> {

public:

  MiniTensorVector(minitensor::Vector<T, N> const & v);

  virtual
  ~MiniTensorVector();

  //
  // ROL's interface
  //
  virtual
  void
  set(Vector<T> const & x) final;

  virtual
  void
  plus(Vector<T> const & x) final;

  virtual
  void
  axpy(T const alpha, Vector<T> const & x) final;

  virtual
  void
  scale(T const alpha) final;

  virtual
  T
  dot(Vector<T> const & x) const final;

  virtual
  T
  norm() const final;

  virtual
  ROL::Ptr<Vector<T>>
  clone() const final;

  virtual
  ROL::Ptr<Vector<T>>
  basis(int const i) const final;

  virtual
  int
  dimension() const final;

  virtual
  void
  applyUnary(Elementwise::UnaryFunction<T> const & f) final;

  virtual
  void
  applyBinary(Elementwise::BinaryFunction<T> const & f, Vector<T> const & x) final;

  virtual
  T
  reduce(Elementwise::ReductionOp<T> const & r) const final;

  //
  // Utilities
  //
  void
  set(minitensor::Vector<T, N> const & x);

  minitensor::Vector<T, N>
  getVector() const;

  minitensor::Vector<T, N>
  getVector();

  friend
  std::istream &
  operator>>(std::istream & is, MiniTensorVector & u)
  {
    is >> u.vector_;
    return is;
  }

  friend
  std::ostream &
  operator<<(std::ostream & os, MiniTensorVector const & u)
  {
    os << u.vector_;
    return os;
  }

private:

  minitensor::Vector<T, N>
  vector_;
}; // class MiniTensorVector

///
/// Covert from ROL to MiniTensor
///
template<typename T, minitensor::Index N>
minitensor::Vector<T, N>
MTfromROL(Vector<T> const & x);

//
// Convert from MiniTensor to ROL
//
template<typename T, minitensor::Index N>
void
MTtoROL(minitensor::Vector<T, N> const & xval, Vector<T> & x);

} // namespace ROL

#include "ROL_MiniTensor_Vector_Def.hpp"

#endif // ROL_MiniTensor_Vector_hpp
