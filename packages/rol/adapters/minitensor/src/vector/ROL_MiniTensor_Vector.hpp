// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
