// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
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

#if defined(ENABLE_ROL)
#if !defined(Intrepid2_MiniTensor_ROL_Vector_h)
#define Intrepid2_MiniTensor_ROL_Vector_h

#include <Intrepid2_MiniTensor.h>
#include "ROL_Vector.hpp"

namespace ROL
{

template <typename T, Index N>
class MiniTensorVector : public Vector<T> {

  using uint = Intrepid2::Index;

private:

  Intrepid2::Vector<T, N> &
  vector_;

public:

  MiniTensorVector(Intrepid2::Vector<T, N> & v) : vector_(v)
  {
    return;
  }

  void
  set(Vector<T> & x)
  {
    assert(dimension() == x.dimension());

    MiniTensorVector<T, N> const &
    ex = Teuchos::dyn_cast<MiniTensorVector<T, N> const>(x);

    Intrepid2::Vector<T, N> const &
    xval = ex.getVector();

    vector_ = xval;
  }

  void
  plus(Vector<T> & x)
  {
    assert(dimension() == x.dimension());

    MiniTensorVector<T, N> const &
    ex = Teuchos::dyn_cast<MiniTensorVector<T, N> const>(x);

    Intrepid2::Vector<T, N> const &
    xval = ex.getVector();

    auto const
    dim = xval.get_dimension();

    for (auto const i = 0; i < dim; ++i) {
      vector_(i) += xval(i);
    }
  }

  void
  axpy(T const alpha, Vector<T> const & x)
  {
    assert(dimension() == x.dimension());

    MiniTensorVector<T, N> const &
    ex = Teuchos::dyn_cast<MiniTensorVector<T, N> const>(x);

    Intrepid2::Vector<T, N> const &
    xval = ex.getVector();

    auto const
    dim = xval.get_dimension();

    assert(vector_.get_dimension() == dim);

    for (auto const i = 0; i < dim); ++i) {
      vector_(i) += alpha * xval(i);
    }
  }

  void
  scale(T const alpha)
  {
    auto const
    dim = vector_.get_dimension();

    for (auto const i = 0; i < dim); ++i) {
      vector_(i) *= alpha;
    }
  }

  virtual
  T
  dot(Vector<T> const & x) const
  {
    MiniTensorVector<T, N> const &
    ex = Teuchos::dyn_cast<MiniTensorVector<T, N> const>(x);

    Intrepid2::Vector<T, N> const &
    xval = ex.getVector();

    return Intrepid2::dot(vector_, xval);
  }

  T
  norm() const
  {
    return Intrepid2::norm(vector_);
  }

  virtual
  Teuchos::RCP<Vector<T>>
  clone() const
  {
    auto const
    dim = vector_.get_dimension();

    auto
    p_mt_vector = Teuchos::rcp(new Intrepid2::Vector<T, N>(dim));

    Teuchos::RCP<MiniTensorVector>
    e = Teuchos::rcp(new MiniTensorVector(*p_mt_vector));

    return e;
  }

  Intrepid2::Vector<T, N>
  getVector() const {
    return vector_;
  }

  Intrepid2::Vector<T, N>
  getVector()
  {
    return vector_;
  }

  Teuchos::RCP<Vector<T>>
  basis(int const i) const
  {
    auto const
    dim = vector_.get_dimension();

    auto
    p_mt_vector = Teuchos::rcp(new Intrepid2::Vector<T, N>(dim));

    p_mt_vector->fill(Intrepid2::ZEROS);

    (*p_mt_vector)(i) = 1.0;

    Teuchos::RCP<MiniTensorVector>
    e = Teuchos::rcp(new MiniTensorVector(*p_mt_vector));

    return e;
  }

  int
  dimension() const
  {
    return static_cast<int>(vector_->get_dimension());
  }

  void
  applyUnary(Elementwise::UnaryFunction<T> const & f)
  {
    auto const
    dim  = vector_->get_dimension();

    for(auto i{0}; i < dim; ++i) {
      vector_(i) = f.apply(vector_(i));
    }
  }

  void
  applyBinary(Elementwise::BinaryFunction<T> const & f, Vector<T> const & x )
  {
    MiniTensorVector<T, N> const &
    ex = Teuchos::dyn_cast<MiniTensorVector<T, N> const>(x);

    Intrepid2::Vector<T, N> const &
    xval = ex.getVector();

    auto const
    dim  = vector_->get_dimension();

    for(auto i{0}; i < dim; ++i) {
      vector_(i) = f.apply(vector_(i), xval(i));
    }
  }

  T
  reduce(Elementwise::ReductionOp<T> & r ) const
  {
    T
    result = r.initialValue();

    auto const
    dim  = vector_->get_dimension();

    for(auto i{0}; i < dim; ++i) {
      r.reduce(vector_(i), result);
    }

    return result;
  }

}; // class MiniTensorVector
}

#include "Intrepid2_MiniTensor_ROL_Vector.t.h"

#endif // Intrepid2_MiniTensor_ROL_Vector_h
#endif // ENABLE_ROL
