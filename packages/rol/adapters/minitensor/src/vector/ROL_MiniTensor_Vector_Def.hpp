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

namespace ROL {

//
// Covert from ROL to MiniTensor
//
template<typename T, minitensor::Index N>
minitensor::Vector<T, N>
MTfromROL(Vector<T> const & x)
{
  MiniTensorVector<T, N> const &
  xe = dynamic_cast<MiniTensorVector<T, N> const&>(x);

  minitensor::Vector<T, N> const &
  xval = xe.getVector();

  return xval;
}

//
// Convert from MiniTensor to ROL
//
template<typename T, minitensor::Index N>
void
MTtoROL(minitensor::Vector<T, N> const & xval, Vector<T> & x)
{
  MiniTensorVector<T, N> &
  xe = dynamic_cast<MiniTensorVector<T, N>&>(x);

  xe.set(xval);
}

//
//
//
template<typename T, minitensor::Index N>
MiniTensorVector<T, N>::
MiniTensorVector(minitensor::Vector<T, N> const & v) : vector_(v)
{
  return;
}

//
//
//
template<typename T, minitensor::Index N>
MiniTensorVector<T, N>::
~MiniTensorVector()
{
  return;
}

//
//
//
template<typename T, minitensor::Index N>
void
MiniTensorVector<T, N>::
set(Vector<T> const & x)
{
  vector_ = MTfromROL<T, N>(x);
}

//
//
//
template<typename T, minitensor::Index N>
void
MiniTensorVector<T, N>::
set(minitensor::Vector<T, N> const & x)
{
  vector_ = x;
}

//
//
//
template<typename T, minitensor::Index N>
void
MiniTensorVector<T, N>::
plus(Vector<T> const & x)
{
  minitensor::Vector<T, N> const
  xval = MTfromROL<T, N>(x);

  auto const
  dim = xval.get_dimension();

  assert(vector_.get_dimension() == dim);

  for (minitensor::Index i{0}; i < dim; ++i) {
    vector_(i) += xval(i);
  }
}

//
//
//
template<typename T, minitensor::Index N>
void
MiniTensorVector<T, N>::
axpy(T const alpha, Vector<T> const & x)
{
  minitensor::Vector<T, N> const
  xval = MTfromROL<T, N>(x);

  auto const
  dim = xval.get_dimension();

  assert(vector_.get_dimension() == dim);

  for (minitensor::Index i{0}; i < dim; ++i) {
    vector_(i) += alpha * xval(i);
  }
}

//
//
//
template<typename T, minitensor::Index N>
void
MiniTensorVector<T, N>::
scale(T const alpha)
{
  auto const
  dim = vector_.get_dimension();

  for (minitensor::Index i{0}; i < dim; ++i) {
    vector_(i) *= alpha;
  }
}

//
//
//
template<typename T, minitensor::Index N>
T
MiniTensorVector<T, N>::
dot(Vector<T> const & x) const
{
  minitensor::Vector<T, N> const
  xval = MTfromROL<T, N>(x);

  return minitensor::dot(vector_, xval);
}

//
//
//
template<typename T, minitensor::Index N>
T
MiniTensorVector<T, N>::
norm() const
{
  return minitensor::norm(vector_);
}

//
//
//
template<typename T, minitensor::Index N>
ROL::Ptr<Vector<T>>
MiniTensorVector<T, N>::
clone() const
{
  auto const
  dim = vector_.get_dimension();

  minitensor::Vector<T, N>
  val(dim, minitensor::Filler::ZEROS);

  ROL::Ptr<MiniTensorVector>
  e = ROL::makePtr<MiniTensorVector>(val);

  return e;
}

//
//
//
template<typename T, minitensor::Index N>
minitensor::Vector<T, N>
MiniTensorVector<T, N>::
getVector() const
{
  return vector_;
}

//
//
//
template<typename T, minitensor::Index N>
minitensor::Vector<T, N>
MiniTensorVector<T, N>::
getVector()
{
  return vector_;
}

//
//
//
template<typename T, minitensor::Index N>
ROL::Ptr<Vector<T>>
MiniTensorVector<T, N>::
basis(int const i) const
{
  auto const
  dim = vector_.get_dimension();

  minitensor::Vector<T, N>
  val(dim, minitensor::Filler::ZEROS);

  val(i) = 1.0;

  ROL::Ptr<MiniTensorVector>
  e = ROL::makePtr<MiniTensorVector>(val);

  return e;
}

//
//
//
template<typename T, minitensor::Index N>
int
MiniTensorVector<T, N>::
dimension() const
{
  return static_cast<int>(vector_.get_dimension());
}

//
//
//
template<typename T, minitensor::Index N>
void
MiniTensorVector<T, N>::
applyUnary(Elementwise::UnaryFunction<T> const & f)
{
  auto const
  dim = vector_.get_dimension();

  for(minitensor::Index i{0}; i < dim; ++i) {
    vector_(i) = f.apply(vector_(i));
  }
}

//
//
//
template<typename T, minitensor::Index N>
void
MiniTensorVector<T, N>::
applyBinary(Elementwise::BinaryFunction<T> const & f, Vector<T> const & x)
{
  minitensor::Vector<T, N> const
  xval = MTfromROL<T, N>(x);

  auto const
  dim  = vector_.get_dimension();

  for(minitensor::Index i{0}; i < dim; ++i) {
    vector_(i) = f.apply(vector_(i), xval(i));
  }
}

//
//
//
template<typename T, minitensor::Index N>
T
MiniTensorVector<T, N>::
reduce(Elementwise::ReductionOp<T> const & r) const
{
  T
  result = r.initialValue();

  auto const
  dim = vector_.get_dimension();

  for(minitensor::Index i{0}; i < dim; ++i) {
    r.reduce(vector_(i), result);
  }

  return result;
}

} // namespace ROL
