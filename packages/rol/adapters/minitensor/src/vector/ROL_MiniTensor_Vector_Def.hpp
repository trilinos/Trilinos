// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
