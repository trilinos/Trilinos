// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

namespace ROL {

using Index = minitensor::Index;

//
//
//
template<typename MSEC, typename S, Index M, Index N>
MiniTensor_EqualityConstraint<MSEC, S, M, N>::
MiniTensor_EqualityConstraint(MSEC & msec) : minisolver_ec_(msec)
{
  return;
}

//
//
//
template<typename MSEC, typename S, Index M, Index N>
MiniTensor_EqualityConstraint<MSEC, S, M, N>::
~MiniTensor_EqualityConstraint()
{
  return;
}

//
//
//
template<typename MSEC, typename S, Index M, Index N>
void
MiniTensor_EqualityConstraint<MSEC, S, M, N>::
value(Vector<S> & c, Vector<S> const & x, S &)
{
  minitensor::Vector<S, N> const
  xval = MTfromROL<S, N>(x);

  minitensor::Vector<S, M> const
  cval = minisolver_ec_.value(xval);

  MTtoROL<S, M>(cval, c);
}

//
//
//
template<typename MSEC, typename S, Index M, Index N>
void
MiniTensor_EqualityConstraint<MSEC, S, M, N>::
applyJacobian(Vector<S> & jv, Vector<S> const & v, Vector<S> const & x, S &)
{
  minitensor::Vector<S, N> const
  xval = MTfromROL<S, N>(x);

  minitensor::Vector<S, N> const
  vval = MTfromROL<S, N>(v);

  minitensor::Matrix<S, M, N> const
  J = minisolver_ec_.gradient(xval);

  minitensor::Vector<S, M> const
  jvval = minitensor::dot(J, vval);

  MTtoROL<S, M>(jvval, jv);
}

//
//
//
template<typename MSEC, typename S, Index M, Index N>
void
MiniTensor_EqualityConstraint<MSEC, S, M, N>::
applyAdjointJacobian(Vector<S> & ajv, Vector<S> const & v, Vector<S> const & x, S &)
{
  minitensor::Vector<S, N> const
  xval = MTfromROL<S, N>(x);

  minitensor::Vector<S, M> const
  vval = MTfromROL<S, M>(v);

  minitensor::Matrix<S, M, N> const
  J = minisolver_ec_.gradient(xval);

  minitensor::Vector<S, N> const
  ajvval = minitensor::dot(vval, J);

  MTtoROL<S, N>(ajvval, ajv);
}

} // namespace ROL
