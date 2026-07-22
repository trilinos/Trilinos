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
template<typename MSFN, typename S, Index M>
MiniTensor_Objective<MSFN, S, M>::
MiniTensor_Objective(MSFN & msfn) : minisolver_fn_(msfn)
{
  return;
}

//
//
//
template<typename MSFN, typename S, Index M>
MiniTensor_Objective<MSFN, S, M>::
~MiniTensor_Objective()
{
  return;
}

//
//
//
template<typename MSFN, typename S, Index M>
S
MiniTensor_Objective<MSFN, S, M>::
value(Vector<S> const & x, S &)
{
  minitensor::Vector<S, M> const
  xval = MTfromROL<S, M>(x);

  return minisolver_fn_.value(xval);
}

//
//
//
template<typename MSFN, typename S, Index M>
void
MiniTensor_Objective<MSFN, S, M>::
gradient(Vector<S> & g, Vector<S> const & x, S &)
{
  minitensor::Vector<S, M> const
  xval = MTfromROL<S, M>(x);

  minitensor::Vector<S, M> const
  gval = minisolver_fn_.gradient(xval);

  MTtoROL<S, M>(gval, g);
}

//
//
//
template<typename MSFN, typename S, Index M>
void
MiniTensor_Objective<MSFN, S, M>::
hessVec(Vector<S> & hv, Vector<S> const & v, Vector<S> const & x, S &)
{
  minitensor::Vector<S, M> const
  xval = MTfromROL<S, M>(x);

  minitensor::Vector<S, M> const
  vval = MTfromROL<S, M>(v);

  minitensor::Tensor<S, M> const
  H = minisolver_fn_.hessian(xval);

  minitensor::Vector<S, M> const
  hvval = H * vval;

  MTtoROL<S, M>(hvval, hv);
}

//
//
//
template<typename MSFN, typename S, Index M>
void
MiniTensor_Objective<MSFN, S, M>::
invHessVec(Vector<S> & hv, Vector<S> const & v, Vector<S> const & x, S &)
{
  minitensor::Vector<S, M> const
  xval = MTfromROL<S, M>(x);

  minitensor::Vector<S, M> const
  vval = MTfromROL<S, M>(v);

  minitensor::Tensor<S, M> const
  H = minisolver_fn_.hessian(xval);

  minitensor::Vector<S, M> const
  hvval = minitensor::solve(H, vval);

  MTtoROL<S, M>(hvval, hv);
}

} // namespace ROL
