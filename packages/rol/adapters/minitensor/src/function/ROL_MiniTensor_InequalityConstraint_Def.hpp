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

using Index = minitensor::Index;

//
//
//
template<typename MSIC, typename S, Index M, Index N>
MiniTensor_InequalityConstraint<MSIC, S, M, N>::
MiniTensor_InequalityConstraint(MSIC & msec) : minisolver_ic_(msec)
{
  return;
}

//
//
//
template<typename MSIC, typename S, Index M, Index N>
MiniTensor_InequalityConstraint<MSIC, S, M, N>::
~MiniTensor_InequalityConstraint()
{
  return;
}

//
//
//
template<typename MSIC, typename S, Index M, Index N>
void
MiniTensor_InequalityConstraint<MSIC, S, M, N>::
value(Vector<S> & c, Vector<S> const & x, S &)
{
  minitensor::Vector<S, N> const
  xval = MTfromROL<S, N>(x);

  minitensor::Vector<S, M> const
  cval = minisolver_ic_.value(xval);

  MTtoROL<S, M>(cval, c);
}

//
//
//
template<typename MSIC, typename S, Index M, Index N>
void
MiniTensor_InequalityConstraint<MSIC, S, M, N>::
applyJacobian(Vector<S> & jv, Vector<S> const & v, Vector<S> const & x, S &)
{
  minitensor::Vector<S, N> const
  xval = MTfromROL<S, N>(x);

  minitensor::Vector<S, N> const
  vval = MTfromROL<S, N>(v);

  minitensor::Matrix<S, M, N> const
  J = minisolver_ic_.gradient(xval);

  minitensor::Vector<S, M> const
  jvval = minitensor::dot(J, vval);

  MTtoROL<S, M>(jvval, jv);
}

//
//
//
template<typename MSIC, typename S, Index M, Index N>
void
MiniTensor_InequalityConstraint<MSIC, S, M, N>::
applyAdjointJacobian(Vector<S> & ajv, Vector<S> const & v, Vector<S> const & x, S &)
{
  minitensor::Vector<S, N> const
  xval = MTfromROL<S, N>(x);

  minitensor::Vector<S, M> const
  vval = MTfromROL<S, M>(v);

  minitensor::Matrix<S, M, N> const
  J = minisolver_ic_.gradient(xval);

  minitensor::Vector<S, N> const
  ajvval = minitensor::dot(vval, J);

  MTtoROL<S, N>(ajvval, ajv);
}

} // namespace ROL
