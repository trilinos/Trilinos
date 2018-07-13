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
