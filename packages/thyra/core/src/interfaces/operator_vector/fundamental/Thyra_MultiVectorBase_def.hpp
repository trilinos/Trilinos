// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_MULTI_VECTOR_BASE_HPP
#define THYRA_MULTI_VECTOR_BASE_HPP

#include "Thyra_MultiVectorBase_decl.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"


namespace Thyra {


// Provide access to the columns as VectorBase objects


template<class Scalar>
RCP<const VectorBase<Scalar> >
MultiVectorBase<Scalar>::colImpl(Ordinal j) const
{
  return const_cast<MultiVectorBase*>(this)->nonconstColImpl(j);
}


// Overridden methods from LinearOpBase


template<class Scalar>
RCP<const LinearOpBase<Scalar> >
MultiVectorBase<Scalar>::clone() const
{
  return this->clone_mv();
}

#ifndef THYRA_HIDE_DEPRECATED_CODE
//
// Deprecated public function
//


template<class Scalar>
void MultiVectorBase<Scalar>::applyOp(
  const RTOpPack::RTOpT<Scalar> &primary_op,
  const int num_multi_vecs,
  const MultiVectorBase<Scalar>*const multi_vecs_in[],
  const int num_targ_multi_vecs,
  MultiVectorBase<Scalar>*const targ_multi_vecs_inout[],
  RTOpPack::ReductTarget*const reduct_objs_inout[],
  const Ordinal primary_global_offset
  ) const
{
  Array<Ptr<const MultiVectorBase<Scalar> > > multi_vecs;
  for (int k = 0; k < num_multi_vecs; ++k)
    multi_vecs.push_back(Teuchos::ptr(multi_vecs_in[k]));
  Array<Ptr<MultiVectorBase<Scalar> > > targ_multi_vecs;
  for (int k = 0; k < num_targ_multi_vecs; ++k)
    targ_multi_vecs.push_back(Teuchos::ptr(targ_multi_vecs_inout[k]));
  Array<Ptr<RTOpPack::ReductTarget> > reduct_objs;
  if (reduct_objs_inout) {
    const int secondary_sub_dim = ( num_multi_vecs
      ? multi_vecs[0]->domain() : targ_multi_vecs[0]->domain()
      )->dim();
    for (int k = 0; k < secondary_sub_dim; ++k)
      reduct_objs.push_back(Teuchos::ptr(reduct_objs_inout[k]));
  }
  mvMultiReductApplyOpImpl(
    primary_op,
    multi_vecs, targ_multi_vecs,
    reduct_objs, primary_global_offset);
}


template<class Scalar>
void MultiVectorBase<Scalar>::applyOp(
  const RTOpPack::RTOpT<Scalar> &primary_op,
  const RTOpPack::RTOpT<Scalar> &secondary_op,
  const int num_multi_vecs,
  const MultiVectorBase<Scalar>*const multi_vecs_in[],
  const int num_targ_multi_vecs,
  MultiVectorBase<Scalar>*const targ_multi_vecs_inout[],
  RTOpPack::ReductTarget* reduct_obj,
  const Ordinal primary_global_offset
  ) const
{
  Array<Ptr<const MultiVectorBase<Scalar> > > multi_vecs;
  for (int k = 0; k < num_multi_vecs; ++k)
    multi_vecs.push_back(Teuchos::ptr(multi_vecs_in[k]));
  Array<Ptr<MultiVectorBase<Scalar> > > targ_multi_vecs;
  for (int k = 0; k < num_targ_multi_vecs; ++k)
    targ_multi_vecs.push_back(Teuchos::ptr(targ_multi_vecs_inout[k]));
  mvSingleReductApplyOpImpl(
    primary_op, secondary_op,
    multi_vecs, targ_multi_vecs,
    Teuchos::ptr(reduct_obj),
    primary_global_offset);
}
#endif // THYRA_HIDE_DEPRECATED_CODE

} // end namespace Thyra


#endif // THYRA_MULTI_VECTOR_BASE_HPP
