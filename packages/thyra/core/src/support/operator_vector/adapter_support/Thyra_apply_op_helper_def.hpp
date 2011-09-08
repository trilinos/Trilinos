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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_APPLY_OP_HELPER_HPP
#define THYRA_APPLY_OP_HELPER_HPP

#include "Thyra_apply_op_helper_decl.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_AssertOp.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_as.hpp"


template<class Scalar>
void Thyra::apply_op_validate_input(
  const std::string &func_name,
  const VectorSpaceBase<Scalar> &space,
  const RTOpPack::RTOpT<Scalar> &op,
  const ArrayView<const Ptr<const VectorBase<Scalar> > > &vecs,
  const ArrayView<const Ptr<VectorBase<Scalar> > > &targ_vecs,
  const Ptr<RTOpPack::ReductTarget> &reduct_obj,
  const Ordinal global_offset_in
  )
{
  const int num_vecs = vecs.size();
  const int num_targ_vecs = targ_vecs.size();
  TEST_FOR_EXCEPTION(
    global_offset_in < 0, std::logic_error
    ,func_name << " : Error! global_offset_in = "
    <<global_offset_in<<" is not valid" );
  for (int k = 0; k < num_vecs; ++k)
    THYRA_ASSERT_VEC_SPACES(func_name,space,*vecs[k]->space());
  for (int k = 0; k < num_targ_vecs; ++k)
    THYRA_ASSERT_VEC_SPACES(func_name,space,*targ_vecs[k]->space());
}


template<class Scalar>
void Thyra::apply_op_validate_input(
  const std::string &func_name,
  const VectorSpaceBase<Scalar> &domain,
  const VectorSpaceBase<Scalar> &range,
  const RTOpPack::RTOpT<Scalar> &primary_op,
  const ArrayView<const Ptr<const MultiVectorBase<Scalar> > > &multi_vecs,
  const ArrayView<const Ptr<MultiVectorBase<Scalar> > > &targ_multi_vecs,
  const ArrayView<const Ptr<RTOpPack::ReductTarget> > &reduct_objs,
  const Ordinal primary_global_offset_in
  )
{
  using Teuchos::as;
  // Validate primary range arguments
  TEST_FOR_EXCEPTION(
    primary_global_offset_in < 0, std::logic_error
    ,func_name << " : Error! primary_global_offset_in = "
    <<primary_global_offset_in<<" is not valid" );
  // Validate secondary domain arguments
  // Validate spaces
  for (int k = 0; k < multi_vecs.size(); ++k) {
    THYRA_ASSERT_VEC_SPACES(func_name,domain,*multi_vecs[k]->domain());
    THYRA_ASSERT_VEC_SPACES(func_name,range,*multi_vecs[k]->range());
  }
  for (int k = 0; k < targ_multi_vecs.size(); ++k) {
    THYRA_ASSERT_VEC_SPACES(func_name,domain,*targ_multi_vecs[k]->domain());
    THYRA_ASSERT_VEC_SPACES(func_name,range,*targ_multi_vecs[k]->range());
  }
}


//
// Explicit instant macro
//

#define THYRA_APPLY_OP_HELPER_INSTANT(SCALAR) \
   \
  template void apply_op_validate_input( \
    const std::string &func_name, \
    const VectorSpaceBase<SCALAR > &space, \
    const RTOpPack::RTOpT<SCALAR > &op, \
    const ArrayView<const Ptr<const VectorBase<SCALAR > > > &vecs, \
    const ArrayView<const Ptr<VectorBase<SCALAR > > > &targ_vecs, \
    const Ptr<RTOpPack::ReductTarget> &reduct_obj, \
    const Ordinal global_offset_in \
    ); \
   \
  template void apply_op_validate_input( \
    const std::string &func_name, \
    const VectorSpaceBase<SCALAR > &domain, \
    const VectorSpaceBase<SCALAR > &range, \
    const RTOpPack::RTOpT<SCALAR > &primary_op, \
    const ArrayView<const Ptr<const MultiVectorBase<SCALAR > > > &multi_vecs, \
    const ArrayView<const Ptr<MultiVectorBase<SCALAR > > > &targ_multi_vecs, \
    const ArrayView<const Ptr<RTOpPack::ReductTarget> > &reduct_objs, \
    const Ordinal primary_global_offset_in \
    ); \



#endif // THYRA_APPLY_OP_HELPER_HPP
