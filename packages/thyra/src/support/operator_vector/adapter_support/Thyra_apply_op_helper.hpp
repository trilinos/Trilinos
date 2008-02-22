// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
  const Index first_ele_offset_in,
  const Index sub_dim_in,
  const Index global_offset_in
  )
{
  const int num_vecs = vecs.size();
  const int num_targ_vecs = targ_vecs.size();
  const Index
    dim = space.dim();
  TEST_FOR_EXCEPTION(
    global_offset_in < 0, std::logic_error
    ,func_name << " : Error! global_offset_in = "
    <<global_offset_in<<" is not valid" );
  TEST_FOR_EXCEPTION(
    first_ele_offset_in+1 > dim, std::logic_error
    ,func_name << " : Error! first_ele_offset_in = "
    <<first_ele_offset_in<<" is not compatible with space.dim() = " << dim );
  TEST_FOR_EXCEPTION(
    (sub_dim_in > 0 && sub_dim_in > dim-first_ele_offset_in), std::logic_error
    ,func_name << " : Error! first_ele_offset_in = "
    <<first_ele_offset_in<<" and sub_dim_in = "<<sub_dim_in
    <<" is not compatible with space.dim() = " << dim );
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
  const Index primary_first_ele_offset_in,
  const Index primary_sub_dim_in,
  const Index primary_global_offset_in,
  const Index secondary_first_ele_offset_in,
  const Index secondary_sub_dim_in
  )
{
  using Teuchos::as;
  // Validate primary range arguments
  const Index
    range_dim = range.dim();
  TEST_FOR_EXCEPTION(
    primary_global_offset_in < 0, std::logic_error
    ,func_name << " : Error! primary_global_offset_in = "
    <<primary_global_offset_in<<" is not valid" );
  TEST_FOR_EXCEPTION(
    primary_first_ele_offset_in < 0 || range_dim < primary_first_ele_offset_in+1, std::logic_error
    ,func_name << " : Error! primary_first_ele_offset_in = "
    <<primary_first_ele_offset_in<<" is not compatible with range.dim() = " << range_dim );
  TEST_FOR_EXCEPTION(
    (primary_sub_dim_in > 0 && primary_sub_dim_in > range_dim-primary_first_ele_offset_in)
    , std::logic_error
    ,func_name << " : Error! primary_first_ele_offset_in = "
    <<primary_first_ele_offset_in<<" and primary_sub_dim_in = "<<primary_sub_dim_in
    <<" are not compatible with range.dim() = " << range_dim );
  // Validate secondary domain arguments
  const Index
    domain_dim = domain.dim();
  TEST_FOR_EXCEPTION(
    secondary_first_ele_offset_in < 0 || domain_dim < secondary_first_ele_offset_in+1, std::logic_error
    ,func_name << " : Error! secondary_first_ele_offset_in = "
    <<secondary_first_ele_offset_in<<" is not compatible with domain.dim() = " << domain_dim );
  TEST_FOR_EXCEPTION(
    (secondary_sub_dim_in > 0 && secondary_sub_dim_in > domain_dim-secondary_first_ele_offset_in)
    , std::logic_error
    ,func_name << " : Error! secondary_first_ele_offset_in = "
    <<secondary_first_ele_offset_in<<" and secondary_sub_dim_in = "<<secondary_sub_dim_in
    <<" are not compatible with domain.dim() = " << domain_dim );
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


#endif // THYRA_APPLY_OP_HELPER_HPP
