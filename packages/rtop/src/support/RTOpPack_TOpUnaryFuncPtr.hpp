// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef RTOPPACK_UNARY_FUNC_PTR_HPP
#define RTOPPACK_UNARY_FUNC_PTR_HPP

#include "RTOpPack_TOpUnaryFuncPtrDecl.hpp"

namespace RTOpPack {

template<class Scalar>
TOpUnaryFuncPtr<Scalar>::TOpUnaryFuncPtr()
  :RTOpT<Scalar>("TOpUnaryFuncPtr")
{
  set_initialized();
}

template<class Scalar>
TOpUnaryFuncPtr<Scalar>::TOpUnaryFuncPtr(
  unary_func_ptr_t        unary_func_ptr
  ,const std::string      &op_name
  )
  :RTOpT<Scalar>("TOpUnaryFuncPtr")
{
  initialize(unary_func_ptr,op_name);
}

template<class Scalar>
void TOpUnaryFuncPtr<Scalar>::initialize(
  unary_func_ptr_t        unary_func_ptr
  ,const std::string      &op_name
  )
{
  TEST_FOR_EXCEPTION( unary_func_ptr==NULL, std::invalid_argument, "Error!" );
  unary_func_ptr_ = unary_func_ptr;
  op_name_ = op_name;
}

template<class Scalar>
void TOpUnaryFuncPtr<Scalar>::set_initialized(
  unary_func_ptr_t    *unary_func_ptr
  ,std::string        *op_name
  )
{
  if(unary_func_ptr) *unary_func_ptr = unary_func_ptr_;
  if(op_name) *op_name = op_name_;

  unary_func_ptr_ = NULL;
  op_name_ = "uninitialized()";
}

// Overridden from RTOpT

template<class Scalar>
const char* TOpUnaryFuncPtr<Scalar>::op_name() const
{
  return op_name_.c_str();
}

template<class Scalar>
void TOpUnaryFuncPtr<Scalar>::apply_op(
  const int   num_vecs,       const ConstSubVectorView<Scalar>         sub_vecs[]
  ,const int  num_targ_vecs,  const SubVectorView<Scalar>  targ_sub_vecs[]
  ,ReductTarget *reduct_obj
  ) const
{
  TEST_FOR_EXCEPTION( num_vecs != 1 || sub_vecs == NULL, std::invalid_argument, "Error!" );
  TEST_FOR_EXCEPTION( num_targ_vecs != 1 || targ_sub_vecs == NULL, std::invalid_argument, "Error!" );
  TEST_FOR_EXCEPTION( reduct_obj != NULL, std::invalid_argument, "Error!" );
  TEST_FOR_EXCEPTION( sub_vecs[0].stride() != 1, std::invalid_argument, "Error, can't handle non-unit strides here!" );
  TEST_FOR_EXCEPTION( targ_sub_vecs[0].stride() != 1, std::invalid_argument, "Error, can't handle non-unit strides here!" );
  TEST_FOR_EXCEPTION( sub_vecs[0].subDim() != targ_sub_vecs[0].subDim(), std::invalid_argument, "Error!" );
  TEST_FOR_EXCEPTION( sub_vecs[0].globalOffset() != targ_sub_vecs[0].globalOffset(), std::invalid_argument, "Error!" );

  unary_func_ptr_( sub_vecs[0].values(), sub_vecs[0].subDim(), targ_sub_vecs[0].values() );

}

} // end namespace RTOpPack

#endif // RTOPPACK_UNARY_FUNC_PTR_HPP
