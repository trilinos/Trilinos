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
