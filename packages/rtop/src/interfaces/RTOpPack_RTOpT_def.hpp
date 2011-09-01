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

#ifndef RTOPPACK_RTOP_NEW_T_HPP
#define RTOPPACK_RTOP_NEW_T_HPP

#include "RTOpPack_RTOpT_decl.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_ScalarTraits.hpp"


namespace RTOpPack {


// Protected functions to be overridden by subclasses


template<class Scalar>
void RTOpT<Scalar>::get_reduct_type_num_entries_impl(
  const Ptr<int> &num_values,
  const Ptr<int> &num_indexes,
  const Ptr<int> &num_chars
  ) const
{
  *num_values  = 0;
  *num_indexes = 0;
  *num_chars   = 0;
}


template<class Scalar>
Teuchos::RCP<ReductTarget>
RTOpT<Scalar>::reduct_obj_create_impl() const
{
  return Teuchos::null;
}


template<class Scalar>
void RTOpT<Scalar>::reduce_reduct_objs_impl(
  const ReductTarget& in_reduct_obj, const Ptr<ReductTarget>& inout_reduct_obj
  ) const
{
  throwNoReductError();
}


template<class Scalar>
void RTOpT<Scalar>::reduct_obj_reinit_impl(
  const Ptr<ReductTarget> &reduct_obj ) const
{
  throwNoReductError();
}


template<class Scalar>
void RTOpT<Scalar>::extract_reduct_obj_state_impl(
  const ReductTarget &reduct_obj,
  const ArrayView<primitive_value_type> &value_data,
  const ArrayView<index_type> &index_data,
  const ArrayView<char_type> &char_data
  ) const
{
  throwNoReductError();
}


template<class Scalar>
void RTOpT<Scalar>::load_reduct_obj_state_impl(
  const ArrayView<const primitive_value_type> &value_data,
  const ArrayView<const index_type> &index_data,
  const ArrayView<const char_type> &char_data,
  const Ptr<ReductTarget> &reduct_obj
  ) const
{
  throwNoReductError();
}


template<class Scalar>
std::string RTOpT<Scalar>::op_name_impl() const
{
  return op_name_;
}


template<class Scalar>
bool RTOpT<Scalar>::coord_invariant_impl() const
{
  return true;
}


template<class Scalar>
Range1D RTOpT<Scalar>::range_impl() const
{
  return Range1D();
}


// Nonvirtual protected functions


template<class Scalar>
RTOpT<Scalar>::RTOpT( const std::string &op_name_base )
{
  setOpNameBase(op_name_base);
}


template<class Scalar>
void RTOpT<Scalar>::setOpNameBase( const std::string &op_name_base )
{
  op_name_ = op_name_base+"<"+ScalarTraits<Scalar>::name()+">";
}


// private


template<class Scalar>
void RTOpT<Scalar>::throwNoReductError() const
{
  TEST_FOR_EXCEPTION(true, std::logic_error,
    "Error, no reduction is defined for concrete reduction op \'"
    << this->description() << "\'!" );
}



} // end namespace RTOpPack


#endif // RTOPPACK_RTOP_NEW_T_HPP
