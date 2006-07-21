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

#ifndef RTOPPACK_RTOP_NEW_T_HPP
#define RTOPPACK_RTOP_NEW_T_HPP

#include "RTOpPack_RTOpTDecl.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace RTOpPack {

template<class Scalar>
RTOpT<Scalar>::RTOpT( const std::string &op_name_base )
  :op_name_(op_name_base + std::string(Teuchos::ScalarTraits<Scalar>::name()))
{}

// Reduction object functions

template<class Scalar>
void RTOpT<Scalar>::get_reduct_type_num_entries(
  int*   num_values
  ,int*  num_indexes
  ,int*  num_chars
  ) const
{
#ifdef RTOp_DEBUG
  TEST_FOR_EXCEPTION( !num_values, std::logic_error, "Error!" );
  TEST_FOR_EXCEPTION( !num_indexes, std::logic_error, "Error!"  );
  TEST_FOR_EXCEPTION( !num_chars, std::logic_error, "Error!"  );
#endif
  *num_values  = 0;
  *num_indexes = 0;
  *num_chars   = 0;
}

template<class Scalar>
Teuchos::RefCountPtr<ReductTarget>
RTOpT<Scalar>::reduct_obj_create() const
{
  return Teuchos::null;
}

template<class Scalar>
void RTOpT<Scalar>::reduce_reduct_objs(
  const ReductTarget& in_reduct_obj, ReductTarget* inout_reduct_obj
  ) const
{
  TEST_FOR_EXCEPTION(true,std::logic_error,"Error, should not call!");
}

template<class Scalar>
void RTOpT<Scalar>::reduct_obj_reinit( ReductTarget* reduct_obj ) const
{
  TEST_FOR_EXCEPTION(true,std::logic_error,"Error, should not call!");
}

template<class Scalar>
void RTOpT<Scalar>::extract_reduct_obj_state(
  const ReductTarget     &reduct_obj
  ,int                      num_values
  ,primitive_value_type     value_data[]
  ,int                      num_indexes
  ,index_type               index_data[]
  ,int                      num_chars
  ,char_type                char_data[]
  ) const
{
  TEST_FOR_EXCEPTION(true,std::logic_error,"Error, should not call!");
}

template<class Scalar>
void RTOpT<Scalar>::load_reduct_obj_state(
  int                            num_values
  ,const primitive_value_type    value_data[]
  ,int                           num_indexes
  ,const index_type              index_data[]
  ,int                           num_chars
  ,const char_type               char_data[]
  ,ReductTarget               *reduct_obj
  ) const
{
  TEST_FOR_EXCEPTION(true,std::logic_error,"Error, should not call!");
}

// Operator functions

template<class Scalar>
const char* RTOpT<Scalar>::op_name() const
{
  return op_name_.c_str();
}

template<class Scalar>
RTOpT<Scalar>& RTOpT<Scalar>::operator=(const RTOpT<Scalar>& op)
{
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  int num_values = 0, num_indexes = 0, num_chars = 0;
  op.get_op_type_num_entries( &num_values, &num_indexes, &num_chars );
  Workspace<primitive_value_type> value_data(wss,num_values,false);
  Workspace<index_type>           index_data(wss,num_indexes,false);
  Workspace<char_type>            char_data(wss,num_chars,false);
  op.extract_op_state(
    num_values,   num_values  ? &value_data[0] : NULL
    ,num_indexes, num_indexes ? &index_data[0] : NULL
    ,num_chars,   num_chars   ? &char_data[0]  : NULL
    );
  this->load_op_state(
    num_values,   num_values  ? &value_data[0] : NULL
    ,num_indexes, num_indexes ? &index_data[0] : NULL
    ,num_chars,   num_chars   ? &char_data[0]  : NULL
    );
  return *this;
}

template<class Scalar>
void RTOpT<Scalar>::get_op_type_num_entries(
  int*  num_values
  ,int* num_indexes
  ,int* num_chars
  ) const
{
#ifdef RTOp_DEBUG
  TEST_FOR_EXCEPTION( !num_values, std::logic_error, "Error!" );
  TEST_FOR_EXCEPTION( !num_indexes, std::logic_error, "Error!"  );
  TEST_FOR_EXCEPTION( !num_chars, std::logic_error, "Error!"  );
#endif
  *num_values  = 0;
  *num_indexes = 0;
  *num_chars   = 0;
}

template<class Scalar>
void RTOpT<Scalar>::extract_op_state(
  int                             num_values
  ,primitive_value_type           value_data[]
  ,int                            num_indexes
  ,index_type                     index_data[]
  ,int                            num_chars
  ,char_type                      char_data[]
  ) const
{
  TEST_FOR_EXCEPTION(true,std::logic_error,"Error, should not call!");
}

template<class Scalar>
void RTOpT<Scalar>::load_op_state(
  int                           num_values
  ,const primitive_value_type   value_data[]
  ,int                          num_indexes
  ,const index_type             index_data[]
  ,int                          num_chars
  ,const char_type              char_data[]
  )
{
  TEST_FOR_EXCEPTION(true,std::logic_error,"Error, should not call!");
}

template<class Scalar>
bool RTOpT<Scalar>::coord_invariant() const
{
  return true;
}

} // end namespace RTOpPack

#endif // RTOPPACK_RTOP_NEW_T_HPP
