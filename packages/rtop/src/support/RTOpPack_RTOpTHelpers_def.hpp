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

#ifndef RTOPPACK_RTOP_T_HELPERS_DEF_HPP
#define RTOPPACK_RTOP_T_HELPERS_DEF_HPP


#include "RTOpPack_RTOpTHelpers_decl.hpp"


namespace RTOpPack {


//
// DefaultReductTarget
//


template<class ConcreteReductObj>
std::string DefaultReductTarget<ConcreteReductObj>::description() const
{
  std::ostringstream oss;
  oss
    << "RTOpPack::DefaultReductTarget<"
    <<TypeNameTraits<ConcreteReductObj>::name()<<">"
    << "{concreteReductObj="<<concreteReductObj_<<"}";
  return oss.str();
}


//
// ROpScalarReductionBase
//


template<class Scalar, class ConcreteReductObj>
void ROpScalarReductionBase<Scalar, ConcreteReductObj>::reduce_reduct_objs(
  const ReductTarget& in_reduct_obj, ReductTarget* inout_reduct_obj
  ) const
{
  TEST_FOR_EXCEPT(true); // ToDo: This function should not even exist!
}


template<class Scalar, class ConcreteReductObj>
void ROpScalarReductionBase<Scalar, ConcreteReductObj>::load_reduct_obj_state(
  int num_values,
  const primitive_value_type value_data[],
  int num_indexes,
  const index_type index_data[],
  int num_chars,
  const char_type char_data[],
  ReductTarget *reduct_obj
  ) const
{

  using Teuchos::arrayView;
  typedef ScalarTraits<Scalar> ST;
  typedef PrimitiveTypeTraits<Scalar, ConcreteReductObj> PTT;

#ifdef RTOPPACK_RTOPT_HELPER_DUMP_OUTPUT
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab tab(out);
  if(rtop_helpers_dump_all) {
    *out << "\nEntering ROpScalarReductionBase::load_reduct_obj_state(...) ...\n"
         << "\nOn input:\n";
    Teuchos::OSTab tab(out);
    *out << "op = " << this->description() << "\n";
    *out << "num_values = " << num_values << "\n";
    if(num_values) {
      *out <<"value_data[] = { ";
      for( int i = 0; i < num_values-1; ++i )
        *out << value_data[i] << ", ";
      *out << value_data[num_values-1] << " }\n";
    }
    *out << "num_indexes = " << num_indexes << "\n";
    if(num_indexes) {
      *out <<"index_data[] = { ";
      for( int i = 0; i < num_indexes-1; ++i )
        *out << index_data[i] << ", ";
      *out << index_data[num_indexes-1] << " }\n";
    }
    *out << "num_chars = " << num_chars << "\n";
    if(num_chars) {
      *out <<"char_data[] = { ";
      for( int i = 0; i < num_chars-1; ++i )
        *out << char_data[i] << ", ";
      *out << char_data[num_chars-1] << " }\n";
    }
  }
#endif // RTOPPACK_RTOPT_HELPER_DUMP_OUTPUT

  ConcreteReductObj concrete_reduct_obj;

  PTT::loadPrimitiveObjs(
    arrayView(value_data, num_values), arrayView(index_data, num_indexes),
    arrayView(char_data, num_chars),
    Teuchos::outArg(concrete_reduct_obj) );

  this->setRawVal( concrete_reduct_obj, reduct_obj );

#ifdef RTOPPACK_RTOPT_HELPER_DUMP_OUTPUT
  if(rtop_helpers_dump_all) {
    *out << "\nOn output:\n";
    Teuchos::OSTab tab(out);
    *out << "val = " << val << "\n";
    *out << "reduct_op = " << Teuchos::describe(*reduct_obj, Teuchos::VERB_EXTREME);
  }
#endif // RTOPPACK_RTOPT_HELPER_DUMP_OUTPUT

}


//
// ROpScalarReductionWithOpBase
//


} // namespace RTOpPack


template<class Scalar>
void RTOpPack::validate_apply_op(
  const RTOpT<Scalar> &op,
  const int allowed_num_sub_vecs,
  const int allowed_num_targ_sub_vecs,
  const bool expect_reduct_obj,
  const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
  const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
  const Ptr<const ReductTarget> &reduct_obj
  )
{

  using Teuchos::as;

  const int num_sub_vecs = sub_vecs.size();
  const int num_targ_sub_vecs = targ_sub_vecs.size();

  const std::string op_name_str = "op.name() = " + op.op_name();

  if (allowed_num_sub_vecs >= 0) {
    TEST_FOR_EXCEPTION( allowed_num_sub_vecs != as<int>(sub_vecs.size()),
      InvalidNumVecs,
      op_name_str<<": Error, sub_vecs.size()="<<sub_vecs.size()
      <<" != allowed_num_sub_vecs="<<allowed_num_sub_vecs<<"!" );
  }

  if (allowed_num_targ_sub_vecs >= 0) {
    TEST_FOR_EXCEPTION( allowed_num_targ_sub_vecs != as<int>(targ_sub_vecs.size()),
      InvalidNumTargVecs,
      op_name_str<<": Error, targ_sub_vecs.size()="<<targ_sub_vecs.size()
      <<" != allowed_num_targ_sub_vecs="<<allowed_num_targ_sub_vecs<<"!" );
  }
  
  TEST_FOR_EXCEPTION( sub_vecs.size() == 0 && targ_sub_vecs.size() == 0,
    InvalidNumVecs,
    op_name_str<<": Error, apply_op(...) must be passed some vectors!"
    );
  
  const index_type subDim =
    (sub_vecs.size() ? sub_vecs[0].subDim() : targ_sub_vecs[0].subDim());
  
  const index_type globalOffset =
    (sub_vecs.size() ? sub_vecs[0].globalOffset() : targ_sub_vecs[0].globalOffset());
  
  for (int k = 0; k < num_sub_vecs; ++k ) {
    TEST_FOR_EXCEPTION(
      sub_vecs[k].subDim() != subDim || sub_vecs[k].globalOffset() != globalOffset,
      IncompatibleVecs,
      op_name_str<<": Error, sub_vec["<<k<<"] "
      "(subDim="<<sub_vecs[k].subDim()<<",globalOffset="<<sub_vecs[k].globalOffset()<<")"
      " is not compatible with (subDim="<<subDim<<",globalOffset="<<globalOffset<<")!"
      );
  }
  
  for (int k = 0; k < num_targ_sub_vecs; ++k ) {
    TEST_FOR_EXCEPTION(
      targ_sub_vecs[k].subDim() != subDim || targ_sub_vecs[k].globalOffset() != globalOffset,
      IncompatibleVecs,
      op_name_str<<": Error, sub_vec["<<k<<"] "
      "(subDim="<<targ_sub_vecs[k].subDim()<<",globalOffset="<<targ_sub_vecs[k].globalOffset()<<")"
      " is not compatible with (subDim="<<subDim<<",globalOffset="<<globalOffset<<")!"
      );
  }
  
  if (expect_reduct_obj) {
    
    TEST_FOR_EXCEPTION( is_null(reduct_obj),
      IncompatibleReductObj,
      op_name_str<<": Error, expected a reduction target object!"
      );
    
    const RCP<ReductTarget> dummy_reduct_obj = op.reduct_obj_create();
    
    const std::type_info &reduct_obj_type = typeid(*reduct_obj);
    const std::type_info &dummy_reduct_obj_type = typeid(*dummy_reduct_obj);
    
    TEST_FOR_EXCEPTION( reduct_obj_type != dummy_reduct_obj_type,
      IncompatibleReductObj,
      op_name_str<<": Error, the type of the input reduct_obj = "
      <<Teuchos::demangleName(reduct_obj_type.name())<<" is not equal to"
      " the expected type "<<Teuchos::demangleName(dummy_reduct_obj_type.name())
      <<"!"
      );

  }
      
}



//
// Explicit Instantiation Macro
//


#define RTOPPACK_RTOPT_HELPERS_DEFAULTREDUCTTARGET_INSTANT(SCALAR) \
  \
  template std::string DefaultReductTarget<SCALAR >::description() const;


#define RTOPPACK_RTOPT_HELPERS_ROPSCALARREDUCTIONBASE_INSTANT( \
  SCALAR, REDUCTSCALAR \
  ) \
  \
  template void ROpScalarReductionBase<SCALAR, REDUCTSCALAR >::reduce_reduct_objs( \
    const ReductTarget& in_reduct_obj, ReductTarget* inout_reduct_obj \
    ) const; \
  \
  template void ROpScalarReductionBase<SCALAR, REDUCTSCALAR >::load_reduct_obj_state( \
    int num_values, \
    const primitive_value_type value_data[], \
    int num_indexes, \
    const index_type index_data[], \
    int num_chars, \
    const char_type char_data[], \
    ReductTarget *reduct_obj \
    ) const;


#define RTOPPACK_RTOPT_HELPERS_ROPSCALARREDUCTIONBASE_INDEX_INSTANT( \
  SCALAR \
  ) \
  \
  RTOPPACK_RTOPT_HELPERS_ROPSCALARREDUCTIONBASE_INSTANT(SCALAR, index_type)


#define RTOPPACK_RTOPT_HELPERS_ROPSCALARREDUCTIONBASE_SCALARINDEX_INSTANT( \
  SCALAR \
  ) \
  \
  RTOPPACK_RTOPT_HELPERS_ROPSCALARREDUCTIONBASE_INSTANT(SCALAR, ScalarIndex<SCALAR >)


#define RTOPPACK_RTOPT_HELPERS_ROPSCALARREDUCTIONBASE_SUBVECTORVIEW_INSTANT( \
  SCALAR \
  ) \
  \
  RTOPPACK_RTOPT_HELPERS_ROPSCALARREDUCTIONBASE_INSTANT(SCALAR, SubVectorView<SCALAR >)


#define RTOPPACK_RTOPT_HELPERS_INSTANT_SCALAR(SCALAR) \
  \
  RTOPPACK_RTOPT_HELPERS_DEFAULTREDUCTTARGET_INSTANT(SCALAR) \
  \
  RTOPPACK_RTOPT_HELPERS_DEFAULTREDUCTTARGET_INSTANT(ScalarIndex<SCALAR >) \
  \
  RTOPPACK_RTOPT_HELPERS_DEFAULTREDUCTTARGET_INSTANT(SubVectorView<SCALAR >) \
  \
  RTOPPACK_RTOPT_HELPERS_ROPSCALARREDUCTIONBASE_INSTANT(SCALAR, SCALAR) \
  \
  template void validate_apply_op<SCALAR >( \
    const RTOpT<SCALAR > &op, \
    const int allowed_num_sub_vecs, \
    const int allowed_num_targ_sub_vecs, \
    const bool expect_reduct_obj, \
    const ArrayView<const ConstSubVectorView<SCALAR > > &sub_vecs, \
    const ArrayView<const SubVectorView<SCALAR > > &targ_sub_vecs, \
    const Ptr<const ReductTarget> &reduct_obj \
    ); \
  \



#endif // RTOPPACK_RTOP_T_HELPERS_DEF_HPP
