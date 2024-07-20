// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    TEUCHOS_TEST_FOR_EXCEPTION( allowed_num_sub_vecs != as<int>(sub_vecs.size()),
      InvalidNumVecs,
      op_name_str<<": Error, sub_vecs.size()="<<sub_vecs.size()
      <<" != allowed_num_sub_vecs="<<allowed_num_sub_vecs<<"!" );
  }

  if (allowed_num_targ_sub_vecs >= 0) {
    TEUCHOS_TEST_FOR_EXCEPTION( allowed_num_targ_sub_vecs != as<int>(targ_sub_vecs.size()),
      InvalidNumTargVecs,
      op_name_str<<": Error, targ_sub_vecs.size()="<<targ_sub_vecs.size()
      <<" != allowed_num_targ_sub_vecs="<<allowed_num_targ_sub_vecs<<"!" );
  }

  TEUCHOS_TEST_FOR_EXCEPTION( sub_vecs.size() == 0 && targ_sub_vecs.size() == 0,
    InvalidNumVecs,
    op_name_str<<": Error, apply_op(...) must be passed some vectors!"
    );

  const index_type subDim =
    (sub_vecs.size() ? sub_vecs[0].subDim() : targ_sub_vecs[0].subDim());

  const index_type globalOffset =
    (sub_vecs.size() ? sub_vecs[0].globalOffset() : targ_sub_vecs[0].globalOffset());

  for (int k = 0; k < num_sub_vecs; ++k ) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      sub_vecs[k].subDim() != subDim || sub_vecs[k].globalOffset() != globalOffset,
      IncompatibleVecs,
      op_name_str<<": Error, sub_vec["<<k<<"] "
      "(subDim="<<sub_vecs[k].subDim()<<",globalOffset="<<sub_vecs[k].globalOffset()<<")"
      " is not compatible with (subDim="<<subDim<<",globalOffset="<<globalOffset<<")!"
      );
  }

  for (int k = 0; k < num_targ_sub_vecs; ++k ) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      targ_sub_vecs[k].subDim() != subDim || targ_sub_vecs[k].globalOffset() != globalOffset,
      IncompatibleVecs,
      op_name_str<<": Error, sub_vec["<<k<<"] "
      "(subDim="<<targ_sub_vecs[k].subDim()<<",globalOffset="<<targ_sub_vecs[k].globalOffset()<<")"
      " is not compatible with (subDim="<<subDim<<",globalOffset="<<globalOffset<<")!"
      );
  }

  if (expect_reduct_obj) {

    TEUCHOS_TEST_FOR_EXCEPTION( is_null(reduct_obj),
      IncompatibleReductObj,
      op_name_str<<": Error, expected a reduction target object!"
      );

    const RCP<ReductTarget> dummy_reduct_obj = op.reduct_obj_create();

    // mfh 22 Sep 2015: Clang 3.6 warns when typeid's expression may
    // have side effects.  In this case, Ptr::operator* may throw in a
    // debug build if the pointer inside is NULL.  (Throwing an
    // exception counts as a side effect here.)
    //
    // const std::type_info &reduct_obj_type = typeid(*reduct_obj);
    // const std::type_info &dummy_reduct_obj_type = typeid(*dummy_reduct_obj);

    const ReductTarget& reduct_obj_thing = *reduct_obj;
    const std::type_info& reduct_obj_type = typeid (reduct_obj_thing);
    ReductTarget& dummy_reduct_obj_thing = *dummy_reduct_obj;
    const std::type_info& dummy_reduct_obj_type = typeid (dummy_reduct_obj_thing);

    TEUCHOS_TEST_FOR_EXCEPTION( reduct_obj_type != dummy_reduct_obj_type,
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


#define RTOPPACK_RTOPT_HELPERS_INSTANT_SCALAR(SCALAR) \
  \
  RTOPPACK_RTOPT_HELPERS_DEFAULTREDUCTTARGET_INSTANT(SCALAR) \
  \
  RTOPPACK_RTOPT_HELPERS_DEFAULTREDUCTTARGET_INSTANT(ScalarIndex<SCALAR >) \
  \
  RTOPPACK_RTOPT_HELPERS_DEFAULTREDUCTTARGET_INSTANT(SubVectorView<SCALAR >) \
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
