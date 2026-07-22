// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_RTOP_T_HELPERS_DECL_HPP
#define RTOPPACK_RTOP_T_HELPERS_DECL_HPP


//#define RTOPPACK_RTOPT_HELPER_DUMP_OUTPUT


#include <typeinfo>


#include "RTOpPack_RTOpT.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_TypeNameTraits.hpp"


#ifdef RTOPPACK_RTOPT_HELPER_DUMP_OUTPUT
#  include "Teuchos_VerboseObject.hpp"
namespace RTOpPack { extern bool rtop_helpers_dump_all; }
#endif // RTOPPACK_RTOPT_HELPER_DUMP_OUTPUT


namespace RTOpPack {


/** \brief Simple struct for a Scalar and an Ordinal object.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class Scalar>
struct ScalarIndex {
  /** \brief. */
  Scalar scalar;
  /** \brief. */
  Ordinal  index;
  /** \brief. */
  ScalarIndex( const Scalar &_scalar, const Ordinal &_index )
    : scalar(_scalar), index(_index)
    {}
  /** \brief. */
  ScalarIndex()
    : scalar(ScalarTraits<Scalar>::zero()), index(-1)
    {}
};


/** \brief .
 *
 * \relates ScalarIndex
 */
template<class Scalar>
std::ostream& operator<<(std::ostream &out, const ScalarIndex<Scalar> &scalarIndex)
{
  out << "{"<<scalarIndex.scalar<<","<<scalarIndex.index<<"}";
  return out;
}


/** \brief Partial specialization of <tt>PrimitiveTypeTraits</tt> for
 * <tt>ScalarIndex</tt>.
 */
template <class Scalar>
class PrimitiveTypeTraits<Scalar, ScalarIndex<Scalar> > {
public:
  /** \brief . */
  typedef PrimitiveTypeTraits<Scalar,Scalar> ScalarPrimitiveTypeTraits;
  /** \brief . */
  typedef typename ScalarPrimitiveTypeTraits::primitiveType primitiveType;
  /** \brief . */
  static int numPrimitiveObjs() { return ScalarPrimitiveTypeTraits::numPrimitiveObjs(); }
  /** \brief . */
  static int numIndexObjs() { return 1; }
  /** \brief . */
  static int numCharObjs() { return 0; }
  /** \brief . */
  static void extractPrimitiveObjs(
    const ScalarIndex<Scalar> &obj,
    const ArrayView<primitiveType> &primitiveObjs,
    const ArrayView<index_type> &indexObjs,
    const ArrayView<char> &charObjs
    )
    {
      assertInput(primitiveObjs, indexObjs, charObjs);
      ScalarPrimitiveTypeTraits::extractPrimitiveObjs(
        obj.scalar, primitiveObjs, Teuchos::null, Teuchos::null );
      indexObjs[0] = obj.index;
    }
  /** \brief . */
  static void loadPrimitiveObjs(
    const ArrayView<const primitiveType> &primitiveObjs,
    const ArrayView<const index_type> &indexObjs,
    const ArrayView<const char> &charObjs,
    const Ptr<ScalarIndex<Scalar> > &obj
    )
    {
      assertInput(primitiveObjs, indexObjs, charObjs);
      ScalarPrimitiveTypeTraits::loadPrimitiveObjs(
        primitiveObjs, Teuchos::null, Teuchos::null,
        Teuchos::outArg(obj->scalar) );
      obj->index = indexObjs[0];
    }
private:
  static void assertInput(
    const ArrayView<const primitiveType> &primitiveObjs,
    const ArrayView<const index_type> &indexObjs,
    const ArrayView<const char> &charObjs
    )
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPT(
        primitiveObjs.size()!=ScalarPrimitiveTypeTraits::numPrimitiveObjs()
        || indexObjs.size()!=1
        || charObjs.size()!=0 );
#else
      (void)primitiveObjs;
      (void)indexObjs;
      (void)charObjs;
#endif
    }
};


/** \brief Simple <tt>ReductTarget</tt> subclass for simple scalar objects.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class ConcreteReductObj>
class DefaultReductTarget : public ReductTarget {
public:
  /** \brief. */
  DefaultReductTarget( const ConcreteReductObj &concreteReductObj )
    : concreteReductObj_(concreteReductObj)
    {}
  /** \brief. */
  void set( const ConcreteReductObj &concreteReductObj )
    { concreteReductObj_ = concreteReductObj; }
  /** \brief. */
  const ConcreteReductObj& get() const
    { return concreteReductObj_; }
  /** \brief. */
  std::string description() const;
private:
  ConcreteReductObj concreteReductObj_;
};


/** \brief Nonmember constructor.
 *
 * \relates DefaultReductTarget
 */
template<class ConcreteReductObj>
const RCP<DefaultReductTarget<ConcreteReductObj> >
defaultReductTarget( const ConcreteReductObj &concreteReductObj )
{
  return Teuchos::rcp(
    new DefaultReductTarget<ConcreteReductObj>(concreteReductObj));
}


/** \brief Validate the input to an apply_op(...) function.
 *
 * \param op [in] The RTOpT object we are validating apply_op(...) input for.
 *
 * \param allowed_num_sub_vecs [in] The allowed number of subvectors for
 * sub_vecs.size().  If <tt>allowed_num_sub_vecs < 0</tt> then this number is
 * not valided.
 *
 * \param allowed_num_targ_sub_vecs [in] The allowed number of subvectors for
 * targ_sub_vecs.size().  If <tt>allowed_num_targ_sub_vecs < 0</tt> then this
 * number is not valided.
 *
 * \param expect_reduct_obj [in] Determines if reduct_obj must be present or
 * not and the type will be validated as well.
 *
 * \param sub_vecs [in] Input to apply_op(...) being validated
 *
 * \param targ_sub_vecs [in] Input to apply_op(...) being validated, not
 * modified here
 *
 * \param reduct_obj_in [in] Input to apply_op(...) being validated
 */
template<class Scalar>
void validate_apply_op(
  const RTOpT<Scalar> &op,
  const int allowed_num_sub_vecs,
  const int allowed_num_targ_sub_vecs,
  const bool expect_reduct_obj,
  const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
  const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
  const Ptr<const ReductTarget> &reduct_obj
  );


//
// Reduction object operator support
//


/** \brief. */
enum EBasicReductTypes { REDUCT_TYPE_SUM, REDUCT_TYPE_MAX, REDUCT_TYPE_MIN };


/** \brief. */
template<class ConcreteReductObj, int ReductionType>
class BasicReductObjReductionOp {
public:
  /** \brief . */
  inline void operator()(const ConcreteReductObj& in_reduct, ConcreteReductObj& inout_reduct) const
    {
      return in_reduct.this_reduction_type_needs_a_specialization();
    }
};


/** \brief. */
template<class ConcreteReductObj>
class BasicReductObjReductionOp<ConcreteReductObj, REDUCT_TYPE_SUM> {
public:
  /** \brief . */
  inline void operator()(const ConcreteReductObj& in_reduct, ConcreteReductObj& inout_reduct) const
    {
      inout_reduct += in_reduct;
    }
};


/** \brief. */
template<class ConcreteReductObj>
class BasicReductObjReductionOp<ConcreteReductObj, REDUCT_TYPE_MAX> {
public:
  /** \brief . */
  inline void operator()(const ConcreteReductObj& in_reduct, ConcreteReductObj& inout_reduct) const
    {
      inout_reduct = std::max(inout_reduct, in_reduct);
    }
};


/** \brief. */
template<class ConcreteReductObj>
class BasicReductObjReductionOp<ConcreteReductObj, REDUCT_TYPE_MIN> {
public:
  /** \brief . */
  inline void operator()(const ConcreteReductObj& in_reduct, ConcreteReductObj& inout_reduct) const
    {
      inout_reduct = std::min(inout_reduct, in_reduct);
    }
};


/** \brief Null reduction object reduction operator. */
template<class Scalar>
class SumScalarReductObjReduction {
public:
  /** \brief . */
  inline void operator()(const Scalar& in_reduct, Scalar& inout_reduct) const
    {
      inout_reduct += in_reduct;
    }
};
// 2008/07/03: rabart: Above: I have broken from the Thyra guideline of
// passing in-out arguments as const Ptr<Type>& and used raw non-const
// reference Type& instead to allow the user function to be more readable.


/** \brief . */
template<class Scalar, class ConcreteReductObj, class ReductObjReduction>
class ROpScalarReductionWithOpBase : public RTOpT<Scalar>
{
public:

  /** \brief . */
  using RTOpT<Scalar>::apply_op;

  /** \brief . */
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;

  /** \brief . */
  ROpScalarReductionWithOpBase(
    const ConcreteReductObj &initReductObjValue_in = ScalarTraits<Scalar>::zero(),
    ReductObjReduction reductObjReduction_in = ReductObjReduction()
    )
    : initReductObjValue_(initReductObjValue_in),
      reductObjReduction_(reductObjReduction_in)
    {}

  /** \brief . */
  const ConcreteReductObj& getRawVal( const ReductTarget &reduct_obj ) const
    {
      using Teuchos::dyn_cast;
      return dyn_cast<const DefaultReductTarget<ConcreteReductObj> >(reduct_obj).get();
    }

  /** \brief . */
  void setRawVal( const ConcreteReductObj &rawVal,
    const Ptr<ReductTarget> &reduct_obj
    ) const
    {
      using Teuchos::dyn_cast;
      dyn_cast<DefaultReductTarget<ConcreteReductObj> >(*reduct_obj).set(rawVal);
    }

  /** \brief . */
  ConcreteReductObj operator()(const ReductTarget& reduct_obj ) const
    {
      return this->getRawVal(reduct_obj);
    }

  /** @name Overridden from RTOpT */
  //@{

  /** \brief . */
  void get_reduct_type_num_entries_impl(
    const Ptr<int> &num_values,
    const Ptr<int> &num_indexes,
    const Ptr<int> &num_chars
    ) const
    {
      typedef PrimitiveTypeTraits<Scalar, ConcreteReductObj> PTT;
      *num_values = PTT::numPrimitiveObjs();
      *num_indexes = PTT::numIndexObjs();
      *num_chars = PTT::numCharObjs();
    }

  /** \brief . */
  Teuchos::RCP<ReductTarget> reduct_obj_create_impl() const
    {
      return Teuchos::rcp(
        new DefaultReductTarget<ConcreteReductObj>(initReductObjValue()));
    }

  /** \brief . */
  virtual void reduce_reduct_objs_impl(
    const ReductTarget& in_reduct_obj, const Ptr<ReductTarget>& inout_reduct_obj
    ) const
    {
      const ConcreteReductObj scalar_in_reduct_obj = this->getRawVal(in_reduct_obj);
      ConcreteReductObj scalar_inout_reduct_obj = this->getRawVal(*inout_reduct_obj);
      reductObjReduction_(scalar_in_reduct_obj, scalar_inout_reduct_obj);
      this->setRawVal( scalar_inout_reduct_obj, inout_reduct_obj );
    }

  /** \brief . */
  void reduct_obj_reinit_impl( const Ptr<ReductTarget> &reduct_obj ) const
    {
      setRawVal( initReductObjValue(), reduct_obj );
    }

  /** \brief . */
  void extract_reduct_obj_state_impl(
    const ReductTarget &reduct_obj,
    const ArrayView<primitive_value_type> &value_data,
    const ArrayView<index_type> &index_data,
    const ArrayView<char_type> &char_data
    ) const
    {
      typedef PrimitiveTypeTraits<Scalar, ConcreteReductObj> PTT;
      PTT::extractPrimitiveObjs( getRawVal(reduct_obj),
        value_data, index_data, char_data );
    }

  /** \brief . */
  void load_reduct_obj_state_impl(
    const ArrayView<const primitive_value_type> &value_data,
    const ArrayView<const index_type> &index_data,
    const ArrayView<const char_type> &char_data,
    const Ptr<ReductTarget> &reduct_obj
    ) const
    {
      typedef PrimitiveTypeTraits<Scalar, ConcreteReductObj> PTT;
      ConcreteReductObj concrete_reduct_obj;
      PTT::loadPrimitiveObjs( value_data, index_data, char_data,
        Teuchos::outArg(concrete_reduct_obj) );
      this->setRawVal( concrete_reduct_obj, reduct_obj );
    }

  //@}

protected:

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ConcreteReductObj, initReductObjValue );

private:

  ReductObjReduction reductObjReduction_;

};


//
// ROp 1 vector scalar reduction
//


/** \brief Base class for scalar reduction RTOps with one input vector. */
template<class Scalar, class ConcreteReductObj, class EleWiseReduction,
  class ReductObjReduction = SumScalarReductObjReduction<ConcreteReductObj> >
class ROp_1_ScalarReduction
  : public ROpScalarReductionWithOpBase<Scalar, ConcreteReductObj, ReductObjReduction>
{
public:

  /** \brief . */
  typedef ROpScalarReductionWithOpBase<Scalar, ConcreteReductObj, ReductObjReduction> base_t;

  /** \brief . */
  ROp_1_ScalarReduction(
    const ConcreteReductObj &initReductObjValue_in = ConcreteReductObj(),
    EleWiseReduction eleWiseReduction_in = EleWiseReduction(),
    ReductObjReduction reductObjReduction_in = ReductObjReduction()
    )
    : base_t(initReductObjValue_in, reductObjReduction_in),
      eleWiseReduction_(eleWiseReduction_in)
    {}

  /** @name Overridden from RTOpT */
  //@{

  /** \brief . */
  void apply_op_impl(
    const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
    const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
    const Ptr<ReductTarget> &reduct_obj_inout
    ) const
    {
      typedef typename Teuchos::ArrayRCP<const Scalar>::iterator const_iter_t;
      using Teuchos::dyn_cast;

#ifdef TEUCHOS_DEBUG
      validate_apply_op<Scalar>(*this, 1, 0, true,
        sub_vecs, targ_sub_vecs, reduct_obj_inout.getConst());
#else
      (void)targ_sub_vecs;
#endif

      DefaultReductTarget<ConcreteReductObj> &reduct_obj =
        dyn_cast<DefaultReductTarget<ConcreteReductObj> >(*reduct_obj_inout); 
      ConcreteReductObj reduct = reduct_obj.get();
      
      const RTOpPack::index_type subDim = sub_vecs[0].subDim();

      const_iter_t v0_val = sub_vecs[0].values().begin();
      const ptrdiff_t v0_s = sub_vecs[0].stride();

      if ( v0_s == 1 ) {
        for( Teuchos_Ordinal i = 0; i < subDim; ++i )
          eleWiseReduction_( *v0_val++, reduct);
      }
      else {
        for( Teuchos_Ordinal i = 0; i < subDim; ++i, v0_val += v0_s )
          eleWiseReduction_( *v0_val, reduct);
      }
      
      reduct_obj.set(reduct);
      
    }
  
  //@}
  
private:

  EleWiseReduction eleWiseReduction_;

};


/** \brief. */
#define RTOP_ROP_1_REDUCT_SCALAR_CUSTOM_DEFAULT( ROP_CLASS_NAME, REDUCT_SCALAR, \
  BASIC_REDUCT_TYPE_ENUM, CUSTOM_DEFAULT \
  ) \
  \
  template<class Scalar, class ReductScalar> \
  class ROP_CLASS_NAME ## EleWiseReduction \
  { \
  public: \
    inline void operator()( \
      const Scalar &v0, \
      ReductScalar &reduct \
      ) const; \
  }; \
  \
  \
  template<class Scalar> \
  class ROP_CLASS_NAME  \
    : public RTOpPack::ROp_1_ScalarReduction< \
        Scalar, \
        REDUCT_SCALAR, \
        ROP_CLASS_NAME ## EleWiseReduction<Scalar, REDUCT_SCALAR >, \
        RTOpPack::BasicReductObjReductionOp<REDUCT_SCALAR, BASIC_REDUCT_TYPE_ENUM> > \
  { \
    typedef RTOpPack::ROp_1_ScalarReduction< \
      Scalar, \
      REDUCT_SCALAR, \
      ROP_CLASS_NAME ## EleWiseReduction<Scalar, REDUCT_SCALAR >, \
      RTOpPack::BasicReductObjReductionOp<REDUCT_SCALAR, BASIC_REDUCT_TYPE_ENUM> > \
      base_t; \
  public: \
    ROP_CLASS_NAME() \
      : base_t(CUSTOM_DEFAULT) \
      { \
        this->setOpNameBase( #ROP_CLASS_NAME ); \
      } \
  }; \
  \
  \
  template<class Scalar, class ReductScalar> \
  void ROP_CLASS_NAME ## EleWiseReduction<Scalar, ReductScalar>::operator()( \
    const Scalar &v0, ReductScalar &reduct \
    ) const


/** \brief. */
#define RTOP_ROP_1_REDUCT_SCALAR( ROP_CLASS_NAME, REDUCT_SCALAR, \
  BASIC_REDUCT_TYPE_ENUM \
  ) \
  RTOP_ROP_1_REDUCT_SCALAR_CUSTOM_DEFAULT(ROP_CLASS_NAME, REDUCT_SCALAR, \
    BASIC_REDUCT_TYPE_ENUM, Teuchos::ScalarTraits<REDUCT_SCALAR >::zero() )


//
// ROp 1 coordinate-variant vector scalar reduction
//


/** \brief Base class for coordinate-variant scalar reduction RTOps with one
 * input vector. */
template<
  class Scalar,
  class ReductScalar,
  class EleWiseReduction,
  class ReductObjReduction = SumScalarReductObjReduction<ReductScalar>
  >
class ROp_1_CoordVariantScalarReduction
  : public ROpScalarReductionWithOpBase<Scalar, ReductScalar, ReductObjReduction>
{
public:

  /** \name Public members */
  //@{

  /** \brief . */
  typedef ROpScalarReductionWithOpBase<Scalar, ReductScalar, ReductObjReduction> base_t;

  /** \brief . */
  ROp_1_CoordVariantScalarReduction(
    const ReductScalar &initReductObjValue_in = ReductScalar(),
    EleWiseReduction eleWiseReduction_in = EleWiseReduction(),
    ReductObjReduction reductObjReduction_in = ReductObjReduction()
    )
    : base_t(initReductObjValue_in, reductObjReduction_in),
      eleWiseReduction_(eleWiseReduction_in)
    {}

  /** \brief . */
  void setEleWiseReduction(EleWiseReduction eleWiseReduction_in)
    { eleWiseReduction_ = eleWiseReduction_in; }

  /** \brief . */
  const EleWiseReduction& getEleWiseReduction() const
    { return eleWiseReduction_; }
  
  //@}

  /** @name Overridden from RTOpT */
  //@{

  /** \brief This RTOp is NOT coordinate invariant! . */
  bool coord_invariant_impl() const { return false; }

  /** \brief . */
  void apply_op_impl(
    const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
    const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
    const Ptr<ReductTarget> &reduct_obj_inout
    ) const
    {
      typedef typename Teuchos::ArrayRCP<const Scalar>::iterator const_iter_t;
      using Teuchos::dyn_cast;

#ifdef TEUCHOS_DEBUG
      validate_apply_op<Scalar>(*this, 1, 0, true,
        sub_vecs, targ_sub_vecs, reduct_obj_inout);
#else
      (void)targ_sub_vecs;
#endif

      DefaultReductTarget<ReductScalar> &reduct_obj =
        dyn_cast<DefaultReductTarget<ReductScalar> >(*reduct_obj_inout); 
      ReductScalar reduct = reduct_obj.get();
      
      const RTOpPack::index_type subDim = sub_vecs[0].subDim();

      const_iter_t v0_val = sub_vecs[0].values().begin();
      const ptrdiff_t v0_s = sub_vecs[0].stride();

      RTOpPack::index_type global_i = sub_vecs[0].globalOffset();

      if ( v0_s == 1 ) {
        for( Teuchos_Ordinal i = 0; i < subDim; ++i, ++global_i )
          eleWiseReduction_( global_i, *v0_val++, reduct);
      }
      else {
        for( Teuchos_Ordinal i = 0; i < subDim; ++i, v0_val += v0_s, ++global_i )
          eleWiseReduction_( global_i, *v0_val, reduct);
      }
      
      reduct_obj.set(reduct);
      
    }
  
  //@}
  
private:

  EleWiseReduction eleWiseReduction_;

};


//
// ROp 2 vector scalar reduction
//


/** \brief Base class for scalar reduction RTOps with two input vectors. */
template<
  class Scalar,
  class ReductScalar,
  class EleWiseReduction,
  class ReductObjReduction = SumScalarReductObjReduction<ReductScalar>
  >
class ROp_2_ScalarReduction
  : public ROpScalarReductionWithOpBase<Scalar, ReductScalar, ReductObjReduction>
{
public:

  /** \brief . */
  typedef ROpScalarReductionWithOpBase<Scalar, ReductScalar, ReductObjReduction>
    base_t;

  /** \brief . */
  ROp_2_ScalarReduction(
    const ReductScalar &initReductObjValue_in = ReductScalar(),
    EleWiseReduction eleWiseReduction_in = EleWiseReduction(),
    ReductObjReduction reductObjReduction_in = ReductObjReduction()
    )
    : base_t(initReductObjValue_in, reductObjReduction_in),
      eleWiseReduction_(eleWiseReduction_in)
    {}

  /** @name Overridden from RTOpT */
  //@{

  /** \brief . */
  void apply_op_impl(
    const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
    const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
    const Ptr<ReductTarget> &reduct_obj_inout
    ) const
    {
      typedef typename Teuchos::ArrayRCP<const Scalar>::iterator const_iter_t;
      using Teuchos::dyn_cast;

#ifdef TEUCHOS_DEBUG
      validate_apply_op<Scalar>(*this, 2, 0, true,
        sub_vecs, targ_sub_vecs, reduct_obj_inout.getConst());
#else
      (void)targ_sub_vecs;
#endif

      DefaultReductTarget<Scalar> &reduct_obj =
        dyn_cast<DefaultReductTarget<Scalar> >(*reduct_obj_inout); 
      Scalar reduct = reduct_obj.get();

      const RTOpPack::index_type subDim = sub_vecs[0].subDim();

      const_iter_t v0_val = sub_vecs[0].values().begin();
      const ptrdiff_t v0_s = sub_vecs[0].stride();
      const_iter_t v1_val = sub_vecs[1].values().begin();
      const ptrdiff_t v1_s = sub_vecs[1].stride();

      if( v0_s == 1 && v1_s == 1 ) {
        for( Teuchos_Ordinal i = 0; i < subDim; ++i )
          eleWiseReduction_( *v0_val++, *v1_val++, reduct);
      }
      else {
        for( Teuchos_Ordinal i = 0; i < subDim; ++i, v0_val += v0_s, v1_val += v1_s )
          eleWiseReduction_( *v0_val, *v1_val, reduct);
      }

      reduct_obj.set(reduct);

    }

  //@}

private:

  EleWiseReduction eleWiseReduction_;

};


/** \brief Declare and define a concreate reduction RTOp that accepts to
 * vector arguments.
 */
#define RTOP_ROP_2_REDUCT_SCALAR( ROP_CLASS_NAME, REDUCT_SCALAR, \
  BASIC_REDUCT_TYPE_ENUM \
  ) \
  \
  template<class Scalar, class ReductScalar> \
  class ROP_CLASS_NAME ## EleWiseReduction \
  { \
  public: \
    inline void operator()(const Scalar &v0, \
      const Scalar &v1, \
      ReductScalar &reduct \
      ) const; \
  }; \
  \
  \
  template<class Scalar> \
  class ROP_CLASS_NAME  \
    : public RTOpPack::ROp_2_ScalarReduction< \
        Scalar, \
        REDUCT_SCALAR, \
        ROP_CLASS_NAME ## EleWiseReduction<Scalar, REDUCT_SCALAR >, \
        RTOpPack::BasicReductObjReductionOp<REDUCT_SCALAR, BASIC_REDUCT_TYPE_ENUM> > \
  { \
  public: \
    ROP_CLASS_NAME() \
      { \
        this->setOpNameBase( #ROP_CLASS_NAME ); \
        this->initReductObjValue(ScalarTraits<REDUCT_SCALAR >::zero()); \
      } \
  }; \
  \
  template<class Scalar, class ReductScalar> \
  void ROP_CLASS_NAME ## EleWiseReduction<Scalar, ReductScalar>::operator()( \
    const Scalar &v0, const Scalar &v1, ReductScalar &reduct) const


//
// TOp 0 to 1 vector transformation
//


/** \brief Base class for transformations for 0 input and 1 output vector. */
template<class Scalar, class EleWiseTransformation>
class TOp_0_1_Base : public RTOpT<Scalar>
{
public:

  /** \brief . */
  TOp_0_1_Base(
    EleWiseTransformation eleWiseTransformation = EleWiseTransformation()
    )
    : eleWiseTransformation_(eleWiseTransformation)
    {}

  /** @name Overridden from RTOpT */
  //@{

  /** \brief . */
  void apply_op_impl(
    const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
    const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
    const Ptr<ReductTarget> &reduct_obj_inout
    ) const
    {
      typedef typename Teuchos::ArrayRCP<Scalar>::iterator iter_t;

#ifdef TEUCHOS_DEBUG
      validate_apply_op<Scalar>(*this, 0, 1, false,
        sub_vecs, targ_sub_vecs, reduct_obj_inout.getConst());
#else
      (void)sub_vecs;
      (void)reduct_obj_inout;
#endif
      
      const RTOpPack::index_type subDim = targ_sub_vecs[0].subDim();

      iter_t z0_val = targ_sub_vecs[0].values().begin();
      const ptrdiff_t z0_s = targ_sub_vecs[0].stride();

      if ( z0_s == 1 ) {
        for( Teuchos_Ordinal i = 0; i < subDim; ++i )
          eleWiseTransformation_( *z0_val++);
      }
      else {
        for( Teuchos_Ordinal i = 0; i < subDim; ++i, z0_val += z0_s )
          eleWiseTransformation_( *z0_val);
      }
      
    }
  
  //@}
  
private:

  EleWiseTransformation eleWiseTransformation_;

};


/** \brief Base class for coordinate variant transformations for 0 input and 1
 * output vector. */
template<class Scalar, class EleWiseTransformation>
class TOp_0_1_CoordVariantBase : public RTOpT<Scalar>
{
public:

  /** \brief . */
  TOp_0_1_CoordVariantBase(
    EleWiseTransformation eleWiseTransformation = EleWiseTransformation()
    )
    : eleWiseTransformation_(eleWiseTransformation)
    {}

  /** \brief . */
  void setEleWiseTransformation(EleWiseTransformation eleWiseTransformation)
    {
      eleWiseTransformation_ = eleWiseTransformation;
    }

  /** \brief . */
  const EleWiseTransformation& getEleWiseTransformation() const
    {
      return eleWiseTransformation_;
    }
  
  /** @name Overridden from RTOpT */
  //@{

  /** \brief This RTOp is NOT coordinate invariant! . */
  bool coord_invariant_impl() const { return false; }

  /** \brief . */
  void apply_op_impl(
    const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
    const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
    const Ptr<ReductTarget> &reduct_obj_inout
    ) const
    {
      typedef typename Teuchos::ArrayRCP<Scalar>::iterator iter_t;

#ifdef TEUCHOS_DEBUG
      validate_apply_op<Scalar>(*this, 0, 1, false,
        sub_vecs, targ_sub_vecs, reduct_obj_inout.getConst());
#else
      (void)sub_vecs;
      (void)reduct_obj_inout;
#endif
      
      const RTOpPack::index_type subDim = targ_sub_vecs[0].subDim();

      iter_t z0_val = targ_sub_vecs[0].values().begin();
      const ptrdiff_t z0_s = targ_sub_vecs[0].stride();

      RTOpPack::index_type global_i = targ_sub_vecs[0].globalOffset();

      if ( z0_s == 1 ) {
        for( Teuchos_Ordinal i = 0; i < subDim; ++i, ++global_i )
          eleWiseTransformation_(global_i, *z0_val++);
      }
      else {
        for( Teuchos_Ordinal i = 0; i < subDim; ++i, z0_val += z0_s, ++global_i )
          eleWiseTransformation_(global_i, *z0_val);
      }
      
    }
  
  //@}
  
private:

  EleWiseTransformation eleWiseTransformation_;

};


//
// TOp 1 to 1 vector transformation
//


/** \brief Base class for transformations for 1 input and 1 output vector. */
template<class Scalar, class EleWiseTransformation>
class TOp_1_1_Base : public RTOpT<Scalar>
{
public:

  /** \brief . */
  TOp_1_1_Base(
    EleWiseTransformation eleWiseTransformation = EleWiseTransformation()
    )
    : eleWiseTransformation_(eleWiseTransformation)
    {}

  /** @name Overridden from RTOpT */
  //@{

  /** \brief . */
  void apply_op_impl(
    const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
    const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
    const Ptr<ReductTarget> &reduct_obj_inout
    ) const
    {
      typedef typename Teuchos::ArrayRCP<const Scalar>::iterator const_iter_t;
      typedef typename Teuchos::ArrayRCP<Scalar>::iterator iter_t;

#ifdef TEUCHOS_DEBUG
      validate_apply_op<Scalar>(*this, 1, 1, false,
        sub_vecs, targ_sub_vecs, reduct_obj_inout.getConst());
#else
      (void)reduct_obj_inout;
#endif
      
      const RTOpPack::index_type subDim = sub_vecs[0].subDim();

      const_iter_t v0_val = sub_vecs[0].values().begin();
      const ptrdiff_t v0_s = sub_vecs[0].stride();

      iter_t z0_val = targ_sub_vecs[0].values().begin();
      const ptrdiff_t z0_s = targ_sub_vecs[0].stride();

      if ( v0_s == 1 && z0_s == 1 ) {
        for( Teuchos_Ordinal i = 0; i < subDim; ++i )
          eleWiseTransformation_( *v0_val++, *z0_val++);
      }
      else {
        for( Teuchos_Ordinal i = 0; i < subDim; ++i, v0_val += v0_s, z0_val += z0_s )
          eleWiseTransformation_( *v0_val, *z0_val);
      }
      
    }
  
  //@}
  
private:

  EleWiseTransformation eleWiseTransformation_;

};


/** \brief. */
#define RTOP_TOP_1_1( TOP_CLASS_NAME ) \
  \
  template<class Scalar> \
  class TOP_CLASS_NAME ## EleWiseTransformation \
  { \
  public: \
    inline void operator()( const Scalar &v0, Scalar &z0 ) const; \
  }; \
  \
  \
  template<class Scalar> \
  class TOP_CLASS_NAME  \
    : public RTOpPack::TOp_1_1_Base< Scalar, \
        TOP_CLASS_NAME ## EleWiseTransformation<Scalar> > \
  { \
  public: \
    TOP_CLASS_NAME() \
      { \
        this->setOpNameBase( #TOP_CLASS_NAME ); \
      } \
  }; \
  \
  \
  template<class Scalar> \
  void TOP_CLASS_NAME ## EleWiseTransformation<Scalar>::operator()( \
    const Scalar &v0, Scalar &z0 \
    ) const


//
// TOp 2 to 1 vector transformation
//


/** \brief Base class for transformations for 2 input and 1 output vector. */
template<class Scalar, class EleWiseTransformation>
class TOp_2_1_Base : public RTOpT<Scalar>
{
public:

  /** \brief . */
  TOp_2_1_Base(
    EleWiseTransformation eleWiseTransformation = EleWiseTransformation()
    )
    : eleWiseTransformation_(eleWiseTransformation)
    {}

  /** @name Overridden from RTOpT */
  //@{

  /** \brief . */
  void apply_op_impl(
    const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
    const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
    const Ptr<ReductTarget> &reduct_obj_inout
    ) const
    {
      typedef typename Teuchos::ArrayRCP<const Scalar>::iterator const_iter_t;
      typedef typename Teuchos::ArrayRCP<Scalar>::iterator iter_t;

#ifdef TEUCHOS_DEBUG
      validate_apply_op<Scalar>(*this, 2, 1, false,
        sub_vecs, targ_sub_vecs, reduct_obj_inout);
#else
      (void)reduct_obj_inout;
#endif
      
      const RTOpPack::index_type subDim = sub_vecs[0].subDim();

      const_iter_t v0_val = sub_vecs[0].values().begin();
      const ptrdiff_t v0_s = sub_vecs[0].stride();

      const_iter_t v1_val = sub_vecs[1].values().begin();
      const ptrdiff_t v1_s = sub_vecs[1].stride();

      iter_t z0_val = targ_sub_vecs[0].values().begin();
      const ptrdiff_t z0_s = targ_sub_vecs[0].stride();

      if ( v0_s == 1 && v1_s == 1 && z0_s == 1 ) {
        for( Teuchos_Ordinal i = 0; i < subDim; ++i )
          eleWiseTransformation_( *v0_val++, *v1_val++, *z0_val++ );
      }
      else {
        for(
          Teuchos_Ordinal i = 0;
          i < subDim;
          ++i, v0_val += v0_s, v1_val += v1_s, z0_val += z0_s
          )
        {
          eleWiseTransformation_( *v0_val, *v1_val, *z0_val );
        }
      }
      
    }
  

  //@}

private:

  EleWiseTransformation eleWiseTransformation_;

};

/** \brief Base class for transformations for 3 input and 1 output vector. */
template<class Scalar, class EleWiseTransformation>
class TOp_3_1_Base : public RTOpT<Scalar>
{
public:

  /** \brief . */
  TOp_3_1_Base(
    EleWiseTransformation eleWiseTransformation = EleWiseTransformation()
    )
    : eleWiseTransformation_(eleWiseTransformation)
    {}

  /** @name Overridden from RTOpT */
  //@{

  /** \brief . */
  void apply_op_impl(
    const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
    const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
    const Ptr<ReductTarget> &reduct_obj_inout
    ) const
    {
      typedef typename Teuchos::ArrayRCP<const Scalar>::iterator const_iter_t;
      typedef typename Teuchos::ArrayRCP<Scalar>::iterator iter_t;

#ifdef TEUCHOS_DEBUG
      validate_apply_op<Scalar>(*this, 3, 1, false,
        sub_vecs, targ_sub_vecs, reduct_obj_inout);
#endif

      const RTOpPack::index_type subDim = sub_vecs[0].subDim();

      const_iter_t v0_val = sub_vecs[0].values().begin();
      const ptrdiff_t v0_s = sub_vecs[0].stride();

      const_iter_t v1_val = sub_vecs[1].values().begin();
      const ptrdiff_t v1_s = sub_vecs[1].stride();

      const_iter_t v2_val = sub_vecs[2].values().begin();
      const ptrdiff_t v2_s = sub_vecs[2].stride();

      iter_t z0_val = targ_sub_vecs[0].values().begin();
      const ptrdiff_t z0_s = targ_sub_vecs[0].stride();

      if ( v0_s == 1 && v1_s == 1 && v2_s == 1 && z0_s == 1 ) {
        for( Teuchos_Ordinal i = 0; i < subDim; ++i )
          eleWiseTransformation_( *v0_val++, *v1_val++, *v2_val++, *z0_val++ );
      }
      else {
        for(
          Teuchos_Ordinal i = 0;
          i < subDim;
          ++i, v0_val += v0_s, v1_val += v1_s, v2_val += v2_s, z0_val += z0_s
          )
        {
          eleWiseTransformation_( *v0_val, *v1_val, *v2_val, *z0_val );
        }
      }

    }

  //@}
  
private:

  EleWiseTransformation eleWiseTransformation_;

};


} // namespace RTOpPack


#endif // RTOPPACK_RTOP_T_HELPERS_DECL_HPP
