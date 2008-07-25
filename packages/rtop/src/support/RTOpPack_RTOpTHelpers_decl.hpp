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


/** \brief Simple struct for a Scalar and an Index object.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class Scalar>
struct ScalarIndex {
  /** \brief. */
  Scalar scalar;
  /** \brief. */
  Index  index;
  /** \brief. */
  ScalarIndex( const Scalar &_scalar, const Index &_index )
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
      TEST_FOR_EXCEPT(
        primitiveObjs.size()!=ScalarPrimitiveTypeTraits::numPrimitiveObjs()
        || indexObjs.size()!=1
        || charObjs.size()!=0 );
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


/** \brief Simple base class for all reduction operators that return a simple
 * scalar reduction object.
 *
 * Subclasses have to minimally define <tt>apply_op()</tt>.
 * Subclasses should also override <tt>reduce_reduct_objs()</tt> if
 * the reduction is not a simple summation.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class Scalar, class ConcreteReductObj = Scalar>
class ROpScalarReductionBase : virtual public RTOpT<Scalar> {
public:

  /** \brief . */
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;

  /** \brief . */
  ROpScalarReductionBase(
    const ConcreteReductObj &initReductObjValue = ConcreteReductObj()
    )
    :RTOpT<Scalar>(""), initReductObjValue_(initReductObjValue) 
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

  /** @name Overridden from RTOpT */
  //@{

  /** \brief . */
  void get_reduct_type_num_entries(
    int* num_values,
    int* num_indexes,
    int* num_chars
    ) const
    {
      typedef PrimitiveTypeTraits<Scalar, ConcreteReductObj> PTT;
      *num_values = PTT::numPrimitiveObjs();
      *num_indexes = PTT::numIndexObjs();
      *num_chars = PTT::numCharObjs();
    }

  /** \brief . */
  Teuchos::RCP<ReductTarget> reduct_obj_create() const
    {
      return Teuchos::rcp(
        new DefaultReductTarget<ConcreteReductObj>(initReductObjValue()));
    }

  /** \brief Default implementation here is for a sum. */
  void reduce_reduct_objs(
    const ReductTarget& in_reduct_obj, ReductTarget* inout_reduct_obj ) const;

  /** \brief . */
  void reduct_obj_reinit( ReductTarget* reduct_obj ) const
    {
      setRawVal( initReductObjValue(), reduct_obj );
    }

  /** \brief . */
  void extract_reduct_obj_state(
    const ReductTarget &reduct_obj,
    int num_values,
    primitive_value_type value_data[],
    int num_indexes,
    index_type index_data[],
    int num_chars,
    char_type char_data[]
    ) const
    {
      using Teuchos::arrayView;
      typedef PrimitiveTypeTraits<Scalar, ConcreteReductObj> PTT;
      PTT::extractPrimitiveObjs( getRawVal(reduct_obj),
        arrayView(value_data, num_values), arrayView(index_data, num_indexes),
        arrayView(char_data, num_chars) );
    }

  /** \brief . */
  void load_reduct_obj_state(
    int num_values,
    const primitive_value_type value_data[],
    int num_indexes,
    const index_type index_data[],
    int num_chars,
    const char_type char_data[],
    ReductTarget *reduct_obj
    ) const;

  //@}

  /** \name Deprecated. */

  //@{

  /** \brief Deprecated. */
  void setRawVal( const ConcreteReductObj &rawVal, ReductTarget *reduct_obj ) const
    {
      setRawVal(rawVal, Teuchos::ptr(reduct_obj));
    }

  //@}

protected:

  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ConcreteReductObj, initReductObjValue );

};


/** \brief Simple base class for all transformation operators that
 * use a single piece of Scalar data.
 *
 * Subclasses have to minimally define <tt>apply_op()</tt>.
 * Subclasses should also define access functions for changing the
 * scalar data.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class Scalar>
class ROpScalarTransformationBase : virtual public RTOpT<Scalar> {
public:
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;
  /** \brief . */
  ROpScalarTransformationBase( const Scalar &scalarData = ScalarTraits<Scalar>::zero() )
    :RTOpT<Scalar>(""), scalarData_(scalarData)
    {}
  /** @name Overridden from RTOpT */
  //@{
  /** \brief . */
  void get_op_type_num_entries(
    int*  num_values
    ,int* num_indexes
    ,int* num_chars
    ) const
    {
      *num_values = num_values_;
      *num_indexes = 0;
      *num_chars = 0;
    }
  /** \brief . */
  void extract_op_state(
    int                             num_values
    ,primitive_value_type           value_data[]
    ,int                            num_indexes
    ,index_type                     index_data[]
    ,int                            num_chars
    ,char_type                      char_data[]
    ) const
    {
      TEST_FOR_EXCEPT(true); // ToDo: Remove this function for now
    }
  /** \brief . */
  void load_op_state(
    int                           num_values
    ,const primitive_value_type   value_data[]
    ,int                          num_indexes
    ,const index_type             index_data[]
    ,int                          num_chars
    ,const char_type              char_data[]
    )
    {
      TEST_FOR_EXCEPT(true); // ToDo: Remove this function!
    }
  //@}
protected:
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( Scalar, scalarData );
private:
  static const int num_values_;
};


template<class Scalar>
const int ROpScalarTransformationBase<Scalar>::num_values_=PrimitiveTypeTraits<Scalar,Scalar>::numPrimitiveObjs();


/** \brief Simple base class for all transformation operators that
 * use a pair of Scalar data members.
 *
 * Subclasses have to minimally define <tt>apply_op()</tt>.
 * Subclasses should also define access functions for changing the
 * scalar data.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class Scalar>
class ROpScalarScalarTransformationBase : virtual public RTOpT<Scalar> {
public:
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;
  /** \brief . */
  ROpScalarScalarTransformationBase(
    const Scalar &scalarData1 = ScalarTraits<Scalar>::zero()
    ,const Scalar &scalarData2 = ScalarTraits<Scalar>::zero()
    )
    :RTOpT<Scalar>(""), scalarData1_(scalarData1), scalarData2_(scalarData2)
    {}
  /** @name Overridden from RTOpT */
  //@{
  /** \brief . */
  void get_op_type_num_entries(
    int*  num_values
    ,int* num_indexes
    ,int* num_chars
    ) const
    {
      *num_values = num_values_;
      *num_indexes = 0;
      *num_chars = 0;
    }
  /** \brief . */
  void extract_op_state(
    int                             num_values
    ,primitive_value_type           value_data[]
    ,int                            num_indexes
    ,index_type                     index_data[]
    ,int                            num_chars
    ,char_type                      char_data[]
    ) const
    {
      TEST_FOR_EXCEPT(true); // ToDo: Remove this function!
    }
  /** \brief . */
  void load_op_state(
    int                           num_values
    ,const primitive_value_type   value_data[]
    ,int                          num_indexes
    ,const index_type             index_data[]
    ,int                          num_chars
    ,const char_type              char_data[]
    )
    {
      TEST_FOR_EXCEPT(true); // ToDo: Remove this function!
    }
  //@}
protected:
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( Scalar, scalarData1 );
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( Scalar, scalarData2 );
private:
  static const int num_values_;
};


template<class Scalar>
const int ROpScalarScalarTransformationBase<Scalar>::num_values_=2*PrimitiveTypeTraits<Scalar,Scalar>::numPrimitiveObjs();


} // namespace RTOpPack


/** \brief Use within an apply_op(...) function implementation where num_vecs==0, num_targ_vecs==1.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
#define RTOP_APPLY_OP_0_1( NUM_VECS, SUB_VECS, NUM_TARG_VECS, TARG_SUB_VECS ) \
  typedef typename Teuchos::ArrayRCP<Scalar>::iterator iter_t; \
  TEST_FOR_EXCEPTION( \
    (NUM_VECS)!=0 || (SUB_VECS)!=NULL \
    ,RTOpPack::InvalidNumVecs \
    ,"Error, num_vecs="<<(NUM_VECS)<<" not allowed, only num_vecs==0, sub_vecs==NULL" \
    ); \
  TEST_FOR_EXCEPTION( \
    (NUM_TARG_VECS)!=1 || (TARG_SUB_VECS)==NULL \
    ,RTOpPack::InvalidNumTargVecs \
    ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==1, targ_sub_vecs!=NULL" \
    ); \
  const RTOpPack::index_type subDim = (TARG_SUB_VECS)[0].subDim(); \
  const RTOpPack::index_type globalOffset = (TARG_SUB_VECS)[0].globalOffset(); \
  TEST_FOR_EXCEPT(globalOffset<0); \
  iter_t z0_val = (TARG_SUB_VECS)[0].values().begin(); \
  const ptrdiff_t z0_s = (TARG_SUB_VECS)[0].stride()


/** \brief Use within an apply_op(...) function implementation where num_vecs==1, num_targ_vecs==1.
 */
#define RTOP_APPLY_OP_1_1( NUM_VECS, SUB_VECS, NUM_TARG_VECS, TARG_SUB_VECS ) \
  typedef typename Teuchos::ArrayRCP<Scalar>::iterator iter_t; \
  typedef typename Teuchos::ArrayRCP<const Scalar>::iterator const_iter_t; \
  TEST_FOR_EXCEPTION( \
    (NUM_VECS)!=1 || (SUB_VECS)==NULL \
    ,RTOpPack::InvalidNumVecs \
    ,"Error, num_vecs="<<(NUM_VECS)<<" not allowed, only num_vecs==1, sub_vecs!=NULL" \
    ); \
  TEST_FOR_EXCEPTION( \
    (NUM_TARG_VECS)!=1 || (TARG_SUB_VECS)==NULL \
    ,RTOpPack::InvalidNumTargVecs \
    ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==1, targ_sub_vecs!=NULL" \
    ); \
  TEST_FOR_EXCEPTION( \
    (SUB_VECS)[0].subDim() != (TARG_SUB_VECS)[0].subDim() || \
    (SUB_VECS)[0].globalOffset() != (TARG_SUB_VECS)[0].globalOffset() \
    ,IncompatibleVecs \
    ,"Error, sub_vec[0] (subDim="<<(SUB_VECS)[0].subDim()<<",globalOffset="<<(SUB_VECS)[0].globalOffset()<<")" \
    " is not compatible with targ_sub_vec[0] (subDim="<<(TARG_SUB_VECS)[0].subDim()<<",globalOffset="<<(TARG_SUB_VECS)[0].globalOffset()<<")" \
    ); \
  const RTOpPack::index_type   subDim  = (SUB_VECS)[0].subDim(); \
 const RTOpPack::index_type globalOffset = (SUB_VECS)[0].globalOffset(); \
  TEST_FOR_EXCEPT(globalOffset<0); \
  const_iter_t v0_val = (SUB_VECS)[0].values().begin(); \
  const ptrdiff_t v0_s = (SUB_VECS)[0].stride(); \
  iter_t z0_val = (TARG_SUB_VECS)[0].values().begin(); \
  const ptrdiff_t z0_s = (TARG_SUB_VECS)[0].stride()


/** \brief Use within an apply_op(...) function implementation where num_vecs==2, num_targ_vecs==1.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
#define RTOP_APPLY_OP_2_1( NUM_VECS, SUB_VECS, NUM_TARG_VECS, TARG_SUB_VECS ) \
  typedef typename Teuchos::ArrayRCP<Scalar>::iterator iter_t; \
  typedef typename Teuchos::ArrayRCP<const Scalar>::iterator const_iter_t; \
  TEST_FOR_EXCEPTION( \
    (NUM_VECS)!=2 || (SUB_VECS)==NULL \
    ,RTOpPack::InvalidNumVecs \
    ,"Error, num_vecs="<<(NUM_VECS)<<" not allowed, only num_vecs==2, sub_vecs!=NULL" \
    ); \
  TEST_FOR_EXCEPTION( \
    (NUM_TARG_VECS)!=1 || (TARG_SUB_VECS)==NULL \
    ,RTOpPack::InvalidNumTargVecs \
    ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==1, targ_sub_vecs!=NULL" \
    ); \
  TEST_FOR_EXCEPTION( \
    (SUB_VECS)[0].subDim() != (SUB_VECS)[1].subDim() || \
    (SUB_VECS)[0].subDim() != (TARG_SUB_VECS)[0].subDim() || \
    (SUB_VECS)[0].globalOffset() != (SUB_VECS)[1].globalOffset() || \
    (SUB_VECS)[0].globalOffset() != (TARG_SUB_VECS)[0].globalOffset() \
    ,IncompatibleVecs \
    ,"Error, sub_vec[0] (subDim="<<(SUB_VECS)[0].subDim()<<",globalOffset="<<(SUB_VECS)[0].globalOffset()<<")," \
    " sub_vec[1] (subDim="<<(SUB_VECS)[1].subDim()<<",globalOffset="<<(SUB_VECS)[1].globalOffset()<<")," \
    " and targ_sub_vec[0] (subDim="<<(TARG_SUB_VECS)[0].subDim()<<",globalOffset="<<(TARG_SUB_VECS)[0].globalOffset()<<")" \
    " are not compatible." \
    ); \
  const RTOpPack::index_type subDim = (SUB_VECS)[0].subDim(); \
  const RTOpPack::index_type globalOffset = (SUB_VECS)[0].globalOffset(); \
  TEST_FOR_EXCEPT(globalOffset<0); \
  const_iter_t v0_val = (SUB_VECS)[0].values().begin(); \
  const ptrdiff_t v0_s = (SUB_VECS)[0].stride(); \
  const_iter_t v1_val = (SUB_VECS)[1].values().begin(); \
  const ptrdiff_t v1_s = (SUB_VECS)[1].stride(); \
  iter_t z0_val = (TARG_SUB_VECS)[0].values().begin(); \
  const ptrdiff_t z0_s = (TARG_SUB_VECS)[0].stride()















//
// New stuff
//


namespace RTOpPack {


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
class ROpScalarReductionWithOpBase
  : public ROpScalarReductionBase<Scalar, ConcreteReductObj>
{
public:

  /** \brief . */
  ROpScalarReductionWithOpBase(
    const ConcreteReductObj &initReductObjValue = ScalarTraits<Scalar>::zero(),
    ReductObjReduction reductObjReduction = ReductObjReduction()
    )
    : RTOpT<Scalar>(""),
      ROpScalarReductionBase<Scalar, ConcreteReductObj>(initReductObjValue),
      reductObjReduction_(reductObjReduction)
    {}

  /** \brief . */
  ConcreteReductObj operator()(const ReductTarget& reduct_obj ) const
    { return this->getRawVal(reduct_obj); }

  /** @name Overridden from RTOpT */
  //@{

  /** \brief . */
  virtual void reduct_reduct_objs(
    const ReductTarget& in_reduct_obj, const Ptr<ReductTarget>& inout_reduct_obj
    ) const
    {
      const ConcreteReductObj scalar_in_reduct_obj = this->getRawVal(in_reduct_obj);
      ConcreteReductObj scalar_inout_reduct_obj = this->getRawVal(*inout_reduct_obj);
      reductObjReduction_(scalar_in_reduct_obj, scalar_inout_reduct_obj);
      this->setRawVal( scalar_inout_reduct_obj, inout_reduct_obj );
    }

  /** \brief Deprecated. */
  virtual void reduce_reduct_objs(
    const ReductTarget& in_reduct_obj, ReductTarget* inout_reduct_obj
    ) const
    {
      reduct_reduct_objs( in_reduct_obj, Teuchos::ptr(inout_reduct_obj) );
    }

  /** \brief Deprecated. */
  virtual void apply_op(
    const int num_vecs, const ConstSubVectorView<Scalar> sub_vecs[],
    const int num_targ_vecs, const SubVectorView<Scalar> targ_sub_vecs[],
    ReductTarget *reduct_obj
    ) const
    {
      RTOpT<Scalar>::apply_op(
        Teuchos::arrayView(sub_vecs, num_vecs),
        Teuchos::arrayView(targ_sub_vecs, num_targ_vecs),
        Teuchos::ptr(reduct_obj)
        );
    }

  //@}

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
    const ConcreteReductObj &initReductObjValue = ConcreteReductObj(),
    EleWiseReduction eleWiseReduction = EleWiseReduction(),
    ReductObjReduction reductObjReduction = ReductObjReduction()
    )
    : RTOpT<Scalar>(""),
      base_t(initReductObjValue, reductObjReduction),
      eleWiseReduction_(eleWiseReduction)
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
      typedef ScalarTraits<Scalar> ST;

#ifdef TEUCHOS_DEBUG
      validate_apply_op<Scalar>(*this, 1, 0, true,
        sub_vecs, targ_sub_vecs, reduct_obj_inout);
#endif

      DefaultReductTarget<ConcreteReductObj> &reduct_obj =
        dyn_cast<DefaultReductTarget<ConcreteReductObj> >(*reduct_obj_inout); 
      ConcreteReductObj reduct = reduct_obj.get();
      
      const RTOpPack::index_type subDim = sub_vecs[0].subDim();

      const_iter_t v0_val = sub_vecs[0].values().begin();
      const ptrdiff_t v0_s = sub_vecs[0].stride();

      if ( v0_s == 1 ) {
        for( Teuchos_Index i = 0; i < subDim; ++i )
          eleWiseReduction_( *v0_val++, reduct);
      }
      else {
        for( Teuchos_Index i = 0; i < subDim; ++i, v0_val += v0_s )
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
      : RTOpPack::RTOpT<Scalar>( #ROP_CLASS_NAME ), \
        base_t(CUSTOM_DEFAULT) \
      {} \
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

  /** \brief . */
  typedef ROpScalarReductionWithOpBase<Scalar, ReductScalar, ReductObjReduction> base_t;

  /** \brief . */
  ROp_1_CoordVariantScalarReduction(
    const ReductScalar &initReductObjValue = ReductScalar(),
    EleWiseReduction eleWiseReduction = EleWiseReduction(),
    ReductObjReduction reductObjReduction = ReductObjReduction()
    )
    : RTOpT<Scalar>(""),
      base_t(initReductObjValue, reductObjReduction),
      eleWiseReduction_(eleWiseReduction)
    {}

  /** \brief . */
  void setEleWiseReduction(EleWiseReduction eleWiseReduction)
    { eleWiseReduction_ = eleWiseReduction; }

  /** @name Overridden from RTOpT */
  //@{

  /** \brief This RTOp is NOT coordinate invariant! . */
  bool coord_invariant() const { return false; }

  /** \brief . */
  void apply_op_impl(
    const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
    const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
    const Ptr<ReductTarget> &reduct_obj_inout
    ) const
    {
      typedef typename Teuchos::ArrayRCP<const Scalar>::iterator const_iter_t;
      using Teuchos::dyn_cast;
      typedef ScalarTraits<Scalar> ST;

#ifdef TEUCHOS_DEBUG
      validate_apply_op<Scalar>(*this, 1, 0, true,
        sub_vecs, targ_sub_vecs, reduct_obj_inout);
#endif

      DefaultReductTarget<ReductScalar> &reduct_obj =
        dyn_cast<DefaultReductTarget<ReductScalar> >(*reduct_obj_inout); 
      ReductScalar reduct = reduct_obj.get();
      
      const RTOpPack::index_type subDim = sub_vecs[0].subDim();

      const_iter_t v0_val = sub_vecs[0].values().begin();
      const ptrdiff_t v0_s = sub_vecs[0].stride();

      RTOpPack::index_type global_i = sub_vecs[0].globalOffset();

      if ( v0_s == 1 ) {
        for( Teuchos_Index i = 0; i < subDim; ++i, ++global_i )
          eleWiseReduction_( global_i, *v0_val++, reduct);
      }
      else {
        for( Teuchos_Index i = 0; i < subDim; ++i, v0_val += v0_s, ++global_i )
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
    const ReductScalar &initReductObjValue = ReductScalar(),
    EleWiseReduction eleWiseReduction = EleWiseReduction(),
    ReductObjReduction reductObjReduction = ReductObjReduction()
    )
    : RTOpT<Scalar>(""),
      base_t(initReductObjValue, reductObjReduction),
      eleWiseReduction_(eleWiseReduction)
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
      typedef ScalarTraits<Scalar> ST;

#ifdef TEUCHOS_DEBUG
      validate_apply_op<Scalar>(*this, 2, 0, true,
        sub_vecs, targ_sub_vecs, reduct_obj_inout);
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
        for( Teuchos_Index i = 0; i < subDim; ++i )
          eleWiseReduction_( *v0_val++, *v1_val++, reduct);
      }
      else {
        for( Teuchos_Index i = 0; i < subDim; ++i, v0_val += v0_s, v1_val += v1_s )
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
      : RTOpPack::RTOpT<Scalar>( #ROP_CLASS_NAME ) \
      { \
        initReductObjValue(ScalarTraits<REDUCT_SCALAR >::zero()); \
      } \
  }; \
  \
  template<class Scalar, class ReductScalar> \
  void ROP_CLASS_NAME ## EleWiseReduction<Scalar, ReductScalar>::operator()( \
    const Scalar &v0, const Scalar &v1, ReductScalar &reduct) const


} // namespace RTOpPack


#endif // RTOPPACK_RTOP_T_HELPERS_DECL_HPP
