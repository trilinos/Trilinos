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
};


/** \brief Simple <tt>ReductTarget</tt> subclass for simple scalar objects.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class Scalar>
class ReductTargetScalar : public ReductTarget {
public:
  /** \brief. */
  ReductTargetScalar( const Scalar &scalar = ScalarTraits<Scalar>::zero() )
    : scalar_(scalar)
    {}
  /** \brief. */
  void set( const Scalar &scalar ) { scalar_ = scalar; }
  /** \brief. */
  const Scalar& get() const { return scalar_; }
  /** \brief. */
  std::string description() const;
private:
  Scalar scalar_;
};


/** \brief Simple <tt>ReductTarget</tt> subclass for <tt>Scalar,Index</tt> objects.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class Scalar>
class ReductTargetScalarIndex : public ReductTarget {
public:
  /** \brief. */
  ReductTargetScalarIndex()
    :scalarIndex_(ScalarTraits<Scalar>::zero(), ScalarTraits<Index>::zero())
    {}
  /** \brief. */
  ReductTargetScalarIndex(const Scalar &scalar, const Index &index)
    :scalarIndex_(scalar, index)
    {}
  /** \brief. */
  ReductTargetScalarIndex(const ScalarIndex<Scalar> &scalarIndex)
    :scalarIndex_(scalarIndex)
    {}
  /** \brief. */
  void set( const ScalarIndex<Scalar> &scalarIndex ) { scalarIndex_ = scalarIndex; }
  /** \brief. */
  const ScalarIndex<Scalar>& get() const { return scalarIndex_; }
private:
  ScalarIndex<Scalar> scalarIndex_;
};


template<bool isComplex, bool isScalarReductScalar, class Scalar, class ReductScalar>
class ROpScalarReductionBaseRawValSetter;

/** \brief Simple base class for all reduction operators that return a simple
 * scalar reduction object.
 *
 * Subclasses have to minimally define <tt>apply_op()</tt>.
 * Subclasses should also override <tt>reduce_reduct_objs()</tt> if
 * the reduction is not a simple summation.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class Scalar, class ReductScalar = Scalar>
class ROpScalarReductionBase : virtual public RTOpT<Scalar> {
public:
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;
  /** \brief . */
  ROpScalarReductionBase(
    const ReductScalar &initReductObjValue = ScalarTraits<ReductScalar>::zero()
    )
    :RTOpT<Scalar>(""), initReductObjValue_(initReductObjValue) 
    {}
  /** \brief . */
  const ReductScalar& getRawVal( const ReductTarget &reduct_obj ) const
    {
      using Teuchos::dyn_cast;
      return dyn_cast<const ReductTargetScalar<ReductScalar> >(reduct_obj).get();
    }
  /** \brief . */
  void setRawVal( const ReductScalar &rawVal,
    const Ptr<ReductTarget> &reduct_obj
    ) const
    {
      using Teuchos::dyn_cast;
      dyn_cast<ReductTargetScalar<ReductScalar> >(*reduct_obj).set(rawVal);
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
      *num_values = Teuchos::PrimitiveTypeTraits<Scalar>::numPrimitiveObjs();
      *num_indexes = 0;
      *num_chars = 0;
    }
  /** \brief . */
  Teuchos::RCP<ReductTarget> reduct_obj_create() const
    {
      return Teuchos::rcp(
        new ReductTargetScalar<ReductScalar>(initReductObjValue()));
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
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(
        num_values==0 || value_data==NULL || num_indexes!=0 || index_data!=NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      Teuchos::PrimitiveTypeTraits<Scalar>::extractPrimitiveObjs( getRawVal(reduct_obj), num_values, value_data );
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
  void setRawVal( const ReductScalar &rawVal, ReductTarget *reduct_obj ) const
    {
      setRawVal(rawVal, Teuchos::ptr(reduct_obj));
    }
  //@}
protected:
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ReductScalar, initReductObjValue );
};


template<class Scalar, class ReductScalar>
class ROpScalarReductionBaseRawValSetter<true,false,Scalar,ReductScalar> {
public:
  static void setRawVal(
    const ROpScalarReductionBase<Scalar,ReductScalar> &rtop
    ,const Scalar &rawVal, ReductTarget *reduct_obj
    )
    { rtop.setRawVal(ScalarTraits<Scalar>::real(rawVal),reduct_obj); }
};


template<class Scalar, class ReductScalar>
class ROpScalarReductionBaseRawValSetter<true,true,Scalar,ReductScalar> {
public:
  static void setRawVal(
    const ROpScalarReductionBase<Scalar,ReductScalar> &rtop
    ,const Scalar &rawVal, ReductTarget *reduct_obj
    )
    { rtop.setRawVal(rawVal,reduct_obj); }
};


template<bool isScalarReductScalar, class Scalar, class ReductScalar>
class ROpScalarReductionBaseRawValSetter<false,isScalarReductScalar,Scalar,ReductScalar> {
public:
  static void setRawVal(
    const ROpScalarReductionBase<Scalar,ReductScalar> &rtop
    ,const Scalar &rawVal, ReductTarget *reduct_obj
    )
    { rtop.setRawVal(rawVal,reduct_obj); }
};


/** \brief Base class for all reduction operators that return a
 * <tt>ScalarIndex</tt> reduction object.
 *
 * Subclasses have to minimally define <tt>apply_op()</tt> and
 * <tt>reduce_reduct_objs()</tt>.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class Scalar>
class ROpScalarIndexReductionBase : virtual public RTOpT<Scalar> {
public:
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;
  /** \brief . */
  ROpScalarIndexReductionBase(
    const Scalar &initScalarReductObjValue = ScalarTraits<Scalar>::zero()
    ,const Index  &initIndexReductObjValue = ScalarTraits<Index>::zero()
    )
    :RTOpT<Scalar>("")
    ,initScalarReductObjValue_(initScalarReductObjValue)
    ,initIndexReductObjValue_(initIndexReductObjValue)
    {}
  /** \brief . */
  const ScalarIndex<Scalar>& getRawVal( const ReductTarget &reduct_obj ) const
    {
      using Teuchos::dyn_cast;
      return dyn_cast<const ReductTargetScalarIndex<Scalar> >(reduct_obj).get();
    }
  /** \brief . */
  void setRawVal( const ScalarIndex<Scalar> &rawVal, ReductTarget *reduct_obj ) const
    {
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION( reduct_obj==NULL, std::invalid_argument, "Error!" );
#endif
      using Teuchos::dyn_cast;
      dyn_cast<ReductTargetScalarIndex<Scalar> >(*reduct_obj).set(rawVal);
    }
  /** @name Overridden from RTOpT */
  //@{
  /** \brief . */
  void get_reduct_type_num_entries(
    int*   num_values
    ,int*  num_indexes
    ,int*  num_chars
    ) const
    {
      *num_values = num_values_;
      *num_indexes = 1;
      *num_chars = 0;
    }
  /** \brief . */
  Teuchos::RCP<ReductTarget> reduct_obj_create() const
    {
      return Teuchos::rcp(
        new ReductTargetScalarIndex<Scalar>(
          ScalarIndex<Scalar>(initScalarReductObjValue(),initIndexReductObjValue())
          )
        );
    }
  /** \brief . */
  void reduct_obj_reinit( ReductTarget* reduct_obj ) const
    {
      setRawVal( ScalarIndex<Scalar>(initScalarReductObjValue(),initIndexReductObjValue()), reduct_obj );
    }
  /** \brief . */
  void extract_reduct_obj_state(
    const ReductTarget        &reduct_obj
    ,int                      num_values
    ,primitive_value_type     value_data[]
    ,int                      num_indexes
    ,index_type               index_data[]
    ,int                      num_chars
    ,char_type                char_data[]
    ) const
    {
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(
        num_values!=num_values_ || value_data==NULL || num_indexes!=1 || index_data==NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      const ScalarIndex<Scalar> &scalarIndex = getRawVal(reduct_obj);
      Teuchos::PrimitiveTypeTraits<Scalar>::extractPrimitiveObjs( scalarIndex.scalar, num_values, value_data );
      index_data[0] = scalarIndex.index;
    }
  /** \brief . */
  void load_reduct_obj_state(
    int                            num_values
    ,const primitive_value_type    value_data[]
    ,int                           num_indexes
    ,const index_type              index_data[]
    ,int                           num_chars
    ,const char_type               char_data[]
    ,ReductTarget                  *reduct_obj
    ) const
    {
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(
        num_values!=num_values_ || value_data==NULL || num_indexes!=1 || index_data==NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      Scalar val = ScalarTraits<Scalar>::nan();
      Teuchos::PrimitiveTypeTraits<Scalar>::loadPrimitiveObjs( num_values, value_data, &val );
      setRawVal( ScalarIndex<Scalar>(val,index_data[0]), reduct_obj );
    }
  //@}
protected:
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( Scalar, initScalarReductObjValue );
  STANDARD_MEMBER_COMPOSITION_MEMBERS( Index, initIndexReductObjValue );
private:
  static const int num_values_;
};


template<class Scalar>
const int ROpScalarIndexReductionBase<Scalar>::num_values_=Teuchos::PrimitiveTypeTraits<Scalar>::numPrimitiveObjs();


/** \brief Simple base class for all reduction operators that return a simple
 * index reduction object.
 *
 * Subclasses have to minimally define <tt>apply_op()</tt>.
 * Subclasses should also override <tt>reduce_reduct_objs()</tt> if
 * the reduction is not a simple summation.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class Scalar>
class ROpIndexReductionBase : virtual public RTOpT<Scalar> {
public:
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;
  /** \brief . */
  ROpIndexReductionBase( const index_type &initReductObjValue = ScalarTraits<index_type>::zero() )
    :RTOpT<Scalar>(""), initReductObjValue_(initReductObjValue) 
    {}
  /** \brief . */
  index_type getRawVal( const ReductTarget &reduct_obj ) const
    {
      using Teuchos::dyn_cast;
      return dyn_cast<const ReductTargetScalar<index_type> >(reduct_obj).get();
    }
  /** \brief . */
  void setRawVal( const index_type &rawVal, ReductTarget *reduct_obj ) const
    {
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION( reduct_obj==NULL, std::invalid_argument, "Error!" );
#endif
      using Teuchos::dyn_cast;
      dyn_cast<ReductTargetScalar<index_type> >(*reduct_obj).set(rawVal);
    }
  /** @name Overridden from RTOpT */
  //@{
  /** \brief . */
  void get_reduct_type_num_entries(
    int*   num_values
    ,int*  num_indexes
    ,int*  num_chars
    ) const
    {
      *num_values = 0;
      *num_indexes = 1;
      *num_chars = 0;
    }
  /** \brief . */
  Teuchos::RCP<ReductTarget> reduct_obj_create() const
    {
      return Teuchos::rcp(new ReductTargetScalar<index_type>(initReductObjValue()));
    }
  /// Default implementation here is for a sum
  void reduce_reduct_objs(
    const ReductTarget& _in_reduct_obj, ReductTarget* _inout_reduct_obj
    ) const
    {
      using Teuchos::dyn_cast;
      const ReductTargetScalar<index_type> &in_reduct_obj    = dyn_cast<const ReductTargetScalar<index_type> >(_in_reduct_obj); 
      ReductTargetScalar<index_type>       &inout_reduct_obj = dyn_cast<ReductTargetScalar<index_type> >(*_inout_reduct_obj); 
      inout_reduct_obj.set( inout_reduct_obj.get() + in_reduct_obj.get() );
    }
  /** \brief . */
  void reduct_obj_reinit( ReductTarget* reduct_obj ) const
    {
      setRawVal( initReductObjValue(), reduct_obj );
    }
  /** \brief . */
  void extract_reduct_obj_state(
    const ReductTarget        &reduct_obj
    ,int                      num_values
    ,primitive_value_type     value_data[]
    ,int                      num_indexes
    ,index_type               index_data[]
    ,int                      num_chars
    ,char_type                char_data[]
    ) const
    {
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(
        num_values!=0 || value_data!=NULL || num_indexes!=1 || index_data==NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      index_data[0] = getRawVal(reduct_obj);
    }
  /** \brief . */
  void load_reduct_obj_state(
    int                            num_values
    ,const primitive_value_type    value_data[]
    ,int                           num_indexes
    ,const index_type              index_data[]
    ,int                           num_chars
    ,const char_type               char_data[]
    ,ReductTarget                  *reduct_obj
    ) const
    {
      using Teuchos::dyn_cast;
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(
        num_values!=0 || value_data!=NULL || num_indexes!=1 || index_data==NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      setRawVal( index_data[0], reduct_obj );
    }
  //@}
protected:
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( index_type, initReductObjValue );
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
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(
        num_values!=num_values_ || value_data==NULL || num_indexes!=0 || index_data!=NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      Teuchos::PrimitiveTypeTraits<Scalar>::extractPrimitiveObjs( scalarData_, num_values, value_data );
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
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(
        num_values!=num_values_ || value_data==NULL || num_indexes!=0 || index_data!=NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      Teuchos::PrimitiveTypeTraits<Scalar>::loadPrimitiveObjs( num_values, value_data, &scalarData_ );
    }
  //@}
protected:
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( Scalar, scalarData );
private:
  static const int num_values_;
};


template<class Scalar>
const int ROpScalarTransformationBase<Scalar>::num_values_=Teuchos::PrimitiveTypeTraits<Scalar>::numPrimitiveObjs();


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
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(
        num_values!=num_values_ || value_data==NULL || num_indexes!=0 || index_data!=NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      Teuchos::PrimitiveTypeTraits<Scalar>::extractPrimitiveObjs( scalarData1_, num_values/2, value_data );
      Teuchos::PrimitiveTypeTraits<Scalar>::extractPrimitiveObjs( scalarData2_, num_values/2, value_data + num_values/2 );
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
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(
        num_values!=num_values_ || value_data==NULL || num_indexes!=0 || index_data!=NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      Teuchos::PrimitiveTypeTraits<Scalar>::loadPrimitiveObjs( num_values/2, value_data, &scalarData1_ );
      Teuchos::PrimitiveTypeTraits<Scalar>::loadPrimitiveObjs( num_values/2, value_data+num_values/2, &scalarData2_ );
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
const int ROpScalarScalarTransformationBase<Scalar>::num_values_=2*Teuchos::PrimitiveTypeTraits<Scalar>::numPrimitiveObjs();


/** \brief Do a transformation and reduce to a bool. Needed for the NVector
 * adapters for the SUNDIALS interface.
 *
 * \author K. Long
 */
template<class Scalar>
class RTOpBoolReduceAndTransform 
  : public ROpIndexReductionBase<Scalar>,
    public ROpScalarTransformationBase<Scalar>
{
public:
  /** \brief . */
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;
  /** \brief . */
  RTOpBoolReduceAndTransform()
    : RTOpT<Scalar>(""), 
      ROpIndexReductionBase<Scalar>(1),
      ROpScalarTransformationBase<Scalar>() 
    {;}
  
  /** \brief . */
  virtual ~RTOpBoolReduceAndTransform(){;}
  
  /** \brief . */
  index_type operator()(const ReductTarget& reduct_obj ) const 
    { return this->getRawVal(reduct_obj); }
  /** \brief Default implementation here is for a logical AND. */
  void reduce_reduct_objs(const ReductTarget& in_reduct_obj, 
                          ReductTarget* inout_reduct_obj) const
    {
      const index_type in_val    = this->getRawVal(in_reduct_obj);
      const index_type inout_val = this->getRawVal(*inout_reduct_obj);
      this->setRawVal( in_val && inout_val, inout_reduct_obj );
    }
};


} // namespace RTOpPack


/** \brief Use within an apply_op(...) function implementation where num_vecs==1, num_targ_vecs==0.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
#define RTOP_APPLY_OP_1_0( NUM_VECS, SUB_VECS, NUM_TARG_VECS, TARG_SUB_VECS ) \
  typedef typename Teuchos::ArrayRCP<const Scalar>::iterator const_iter_t; \
  TEST_FOR_EXCEPTION( \
    (NUM_VECS)!=1 || (SUB_VECS)==NULL \
    ,RTOpPack::InvalidNumVecs \
    ,"Error, num_vecs="<<(NUM_VECS)<<" not allowed, only num_vecs==1, sub_vecs!=NULL" \
    ); \
  TEST_FOR_EXCEPTION( \
    (NUM_TARG_VECS)!=0 || (TARG_SUB_VECS)!=NULL \
    ,RTOpPack::InvalidNumTargVecs \
    ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==0, targ_sub_vecs==NULL" \
    ); \
  const RTOpPack::index_type subDim  = (SUB_VECS)[0].subDim(); \
  const RTOpPack::index_type globalOffset = (SUB_VECS)[0].globalOffset(); \
  TEST_FOR_EXCEPT(globalOffset<0); \
  const_iter_t v0_val = (SUB_VECS)[0].values().begin(); \
  const ptrdiff_t v0_s = (SUB_VECS)[0].stride()


/** \brief Use within an apply_op(...) function implementation where num_vecs==2, num_targ_vecs==0.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
#define RTOP_APPLY_OP_2_0( NUM_VECS, SUB_VECS, NUM_TARG_VECS, TARG_SUB_VECS ) \
  typedef typename Teuchos::ArrayRCP<const Scalar>::iterator const_iter_t; \
  TEST_FOR_EXCEPTION( \
    (NUM_VECS)!=2 || (SUB_VECS)==NULL \
    ,RTOpPack::InvalidNumVecs \
    ,"Error, num_vecs="<<(NUM_VECS)<<" not allowed, only num_vecs==1, sub_vecs!=NULL" \
    ); \
  TEST_FOR_EXCEPTION( \
    (NUM_TARG_VECS)!=0 || (TARG_SUB_VECS)!=NULL \
    ,RTOpPack::InvalidNumTargVecs \
    ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==0, targ_sub_vecs==NULL" \
    ); \
  TEST_FOR_EXCEPTION( \
    (SUB_VECS)[0].subDim() != (SUB_VECS)[1].subDim() || \
    (SUB_VECS)[0].globalOffset() != (SUB_VECS)[1].globalOffset() \
    ,IncompatibleVecs \
    ,"Error, sub_vec[0] (subDim="<<(SUB_VECS)[0].subDim()<<",globalOffset="<<(SUB_VECS)[0].globalOffset()<<")" \
    " is not compatible with sub_vec[1] (subDim="<<(SUB_VECS)[1].subDim()<<",globalOffset="<<(SUB_VECS)[1].globalOffset()<<")" \
    ); \
  const RTOpPack::index_type subDim  = (SUB_VECS)[0].subDim(); \
  const RTOpPack::index_type globalOffset = (SUB_VECS)[0].globalOffset(); \
  TEST_FOR_EXCEPT(globalOffset<0); \
  const_iter_t v0_val = (SUB_VECS)[0].values().begin(); \
  const ptrdiff_t v0_s = (SUB_VECS)[0].stride(); \
  const_iter_t v1_val = (SUB_VECS)[1].values().begin(); \
  const ptrdiff_t v1_s = (SUB_VECS)[1].stride()


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


/** \brief Use within an apply_op(...) function implementation where
 * num_vecs==3, num_targ_vecs==0.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
#define RTOP_APPLY_OP_3_0( NUM_VECS, SUB_VECS, NUM_TARG_VECS, TARG_SUB_VECS ) \
  typedef typename Teuchos::ArrayRCP<Scalar>::iterator iter_t; \
  typedef typename Teuchos::ArrayRCP<const Scalar>::iterator const_iter_t; \
  TEST_FOR_EXCEPTION( \
    (NUM_VECS)!=3 || (SUB_VECS)==NULL \
    ,RTOpPack::InvalidNumVecs \
    ,"Error, num_vecs="<<(NUM_VECS)<<" not allowed, only num_vecs==3, sub_vecs!=NULL" \
    ); \
  TEST_FOR_EXCEPTION( \
    (NUM_TARG_VECS)!=0 || (TARG_SUB_VECS)!=NULL \
    ,RTOpPack::InvalidNumTargVecs \
    ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==0, targ_sub_vecs==NULL" \
    ); \
  TEST_FOR_EXCEPTION( \
    (SUB_VECS)[0].subDim() != (SUB_VECS)[1].subDim() \
    || (SUB_VECS)[0].subDim() != (SUB_VECS)[2].subDim() \
    ||(SUB_VECS)[0].globalOffset() != (SUB_VECS)[1].globalOffset() \
    ||(SUB_VECS)[0].globalOffset() != (SUB_VECS)[1].globalOffset() \
    ,IncompatibleVecs \
    ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==0, targ_sub_vecs==NULL" \
    ); \
  const RTOpPack::index_type subDim = (SUB_VECS)[0].subDim(); \
  const_iter_t v0_val = (SUB_VECS)[0].values().begin(); \
  const ptrdiff_t v0_s = (SUB_VECS)[0].stride(); \
  const_iter_t v1_val = (SUB_VECS)[1].values().begin(); \
  const ptrdiff_t v1_s = (SUB_VECS)[1].stride(); \
  const_iter_t v2_val = (SUB_VECS)[2].values().begin(); \
  const ptrdiff_t v2_s = (SUB_VECS)[2].stride();






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
template<class ReductScalar, int ReductionType>
class BasicReductObjReductionOp {
public:
  /** \brief . */
  inline void operator()(const ReductScalar& in_reduct, ReductScalar& inout_reduct) const
    {
      return in_reduct.this_reduction_type_needs_a_specialization();
    }
};


/** \brief. */
template<class ReductScalar>
class BasicReductObjReductionOp<ReductScalar, REDUCT_TYPE_SUM> {
public:
  /** \brief . */
  inline void operator()(const ReductScalar& in_reduct, ReductScalar& inout_reduct) const
    {
      inout_reduct += in_reduct;
    }
};


/** \brief. */
template<class ReductScalar>
class BasicReductObjReductionOp<ReductScalar, REDUCT_TYPE_MAX> {
public:
  /** \brief . */
  inline void operator()(const ReductScalar& in_reduct, ReductScalar& inout_reduct) const
    {
      inout_reduct = std::max(inout_reduct, in_reduct);
    }
};


/** \brief. */
template<class ReductScalar>
class BasicReductObjReductionOp<ReductScalar, REDUCT_TYPE_MIN> {
public:
  /** \brief . */
  inline void operator()(const ReductScalar& in_reduct, ReductScalar& inout_reduct) const
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
template<class Scalar, class ReductScalar, class ReductObjReduction>
class ROpScalarReductionWithOpBase
  : public ROpScalarReductionBase<Scalar, ReductScalar>
{
public:

  /** \brief . */
  using RTOpT<Scalar>::apply_op;

  /** \brief . */
  ROpScalarReductionWithOpBase(
    const ReductScalar &initReductObjValue = ScalarTraits<Scalar>::zero(),
    ReductObjReduction reductObjReduction = ReductObjReduction()
    )
    : RTOpT<Scalar>(""),
      ROpScalarReductionBase<Scalar, ReductScalar>(initReductObjValue),
      reductObjReduction_(reductObjReduction)
    {}

  /** \brief . */
  ReductScalar operator()(const ReductTarget& reduct_obj ) const
    { return this->getRawVal(reduct_obj); }

  /** @name Overridden from RTOpT */
  //@{

  /** \brief . */
  virtual void reduct_reduct_objs(
    const ReductTarget& in_reduct_obj, const Ptr<ReductTarget>& inout_reduct_obj
    ) const
    {
      const ReductScalar scalar_in_reduct_obj = this->getRawVal(in_reduct_obj);
      ReductScalar scalar_inout_reduct_obj = this->getRawVal(*inout_reduct_obj);
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
      apply_op(
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
template<class Scalar, class ReductScalar, class EleWiseReduction,
  class ReductObjReduction = SumScalarReductObjReduction<ReductScalar> >
class ROp_1_ScalarReduction
  : public ROpScalarReductionWithOpBase<Scalar, ReductScalar, ReductObjReduction>
{
public:

  /** \brief . */
  using RTOpT<Scalar>::apply_op;

  /** \brief . */
  typedef ROpScalarReductionWithOpBase<Scalar, ReductScalar, ReductObjReduction> base_t;

  /** \brief . */
  ROp_1_ScalarReduction(
    const ReductScalar &initReductObjValue = ScalarTraits<Scalar>::zero(),
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
  void apply_op(
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

      ReductTargetScalar<ReductScalar> &reduct_obj =
        dyn_cast<ReductTargetScalar<ReductScalar> >(*reduct_obj_inout); 
      ReductScalar reduct = reduct_obj.get();
      
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
    : public ROp_1_ScalarReduction< \
        Scalar, \
        REDUCT_SCALAR, \
        ROP_CLASS_NAME ## EleWiseReduction<Scalar, REDUCT_SCALAR>, \
        BasicReductObjReductionOp<REDUCT_SCALAR, BASIC_REDUCT_TYPE_ENUM> > \
  { \
    typedef ROp_1_ScalarReduction< \
      Scalar, \
      REDUCT_SCALAR, \
      ROP_CLASS_NAME ## EleWiseReduction<Scalar, REDUCT_SCALAR>, \
      BasicReductObjReductionOp<REDUCT_SCALAR, BASIC_REDUCT_TYPE_ENUM> > \
      base_t; \
  public: \
    ROP_CLASS_NAME() \
      : RTOpT<Scalar>( #ROP_CLASS_NAME ), \
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
    BASIC_REDUCT_TYPE_ENUM, Teuchos::ScalarTraits<REDUCT_SCALAR>::zero() )


//
// ROp 2 vector scalar reduction
//


/** \brief Base class for scalar reduction RTOps with two input vectors. */
template<class Scalar, class ReductScalar, class EleWiseReduction,
  class ReductObjReduction = SumScalarReductObjReduction<ReductScalar> >
class ROp_2_ScalarReduction
  : public ROpScalarReductionWithOpBase<Scalar, ReductScalar, ReductObjReduction>
{
public:

  /** \brief . */
  using RTOpT<Scalar>::apply_op;

  /** \brief . */
  typedef ROpScalarReductionWithOpBase<Scalar, ReductScalar, ReductObjReduction> base_t;

  /** \brief . */
  ROp_2_ScalarReduction(
    EleWiseReduction eleWiseReduction = EleWiseReduction(),
    ReductObjReduction reductObjReduction = ReductObjReduction()
    )
    : RTOpT<Scalar>(""),
      base_t(ScalarTraits<ReductScalar>::zero(), reductObjReduction),
      eleWiseReduction_(eleWiseReduction)
    {}

  /** @name Overridden from RTOpT */
  //@{

  /** \brief . */
  void apply_op(
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

      ReductTargetScalar<Scalar> &reduct_obj =
        dyn_cast<ReductTargetScalar<Scalar> >(*reduct_obj_inout); 
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
    : public ROp_2_ScalarReduction< \
        Scalar, \
        REDUCT_SCALAR, \
        ROP_CLASS_NAME ## EleWiseReduction<Scalar, REDUCT_SCALAR>, \
        BasicReductObjReductionOp<REDUCT_SCALAR, BASIC_REDUCT_TYPE_ENUM> > \
  { \
  public: \
    ROP_CLASS_NAME() \
      : RTOpT<Scalar>( #ROP_CLASS_NAME ) \
      {} \
  }; \
  \
  template<class Scalar, class ReductScalar> \
  void ROP_CLASS_NAME ## EleWiseReduction<Scalar, ReductScalar>::operator()( \
    const Scalar &v0, const Scalar &v1, ReductScalar &reduct) const


} // namespace RTOpPack


#endif // RTOPPACK_RTOP_T_HELPERS_DECL_HPP
