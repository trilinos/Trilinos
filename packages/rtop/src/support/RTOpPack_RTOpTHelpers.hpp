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

#ifndef RTOPPACK_RTOP_NEW_T_HELPERS_HPP
#define RTOPPACK_RTOP_NEW_T_HELPERS_HPP

//#define RTOPPACK_RTOPT_HELPER_DUMP_OUTPUT

#include "RTOpPack_RTOpT.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_dyn_cast.hpp"

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
  Scalar scalar;
  Index  index;
  ScalarIndex( const Scalar &_scalar, const Index &_index ) : scalar(_scalar), index(_index) {}
};

/** \brief Simple <tt>ReductTarget</tt> subclass for simple scalar objects.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class Scalar>
class ReductTargetScalar : public ReductTarget {
public:
  ReductTargetScalar( const Scalar &scalar = Teuchos::ScalarTraits<Scalar>::zero() ) : scalar_(scalar) {}
  void set( const Scalar &scalar ) { scalar_ = scalar; }
  const Scalar& get() const { return scalar_; }
  std::string description() const
    {
      std::ostringstream oss;
      oss << "RTOpPack::ReductTargetScalar<"<<Teuchos::ScalarTraits<Scalar>::name()<<">{scalar="<<scalar_<<"}";
      return oss.str();
    }
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
  ReductTargetScalarIndex(
    const Scalar &scalar = Teuchos::ScalarTraits<Scalar>::zero()
    ,const Index &index  = Teuchos::ScalarTraits<Index>::zero()
    )
    :scalarIndex_(scalar,index)
    {}
  ReductTargetScalarIndex(
    const ScalarIndex<Scalar> &scalarIndex = ScalarIndex<Scalar>(Teuchos::ScalarTraits<Scalar>::zero(),Teuchos::ScalarTraits<Index>::zero())
    )
    :scalarIndex_(scalarIndex)
    {}
  void set( const ScalarIndex<Scalar> &scalarIndex ) { scalarIndex_ = scalarIndex; }
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
  typedef typename RTOpT<Scalar>::primitive_value_type        primitive_value_type;
  /** \brief . */
  ROpScalarReductionBase(
    const ReductScalar &initReductObjValue = Teuchos::ScalarTraits<ReductScalar>::zero()
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
  void setRawVal( const ReductScalar &rawVal, ReductTarget *reduct_obj ) const
    {
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION( reduct_obj==NULL, std::invalid_argument, "Error!" );
#endif
      using Teuchos::dyn_cast;
      dyn_cast<ReductTargetScalar<ReductScalar> >(*reduct_obj).set(rawVal);
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
      *num_values = Teuchos::PrimitiveTypeTraits<Scalar>::numPrimitiveObjs();
      *num_indexes = 0;
      *num_chars = 0;
    }
  /** \brief . */
  Teuchos::RefCountPtr<ReductTarget> reduct_obj_create() const
    {
      return Teuchos::rcp(new ReductTargetScalar<ReductScalar>(initReductObjValue()));
    }
  /** \brief Default implementation here is for a sum. */
  void reduce_reduct_objs(
    const ReductTarget& in_reduct_obj, ReductTarget* inout_reduct_obj
    ) const
    {
      const ReductScalar in_val    = getRawVal(in_reduct_obj);
      const ReductScalar inout_val = getRawVal(*inout_reduct_obj);
#ifdef RTOPPACK_RTOPT_HELPER_DUMP_OUTPUT
      Teuchos::RefCountPtr<Teuchos::FancyOStream>
        out = Teuchos::VerboseObjectBase::getDefaultOStream();
      Teuchos::OSTab tab(out);
      if(rtop_helpers_dump_all) {
        *out << "\nEntering RTOpPack::ROpScalarReductionBase::reduce_reduct_objs(...) ...\n";
        *out
          << "\nop = " << this->description() << "\n"
          << "in_reduct_obj = " << Teuchos::describe(in_reduct_obj,Teuchos::VERB_EXTREME)
          << "in_val = " << in_val << "\n"
          << "inout_reduct_obj (before reduction) = "
          <<     Teuchos::describe(*inout_reduct_obj,Teuchos::VERB_EXTREME)
          << "inout_val (before reduction) = " << inout_val << "\n";
      }
#endif // RTOPPACK_RTOPT_HELPER_DUMP_OUTPUT
      setRawVal( in_val + inout_val, inout_reduct_obj );
#ifdef RTOPPACK_RTOPT_HELPER_DUMP_OUTPUT
      if(rtop_helpers_dump_all) {
        *out
          << "\ninout_reduct_obj (after reduction) = "
          << Teuchos::describe(*inout_reduct_obj,Teuchos::VERB_EXTREME);
      }
#endif // RTOPPACK_RTOPT_HELPER_DUMP_OUTPUT
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
        num_values==0 || value_data==NULL || num_indexes!=0 || index_data!=NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      Teuchos::PrimitiveTypeTraits<Scalar>::extractPrimitiveObjs( getRawVal(reduct_obj), num_values, value_data );
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
      typedef Teuchos::ScalarTraits<Scalar> ST;
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(
        num_values==0 || value_data==NULL || num_indexes!=0 || index_data!=NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
#ifdef RTOPPACK_RTOPT_HELPER_DUMP_OUTPUT
      Teuchos::RefCountPtr<Teuchos::FancyOStream>
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
      Scalar val = ST::nan();
      Teuchos::PrimitiveTypeTraits<Scalar>::loadPrimitiveObjs( num_values, value_data, &val );
      ROpScalarReductionBaseRawValSetter<ST::isComplex,sizeof(Scalar)==sizeof(ReductScalar),Scalar,ReductScalar>::setRawVal( *this, val, reduct_obj );
#ifdef RTOPPACK_RTOPT_HELPER_DUMP_OUTPUT
      if(rtop_helpers_dump_all) {
        *out << "\nOn output:\n";
        Teuchos::OSTab tab(out);
        *out << "val = " << val << "\n";
        *out << "reduct_op = " << Teuchos::describe(*reduct_obj,Teuchos::VERB_EXTREME);
      }
#endif // RTOPPACK_RTOPT_HELPER_DUMP_OUTPUT
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
    { rtop.setRawVal(Teuchos::ScalarTraits<Scalar>::real(rawVal),reduct_obj); }
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
    const Scalar &initScalarReductObjValue = Teuchos::ScalarTraits<Scalar>::zero()
    ,const Index  &initIndexReductObjValue = Teuchos::ScalarTraits<Index>::zero()
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
  Teuchos::RefCountPtr<ReductTarget> reduct_obj_create() const
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
      Scalar val = Teuchos::ScalarTraits<Scalar>::nan();
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
  ROpIndexReductionBase( const index_type &initReductObjValue = Teuchos::ScalarTraits<index_type>::zero() )
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
  Teuchos::RefCountPtr<ReductTarget> reduct_obj_create() const
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
  ROpScalarTransformationBase( const Scalar &scalarData = Teuchos::ScalarTraits<Scalar>::zero() )
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
    const Scalar &scalarData1 = Teuchos::ScalarTraits<Scalar>::zero()
    ,const Scalar &scalarData2 = Teuchos::ScalarTraits<Scalar>::zero()
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
  const RTOpPack::index_type   subDim  = (SUB_VECS)[0].subDim(); \
  const RTOpPack::index_type   globalOffset = (SUB_VECS)[0].globalOffset(); \
  TEST_FOR_EXCEPT(globalOffset<0); \
  const Scalar                 *v0_val = (SUB_VECS)[0].values(); \
  const ptrdiff_t              v0_s    = (SUB_VECS)[0].stride()

/** \brief Use within an apply_op(...) function implementation where num_vecs==2, num_targ_vecs==0.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
#define RTOP_APPLY_OP_2_0( NUM_VECS, SUB_VECS, NUM_TARG_VECS, TARG_SUB_VECS ) \
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
  const RTOpPack::index_type   subDim  = (SUB_VECS)[0].subDim(); \
  const RTOpPack::index_type   globalOffset = (SUB_VECS)[0].globalOffset(); \
  TEST_FOR_EXCEPT(globalOffset<0); \
  const Scalar                 *v0_val = (SUB_VECS)[0].values(); \
  const ptrdiff_t              v0_s    = (SUB_VECS)[0].stride(); \
  const Scalar                 *v1_val = (SUB_VECS)[1].values(); \
  const ptrdiff_t              v1_s    = (SUB_VECS)[1].stride()

/** \brief Use within an apply_op(...) function implementation where num_vecs==0, num_targ_vecs==1.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
#define RTOP_APPLY_OP_0_1( NUM_VECS, SUB_VECS, NUM_TARG_VECS, TARG_SUB_VECS ) \
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
  const RTOpPack::index_type   subDim  = (TARG_SUB_VECS)[0].subDim(); \
  const RTOpPack::index_type   globalOffset = (TARG_SUB_VECS)[0].globalOffset(); \
  TEST_FOR_EXCEPT(globalOffset<0); \
  Scalar                       *z0_val = (TARG_SUB_VECS)[0].values(); \
  const ptrdiff_t              z0_s    = (TARG_SUB_VECS)[0].stride()

/** \brief Use within an apply_op(...) function implementation where num_vecs==1, num_targ_vecs==1.
 */
#define RTOP_APPLY_OP_1_1( NUM_VECS, SUB_VECS, NUM_TARG_VECS, TARG_SUB_VECS ) \
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
  const RTOpPack::index_type   globalOffset = (SUB_VECS)[0].globalOffset(); \
  TEST_FOR_EXCEPT(globalOffset<0); \
  const Scalar                 *v0_val = (SUB_VECS)[0].values(); \
  const ptrdiff_t              v0_s    = (SUB_VECS)[0].stride(); \
  Scalar                       *z0_val = (TARG_SUB_VECS)[0].values(); \
  const ptrdiff_t              z0_s    = (TARG_SUB_VECS)[0].stride()

/** \brief Use within an apply_op(...) function implementation where num_vecs==2, num_targ_vecs==1.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
#define RTOP_APPLY_OP_2_1( NUM_VECS, SUB_VECS, NUM_TARG_VECS, TARG_SUB_VECS ) \
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
  const RTOpPack::index_type   subDim  = (SUB_VECS)[0].subDim(); \
  const RTOpPack::index_type   globalOffset = (SUB_VECS)[0].globalOffset(); \
  TEST_FOR_EXCEPT(globalOffset<0); \
  const Scalar                 *v0_val = (SUB_VECS)[0].values(); \
  const ptrdiff_t              v0_s    = (SUB_VECS)[0].stride(); \
  const Scalar                 *v1_val = (SUB_VECS)[1].values(); \
  const ptrdiff_t              v1_s    = (SUB_VECS)[1].stride(); \
  Scalar                       *z0_val = (TARG_SUB_VECS)[0].values(); \
  const ptrdiff_t              z0_s    = (TARG_SUB_VECS)[0].stride()


/** \brief Use within an apply_op(...) function implementation where num_vecs==3, num_targ_vecs==0.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
#define RTOP_APPLY_OP_3_0( NUM_VECS, SUB_VECS, NUM_TARG_VECS, TARG_SUB_VECS ) \
  TEST_FOR_EXCEPTION(                                                   \
                     (NUM_VECS)!=3 || (SUB_VECS)==NULL                  \
                     ,RTOpPack::InvalidNumVecs                          \
                     ,"Error, num_vecs="<<(NUM_VECS)<<" not allowed, only num_vecs==3, sub_vecs!=NULL" \
                     );                                                 \
  TEST_FOR_EXCEPTION(                                                   \
                     (NUM_TARG_VECS)!=0 || (TARG_SUB_VECS)!=NULL        \
                     ,RTOpPack::InvalidNumTargVecs                      \
                     ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==0, targ_sub_vecs==NULL" \
                     );                                                 \
  TEST_FOR_EXCEPTION(                                                   \
                     (SUB_VECS)[0].subDim() != (SUB_VECS)[1].subDim()   \
                     || (SUB_VECS)[0].subDim() != (SUB_VECS)[2].subDim() \
                     ||(SUB_VECS)[0].globalOffset() != (SUB_VECS)[1].globalOffset() \
                     ||(SUB_VECS)[0].globalOffset() != (SUB_VECS)[1].globalOffset() \
                     ,IncompatibleVecs                                  \
                     ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==0, targ_sub_vecs==NULL" \
                     );                                                 \
  const RTOpPack::index_type   subDim  = (SUB_VECS)[0].subDim();        \
  const Scalar                 *v0_val = (SUB_VECS)[0].values();        \
  const ptrdiff_t              v0_s    = (SUB_VECS)[0].stride();        \
  const Scalar                 *v1_val = (SUB_VECS)[1].values();        \
  const ptrdiff_t              v1_s    = (SUB_VECS)[1].stride();        \
  const Scalar                 *v2_val = (SUB_VECS)[2].values();        \
  const ptrdiff_t              v2_s    = (SUB_VECS)[2].stride();



#endif // RTOPPACK_RTOP_NEW_T_HELPERS_HPP
