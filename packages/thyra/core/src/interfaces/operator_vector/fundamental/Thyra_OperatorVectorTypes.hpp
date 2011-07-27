// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_OPERATOR_VECTOR_TYPES_HPP
#define THYRA_OPERATOR_VECTOR_TYPES_HPP

#include "Thyra_ConfigDefs.hpp"
#include "RTOpPack_Types.hpp"
#include "Teuchos_Range1D.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TypeNameTraits.hpp"


namespace Thyra {

// Using declarations from Teuchos

/** \brief .
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
using Teuchos::Ptr;
/** \brief .
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
using Teuchos::RCP;
/** \brief .
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
using Teuchos::Array;
/** \brief .
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
using Teuchos::ArrayView;
/** \brief .
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
using Teuchos::ArrayRCP;
/** \brief .
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
using Teuchos::Tuple;
/** \brief .
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
typedef Teuchos::Range1D Range1D;
/** \brief .
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
using Teuchos::FancyOStream;
/** \brief .
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
using Teuchos::ParameterList;
/** \brief .
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
using Teuchos::ScalarTraits;
/** \brief .
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
using Teuchos::typeName;
/** \brief .
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
using Teuchos::TypeNameTraits;
/** \brief Type for the dimension of a vector space.
`*
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
typedef Teuchos::Ordinal Ordinal;

/** \brief Deprecated: Use Ordinal instead!
`*
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
THYRA_DEPRECATED typedef Ordinal Index;


/** \brief Enumeration for determining how a linear operator is applied.
`*
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
enum EConj {
  NONCONJ_ELE     ///< Use the linear operator with non-conjugate elements.
  ,CONJ_ELE       ///< Use the linear operator with conjugate elements.
};


/** \brief Return a string name for a <tt>EOpTransp</tt> value.
`*
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
inline
const char* toString(EConj conj)
{
  switch(conj) {
    case NONCONJ_ELE:    return "NONCONJ_ELE";
    case CONJ_ELE:       return "CONJ_ELE";
    default: TEST_FOR_EXCEPT(true);
  }
  return "BAD"; // Should never be called!
}


/** \brief Enumeration for determining how a linear operator is applied.
`*
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
enum EOpTransp {
  /** \brief Use the non-transposed operator. */
  NOTRANS,
  /** \brief Use the non-transposed operator with complex-conjugate elements
   * (same as <tt>NOTRANS</tt> for real scalar types).
   */
  CONJ,
  /** \brief Use the transposed operator. */ 
  TRANS,
  /** \brief Use the transposed operator with complex-conjugate clements (same
   * as <tt>TRANS</tt> for real scalar types).
   */
  CONJTRANS
};


/** \brief Deprecated (use EOpTransp instead).
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
THYRA_DEPRECATED typedef EOpTransp ETransp;


/** \brief Return a string name for a <tt>EOpTransp</tt> value.
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
inline
const char* toString(EOpTransp transp)
{
  switch(transp) {
    case NOTRANS:    return "NOTRANS";
    case CONJ:       return "CONJ";
    case TRANS:      return "TRANS";
    case CONJTRANS:  return "CONJTRANS";
    default: TEST_FOR_EXCEPT(true);
  }
  return "BAD"; // Should never be called!
}


/** \brief Return <tt>NOTRANS</tt> or <tt>TRANS</tt> for real scalar valued
 * operators and this also is used for determining structural transpose.
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
inline
EOpTransp real_trans(EOpTransp transp)
{
  switch(transp) {
    case NOTRANS:    return NOTRANS;
    case CONJ:       return NOTRANS;
    case TRANS:      return TRANS;
    case CONJTRANS:  return TRANS;
    default: TEST_FOR_EXCEPT(true);
  }
  return NOTRANS; // Will never be called!
}


/** \brief Perform a not operation on an EOpTransp value
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
inline 
EOpTransp not_trans( EOpTransp transp )
{
  switch(transp) {
    case NOTRANS:    return TRANS;
    case CONJ:       return CONJTRANS;
    case TRANS:      return CONJ;
    case CONJTRANS:  return NOTRANS;
    default: TEST_FOR_EXCEPT(true);
  }
  return NOTRANS; // Will never be called!
}


/** \brief Combine two transpose arguments
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
inline
EOpTransp trans_trans( EOpTransp trans1, EOpTransp trans2 )
{
  if( trans1 == trans2 )
    return NOTRANS;
  if( trans1 == NOTRANS )
    return trans2;
  if( trans2 == NOTRANS )
    return trans1;
  if( ( trans1 == CONJ && trans2 == TRANS ) || ( trans2 == CONJ && trans1 == TRANS ) )
    return CONJTRANS;
  if( ( trans1 == TRANS && trans2 == CONJTRANS ) || ( trans2 == TRANS && trans1 == CONJTRANS ) )
    return CONJ;
  if( ( trans1 == CONJ && trans2 == CONJTRANS ) || ( trans2 == CONJ && trans1 == CONJTRANS ) )
    return TRANS;
  else
    TEST_FOR_EXCEPT(true);
  return NOTRANS; // Will never be executed!
}


/** \brief Convert from <tt>EOpTransp</tt> to <tt>EConj</tt>.
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
inline
EConj transToConj( EOpTransp trans )
{
  switch(trans) {
    case NOTRANS:    return NONCONJ_ELE;
    case CONJ:       return CONJ_ELE;
    case TRANS:      return NONCONJ_ELE;
    case CONJTRANS:  return CONJ_ELE;
    default: TEST_FOR_EXCEPT(true);
  }
  return NONCONJ_ELE; // Will never be called!
}

/** \brief Convert from <tt>EConj</tt> to <tt>EOpTransp</tt> for forward apply.
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
inline
EOpTransp applyConjToTrans( EConj conj )
{
  switch(conj) {
    case NONCONJ_ELE: return NOTRANS;
    case CONJ_ELE:    return CONJ;
    default: TEST_FOR_EXCEPT(true);
  }
  return NOTRANS; // Will never be called!
}


/** \brief Convert from <tt>EConj</tt> to <tt>EOpTransp</tt> for forward apply.
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
inline
EOpTransp applyTransposeConjToTrans( EConj conj )
{
  switch(conj) {
    case NONCONJ_ELE: return TRANS;
    case CONJ_ELE:    return CONJTRANS;
    default: TEST_FOR_EXCEPT(true);
  }
  return NOTRANS; // Will never be called!
}

/** \brief Determines if a view is a direct view of data or a detached copy of
 * data.
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
enum EViewType {
  VIEW_TYPE_DIRECT   ///< The view is a direct view of data and no copies are made.
  ,VIEW_TYPE_DETACHED  ///< The view is a detached copy of the data.
};


/** \brief Determine if data is unit stride or non-unit stride.
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
enum EStrideType {
  /// The stride between elements in an array is one.
  STRIDE_TYPE_UNIT,
  /// The stride between elements in an array is greater than or equal to one.
  STRIDE_TYPE_NONUNIT
};


namespace Exceptions {


/** \brief Thrown if any member functions are called before initialize() has
 * been called.
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
class UnInitialized : public std::logic_error
{public: UnInitialized(const std::string& what_arg)
   : std::logic_error(what_arg) {}};


/** \brief Thrown if vector spaces are incompatible.
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
class IncompatibleVectorSpaces : public std::logic_error
{public:
  IncompatibleVectorSpaces(const std::string& what_arg)
    : std::logic_error(what_arg) {}
};


/** \biref Thrown if the argument <tt>M_trans</tt> is not supported.
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
class OpNotSupported : public std::logic_error
{public: OpNotSupported(const std::string& what_arg)
   : std::logic_error(what_arg) {}};


} // namespace Exceptions


// Fundamental ANA operator/vector interface classes


template<class Scalar> class VectorSpaceFactoryBase;
template<class Scalar> class VectorSpaceBase;
template<class Scalar> class LinearOpBase;
template<class Scalar> class MultiVectorBase;
template<class Scalar> class VectorBase;


} // end namespace Thyra


#endif // THYRA_OPERATOR_VECTOR_TYPES_HPP
