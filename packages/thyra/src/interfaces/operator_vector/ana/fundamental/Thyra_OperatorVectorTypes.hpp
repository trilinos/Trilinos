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

#include "RTOpPack_Types.hpp"
#include "Teuchos_Range1D.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_TypeNameTraits.hpp"

namespace Thyra {

using Teuchos::typeName;
using Teuchos::TypeNameTraits;

/** \defgroup Thyra_Op_Vec_BasicTypes_grp Basic Thyra types.
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
//@{

//
// Basic types
//

/// Type for the dimension of a vector space
typedef Teuchos::Range1D::Index  Index;

/// Type for a range of indices
typedef Teuchos::Range1D   Range1D;

/** \brief Enumeration for determining how a linear operator is applied.
 */
enum EConj {
  NONCONJ_ELE     ///< Use the linear operator with non-conjugate elements.
  ,CONJ_ELE       ///< Use the linear operator with conjugate elements.
};

/** \brief Return a string name for a <tt>ETransp</tt> value.
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

/** \brief Enumeration for determining how a linear operator (with a single
    scalar * type) is applied.
*/
enum ETransp {
  NOTRANS     ///< Use the non-transposed operator
  ,CONJ       ///< Use the non-transposed operator with complex-conjugate elements (same as <tt>NOTRANS</tt> for real scalar types)
  ,TRANS      ///< Use the transposed operator 
  ,CONJTRANS  ///< Use the transposed operator with complex-conjugate clements (same as <tt>TRANS</tt> for real scalar types)
};

/** \brief Return a string name for a <tt>ETransp</tt> value.
 */
inline
const char* toString(ETransp transp)
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
 */
inline
ETransp real_trans(ETransp transp)
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

/** \brief Perform a not operation on an ETransp value
 */
inline 
ETransp not_trans( ETransp transp )
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
 */
inline
ETransp trans_trans( ETransp trans1, ETransp trans2 )
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

/** \brief Convert from <tt>ETransp</tt> to <tt>EConj</tt>.
 */
inline
EConj transToConj( ETransp trans )
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

/** \brief Convert from <tt>EConj</tt> to <tt>ETransp</tt> for forward apply.
 */
inline
ETransp applyConjToTrans( EConj conj ) {
  switch(conj) {
    case NONCONJ_ELE: return NOTRANS;
    case CONJ_ELE:    return CONJ;
    default: TEST_FOR_EXCEPT(true);
  }
  return NOTRANS; // Will never be called!
}

/** \brief Convert from <tt>EConj</tt> to <tt>ETransp</tt> for forward apply.
 */
inline
ETransp applyTransposeConjToTrans( EConj conj ) {
  switch(conj) {
    case NONCONJ_ELE: return TRANS;
    case CONJ_ELE:    return CONJTRANS;
    default: TEST_FOR_EXCEPT(true);
  }
  return NOTRANS; // Will never be called!
}

/** \brief Determines if a view is a direct view of data or a detached copy of
 * data.
 */
enum EViewType {
  VIEW_TYPE_DIRECT   ///< The view is a direct view of data and no copies are made.
  ,VIEW_TYPE_DETACHED  ///< The view is a detached copy of the data.
};

/** \brief Determine if data is unit stride or non-unit stride. */
enum EStrideType {
  STRIDE_TYPE_UNIT      ///< The stride between elements in an array is one.
  ,STRIDE_TYPE_NONUNIT   ///< The stride between elements in an array is greater than or equal to one.
};

//@}

namespace Exceptions {

/** \defgroup Thyra_Op_Vec_Exceptions_grp Basic Thyra exception types.
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
//@{

/// Thrown if any member functions are called before initialize() has been called.
class UnInitialized : public std::logic_error
{public: UnInitialized(const std::string& what_arg) : std::logic_error(what_arg) {}};

/// Thrown if vector spaces are incompatible
class IncompatibleVectorSpaces : public std::logic_error
{public:
  IncompatibleVectorSpaces(const std::string& what_arg) : std::logic_error(what_arg) {}
//	IncompatibleVectorSpaces(const IncompatibleVectorSpaces& ivs) : std::logic_error(ivs.what()) {}
};

/// Thrown if the argument <tt>M_trans</tt> is not supported,
class OpNotSupported : public std::logic_error
{public: OpNotSupported(const std::string& what_arg) : std::logic_error(what_arg) {}};

//@}

} // namespace Exceptions

// Fundamental ANA operator/vector interface classes

template<class Scalar>                                        class VectorSpaceFactoryBase;
template<class Scalar>                                        class VectorSpaceBase;
template<class RangeScalar, class DomainScalar = RangeScalar> class LinearOpBase;
template<class Scalar>                                        class MultiVectorBase;
template<class Scalar>                                        class VectorBase;

} // end namespace Thyra

#endif // THYRA_OPERATOR_VECTOR_TYPES_HPP
