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

#ifndef TSF_CORE_OPERATOR_VECTOR_TYPES_HPP
#define TSF_CORE_OPERATOR_VECTOR_TYPES_HPP

#include "RTOpPack_Types.hpp"
#include "Thyra_Range1D.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace Thyra {

/** \defgroup Thyra_Op_Vec_BasicTypes_grp Basic Thyra types.
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
//@{

//
// Basic types
//

/// Type for the dimension of a vector space
typedef RTOp_index_type  Index;

/// Type for a range of indices
typedef RangePack::Range1D   Range1D;

/** \brief Enumeration for determining how an operator is applied.
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

//
// Thyra
//

// Core abstract interface classes

template<class Scalar> class VectorSpaceFactoryBase;
template<class Scalar> class VectorSpaceBase;
template<class Scalar> class VectorBase;
template<class Scalar> class LinearOpBase;
template<class Scalar> class MultiVectorBase;

// Basic node support subclasses and interfaces

template<class Scalar> class ScalarProdBase;
template<class Scalar> class ScalarProdVectorSpaceBase;
template<class Scalar> class EuclideanLinearOpBase;
template<class Scalar> class SerialVectorSpaceBase;
template<class Scalar> class SerialVectorBase;

// Basic concrete support subclasses

template<class Scalar> class EuclideanScalarProd;
template<class Scalar> class LinearOpScalarProd;
template<class Scalar> class SerialVectorSpaceFactoryStd;
template<class Scalar> class SerialVectorSpaceStd;
template<class Scalar> class SerialVectorStd;
template<class Scalar> class SerialMultiVectorStd;
template<class Scalar> class MultiVectorCols;
template<class Scalar> class VectorMultiVector;

} // end namespace Thyra

#endif // TSF_CORE_OPERATOR_VECTOR_TYPES_HPP
