// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_SERIAL_VECTOR_SPACE_CONVERTED_BASE_HPP
#define THYRA_SERIAL_VECTOR_SPACE_CONVERTED_BASE_HPP

#include "Thyra_VectorSpaceConverterBase.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"


namespace Thyra {


/** \brief Node base class for converting serial multi-vectors (and vectors)
 * from one scalar type to another.
 *
 * This node base class defines the function <tt>convert()</tt> for all serial
 * vectors.  A concrete subclass is created by deriving from this interface
 * and then defining the function <tt>createVectorSpace()</tt>.
 *
 * \ingroup Thyra_Op_Vec_adapters_Spmd_support_grp
 */
template<class ScalarFrom, class ScalarTo>
class SerialVectorSpaceConverterBase
  : virtual public VectorSpaceConverterBase<ScalarFrom,ScalarTo>
{
public:

  /** @name Overridden from VectorSpaceConverterBase */
  //@{

  /** \brief . */
  virtual void convert(
    const MultiVectorBase<ScalarFrom> &mv_from,
    MultiVectorBase<ScalarTo> *mv_to
    ) const;

  //@}
  
};


// Implementation


template<class ScalarFrom, class ScalarTo>
void SerialVectorSpaceConverterBase<ScalarFrom,ScalarTo>::convert(
  const MultiVectorBase<ScalarFrom> &mv_from,
  MultiVectorBase<ScalarTo> *mv_to
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(mv_to==NULL);
#endif
  ConstDetachedMultiVectorView<ScalarFrom> emv_from(mv_from);
  DetachedMultiVectorView<ScalarTo> emv_to(*mv_to);
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(emv_from.subDim() != emv_to.subDim());
  TEUCHOS_TEST_FOR_EXCEPT(emv_from.numSubCols() != emv_to.numSubCols());
#endif
  for( Ordinal j = 0; j < emv_from.numSubCols(); ++j ) {
    for( Ordinal i = 0; i < emv_from.subDim(); ++i ) {
      emv_to(i,j) = emv_from(i,j); // ToDo: Make this faster using optimized copy functions?
    }
  }
}


} // namespace Thyra


#endif // THYRA_SERIAL_VECTOR_SPACE_CONVERTED_BASE_HPP
