// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_SERIAL_VECTOR_SPACE_CONVERTED_STD_DECL_HPP
#define THYRA_SERIAL_VECTOR_SPACE_CONVERTED_STD_DECL_HPP

#include "Thyra_SerialVectorSpaceConverterBase.hpp"


namespace Thyra {


/** \brief Concrete subclass for a converter subclass for converting serial
 * multi-vectors and vectors.
 *
 * While this concrete subclass creates concrete vector spaces of type
 * <tt>DefaultSerialVectorSpace</tt>, it should be usable with any serial vector
 * space type and therefore this subclass is more general then it may appear
 * at first.
 *
 * \ingroup Thyra_Op_Vec_adapters_Spmd_concrete_std_grp
 */
template<class ScalarFrom, class ScalarTo>
class DefaultSerialVectorSpaceConverter : virtual public SerialVectorSpaceConverterBase<ScalarFrom,ScalarTo> {
public:

  /** @name Overridden from VectorSpaceConverterBase */
  //@{

  /** \brief . */
  virtual
  Teuchos::RCP<const VectorSpaceBase<ScalarTo> >
  createVectorSpaceTo(
    const VectorSpaceBase<ScalarFrom> &vecSpc
    ) const;

  /** \brief . */
  virtual
  Teuchos::RCP<const VectorSpaceBase<ScalarFrom> >
  createVectorSpaceFrom(
    const VectorSpaceBase<ScalarTo> &vecSpc
    ) const;

  //@}
  
};


// Implementation


template<class ScalarFrom, class ScalarTo>
Teuchos::RCP<const VectorSpaceBase<ScalarTo> >
DefaultSerialVectorSpaceConverter<ScalarFrom,ScalarTo>::createVectorSpaceTo(
  const VectorSpaceBase<ScalarFrom>  &vecSpc
  ) const
{
  return defaultSpmdVectorSpace<ScalarTo>(vecSpc.dim());
}


template<class ScalarFrom, class ScalarTo>
Teuchos::RCP<const VectorSpaceBase<ScalarFrom> >
DefaultSerialVectorSpaceConverter<ScalarFrom,ScalarTo>::createVectorSpaceFrom(
  const VectorSpaceBase<ScalarTo>  &vecSpc
  ) const
{
  return defaultSpmdVectorSpace<ScalarFrom>(vecSpc.dim());
}


} // namespace Thyra


#endif // THYRA_SERIAL_VECTOR_SPACE_CONVERTED_STD_DECL_HPP
