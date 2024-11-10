// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_EPETRAEXT_ADD_TRANSFORMER_HPP
#define THYRA_EPETRAEXT_ADD_TRANSFORMER_HPP

#include "Thyra_LinearOpTransformerBase.hpp"


namespace Thyra {


/** \brief Transformer subclass for adding Epetra/Thyra operators using
 * EpetraExt::MatrixMatrix.
 *
 * \ingroup EpetraExt_Thyra_Op_Vec_adapters_grp
 */
class EpetraExtAddTransformer : public LinearOpTransformerBase<double>
{
public:

  /** \name Overridden from LinearOpTransformerBase. */
  //@{

  /** \brief . */
  virtual bool isCompatible(const LinearOpBase<double> &op_in) const;

  /** \brief . */
  virtual RCP<LinearOpBase<double> > createOutputOp() const;

  /** \brief . */
  virtual void transform(
    const LinearOpBase<double> &op_in,
    const Ptr<LinearOpBase<double> > &op_inout
    ) const;

  //@}

private:
  
};


/** \brief Nonmember constructor.
 *
 * \relates EpetraExtAddTransformer
 */
inline
RCP<EpetraExtAddTransformer>
epetraExtAddTransformer()
{
  return Teuchos::rcp(new EpetraExtAddTransformer());
}


} // namespace Thyra


#endif	// THYRA_EPETRAEXT_ADD_TRANSFORMER_HPP

#if defined(Thyra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ThyraEpetraExtAdapters package is deprecated"
#endif
#endif

