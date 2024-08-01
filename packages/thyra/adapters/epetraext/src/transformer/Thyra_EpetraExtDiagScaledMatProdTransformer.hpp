// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_EPETRAEXT_DIAG_SCALED_MAT_PROD_TRANSFORMER_HPP
#define THYRA_EPETRAEXT_DIAG_SCALED_MAT_PROD_TRANSFORMER_HPP


#include "Thyra_LinearOpTransformerBase.hpp"


namespace Thyra {


/** \brief Transformer subclass for diagonally scaling and multiplying
 * Epetra/Thyra operators.
 *
 * \ingroup EpetraExt_Thyra_Op_Vec_adapters_grp
 */
class EpetraExtDiagScaledMatProdTransformer : public LinearOpTransformerBase<double>
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
 * \relates EpetraExtDiagScaledMatProdTransformer
 */
inline
RCP<EpetraExtDiagScaledMatProdTransformer>
epetraExtDiagScaledMatProdTransformer()
{
  return Teuchos::rcp(new EpetraExtDiagScaledMatProdTransformer());
}


} // namespace Thyra


#endif	// THYRA_EPETRAEXT_DIAG_SCALED_MAT_PROD_TRANSFORMER_HPP

#if defined(Thyra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ThyraEpetraExtAdapters package is deprecated"
#endif
#endif

