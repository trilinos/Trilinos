// ***********************************************************************
// 
//               Thyra: Trilinos Solver Framework Core
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

#ifndef THYRA_EPETRA_OPERATOR_VIEW_EXTRACTOR_STD_HPP
#define THYRA_EPETRA_OPERATOR_VIEW_EXTRACTOR_STD_HPP

#include "Thyra_EpetraOperatorViewExtractorBase.hpp"

namespace Thyra {

/** \brief Standard strategy subclass for extracting an
 * <tt>Epetra_Operator</tt> view out of a <tt>Thyra::LinearOpBase<double></tt>
 * object by dynamic casting to the <tt>EpetraLinearOpBase</tt> interface.
 *
 * ToDo: Finish documentation!
 */
class EpetraOperatorViewExtractorStd : virtual public EpetraOperatorViewExtractorBase
{
public:

  /** \name Overridden from EpetraOperatorViewExtractorBase. */
  //@{

  /** \brief . */
  bool isCompatible( const LinearOpBase<double> &fwdOp ) const;
  /** \brief . */
  void getEpetraOpView(
    const Teuchos::RCP<LinearOpBase<double> >   &fwdOp
    ,Teuchos::RCP<Epetra_Operator>              *epetraOp
    ,EOpTransp                                            *epetraOpTransp
    ,EApplyEpetraOpAs                                   *epetraOpApplyAs
    ,EAdjointEpetraOp                                   *epetraOpAdjointSupport
    ,double                                             *epetraOpScalar
    ) const;
  /** \brief . */
  void getEpetraOpView(
    const Teuchos::RCP<const LinearOpBase<double> >   &fwdOp
    ,Teuchos::RCP<const Epetra_Operator>              *epetraOp
    ,EOpTransp                                                  *epetraOpTransp
    ,EApplyEpetraOpAs                                         *epetraOpApplyAs
    ,EAdjointEpetraOp                                         *epetraOpAdjointSupport
    ,double                                                   *epetraOpScalar
    ) const;

  //@}

};

} // namespace Thyra

#endif // THYRA_EPETRA_OPERATOR_VIEW_EXTRACTOR_STD_HPP
