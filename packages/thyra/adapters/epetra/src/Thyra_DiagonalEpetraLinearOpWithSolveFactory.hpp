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

#ifndef THYRA_DIAGONAL_EPETRA_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
#define THYRA_DIAGONAL_EPETRA_LINEAR_OP_WITH_SOLVE_FACTORY_HPP


#include "Thyra_LinearOpWithSolveFactoryBase.hpp"


namespace Thyra {


/** \brief Create a DefaultDiagonalLinearOpWithSolve out of a diagonal
 * Epetra_RowMatrix object.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
class DiagonalEpetraLinearOpWithSolveFactory : public LinearOpWithSolveFactoryBase<double> {
public:

  /** @name Overridden from LinearOpWithSolveFactoryBase */
  //@{

  /** \brief . */
  bool isCompatible( const LinearOpSourceBase<double> &fwdOpSrc ) const;

  /** \brief . */
  Teuchos::RCP<LinearOpWithSolveBase<double> > createOp() const;

  /** \brief . */
  void initializeOp(
    const Teuchos::RCP<const LinearOpSourceBase<double> > &fwdOpSrc,
    LinearOpWithSolveBase<double> *Op,
    const ESupportSolveUse supportSolveUse
    ) const;

  void uninitializeOp(
    LinearOpWithSolveBase<double> *Op,
    Teuchos::RCP<const LinearOpSourceBase<double> > *fwdOpSrc,
    Teuchos::RCP<const PreconditionerBase<double> > *prec,
    Teuchos::RCP<const LinearOpSourceBase<double> > *approxFwdOpSrc,
    ESupportSolveUse *supportSolveUse
    ) const;

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

};


} // namespace Thyra


#endif // THYRA_DIAGONAL_EPETRA_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
