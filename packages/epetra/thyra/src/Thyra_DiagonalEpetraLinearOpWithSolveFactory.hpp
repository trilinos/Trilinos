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

#ifndef __sun

#include "Thyra_LinearOpWithSolveFactoryBase.hpp"

namespace Thyra {

/** \brief. */
class DiagonalEpetraLinearOpWithSolveFactory : public LinearOpWithSolveFactoryBase<double> {
public:

  /** @name Overridden from LinearOpWithSolveFactoryBase */
  //@{

  /** \brief . */
  bool isCompatible( const LinearOpBase<double> &fwdOp ) const;

  /** \brief . */
  Teuchos::RefCountPtr<LinearOpWithSolveBase<double> > createOp() const;

  /** \brief . */
  void initializeOp(
    const Teuchos::RefCountPtr<const LinearOpBase<double> >    &fwdOp
    ,LinearOpWithSolveBase<double>                             *Op
    ) const;

  void uninitializeOp(
    LinearOpWithSolveBase<double>                        *Op
    ,Teuchos::RefCountPtr<const LinearOpBase<double> >   *fwdOp
    ,Teuchos::RefCountPtr<const LinearOpBase<double > >  *precOp
    ,EPreconditionerInputType                            *precOpType
    ) const;

  //@}

};

//@}

} // namespace Thyra

#endif // __sun

#endif // THYRA_DIAGONAL_EPETRA_LINEAR_OP_WITH_SOLVE_FACTORY_HPP
