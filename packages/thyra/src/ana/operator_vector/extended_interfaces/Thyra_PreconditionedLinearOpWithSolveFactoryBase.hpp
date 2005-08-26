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

#ifndef THYRA_PRECONDITIONED_LINEAR_OP_WITH_SOLVE_FACTORY_BASE_DECL_HPP
#define THYRA_PRECONDITIONED_LINEAR_OP_WITH_SOLVE_FACTORY_BASE_DECL_HPP

#include "Thyra_LinearOpWithSolveFactoryBase.hpp"

namespace Thyra {

/** \brief Factory interface for creating <tt>LinearOpWithSolveBase</tt>
 * objects from <tt>LinearOpBase</tt> objects for the operator and an optional
 * preconditioner.
 *
 * ToDo: Finish documation!
 *
 * \ingroup Thyra_Op_Vec_Interoperability_Extended_Interfaces_grp
 */
template <class RangeScalar, class DomainScalar = RangeScalar>
class PreconditionedLinearOpWithSolveFactoryBase
  : virtual public LinearOpWithSolveFactoryBase<RangeScalar,DomainScalar> 
{
public:

  /** @name Pure virtual functions that must be overridden in subclasses */
  //@{

  /** \brief Initialize a pre-created <tt>LinearOpWithSolveBase</tt> object
   * given a "compatible" <tt>LinearOpBase</tt> object and an optional
   * preconditioner.
   *
   * ToDo: Finish documation!
   */
  virtual void initializePreconditionedOp(
    const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >     &fwdOp
    ,const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &precOp
    ,const EPreconditionerInputType                                               precOpType
    ,LinearOpWithSolveBase<RangeScalar,DomainScalar>                              *Op
    ) const = 0;

  /** \brief Uninitialize a <tt>LinearOpWithSolveBase</tt> object and return its
   * remembered forward linear operator and preconditioner objects.
   *
   * ToDo: Finish documation!
   */
  virtual void uninitializePreconditionedOp(
    LinearOpWithSolveBase<RangeScalar,DomainScalar>                       *Op
    ,Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >  *fwdOp       = NULL
    ,Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >  *precOp      = NULL
    ,EPreconditionerInputType                                             *precOpType  = NULL
    ) const = 0;
  
  //@}

};

//@}

} // namespace Thyra

#endif // THYRA_PRECONDITIONED_LINEAR_OP_WITH_SOLVE_FACTORY_BASE_DECL_HPP
