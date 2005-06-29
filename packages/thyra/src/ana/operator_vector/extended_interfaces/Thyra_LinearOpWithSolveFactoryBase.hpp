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

#ifndef THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_BASE_DECL_HPP
#define THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_BASE_DECL_HPP

#include "Thyra_LinearOpWithSolveBaseDecl.hpp"

namespace Thyra {

/** \brief Factory interface for creating <tt>LinearOpWithSolveBase</tt> objects. */
template <class RangeScalar, class DomainScalar = RangeScalar>
class LinearOpWithSolveFactoryBase {
public:

  /** \brief . */
  virtual ~LinearOpWithSolveFactoryBase() {}

  /** @name Pure virtual functions that must be overridden in subclasses */
  //@{

  /** \brief Create an (uninitialized) <tt>LinearOpWithSolveBase</tt> object to be
   * initialized later.
   */
  virtual Teuchos::RefCountPtr<LinearOpWithSolveBase<RangeScalar,DomainScalar> > createOp() const = 0;

  /** \brief Initalize a precreated <tt>LinearOpWithSolveBase</tt> object
   * given a "compatible" <tt>LinearOpBase</tt>.
   *
   * The assumption is that the output <tt>Op</tt> object may keep a memory to
   * the input <tt>fwdOp</tt> object but the factory itself may not.
   */
  virtual void initializeOp(
    const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &fwdOp
    ,LinearOpWithSolveBase<RangeScalar,DomainScalar>                             *Op
    ) const = 0;

  //@}

};

//@}

} // namespace Thyra

#endif // THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_BASE_DECL_HPP
