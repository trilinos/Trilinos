/*
// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
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
*/

#ifndef THYRA_AMESOS_LINEAR_OP_WITH_SOLVE_HPP
#define THYRA_AMESOS_LINEAR_OP_WITH_SOLVE_HPP

#include "Thyra_SingleRhsLinearOpWithSolveBase.hpp"
#include "Thyra_EpetraLinearOp.hpp"

namespace Thyra {

/** \brief Concrete <tt>LinearOpWithSolveBase</tt> subclass that adapts any
 * <tt>Amesos_Solver</tt> object.
 *
 * ToDo: Finish documentation!
 */
class AmesosLinearOpWithSolve
  : virtual public LinearOpWithSolveBase<double>               // Public interface
  , virtual protected SingleRhsLinearOpWithSolveBase<double>   // Implementation detail
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief. */
  AmesosLinearOpWithSolve();

  /** @name Overridden from LinearOpBase */
  //@{
  /** \brief. */
  Teuchos::RefCountPtr< const VectorSpaceBase<double> > range() const;
  /** \brief. */
  Teuchos::RefCountPtr< const VectorSpaceBase<double> > domain() const;
  /** \brief. */
  Teuchos::RefCountPtr<const LinearOpBase<double> > clone() const;
  //@}

protected:

  /** @name Overridden from SingleScalarLinearOpBase */
  //@{
  /** \brief . */
  bool opSupported(ETransp M_trans) const;
  //@}

  /** @name Overridden from SingleRhsLinearOpBase */
  //@{
  /** \brief . */
  void apply(
    const ETransp                M_trans
    ,const VectorBase<double>    &x
    ,VectorBase<double>          *y
    ,const double                alpha
    ,const double                beta
    ) const;
  //@}

  /** @name Overridden from SingleScalarLinearOpWithSolveBase */
  //@{
  /** \brief . */
  bool solveSupportsTrans(ETransp M_trans) const;
  /** \brief . */
  bool solveSupportsSolveTolType(ETransp M_trans, ESolveTolType solveTolType) const;
  //@}

  /** @name Overridden from SingleRhsLinearOpWithSolveBase */
  //@{
  /** \brief . */
  SolveStatus<double> solve(
    const ETransp                         M_trans
    ,const VectorBase<double>             &b
    ,VectorBase<double>                   *x
    ,const SolveCriteria<double>          *solveCriteria
    ) const;
  //@}

private:

  EpetraLinearOp  epetraOp_;

  void assertInitialized() const;

};

} // namespace Thyra

#endif	// THYRA_AMESOS_LINEAR_OP_WITH_SOLVE_HPP
