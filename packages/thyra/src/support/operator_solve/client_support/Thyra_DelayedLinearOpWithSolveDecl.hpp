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


#ifndef THYRA_DELAYED_LINEAR_OP_WITH_SOLVE_DECL_HPP
#define THYRA_DELAYED_LINEAR_OP_WITH_SOLVE_DECL_HPP


#include "Thyra_SingleScalarLinearOpWithSolveBase.hpp"


namespace Thyra {


/** \brief Delayed linear solver construction LinearOpWithSolveBase decorator
 * class.
 *
 * This simple class takes a LinearOpWithSolveFactoryBase object and the
 * arguments for a call to LinearOpWithSolveBase::initializeOp() (or one of
 * its related functions) and then waits until a solve is called before it
 * constructs the preconditioner.
 *
 * This class should never be created directly by a client.  Instead, it
 * should be created by a DelayedLinearOpWithSolveFactory object.
 */
template <class Scalar>
class DelayedLinearOpWithSolve
  : virtual public LinearOpWithSolveBase<Scalar> // Interface
  , virtual protected SingleScalarLinearOpWithSolveBase<Scalar> // Implementation
{
public:

  /** \name Constructor/Initializers */
  //@{

  /** \brief . */
  DelayedLinearOpWithSolve();

  /** \brief . */
  void initialize(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    const RCP<const PreconditionerBase<Scalar> > &prec,
    const RCP<const LinearOpSourceBase<Scalar> > &approxFwdOpSrc,
    const ESupportSolveUse supportSolveUse,
    const RCP<LinearOpWithSolveFactoryBase<Scalar> > &lowsf
    );

  /** \brief . */
  RCP<const LinearOpSourceBase<Scalar> > getFwdOpSrc() const;

  /** \brief . */
  RCP<const PreconditionerBase<Scalar> > getPrec() const;
 
  /** \brief . */
  RCP<const LinearOpSourceBase<Scalar> > getApproxFwdOpSrc() const;

  /** \brief . */
  ESupportSolveUse getSupportSolveUse() const;

  //@}

  /** \name Overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief . */
  RCP< const VectorSpaceBase<Scalar> > range() const;
  /** \brief . */
  RCP< const VectorSpaceBase<Scalar> > domain() const;
  /** \brief . */
  RCP<const LinearOpBase<Scalar> > clone() const;

  //@}

protected:

  /** \brief Overridden from Teuchos::VerboseObjectBase */
  //@{

  /** \brief . */
  void informUpdatedVerbosityState() const;

  //@}

  /** @name Overridden from SingleScalarLinearOpBase  */
  //@{

  /** \brief . */
  virtual bool opSupported(EOpTransp M_trans) const;
  /** \brief . */
  virtual void apply(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &X,
    MultiVectorBase<Scalar> *Y,
    const Scalar alpha,
    const Scalar beta
    ) const;

  //@}

  /** @name Overridden from SingleScalarLinearOpWithSolveBase. */
  //@{

  /** \brief . */
  virtual bool solveSupportsTrans(EOpTransp M_trans) const;
  /** \brief . */
  virtual bool solveSupportsSolveMeasureType(
    EOpTransp M_trans, const SolveMeasureType& solveMeasureType
    ) const;
  /** \brief . */
  virtual void solve(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &B,
    MultiVectorBase<Scalar> *X,
    const int numBlocks,
    const BlockSolveCriteria<Scalar> blockSolveCriteria[],
    SolveStatus<Scalar> blockSolveStatus[]
    ) const;
  
  //@}

private:

  // //////////////////////
  // Private data members

  RCP<const LinearOpSourceBase<Scalar> > fwdOpSrc_;
  RCP<const PreconditionerBase<Scalar> > prec_;
  RCP<const LinearOpSourceBase<Scalar> > approxFwdOpSrc_;
  ESupportSolveUse supportSolveUse_;
  RCP<LinearOpWithSolveFactoryBase<Scalar> > lowsf_;

  RCP<const LinearOpBase<Scalar> > fwdOp_;

  mutable bool lows_is_valid_;
  mutable RCP<LinearOpWithSolveBase<Scalar> > lows_;

  // ///////////////////////////
  // Private member functions

  void updateSolver() const;

};


} // namespace Thyra


#endif // THYRA_DELAYED_LINEAR_OP_WITH_SOLVE_DECL_HPP
