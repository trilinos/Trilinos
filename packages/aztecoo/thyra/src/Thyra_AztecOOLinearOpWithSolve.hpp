/*@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
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
//@HEADER
*/

#ifndef THYRA_AZTECOO_LINEAR_OP_WITH_SOLVE_HPP
#define THYRA_AZTECOO_LINEAR_OP_WITH_SOLVE_HPP

#include "Thyra_SingleRhsLinearOpWithSolveBase.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "AztecOO.h"

namespace Thyra {

/** \brief Concrete <tt>LinearOpWithSolveBase</tt> subclass in terms of
 * <tt>AztecOO</tt>.
 *
 * This subclass is designed to be very flexible and handle a number of
 * different use cases.  It supports forward and optionally adjoint
 * (transpose) solves.  I can support inexact solves based on a residual norm
 * tolerance or just allow for a default (i.e. tight) linear solve tolerance.
 * Currently, this subclass does not support inexact solves by specifying a
 * tolerance on the estimate of the solution error and it is unlikely that
 * this subclass with ever support this mode.
 *
 * This subclass is not designed to be used directly by users but instead by
 * subclasses of <tt>LinearOpWithSolveFactoryBase</tt>.  One standard
 * implementation that is fairly flexible (and will be make more flexible in
 * the future) is <tt>AztecOOLinearOpWithSolveFactory</tt>.
 *
 * This subclass allows for user-defined preconditioners or for built-in aztec
 * preconditioners.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup AztecOO_Thyra_adapters_grp
 */
class AztecOOLinearOpWithSolve
  : virtual public LinearOpWithSolveBase<double>               // Public interface
  , virtual protected SingleRhsLinearOpWithSolveBase<double>   // Implementation detail
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

	/** Construct uninitialized but with default option values.
	 *
	 * Note, these defaults where taken from
	 * NOX::EpetraNew::LinearSystemAztecOO::applyJacobianInverse(...) on
	 * 2005/08/15.
	 */
 	AztecOOLinearOpWithSolve(
	 	const int       fwdDefaultMaxIterations       = 400
		,const double   fwdDefaultTol                 = 1e-6
	 	,const int      adjDefaultMaxIterations       = 400
		,const double   adjDefaultTol                 = 1e-6
		);

	/** \brief The default maximum number of iterations for forward solves. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, fwdDefaultMaxIterations )
  /** \brief The default solution tolerance on the residual for forward solves. */
 	STANDARD_MEMBER_COMPOSITION_MEMBERS( double, fwdDefaultTol )
	/** \brief The default maximum number of iterations for adjoint solves. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, adjDefaultMaxIterations )
  /** \brief The default solution tolerance on the residual for adjoint solves. */
 	STANDARD_MEMBER_COMPOSITION_MEMBERS( double, adjDefaultTol )

  /** \brief Sets up this object.
   *
   * \param  fwdOp   [in] The forward operator object that defines this objects
   *                 <tt>LinearOpBase</tt> interface.  This also should be the
   *                 exact same object that is passed in through a <tt>LinearOpWithSolveFactoryBase</tt>
   *                 interface.
   * \param  precOp  [in] The original abstract preconditioner object that was passed through the
   *                 <tt>PreconditionedLinearOpWithSolveFactoryBase</tt> interface.  This object
   *                 is not used for anything and can be set as <tt>precOp.get()==NULL</tt>.
   * \param  precOpType
   *                 [in] The original argument passed through the 
   *                 <tt>PreconditionedLinearOpWithSolveFactoryBase</tt> interface.  This value
   *                 is not used for anything can be set to any value the client desires.
   * \param  aztecFwdSolver
   *                 [in] The <tt>AztecOO</tt> object used to perform forward solves.  This object must be
   *                 be ready to call <tt>aztecFwdSolver->SetRHS()</tt> and <tt>aztecFwdSolver->SetLHS()</tt>
   *                 and then call <tt>aztecFwdSolver->Solve()</tt>.
   * \param  allowInexactFwdSolve
   *                 [in] Determines if <tt>this->solveSupportsSolveTolType(NOTRANS,SOLVE_TOL_REL_RESIDUAL_NORM)</tt>
   *                 returns true or not.  With the current design, an inexact forward solve can not be supported
   *                 if there is left scaling or a left preconditioner aggregated with <tt>*aztecFwdOp</tt>.
   * \param  aztecAdjSolver
   *                 [in] The <tt>AztecOO</tt> object used to perform adjoint solves.  This object must be
   *                 be ready to call <tt>aztecAdjSolver->SetRHS()</tt> and <tt>aztecAdjSolver->SetLHS()</tt>
   *                 and then call <tt>aztecAdjSolver->Solve()</tt>.
   * \param  allowInexactAdjSolve
   *                 [in] Determines if <tt>this->solveSupportsSolveTolType(TRANS,SOLVE_TOL_REL_RESIDUAL_NORM)</tt>
   *                 returns true or not.  With the current design, an inexact forward solve can not be supported
   *                 if there is left scaling or a left preconditioner aggregated with <tt>*aztecFwdOp</tt>.
   * \param  linearSystemTransformer
   *                 [in] This is a transformation object that is called to pre-preprocess the linear problem
   *                 before a forward and adjoint linear solver and post-process the linear problem after
   *                 forward and adjoint linear solve.  This abstract object is used to deal with scaling
   *                 and aggregated preconditioners.  It is what makes this implementation fairly flexible.
	 *
   * <b>Preconditions:</b><ul>
   * <li><tt>fwdOp.get()!=NULL</tt>
   * <li><tt>fwdFwdSolver.get()!=NULL</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->range() == fwdOp->range()</tt>
   * <li><tt>this->domain() == fwdOp->domain()</tt>
   * <li><tt>this->opSupports(M_trans) == opSupports(*fwdOp,M_trans)</tt>
   * <li><tt>this->solveSupportsTrans(M_trans) == (aztecAdjSolver.get()!=NULL)</tt>
   * <li><tt>this->solveSupportsSolveTolType([NOTRANS,CONJ],SOLVE_TOL_REL_RESIDUAL_NORM) == allowInexactFwdSolve</tt>
   * <li><tt>this->solveSupportsSolveTolType([NOTRANS,CONJ],SOLVE_TOL_REL_ERR_NORM) == false</tt>
   * <li><tt>this->solveSupportsSolveTolType([TRANS,CONJTRANS],SOLVE_TOL_REL_RESIDUAL_NORM) == (aztecAdjSolver.get()!=NULL&&allowInexactAdjSolve)</tt>
   * <li><tt>this->solveSupportsSolveTolType([TRANS,CONJTRANS],SOLVE_TOL_REL_ERR_NORM) == false</tt>
   * </ul>
   *
   * ToDo: Finish documentation!
	 */
  void initialize(
    const Teuchos::RefCountPtr<const LinearOpBase<double> >                 &fwdOp
    ,const Teuchos::RefCountPtr<const LinearOpBase<double> >                &precOp
    ,const EPreconditionerInputType                                         precOpType
    ,const Teuchos::RefCountPtr<AztecOO>                                    &aztecFwdSolver
    ,const bool                                                             allowInexactFwdSolve       = false
    ,const Teuchos::RefCountPtr<AztecOO>                                    &aztecAdjSolver            = Teuchos::null
    ,const bool                                                             allowInexactAdjSolve       = false
    ,const double                                                           aztecSolverScalar          = 1.0
    );

  /** \brief Extract the forward <tt>LinearOpBase<double></tt> object so that it can be modified.
   * 
   * <b>Postconditions:</b><ul>
   * <li><tt>return.get()</tt> is the same as <tt>this->get_fwdOp().get()</tt> before call.
   * <li><tt><tt>this->get_fwdOp().get()==NULL</tt>
   * </ul>
   */
  Teuchos::RefCountPtr<const LinearOpBase<double> > extract_fwdOp();

  /** \brief Extract the original preconditioner.
   * 
   * <b>Postconditions:</b><ul>
   * <li><tt>return.get()</tt> is the same as <tt>this->get_precOp().get()</tt> before call.
   * <li><tt><tt>this->get_precOp().get()==NULL</tt>
   * </ul>
   */
  Teuchos::RefCountPtr<const LinearOpBase<double> > extract_precOp();

  /** \brief Extract the original preconditioner type.
   */
  EPreconditionerInputType extract_precOpType();
	
	/** \brief Uninitialize. */
	void uninitialize(
    Teuchos::RefCountPtr<const LinearOpBase<double> >                 *fwdOp                     = NULL
    ,Teuchos::RefCountPtr<const LinearOpBase<double> >                *precOp                    = NULL
    ,EPreconditionerInputType                                         *precOpType                = NULL
    ,Teuchos::RefCountPtr<AztecOO>                                    *aztecFwdSolver            = NULL
    ,bool                                                             *allowInexactFwdSolve      = NULL
    ,Teuchos::RefCountPtr<AztecOO>                                    *aztecAdjSolver            = NULL
    ,bool                                                             *allowInexactAdjSolve      = NULL
    ,double                                                           *aztecSolverScalar         = NULL
    );

  //@}

  /** @name Overridden from LinearOpBase */
  //@{
  /** \brief. */
  Teuchos::RefCountPtr< const VectorSpaceBase<double> > range() const;
  /** \brief. */
  Teuchos::RefCountPtr< const VectorSpaceBase<double> > domain() const;
  /** \brief. */
  Teuchos::RefCountPtr<const LinearOpBase<double> > clone() const;
  //@}

  /** @name Overridden from Teuchos::Describable */
  //@{
  /** \brief . */
  std::string description() const;
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
  /** \brief . */
  int defaultSolveMaxIterations(ETransp M_trans, ESolveTolType solveTolType) const;
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
  
  Teuchos::RefCountPtr<const LinearOpBase<double> >                fwdOp_;
  Teuchos::RefCountPtr<const LinearOpBase<double> >                precOp_;
  EPreconditionerInputType                                         precOpType_;
  Teuchos::RefCountPtr<AztecOO>                                    aztecFwdSolver_;
  bool                                                             allowInexactFwdSolve_;
  Teuchos::RefCountPtr<AztecOO>                                    aztecAdjSolver_;
  bool                                                             allowInexactAdjSolve_;
  double                                                           aztecSolverScalar_;
                                                     
  void assertInitialized() const;
  
};

} // namespace Thyra

#endif	// THYRA_AZTECOO_LINEAR_OP_WITH_SOLVE_HPP
