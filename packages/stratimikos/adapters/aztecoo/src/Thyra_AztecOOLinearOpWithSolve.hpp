// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_AZTECOO_LINEAR_OP_WITH_SOLVE_HPP
#define THYRA_AZTECOO_LINEAR_OP_WITH_SOLVE_HPP

#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_LinearOpSourceBase.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_PreconditionerBase.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "AztecOO.h"


namespace Thyra {


/** \brief Concrete <tt>LinearOpWithSolveBase</tt> subclass implemented using
 * <tt>AztecOO</tt>.
 *
 * This subclass is designed to be very flexible and handle a number of
 * different use cases.  It supports forward and optionally adjoint
 * (transpose) solves.  I can support inexact solves based on a relative
 * residual norm tolerance or just allow for a default (i.e. tight) linear
 * solve tolerance.
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
class AztecOOLinearOpWithSolve : virtual public LinearOpWithSolveBase<double>
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
     const int fwdDefaultMaxIterations = 400,
     const double fwdDefaultTol = 1e-6,
     const int adjDefaultMaxIterations = 400,
     const double adjDefaultTol = 1e-6,
     const bool outputEveryRhs = false
     );

  /** \brief The default maximum number of iterations for forward solves. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, fwdDefaultMaxIterations );
  /** \brief The default solution tolerance on the residual for forward solves. */
   STANDARD_MEMBER_COMPOSITION_MEMBERS( double, fwdDefaultTol );
  /** \brief The default maximum number of iterations for adjoint solves. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, adjDefaultMaxIterations );
  /** \brief The default solution tolerance on the residual for adjoint solves. */
   STANDARD_MEMBER_COMPOSITION_MEMBERS( double, adjDefaultTol );
  /** \brief Determine if output for every RHS will be printed or not. */
   STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, outputEveryRhs );

  /** \brief Sets up this object.
   *
   * \param fwdOp [in] The forward operator object that defines this objects
   * <tt>LinearOpBase</tt> interface.  interface.
   *
   * \param fwdOpSrc [in] The source for the forward operator object
   * <tt>fwdOp</tt>.  This also should be the exact same object that is passed
   * in through a <tt>LinearOpWithSolveFactoryBase</tt> interface.
   *
   * \param prec [in] The original abstract preconditioner object that was
   * passed through the <tt>LinearOpWithSolveFactoryBase</tt> interface.  This
   * object is not used for anything and can be set as
   * <tt>prec==Teuchos::null</tt>.
   *
   * \param isExternalPrec [in] True if the precondition was created
   * externally from the <tt>LinearOpWithSolveFactoryBase</tt> object, false
   * otherwise.
   *
   * \param approxFwdOpSrc [in] The source for the original abstract
   * approximate forward operator object that was passed through the
   * <tt>LinearOpWithSolveFactoryBase</tt> interface.  This object is not used
   * for anything and can be set as <tt>approxFwdOpSrc==Teuchos::null</tt>.
   *
   * \param aztecFwdSolver [in] The <tt>AztecOO</tt> object used to perform
   * forward solves.  This object must be be ready to call
   * <tt>aztecFwdSolver->SetRHS()</tt> and <tt>aztecFwdSolver->SetLHS()</tt>
   * and then call <tt>aztecFwdSolver->Solve()</tt>.
   *
   * \param allowInexactFwdSolve [in] Determines if
   * <tt>this->solveSupportsSolveTolType(NOTRANS,SOLVE_TOL_REL_RESIDUAL_NORM)</tt>
   * returns true or not.  With the current design, an inexact forward solve
   * can not be supported if there is left scaling or a left preconditioner
   * aggregated with <tt>*aztecFwdOp</tt>.
   *
   * \param aztecAdjSolver [in] The <tt>AztecOO</tt> object used to perform
   * adjoint solves.  This object must be be ready to call
   * <tt>aztecAdjSolver->SetRHS()</tt> and <tt>aztecAdjSolver->SetLHS()</tt>
   * and then call <tt>aztecAdjSolver->Solve()</tt>.
   *
   * \param allowInexactAdjSolve [in] Determines if
   * <tt>this->solveSupportsSolveTolType(TRANS,SOLVE_TOL_REL_RESIDUAL_NORM)</tt>
   * returns true or not.  With the current design, an inexact forward solve
   * can not be supported if there is left scaling or a left preconditioner
   * aggregated with <tt>*aztecFwdOp</tt>.
   *
   * \param linearSystemTransformer [in] This is a transformation object that
   * is called to pre-preprocess the linear problem before a forward and
   * adjoint linear solver and post-process the linear problem after forward
   * and adjoint linear solve.  This abstract object is used to deal with
   * scaling and aggregated preconditioners.  It is what makes this
   * implementation fairly flexible.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>fwdOp.get()!=NULL</tt>
   * <li><tt>fwdOpSrc.get()!=NULL</tt>
   * <li><tt>fwdFwdSolver.get()!=NULL</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->range() == fwdOp->range()</tt>
   * <li><tt>this->domain() == fwdOp->domain()</tt>
   * <li><tt>this->opSupports(M_trans) == opSupports(*fwdOp,M_trans)</tt>
   * <li><tt>this->solveSupportsTrans(M_trans) == (aztecAdjSolver.get()!=NULL)</tt>
   * <li><tt>this->solveSupportsSolveTolType([NOTRANS,CONJ],
             SolveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MESURE_NORM_RHS))
   *         == allowInexactFwdSolve</tt>
   * <li><tt>this->solveSupportsSolveTolType([TRANS,CONJTRANS],
             SolveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MESURE_NORM_RHS))
   *         == (aztecAdjSolver.get()!=NULL&&allowInexactAdjSolve)</tt>
   * </ul>
   *
   * ToDo: Finish documentation!
   */
  void initialize(
    const RCP<const LinearOpBase<double> > &fwdOp,
    const RCP<const LinearOpSourceBase<double> > &fwdOpSrc,
    const RCP<const PreconditionerBase<double> > &prec,
    const bool isExternalPrec,
    const RCP<const LinearOpSourceBase<double> > &approxFwdOpSrc,
    const RCP<AztecOO> &aztecFwdSolver,
    const bool allowInexactFwdSolve = false,
    const RCP<AztecOO> &aztecAdjSolver = Teuchos::null,
    const bool allowInexactAdjSolve = false,
    const double aztecSolverScalar = 1.0
    );
  
  /** \brief Extract the forward <tt>LinearOpBase<double></tt> object so that
   * it can be modified.
   */
  RCP<const LinearOpSourceBase<double> > extract_fwdOpSrc();

  /** \brief Extract the preconditioner.
   */
  RCP<const PreconditionerBase<double> > extract_prec();

  /** \brief Determine if the preconditioner was external or not.
   */
  bool isExternalPrec() const;


  /** \brief Extract the approximate forward <tt>LinearOpBase<double></tt>
   * object used to build the preconditioner.
   */
  RCP<const LinearOpSourceBase<double> > extract_approxFwdOpSrc();
  
  /** \brief Uninitialize. */
  void uninitialize(
    RCP<const LinearOpBase<double> > *fwdOp = NULL,
    RCP<const LinearOpSourceBase<double> > *fwdOpSrc = NULL,
    RCP<const PreconditionerBase<double> > *prec = NULL,
    bool *isExternalPrec = NULL,
    RCP<const LinearOpSourceBase<double> > *approxFwdOpSrc = NULL,
    RCP<AztecOO> *aztecFwdSolver = NULL,
    bool *allowInexactFwdSolve = NULL,
    RCP<AztecOO> *aztecAdjSolver = NULL,
    bool *allowInexactAdjSolve = NULL,
    double *aztecSolverScalar = NULL
    );

  //@}

  /** @name Overridden from LinearOpBase */
  //@{
  /** \brief. */
  RCP< const VectorSpaceBase<double> > range() const;
  /** \brief. */
  RCP< const VectorSpaceBase<double> > domain() const;
  /** \brief. */
  RCP<const LinearOpBase<double> > clone() const;
  //@}

  /** @name Overridden from Teuchos::Describable */
  //@{
  /** \brief . */
  std::string description() const;
  /** \brief . */
  void describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;
  //@}

protected:

  /** @name Overridden from LinearOpBase  */
  //@{
  /** \brief . */
  virtual bool opSupportedImpl(EOpTransp M_trans) const;
  /** \brief . */
  virtual void applyImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<double> &X,
    const Ptr<MultiVectorBase<double> > &Y,
    const double alpha,
    const double beta
    ) const;
  //@}

  /** @name Overridden from LinearOpWithSolveBase. */
  //@{
  /** \brief . */
  virtual bool solveSupportsImpl(EOpTransp M_trans) const;
  /** \brief . */
  virtual bool solveSupportsSolveMeasureTypeImpl(
    EOpTransp M_trans, const SolveMeasureType& solveMeasureType
    ) const;
  /** \brief . */
  SolveStatus<double> solveImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<double> &B,
    const Ptr<MultiVectorBase<double> > &X,
    const Ptr<const SolveCriteria<double> > solveCriteria
    ) const;
  //@}
  
private:
  
  RCP<const LinearOpBase<double> > fwdOp_;
  RCP<const LinearOpSourceBase<double> > fwdOpSrc_;
  RCP<const PreconditionerBase<double> > prec_;
  bool isExternalPrec_;
  RCP<const LinearOpSourceBase<double> > approxFwdOpSrc_;
  RCP<AztecOO> aztecFwdSolver_;
  bool allowInexactFwdSolve_;
  RCP<AztecOO> aztecAdjSolver_;
  bool allowInexactAdjSolve_;
  double aztecSolverScalar_;
  
  void assertInitialized() const;

};


} // namespace Thyra

#endif	// THYRA_AZTECOO_LINEAR_OP_WITH_SOLVE_HPP
