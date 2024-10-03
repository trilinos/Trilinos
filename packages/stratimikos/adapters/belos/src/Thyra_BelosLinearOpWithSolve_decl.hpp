// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_BELOS_LINEAR_OP_WITH_SOLVE_DECL_HPP
#define THYRA_BELOS_LINEAR_OP_WITH_SOLVE_DECL_HPP

#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_LinearOpSourceBase.hpp"
#include "BelosSolverManager.hpp"
#include "BelosThyraAdapter.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"


namespace Thyra {


/** \brief Concrete <tt>LinearOpWithSolveBase</tt> subclass in terms of
 * <tt>Belos</tt>.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Belos_Thyra_adapters_grp
 */
template<class Scalar>
class BelosLinearOpWithSolve : virtual public LinearOpWithSolveBase<Scalar>
{
public:

  /** @name Public typedefs */
  //@{

  /** \brief. */
  typedef MultiVectorBase<Scalar>    MV_t;
  /** \brief. */
  typedef LinearOpBase<Scalar>       LO_t;

  //@}

  /** @name Constructors/initializers/accessors */
  //@{

  /// Construct to unintialize.
  BelosLinearOpWithSolve();

  /** \brief Initializes given precreated solver objects.
   *
   * \param lp [in] The linear problem that was used to initialize the
   * iterative solver.  The RHS and LHS arguments are set on this object to
   * solve a linear system.
   *
   * \param solverPL [in] Parameter list that is used by the iterative solver.
   *
   * \param iterativeSolver [in] The iterative solver manager that will be
   * used to solve for linear systems.  This has links to <tt>*lp</tt>,
   * <tt>*solverPL</tt> already embedded.
   *
   * \param fwdOpSrc [in] The source for the forward operator object defining
   * the linear system.  This object is not used here, it is just being
   * "remembered" so that it can be extracted by
   * <tt>BelosLinearOpWithSolveFactory::unitializeOp()</tt>.
   *
   * \param prec [in] The preconditioner object that was used to get the
   * preconditioners set in <tt>*lp</tt> This object is not used here, it is
   * just being "remembered" so that it can be extracted by
   * <tt>BelosLinearOpWithSolveFactory::unitializeOp()</tt>.
   *
   * \param isExternalPrec [in] Determines if the preconditioner was set by an
   * external client or was created internally by the
   * <tt>BelosLinearOpWithSolveFactory</tt> object.  This is not used here, it
   * is just being "remembered" so that it can be used in the logic for
   * <tt>BelosLinearOpWithSolveFactory::unitializeOp()</tt>.
   *
   * \param approxFwdOpSrc [in] The external approximate forward operator
   * object that was used to create the internal preconditioner.  This object
   * is not used here, it is just being "remembered" so that it can be
   * extracted by <tt>BelosLinearOpWithSolveFactory::unitializeOp()</tt>.
   *
   * \param supportSolveUse [in] Argument passed to
   * <tt>BelosLinearOpWithSolveFactory</tt> that is being remembered here to
   * be passed back to <tt>BelosLinearOpWithSolveFactory::unitializeOp()</tt>.
   *
   * ToDo: Finish documentation!
   */
  void initialize(
    const RCP<Belos::LinearProblem<Scalar,MV_t,LO_t> > &lp,
    const RCP<Teuchos::ParameterList> &solverPL,
    const RCP<Belos::SolverManager<Scalar,MV_t,LO_t> > &iterativeSolver,
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    const RCP<const PreconditionerBase<Scalar> > &prec,
    const bool isExternalPrec,
    const RCP<const LinearOpSourceBase<Scalar> > &approxFwdOpSrc,
    const ESupportSolveUse &supportSolveUse,
    const int convergenceTestFrequency
    );

  /** \brief . */
  RCP<const LinearOpSourceBase<Scalar> > extract_fwdOpSrc();

  /** \brief . */
  RCP<const PreconditionerBase<Scalar> > extract_prec();

  /** \brief . */
  bool isExternalPrec() const;

  /** \brief . */
  RCP<const LinearOpSourceBase<Scalar> > extract_approxFwdOpSrc();

  /** \brief . */
  ESupportSolveUse supportSolveUse() const;

  /** \brief Uninitializes and returns stored quantities.
   *
   * ToDo: Finish documentation!
   */
  void uninitialize(
    RCP<Belos::LinearProblem<Scalar,MV_t,LO_t> > *lp = NULL,
    RCP<Teuchos::ParameterList> *solverPL = NULL,
    RCP<Belos::SolverManager<Scalar,MV_t,LO_t> > *iterativeSolver = NULL,
    RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc = NULL,
    RCP<const PreconditionerBase<Scalar> > *prec = NULL,
    bool *isExternalPrec = NULL,
    RCP<const LinearOpSourceBase<Scalar> > *approxFwdOpSrc = NULL,
    ESupportSolveUse *supportSolveUse = NULL
    );

  //@}

  /** @name Overridden from LinearOpBase */
  //@{
  /** \brief. */
  RCP< const VectorSpaceBase<Scalar> > range() const;
  /** \brief. */
  RCP< const VectorSpaceBase<Scalar> > domain() const;
  /** \brief. */
  RCP<const LinearOpBase<Scalar> > clone() const;
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

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  RCP<Teuchos::ParameterList> getNonconstParameterList();
  /** \brief . */
  RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  RCP<const Teuchos::ParameterList> getValidParameters() const;

  //@}

protected:

  /** @name Overridden from LinearOpBase  */
  //@{
  /** \brief . */
  virtual bool opSupportedImpl(EOpTransp M_trans) const;
  /** \brief . */
  virtual void applyImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &X,
    const Ptr<MultiVectorBase<Scalar> > &Y,
    const Scalar alpha,
    const Scalar beta
    ) const;
  //@}

  /** @name Overridden from LinearOpWithSolveBase. */
  //@{
  /** \brief . */
  virtual bool solveSupportsImpl(EOpTransp M_trans) const;
  /** \brief . */
  virtual bool solveSupportsNewImpl(EOpTransp transp,
    const Ptr<const SolveCriteria<Scalar> > solveCriteria) const;
  /** \brief . */
  virtual bool solveSupportsSolveMeasureTypeImpl(
    EOpTransp M_trans, const SolveMeasureType& solveMeasureType
    ) const;
  /** \brief . */
  virtual SolveStatus<Scalar> solveImpl(
    const EOpTransp transp,
    const MultiVectorBase<Scalar> &B,
    const Ptr<MultiVectorBase<Scalar> > &X,
    const Ptr<const SolveCriteria<Scalar> > solveCriteria
    ) const;
  //@}
  
private:
  
  // ///////////////////////////////
  // Private data members


  RCP<Belos::LinearProblem<Scalar,MV_t,LO_t> > lp_;
  RCP<Teuchos::ParameterList> solverPL_;
  RCP<Belos::SolverManager<Scalar,MV_t,LO_t> > iterativeSolver_;
  int convergenceTestFrequency_;

  RCP<const LinearOpSourceBase<Scalar> > fwdOpSrc_;
  RCP<const PreconditionerBase<Scalar> > prec_;
  bool isExternalPrec_;
  RCP<const LinearOpSourceBase<Scalar> > approxFwdOpSrc_;
  ESupportSolveUse supportSolveUse_;

  typename Teuchos::ScalarTraits<Scalar>::magnitudeType defaultTol_;
                                                     
  void assertInitialized() const;

  std::string label_, filenameLHS_, filenameRHS_;
  mutable int counter_;
  
};


} // namespace Thyra


#endif	// THYRA_BELOS_LINEAR_OP_WITH_SOLVE_DECL_HPP
