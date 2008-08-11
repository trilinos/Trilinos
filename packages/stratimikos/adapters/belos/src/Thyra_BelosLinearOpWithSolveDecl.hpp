
#ifndef THYRA_BELOS_LINEAR_OP_WITH_SOLVE_DECL_HPP
#define THYRA_BELOS_LINEAR_OP_WITH_SOLVE_DECL_HPP

#include "Thyra_SingleRhsLinearOpWithSolveBase.hpp"
#include "BelosSolverManager.hpp"
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
class BelosLinearOpWithSolve
  : virtual public LinearOpWithSolveBase<Scalar>                  // Public interface
  , virtual protected SingleScalarLinearOpWithSolveBase<Scalar>   // Implementation detail
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

  /// Calls <tt>initialize()</tt>
  BelosLinearOpWithSolve(
    const Teuchos::RCP<Belos::LinearProblem<Scalar,MV_t,LO_t> >         &lp
    ,const Teuchos::RCP<Teuchos::ParameterList>                         &solverPL
    ,const Teuchos::RCP<Belos::SolverManager<Scalar,MV_t,LO_t> >        &iterativeSolver
    ,const Teuchos::RCP<const LinearOpSourceBase<Scalar> >      &fwdOpSrc
    ,const Teuchos::RCP<const PreconditionerBase<Scalar> >             &prec
    ,const bool                                                                 isExternalPrec
    ,const Teuchos::RCP<const LinearOpSourceBase<Scalar> >              &approxFwdOpSrc
    ,const ESupportSolveUse                                                     &supportSolveUse
    );

  /** \brief Initializes given precreated solver objects.
   *
   * \param lp   [in] The linear problem that was used to initialize the iterative solver.
   *             The RHS and LHS arguments are set on this object to solve a linear system.
   * \param solverPL
   *             [in] Parameter list that is used by the iterative solver.
   * \param iterativeSolver
   *             [in] The iterative solver manager that will be used to solve for linear systems.  This has
   *             links to <tt>*lp</tt>, <tt>*solverPL</tt> already embedded.
   * \param fwdOpSrc
   *             [in] The source for the forward operator object defining the linear system.
   *             This object is not used here, it is just being "remembered" so that it can be extracted
   *             by <tt>BelosLinearOpWithSolveFactory::unitializeOp()</tt>.
   * \param prec
   *             [in] The preconditioner object that was used to get the precondtioners set in <tt>*lp</tt>
   *             This object is not used here, it is just being "remembered" so that it can be extracted
   *             by <tt>BelosLinearOpWithSolveFactory::unitializeOp()</tt>.
   * \param isExternalPrec
   *             [in] Determines if the preconditioner was set by an external client or was created internally
   *             by the <tt>BelosLinearOpWithSolveFactory</tt> object.
   *             This is not used here, it is just being "remembered" so that it can be used in the logic for
   *             <tt>BelosLinearOpWithSolveFactory::unitializeOp()</tt>.
   * \param approxFwdOpSrc
   *             [in] The external approximate forward operator object that was used to create the internal
   *             preconditioner.
   *             This object is not used here, it is just being "remembered" so that it can be extracted
   *             by <tt>BelosLinearOpWithSolveFactory::unitializeOp()</tt>.
   * \param supportSolveUse
   *             [in] Argument passed to <tt>BelosLinearOpWithSolveFactory</tt> that is being remembered here
   *             to be passed back to <tt>BelosLinearOpWithSolveFactory::unitializeOp()</tt>.
   *
   * ToDo: Finish documentation!
   */
  void initialize(
    const Teuchos::RCP<Belos::LinearProblem<Scalar,MV_t,LO_t> >         &lp
    ,const Teuchos::RCP<Teuchos::ParameterList>                         &solverPL
    ,const Teuchos::RCP<Belos::SolverManager<Scalar,MV_t,LO_t> >        &iterativeSolver
    ,const Teuchos::RCP<const LinearOpSourceBase<Scalar> >              &fwdOpSrc
    ,const Teuchos::RCP<const PreconditionerBase<Scalar> >              &prec
    ,const bool                                                         isExternalPrec
    ,const Teuchos::RCP<const LinearOpSourceBase<Scalar> >              &approxFwdOpSrc
    ,const ESupportSolveUse                                             &supportSolveUse
    );

  /** \brief . */
  Teuchos::RCP<const LinearOpSourceBase<Scalar> > extract_fwdOpSrc();

  /** \brief . */
  Teuchos::RCP<const PreconditionerBase<Scalar> > extract_prec();

  /** \brief . */
  bool isExternalPrec() const;

  /** \brief . */
  Teuchos::RCP<const LinearOpSourceBase<Scalar> > extract_approxFwdOpSrc();

  /** \brief . */
  ESupportSolveUse supportSolveUse() const;

  /** \brief Uninitializes and returns stored quantities.
   *
   * ToDo: Finish documentation!
   */
  void uninitialize(
    Teuchos::RCP<Belos::LinearProblem<Scalar,MV_t,LO_t> >         *lp                        = NULL
    ,Teuchos::RCP<Teuchos::ParameterList>                         *solverPL                   = NULL
    ,Teuchos::RCP<Belos::SolverManager<Scalar,MV_t,LO_t> >        *iterativeSolver           = NULL
    ,Teuchos::RCP<const LinearOpSourceBase<Scalar> >              *fwdOpSrc                  = NULL
    ,Teuchos::RCP<const PreconditionerBase<Scalar> >              *prec                      = NULL
    ,bool                                                         *isExternalPrec            = NULL
    ,Teuchos::RCP<const LinearOpSourceBase<Scalar> >              *approxFwdOpSrc            = NULL
    ,ESupportSolveUse                                             *supportSolveUse           = NULL
    );

  //@}

  /** @name Overridden from LinearOpBase */
  //@{
  /** \brief. */
  Teuchos::RCP< const VectorSpaceBase<Scalar> > range() const;
  /** \brief. */
  Teuchos::RCP< const VectorSpaceBase<Scalar> > domain() const;
  /** \brief. */
  Teuchos::RCP<const LinearOpBase<Scalar> > clone() const;
  //@}

  /** @name Overridden from Teuchos::Describable */
  //@{
  /** \brief . */
  std::string description() const;
  /** \brief . */
  void describe(
    Teuchos::FancyOStream                &out
    ,const Teuchos::EVerbosityLevel      verbLevel
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

protected:

  /** @name Overridden from SingleScalarLinearOpBase */
  //@{
  /** \brief . */
  bool opSupported(EOpTransp M_trans) const;
  /** \brief . */
  void apply(
    const EOpTransp                     M_trans
    ,const MultiVectorBase<Scalar>    &X
    ,MultiVectorBase<Scalar>          *Y
    ,const Scalar                     alpha
    ,const Scalar                     beta
    ) const;
  //@}

  /** @name Overridden from SingleScalarLinearOpWithSolveBase */
  //@{
  /** \brief . */
  bool solveSupportsTrans(EOpTransp M_trans) const;
  /** \brief . */
  bool solveSupportsSolveMeasureType(EOpTransp M_trans, const SolveMeasureType& solveMeasureType) const;
  /** \brief . */
  void solve(
    const EOpTransp                         M_trans
    ,const MultiVectorBase<Scalar>        &B
    ,MultiVectorBase<Scalar>              *X
    ,const int                            numBlocks
    ,const BlockSolveCriteria<Scalar>     blockSolveCriteria[]
    ,SolveStatus<Scalar>                  blockSolveStatus[]
    ) const;
  //@}
  
private:
  
  // ///////////////////////////////
  // Private data members


  Teuchos::RCP<Belos::LinearProblem<Scalar,MV_t,LO_t> >           lp_;
  Teuchos::RCP<Teuchos::ParameterList>                            solverPL_;
  Teuchos::RCP<Belos::SolverManager<Scalar,MV_t,LO_t> >           iterativeSolver_;

  Teuchos::RCP<const LinearOpSourceBase<Scalar> >                 fwdOpSrc_;
  Teuchos::RCP<const PreconditionerBase<Scalar> >                 prec_;
  bool                                                            isExternalPrec_;
  Teuchos::RCP<const LinearOpSourceBase<Scalar> >                 approxFwdOpSrc_;
  ESupportSolveUse                                                supportSolveUse_;

  typename Teuchos::ScalarTraits<Scalar>::magnitudeType           defaultTol_;
                                                     
  void assertInitialized() const;
  
};

} // namespace Thyra

#endif	// THYRA_BELOS_LINEAR_OP_WITH_SOLVE_DECL_HPP
