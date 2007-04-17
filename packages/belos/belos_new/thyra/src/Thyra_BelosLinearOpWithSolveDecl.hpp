
#ifndef THYRA_BELOS_LINEAR_OP_WITH_SOLVE_DECL_HPP
#define THYRA_BELOS_LINEAR_OP_WITH_SOLVE_DECL_HPP

#include "Thyra_SingleRhsLinearOpWithSolveBase.hpp"
#include "BelosIterativeSolver.hpp"
#include "BelosStatusTestResNorm.hpp"
#include "BelosOutputManager.hpp"
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
    const Teuchos::RefCountPtr<Belos::LinearProblem<Scalar,MV_t,LO_t> >         &lp
    ,const bool                                                                 adjustableBlockSize
    ,const int                                                                  maxNumberOfKrylovVectors
    ,const Teuchos::RefCountPtr<Teuchos::ParameterList>                         &gmresPL
    ,const Teuchos::RefCountPtr<Belos::StatusTestResNorm<Scalar,MV_t,LO_t> >    &resNormST
    ,const Teuchos::RefCountPtr<Belos::IterativeSolver<Scalar,MV_t,LO_t> >      &iterativeSolver
    ,const Teuchos::RefCountPtr<Belos::OutputManager<Scalar> >                  &outputManager
    ,const Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> >              &fwdOpSrc
    ,const Teuchos::RefCountPtr<const PreconditionerBase<Scalar> >              &prec
    ,const bool                                                                 isExternalPrec
    ,const Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> >              &approxFwdOpSrc
    ,const ESupportSolveUse                                                     &supportSolveUse
    );

  /** \brief Initializes given precreated solver objects.
   *
   * \param lp   [in] The linear problem that was used to initialize the iterative solver.
   *             The RHS and LHS arguments are set on this object to solve a linear system.
   * \param adjustableBlockSize
   *             [in] If <tt>true</tt>, then the block size in <tt>*lp</tt> will be adjusted according
   *             to the linear system being solved automatically.
   * \param maxNumberOfKrylovVectors
   *             [in] Total number of Krylov vectors that can be stored and manipulated.  This more-or-less
   *             bounds the total amount of storage that the algorithm can use.
   * \param gmresPL
   *             [in] Parameter list that is used for GMRES if GMRES is used.
   * \param resNormST
   *             [in] The residual norm status test that was used to initialize the iterative solver.
   *             This is accessed to set the relative tolerance.
   * \param iterativeSolver
   *             [in] The iterative solver that will be used to solve for linear systems.  This has
   *             links to <tt>*lp</tt>, <tt>*resNormST</tt> and <tt>*outputManager</tt> already
   *             embedded.
   * \param outputManager
   *             [in] The output manager that was passed to the iterative solver.
   *             This is used to reset the output level and the outptu stream on a solve by solve basis.
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
    const Teuchos::RefCountPtr<Belos::LinearProblem<Scalar,MV_t,LO_t> >         &lp
    ,const bool                                                                 adjustableBlockSize
    ,const int                                                                  maxNumberOfKrylovVectors
    ,const Teuchos::RefCountPtr<Teuchos::ParameterList>                         &gmresPL
    ,const Teuchos::RefCountPtr<Belos::StatusTestResNorm<Scalar,MV_t,LO_t> >    &resNormST
    ,const Teuchos::RefCountPtr<Belos::IterativeSolver<Scalar,MV_t,LO_t> >      &iterativeSolver
    ,const Teuchos::RefCountPtr<Belos::OutputManager<Scalar> >                  &outputManager
    ,const Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> >              &fwdOpSrc
    ,const Teuchos::RefCountPtr<const PreconditionerBase<Scalar> >              &prec
    ,const bool                                                                 isExternalPrec
    ,const Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> >              &approxFwdOpSrc
    ,const ESupportSolveUse                                                     &supportSolveUse
    );

  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> > extract_fwdOpSrc();

  /** \brief . */
  Teuchos::RefCountPtr<const PreconditionerBase<Scalar> > extract_prec();

  /** \brief . */
  bool isExternalPrec() const;

  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> > extract_approxFwdOpSrc();

  /** \brief . */
  ESupportSolveUse supportSolveUse() const;

  /** \brief Uninitializes and returns stored quantities.
   *
   * ToDo: Finish documentation!
   */
  void uninitialize(
    Teuchos::RefCountPtr<Belos::LinearProblem<Scalar,MV_t,LO_t> >         *lp                        = NULL
    ,bool                                                                 *adjustableBlockSize       = NULL
    ,int                                                                  *maxNumberOfKrylovVectors  = NULL
    ,Teuchos::RefCountPtr<Teuchos::ParameterList>                         *gmresPL                   = NULL
    ,Teuchos::RefCountPtr<Belos::StatusTestResNorm<Scalar,MV_t,LO_t> >    *resNormST                 = NULL
    ,Teuchos::RefCountPtr<Belos::IterativeSolver<Scalar,MV_t,LO_t> >      *iterativeSolver           = NULL
    ,Teuchos::RefCountPtr<Belos::OutputManager<Scalar> >                  *outputManager             = NULL
    ,Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> >              *fwdOpSrc                  = NULL
    ,Teuchos::RefCountPtr<const PreconditionerBase<Scalar> >              *prec                      = NULL
    ,bool                                                                 *isExternalPrec            = NULL
    ,Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> >              *approxFwdOpSrc            = NULL
    ,ESupportSolveUse                                                     *supportSolveUse           = NULL
    );

  //@}

  /** @name Overridden from LinearOpBase */
  //@{
  /** \brief. */
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > range() const;
  /** \brief. */
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > domain() const;
  /** \brief. */
  Teuchos::RefCountPtr<const LinearOpBase<Scalar> > clone() const;
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
  void setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> getParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getValidParameters() const;

  //@}

protected:

  /** @name Overridden from SingleScalarLinearOpBase */
  //@{
  /** \brief . */
  bool opSupported(ETransp M_trans) const;
  /** \brief . */
  void apply(
    const ETransp                     M_trans
    ,const MultiVectorBase<Scalar>    &X
    ,MultiVectorBase<Scalar>          *Y
    ,const Scalar                     alpha
    ,const Scalar                     beta
    ) const;
  //@}

  /** @name Overridden from SingleScalarLinearOpWithSolveBase */
  //@{
  /** \brief . */
  bool solveSupportsTrans(ETransp M_trans) const;
  /** \brief . */
  bool solveSupportsSolveMeasureType(ETransp M_trans, const SolveMeasureType& solveMeasureType) const;
  /** \brief . */
  void solve(
    const ETransp                         M_trans
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

  typedef Belos::StatusTestResNorm<Scalar,MV_t,LO_t>                      StatusTestResNorm_t;

  Teuchos::RefCountPtr<Belos::LinearProblem<Scalar,MV_t,LO_t> >           lp_;
  bool                                                                    adjustableBlockSize_;
  int                                                                     maxNumberOfKrylovVectors_;
  Teuchos::RefCountPtr<Teuchos::ParameterList>                            gmresPL_;
  Teuchos::RefCountPtr<StatusTestResNorm_t>                               resNormST_;
  Teuchos::RefCountPtr<Belos::IterativeSolver<Scalar,MV_t,LO_t> >         iterativeSolver_;
  Teuchos::RefCountPtr<Belos::OutputManager<Scalar> >                     outputManager_;

  Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> >                 fwdOpSrc_;
  Teuchos::RefCountPtr<const PreconditionerBase<Scalar> >                 prec_;
  bool                                                                    isExternalPrec_;
  Teuchos::RefCountPtr<const LinearOpSourceBase<Scalar> >                 approxFwdOpSrc_;
  ESupportSolveUse                                                        supportSolveUse_;

  typename Teuchos::ScalarTraits<Scalar>::magnitudeType                   defaultTol_;
                                                     
  void assertInitialized() const;
  
};

} // namespace Thyra

#endif	// THYRA_BELOS_LINEAR_OP_WITH_SOLVE_DECL_HPP
