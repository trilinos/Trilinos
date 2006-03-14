
#ifndef THYRA_BELOS_LINEAR_OP_WITH_SOLVE_DECL_HPP
#define THYRA_BELOS_LINEAR_OP_WITH_SOLVE_DECL_HPP

#include "Thyra_SingleRhsLinearOpWithSolveBase.hpp"
#include "BelosIterativeSolver.hpp"
#include "BelosStatusTestResNorm.hpp"
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
    ,const Teuchos::RefCountPtr<Belos::StatusTestResNorm<Scalar,MV_t,LO_t> >    &resNormST
    ,const Teuchos::RefCountPtr<Belos::IterativeSolver<Scalar,MV_t,LO_t> >      &iterativeSolver
    );

  /** \brief Initializes given precreated solver objects ...
   *
   */
  void initialize(
    const Teuchos::RefCountPtr<Belos::LinearProblem<Scalar,MV_t,LO_t> >         &lp
    ,const Teuchos::RefCountPtr<Belos::StatusTestResNorm<Scalar,MV_t,LO_t> >    &resNormST
    ,const Teuchos::RefCountPtr<Belos::IterativeSolver<Scalar,MV_t,LO_t> >      &iterativeSolver
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
  bool solveSupportsSolveTolType(ETransp M_trans, ESolveTolType solveTolType) const;
  /** \brief . */
  int defaultSolveMaxIterations(ETransp M_trans, ESolveTolType solveTolType) const;
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

  typedef Belos::StatusTestResNorm<Scalar,MV_t,LO_t>                      StatusTestResNorm_t;

  Teuchos::RefCountPtr<Belos::LinearProblem<Scalar,MV_t,LO_t> >           lp_;
  Teuchos::RefCountPtr<StatusTestResNorm_t>                               resNormST_;
  Teuchos::RefCountPtr<Belos::IterativeSolver<Scalar,MV_t,LO_t> >         iterativeSolver_;

  typename Teuchos::ScalarTraits<Scalar>::magnitudeType                   defaultTol_;
                                                     
  void assertInitialized() const;
  
};

} // namespace Thyra

#endif	// THYRA_BELOS_LINEAR_OP_WITH_SOLVE_DECL_HPP
