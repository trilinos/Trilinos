//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2026 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluator_decl_hpp
#define Tempus_PhiEvaluator_decl_hpp

#include "Tempus_config.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"

#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_ModelEvaluator.hpp"

#include "Teuchos_DefaultComm.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Teuchos_CommHelpers.hpp"

namespace Tempus {

/**
 * @brief Selects which cached operators a linearization update must provide.
 *
 * `ONLY_MASS` prepares the mass operator \f$M\f$ and its inverse action;
 * `JACOBIAN_AND_MASS` additionally prepares the implicit residual Jacobian
 * \f$J_{\mathrm{impl}}\f$.  Operations involving the exponential operator
 * require the latter mode, while mass application and mass solves require
 * only the former.
 */
enum class PhiInitialization {
    ONLY_MASS,          ///< Assemble or require only \f$M\f$ and \f$M^{-1}\f$.
    JACOBIAN_AND_MASS,  ///< Also assemble \f$J_{\mathrm{impl}}\f$.
};

/**
 * @brief Owns the mass and implicit-Jacobian operators used by a PhiEvaluator.
 *
 * The supplied `Thyra::ModelEvaluator<Scalar>` is queried at a caller-selected
 * linearization point to assemble \f$M\f$ and \f$J_{\mathrm{impl}}\f$.  The
 * exponential operator constructed by this class is
 * \f$L=-\Delta t M^{-1}J_{\mathrm{impl}}\f$.  `Scalar` is the scalar type of
 * the Thyra vector and operator spaces.  A requested lumped mass matrix is a
 * row-sum diagonal approximation; otherwise inverse mass actions are supplied
 * by the `LinearOpWithSolveFactory` of the model.
 */
template <class Scalar>
class PhiLinearSolver {
 public:
  /**
   * @brief Associates this solver with an application model.
   *
   * @param appModel Const `Teuchos::RCP` to the model that creates and fills
   *   the mass and Jacobian operators; it must remain valid while this solver
   *   is used.
   * @param lumpMass When `true`, request the row-sum diagonal mass
   *   approximation on the next mass assembly.
   */
  PhiLinearSolver(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> appModel, bool lumpMass=false)
    : appModel_(appModel), lumpMass_(lumpMass) {
  }

  ~PhiLinearSolver() {}

  /** @brief Selects full or row-sum-lumped mass actions and invalidates cached mass operators when the selection changes.
   * @param lump Boolean selecting row-sum lumping when `true`.
   */
  void setLumpMassMatrix(const bool lump);

  /** @brief Returns `true` when both the cached `Thyra::LinearOpBase<Scalar>` mass operator and its inverse action are available. */
  bool massInitialized() const;

  /** @brief Releases cached mass, inverse-mass, and Jacobian operators; a subsequent operation must reassemble the operators it needs. */
  void clearMemory();

  /**
   * @brief Stores settings for `computeJacobianSpectrumBounds`.
   *
   * The `Teuchos::ParameterList` uses the "Eigensolver" entries from
   * `PhiEvaluator::getValidParametersBasic()`, including iterative-solver
   * controls and the dense fallback dimension.
   * @param pl Nonnull nonconst `Teuchos::RCP` to the eigensolver settings; set
   *   this before requesting spectrum bounds.
   */
  void setEigensolverParams(Teuchos::RCP<Teuchos::ParameterList> pl);

  /**
   * @brief Verifies that the model and the operators required by `mode` are cached.
   *
   * This inexpensive check does not assemble or repair any operator.
   * @param mode `PhiInitialization` level required by the pending operation.
   * @throws std::logic_error if the model, mass operators, or requested
   *   Jacobian operator is unavailable.
   */
  void checkInitialized(const PhiInitialization& mode) const;

  /** @brief Assembles \f$M\f$ at `inArgs` and caches its `Thyra::LinearOpBase<Scalar>` representation and inverse action. */
  void computeMassMatrix(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs);

  /** @brief Computes \f$Mf\f$ using the cached full or lumped mass operator.
   * @param Mf Writable `Thyra::VectorBase<Scalar>` in the mass range.
   * @param f Const `Thyra::VectorBase<Scalar>` in the mass domain.
   */
  void applyMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> Mf,
                  const Teuchos::RCP<const Thyra::VectorBase<Scalar>> f) const;

  /** @brief Computes the cached inverse-mass action \f$f=M^{-1}Mf\f$.
   * @param f Writable `Thyra::VectorBase<Scalar>` in the mass domain.
   * @param Mf Const mass-weighted `Thyra::VectorBase<Scalar>` in the mass range.
   */
  void solveMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> f,
                  const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mf) const;

  /** @brief Assembles and caches \f$J_{\mathrm{impl}}\f$ at `inArgs`; a mass operator must already be available. */
  void computeJacobian(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs);

  /** @brief Computes \f$J_{\mathrm{impl}}f\f$ using the cached implicit residual Jacobian.
   * @param Jf Writable `Thyra::VectorBase<Scalar>` in the Jacobian range.
   * @param f Const `Thyra::VectorBase<Scalar>` in the Jacobian domain.
   */
  void applyJacobian(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> Jf,
                      const Teuchos::RCP<const Thyra::VectorBase<Scalar>> f) const;

  /**
   * @brief Builds the linear operator \f$L=-dt M^{-1}J_{\mathrm{impl}}\f$.
   * @param dt Scalar time increment used to scale the operator.
   * @return Const `Teuchos::RCP` to the composed Thyra operator.
   * @pre Mass and Jacobian operators have been assembled at the intended
   *   linearization point.
   */
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> buildL(const Scalar dt) const;

  /**
   * @brief Builds the block extension used to evaluate a linear combination of phi functions.
   *
   * For `p = rhs_B.size() - 1`, returns
   * \f$\widetilde A=[L\ B;0\ K]\f$, where \f$L=-dt M^{-1}J_{\mathrm{impl}}\f$,
   * \f$K\in\mathbb{R}^{p\times p}\f$ has ones on its superdiagonal, and
   * \f$B\in V^{N\times p}\f$ contains `rhs_B[p]` through `rhs_B[1]` in
   * reverse column order.  Null right-hand-side RCPs produce zero columns.
   * @param dt Scalar time increment.
   * @param rhs_B Array view with `p + 1` mass-solved vectors; `p` must be at
   *   least one.
   * @return Const `Teuchos::RCP` to the two-block Thyra operator on
   *   \f$V^N\times\mathbb{R}^p\f$.
   */
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> buildATilde(
      const Scalar dt,
      const Teuchos::ArrayView<const Teuchos::RCP<const Thyra::VectorBase<Scalar>>> &rhs_B) const;

  /**
   * @brief Builds the extended initial vector \f$[x_0;e_p]\f$ for `buildATilde`.
   *
   * The first product block receives `x0` (or zeros for a null RCP); the
   * auxiliary block has dimension `p` and is zero except for its final entry,
   * which is one.
   * @param space Product `Thyra::VectorSpaceBase<Scalar>` returned by the
   *   extended operator domain.
   * @param x0 Const state-space `Thyra::VectorBase<Scalar>` for the first
   *   block, or null for a zero block.
   * @return Writable two-block `Thyra::ProductVectorBase<Scalar>`.
   */
  Teuchos::RCP<Thyra::ProductVectorBase<Scalar>> buildv(
      const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> space,
      const Teuchos::RCP<const Thyra::VectorBase<Scalar>> x0) const;

  /**
   * @brief Estimates zero-inclusive bounds of the spectrum of \f$-M^{-1}J_{\mathrm{impl}}\f$.
   *
   * For systems at or below the configured dense fallback dimension, the
   * method forms a replicated dense operator and computes all eigenvalues with
   * LAPACK.  Larger systems use Anasazi block Krylov-Schur Ritz values, which
   * may be unconverged and need not bound the full spectrum.  The results are
   * unscaled by `dt` and are intended for evaluator tuning such as a Leja
   * ellipse, not as certified spectral bounds.
   * @param a Output `double` set to the lesser of zero and the smallest
   *   observed real part.
   * @param b Output `double` set to the greater of zero and the largest
   *   observed real part.
   * @param c Output `double` set to the largest observed absolute imaginary
   *   part.
   */
  void computeJacobianSpectrumBounds(double& a, double& b, double& c);

  /** @brief Assembles \f$alpha M + beta J_{\mathrm{impl}}\f$ at `inArgs` and solves the resulting linear system.
   * @param inArgs Model-evaluation inputs used for the assembly.
   * @param x Writable `Thyra::VectorBase<Scalar>` receiving the solution.
   * @param Mf Const right-hand-side `Thyra::VectorBase<Scalar>`.
   * @param alpha Scalar multiplying the mass operator.
   * @param beta Scalar multiplying the Jacobian operator.
   * @return Thyra status reported by the linear solve.
   * @note This code path does not support mass lumping and is currently unused.
   */
  Thyra::SolveStatus<Scalar> assembleAndsolveMpJ(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
                                                 const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
                                                 const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mf,
                                                 Scalar alpha = 1., Scalar beta = 0.) const;

  /** @brief Solves \f$(alpha M + beta J_{\mathrm{impl}})x=Mf\f$ with cached operators.
   * @param x Writable `Thyra::VectorBase<Scalar>` receiving the solution.
   * @param Mf Const right-hand-side `Thyra::VectorBase<Scalar>`.
   * @param alpha Scalar multiplying the cached mass operator.
   * @param beta Scalar multiplying the cached Jacobian operator.
   * @return Thyra status reported by the linear solve.
   * @pre Mass and Jacobian operators have been assembled.
   * @note This code path does not currently support pre-conditioning and uses a default iterative solver.
   */
  Thyra::SolveStatus<Scalar> solveMpJ(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
                                      const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mf,
                                      Scalar alpha=1., Scalar beta=0.) const;

 private:
  /// Const `Teuchos::RCP` to the application model that assembles operators.
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > appModel_;
  /// `bool` selecting row-sum lumping for future mass assemblies.
  bool lumpMass_;

  /// Cached full or row-sum-lumped `Thyra::LinearOpBase<Scalar>` mass operator.
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> massMatrix_;
  /// Cached `Thyra::LinearOpBase<Scalar>` implementing the corresponding (full or lumped) inverse-mass action.
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> inverseMassMatrix_;

  /// Cached implicit residual `Thyra::LinearOpBase<Scalar>` Jacobian at the latest requested linearization point.
  Teuchos::RCP<Thyra::LinearOpBase<Scalar>> jacobianMatrix_;

  /// Previous Krylov-Schur eigenvectors for warm-starting spectrum estimates of \f$-M^{-1}J_{\mathrm{impl}}\f$.
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar>> evecsMinvJ_;

  /// Eigensolver settings forwarded from the owning `PhiEvaluator`.
  Teuchos::RCP<Teuchos::ParameterList> eigensolverPL_;
  /** @brief Builds the `p` by `p` superdiagonal shift operator used by the block extension. */
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> buildK(const Thyra::Ordinal max_phi_order) const;
  /** @brief Builds the `N` by `p` extension block from mass-solved phi right-hand sides. */
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> buildB(
      const Teuchos::ArrayView<const Teuchos::RCP<const Thyra::VectorBase<Scalar>>> &rhs_B) const;
};


/**
 * @brief Abstract mass-aware interface for actions of phi functions in exponential integrators.
 *
 * At a selected model linearization point, the base implementation uses
 * \f$L=-dtM^{-1}J_{\mathrm{impl}}\f$.  Positive-order right-hand sides use
 * the mass-aware interface from the exponential-integrator formulation:
 * callers provide `Mrhs_b = M b` and the evaluator obtains `b` internally.
 * Derived classes can implement `computeLinOpPhi` for their approximation method.
 *
 * @tparam Scalar Scalar type used by the model, Thyra spaces, and operators.
 */
template <class Scalar>
class PhiEvaluator
  : virtual public Teuchos::Describable,
    virtual public Teuchos::VerboseObject<PhiEvaluator<Scalar> > {
 public:
  /** @brief Constructs an unconfigured evaluator named "Phi Evaluator". */
  PhiEvaluator();

  /** @brief Constructs an unconfigured evaluator with a caller-supplied name.
   * @param name Descriptive `std::string` used in diagnostics and parameters.
   */
  PhiEvaluator(std::string name);

  ~PhiEvaluator() {}

  /// \name Basic PhiEvaluator Methods
  //@{
  /// @brief Returns this evaluator's descriptive `std::string` name.
  std::string getName() const { return name_; }

  /** @brief Replaces this evaluator's descriptive name.
   * @param name `std::string` used in descriptions and parameter-list names.
   */
  void setName(std::string name) { name_ = name; }

  /** @brief Returns the valid `Teuchos::ParameterList` schema for this evaluator, including derived-class settings. */
  virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  /**
   * @brief Returns the base evaluator parameter schema.
   *
   * "Lump Mass Matrix" selects the row-sum diagonal mass approximation.
   * "Constant Mass Matrix" permits reuse of the cached mass operator across
   * linearization updates; callers must set it consistently with the model.
   * The "Eigensolver" sublist controls spectrum estimates used by evaluators
   * that adapt to the current Jacobian.
   */
  Teuchos::RCP<Teuchos::ParameterList> getValidParametersBasic() const;

  /** @brief Returns a non-const `Teuchos::ParameterList` view of the valid parameter schema. */
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();

  /**
   * @brief Validates and applies base evaluator settings from `pl`.
   *
   * The mass-lumping and mass-caching entries control subsequent
   * `setLinearizationPoint` calls.  The "Eigensolver" sublist is copied and
   * forwarded to the `PhiLinearSolver` when one is available.
   * @param pl Nonnull `Teuchos::ParameterList` containing this evaluator's
   *   supported settings.
   */
  virtual void setPhiEvaluatorValues(Teuchos::RCP<Teuchos::ParameterList> pl);

  /**
   * @brief Marks the evaluator ready after verifying that a model and linear solver exist.
   *
   * This inexpensive operation only checks that model and internal phiSolver have been initialized
   * but does not assemble \f$M\f$ or \f$J_{\mathrm{impl}}\f$;
   * call `setLinearizationPoint` before an operation that requires them.
   */
  void initialize();

  /// @brief Returns whether `initialize()` has completed successfully.
  bool isInitialized() const { return isInitialized_; }

  /** @brief Throws `std::logic_error` unless `initialize()` has completed. */
  void checkInitialized() const;
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
  virtual std::string description() const;
  virtual void describe(Teuchos::FancyOStream& out,
                        const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  /** @brief Selects row-sum mass lumping for subsequent mass assemblies.
   * @param lumpMassMatrix Boolean selecting the lumped approximation.
   */
  void setLumpMassMatrix(bool lumpMassMatrix);
  /** @brief Selects whether a cached mass operator may be reused across linearization updates.
   * @param constantMassMatrix Boolean indicating that \f$M\f$ is constant in inArgs
   * and only depends on the model (i.e. does not need to be rebuilt if already present).
   * @note Switching from `true` to `false` clears cached mass and Jacobian
   *   operators; no model property is independently verified.
   */
  void setConstantMassMatrix(bool constantMassMatrix);

  /**
   * @brief Associates an application model and creates its `PhiLinearSolver`.
   * @param appModel Const `Teuchos::RCP` to the `Thyra::ModelEvaluator<Scalar>`
   *   from which \f$M\f$ and \f$J_{\mathrm{impl}}\f$ are assembled.
   * @note Call `initialize()` after setting the model.
   */
  void setModel(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > appModel);

  /**
   * @brief Updates cached operators for a model linearization point.
   *
   * The state, time, and other model inputs are copied from `inArgs`.  The
   * mass operator is assembled unless it is cached under the caller-declared
   * constant-mass setting.  `JACOBIAN_AND_MASS` also reassembles
   * \f$J_{\mathrm{impl}}\f$; use it before `computePhi`, `computePhis`, or
   * spectrum adaptation.  `ONLY_MASS` is sufficient for mass application or
   * inverse-mass actions.
   * @param inArgs `Thyra::ModelEvaluatorBase::InArgs<Scalar>` defining the
   *   intended linearization point.
   * @param mode Requested `PhiInitialization` level; the default assembles
   *   both mass and Jacobian.
   */
  void setLinearizationPoint(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
                             const PhiInitialization& mode = PhiInitialization::JACOBIAN_AND_MASS);

  /** @brief Adapts method-specific settings to the currently cached Jacobian.
   *
   * Call after `setLinearizationPoint` with `JACOBIAN_AND_MASS`.  The base
   * implementation is a no-op; derived evaluators may estimate the spectrum
   * or adjust approximation parameters.
   */
  virtual void adaptEvaluator() { this->isInitialized(); };

  /**
   * @brief Computes one mass-aware phi-function action at the cached linearization point.
   *
   * With \f$L=-c \Delta t M^{-1}J_{\mathrm{impl}}\f$, the base implementation applies
   * the cached inverse mass action to `Mrhs_b` before evaluating
   * \f$\varphi_k(L)b\f$.
   * @param x Writable `Thyra::VectorBase<Scalar>` receiving the result.
   * @param phi_order Nonnegative `int` phi-function index `k`.
   * @param cdt Scalar multiplier  \f$ c \Delta t\f$ for the cached mass-inverted Jacobian.
   * @param Mrhs_b Const state-space `Thyra::VectorBase<Scalar>` holding the
   *   mass-containing right-hand side \f$M b\f$; the base implementation applies the
   *   inverse mass action for every order.
   * @return `Thyra::SolveStatus<Scalar>` returned by the selected evaluator.
   * @pre The evaluator is initialized and mass and Jacobian operators have
   *   been assembled for the intended linearization point.
   */
  virtual Thyra::SolveStatus<Scalar> computePhi(
      const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
      const int phi_order, const Scalar cdt,
      const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &Mrhs_b);

  /**
   * @brief Computes a linear combination of mass-aware phi-function actions.
   *
   * For `p = Mrhs_B.size() - 1`, computes
   * \f$\sum_{s=0}^{p}\varphi_s(-c \Delta t M^{-1}J_{\mathrm{impl}})b_s\f$.
   * Each positive entry `Mrhs_B[s]` represents \f$M b_s\f$.
   * The array has one entry per order from zero through `p`; a null RCP
   * denotes a zero vector.
   * @param x Writable `Thyra::VectorBase<Scalar>` receiving the sum.
   * @param cdt Scalar multiplier  \f$ c \Delta t\f$ for the cached mass-inverted Jacobian.
   * @param Mrhs_B Nonempty `Teuchos::ArrayView` of const vector RCPs, indexed
   *   by phi order represents \f$M b_s\f$.
   * @return `Thyra::SolveStatus<Scalar>` returned by the extended evaluation.
   * @pre The evaluator is initialized and mass and Jacobian operators have
   *   been assembled for the intended linearization point.
   */
  virtual Thyra::SolveStatus<Scalar> computePhis(
      const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
      const Scalar cdt,
      const Teuchos::ArrayView<const Teuchos::RCP<const Thyra::VectorBase<Scalar>>> &Mrhs_B);

  /** @brief Applies the cached full or lumped mass operator.
   * @param Mf Writable `Thyra::VectorBase<Scalar>` receiving \f$Mf\f$.
   * @param f Const state-space `Thyra::VectorBase<Scalar>`.
   * @pre A mass operator has been assembled.
   */
  void applyMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> Mf,
                 const Teuchos::RCP<const Thyra::VectorBase<Scalar>> f) const;

  /** @brief Applies the cached inverse-mass action.
   * @param f Writable `Thyra::VectorBase<Scalar>` receiving \f$M^{-1}Mf\f$.
   * @param Mf Const mass-containing `Thyra::VectorBase<Scalar>`.
   * @pre A mass operator has been assembled.
   */
  void solveMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> f,
                 const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mf) const;

  /** @brief Applies the cached implicit residual Jacobian.
   * @param MJf Writable `Thyra::VectorBase<Scalar>` receiving
   *   \f$J_{\mathrm{impl}}f\f$.
   * @param f Const state-space `Thyra::VectorBase<Scalar>`.
   * @pre Mass and Jacobian operators have been assembled.
   */
  void applyJacobian(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> MJf,
                     const Teuchos::RCP<const Thyra::VectorBase<Scalar>> f) const;

 protected:
  /// Descriptive `std::string` used by parameter lists and diagnostics.
  std::string name_;
  /// `bool` selecting row-sum mass lumping for the managed linear solver.
  bool lumpMassMatrix_;
  /// Caller-declared `bool` permitting reuse of a cached mass operator.
  bool constantMassMatrix_;
  /// `bool` selecting the block extension for positive-order single-RHS calls.
  bool useAtildeForSingleRHS_;

  bool isInitialized_;  ///< Bool if PhiEvaluator is initialized.

  /// Const `Teuchos::RCP` to the model used to define the cached operators.
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> appModel_;
  /// `Teuchos::RCP` to the mass/Jacobian cache and operator builder.
  Teuchos::RCP<Tempus::PhiLinearSolver<Scalar>> phiLinSolv_;

  /// `Teuchos::ParameterList` forwarded to spectrum estimation when configured.
  Teuchos::RCP<Teuchos::ParameterList> eigensolverPL_;

  /// Owned `InArgs` copy describing the most recently requested linearization point.
  Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs_lin_;

  /**
   * @brief Evaluates a method-specific phi action in place.
   *
   * The base implementations of `computePhi` and `computePhis` construct an
   * ordinary or extended operator and use this hook to overwrite `v` with
   * \f$\varphi_{\mathrm{phi\_order}}(L)v\f$.  The extended route always
   * requests order zero because it evaluates a matrix exponential.
   * @param phi_order Nonnegative `int` phi index required by the method.
   * @param L Const `Teuchos::RCP` to the already scaled Thyra operator.
   * @param v Writable `Thyra::VectorBase<Scalar>` input/output vector in the
   *   operator domain/range.
   * @param cdt Scalar timestep factor retained for methods that need it for
   *   approximation scaling; `L` already contains this factor in the base path.
   * @return Method-specific `Thyra::SolveStatus<Scalar>`.
   */
  virtual Thyra::SolveStatus<Scalar> computeLinOpPhi(
      const int phi_order,
      const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> L,
      const Teuchos::Ptr<Thyra::VectorBase<Scalar>> v,
      const Scalar cdt=1.0
    ) = 0;
};

/**
 * @brief Materializes a Thyra linear operator as a replicated serial dense matrix.
 *
 * The returned `Teuchos::SerialDenseMatrix<int, Scalar>` has
 * `lop->range()->dim()` rows and `lop->domain()->dim()` columns.  Every rank
 * applies `lop` to the distributed identity and participates in a sum
 * reduction, so each rank receives the complete dense matrix.  This helper is
 * therefore intended only for small operators, such as the spectrum-estimation
 * dense fallback.
 * @param lop Const `Teuchos::RCP` to an operator with SPMD domain and range
 *   spaces.
 * @return Dense global matrix in column-major logical indexing.
 */
template <class Scalar>
Teuchos::SerialDenseMatrix<int, Scalar>
operatorToDense(Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> lop)
{
  const int numCols = static_cast<int>(lop->domain()->dim());
  const int numRows = static_cast<int>(lop->range()->dim());

  auto rangeSpmd =
    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<Scalar>>(
        lop->range(), true);

  auto domainSpmd =
    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<Scalar>>(
        lop->domain(), true);

  auto comm = rangeSpmd->getComm();
  const int rank = comm.is_null() ? 0 : comm->getRank();
  const int numProcs = comm.is_null() ? 1 : comm->getSize();

  const int rowOffset   = static_cast<int>(rangeSpmd->localOffset());
  const int rowLocalDim = static_cast<int>(rangeSpmd->localSubDim());

  const int colOffset   = static_cast<int>(domainSpmd->localOffset());
  const int colLocalDim = static_cast<int>(domainSpmd->localSubDim());

  Teuchos::RCP<Thyra::MultiVectorBase<Scalar>> Id =
      Thyra::createMembers(lop->domain(), numCols);

  Thyra::assign(Id.ptr(), Scalar(0));

  // Build the distributed identity: each rank writes only its owned domain
  // entries, then the operator application obtains all global columns.
  for (int jLocal = 0; jLocal < colLocalDim; ++jLocal) {
    const int jGlobal = colOffset + jLocal;

    Teuchos::RCP<Thyra::VectorBase<Scalar>> col_j = Id->col(jGlobal);
    Thyra::set_ele(jGlobal, Scalar(1), col_j.ptr());
  }

  Teuchos::RCP<Thyra::MultiVectorBase<Scalar>> Y =
      Thyra::createMembers(lop->range(), numCols);

  Thyra::apply(*lop, Thyra::NOTRANS, *Id, Y.ptr(), Scalar(1), Scalar(0));

  // Store locally owned output rows in global, column-major positions
  // (i + j * numRows) before the collective reconstruction.
  std::vector<Scalar> localDense(numRows * numCols, Scalar(0));
  std::vector<Scalar> globalDense(numRows * numCols, Scalar(0));

  for (int j = 0; j < numCols; ++j) {
    Teuchos::RCP<const Thyra::VectorBase<Scalar>> colY = Y->col(j);

    if (rowLocalDim > 0) {
      Thyra::Range1D localRange(rowOffset, rowOffset + rowLocalDim - 1);

      Thyra::ConstDetachedVectorView<Scalar> localView(*colY, localRange);

      for (int iLocal = 0; iLocal < rowLocalDim; ++iLocal) {
        const int iGlobal = rowOffset + iLocal;
        localDense[iGlobal + j * numRows] = localView[iGlobal];
      }
    }
  }

  // Because each row has one owner, this collective sum reconstructs the full
  // dense matrix on every rank.  A null communicator needs no reduction.
  if (comm.is_null()) {
    globalDense = localDense;
  }
  else {
    Teuchos::reduceAll(
        *comm,
        Teuchos::REDUCE_SUM,
        static_cast<long int>(localDense.size()),
        localDense.data(),
        globalDense.data());
  }

  Teuchos::SerialDenseMatrix<int, Scalar> globalDenseMat(numRows, numCols);

  for (int j = 0; j < numCols; ++j) {
    for (int i = 0; i < numRows; ++i) {
      globalDenseMat(i, j) = globalDense[i + j * numRows];
    }
  }

  return globalDenseMat;
}

/**
 * @brief Computes eigenvalues of a square dense matrix using LAPACK GEEV.
 * @tparam OrdinalType Ordinal type used by the dense matrix and LAPACK.
 * @tparam Scalar Scalar type of matrix entries and eigenvalue components.
 * @param A Square `Teuchos::SerialDenseMatrix` copied before LAPACK overwrites it.
 * @param eigs_re Output `Teuchos::Array<Scalar>` with one real component per
 *   matrix row.
 * @param eigs_im Output `Teuchos::Array<Scalar>` with one imaginary component
 *   per matrix row.
 * @return LAPACK `info` status from GEEV.
 */
template<typename OrdinalType, typename Scalar>
int denseEigenvalues(const Teuchos::SerialDenseMatrix<OrdinalType, Scalar>& A,
                           Teuchos::Array<Scalar>& eigs_re,
                           Teuchos::Array<Scalar>& eigs_im) {
    OrdinalType n = A.numRows();
    Teuchos::SerialDenseMatrix<OrdinalType, Scalar> A_copy(Teuchos::Copy, A);
    Teuchos::LAPACK<OrdinalType, Scalar> lapack;
    // Quick workspace query
    Scalar lworkQuery = 0.0;
    OrdinalType info = 0;
    lapack.GEEV('N', 'N', n, A_copy.values(), A_copy.stride(), eigs_re.getRawPtr(), eigs_im.getRawPtr(),
                nullptr, 1, nullptr, 1, &lworkQuery, -1, &info);
    // Solve
    OrdinalType lwork = static_cast<OrdinalType>(lworkQuery);
    Teuchos::Array<Scalar> work(lwork);
    lapack.GEEV('N', 'N', n, A_copy.values(), A_copy.stride(), eigs_re.getRawPtr(), eigs_im.getRawPtr(),
                nullptr, 1, nullptr, 1, work.getRawPtr(), lwork, &info);
    return info;
}

}  // namespace Tempus

#endif  // Tempus_PhiEvaluator_decl_hpp
