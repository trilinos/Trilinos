//@HEADER
// *****************************************************************************
// TODO
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluator_decl_hpp
#define Tempus_PhiEvaluator_decl_hpp

#include "Tempus_config.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_ModelEvaluator.hpp"

namespace Tempus {

/** \brief PhiLinearSolver
 * Uses the ModelEvaluator to compute and hold the Mass matrix and Jacobian
 */
template <class Scalar>
class PhiLinearSolver {
 public:
  PhiLinearSolver(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> appModel, bool lumpMass=false)
    : appModel_(appModel), lumpMass_(lumpMass), isInitialized_(false) {
  }

  ~PhiLinearSolver() {}

  void setLumpMassMatrix(const bool lump);

  /** \brief Initialize PhiSolver
   *
   *  This function will check if mass matrix and Jacobian have been computed.
   *  This function does not make member data consistent, but just checks it.
   *  This ensures it is inexpensive.
   */
  void initialize();
  /// Return if PhiSolver is initialized.
  bool isInitialized() const { return isInitialized_; }
  void checkInitialized() const;

  void computeMassMatrix(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs);
  void applyMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> Mf, const Teuchos::RCP<const Thyra::VectorBase<Scalar>> f) const;
  void solveMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> f, const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mf) const;

  void computeJacobian(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs);
  void applyJacobian(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> Jf, const Teuchos::RCP<const Thyra::VectorBase<Scalar>> f) const;

  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> buildL(const Scalar dt) const;

  // TODO: make that one public function
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> buildATilde(const Scalar dt);
  void buildK(const Thyra::Ordinal n);
  void buildb(const Teuchos::ArrayView<const Teuchos::RCP<const Thyra::VectorBase<Scalar>>> &rhs_B);
  Teuchos::RCP<Thyra::ProductVectorBase<Scalar>> buildv(const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> space,
							const Teuchos::RCP<const Thyra::VectorBase<Scalar>> x0) const;
  Teuchos::Tuple<Scalar, 3> computeJacobianSpectrumBounds(int inev);

  // Solve Mass plus Jacobian, for given inArgs (recompute matrices from from ModelEvaluator)
  Thyra::SolveStatus<Scalar> assembleAndsolveMpJ(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
                                                 const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
                                                 const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mf,
                                                 Scalar alpha = 1., Scalar beta = 0.) const;

  // Solve Mass plus Jacobian, for precomputed Mass and Jacobian
  Thyra::SolveStatus<Scalar> solveMpJ(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
                                      const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mf,
                                      Scalar alpha=1., Scalar beta=0.) const;

 private:
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > appModel_;
  bool lumpMass_;
  bool isInitialized_;

  Teuchos::RCP<Thyra::LinearOpBase<Scalar>> fullMassMatrix_;
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> lumpedMassMatrix_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > lumpedMassDiagonal_;

  // massMatrix_ and inverseMassMatrix_ are either lumped, or not, depending on bool lumpMass_
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> massMatrix_;
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> inverseMassMatrix_;

  Teuchos::RCP<Thyra::LinearOpBase<Scalar>> jacobianMatrix_;

  // internal variables for extended matrix strategy
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> expMassMatrix_;
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> Atilde_;
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> KMatrix_;
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> bMatrix_;

  // internal storage for eigenvectors of the M^{-1}*J LinOp
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar>> evecsMinvJ_;
};


/** \brief PhiEvaluator evaluates
 *
 *  \f$[x = \varphi_k(J) b]\f$, where
 *
 *   - b is a right hand side vector
 *   - J is a linear operator
 */
template <class Scalar>
class PhiEvaluator
  : virtual public Teuchos::Describable,
    virtual public Teuchos::VerboseObject<PhiEvaluator<Scalar> > {
 public:
  PhiEvaluator();

  PhiEvaluator(
      std::string name);

  ~PhiEvaluator() {}

  /// \name Basic PhiEvaluator Methods
  //@{

  /// Make a shallow copy of PhiEvaluator (i.e., only RCPs).
  void copy(Teuchos::RCP<const PhiEvaluator<Scalar> > sh);
  //@}

  /// \name Accessor methods
  //@{
  /// Get this PhiEvaluator's name
  std::string getName() const { return name_; }

  /// Set this PhiEvaluator's name
  void setName(std::string name) { name_ = name; }

  /// Return a valid ParameterList with current settings.
  virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  /// Return a valid ParameterList with basic settings.
  Teuchos::RCP<Teuchos::ParameterList> getValidParametersBasic() const;

  /// Return a valid non-const ParameterList with current settings.
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();

  /// Set the parameters from a ParameterList
  virtual void setPhiEvaluatorValues(Teuchos::RCP<Teuchos::ParameterList> pl);

  /// \name Overridden from Teuchos::Describable
  //@{
  virtual std::string description() const;
  virtual void describe(Teuchos::FancyOStream& out,
                        const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  /** \brief Initialize PhiEvaluator
   *
   *  This function will check if all member data is initialized
   *  and is consistent.  This function does not make member data
   *  consistent, but just checks it.  This ensures it is inexpensive.
   */
  void initialize();

  /// Return if PhiEvaluator is initialized.
  bool isInitialized() const { return isInitialized_; }

  void checkInitialized() const;

  void setLumpMassMatrix(bool lumpMassMatrix);

  /// set the ModelEvaluator
  void setModel(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > appModel);

  /** \brief   Set the linearization point for the Jacobian calculation
   *
   *  The linearization point x and time t are taken from inArgs.
   *  This computes the Mass and Jacobian matrix for future use.
   */
  void setLinearizationPoint(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs);

  /** \brief   Adapt internal PhiEvaluator hyperparameters to the current Jacobian.
   *
   *  Called after setLinearizationPoint.  Default implementation is a null-op.
   *  This method should be overridden by inheriting PhiEvaluator implementations when
   *  there are hyperparameters of the method that require analysis of the Jacobian.
   */
  virtual void adaptEvaluator() { this->isInitialized(); };

  /** \brief  Compute the Phi_k function of cdt*Jacobian for right hand side Mrhs_b
   *
   *  phi_order is the index of the phi-function Phi_k.
   *  For an implicit model, the right hand side contains the mass matrix M,
   *  which is solved as part of this method.
   */
  virtual Thyra::SolveStatus<Scalar> computePhi(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
						const int phi_order, const Scalar cdt,
						const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &Mrhs_b);

  /** \brief  Compute the Phi function of cdt times Jacobian for a linear combination with right hand side vectors Mrhs_B
   *
   *  The vectors in Mrhs_B are at the index of the vector corresponding to the phi_order of the
   *  respective Phi function, Mrhs_b[0] is the rhs for the matrix exponential, Mrhs_b[1] is the
   *  right-hand side for the phi_1 function. A Teuchos::null RCP is interpreted as a zero vector.
   *  For an implicit model, the right hand side contains a multiplication with the mass matrix M,
   *  which is solved as part of this method.
   */
  virtual Thyra::SolveStatus<Scalar> computePhis(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
                                                 const Scalar cdt,
                                                 const Teuchos::ArrayView<const Teuchos::RCP<const Thyra::VectorBase<Scalar>>> &Mrhs_B);

  // Multiply the mass matrix (lumped or not) with right hand side f
  void applyMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> Mf, const Teuchos::RCP<const Thyra::VectorBase<Scalar>> f) const;

  // Invert the mass matrix (lumped or not) with right hand side Mf
  void solveMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> f, const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mf) const;

  // Multiply the MassJacobian matrix with right hand side Mf
  void applyJacobian(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> MJf, const Teuchos::RCP<const Thyra::VectorBase<Scalar>> f) const;

 protected:
  std::string name_;
  bool lumpMassMatrix_;
  bool useAtildeForSingleRHS_;

  bool isInitialized_;  ///< Bool if PhiEvaluator is initialized.

  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> appModel_;
  Teuchos::RCP<Tempus::PhiLinearSolver<Scalar>> phiLinSolv_;

  //mutable
  Teuchos::RCP<const Thyra::ModelEvaluatorBase::InArgs<Scalar>> inArgs_lin_;


  /** \brief  Internal method for a LinOp, used for default impl. of computePhi/computePhis
   *
   *  Computes v := phi_{phi_order}(L)v in place (overwriting the rhs with the result)
   */
  virtual Thyra::SolveStatus<Scalar> computeLinOpPhi(const int phi_order,
						     const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> L,
						     const Teuchos::Ptr<Thyra::VectorBase<Scalar>> v,
						     const Scalar cdt=1.0
						     ) = 0;
};

}  // namespace Tempus

#endif  // Tempus_PhiEvaluator_decl_hpp
