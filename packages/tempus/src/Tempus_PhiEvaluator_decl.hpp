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
#include "Thyra_ModelEvaluator.hpp"

namespace Tempus {

/** \brief PhiLinearSolver
 * Uses the ModelEvaluator to compute and hold the Mass matrix and Jacobian
 */
template <class Scalar>
class PhiLinearSolver {
 public:
  PhiLinearSolver(const Teuchos::RCP<const Thyra::ModelEvaluator<double>> appModel, bool lumpMass=false)
    : appModel_(appModel), lumpMass_(lumpMass) {
  }

  ~PhiLinearSolver() {}

  void computeMassMatrix(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs);
  void applyMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> Mf, const Teuchos::RCP<const Thyra::VectorBase<Scalar>> f);
  void solveMass(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> f, const Teuchos::RCP<const Thyra::VectorBase<Scalar>> Mx);

  void computeJacobian(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs);


 private:
  Teuchos::RCP<const Thyra::ModelEvaluator<double>> appModel_;
  const bool lumpMass_;

  mutable Teuchos::RCP<Thyra::LinearOpBase<Scalar>> fullMassMatrix_;
  mutable Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> lumpMassMatrix_;
  mutable Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> invMassMatrix_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> prototypeInArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> prototypeOutArgs_;

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

  /// Make a shallow copy of PhiEvaluator (i.e., only RCPs to states and
  /// interpolator).
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
  void initialize() const;

  /// Return if PhiEvaluator is initialized.
  bool isInitialized() { return isInitialized_; }

  void checkInitialized();

  /// set the ModelEvaluator
  void setModel(const Teuchos::RCP<const Thyra::ModelEvaluator<double> > appModel);

  /// Set the linearization point for the Jacobian calculation
  virtual void setLinearizationPoint(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs);

  /// compute the Phi_k function of cdt times Jacobian for right hand side rhs_b
  virtual void computePhi(const Teuchos::Ptr<Thyra::VectorBase<Scalar>>,
                          int k, Scalar cdt, const Teuchos::RCP<const Thyra::VectorBase<Scalar>> rhs_b) = 0;

 protected:
  std::string name_;

  mutable bool isInitialized_;  ///< Bool if PhiEvaluator is initialized.

  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> appModel_;
  Teuchos::RCP<const Tempus::PhiLinearSolver<Scalar>> phiLinSolv_;
};

}  // namespace Tempus

#endif  // Tempus_PhiEvaluator_decl_hpp
