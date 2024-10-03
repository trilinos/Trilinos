//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_ModelEvaluatorPairPartIMEX_StaggeredFSA_decl_hpp
#define Tempus_ModelEvaluatorPairPartIMEX_StaggeredFSA_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_SensitivityModelEvaluatorBase.hpp"
#include "Tempus_WrapperModelEvaluatorPairPartIMEX_Basic.hpp"
#include "Tempus_StaggeredForwardSensitivityModelEvaluator.hpp"

#include "Thyra_ProductMultiVectorBase.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"

namespace Tempus {

/** \brief Specialization of IMEX-Part ME for "combined" FSA method.
 *
 * This specializes the implementation of several parts of
 * WrapperModelEvaluatorPairPartIMEX_Basic for forward-sensitivity analysis
 * with StaggeredForwardSensitivityModelEvaluator.  It deals with the product
 * structure of the sensitivity solution vectors and incorporates the
 * sensitivity of the implicit term with respect to the explicit-only vector.
 */
template <typename Scalar>
class WrapperModelEvaluatorPairPartIMEX_StaggeredFSA
  : public SensitivityModelEvaluatorBase<Scalar>,
    public WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> {
 public:
  /// Constructor
  WrapperModelEvaluatorPairPartIMEX_StaggeredFSA(
      const Teuchos::RCP<
          const WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >& forwardModel,
      const bool is_pseudotransient,
      const Teuchos::RCP<const Teuchos::ParameterList>& pList = Teuchos::null);

  /// Destructor
  virtual ~WrapperModelEvaluatorPairPartIMEX_StaggeredFSA() {}

  /// Initialize after setting member data.
  virtual void initialize();

  /// \name Methods that apply to both explicit and implicit terms.
  //@{

  virtual Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(
      int i) const;

  //@}

  //@{ \name Accessors

  /// Extract IMEX vector from a full solution vector
  virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getIMEXVector(
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& full) const;

  /// Extract IMEX vector for reading
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar> > getIMEXVector(
      const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& full) const;

  /// Extract explicit-only vector from a full solution vector
  virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getExplicitOnlyVector(
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& full) const;

  /// Extract explicit-only vector for reading
  virtual Teuchos::RCP<const Thyra::VectorBase<Scalar> > getExplicitOnlyVector(
      const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& full) const;

  //@}

  /// \name Overridden from Tempus::SensitivityModelEvaluatorBase
  //@{

  /// Get the underlying forward model
  virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getForwardModel()
      const;

  /// Set solution history from forward state evaluation (for interpolation)
  virtual void setForwardSolutionHistory(
      const Teuchos::RCP<const Tempus::SolutionHistory<Scalar> >& sh);

  /// Set solution state from forward state evaluation (for frozen state)
  virtual void setForwardSolutionState(
      const Teuchos::RCP<const Tempus::SolutionState<Scalar> >& s);

  /// Set the solver of the underlying model if you want to reuse it
  virtual void setSolver(
      const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
      const bool force_W_update);

  //@}

  /// \name Overridden from Thyra::StateFuncModelEvaluatorBase
  //@{

  virtual Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  virtual void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const;

  //@}

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

 protected:
  /// Build implicit x and end explicit y states from forward_state_
  void buildIMEXStates() const;

  typedef WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> Base;
  typedef Thyra::DefaultMultiVectorProductVectorSpace<Scalar> DMVPVS;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
  typedef Thyra::ProductMultiVectorBase<Scalar> PMVB;
  typedef StaggeredForwardSensitivityModelEvaluator<Scalar> FSAME;

  Teuchos::RCP<const WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
      forwardModel_;
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > appExplicitModel_;
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > appImplicitModel_;
  Teuchos::RCP<FSAME> fsaExplicitModel_;
  Teuchos::RCP<FSAME> fsaImplicitModel_;

  bool use_dfdp_as_tangent_;
  int y_tangent_index_;

  Teuchos::RCP<const DMVPVS> explicit_dydp_prod_space_;
  Teuchos::RCP<const DMVPVS> imex_dxdp_prod_space_;

  Teuchos::RCP<const Tempus::SolutionHistory<Scalar> > sh_;
  mutable Scalar t_interp_;
  mutable Teuchos::RCP<const Tempus::SolutionState<Scalar> > forward_state_;
  mutable Teuchos::RCP<Tempus::SolutionState<Scalar> > nc_forward_state_;
  mutable Teuchos::RCP<const Tempus::SolutionState<Scalar> > explicit_y_state_;
  mutable Teuchos::RCP<const Tempus::SolutionState<Scalar> > implicit_x_state_;

  mutable Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > my_dfdp_mv_;
  mutable Teuchos::RCP<Thyra::LinearOpBase<Scalar> > my_dfdp_op_;
};

}  // namespace Tempus

#endif  // Tempus_ModelEvaluatorPairPartIMEX_Basic_decl_hpp
