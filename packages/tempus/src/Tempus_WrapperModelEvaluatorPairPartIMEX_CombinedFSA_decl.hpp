//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_ModelEvaluatorPairPartIMEX_CombinedFSA_decl_hpp
#define Tempus_ModelEvaluatorPairPartIMEX_CombinedFSA_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_SensitivityModelEvaluatorBase.hpp"
#include "Tempus_WrapperModelEvaluatorPairPartIMEX_Basic.hpp"
#include "Tempus_CombinedForwardSensitivityModelEvaluator.hpp"

#include "Thyra_ProductMultiVectorBase.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"

namespace Tempus {

/** \brief Specialization of IMEX-Part ME for "combined" FSA method.
 *
 * This specializes the implementation of several parts of
 * WrapperModelEvaluatorPairPartIMEX_Basic for forward-sensitivity analysis
 * with CombinedForwardSensitivityModelEvaluator.  It deals with the product
 * structure of the sensitivity solution vectors and incorporates the
 * sensitivity of the implicit term with respect to the explicit-only vector.
 */
template <typename Scalar>
class WrapperModelEvaluatorPairPartIMEX_CombinedFSA
  : public SensitivityModelEvaluatorBase<Scalar>,
    public WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> {
 public:
  /// Constructor
  WrapperModelEvaluatorPairPartIMEX_CombinedFSA(
      const Teuchos::RCP<
          const WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >& forwardModel,
      const Teuchos::RCP<const Teuchos::ParameterList>& pList = Teuchos::null);

  /// Destructor
  virtual ~WrapperModelEvaluatorPairPartIMEX_CombinedFSA() {}

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
  typedef WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> Base;
  typedef Thyra::DefaultMultiVectorProductVectorSpace<Scalar> DMVPVS;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
  typedef Thyra::ProductMultiVectorBase<Scalar> PMVB;
  typedef CombinedForwardSensitivityModelEvaluator<Scalar> FSAME;

  Teuchos::RCP<const WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
      forwardModel_;
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > appExplicitModel_;
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > appImplicitModel_;
  Teuchos::RCP<const FSAME> fsaExplicitModel_;
  Teuchos::RCP<const FSAME> fsaImplicitModel_;

  bool use_dfdp_as_tangent_;
  int y_tangent_index_;

  Teuchos::RCP<const DMVPVS> explicit_y_dydp_prod_space_;
  Teuchos::RCP<const DMVPVS> explicit_dydp_prod_space_;
  Teuchos::RCP<const DMVPVS> imex_x_dxdp_prod_space_;

  mutable Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > my_dfdp_mv_;
  mutable Teuchos::RCP<Thyra::LinearOpBase<Scalar> > my_dfdp_op_;
};

}  // namespace Tempus

#endif  // Tempus_ModelEvaluatorPairPartIMEX_Basic_decl_hpp
