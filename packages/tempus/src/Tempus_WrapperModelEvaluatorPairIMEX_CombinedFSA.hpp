//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_ModelEvaluatorPairIMEX_CombinedFSA_hpp
#define Tempus_ModelEvaluatorPairIMEX_CombinedFSA_hpp

#include "Tempus_config.hpp"
#include "Tempus_SensitivityModelEvaluatorBase.hpp"
#include "Tempus_WrapperModelEvaluatorPairIMEX_Basic.hpp"
#include "Tempus_CombinedForwardSensitivityModelEvaluator.hpp"

namespace Tempus {

/** \brief Specialization of IMEX ME for "combined" FSA method.
 *
 * For the combined forward sensitivitymethod, the implementation found in
 * WrapperModelEvaluatorPairIMEX_Basic works just fine.  We go ahead and
 * create a specialized class to follow the pattern of other methods and also
 * handle the wrapping of the underlying MEs.
 */
template <typename Scalar>
class WrapperModelEvaluatorPairIMEX_CombinedFSA
  : public SensitivityModelEvaluatorBase<Scalar>,
    public WrapperModelEvaluatorPairIMEX_Basic<Scalar> {
 public:
  /// Constructor
  WrapperModelEvaluatorPairIMEX_CombinedFSA(
      const Teuchos::RCP<const WrapperModelEvaluatorPairIMEX_Basic<Scalar> >&
          forwardModel,
      const Teuchos::RCP<const Teuchos::ParameterList>& pList = Teuchos::null)
  {
    forwardModel_     = forwardModel;
    appExplicitModel_ = forwardModel_->getExplicitModel();
    appImplicitModel_ = forwardModel_->getImplicitModel();
    fsaExplicitModel_ = rcp(new FSAME(appExplicitModel_, appExplicitModel_,
                                      appExplicitModel_, pList));
    fsaImplicitModel_ = rcp(new FSAME(appImplicitModel_, appImplicitModel_,
                                      appImplicitModel_, pList));
    Base::setup(fsaExplicitModel_, fsaImplicitModel_);
  }

  /// Destructor
  virtual ~WrapperModelEvaluatorPairIMEX_CombinedFSA() {}

  /// \name Overridden from Tempus::SensitivityModelEvaluatorBase
  //@{

  /// Get the underlying forward model
  virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getForwardModel()
      const
  {
    return forwardModel_;
  }

  //@}

 private:
  /// Default constructor - not allowed
  WrapperModelEvaluatorPairIMEX_CombinedFSA() {}

 protected:
  typedef WrapperModelEvaluatorPairIMEX_Basic<Scalar> Base;
  typedef CombinedForwardSensitivityModelEvaluator<Scalar> FSAME;

  Teuchos::RCP<const WrapperModelEvaluatorPairIMEX_Basic<Scalar> >
      forwardModel_;
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > appExplicitModel_;
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > appImplicitModel_;
  Teuchos::RCP<FSAME> fsaExplicitModel_;
  Teuchos::RCP<FSAME> fsaImplicitModel_;
};

}  // namespace Tempus

#endif  // Tempus_ModelEvaluatorPairIMEX_CombinedFSA_hpp
