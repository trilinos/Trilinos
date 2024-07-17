//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_ModelEvaluatorPairPartIMEX_Basic_impl_hpp
#define Tempus_ModelEvaluatorPairPartIMEX_Basic_impl_hpp

#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"

namespace Tempus {

template <typename Scalar>
WrapperModelEvaluatorPairPartIMEX_Basic<
    Scalar>::WrapperModelEvaluatorPairPartIMEX_Basic()
  : timeDer_(Teuchos::null),
    numExplicitOnlyBlocks_(0),
    parameterIndex_(-1),
    useImplicitModel_(false)
{
}

template <typename Scalar>
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::
    WrapperModelEvaluatorPairPartIMEX_Basic(
        const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& explicitModel,
        const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& implicitModel,
        int numExplicitOnlyBlocks, int parameterIndex)
  : timeDer_(Teuchos::null),
    numExplicitOnlyBlocks_(numExplicitOnlyBlocks),
    parameterIndex_(parameterIndex),
    useImplicitModel_(false)
{
  setExplicitModel(explicitModel);
  setImplicitModel(implicitModel);
  setParameterIndex();
  initialize();
}

template <typename Scalar>
void WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::setup(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& explicitModel,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& implicitModel,
    int numExplicitOnlyBlocks, int parameterIndex)
{
  timeDer_               = Teuchos::null;
  numExplicitOnlyBlocks_ = numExplicitOnlyBlocks;
  parameterIndex_        = parameterIndex;
  useImplicitModel_      = false;
  setExplicitModel(explicitModel);
  setImplicitModel(implicitModel);
  setParameterIndex();
  initialize();
}

template <typename Scalar>
void WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::initialize()
{
  using Teuchos::RCP;

  useImplicitModel_       = true;
  wrapperImplicitInArgs_  = this->createInArgs();
  wrapperImplicitOutArgs_ = this->createOutArgs();
  useImplicitModel_       = false;

  RCP<const Thyra::VectorBase<Scalar> > z =
      explicitModel_->getNominalValues().get_x();

  // A Thyra::VectorSpace requirement
  TEUCHOS_TEST_FOR_EXCEPTION(
      !(getIMEXVector(z)->space()->isCompatible(
          *(implicitModel_->get_x_space()))),
      std::logic_error,
      "Error - WrapperModelEvaluatorPairIMEX_Basic::initialize()\n"
          << "  Explicit and Implicit vector x spaces are incompatible!\n"
          << "  Explicit vector x space = "
          << *(getIMEXVector(z)->space())
          << "\n  Implicit vector x space = "
          << *(implicitModel_->get_x_space()) << "\n");

  // Valid number of blocks?
  RCP<const Thyra::ProductVectorBase<Scalar> > zPVector =
      Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Scalar> >(z);
  TEUCHOS_TEST_FOR_EXCEPTION(
      zPVector == Teuchos::null, std::logic_error,
      "Error - WrapperModelEvaluatorPairPartIMEX_Basic::initialize()\n"
      "  was given a VectorBase that could not be cast to a\n"
      "  ProductVectorBase!\n");

  int numBlocks = zPVector->productSpace()->numBlocks();

  TEUCHOS_TEST_FOR_EXCEPTION(
      !(0 <= numExplicitOnlyBlocks_ && numExplicitOnlyBlocks_ < numBlocks),
      std::logic_error,
      "Error - WrapperModelEvaluatorPairPartIMEX_Basic::initialize()\n"
          << "Invalid number of explicit-only blocks = "
          << numExplicitOnlyBlocks_
          << "\nShould be in the interval [0, numBlocks) = [0, "
          << numBlocks << ")\n");
}

template <typename Scalar>
void WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::setAppModel(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& /* me */)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Error - WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::setAppModel\n"
      "  should not be used.  One should instead use setExplicitModel,\n"
      "  setImplicitModel, or create a new "
      "WrapperModelEvaluatorPairPartIMEX.\n");
}

template <typename Scalar>
Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::getAppModel() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Error - WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::getAppModel\n"
      "  should not be used.  One should instead use getExplicitModel,\n"
      "  getImplicitModel, or directly use WrapperModelEvaluatorPairPartIMEX\n"
      "  object.\n");
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::get_x_space() const
{
  if (useImplicitModel_ == true) return implicitModel_->get_x_space();

  return explicitModel_->get_x_space();
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::get_g_space(int i) const
{
  if (useImplicitModel_ == true) return implicitModel_->get_g_space(i);

  return explicitModel_->get_g_space(i);
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::get_p_space(int i) const
{
  if (useImplicitModel_ == true) return implicitModel_->get_p_space(i);

  return explicitModel_->get_p_space(i);
}

template <typename Scalar>
void WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::setImplicitModel(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model)
{
  implicitModel_ = model;
}

template <typename Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::getIMEXVector(
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& full) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  Teuchos::RCP<Thyra::VectorBase<Scalar> > vector;
  if (full == Teuchos::null) {
    vector = Teuchos::null;
  }
  else if (numExplicitOnlyBlocks_ == 0) {
    vector = full;
  }
  else {
    RCP<Thyra::ProductVectorBase<Scalar> > blk_full =
        rcp_dynamic_cast<Thyra::ProductVectorBase<Scalar> >(full);
    TEUCHOS_TEST_FOR_EXCEPTION(
        blk_full == Teuchos::null, std::logic_error,
        "Error - WrapperModelEvaluatorPairPartIMEX_Basic::getIMEXVector()\n"
        "  was given a VectorBase that could not be cast to a\n"
        "  ProductVectorBase!\n");
    int numBlocks = blk_full->productSpace()->numBlocks();

    // special case where the implicit terms are not blocked
    if (numBlocks == numExplicitOnlyBlocks_ + 1)
      vector = blk_full->getNonconstVectorBlock(numExplicitOnlyBlocks_);

    TEUCHOS_TEST_FOR_EXCEPTION(
        !(numExplicitOnlyBlocks_ == 0 || full == Teuchos::null ||
          numBlocks == numExplicitOnlyBlocks_ + 1),
        std::logic_error, "Error - Invalid values!\n");
  }

  return vector;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::getIMEXVector(
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& full) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  Teuchos::RCP<const Thyra::VectorBase<Scalar> > vector;
  if (full == Teuchos::null) {
    vector = Teuchos::null;
  }
  else if (numExplicitOnlyBlocks_ == 0) {
    vector = full;
  }
  else {
    // special case where the implicit terms are not blocked

    RCP<const Thyra::ProductVectorBase<Scalar> > blk_full =
        rcp_dynamic_cast<const Thyra::ProductVectorBase<Scalar> >(full);
    TEUCHOS_TEST_FOR_EXCEPTION(
        blk_full == Teuchos::null, std::logic_error,
        "Error - WrapperModelEvaluatorPairPartIMEX_Basic::getIMEXVector()\n"
        "  was given a VectorBase that could not be cast to a\n"
        "  ProductVectorBase!\n");
    int numBlocks = blk_full->productSpace()->numBlocks();

    // special case where the implicit terms are not blocked
    if (numBlocks == numExplicitOnlyBlocks_ + 1)
      vector = blk_full->getVectorBlock(numExplicitOnlyBlocks_);

    TEUCHOS_TEST_FOR_EXCEPTION(
        !(numExplicitOnlyBlocks_ == 0 || full == Teuchos::null ||
          numBlocks == numExplicitOnlyBlocks_ + 1),
        std::logic_error, "Error - Invalid values!\n");
  }

  return vector;
}

template <typename Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::getExplicitOnlyVector(
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& full) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  Teuchos::RCP<Thyra::VectorBase<Scalar> > vector;
  if (numExplicitOnlyBlocks_ == 0 || full == Teuchos::null) {
    vector = Teuchos::null;
  }
  else if (numExplicitOnlyBlocks_ == 1) {
    // special case where the explicit terms are not blocked

    RCP<Thyra::ProductVectorBase<Scalar> > blk_full =
        rcp_dynamic_cast<Thyra::ProductVectorBase<Scalar> >(full);
    TEUCHOS_TEST_FOR_EXCEPTION(
        blk_full == Teuchos::null, std::logic_error,
        "Error - WrapperModelEvaluatorPairPartIMEX_Basic::getExplicitOnlyVector\n"
            << "  given a VectorBase that could not be cast to a ProductVectorBase!\n"
            << "  full = " << *full << "\n");

    vector = blk_full->getNonconstVectorBlock(0);
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
      !((numExplicitOnlyBlocks_ == 0 || full == Teuchos::null) ||
        (numExplicitOnlyBlocks_ == 1)),
      std::logic_error, "Error - Invalid values!\n");

  return vector;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::getExplicitOnlyVector(
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& full) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  RCP<const Thyra::VectorBase<Scalar> > vector;
  if (numExplicitOnlyBlocks_ == 0 || full == Teuchos::null) {
    vector = Teuchos::null;
  }
  else if (numExplicitOnlyBlocks_ == 1) {
    // special case where the explicit terms are not blocked

    RCP<const Thyra::ProductVectorBase<Scalar> > blk_full =
        rcp_dynamic_cast<const Thyra::ProductVectorBase<Scalar> >(full);
    TEUCHOS_TEST_FOR_EXCEPTION(
        blk_full == Teuchos::null, std::logic_error,
        "Error - WrapperModelEvaluatorPairPartIMEX_Basic::getExplicitOnlyVector\n"
            << "  given a VectorBase that could not be cast to a ProductVectorBase!\n"
            << "  full = " << *full << "\n");

    vector = blk_full->getVectorBlock(0);
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
      !((numExplicitOnlyBlocks_ == 0 || full == Teuchos::null) ||
        (numExplicitOnlyBlocks_ == 1)),
      std::logic_error, "Error - Invalid values!\n");

  return vector;
}

template <typename Scalar>
void WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::setParameterIndex(
    int parameterIndex)
{
  if (implicitModel_->Np() == 0) {
    if (parameterIndex >= 0) {
      Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
      out->setOutputToRootOnly(0);
      Teuchos::OSTab ostab(
          out, 1,
          "WrapperModelEvaluatorPairPartIMEX_Basic::setParameterIndex()");
      *out << "Warning -- \n"
           << "  Invalid input parameter index = " << parameterIndex << "\n"
           << "  Should not be set since Np = " << implicitModel_->Np() << "\n"
           << "  Not setting parameter index." << std::endl;
    }
    if (parameterIndex_ >= 0) {
      Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
      out->setOutputToRootOnly(0);
      Teuchos::OSTab ostab(
          out, 1,
          "WrapperModelEvaluatorPairPartIMEX_Basic::setParameterIndex()");
      *out << "Warning -- \n"
           << "  Invalid parameter index = " << parameterIndex_ << "\n"
           << "  Should not be set since Np = " << implicitModel_->Np() << "\n"
           << "  Resetting to parameter index to unset state." << std::endl;
      parameterIndex_ = -1;
    }
  }
  else {
    if (parameterIndex >= 0) {
      parameterIndex_ = parameterIndex;
    }
    else if (parameterIndex_ < 0) {
      parameterIndex_ = 0;
      for (int i = 0; i < implicitModel_->Np(); i++) {
        if ((*implicitModel_->get_p_names(i))[0] == "EXPLICIT_ONLY_VECTOR") {
          parameterIndex_ = i;
          break;
        }
      }
    }
  }

  return;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::get_f_space() const
{
  if (useImplicitModel_ == true) return implicitModel_->get_f_space();

  return explicitModel_->get_f_space();
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::getNominalValues() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgsSetup<Scalar> inArgs = this->createInArgs();
  return std::move(inArgs);
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::createInArgs() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar> implicitInArgs = implicitModel_->getNominalValues();
  MEB::InArgs<Scalar> explicitInArgs = explicitModel_->getNominalValues();
  const int np                       = std::max(implicitInArgs.Np(), explicitInArgs.Np());
  if (useImplicitModel_ == true) {
    MEB::InArgsSetup<Scalar> inArgs(implicitInArgs);
    inArgs.setModelEvalDescription(this->description());
    inArgs.set_Np(np);
    return std::move(inArgs);
  }

  MEB::InArgsSetup<Scalar> inArgs(explicitInArgs);
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(np);
  return std::move(inArgs);
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgs<Scalar> implicitOutArgs = implicitModel_->createOutArgs();
  MEB::OutArgs<Scalar> explicitOutArgs = explicitModel_->createOutArgs();
  const int np                         = std::max(implicitOutArgs.Np(), explicitOutArgs.Np());
  if (useImplicitModel_ == true) {
    MEB::OutArgsSetup<Scalar> outArgs(implicitOutArgs);
    outArgs.setModelEvalDescription(this->description());
    outArgs.set_Np_Ng(np, implicitOutArgs.Ng());
    return std::move(outArgs);
  }

  MEB::OutArgsSetup<Scalar> outArgs(explicitOutArgs);
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(np, explicitOutArgs.Ng());
  return std::move(outArgs);
}

template <typename Scalar>
void WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::RCP;

  RCP<const Thyra::VectorBase<Scalar> > x = inArgs.get_x();
  RCP<Thyra::VectorBase<Scalar> > x_dot =
      Thyra::createMember(implicitModel_->get_x_space());
  timeDer_->compute(x, x_dot);

  MEB::InArgs<Scalar> appImplicitInArgs(wrapperImplicitInArgs_);
  MEB::OutArgs<Scalar> appImplicitOutArgs(wrapperImplicitOutArgs_);
  appImplicitInArgs.set_x(x);
  appImplicitInArgs.set_x_dot(x_dot);
  for (int i = 0; i < implicitModel_->Np(); ++i) {
    // Copy over parameters except for the parameter for explicit-only vector!
    if ((inArgs.get_p(i) != Teuchos::null) && (i != parameterIndex_))
      appImplicitInArgs.set_p(i, inArgs.get_p(i));
  }

  appImplicitOutArgs.set_f(outArgs.get_f());
  appImplicitOutArgs.set_W_op(outArgs.get_W_op());

  implicitModel_->evalModel(appImplicitInArgs, appImplicitOutArgs);
}

}  // end namespace Tempus

#endif  // Tempus_ModelEvaluatorPairPartIMEX_Basic_impl_hpp
