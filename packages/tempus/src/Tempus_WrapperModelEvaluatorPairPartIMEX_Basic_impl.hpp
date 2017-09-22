// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_ModelEvaluatorPairPartIMEX_Basic_impl_hpp
#define Tempus_ModelEvaluatorPairPartIMEX_Basic_impl_hpp

#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"


namespace Tempus {

template <typename Scalar>
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::
WrapperModelEvaluatorPairPartIMEX_Basic(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& explicitModel,
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& implicitModel,
  int numExplicitOnlyBlocks, int parameterIndex)
  : timeDer_(Teuchos::null), numExplicitOnlyBlocks_(numExplicitOnlyBlocks),
    parameterIndex_(parameterIndex), useImplicitModel_(false)
{
  setExplicitModel(explicitModel);
  setImplicitModel(implicitModel);
  setParameterIndex();
  initialize();
}

template <typename Scalar>
void
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::
initialize()
{
  using Teuchos::RCP;

  useImplicitModel_ = true;
  wrapperImplicitInArgs_  = this->createInArgs();
  wrapperImplicitOutArgs_ = this->createOutArgs();
  useImplicitModel_ = false;

  RCP<const Thyra::VectorBase<Scalar> > z =
    explicitModel_->getNominalValues().get_x();

  // A Thyra::VectorSpace requirement
  TEUCHOS_TEST_FOR_EXCEPTION( !(getIMEXVector(z)->space()->isCompatible(
                                       *(implicitModel_->get_x_space()))),
    std::logic_error,
    "Error - WrapperModelEvaluatorPairIMEX_Basic::initialize()\n"
    "  Explicit and Implicit vector x spaces are incompatible!\n"
    "  Explicit vector x space = " << *(getIMEXVector(z)->space())    << "\n"
    "  Implicit vector x space = " << *(implicitModel_->get_x_space())<< "\n");

  // Valid number of blocks?
  RCP<const Thyra::ProductVectorBase<Scalar> > zPVector =
   Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Scalar> >(z);
   TEUCHOS_TEST_FOR_EXCEPTION( zPVector == Teuchos::null, std::logic_error,
    "Error - WrapperModelEvaluatorPairPartIMEX_Basic::initialize()\n"
    "  was given a VectorBase that could not be cast to a\n"
    "  ProductVectorBase!\n");

  int numBlocks = zPVector->productSpace()->numBlocks();

  TEUCHOS_TEST_FOR_EXCEPTION( !(0 <= numExplicitOnlyBlocks_ and
                                   numExplicitOnlyBlocks_ < numBlocks),
    std::logic_error,
    "Error - WrapperModelEvaluatorPairPartIMEX_Basic::initialize()\n"
    "Invalid number of explicit-only blocks = "<<numExplicitOnlyBlocks_<<"\n"
    "Should be in the interval [0, numBlocks) = [0, "<<numBlocks<<")\n");
}

template <typename Scalar>
void
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::
setAppModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > & me)
{
  TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
    "Error - WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::setAppModel\n"
    "  should not be used.  One should instead use setExplicitModel,\n"
    "  setImplicitModel, or create a new WrapperModelEvaluatorPairPartIMEX.\n");
}

template <typename Scalar>
Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::
getAppModel() const
{
  TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
    "Error - WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::getAppModel\n"
    "  should not be used.  One should instead use getExplicitModel,\n"
    "  getImplicitModel, or directly use WrapperModelEvaluatorPairPartIMEX\n"
    "  object.\n");
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::
get_x_space() const
{
  if (useImplicitModel_ == true) return implicitModel_->get_x_space();

  return explicitModel_->get_x_space();
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::
get_g_space(int i) const
{
  if (useImplicitModel_ == true) return implicitModel_->get_g_space(i);

  return explicitModel_->get_g_space(i);
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::
get_p_space(int i) const
{
  if (useImplicitModel_ == true) return implicitModel_->get_p_space(i);

  return explicitModel_->get_p_space(i);
}

template <typename Scalar>
void
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::
setImplicitModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > & model )
{
  implicitModel_ = model;
}

template <typename Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::
getIMEXVector(const Teuchos::RCP<Thyra::VectorBase<Scalar> > & full) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  if(numExplicitOnlyBlocks_==0)
    return full;

  RCP<Thyra::ProductVectorBase<Scalar> > blk_full =
    rcp_dynamic_cast<Thyra::ProductVectorBase<Scalar> >(full);
  TEUCHOS_TEST_FOR_EXCEPTION( blk_full == Teuchos::null, std::logic_error,
    "Error - WrapperModelEvaluatorPairPartIMEX_Basic::getIMEXVector()\n"
    "  was given a VectorBase that could not be cast to a\n"
    "  ProductVectorBase!\n");
  int numBlocks = blk_full->productSpace()->numBlocks();

  // special case where the implicit terms are not blocked
  if(numBlocks==numExplicitOnlyBlocks_+1)
    return blk_full->getNonconstVectorBlock(numExplicitOnlyBlocks_);

  TEUCHOS_ASSERT(false);
  return Teuchos::null;
}

template <typename Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::
getExplicitOnlyVector(
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & full) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  if(numExplicitOnlyBlocks_==0)
    return Teuchos::null;

  RCP<Thyra::ProductVectorBase<Scalar> > blk_full =
    rcp_dynamic_cast<Thyra::ProductVectorBase<Scalar> >(full);
  TEUCHOS_TEST_FOR_EXCEPTION( blk_full == Teuchos::null, std::logic_error,
    "Error - WrapperModelEvaluatorPairPartIMEX_Basic::getExplicitOnlyVector()\n"
    "  was given a VectorBase that could not be cast to a ProductVectorBase!\n"
    "  full = " << *full << "\n");

  // special case where the explicit terms are not blocked
  if(numExplicitOnlyBlocks_==1)
    return blk_full->getNonconstVectorBlock(0);

  TEUCHOS_ASSERT(false);
  return Teuchos::null;

}

template <typename Scalar>
void
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::
setParameterIndex(int parameterIndex)
{
  if (implicitModel_->Np() == 0) {
    TEUCHOS_TEST_FOR_EXCEPTION( 0 <= parameterIndex_, std::logic_error,
      "Error - WrapperModelEvaluatorPairPartIMEX_Basic::setParameterIndex()\n"
      "  Invalid parameter index = " << parameterIndex_ << "\n"
      "  Should not be set since Np = "<<implicitModel_->Np()<<")\n");
  } else {
    if(parameterIndex >= 0) {
      parameterIndex_ = parameterIndex;
    } else if (parameterIndex_ < 0) {
      parameterIndex_ = 0;
      for(int i=0; i<implicitModel_->Np(); i++) {
        if((*implicitModel_->get_p_names(i))[0] == "EXPLICIT_ONLY_VECTOR") {
          parameterIndex_ = i;
          break;
        }
      }
    }

    TEUCHOS_TEST_FOR_EXCEPTION( 0 > parameterIndex_ or
                                    parameterIndex_ >= implicitModel_->Np(),
     std::logic_error,
     "Error - WrapperModelEvaluatorPairPartIMEX_Basic::setParameterIndex()\n"
     "  Invalid parameter index = " << parameterIndex_ << "\n"
     "  Should be in the interval [0, Np) = [0, "<<implicitModel_->Np()<<")\n");
  }

  return;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::
get_f_space() const
{
  if (useImplicitModel_ == true) return implicitModel_->get_f_space();

  return explicitModel_->get_f_space();
}


template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::
getNominalValues() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgsSetup<Scalar> inArgs = this->createInArgs();
  return inArgs;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::
createInArgs() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  if (useImplicitModel_ == true) {
    MEB::InArgsSetup<Scalar> inArgs(implicitModel_->getNominalValues());
    inArgs.setModelEvalDescription(this->description());
    return inArgs;
  }

  MEB::InArgsSetup<Scalar> inArgs(explicitModel_->getNominalValues());
  inArgs.setModelEvalDescription(this->description());
  return inArgs;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::
createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  if (useImplicitModel_ == true) {
    MEB::OutArgsSetup<Scalar> outArgs(implicitModel_->createOutArgs());
    outArgs.setModelEvalDescription(this->description());
    return outArgs;
  }

  MEB::OutArgsSetup<Scalar> outArgs(explicitModel_->createOutArgs());
  outArgs.setModelEvalDescription(this->description());
  return outArgs;
}

template <typename Scalar>
void WrapperModelEvaluatorPairPartIMEX_Basic<Scalar>::
evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar>  & inArgs,
              const Thyra::ModelEvaluatorBase::OutArgs<Scalar> & outArgs) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::RCP;

  RCP<const Thyra::VectorBase<Scalar> > x = inArgs.get_x();
  RCP<Thyra::VectorBase<Scalar> >   x_dot =
    Thyra::createMember(implicitModel_->get_x_space());
  timeDer_->compute(x, x_dot);

  MEB::InArgs<Scalar>  appImplicitInArgs (wrapperImplicitInArgs_);
  MEB::OutArgs<Scalar> appImplicitOutArgs(wrapperImplicitOutArgs_);
  appImplicitInArgs.set_x(x);
  appImplicitInArgs.set_x_dot(x_dot);
  for (int i=0; i<implicitModel_->Np(); ++i) {
    // Copy over parameters except for the parameter for explicit-only vector!
    if ((inArgs.get_p(i) != Teuchos::null) and (i != parameterIndex_))
      appImplicitInArgs.set_p(i, inArgs.get_p(i));
  }

  appImplicitOutArgs.set_f(outArgs.get_f());
  appImplicitOutArgs.set_W_op(outArgs.get_W_op());

  implicitModel_->evalModel(appImplicitInArgs,appImplicitOutArgs);
}

} // end namespace Tempus

#endif // Tempus_ModelEvaluatorPairPartIMEX_Basic_impl_hpp
