// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_MODEL_EVALUATOR_IMPL_HPP
#define PANZER_MODEL_EVALUATOR_IMPL_HPP

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_ArrayRCP.hpp"

#include "Panzer_config.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_InArgs.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_ThyraObjContainer.hpp"
#include "Panzer_ThyraObjFactory.hpp"

#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"

#include "Tpetra_CrsMatrix.hpp"

// Constructors/Initializers/Accessors

template<typename Scalar, typename NODE>
panzer::ModelEvaluator<Scalar,NODE>::
ModelEvaluator(const Teuchos::RCP<panzer::FieldManagerBuilder>& fmb,
               const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> >& rLibrary,
               const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> >& lof,
               const std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >& p_names,
               const Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > & solverFactory,
               const Teuchos::RCP<panzer::GlobalData>& global_data,
               bool build_transient_support,
               double t_init)
  : t_init_(t_init)
  , fmb_(fmb)
  , p_names_(p_names)
  , responseLibrary_(rLibrary)
  , global_data_(global_data)
  , build_transient_support_(build_transient_support)
  , lof_(lof)
  , solverFactory_(solverFactory)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::tuple;
  using Thyra::VectorBase;
  using Thyra::createMember;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  TEUCHOS_ASSERT(lof_!=Teuchos::null);

  panzer::AssemblyEngine_TemplateBuilder builder(fmb,lof);
  ae_tm_.buildObjects(builder);

  //
  // Setup parameters
  //
  // this->initializeParameterVector(p_names_,global_data->pl);

  //
  // Setup responses
  //

  //
  // Create the structure for the problem
  //

  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_x);
  if(build_transient_support_) {
    inArgs.setSupports(MEB::IN_ARG_x_dot,true);
    inArgs.setSupports(MEB::IN_ARG_t,true);
    inArgs.setSupports(MEB::IN_ARG_alpha,true);
    inArgs.setSupports(MEB::IN_ARG_beta,true);
  }
  prototypeInArgs_ = inArgs;
  
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(p_init_.size(), g_space_.size());
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);
  prototypeOutArgs_ = outArgs;

  //
  // Build x, f spaces
  //
  
  // dynamic cast to blocked LOF for now
  RCP<const ThyraObjFactory<Scalar> > tof = rcp_dynamic_cast<const ThyraObjFactory<Scalar> >(lof,true);

  x_space_ = tof->getThyraDomainSpace();
  f_space_ = tof->getThyraRangeSpace();

  //
  // Setup nominal values
  //

  MEB::InArgsSetup<Scalar> nomInArgs;
  nomInArgs.setModelEvalDescription(this->description());
  nomInArgs.setSupports(MEB::IN_ARG_x);
  Teuchos::RCP<Thyra::VectorBase<double> > x_nom = Thyra::createMember(x_space_);
  Thyra::assign(x_nom.ptr(),0.0);
  nomInArgs.set_x(x_nom);
  if(build_transient_support_) {
    nomInArgs.setSupports(MEB::IN_ARG_x_dot,true);
    nomInArgs.setSupports(MEB::IN_ARG_t,true);
    nomInArgs.setSupports(MEB::IN_ARG_alpha,true);
    nomInArgs.setSupports(MEB::IN_ARG_beta,true);

    Teuchos::RCP<Thyra::VectorBase<double> > x_dot_nom = Thyra::createMember(x_space_);
    Thyra::assign(x_dot_nom.ptr(),0.0);
    nomInArgs.set_x_dot(x_dot_nom);
    nomInArgs.set_t(t_init_);
    nomInArgs.set_alpha(0.0); // these have no meaning initially!
    nomInArgs.set_beta(0.0);
  }
  nominalValues_ = nomInArgs;
  
}

template<typename Scalar, typename NODE>
panzer::ModelEvaluator<Scalar,NODE>::
ModelEvaluator()
{
  TEUCHOS_ASSERT(false);
}

// Public functions overridden from ModelEvaulator

template<typename Scalar, typename NODE>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
panzer::ModelEvaluator<Scalar,NODE>::get_x_space() const
{
  return x_space_;
}


template<typename Scalar, typename NODE>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
panzer::ModelEvaluator<Scalar,NODE>::get_f_space() const
{
  return f_space_;
}

template<typename Scalar, typename NODE>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
panzer::ModelEvaluator<Scalar,NODE>::get_g_space(int i) const
{
  TEUCHOS_ASSERT(i>=0 && 
		 static_cast<typename std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > >::size_type>(i)<g_space_.size());

  return g_space_[i];
}

template<typename Scalar, typename NODE>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
panzer::ModelEvaluator<Scalar,NODE>::createInArgs() const
{
  return prototypeInArgs_;
}

template<typename Scalar,  typename NODE>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
panzer::ModelEvaluator<Scalar,NODE>::getNominalValues() const
{
  return nominalValues_;
}

// Private functions overridden from ModelEvaulatorDefaultBase


template <typename Scalar, typename NODE>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
panzer::ModelEvaluator<Scalar,NODE>::createOutArgsImpl() const
{
  return prototypeOutArgs_;
}

template <typename Scalar, typename NODE>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> > 
panzer::ModelEvaluator<Scalar,NODE>::
create_W_op() const
{
  Teuchos::RCP<const ThyraObjFactory<Scalar> > tof 
     = Teuchos::rcp_dynamic_cast<const ThyraObjFactory<Scalar> >(lof_,true);

  return tof->getThyraMatrix();
}

template <typename Scalar, typename NODE>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > 
panzer::ModelEvaluator<Scalar,NODE>::
get_W_factory() const
{
  return solverFactory_;
}

template <typename Scalar, typename NODE>
void panzer::ModelEvaluator<Scalar,NODE>::
addNonParameterGlobalEvaluationData(const std::string & key,
                                    const Teuchos::RCP<GlobalEvaluationData> & ged)
{
   nonParamGlobalEvaluationData_.addDataObject(key,ged);
}

template <typename Scalar, typename NODE>
void panzer::ModelEvaluator<Scalar,NODE>::
evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
              const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
   evalModelImpl_basic(inArgs,outArgs); 
}

template <typename Scalar, typename NODE>
void panzer::ModelEvaluator<Scalar,NODE>::
evalModelImpl_basic(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
                    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::Array;
  using Teuchos::tuple;
  using Teuchos::rcp_dynamic_cast;

  typedef Thyra::ModelEvaluatorBase MEB;

  // Transient or steady-state evaluation is determined by the x_dot
  // vector.  If this RCP is null, then we are doing a steady-state
  // fill.
  bool is_transient = false;
  if (inArgs.supports(MEB::IN_ARG_x_dot ))
    is_transient = !Teuchos::is_null(inArgs.get_x_dot());

  // Make sure construction built in transient support
  TEUCHOS_TEST_FOR_EXCEPTION(is_transient && !build_transient_support_, std::runtime_error,
		     "ModelEvaluator was not built with transient support enabled!");

  //
  // Get the output arguments
  //
  const RCP<Thyra::VectorBase<Scalar> > f_out = outArgs.get_f();
  const RCP<Thyra::LinearOpBase<Scalar> > W_out = outArgs.get_W_op();
  bool requiredResponses = required_basic_g(outArgs);

  // see if the user wants us to do anything
  if(Teuchos::is_null(f_out) && Teuchos::is_null(W_out) && !requiredResponses) {
     return;
  }

  // the user requested work from this method
  // keep on moving

  // if neccessary build a ghosted container
  if(Teuchos::is_null(ghostedContainer_)) {
     ghostedContainer_ = lof_->buildGhostedLinearObjContainer();
     lof_->initializeGhostedContainer(panzer::LinearObjContainer::X |
                                      panzer::LinearObjContainer::DxDt |
                                      panzer::LinearObjContainer::F |
                                      panzer::LinearObjContainer::Mat, *ghostedContainer_); 
  }

  //
  // Get the input arguments
  //
  RCP<const Thyra::VectorBase<Scalar> > x_dot; // possibly empty, but otherwise uses x_dot
  const RCP<const Thyra::VectorBase<Scalar> > x = inArgs.get_x();
  panzer::AssemblyEngineInArgs ae_inargs;
  ae_inargs.container_ = lof_->buildLinearObjContainer(); // we use a new global container
  ae_inargs.ghostedContainer_ = ghostedContainer_;        // we can reuse the ghosted container
  ae_inargs.alpha = 0.0;
  ae_inargs.beta = 1.0;
  ae_inargs.evaluate_transient_terms = false;
  if (is_transient) {
    x_dot = inArgs.get_x_dot();
    ae_inargs.alpha = inArgs.get_alpha();
    ae_inargs.beta = inArgs.get_beta();
    ae_inargs.time = inArgs.get_t();
    ae_inargs.evaluate_transient_terms = true;
  }
  ae_inargs.addGlobalEvaluationData(nonParamGlobalEvaluationData_);
  
  // here we are building a container, this operation is fast, simply allocating a struct
  const RCP<panzer::ThyraObjContainer<double> > thGlobalContainer = 
    Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(ae_inargs.container_);

  TEUCHOS_ASSERT(!Teuchos::is_null(thGlobalContainer));

  // Ghosted container objects are zeroed out below only if needed for
  // a particular calculation.  This makes it more efficient than
  // zeroing out all objects in the container here.
  const RCP<panzer::ThyraObjContainer<double> > thGhostedContainer = 
    Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(ae_inargs.ghostedContainer_);
  
  // Set the solution vector (currently all targets require solution).
  // In the future we may move these into the individual cases below.
  // A very subtle (and fragile) point: A non-null pointer in global
  // container triggers export operations during fill.  Also, the
  // introduction of the container is forcing us to cast away const on
  // arguments that should be const.  Another reason to redesign
  // LinearObjContainer layers.
  thGlobalContainer->set_x_th(Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(x));
  if (is_transient)
    thGlobalContainer->set_dxdt_th(Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(x_dot));

  if (!Teuchos::is_null(f_out) && !Teuchos::is_null(W_out)) {

    PANZER_FUNC_TIME_MONITOR("panzer::ModelEvaluator::evalModel(f and J)");

    // Set the targets
    thGlobalContainer->set_f_th(f_out);
    thGlobalContainer->set_A_th(W_out);

    // Zero values in ghosted container objects
    Thyra::assign(thGhostedContainer->get_f_th().ptr(),0.0);
    thGhostedContainer->initializeMatrix(0.0);

    ae_tm_.template getAsObject<panzer::Traits::Jacobian>()->evaluate(ae_inargs);
  }
  else if(!Teuchos::is_null(f_out) && Teuchos::is_null(W_out)) {

    PANZER_FUNC_TIME_MONITOR("panzer::ModelEvaluator::evalModel(f)");

    thGlobalContainer->set_f_th(f_out);

    // Zero values in ghosted container objects
    Thyra::assign(thGhostedContainer->get_f_th().ptr(),0.0);

    ae_tm_.template getAsObject<panzer::Traits::Residual>()->evaluate(ae_inargs);
  }
  else if(Teuchos::is_null(f_out) && !Teuchos::is_null(W_out)) {

    PANZER_FUNC_TIME_MONITOR("panzer::ModelEvaluator::evalModel(J)");

    // this dummy nonsense is needed only for scattering dirichlet conditions
    RCP<Thyra::VectorBase<Scalar> > dummy_f = Thyra::createMember(f_space_);
    thGlobalContainer->set_f_th(dummy_f); 
    thGlobalContainer->set_A_th(W_out);

    // Zero values in ghosted container objects
    thGhostedContainer->initializeMatrix(0.0);

    ae_tm_.template getAsObject<panzer::Traits::Jacobian>()->evaluate(ae_inargs);
  }

  // evaluate responses...uses the stored assembly arguments and containers
  if(requiredResponses)
     evalModelImpl_basic_g(ae_inargs,inArgs,outArgs);

  // Holding a rcp to f produces a seg fault in Rythmos when the next
  // f comes in and the resulting dtor is called.  Need to discuss
  // with Ross.  Clearing all references here works!

  thGlobalContainer->set_x_th(Teuchos::null);
  thGlobalContainer->set_dxdt_th(Teuchos::null);
  thGlobalContainer->set_f_th(Teuchos::null);
  thGlobalContainer->set_A_th(Teuchos::null);

  // forget previous containers
  ae_inargs.container_ = Teuchos::null;
  ae_inargs.ghostedContainer_ = Teuchos::null;
}

template <typename Scalar, typename NODE>
void panzer::ModelEvaluator<Scalar,NODE>::
evalModelImpl_basic_g(panzer::AssemblyEngineInArgs & ae_inargs,
                      const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
                      const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
   // optional sanity check
   // TEUCHOS_ASSERT(required_basic_g(outArgs));

   for(std::size_t i=0;i<g_names_.size();i++) {
      Teuchos::RCP<Thyra::VectorBase<double> > vec = outArgs.get_g(i);
      if(vec!=Teuchos::null) {
        std::string responseName = g_names_[i];
        Teuchos::RCP<panzer::ResponseMESupportBase<panzer::Traits::Residual> > resp 
            = Teuchos::rcp_dynamic_cast<panzer::ResponseMESupportBase<panzer::Traits::Residual> >(responseLibrary_->getResponse<panzer::Traits::Residual>(responseName));
        resp->setVector(vec);
      }
   }

   // evaluator responses
   responseLibrary_->addResponsesToInArgs<panzer::Traits::Residual>(ae_inargs);
   responseLibrary_->evaluate<panzer::Traits::Residual>(ae_inargs);
}

template <typename Scalar, typename NODE>
bool panzer::ModelEvaluator<Scalar,NODE>::
required_basic_g(const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
   // determine if any of the outArgs are not null!
   bool activeGArgs = false;
   for(int i=0;i<outArgs.Ng();i++) 
      activeGArgs |= (outArgs.get_g(i)!=Teuchos::null); 

   return activeGArgs;
}

#endif
