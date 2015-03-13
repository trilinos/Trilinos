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
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_ParameterList_GlobalEvaluationData.hpp"
#include "Panzer_EpetraVector_ReadOnly_GlobalEvaluationData.hpp"
#include "Panzer_LinearObjFactory_Utilities.hpp"

#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_DefaultSpmdVector.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Tpetra_CrsMatrix.hpp"

// Constructors/Initializers/Accessors

template<typename Scalar>
panzer::ModelEvaluator<Scalar>::
ModelEvaluator(const Teuchos::RCP<panzer::FieldManagerBuilder>& fmb,
               const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> >& rLibrary,
               const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> >& lof,
               const std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >& p_names,
               const Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > & solverFactory,
               const Teuchos::RCP<panzer::GlobalData>& global_data,
               bool build_transient_support,
               double t_init)
  : t_init_(t_init)
  , require_in_args_refresh_(true)
  , require_out_args_refresh_(true)
  , responseLibrary_(rLibrary)
  , global_data_(global_data)
  , build_transient_support_(build_transient_support)
  , lof_(lof)
  , solverFactory_(solverFactory)
  , oneTimeDirichletBeta_on_(false)
  , oneTimeDirichletBeta_(0.0)
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
  for(std::size_t i=0;i<p_names.size();i++)
     addParameter(*(p_names[i]));

  //
  // Build x, f spaces
  //
  
  // dynamic cast to blocked LOF for now
  RCP<const ThyraObjFactory<Scalar> > tof = rcp_dynamic_cast<const ThyraObjFactory<Scalar> >(lof,true);

  x_space_ = tof->getThyraDomainSpace();
  f_space_ = tof->getThyraRangeSpace();

  // now that the vector spaces are setup we can allocate the nominal values
  // (i.e. initial conditions)
  initializeNominalValues();
}

template<typename Scalar>
panzer::ModelEvaluator<Scalar>::
ModelEvaluator(const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> >& lof,
               const Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > & solverFactory,
               const Teuchos::RCP<panzer::GlobalData>& global_data,
               bool build_transient_support,double t_init)
  : t_init_(t_init)
  , require_in_args_refresh_(true)
  , require_out_args_refresh_(true)
  , global_data_(global_data)
  , build_transient_support_(build_transient_support)
  , lof_(lof)
  , solverFactory_(solverFactory)
  , oneTimeDirichletBeta_on_(false)
  , oneTimeDirichletBeta_(0.0)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  TEUCHOS_ASSERT(lof_!=Teuchos::null);

  //
  // Build x, f spaces
  //
  
  // dynamic cast to blocked LOF for now
  RCP<const ThyraObjFactory<Scalar> > tof = rcp_dynamic_cast<const ThyraObjFactory<Scalar> >(lof_,true);

  x_space_ = tof->getThyraDomainSpace();
  f_space_ = tof->getThyraRangeSpace();

  // now that the vector spaces are setup we can allocate the nominal values
  // (i.e. initial conditions)
  initializeNominalValues();

  // allocate a response library so that responses can be added, it will be initialized in "setupModel"
  responseLibrary_ = Teuchos::rcp(new panzer::ResponseLibrary<panzer::Traits>());
}

template<typename Scalar>
panzer::ModelEvaluator<Scalar>::
ModelEvaluator()
{
  TEUCHOS_ASSERT(false);
}

// Public functions overridden from ModelEvaulator

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
panzer::ModelEvaluator<Scalar>::get_x_space() const
{
  return x_space_;
}


template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
panzer::ModelEvaluator<Scalar>::get_f_space() const
{
  return f_space_;
}

template<typename Scalar>
Teuchos::RCP<const Teuchos::Array<std::string> > 
panzer::ModelEvaluator<Scalar>::get_p_names(int i) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!(i>=0 && i<Teuchos::as<int>(parameters_.size())),std::runtime_error,
                             "panzer::ModelEvaluator::get_p_names: Requested parameter index out of range.");

  return parameters_[i]->names;
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > 
panzer::ModelEvaluator<Scalar>::get_p_space(int i) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!(i>=0 && i<Teuchos::as<int>(parameters_.size())),std::runtime_error,
                             "panzer::ModelEvaluator::get_p_space: Requested parameter index out of range.");

  return parameters_[i]->space;
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
panzer::ModelEvaluator<Scalar>::get_g_space(int i) const
{
  TEUCHOS_ASSERT(i>=0 && 
                 static_cast<typename std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > >::size_type>(i)<g_space_.size());

  return g_space_[i];
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
panzer::ModelEvaluator<Scalar>::createInArgs() const
{
  return getNominalValues();
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
panzer::ModelEvaluator<Scalar>::getNominalValues() const
{
  if(require_in_args_refresh_) {
    typedef Thyra::ModelEvaluatorBase MEB;

    //
    // Refresh nominal values, this is primarily adding
    // new parameters.
    //
  
    MEB::InArgsSetup<Scalar> nomInArgs;
    nomInArgs = nominalValues_;
    nomInArgs.setSupports(nominalValues_);

    // setup parameter support
    nomInArgs.set_Np(parameters_.size());
    for(std::size_t p=0;p<parameters_.size();p++) 
      nomInArgs.set_p(p,parameters_[p]->initial_value);

    nominalValues_ = nomInArgs;
  }

  // refresh no longer required
  require_in_args_refresh_ = false;

  return nominalValues_;
}

template<typename Scalar>
void
panzer::ModelEvaluator<Scalar>::initializeNominalValues() 
{
  typedef Thyra::ModelEvaluatorBase MEB;

  //
  // Setup nominal values
  //
  
  MEB::InArgsSetup<Scalar> nomInArgs;
  nomInArgs.setModelEvalDescription(this->description());
  nomInArgs.setSupports(MEB::IN_ARG_x);
  Teuchos::RCP<Thyra::VectorBase<Scalar> > x_nom = Thyra::createMember(x_space_);
  Thyra::assign(x_nom.ptr(),0.0);
  nomInArgs.set_x(x_nom);
  if(build_transient_support_) {
    nomInArgs.setSupports(MEB::IN_ARG_x_dot,true);
    nomInArgs.setSupports(MEB::IN_ARG_t,true);
    nomInArgs.setSupports(MEB::IN_ARG_alpha,true);
    nomInArgs.setSupports(MEB::IN_ARG_beta,true);

    Teuchos::RCP<Thyra::VectorBase<Scalar> > x_dot_nom = Thyra::createMember(x_space_);
    Thyra::assign(x_dot_nom.ptr(),0.0);
    nomInArgs.set_x_dot(x_dot_nom);
    nomInArgs.set_t(t_init_);
    nomInArgs.set_alpha(0.0); // these have no meaning initially!
    nomInArgs.set_beta(0.0);
  }

  // setup parameter support
  nomInArgs.set_Np(parameters_.size());
  for(std::size_t p=0;p<parameters_.size();p++) 
    nomInArgs.set_p(p,parameters_[p]->initial_value);

  nominalValues_ = nomInArgs;
}

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
setupModel(const Teuchos::RCP<panzer::WorksetContainer> & wc,
           const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
           const std::vector<panzer::BC> & bcs,
           const panzer::EquationSetFactory & eqset_factory,
           const panzer::BCStrategyFactory& bc_factory,
           const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& volume_cm_factory,
           const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& bc_cm_factory,
           const Teuchos::ParameterList& closure_models,
           const Teuchos::ParameterList& user_data,
           bool writeGraph,const std::string & graphPrefix)
{
  // First: build residual assembly engine
  /////////////////////////////////////////////////////////////////////////////////////////////////

  {
    // 1. build Field manager builder
    /////////////////////////////////////////////////////////////////////////////////////////////////

    Teuchos::RCP<panzer::FieldManagerBuilder> fmb = Teuchos::rcp(new panzer::FieldManagerBuilder);
    fmb->setWorksetContainer(wc);
    fmb->setupVolumeFieldManagers(physicsBlocks,volume_cm_factory,closure_models,*lof_,user_data);
    fmb->setupBCFieldManagers(bcs,physicsBlocks,eqset_factory,bc_cm_factory,bc_factory,closure_models,*lof_,user_data);

    // Print Phalanx DAGs
    if (writeGraph){
      fmb->writeVolumeGraphvizDependencyFiles(graphPrefix, physicsBlocks);
      fmb->writeBCGraphvizDependencyFiles(graphPrefix+"BC_");
    }

    panzer::AssemblyEngine_TemplateBuilder builder(fmb,lof_);
    ae_tm_.buildObjects(builder);
  }

  // Second: build the responses
  /////////////////////////////////////////////////////////////////////////////////////////////////

  responseLibrary_->initialize(wc,lof_->getUniqueGlobalIndexerBase(),lof_);

  buildResponses(physicsBlocks,eqset_factory,volume_cm_factory,closure_models,user_data,writeGraph,graphPrefix+"Responses_");
  buildDistroParamDfDp_RL(wc,physicsBlocks,bcs,eqset_factory,bc_factory,volume_cm_factory,closure_models,user_data,writeGraph,graphPrefix+"Response_DfDp_");
}

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
setupAssemblyInArgs(const Thyra::ModelEvaluatorBase::InArgs<Scalar> & inArgs,
                    panzer::AssemblyEngineInArgs & ae_inargs) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  typedef Thyra::ModelEvaluatorBase MEB;

  // if neccessary build a ghosted container
  if(Teuchos::is_null(ghostedContainer_)) {
     ghostedContainer_ = lof_->buildGhostedLinearObjContainer();
     lof_->initializeGhostedContainer(panzer::LinearObjContainer::X |
                                      panzer::LinearObjContainer::DxDt |
                                      panzer::LinearObjContainer::F |
                                      panzer::LinearObjContainer::Mat, *ghostedContainer_); 
  }

  bool is_transient = false;
  if (inArgs.supports(MEB::IN_ARG_x_dot ))
    is_transient = !Teuchos::is_null(inArgs.get_x_dot());

  const RCP<const Thyra::VectorBase<Scalar> > x = inArgs.get_x();
  RCP<const Thyra::VectorBase<Scalar> > x_dot; // possibly empty, but otherwise uses x_dot

  // Make sure construction built in transient support
  TEUCHOS_TEST_FOR_EXCEPTION(is_transient && !build_transient_support_, std::runtime_error,
                     "ModelEvaluator was not built with transient support enabled!");

  ae_inargs.container_ = lof_->buildLinearObjContainer(); // we use a new global container
  ae_inargs.ghostedContainer_ = ghostedContainer_;        // we can reuse the ghosted container
  ae_inargs.alpha = 0.0;
  ae_inargs.beta = 1.0;
  ae_inargs.evaluate_transient_terms = false;
  if (build_transient_support_) {
    x_dot = inArgs.get_x_dot();
    ae_inargs.alpha = inArgs.get_alpha();
    ae_inargs.beta = inArgs.get_beta();
    ae_inargs.time = inArgs.get_t();
    ae_inargs.evaluate_transient_terms = true;
  }

  // this member is handled in the individual functions
  ae_inargs.apply_dirichlet_beta = false;

  // Set input parameters
  for (int i=0; i<inArgs.Np(); i++) {
    
    RCP<const Thyra::VectorBase<Scalar> > p = inArgs.get_p(i);

    if ( p!=Teuchos::null && !parameters_[i]->is_distributed) {
      // non distributed parameters

      Teuchos::ArrayRCP<const Scalar> p_data;
      rcp_dynamic_cast<const Thyra::SpmdVectorBase<Scalar> >(p,true)->getLocalData(Teuchos::ptrFromRef(p_data));

      for (unsigned int j=0; j < parameters_[i]->scalar_value.size(); j++) {
        parameters_[i]->scalar_value[j].baseValue = p_data[j];
        parameters_[i]->scalar_value[j].family->setRealValueForAllTypes(parameters_[i]->scalar_value[j].baseValue);
      }
    }
    else if ( p!=Teuchos::null && parameters_[i]->is_distributed) {
      // distributed parameters

      std::string key = (*parameters_[i]->names)[0];
      RCP<GlobalEvaluationData> ged = distrParamGlobalEvaluationData_.getDataObject(key);

      TEUCHOS_ASSERT(ged!=Teuchos::null);

      // cast to a LOCPair throwing an exception if the cast doesn't work.
      RCP<LOCPair_GlobalEvaluationData> loc_pair_ged = rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(ged);
      RCP<EpetraVector_ReadOnly_GlobalEvaluationData> ro_ged = rcp_dynamic_cast<EpetraVector_ReadOnly_GlobalEvaluationData>(ged);
      if(loc_pair_ged!=Teuchos::null) {
        // cast to a ThyraObjContainer throwing an exception if the cast doesn't work.
        RCP<ThyraObjContainer<Scalar> > th_ged = rcp_dynamic_cast<ThyraObjContainer<Scalar> >(loc_pair_ged->getGlobalLOC(),true);
        th_ged->set_x_th(Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(p));
      }
      else {
        TEUCHOS_ASSERT(ro_ged!=Teuchos::null);
        ro_ged->setUniqueVector(p);
      }
    }
  }

  ae_inargs.addGlobalEvaluationData(nonParamGlobalEvaluationData_);
  ae_inargs.addGlobalEvaluationData(distrParamGlobalEvaluationData_);

  // here we are building a container, this operation is fast, simply allocating a struct
  const RCP<panzer::ThyraObjContainer<Scalar> > thGlobalContainer = 
    Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<Scalar> >(ae_inargs.container_);

  TEUCHOS_ASSERT(!Teuchos::is_null(thGlobalContainer));

  // Ghosted container objects are zeroed out below only if needed for
  // a particular calculation.  This makes it more efficient than
  // zeroing out all objects in the container here.
  const RCP<panzer::ThyraObjContainer<Scalar> > thGhostedContainer = 
    Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<Scalar> >(ae_inargs.ghostedContainer_);
  
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
}

// Private functions overridden from ModelEvaulatorDefaultBase


template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
panzer::ModelEvaluator<Scalar>::createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  if(require_out_args_refresh_) {
    MEB::OutArgsSetup<Scalar> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.set_Np_Ng(parameters_.size(), g_space_.size());
    outArgs.setSupports(MEB::OUT_ARG_f);
    outArgs.setSupports(MEB::OUT_ARG_W_op);
  
    // add in dg/dx (if appropriate)
    for(std::size_t i=0;i<g_names_.size();i++) {
      typedef panzer::Traits::Jacobian RespEvalT;
  
      // check dg/dx and add it in if appropriate
      Teuchos::RCP<panzer::ResponseBase> respJacBase = responseLibrary_->getResponse<RespEvalT>(g_names_[i]);
      if(respJacBase!=Teuchos::null) {
        // cast is guranteed to succeed because of check in addResponse
        Teuchos::RCP<panzer::ResponseMESupportBase<RespEvalT> > resp 
           = Teuchos::rcp_dynamic_cast<panzer::ResponseMESupportBase<RespEvalT> >(respJacBase);
   
        // class must supppot a derivative 
        if(resp->supportsDerivative())
          outArgs.setSupports(MEB::OUT_ARG_DgDx,i,MEB::DerivativeSupport(MEB::DERIV_MV_GRADIENT_FORM));
      }
    }

    // setup parameter support
    for(std::size_t p=0;p<parameters_.size();p++) {

      if(!parameters_[p]->is_distributed)
        outArgs.setSupports(MEB::OUT_ARG_DfDp,p,MEB::DerivativeSupport(MEB::DERIV_MV_BY_COL));
      else if(parameters_[p]->is_distributed && parameters_[p]->global_indexer!=Teuchos::null)
        outArgs.setSupports(MEB::OUT_ARG_DfDp,p,MEB::DerivativeSupport(MEB::DERIV_LINEAR_OP));
    }
  
    prototypeOutArgs_ = outArgs;
  }

  // we don't need to build it anymore
  require_out_args_refresh_ = false;

  return prototypeOutArgs_;
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> > 
panzer::ModelEvaluator<Scalar>::
create_W_op() const
{
  Teuchos::RCP<const ThyraObjFactory<Scalar> > tof 
     = Teuchos::rcp_dynamic_cast<const ThyraObjFactory<Scalar> >(lof_,true);

  return tof->getThyraMatrix();
}

template <typename Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > 
panzer::ModelEvaluator<Scalar>::
get_W_factory() const
{
  return solverFactory_;
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> > 
panzer::ModelEvaluator<Scalar>::
create_DfDp_op(int p) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  typedef Thyra::ModelEvaluatorBase MEB;

  TEUCHOS_ASSERT(0<=p && p<this->Np());

  // assert that DfDp is supported
  const ParameterObject & po = *parameters_[p];

  if(po.is_distributed && po.global_indexer!=Teuchos::null) {
    TEUCHOS_ASSERT(prototypeOutArgs_.supports(MEB::OUT_ARG_DfDp,p).supports(MEB::DERIV_LINEAR_OP));

    // for a distributed parameter, figure it out from the 
    // response library
    RCP<Response_Residual<Traits::Jacobian> > response_jacobian
      = rcp_dynamic_cast<Response_Residual<Traits::Jacobian> >(po.dfdp_rl->template getResponse<Traits::Jacobian>("RESIDUAL"));

    return response_jacobian->allocateJacobian();
  }
  else if(!po.is_distributed) {
    TEUCHOS_ASSERT(prototypeOutArgs_.supports(MEB::OUT_ARG_DfDp,p).supports(MEB::DERIV_MV_BY_COL));

    // this is a scalar parameter (its easy to create!)
    return Thyra::createMember(*get_f_space());
  }

  // shourld never get here
  TEUCHOS_ASSERT(false);

  return Teuchos::null;
}

template <typename Scalar>
int panzer::ModelEvaluator<Scalar>::
addParameter(const std::string & name)
{
  Teuchos::Array<std::string> tmp_names;
  tmp_names.push_back(name);

  return addParameter(tmp_names);
}

template <typename Scalar>
int panzer::ModelEvaluator<Scalar>::
addParameter(const Teuchos::Array<std::string> & names)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::ptrFromRef;

  int parameter_index = parameters_.size();
  parameters_.push_back(createScalarParameter(names));

  require_in_args_refresh_ = true;
  require_out_args_refresh_ = true;

  return parameter_index;
}

template <typename Scalar>
int panzer::ModelEvaluator<Scalar>::
addDistributedParameter(const std::string & key,
                        const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > & vs,
                        const Teuchos::RCP<GlobalEvaluationData> & ged,
                        const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & initial,
                        const Teuchos::RCP<const UniqueGlobalIndexerBase> & ugi)
{
  distrParamGlobalEvaluationData_.addDataObject(key,ged);

  int parameter_index = parameters_.size();
  parameters_.push_back(createDistributedParameter(key,vs,initial,ugi));

  require_in_args_refresh_ = true;
  require_out_args_refresh_ = true;

  return parameter_index;
}

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
addNonParameterGlobalEvaluationData(const std::string & key,
                                    const Teuchos::RCP<GlobalEvaluationData> & ged)
{
   nonParamGlobalEvaluationData_.addDataObject(key,ged);
}


template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
applyDirichletBCs(const Teuchos::RCP<Thyra::VectorBase<Scalar> > & x,
                  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & f) const
{
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::Array;
  using Teuchos::tuple;
  using Teuchos::rcp_dynamic_cast;

  // if neccessary build a ghosted container
  if(Teuchos::is_null(ghostedContainer_)) {
     ghostedContainer_ = lof_->buildGhostedLinearObjContainer();
     lof_->initializeGhostedContainer(panzer::LinearObjContainer::X |
                                      panzer::LinearObjContainer::F,*ghostedContainer_); 
  }

  panzer::AssemblyEngineInArgs ae_inargs;
  ae_inargs.container_ = lof_->buildLinearObjContainer(); // we use a new global container
  ae_inargs.ghostedContainer_ = ghostedContainer_;        // we can reuse the ghosted container
  ae_inargs.alpha = 0.0;
  ae_inargs.beta = 1.0;
  ae_inargs.evaluate_transient_terms = false;
  ae_inargs.addGlobalEvaluationData(nonParamGlobalEvaluationData_);
  ae_inargs.addGlobalEvaluationData(distrParamGlobalEvaluationData_);

  // this is the tempory target
  lof_->initializeContainer(panzer::LinearObjContainer::F,*ae_inargs.container_); 

  // here we are building a container, this operation is fast, simply allocating a struct
  const RCP<panzer::ThyraObjContainer<Scalar> > thGlobalContainer = 
    Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<Scalar> >(ae_inargs.container_);

  TEUCHOS_ASSERT(!Teuchos::is_null(thGlobalContainer));

  // Ghosted container objects are zeroed out below only if needed for
  // a particular calculation.  This makes it more efficient than
  // zeroing out all objects in the container here.
  const RCP<panzer::ThyraObjContainer<Scalar> > thGhostedContainer = 
    Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<Scalar> >(ae_inargs.ghostedContainer_);
  Thyra::assign(thGhostedContainer->get_f_th().ptr(),0.0);
  
  // Set the solution vector (currently all targets require solution).
  // In the future we may move these into the individual cases below.
  // A very subtle (and fragile) point: A non-null pointer in global
  // container triggers export operations during fill.  Also, the
  // introduction of the container is forcing us to cast away const on
  // arguments that should be const.  Another reason to redesign
  // LinearObjContainer layers.
  thGlobalContainer->set_x_th(x);

  // evaluate dirichlet boundary conditions
  RCP<panzer::LinearObjContainer> counter = ae_tm_.template getAsObject<panzer::Traits::Residual>()->evaluateOnlyDirichletBCs(ae_inargs);

  // allocate the result container
  RCP<panzer::LinearObjContainer> result = lof_->buildLinearObjContainer(); // we use a new global container

  // stuff the evaluate boundary conditions into the f spot of the counter ... the x is already filled
  Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<Scalar> >(counter)->set_f_th(thGlobalContainer->get_f_th());
  
  // stuff the vector that needs applied dirichlet conditions in the the f spot of the result LOC
  Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<Scalar> >(result)->set_f_th(f);
  
  // use the linear object factory to apply the result
  lof_->applyDirichletBCs(*counter,*result);
}

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
              const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  evalModelImpl_basic(inArgs,outArgs); 

  // evaluate responses...uses the stored assembly arguments and containers
  if(required_basic_g(outArgs)) {
    evalModelImpl_basic_g(inArgs,outArgs);

    // evaluate response derivatives 
    if(required_basic_dgdx(outArgs))
      evalModelImpl_basic_dgdx(inArgs,outArgs);
  }

  if(required_basic_dfdp_scalar(outArgs))
    evalModelImpl_basic_dfdp_scalar(inArgs,outArgs);

  if(required_basic_dfdp_distro(outArgs))
    evalModelImpl_basic_dfdp_distro(inArgs,outArgs);
}

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
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

  // see if the user wants us to do anything
  if(Teuchos::is_null(f_out) && Teuchos::is_null(W_out) ) {
     return;
  }

  // setup all the assembly in arguments (this is parameters and
  // x/x_dot). At this point with the exception of the one time dirichlet
  // beta that is all thats neccessary.
  panzer::AssemblyEngineInArgs ae_inargs;
  setupAssemblyInArgs(inArgs,ae_inargs);

  // handle application of the one time dirichlet beta in the
  // assembly engine. Note that this has to be set explicitly
  // each time because this badly breaks encapsulation. Essentially
  // we must work around the model evaluator abstraction!
  if(oneTimeDirichletBeta_on_) {
    ae_inargs.dirichlet_beta = oneTimeDirichletBeta_;
    ae_inargs.apply_dirichlet_beta = true;

    oneTimeDirichletBeta_on_ = false;
  }

  // here we are building a container, this operation is fast, simply allocating a struct
  const RCP<panzer::ThyraObjContainer<Scalar> > thGlobalContainer = 
    Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<Scalar> >(ae_inargs.container_);
  const RCP<panzer::ThyraObjContainer<Scalar> > thGhostedContainer = 
    Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<Scalar> >(ae_inargs.ghostedContainer_);

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

  // HACK: set A to null before calling responses to avoid touching the
  // the Jacobian after it has been properly assembled.  Should be fixed
  // by using a modified version of ae_inargs instead.
  thGlobalContainer->set_A_th(Teuchos::null);

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

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
evalModelImpl_basic_g(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
                      const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  // optional sanity check
  // TEUCHOS_ASSERT(required_basic_g(outArgs));

  // setup all the assembly in arguments (this is parameters and
  // x/x_dot). At this point with the exception of the one time dirichlet
  // beta that is all thats neccessary.
  panzer::AssemblyEngineInArgs ae_inargs;
  setupAssemblyInArgs(inArgs,ae_inargs);

  for(std::size_t i=0;i<g_names_.size();i++) {
    Teuchos::RCP<Thyra::VectorBase<Scalar> > vec = outArgs.get_g(i);
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

template <typename Scalar>
void 
panzer::ModelEvaluator<Scalar>::
evalModelImpl_basic_dgdx(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
                         const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  // optional sanity check
  TEUCHOS_ASSERT(required_basic_dgdx(outArgs));

  for(std::size_t i=0;i<g_names_.size();i++) {
    // get "Vector" out of derivative, if its something else, throw an exception
    MEB::Derivative<Scalar> deriv = outArgs.get_DgDx(i);
    if(deriv.isEmpty())
      continue;

    Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > vec = deriv.getMultiVector();

    if(vec!=Teuchos::null) {

      std::string responseName = g_names_[i];
      Teuchos::RCP<panzer::ResponseMESupportBase<panzer::Traits::Jacobian> > resp 
          = Teuchos::rcp_dynamic_cast<panzer::ResponseMESupportBase<panzer::Traits::Jacobian> >(responseLibrary_->getResponse<panzer::Traits::Jacobian>(responseName));
      resp->setDerivative(vec);
    }
  }

  // setup all the assembly in arguments (this is parameters and
  // x/x_dot). At this point with the exception of the one time dirichlet
  // beta that is all thats neccessary.
  panzer::AssemblyEngineInArgs ae_inargs;
  setupAssemblyInArgs(inArgs,ae_inargs);

  // evaluator responses
  responseLibrary_->addResponsesToInArgs<panzer::Traits::Jacobian>(ae_inargs);
  responseLibrary_->evaluate<panzer::Traits::Jacobian>(ae_inargs);
}

template <typename Scalar>
void 
panzer::ModelEvaluator<Scalar>::
evalModelImpl_basic_dfdp_scalar(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
                                const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   typedef Thyra::ModelEvaluatorBase MEB;

   TEUCHOS_ASSERT(required_basic_dfdp_scalar(outArgs));

   // setup all the assembly in arguments (this is parameters and
   // x/x_dot). At this point with the exception of the one time dirichlet
   // beta that is all thats neccessary.
   panzer::AssemblyEngineInArgs ae_inargs;
   setupAssemblyInArgs(inArgs,ae_inargs);

   // First: Fill the output vectors from the input arguments structure. Put them
   //        in the global evaluation data container so they are correctly communicated.
   ///////////////////////////////////////////////////////////////////////////////////////

   std::vector<std::string> activeParameters;

   int totalParameterCount = 0;
   for(std::size_t i=0; i < parameters_.size(); i++) {
     // skip non-scalar parameters
     if(parameters_[i]->is_distributed)
       continue;

     // have derivatives been requested?
     MEB::Derivative<Scalar> deriv = outArgs.get_DfDp(i);
     if(deriv.isEmpty())
       continue;

     // grab multivector, make sure its the right dimension
     Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > mVec = deriv.getMultiVector();
     TEUCHOS_ASSERT(mVec->domain()->dim()==Teuchos::as<int>(parameters_[i]->scalar_value.size()));

     for (std::size_t j=0; j < parameters_[i]->scalar_value.size(); j++) {

       // build containers for each vector
       RCP<LOCPair_GlobalEvaluationData> loc_pair = Teuchos::rcp(new LOCPair_GlobalEvaluationData(lof_,LinearObjContainer::F));
       RCP<LinearObjContainer> globalContainer = loc_pair->getGlobalLOC();

       // stuff target vector into global container
       RCP<Thyra::VectorBase<Scalar> > vec = mVec->col(j);
       RCP<panzer::ThyraObjContainer<Scalar> > thGlobalContainer = 
         Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<Scalar> >(globalContainer);
       thGlobalContainer->set_f_th(vec);

       // add container into in args object
       std::string name = "PARAMETER_SENSITIVIES: "+(*parameters_[i]->names)[j];
       ae_inargs.addGlobalEvaluationData(name,loc_pair->getGhostedLOC());
       ae_inargs.addGlobalEvaluationData(name+"_pair",loc_pair);

       activeParameters.push_back(name);
       totalParameterCount++;
     }
   }

   // Second: For all parameters that require derivative sensitivities, put in a name 
   //         so that the scatter can realize which sensitivity vectors it needs to fill
   ///////////////////////////////////////////////////////////////////////////////////////

   RCP<GlobalEvaluationData> ged_activeParameters = Teuchos::rcp(new ParameterList_GlobalEvaluationData(activeParameters));
   ae_inargs.addGlobalEvaluationData("PARAMETER_NAMES",ged_activeParameters);

   // Third: Now seed all the parameters in the parameter vector so that derivatives
   //        can be properly computed.
   ///////////////////////////////////////////////////////////////////////////////////////

   int paramIndex = 0;
   for(std::size_t i=0; i < parameters_.size(); i++) {
     // skip non-scalar parameters
     if(parameters_[i]->is_distributed)
       continue;

     // don't modify the parameter if its not needed
     MEB::Derivative<Scalar> deriv = outArgs.get_DfDp(i);
     if(deriv.isEmpty()) {
       // reinitialize values that should not have sensitivities computed (this is a precaution)
       for (unsigned int j=0; j < parameters_[i]->scalar_value.size(); j++) {
         Traits::FadType p = Traits::FadType(totalParameterCount, parameters_[i]->scalar_value[j].baseValue);
         parameters_[i]->scalar_value[j].family->template setValue<panzer::Traits::Tangent>(p);
       }
       continue;
     }
     else {
       // loop over each parameter in the vector, initializing the AD type
       for (unsigned int j=0; j < parameters_[i]->scalar_value.size(); j++) {
         Traits::FadType p = Traits::FadType(totalParameterCount, parameters_[i]->scalar_value[j].baseValue);
         p.fastAccessDx(paramIndex) = 1.0;
         parameters_[i]->scalar_value[j].family->template setValue<panzer::Traits::Tangent>(p);
         paramIndex++;
       }
     }
   }

   // make sure that the total parameter count and the total parameter index match up
   TEUCHOS_ASSERT(paramIndex==totalParameterCount);

   // Fourth: Actually evaluate the residual's sensitivity to the parameters
   ///////////////////////////////////////////////////////////////////////////////////////

   if(totalParameterCount>0) {
     ae_tm_.getAsObject<panzer::Traits::Tangent>()->evaluate(ae_inargs);
   }
}

template <typename Scalar>
void 
panzer::ModelEvaluator<Scalar>::
evalModelImpl_basic_dfdp_distro(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
                                const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::null;

   typedef Thyra::ModelEvaluatorBase MEB;

  TEUCHOS_ASSERT(required_basic_dfdp_distro(outArgs));

  // loop over parameters, and then build a dfdp_rl only if they are distributed
  // and the user has provided the UGI. Note that this may be overly expensive if they
  // don't actually want those sensitivites because memory will be allocated unneccesarily.
  // It would be good to do this "just in time", but for now this is sufficient.
  for(std::size_t p=0;p<parameters_.size();p++) {

    // parameter is not distributed, a different path is
    // taken for those to compute dfdp
    if(!parameters_[p]->is_distributed)
      continue;

    // parameter is distributed but has no global indexer.
    // thus the user doesn't want sensitivities!
    if(parameters_[p]->dfdp_rl==null)
      continue;

    // have derivatives been requested?
    MEB::Derivative<Scalar> deriv = outArgs.get_DfDp(p);
    if(deriv.isEmpty())
      continue;
 
    ResponseLibrary<Traits> & rLibrary = *parameters_[p]->dfdp_rl;

    // get the response and tell it to fill the derivative operator
    RCP<Response_Residual<Traits::Jacobian> > response_jacobian =
      rcp_dynamic_cast<Response_Residual<Traits::Jacobian> >(rLibrary.getResponse<Traits::Jacobian>("RESIDUAL"));
    response_jacobian->setJacobian(deriv.getLinearOp());

    // setup all the assembly in arguments (this is parameters and x/x_dot). 
    // make sure the correct seeding is performed
    panzer::AssemblyEngineInArgs ae_inargs;
    setupAssemblyInArgs(inArgs,ae_inargs);

    ae_inargs.sensitivities_name = (*parameters_[p]->names)[0]; // distributed parameters can only have one name!
    ae_inargs.gather_seeds.push_back(1.0); // this assumes that gather point is always the zero index of
                                           // gather seeds
    rLibrary.addResponsesToInArgs<Traits::Jacobian>(ae_inargs);

    rLibrary.evaluate<Traits::Jacobian>(ae_inargs);
  }
}

template <typename Scalar>
bool panzer::ModelEvaluator<Scalar>::
required_basic_g(const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
   // determine if any of the outArgs are not null!
   bool activeGArgs = false;
   for(int i=0;i<outArgs.Ng();i++) 
      activeGArgs |= (outArgs.get_g(i)!=Teuchos::null); 

   return activeGArgs | required_basic_dgdx(outArgs);
}

template <typename Scalar>
bool panzer::ModelEvaluator<Scalar>::
required_basic_dgdx(const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
   typedef Thyra::ModelEvaluatorBase MEB;

   // determine if any of the outArgs are not null!
   bool activeGArgs = false;
   for(int i=0;i<outArgs.Ng();i++) {
     // no derivatives are supported
     if(outArgs.supports(MEB::OUT_ARG_DgDx,i).none())
       continue;

     // this is basically a redundant computation
     activeGArgs |= (!outArgs.get_DgDx(i).isEmpty());
   }

   return activeGArgs;
}

template <typename Scalar>
bool panzer::ModelEvaluator<Scalar>::
required_basic_dfdp_scalar(const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
   typedef Thyra::ModelEvaluatorBase MEB;

   // determine if any of the outArgs are not null!
   bool activeFPArgs = false;
   for(int i=0;i<outArgs.Np();i++) {
 
     // this is for scalar parameters only
     if(parameters_[i]->is_distributed)
       continue;

     // no derivatives are supported
     if(outArgs.supports(MEB::OUT_ARG_DfDp,i).none())
       continue;

     // this is basically a redundant computation
     activeFPArgs |= (!outArgs.get_DfDp(i).isEmpty());
   }

   return activeFPArgs;
}

template <typename Scalar>
bool panzer::ModelEvaluator<Scalar>::
required_basic_dfdp_distro(const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
   typedef Thyra::ModelEvaluatorBase MEB;

   // determine if any of the outArgs are not null!
   bool activeFPArgs = false;
   for(int i=0;i<outArgs.Np();i++) {
 
     // this is for scalar parameters only
     if(!parameters_[i]->is_distributed)
       continue;

     // no derivatives are supported
     if(outArgs.supports(MEB::OUT_ARG_DfDp,i).none())
       continue;

     // this is basically a redundant computation
     activeFPArgs |= (!outArgs.get_DfDp(i).isEmpty());
   }

   return activeFPArgs;
}

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
buildDistroParamDfDp_RL(
       const Teuchos::RCP<panzer::WorksetContainer> & wc,
       const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
       const std::vector<panzer::BC> & bcs,
       const panzer::EquationSetFactory & eqset_factory,
       const panzer::BCStrategyFactory& bc_factory,
       const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
       const Teuchos::ParameterList& closure_models,
       const Teuchos::ParameterList& user_data,
       const bool write_graphviz_file,
       const std::string& graphviz_file_prefix)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::null;

  // loop over parameters, and then build a dfdp_rl only if they are distributed
  // and the user has provided the UGI. Note that this may be overly expensive if they
  // don't actually want those sensitivites because memory will be allocated unneccesarily.
  // It would be good to do this "just in time", but for now this is sufficient.
  for(std::size_t p=0;p<parameters_.size();p++) {
    // parameter is not distributed, a different path is
    // taken for those to compute dfdp
    if(!parameters_[p]->is_distributed)
      continue;

    // parameter is distributed but has no global indexer.
    // thus the user doesn't want sensitivities!
    if(parameters_[p]->global_indexer==null)
      continue;

    // build the linear object factory that has the correct sizing for
    // the sensitivity matrix (parameter sized domain, residual sized range)
    RCP<const LinearObjFactory<Traits> > param_lof = cloneWithNewDomain(*lof_,parameters_[p]->global_indexer);

    // the user wants global sensitivities, hooray! Build and setup the response library 
    RCP<ResponseLibrary<Traits> > rLibrary
        = Teuchos::rcp(new ResponseLibrary<Traits>(wc,lof_->getUniqueGlobalIndexerBase(),param_lof,true));
    rLibrary->buildResidualResponseEvaluators(physicsBlocks,eqset_factory,bcs,bc_factory,
                                              cm_factory,closure_models,user_data,
                                              write_graphviz_file,graphviz_file_prefix);

    // make sure parameter response library is correct
    parameters_[p]->dfdp_rl = rLibrary;
  }
}

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
setOneTimeDirichletBeta(const Scalar & beta) const
{
  oneTimeDirichletBeta_on_ = true;
  oneTimeDirichletBeta_    = beta;
}

template <typename Scalar>
Teuchos::RCP<typename panzer::ModelEvaluator<Scalar>::ParameterObject> 
panzer::ModelEvaluator<Scalar>::
createScalarParameter(const Teuchos::Array<std::string> & in_names) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  RCP<ParameterObject> paramObj = rcp(new ParameterObject);

  paramObj->names = rcp(new Teuchos::Array<std::string>(in_names));
  paramObj->is_distributed = false;

  // associate vector with the ParamLib
  paramObj->scalar_value = panzer::ParamVec();

  global_data_->pl->fillVector<panzer::Traits::Residual>(*paramObj->names, paramObj->scalar_value);

  // build initial condition vector
  paramObj->space =
    Thyra::locallyReplicatedDefaultSpmdVectorSpace<Scalar>(rcp(new Teuchos::MpiComm<long int>(lof_->getComm().getRawMpiComm())),
                                                                                              paramObj->names->size());
  
  // fill vector with parameter values
  Teuchos::ArrayRCP<Scalar> data;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > initial_value 
      = Thyra::createMember(paramObj->space);
  RCP<Thyra::SpmdVectorBase<Scalar> > vec
      = rcp_dynamic_cast<Thyra::SpmdVectorBase<Scalar> >(initial_value);
  vec->getNonconstLocalData(ptrFromRef(data));
  for (unsigned int i=0; i < paramObj->scalar_value.size(); i++)
    data[i] = paramObj->scalar_value[i].baseValue;

  paramObj->initial_value = initial_value;

  return paramObj;
}

template <typename Scalar>
Teuchos::RCP<typename panzer::ModelEvaluator<Scalar>::ParameterObject> 
panzer::ModelEvaluator<Scalar>::
createDistributedParameter(const std::string & key,
                           const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > & vs,
                           const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & initial,
                           const Teuchos::RCP<const UniqueGlobalIndexerBase> & ugi) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  RCP<ParameterObject> paramObj = rcp(new ParameterObject);

  paramObj->is_distributed = true;
  paramObj->names = rcp(new Teuchos::Array<std::string>());
  paramObj->names->push_back(key);
  paramObj->space = vs;
  paramObj->initial_value = initial;

  paramObj->global_indexer = ugi;

  return paramObj;
}

#endif
