// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_ModelEvaluator_impl_hpp__
#define   __Panzer_ModelEvaluator_impl_hpp__

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_ArrayRCP.hpp"

#include "PanzerDiscFE_config.hpp"
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
#include "Panzer_ParameterLibraryUtilities.hpp"
#include "Panzer_LinearObjFactory_Utilities.hpp"

#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_DefaultSpmdVector.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"

// For writing out residuals/Jacobians
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_TpetraVector.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Tpetra_CrsMatrix.hpp"

// Constructors/Initializers/Accessors

template<typename Scalar>
panzer::ModelEvaluator<Scalar>::
ModelEvaluator(const Teuchos::RCP<panzer::FieldManagerBuilder>& fmb,
               const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> >& rLibrary,
               const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> >& lof,
               const std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >& p_names,
               const std::vector<Teuchos::RCP<Teuchos::Array<double> > >& p_values,
               const Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > & solverFactory,
               const Teuchos::RCP<panzer::GlobalData>& global_data,
               bool build_transient_support,
               double t_init)
  : t_init_(t_init)
  , num_me_parameters_(0)
  , do_fd_dfdp_(false)
  , fd_perturb_size_(1e-7)
  , require_in_args_refresh_(true)
  , require_out_args_refresh_(true)
  , responseLibrary_(rLibrary)
  , global_data_(global_data)
  , build_transient_support_(build_transient_support)
  , lof_(lof)
  , solverFactory_(solverFactory)
  , oneTimeDirichletBeta_on_(false)
  , oneTimeDirichletBeta_(0.0)
  , build_volume_field_managers_(true)
  , build_bc_field_managers_(true)
  , active_evaluation_types_(Sacado::mpl::size<panzer::Traits::EvalTypes>::value, true)
  , write_matrix_count_(0)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::tuple;
  using Thyra::VectorBase;
  using Thyra::createMember;

  TEUCHOS_ASSERT(lof_!=Teuchos::null);

  panzer::AssemblyEngine_TemplateBuilder builder(fmb,lof);
  ae_tm_.buildObjects(builder);

  //
  // Build x, f spaces
  //

  // dynamic cast to blocked LOF for now
  RCP<const ThyraObjFactory<Scalar> > tof = rcp_dynamic_cast<const ThyraObjFactory<Scalar> >(lof,true);

  x_space_ = tof->getThyraDomainSpace();
  f_space_ = tof->getThyraRangeSpace();

  //
  // Setup parameters
  //
  for(std::size_t i=0;i<p_names.size();i++)
     addParameter(*(p_names[i]),*(p_values[i]));

  // now that the vector spaces are setup we can allocate the nominal values
  // (i.e. initial conditions)
  initializeNominalValues();
}

template<typename Scalar>
panzer::ModelEvaluator<Scalar>::
ModelEvaluator(const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> >& lof,
               const Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > & solverFactory,
               const Teuchos::RCP<panzer::GlobalData>& global_data,
               bool build_transient_support,double t_init)
  : t_init_(t_init)
  , num_me_parameters_(0)
  , do_fd_dfdp_(false)
  , fd_perturb_size_(1e-7)
  , require_in_args_refresh_(true)
  , require_out_args_refresh_(true)
  , global_data_(global_data)
  , build_transient_support_(build_transient_support)
  , lof_(lof)
  , solverFactory_(solverFactory)
  , oneTimeDirichletBeta_on_(false)
  , oneTimeDirichletBeta_(0.0)
  , build_volume_field_managers_(true)
  , build_bc_field_managers_(true)
  , active_evaluation_types_(Sacado::mpl::size<panzer::Traits::EvalTypes>::value, true)
  , write_matrix_count_(0)
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
  TEUCHOS_TEST_FOR_EXCEPTION(!(i>=0 && i<num_me_parameters_),std::runtime_error,
                             "panzer::ModelEvaluator::get_p_names: Requested parameter index out of range.");

  if (i < Teuchos::as<int>(parameters_.size()))
    return parameters_[i]->names;
  else if (i < Teuchos::as<int>(parameters_.size()+tangent_space_.size())) {
    Teuchos::RCP<Teuchos::Array<std::string> > names = rcp(new Teuchos::Array<std::string>);
    int param_index = i-parameters_.size();
    std::ostringstream ss;
    ss << "TANGENT VECTOR: " << param_index;
    names->push_back(ss.str());
    return names;
  }
  else if (build_transient_support_ && i < Teuchos::as<int>(parameters_.size()+2*tangent_space_.size())) {
    Teuchos::RCP<Teuchos::Array<std::string> > names = rcp(new Teuchos::Array<std::string>);
    int param_index = i-parameters_.size()-tangent_space_.size();
    std::ostringstream ss;
    ss << "TIME DERIVATIVE TANGENT VECTOR: " << param_index;
    names->push_back(ss.str());
    return names;
  }

  return Teuchos::null;
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
panzer::ModelEvaluator<Scalar>::get_p_space(int i) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!(i>=0 && i<num_me_parameters_),std::runtime_error,
                             "panzer::ModelEvaluator::get_p_space: Requested parameter index out of range.");

  if (i < Teuchos::as<int>(parameters_.size()))
    return parameters_[i]->space;
  else if (i < Teuchos::as<int>(parameters_.size()+tangent_space_.size()))
    return tangent_space_[i-parameters_.size()];
  else if (build_transient_support_ && i < Teuchos::as<int>(parameters_.size()+2*tangent_space_.size()))
    return tangent_space_[i-parameters_.size()-tangent_space_.size()];

  return Teuchos::null;
}

template<typename Scalar>
Teuchos::ArrayView<const std::string>
panzer::ModelEvaluator<Scalar>::get_g_names(int i) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!(i>=0 && i<Teuchos::as<int>(responses_.size())),std::runtime_error,
                             "panzer::ModelEvaluator::get_g_names: Requested response index out of range.");

  return Teuchos::ArrayView<const std::string>(&(responses_[i]->name),1);
}

template<typename Scalar>
const std::string &
panzer::ModelEvaluator<Scalar>::get_g_name(int i) const
{
  TEUCHOS_ASSERT(i>=0 &&
                 static_cast<typename std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > >::size_type>(i)<responses_.size());

  return responses_[i]->name;
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
panzer::ModelEvaluator<Scalar>::get_g_space(int i) const
{
  TEUCHOS_ASSERT(i>=0 &&
                 static_cast<typename std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > >::size_type>(i)<responses_.size());

  return responses_[i]->space;
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
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

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
    nomInArgs.set_Np(num_me_parameters_);
    for(std::size_t p=0;p<parameters_.size();p++) {
      // setup nominal in arguments
      nomInArgs.set_p(p,parameters_[p]->initial_value);

      // We explicitly do not set nominal values for tangent parameters
      // as these are parameters that should be hidden from client code
    }

    nominalValues_ = nomInArgs;
  }

  // refresh no longer required
  require_in_args_refresh_ = false;

  return nominalValues_;
}

template<typename Scalar>
void
panzer::ModelEvaluator<Scalar>::initializeNominalValues() const
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
    nomInArgs.setSupports(MEB::IN_ARG_step_size,true);
    nomInArgs.setSupports(MEB::IN_ARG_stage_number,true);

    Teuchos::RCP<Thyra::VectorBase<Scalar> > x_dot_nom = Thyra::createMember(x_space_);
    Thyra::assign(x_dot_nom.ptr(),0.0);
    nomInArgs.set_x_dot(x_dot_nom);
    nomInArgs.set_t(t_init_);
    nomInArgs.set_alpha(0.0); // these have no meaning initially!
    nomInArgs.set_beta(0.0);
    //TODO: is this needed?
    nomInArgs.set_step_size(0.0);
    nomInArgs.set_stage_number(1.0);
  }

  // setup parameter support -- for each scalar parameter we support the parameter itself and tangent vectors for x, xdot
  nomInArgs.set_Np(num_me_parameters_);
  std::size_t v_index = 0;
  for(std::size_t p=0;p<parameters_.size();p++) {
    nomInArgs.set_p(p,parameters_[p]->initial_value);
    if (!parameters_[p]->is_distributed) {
      Teuchos::RCP<Thyra::VectorBase<Scalar> > v_nom_x    = Thyra::createMember(*tangent_space_[v_index]);
      Thyra::assign(v_nom_x.ptr(),0.0);
      nomInArgs.set_p(v_index+parameters_.size(),v_nom_x);
      if (build_transient_support_) {
        Teuchos::RCP<Thyra::VectorBase<Scalar> > v_nom_xdot = Thyra::createMember(*tangent_space_[v_index]);
        Thyra::assign(v_nom_xdot.ptr(),0.0);
        nomInArgs.set_p(v_index+parameters_.size()+tangent_space_.size(),v_nom_xdot);
      }
      ++v_index;
    }
  }

  nominalValues_ = nomInArgs;
}

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
buildVolumeFieldManagers(const bool value)
{
  build_volume_field_managers_ = value;
}

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
buildBCFieldManagers(const bool value)
{
  build_bc_field_managers_ = value;
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
           bool writeGraph,const std::string & graphPrefix,
           const Teuchos::ParameterList& me_params)
{
  // First: build residual assembly engine
  /////////////////////////////////////////////////////////////////////////////////////////////////
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::ModelEvaluator::setupModel()",setupModel);

  {
    // 1. build Field manager builder
    /////////////////////////////////////////////////////////////////////////////////////////////////

    Teuchos::RCP<panzer::FieldManagerBuilder> fmb;
    {
      PANZER_FUNC_TIME_MONITOR_DIFF("allocate FieldManagerBuilder",allocFMB);
      fmb = Teuchos::rcp(new panzer::FieldManagerBuilder);
      fmb->setActiveEvaluationTypes(active_evaluation_types_);
    }
    {
      PANZER_FUNC_TIME_MONITOR_DIFF("fmb->setWorksetContainer()",setupWorksets);
      fmb->setWorksetContainer(wc);
    }
    if (build_volume_field_managers_) {
      PANZER_FUNC_TIME_MONITOR_DIFF("fmb->setupVolumeFieldManagers()",setupVolumeFieldManagers);
      fmb->setupVolumeFieldManagers(physicsBlocks,volume_cm_factory,closure_models,*lof_,user_data);
    }
    if (build_bc_field_managers_) {
      PANZER_FUNC_TIME_MONITOR_DIFF("fmb->setupBCFieldManagers()",setupBCFieldManagers);
      fmb->setupBCFieldManagers(bcs,physicsBlocks,eqset_factory,bc_cm_factory,bc_factory,closure_models,*lof_,user_data);
    }

    // Print Phalanx DAGs
    if (writeGraph){
      if (build_volume_field_managers_)
        fmb->writeVolumeGraphvizDependencyFiles(graphPrefix, physicsBlocks);
      if (build_bc_field_managers_)
        fmb->writeBCGraphvizDependencyFiles(graphPrefix+"BC_");
    }

    {
      PANZER_FUNC_TIME_MONITOR_DIFF("AssemblyEngine_TemplateBuilder::buildObjects()",AETM_BuildObjects);
      panzer::AssemblyEngine_TemplateBuilder builder(fmb,lof_);
      ae_tm_.buildObjects(builder);
    }
  }

  // Second: build the responses
  /////////////////////////////////////////////////////////////////////////////////////////////////

  {
    PANZER_FUNC_TIME_MONITOR_DIFF("build response library",buildResponses);

    responseLibrary_->initialize(wc,lof_->getRangeGlobalIndexer(),lof_);

    buildResponses(physicsBlocks,eqset_factory,volume_cm_factory,closure_models,user_data,writeGraph,graphPrefix+"Responses_");
    buildDistroParamDfDp_RL(wc,physicsBlocks,bcs,eqset_factory,bc_factory,volume_cm_factory,closure_models,user_data,writeGraph,graphPrefix+"Response_DfDp_");
    buildDistroParamDgDp_RL(wc,physicsBlocks,bcs,eqset_factory,bc_factory,volume_cm_factory,closure_models,user_data,writeGraph,graphPrefix+"Response_DgDp_");

    do_fd_dfdp_ = false;
    fd_perturb_size_ = 1.0e-7;
    if (me_params.isParameter("FD Forward Sensitivities"))
      do_fd_dfdp_ = me_params.get<bool>("FD Forward Sensitivities");
    if (me_params.isParameter("FD Perturbation Size"))
      fd_perturb_size_ = me_params.get<double>("FD Perturbation Size");
  }
}

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
setupAssemblyInArgs(const Thyra::ModelEvaluatorBase::InArgs<Scalar> & inArgs,
                    panzer::AssemblyEngineInArgs & ae_inargs) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_const_cast;
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

  if(Teuchos::is_null(xContainer_))
    xContainer_    = lof_->buildReadOnlyDomainContainer();
  if(Teuchos::is_null(xdotContainer_) && is_transient)
    xdotContainer_ = lof_->buildReadOnlyDomainContainer();

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

    ae_inargs.step_size= inArgs.get_step_size();
    ae_inargs.stage_number = inArgs.get_stage_number();
    ae_inargs.evaluate_transient_terms = true;
  }

  // this member is handled in the individual functions
  ae_inargs.apply_dirichlet_beta = false;

  // Set input parameters
  int num_param_vecs = parameters_.size();
  for (int i=0; i<num_param_vecs; i++) {

    RCP<const Thyra::VectorBase<Scalar> > paramVec = inArgs.get_p(i);
    if ( paramVec!=Teuchos::null && !parameters_[i]->is_distributed) {
      // non distributed parameters

      Teuchos::ArrayRCP<const Scalar> p_data;
      rcp_dynamic_cast<const Thyra::SpmdVectorBase<Scalar> >(paramVec,true)->getLocalData(Teuchos::ptrFromRef(p_data));

      for (unsigned int j=0; j < parameters_[i]->scalar_value.size(); j++) {
        parameters_[i]->scalar_value[j].baseValue = p_data[j];
        parameters_[i]->scalar_value[j].family->setRealValueForAllTypes(parameters_[i]->scalar_value[j].baseValue);
      }
    }
    else if ( paramVec!=Teuchos::null && parameters_[i]->is_distributed) {
      // distributed parameters

      std::string key = (*parameters_[i]->names)[0];
      RCP<GlobalEvaluationData> ged = distrParamGlobalEvaluationData_.getDataObject(key);

      TEUCHOS_ASSERT(ged!=Teuchos::null);

      // cast to a LOCPair throwing an exception if the cast doesn't work.
      RCP<LOCPair_GlobalEvaluationData> loc_pair_ged = rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(ged);
      RCP<ReadOnlyVector_GlobalEvaluationData> ro_ged = rcp_dynamic_cast<ReadOnlyVector_GlobalEvaluationData>(ged);
      if(loc_pair_ged!=Teuchos::null) {
        // cast to a ThyraObjContainer throwing an exception if the cast doesn't work.
        RCP<ThyraObjContainer<Scalar> > th_ged = rcp_dynamic_cast<ThyraObjContainer<Scalar> >(loc_pair_ged->getGlobalLOC(),true);
        th_ged->set_x_th(Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(paramVec));
      }
      else {
        TEUCHOS_ASSERT(ro_ged!=Teuchos::null);
        ro_ged->setOwnedVector(paramVec);
      }
    }
  }

  ae_inargs.addGlobalEvaluationData(distrParamGlobalEvaluationData_);

  // here we are building a container, this operation is fast, simply allocating a struct
  const RCP<panzer::ThyraObjContainer<Scalar> > thGlobalContainer =
    Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<Scalar> >(ae_inargs.container_);

  TEUCHOS_ASSERT(!Teuchos::is_null(thGlobalContainer));

  // Ghosted container objects are zeroed out below only if needed for
  // a particular calculation.  This makes it more efficient than
  // zeroing out all objects in the container here.
  // const RCP<panzer::ThyraObjContainer<Scalar> > thGhostedContainer =
  //   Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<Scalar> >(ae_inargs.ghostedContainer_);

  // Set the solution vector (currently all targets require solution).
  // In the future we may move these into the individual cases below.
  // A very subtle (and fragile) point: A non-null pointer in global
  // container triggers export operations during fill.  Also, the
  // introduction of the container is forcing us to cast away const on
  // arguments that should be const.  Another reason to redesign
  // LinearObjContainer layers.
  thGlobalContainer->set_x_th(Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(x));
  xContainer_->setOwnedVector(x);
  ae_inargs.addGlobalEvaluationData("Solution Gather Container - X",xContainer_);

  if (is_transient) {
    thGlobalContainer->set_dxdt_th(Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(x_dot));
    xdotContainer_->setOwnedVector(x_dot);
    ae_inargs.addGlobalEvaluationData("Solution Gather Container - Xdot",xdotContainer_);
  }

  // Add tangent vectors for x and xdot to GlobalEvaluationData, one for each
  // scalar parameter vector and parameter within that vector.
  // Note:  The keys for the global evaluation data containers for the tangent
  //        vectors are constructed in EquationSet_AddFieldDefaultImpl::
  //        buildAndRegisterGatherAndOrientationEvaluators().
  int vIndex(0);
  for (int i(0); i < num_param_vecs; ++i)
  {
    using std::string;
    using Thyra::ProductVectorBase;
    using Thyra::VectorBase;
    using ROVGED = panzer::ReadOnlyVector_GlobalEvaluationData;
    if (not parameters_[i]->is_distributed)
    {
      auto dxdp = rcp_const_cast<VectorBase<Scalar>>
        (inArgs.get_p(vIndex + num_param_vecs));
      if (not dxdp.is_null())
      {
        // We need to cast away const because the object container requires
        // non-const vectors.
        auto dxdpBlock = rcp_dynamic_cast<ProductVectorBase<Scalar>>(dxdp);
        int numParams(parameters_[i]->scalar_value.size());
        for (int j(0); j < numParams; ++j)
        {
          RCP<ROVGED> dxdpContainer = lof_->buildReadOnlyDomainContainer();
          dxdpContainer->setOwnedVector(dxdpBlock->getNonconstVectorBlock(j));
          string name("X TANGENT GATHER CONTAINER: " +
            (*parameters_[i]->names)[j]);
          ae_inargs.addGlobalEvaluationData(name, dxdpContainer);
        } // end loop over the parameters
      } // end if (not dxdp.is_null())
      if (build_transient_support_)
      {
        // We need to cast away const because the object container requires
        // non-const vectors.
        auto dxdotdp = rcp_const_cast<VectorBase<Scalar>>
          (inArgs.get_p(vIndex + num_param_vecs + tangent_space_.size()));
        if (not dxdotdp.is_null())
        {
          auto dxdotdpBlock =
            rcp_dynamic_cast<ProductVectorBase<Scalar>>(dxdotdp);
          int numParams(parameters_[i]->scalar_value.size());
          for (int j(0); j < numParams; ++j)
          {
            RCP<ROVGED> dxdotdpContainer = lof_->buildReadOnlyDomainContainer();
            dxdotdpContainer->setOwnedVector(
              dxdotdpBlock->getNonconstVectorBlock(j));
            string name("DXDT TANGENT GATHER CONTAINER: " +
              (*parameters_[i]->names)[j]);
            ae_inargs.addGlobalEvaluationData(name, dxdotdpContainer);
          } // end loop over the parameters
        } // end if (not dxdotdp.is_null())
      } // end if (build_transient_support_)
      ++vIndex;
    } // end if (not parameters_[i]->is_distributed)
  } // end loop over the parameter vectors
} // end of setupAssemblyInArgs()

// Private functions overridden from ModelEvaulatorDefaultBase


template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
panzer::ModelEvaluator<Scalar>::createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  if(require_out_args_refresh_) {
    MEB::OutArgsSetup<Scalar> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.set_Np_Ng(num_me_parameters_, responses_.size());
    outArgs.setSupports(MEB::OUT_ARG_f);
    outArgs.setSupports(MEB::OUT_ARG_W_op);

    // add in dg/dx (if appropriate)
    for(std::size_t i=0;i<responses_.size();i++) {
      typedef panzer::Traits::Jacobian RespEvalT;

      // check dg/dx and add it in if appropriate
      Teuchos::RCP<panzer::ResponseBase> respJacBase
          = responseLibrary_->getResponse<RespEvalT>(responses_[i]->name);
      if(respJacBase!=Teuchos::null) {
        // cast is guranteed to succeed because of check in addResponse
        Teuchos::RCP<panzer::ResponseMESupportBase<RespEvalT> > resp
           = Teuchos::rcp_dynamic_cast<panzer::ResponseMESupportBase<RespEvalT> >(respJacBase);

        // class must supppot a derivative
        if(resp->supportsDerivative()) {
          outArgs.setSupports(MEB::OUT_ARG_DgDx,i,MEB::DerivativeSupport(MEB::DERIV_MV_GRADIENT_FORM));


          for(std::size_t p=0;p<parameters_.size();p++) {
            if(parameters_[p]->is_distributed && parameters_[p]->global_indexer!=Teuchos::null)
              outArgs.setSupports(MEB::OUT_ARG_DgDp,i,p,MEB::DerivativeSupport(MEB::DERIV_MV_GRADIENT_FORM));
            if(!parameters_[p]->is_distributed)
              outArgs.setSupports(MEB::OUT_ARG_DgDp,i,p,MEB::DerivativeSupport(MEB::DERIV_MV_JACOBIAN_FORM));
          }
        }
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
  PANZER_FUNC_TIME_MONITOR("panzer::ModelEvaluator::create_W_op");
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

  // The code below uses prototypeOutArgs_, so we need to make sure it is
  // initialized before using it.  This happens through createOutArgs(),
  // however it may be allowable to call create_DfDp_op() before
  // createOutArgs() is called.  Thus we do this here if prototypeOutArgs_
  // needs to be initialized.
  //
  // Previously this was handled in the TEUCHOS_ASSERT below through the call
  // to Np(), however it isn't a good idea to include code in asserts that is
  // required for proper execution (the asserts may be removed in an optimized
  // build, for example).
  if(require_out_args_refresh_) {
    this->createOutArgs();
  }

  TEUCHOS_ASSERT(0<=p && p<Teuchos::as<int>(parameters_.size()));

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
addParameter(const std::string & name,const Scalar & initialValue)
{
  Teuchos::Array<std::string> tmp_names;
  tmp_names.push_back(name);

  Teuchos::Array<Scalar> tmp_values;
  tmp_values.push_back(initialValue);

  return addParameter(tmp_names,tmp_values);
}

template <typename Scalar>
int panzer::ModelEvaluator<Scalar>::
addParameter(const Teuchos::Array<std::string> & names,
             const Teuchos::Array<Scalar> & initialValues)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::ptrFromRef;

  TEUCHOS_ASSERT(names.size()==initialValues.size());

  int parameter_index = parameters_.size();

  // Create parameter object
  RCP<ParameterObject> param = createScalarParameter(names,initialValues);
  parameters_.push_back(param);

  // Create vector space for parameter tangent vector
  RCP< Thyra::VectorSpaceBase<double> > tan_space =
    Thyra::multiVectorProductVectorSpace(x_space_, param->names->size());
  tangent_space_.push_back(tan_space);

  // The number of model evaluator parameters is the number of model parameters (parameters_.size()) plus a tangent
  // vector for each scalar parameter (tangent_space_.size()) plus a tangent vector for xdot for each scalar parameter.
  num_me_parameters_ += 2;
  if (build_transient_support_)
    ++num_me_parameters_;

  require_in_args_refresh_ = true;
  require_out_args_refresh_ = true;
  this->resetDefaultBase();

  return parameter_index;
}

template <typename Scalar>
int panzer::ModelEvaluator<Scalar>::
addDistributedParameter(const std::string & key,
                        const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > & vs,
                        const Teuchos::RCP<GlobalEvaluationData> & ged,
                        const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & initial,
                        const Teuchos::RCP<const GlobalIndexer> & ugi)
{
  distrParamGlobalEvaluationData_.addDataObject(key,ged);

  int parameter_index = parameters_.size();
  parameters_.push_back(createDistributedParameter(key,vs,initial,ugi));
  ++num_me_parameters_;

  require_in_args_refresh_ = true;
  require_out_args_refresh_ = true;
  this->resetDefaultBase();

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
int panzer::ModelEvaluator<Scalar>::
addFlexibleResponse(const std::string & responseName,
            const std::vector<WorksetDescriptor> & wkst_desc,
            const Teuchos::RCP<ResponseMESupportBuilderBase> & builder)
{
   // add a basic response, use x global indexer to define it
   builder->setDerivativeInformation(lof_);

   int respIndex = addResponse(responseName,wkst_desc,*builder);

   // set the builder for this response
   responses_[respIndex]->builder = builder;

   return respIndex;
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
  //TODO: is this really needed?
  ae_inargs.step_size = 0.0;
  ae_inargs.stage_number = 1.0;
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
  RCP<panzer::LinearObjContainer> counter
     = ae_tm_.template getAsObject<panzer::Traits::Residual>()->evaluateOnlyDirichletBCs(ae_inargs);

  // allocate the result container
  RCP<panzer::LinearObjContainer> result = lof_->buildLinearObjContainer(); // we use a new global container

  // stuff the evaluate boundary conditions into the f spot of the counter ... the x is already filled
  Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<Scalar> >(counter)->set_f_th(
        thGlobalContainer->get_f_th());

  // stuff the vector that needs applied dirichlet conditions in the the f spot of the result LOC
  Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<Scalar> >(result)->set_f_th(f);

  // use the linear object factory to apply the result
  lof_->applyDirichletBCs(*counter,*result);
}

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
evalModel_D2gDx2(int respIndex,
                 const Thyra::ModelEvaluatorBase::InArgs<Scalar> & inArgs,
                 const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & delta_x,
                 const Teuchos::RCP<Thyra::VectorBase<Scalar> > & D2gDx2) const
{
#ifdef Panzer_BUILD_HESSIAN_SUPPORT

  // set model parameters from supplied inArgs
  setParameters(inArgs);

  {
    std::string responseName = responses_[respIndex]->name;
    Teuchos::RCP<panzer::ResponseMESupportBase<panzer::Traits::Hessian> > resp
        = Teuchos::rcp_dynamic_cast<panzer::ResponseMESupportBase<panzer::Traits::Hessian> >(
            responseLibrary_->getResponse<panzer::Traits::Hessian>(responseName));
    resp->setDerivative(D2gDx2);
  }

  // setup all the assembly in arguments (this is parameters and
  // x/x_dot). At this point with the exception of the one time dirichlet
  // beta that is all thats neccessary.
  panzer::AssemblyEngineInArgs ae_inargs;
  setupAssemblyInArgs(inArgs,ae_inargs);

  ae_inargs.beta = 1.0;

  auto deltaXContainer = lof_->buildReadOnlyDomainContainer();
  deltaXContainer->setOwnedVector(delta_x);
  ae_inargs.addGlobalEvaluationData("DELTA_Solution Gather Container",deltaXContainer);

  // evaluate responses
  responseLibrary_->addResponsesToInArgs<panzer::Traits::Hessian>(ae_inargs);
  responseLibrary_->evaluate<panzer::Traits::Hessian>(ae_inargs);

  // reset parameters back to nominal values
  resetParameters();
#else
  (void)respIndex;
  (void)inArgs;
  (void)delta_x;
  (void)D2gDx2;
  TEUCHOS_ASSERT(false);
#endif
}

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
evalModel_D2gDxDp(int respIndex,
                  int pIndex,
                  const Thyra::ModelEvaluatorBase::InArgs<Scalar> & inArgs,
                  const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & delta_p,
                  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & D2gDxDp) const
{
#ifdef Panzer_BUILD_HESSIAN_SUPPORT

  // set model parameters from supplied inArgs
  setParameters(inArgs);

  {
    std::string responseName = responses_[respIndex]->name;
    Teuchos::RCP<panzer::ResponseMESupportBase<panzer::Traits::Hessian> > resp
        = Teuchos::rcp_dynamic_cast<panzer::ResponseMESupportBase<panzer::Traits::Hessian> >(
            responseLibrary_->getResponse<panzer::Traits::Hessian>(responseName));
    resp->setDerivative(D2gDxDp);
  }

  // setup all the assembly in arguments (this is parameters and
  // x/x_dot). At this point with the exception of the one time dirichlet
  // beta that is all thats neccessary.
  panzer::AssemblyEngineInArgs ae_inargs;
  setupAssemblyInArgs(inArgs,ae_inargs);

  ae_inargs.beta = 1.0;
  ae_inargs.second_sensitivities_name = (*parameters_[pIndex]->names)[0]; // distributed parameters have one name!

  auto deltaPContainer = parameters_[pIndex]->dfdp_rl->getLinearObjFactory()->buildReadOnlyDomainContainer();
  deltaPContainer->setOwnedVector(delta_p);
  ae_inargs.addGlobalEvaluationData("DELTA_"+(*parameters_[pIndex]->names)[0],deltaPContainer);

  // evaluate responses
  responseLibrary_->addResponsesToInArgs<panzer::Traits::Hessian>(ae_inargs);
  responseLibrary_->evaluate<panzer::Traits::Hessian>(ae_inargs);

  // reset parameters back to nominal values
  resetParameters();
#else
  (void)respIndex;
  (void)pIndex;
  (void)inArgs;
  (void)delta_p;
  (void)D2gDxDp;
  TEUCHOS_ASSERT(false);
#endif
}

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
evalModel_D2gDp2(int respIndex,
                 int pIndex,
                 const Thyra::ModelEvaluatorBase::InArgs<Scalar> & inArgs,
                 const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & delta_p,
                 const Teuchos::RCP<Thyra::VectorBase<Scalar> > & D2gDp2) const
{
#ifdef Panzer_BUILD_HESSIAN_SUPPORT

  // set model parameters from supplied inArgs
  setParameters(inArgs);

  ResponseLibrary<Traits> & rLibrary = *parameters_[pIndex]->dgdp_rl;

  {
    std::string responseName = responses_[respIndex]->name;
    Teuchos::RCP<panzer::ResponseMESupportBase<panzer::Traits::Hessian> > resp
        = Teuchos::rcp_dynamic_cast<panzer::ResponseMESupportBase<panzer::Traits::Hessian> >(
            rLibrary.getResponse<panzer::Traits::Hessian>(responseName));
    resp->setDerivative(D2gDp2);
  }

  // setup all the assembly in arguments (this is parameters and
  // x/x_dot). At this point with the exception of the one time dirichlet
  // beta that is all thats neccessary.
  panzer::AssemblyEngineInArgs ae_inargs;
  setupAssemblyInArgs(inArgs,ae_inargs);

  ae_inargs.gather_seeds.push_back(1.0); // this assumes that gather point is always the zero index of
                                         // gather seeds
  ae_inargs.first_sensitivities_name  = (*parameters_[pIndex]->names)[0]; // distributed parameters have one name!
  ae_inargs.second_sensitivities_name = (*parameters_[pIndex]->names)[0]; // distributed parameters have one name!

  auto deltaPContainer = parameters_[pIndex]->dfdp_rl->getLinearObjFactory()->buildReadOnlyDomainContainer();
  deltaPContainer->setOwnedVector(delta_p);
  ae_inargs.addGlobalEvaluationData("DELTA_"+(*parameters_[pIndex]->names)[0],deltaPContainer);

  // evaluate responses
  rLibrary.addResponsesToInArgs<panzer::Traits::Hessian>(ae_inargs);
  rLibrary.evaluate<panzer::Traits::Hessian>(ae_inargs);

  // reset parameters back to nominal values
  resetParameters();
#else
  (void)respIndex;
  (void)pIndex;
  (void)inArgs;
  (void)delta_p;
  (void)D2gDp2;
  TEUCHOS_ASSERT(false);
#endif
}

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
evalModel_D2gDpDx(int respIndex,
                  int pIndex,
                  const Thyra::ModelEvaluatorBase::InArgs<Scalar> & inArgs,
                  const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & delta_x,
                  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & D2gDpDx) const
{
#ifdef Panzer_BUILD_HESSIAN_SUPPORT

  // set model parameters from supplied inArgs
  setParameters(inArgs);

  ResponseLibrary<Traits> & rLibrary = *parameters_[pIndex]->dgdp_rl;

  {
    std::string responseName = responses_[respIndex]->name;
    Teuchos::RCP<panzer::ResponseMESupportBase<panzer::Traits::Hessian> > resp
        = Teuchos::rcp_dynamic_cast<panzer::ResponseMESupportBase<panzer::Traits::Hessian> >(
            rLibrary.getResponse<panzer::Traits::Hessian>(responseName));
    resp->setDerivative(D2gDpDx);
  }

  // setup all the assembly in arguments (this is parameters and
  // x/x_dot). At this point with the exception of the one time dirichlet
  // beta that is all thats neccessary.
  panzer::AssemblyEngineInArgs ae_inargs;
  setupAssemblyInArgs(inArgs,ae_inargs);

  ae_inargs.gather_seeds.push_back(1.0); // this assumes that gather point is always the zero index of
                                         // gather seeds
  ae_inargs.first_sensitivities_name  = (*parameters_[pIndex]->names)[0]; // distributed parameters have one name!
  ae_inargs.second_sensitivities_name  = "";

  auto deltaXContainer = lof_->buildReadOnlyDomainContainer();
  deltaXContainer->setOwnedVector(delta_x);
  ae_inargs.addGlobalEvaluationData("DELTA_Solution Gather Container",deltaXContainer);

  // evaluate responses
  rLibrary.addResponsesToInArgs<panzer::Traits::Hessian>(ae_inargs);
  rLibrary.evaluate<panzer::Traits::Hessian>(ae_inargs);

  // reset parameters back to nominal values
  resetParameters();
#else
  (void)respIndex;
  (void)pIndex;
  (void)inArgs;
  (void)delta_x;
  (void)D2gDpDx;
  TEUCHOS_ASSERT(false);
#endif
}

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
evalModel_D2fDx2(const Thyra::ModelEvaluatorBase::InArgs<Scalar> & inArgs,
                 const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & delta_x,
                 const Teuchos::RCP<Thyra::LinearOpBase<Scalar> > & D2fDx2) const
{
#ifdef Panzer_BUILD_HESSIAN_SUPPORT

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
  const RCP<Thyra::LinearOpBase<Scalar> > W_out = D2fDx2;

  // setup all the assembly in arguments (this is parameters and
  // x/x_dot). At this point with the exception of the one time dirichlet
  // beta that is all thats neccessary.
  panzer::AssemblyEngineInArgs ae_inargs;
  setupAssemblyInArgs(inArgs,ae_inargs);

  auto deltaXContainer = lof_->buildReadOnlyDomainContainer();
  deltaXContainer->setOwnedVector(delta_x);
  ae_inargs.addGlobalEvaluationData("DELTA_Solution Gather Container",deltaXContainer);

  // set model parameters from supplied inArgs
  setParameters(inArgs);

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

  {
    PANZER_FUNC_TIME_MONITOR("panzer::ModelEvaluator::evalModel(D2fDx2)");

    // this dummy nonsense is needed only for scattering dirichlet conditions
    RCP<Thyra::VectorBase<Scalar> > dummy_f = Thyra::createMember(f_space_);
    thGlobalContainer->set_f_th(dummy_f);
    thGlobalContainer->set_A_th(W_out);

    // Zero values in ghosted container objects
    thGhostedContainer->initializeMatrix(0.0);

    ae_tm_.template getAsObject<panzer::Traits::Hessian>()->evaluate(ae_inargs);
  }

  // HACK: set A to null before calling responses to avoid touching the
  // the Jacobian after it has been properly assembled.  Should be fixed
  // by using a modified version of ae_inargs instead.
  thGlobalContainer->set_A_th(Teuchos::null);

  // TODO: Clearing all references prevented a seg-fault with Rythmos,
  // which is no longer used. Check if it's still needed.
  thGlobalContainer->set_x_th(Teuchos::null);
  thGlobalContainer->set_dxdt_th(Teuchos::null);
  thGlobalContainer->set_f_th(Teuchos::null);
  thGlobalContainer->set_A_th(Teuchos::null);

  // reset parameters back to nominal values
  resetParameters();
#else
  (void)inArgs;
  (void)delta_x;
  (void)D2fDx2;
  TEUCHOS_ASSERT(false);
#endif
}

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
evalModel_D2fDxDp(int pIndex,
                  const Thyra::ModelEvaluatorBase::InArgs<Scalar> & inArgs,
                  const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & delta_p,
                  const Teuchos::RCP<Thyra::LinearOpBase<Scalar> > & D2fDxDp) const
{
#ifdef Panzer_BUILD_HESSIAN_SUPPORT

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
  const RCP<Thyra::LinearOpBase<Scalar> > W_out = D2fDxDp;

  // setup all the assembly in arguments (this is parameters and
  // x/x_dot). At this point with the exception of the one time dirichlet
  // beta that is all thats neccessary.
  panzer::AssemblyEngineInArgs ae_inargs;
  setupAssemblyInArgs(inArgs,ae_inargs);

  ae_inargs.second_sensitivities_name = (*parameters_[pIndex]->names)[0]; // distributed parameters have one name!

  auto deltaPContainer = parameters_[pIndex]->dfdp_rl->getLinearObjFactory()->buildReadOnlyDomainContainer();
  deltaPContainer->setOwnedVector(delta_p);
  ae_inargs.addGlobalEvaluationData("DELTA_"+(*parameters_[pIndex]->names)[0],deltaPContainer);

  // set model parameters from supplied inArgs
  setParameters(inArgs);

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

  {
    PANZER_FUNC_TIME_MONITOR("panzer::ModelEvaluator::evalModel(D2fDxDp)");

    // this dummy nonsense is needed only for scattering dirichlet conditions
    RCP<Thyra::VectorBase<Scalar> > dummy_f = Thyra::createMember(f_space_);
    thGlobalContainer->set_f_th(dummy_f);
    thGlobalContainer->set_A_th(W_out);

    // Zero values in ghosted container objects
    thGhostedContainer->initializeMatrix(0.0);

    ae_tm_.template getAsObject<panzer::Traits::Hessian>()->evaluate(ae_inargs);
  }

  // HACK: set A to null before calling responses to avoid touching the
  // the Jacobian after it has been properly assembled.  Should be fixed
  // by using a modified version of ae_inargs instead.
  thGlobalContainer->set_A_th(Teuchos::null);

  // TODO: Clearing all references prevented a seg-fault with Rythmos,
  // which is no longer used. Check if it's still needed.
  thGlobalContainer->set_x_th(Teuchos::null);
  thGlobalContainer->set_dxdt_th(Teuchos::null);
  thGlobalContainer->set_f_th(Teuchos::null);
  thGlobalContainer->set_A_th(Teuchos::null);

  // reset parameters back to nominal values
  resetParameters();
#else
  (void)pIndex;
  (void)inArgs;
  (void)delta_p;
  (void)D2fDxDp;
  TEUCHOS_ASSERT(false);
#endif
}

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
evalModel_D2fDpDx(int pIndex,
                  const Thyra::ModelEvaluatorBase::InArgs<Scalar> & inArgs,
                  const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & delta_x,
                  const Teuchos::RCP<Thyra::LinearOpBase<Scalar> > & D2fDpDx) const
{
#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::null;

  // parameter is not distributed
  TEUCHOS_ASSERT(parameters_[pIndex]->is_distributed);

  // parameter is distributed but has no global indexer.
  // thus the user doesn't want sensitivities!
  TEUCHOS_ASSERT(parameters_[pIndex]->dfdp_rl!=null);

  ResponseLibrary<Traits> & rLibrary = *parameters_[pIndex]->dfdp_rl;

  // get the response and tell it to fill the derivative operator
  RCP<Response_Residual<Traits::Hessian> > response_hessian =
    rcp_dynamic_cast<Response_Residual<Traits::Hessian> >(rLibrary.getResponse<Traits::Hessian>("RESIDUAL"));
  response_hessian->setHessian(D2fDpDx);

  // setup all the assembly in arguments (this is parameters and x/x_dot).
  // make sure the correct seeding is performed
  panzer::AssemblyEngineInArgs ae_inargs;
  setupAssemblyInArgs(inArgs,ae_inargs);

  auto deltaXContainer = lof_->buildReadOnlyDomainContainer();
  deltaXContainer->setOwnedVector(delta_x);
  ae_inargs.addGlobalEvaluationData("DELTA_Solution Gather Container",deltaXContainer);

  ae_inargs.gather_seeds.push_back(1.0); // this assumes that gather point is always the zero index of
                                         // gather seeds
  ae_inargs.first_sensitivities_name  = (*parameters_[pIndex]->names)[0]; // distributed parameters have one name!
  ae_inargs.second_sensitivities_name  = "";

  rLibrary.addResponsesToInArgs<Traits::Hessian>(ae_inargs);
  rLibrary.evaluate<Traits::Hessian>(ae_inargs);
#else
  (void)pIndex;
  (void)inArgs;
  (void)delta_x;
  (void)D2fDpDx;
  TEUCHOS_ASSERT(false);
#endif
}

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
evalModel_D2fDp2(int pIndex,
                 const Thyra::ModelEvaluatorBase::InArgs<Scalar> & inArgs,
                 const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & delta_p,
                 const Teuchos::RCP<Thyra::LinearOpBase<Scalar> > & D2fDp2) const
{
#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::null;

  // parameter is not distributed
  TEUCHOS_ASSERT(parameters_[pIndex]->is_distributed);

  // parameter is distributed but has no global indexer.
  // thus the user doesn't want sensitivities!
  TEUCHOS_ASSERT(parameters_[pIndex]->dfdp_rl!=null);

  ResponseLibrary<Traits> & rLibrary = *parameters_[pIndex]->dfdp_rl;

  // get the response and tell it to fill the derivative operator
  RCP<Response_Residual<Traits::Hessian> > response_hessian =
    rcp_dynamic_cast<Response_Residual<Traits::Hessian> >(rLibrary.getResponse<Traits::Hessian>("RESIDUAL"));
  response_hessian->setHessian(D2fDp2);

  // setup all the assembly in arguments (this is parameters and x/x_dot).
  // make sure the correct seeding is performed
  panzer::AssemblyEngineInArgs ae_inargs;
  setupAssemblyInArgs(inArgs,ae_inargs);

  auto deltaPContainer = parameters_[pIndex]->dfdp_rl->getLinearObjFactory()->buildReadOnlyDomainContainer();
  deltaPContainer->setOwnedVector(delta_p);
  ae_inargs.addGlobalEvaluationData("DELTA_"+(*parameters_[pIndex]->names)[0],deltaPContainer);

  ae_inargs.gather_seeds.push_back(1.0); // this assumes that gather point is always the zero index of
                                         // gather seeds
  ae_inargs.first_sensitivities_name  = (*parameters_[pIndex]->names)[0]; // distributed parameters have one name!
  ae_inargs.second_sensitivities_name = (*parameters_[pIndex]->names)[0]; // distributed parameters have one name!

  rLibrary.addResponsesToInArgs<Traits::Hessian>(ae_inargs);
  rLibrary.evaluate<Traits::Hessian>(ae_inargs);
#else
  (void)pIndex;
  (void)inArgs;
  (void)delta_p;
  (void)D2fDp2;
  TEUCHOS_ASSERT(false);
#endif
}

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
              const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  evalModelImpl_basic(inArgs,outArgs);

  // evaluate responses...uses the stored assembly arguments and containers
  if(required_basic_g(outArgs))
    evalModelImpl_basic_g(inArgs,outArgs);

  // evaluate response derivatives
  if(required_basic_dgdx(outArgs))
    evalModelImpl_basic_dgdx(inArgs,outArgs);

  // evaluate response derivatives to scalar parameters
  if(required_basic_dgdp_scalar(outArgs))
    evalModelImpl_basic_dgdp_scalar(inArgs,outArgs);

  // evaluate response derivatives to distributed parameters
  if(required_basic_dgdp_distro(outArgs))
    evalModelImpl_basic_dgdp_distro(inArgs,outArgs);

  if(required_basic_dfdp_scalar(outArgs)) {
    if (do_fd_dfdp_)
      evalModelImpl_basic_dfdp_scalar_fd(inArgs,outArgs);
    else
      evalModelImpl_basic_dfdp_scalar(inArgs,outArgs);
  }

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

  // set model parameters from supplied inArgs
  setParameters(inArgs);

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

    // only add auxiliary global data if Jacobian is being formed
    ae_inargs.addGlobalEvaluationData(nonParamGlobalEvaluationData_);

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

    // don't add auxiliary global data if Jacobian is not computed.
    // this leads to zeroing of aux ops in special cases.

    thGlobalContainer->set_f_th(f_out);

    // Zero values in ghosted container objects
    Thyra::assign(thGhostedContainer->get_f_th().ptr(),0.0);

    ae_tm_.template getAsObject<panzer::Traits::Residual>()->evaluate(ae_inargs);
  }
  else if(Teuchos::is_null(f_out) && !Teuchos::is_null(W_out)) {

    PANZER_FUNC_TIME_MONITOR("panzer::ModelEvaluator::evalModel(J)");

    // only add auxiliary global data if Jacobian is being formed
    ae_inargs.addGlobalEvaluationData(nonParamGlobalEvaluationData_);

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

  // TODO: Clearing all references prevented a seg-fault with Rythmos,
  // which is no longer used. Check if it's still needed.
  thGlobalContainer->set_x_th(Teuchos::null);
  thGlobalContainer->set_dxdt_th(Teuchos::null);
  thGlobalContainer->set_f_th(Teuchos::null);
  thGlobalContainer->set_A_th(Teuchos::null);

  // reset parameters back to nominal values
  resetParameters();

  const bool writeToFile = false;
  if (writeToFile && nonnull(W_out)) {
    const auto check_blocked = Teuchos::rcp_dynamic_cast<::Thyra::BlockedLinearOpBase<double> >(W_out,false);
    if (check_blocked) {
      const int numBlocks = check_blocked->productDomain()->numBlocks();
      const int rangeBlocks = check_blocked->productRange()->numBlocks();
      TEUCHOS_ASSERT(numBlocks == rangeBlocks); // not true for optimization
      for (int row=0; row < numBlocks; ++row) {
        for (int col=0; col < numBlocks; ++col) {
          using LO = panzer::LocalOrdinal;
          using GO = panzer::GlobalOrdinal;
          using NodeT = panzer::TpetraNodeType;
          const auto thyraTpetraOperator = Teuchos::rcp_dynamic_cast<::Thyra::TpetraLinearOp<double,LO,GO,NodeT>>(check_blocked->getNonconstBlock(row,col),true);
          const auto tpetraCrsMatrix = Teuchos::rcp_dynamic_cast<Tpetra::CrsMatrix<double,LO,GO,NodeT>>(thyraTpetraOperator->getTpetraOperator(),true);
          tpetraCrsMatrix->print(std::cout);
          std::stringstream ss;
          ss << "W_out_" << write_matrix_count_ << ".rank_" << tpetraCrsMatrix->getMap()->getComm()->getRank() << ".block_" << row << "_" << col << ".txt";
          std::fstream fs(ss.str().c_str(),std::fstream::out|std::fstream::trunc);
          Teuchos::FancyOStream fos(Teuchos::rcpFromRef(fs));
          tpetraCrsMatrix->describe(fos,Teuchos::VERB_EXTREME);
          fs.close();
        }
      }
    }
    else {
      using LO = panzer::LocalOrdinal;
      using GO = panzer::GlobalOrdinal;
      using NodeT = panzer::TpetraNodeType;
      const auto thyraTpetraOperator = Teuchos::rcp_dynamic_cast<::Thyra::TpetraLinearOp<double,LO,GO,NodeT>>(W_out,true);
      const auto tpetraCrsMatrix = Teuchos::rcp_dynamic_cast<Tpetra::CrsMatrix<double,LO,GO,NodeT>>(thyraTpetraOperator->getTpetraOperator(),true);
      tpetraCrsMatrix->print(std::cout);
      std::stringstream ss;
      ss << "W_out_" << write_matrix_count_ << ".rank_" << tpetraCrsMatrix->getMap()->getComm()->getRank() << ".txt";
      std::fstream fs(ss.str().c_str(),std::fstream::out|std::fstream::trunc);
      Teuchos::FancyOStream fos(Teuchos::rcpFromRef(fs));
      tpetraCrsMatrix->describe(fos,Teuchos::VERB_EXTREME);
      fs.close();
    }
    ++write_matrix_count_;
  }

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

  // set model parameters from supplied inArgs
  setParameters(inArgs);

  for(std::size_t i=0;i<responses_.size();i++) {
    Teuchos::RCP<Thyra::VectorBase<Scalar> > vec = outArgs.get_g(i);
    if(vec!=Teuchos::null) {
      std::string responseName = responses_[i]->name;
      Teuchos::RCP<panzer::ResponseMESupportBase<panzer::Traits::Residual> > resp
          = Teuchos::rcp_dynamic_cast<panzer::ResponseMESupportBase<panzer::Traits::Residual> >(
              responseLibrary_->getResponse<panzer::Traits::Residual>(responseName));
      resp->setVector(vec);
    }
  }

  // evaluator responses
  responseLibrary_->addResponsesToInArgs<panzer::Traits::Residual>(ae_inargs);
  responseLibrary_->evaluate<panzer::Traits::Residual>(ae_inargs);

  // reset parameters back to nominal values
  resetParameters();
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

  // set model parameters from supplied inArgs
  setParameters(inArgs);

  for(std::size_t i=0;i<responses_.size();i++) {
    // get "Vector" out of derivative, if its something else, throw an exception
    MEB::Derivative<Scalar> deriv = outArgs.get_DgDx(i);
    if(deriv.isEmpty())
      continue;

    Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > vec = deriv.getMultiVector();

    if(vec!=Teuchos::null) {

      std::string responseName = responses_[i]->name;
      Teuchos::RCP<panzer::ResponseMESupportBase<panzer::Traits::Jacobian> > resp
          = Teuchos::rcp_dynamic_cast<panzer::ResponseMESupportBase<panzer::Traits::Jacobian> >(
              responseLibrary_->getResponse<panzer::Traits::Jacobian>(responseName));
      resp->setDerivative(vec);
    }
  }

  // setup all the assembly in arguments (this is parameters and
  // x/x_dot). At this point with the exception of the one time dirichlet
  // beta that is all thats neccessary.
  panzer::AssemblyEngineInArgs ae_inargs;
  setupAssemblyInArgs(inArgs,ae_inargs);

  // evaluate responses
  responseLibrary_->addResponsesToInArgs<panzer::Traits::Jacobian>(ae_inargs);
  responseLibrary_->evaluate<panzer::Traits::Jacobian>(ae_inargs);

  // reset parameters back to nominal values
  resetParameters();
}

template <typename Scalar>
void
panzer::ModelEvaluator<Scalar>::
evalModelImpl_basic_dgdp_scalar(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
                                const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  typedef Thyra::ModelEvaluatorBase MEB;

  // optional sanity check
  TEUCHOS_ASSERT(required_basic_dgdp_scalar(outArgs));

  // First find all of the active parameters for all responses
  std::vector<std::string> activeParameterNames;
  std::vector<int> activeParameters;
  int totalParameterCount = 0;
  for(std::size_t j=0; j<parameters_.size(); j++) {

    // skip non-scalar parameters
    if(parameters_[j]->is_distributed)
      continue;

    bool is_active = false;
    for(std::size_t i=0;i<responses_.size(); i++) {

      MEB::Derivative<Scalar> deriv = outArgs.get_DgDp(i,j);
      if(deriv.isEmpty())
        continue;

      Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > vec = deriv.getMultiVector();
      if(vec!=Teuchos::null) {
        // get the response and tell it to fill the derivative vector
        std::string responseName = responses_[i]->name;
        RCP<panzer::ResponseMESupportBase<panzer::Traits::Tangent> > resp =
          rcp_dynamic_cast<panzer::ResponseMESupportBase<panzer::Traits::Tangent> >(
            responseLibrary_->getResponse<panzer::Traits::Tangent>(responseName));
        resp->setVector(vec);
        is_active = true;
      }
    }

    if (is_active) {
      for (std::size_t k=0; k<parameters_[j]->scalar_value.size(); k++) {
        std::string name = "PARAMETER_SENSITIVIES: "+(*parameters_[j]->names)[k];
        activeParameterNames.push_back(name);
        totalParameterCount++;
      }
      activeParameters.push_back(j);
    }
  }

  // setup all the assembly in arguments
  panzer::AssemblyEngineInArgs ae_inargs;
  setupAssemblyInArgs(inArgs,ae_inargs);

  // add active parameter names to assembly in-args
  RCP<panzer::GlobalEvaluationData> ged_activeParameters =
    rcp(new panzer::ParameterList_GlobalEvaluationData(activeParameterNames));
  ae_inargs.addGlobalEvaluationData("PARAMETER_NAMES",ged_activeParameters);

  // Initialize Fad components of all active parameters
  int paramIndex = 0;
  for (std::size_t ap=0; ap<activeParameters.size(); ++ap) {
    const int j = activeParameters[ap];
    for (unsigned int k=0; k < parameters_[j]->scalar_value.size(); k++) {
      panzer::Traits::FadType p(totalParameterCount, parameters_[j]->scalar_value[k].baseValue);
      p.fastAccessDx(paramIndex) = 1.0;
      parameters_[j]->scalar_value[k].family->template setValue<panzer::Traits::Tangent>(p);
      paramIndex++;
    }
  }

  // make sure that the total parameter count and the total parameter index match up
  TEUCHOS_ASSERT(paramIndex==totalParameterCount);

  // evaluate response tangent
  if(totalParameterCount>0) {
    responseLibrary_->addResponsesToInArgs<Traits::Tangent>(ae_inargs);
    responseLibrary_->evaluate<Traits::Tangent>(ae_inargs);
  }
}

template <typename Scalar>
void
panzer::ModelEvaluator<Scalar>::
evalModelImpl_basic_dgdp_distro(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
                                const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  // optional sanity check
  TEUCHOS_ASSERT(required_basic_dgdp_distro(outArgs));

  // loop over parameters, and then build a dfdp_rl only if they are distributed
  // and the user has provided the UGI. Note that this may be overly expensive if they
  // don't actually want those sensitivites because memory will be allocated unneccesarily.
  // It would be good to do this "just in time", but for now this is sufficient.
  for(std::size_t p=0;p<parameters_.size();p++) {

    // parameter is not distributed, a different path is
    // taken for those to compute dfdp
    if(!parameters_[p]->is_distributed)
      continue;

    ResponseLibrary<Traits> & rLibrary = *parameters_[p]->dgdp_rl;

    for(std::size_t r=0;r<responses_.size();r++) {
      // have derivatives been requested?
      MEB::Derivative<Scalar> deriv = outArgs.get_DgDp(r,p);
      if(deriv.isEmpty())
        continue;

      Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > vec = deriv.getMultiVector();

      if(vec!=Teuchos::null) {

        // get the response and tell it to fill the derivative vector
        std::string responseName = responses_[r]->name;
        Teuchos::RCP<panzer::ResponseMESupportBase<panzer::Traits::Jacobian> > resp
            = Teuchos::rcp_dynamic_cast<panzer::ResponseMESupportBase<panzer::Traits::Jacobian> >(
                rLibrary.getResponse<panzer::Traits::Jacobian>(responseName));

        resp->setDerivative(vec);
      }
    }

    // setup all the assembly in arguments (this is parameters and x/x_dot).
    // make sure the correct seeding is performed
    panzer::AssemblyEngineInArgs ae_inargs;
    setupAssemblyInArgs(inArgs,ae_inargs);

    ae_inargs.first_sensitivities_name = (*parameters_[p]->names)[0]; // distributed parameters have one name!
    ae_inargs.gather_seeds.push_back(1.0); // this assumes that gather point is always the zero index of
                                           // gather seeds

    // evaluate responses
    rLibrary.addResponsesToInArgs<Traits::Jacobian>(ae_inargs);
    rLibrary.evaluate<Traits::Jacobian>(ae_inargs);
  }
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
       RCP<LOCPair_GlobalEvaluationData> loc_pair
           = Teuchos::rcp(new LOCPair_GlobalEvaluationData(lof_,LinearObjContainer::F));
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

   RCP<GlobalEvaluationData> ged_activeParameters
       = Teuchos::rcp(new ParameterList_GlobalEvaluationData(activeParameters));
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
         Traits::FadType p = Traits::FadType(totalParameterCount,
                                             parameters_[i]->scalar_value[j].baseValue);
         parameters_[i]->scalar_value[j].family->template setValue<panzer::Traits::Tangent>(p);
       }
       continue;
     }
     else {
       // loop over each parameter in the vector, initializing the AD type
       for (unsigned int j=0; j < parameters_[i]->scalar_value.size(); j++) {
         Traits::FadType p = Traits::FadType(totalParameterCount,
                                             parameters_[i]->scalar_value[j].baseValue);
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
     PANZER_FUNC_TIME_MONITOR("panzer::ModelEvaluator::evalModel(df/dp)");
     ae_tm_.getAsObject<panzer::Traits::Tangent>()->evaluate(ae_inargs);
   }
}

template <typename Scalar>
void
panzer::ModelEvaluator<Scalar>::
evalModelImpl_basic_dfdp_scalar_fd(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
                                   const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
   PANZER_FUNC_TIME_MONITOR("panzer::ModelEvaluator::evalModel(df/dp)");

   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   typedef Thyra::ModelEvaluatorBase MEB;

   TEUCHOS_ASSERT(required_basic_dfdp_scalar(outArgs));

   // First evaluate the model without df/dp for the base point
   // Maybe it would be better to set all outArgs and then remove the df/dp ones,
   // but I couldn't get that to work.
   MEB::OutArgs<Scalar> outArgs_base = this->createOutArgs();
   if (outArgs.get_f() == Teuchos::null)
     outArgs_base.set_f(Thyra::createMember(this->get_f_space()));
   else
     outArgs_base.set_f(outArgs.get_f());
   outArgs_base.set_W_op(outArgs.get_W_op());
   this->evalModel(inArgs, outArgs_base);
   RCP<const Thyra::VectorBase<Scalar> > f = outArgs_base.get_f();
   RCP<const Thyra::VectorBase<Scalar> > x = inArgs.get_x();
   RCP<const Thyra::VectorBase<Scalar> > x_dot;
   if (inArgs.supports(MEB::IN_ARG_x_dot))
     x_dot = inArgs.get_x_dot();

   // Create in/out args for FD calculation
   RCP<Thyra::VectorBase<Scalar> > fd = Thyra::createMember(this->get_f_space());
   MEB::OutArgs<Scalar> outArgs_fd = this->createOutArgs();
   outArgs_fd.set_f(fd);

   RCP<Thyra::VectorBase<Scalar> > xd = Thyra::createMember(this->get_x_space());
   RCP<Thyra::VectorBase<Scalar> > xd_dot;
   if (x_dot != Teuchos::null)
     xd_dot = Thyra::createMember(this->get_x_space());
   MEB::InArgs<Scalar> inArgs_fd = this->createInArgs();
   inArgs_fd.setArgs(inArgs);  // This sets all inArgs that we don't override below
   inArgs_fd.set_x(xd);
   if (x_dot != Teuchos::null)
     inArgs_fd.set_x_dot(xd_dot);

   const double h = fd_perturb_size_;
   for(std::size_t i=0; i < parameters_.size(); i++) {

     // skip non-scalar parameters
     if(parameters_[i]->is_distributed)
       continue;

     // have derivatives been requested?
     MEB::Derivative<Scalar> deriv = outArgs.get_DfDp(i);
     if(deriv.isEmpty())
       continue;

     // grab multivector, make sure its the right dimension
     RCP<Thyra::MultiVectorBase<Scalar> > dfdp = deriv.getMultiVector();
     TEUCHOS_ASSERT(dfdp->domain()->dim()==Teuchos::as<int>(parameters_[i]->scalar_value.size()));

     // Get parameter vector and tangent vectors
     RCP<const Thyra::VectorBase<Scalar> > p       = inArgs.get_p(i);
     RCP<const Thyra::VectorBase<Scalar> > dx_v    = inArgs.get_p(i+parameters_.size());
     RCP<const Thyra::MultiVectorBase<Scalar> > dx =
       rcp_dynamic_cast<const Thyra::DefaultMultiVectorProductVector<Scalar> >(dx_v,true)->getMultiVector();
     RCP<const Thyra::VectorBase<Scalar> > dx_dot_v;
     RCP<const Thyra::MultiVectorBase<Scalar> > dx_dot;
     if (x_dot != Teuchos::null) {
       dx_dot_v =inArgs.get_p(i+parameters_.size()+tangent_space_.size());
       dx_dot =
         rcp_dynamic_cast<const Thyra::DefaultMultiVectorProductVector<Scalar> >(dx_dot_v,true)->getMultiVector();
     }

     // Create perturbed parameter vector
     RCP<Thyra::VectorBase<Scalar> > pd = Thyra::createMember(this->get_p_space(i));
     inArgs_fd.set_p(i,pd);

     for (std::size_t j=0; j < parameters_[i]->scalar_value.size(); j++) {

       // Perturb parameter vector
       Thyra::copy(*p, pd.ptr());
       Thyra::set_ele(j, Thyra::get_ele(*p,j)+h, pd.ptr());

       // Perturb state vectors using tangents
       Thyra::V_VpStV(xd.ptr(), *x, h, *(dx)->col(j));
       if (x_dot != Teuchos::null)
         Thyra::V_VpStV(xd_dot.ptr(), *x_dot, h, *(dx_dot)->col(j));

       // Evaluate perturbed residual
       Thyra::assign(fd.ptr(), 0.0);
       this->evalModel(inArgs_fd, outArgs_fd);

       // FD calculation
       Thyra::V_StVpStV(dfdp->col(j).ptr(), 1.0/h, *fd, -1.0/h, *f);

       // Reset parameter back to un-perturbed value
       parameters_[i]->scalar_value[j].family->setRealValueForAllTypes(Thyra::get_ele(*p,j));

     }
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

    ae_inargs.first_sensitivities_name = (*parameters_[p]->names)[0]; // distributed parameters have one name!
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
required_basic_dgdp_scalar(const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
   typedef Thyra::ModelEvaluatorBase MEB;

   // determine if any of the outArgs are not null!
   bool activeGArgs = false;
   for(int i=0;i<outArgs.Ng();i++) {
     for(int p=0;p<Teuchos::as<int>(parameters_.size());p++) {

       // only look at scalar parameters
       if(parameters_[p]->is_distributed)
         continue;

       // no derivatives are supported
       if(outArgs.supports(MEB::OUT_ARG_DgDp,i,p).none())
         continue;

       activeGArgs |= (!outArgs.get_DgDp(i,p).isEmpty());
     }
   }

   return activeGArgs;
}

template <typename Scalar>
bool panzer::ModelEvaluator<Scalar>::
required_basic_dgdp_distro(const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
   typedef Thyra::ModelEvaluatorBase MEB;

   // determine if any of the outArgs are not null!
   bool activeGArgs = false;
   for(int i=0;i<outArgs.Ng();i++) {
     for(int p=0;p<Teuchos::as<int>(parameters_.size());p++) {

       // only look at distributed parameters
       if(!parameters_[p]->is_distributed)
         continue;

       // no derivatives are supported
       if(outArgs.supports(MEB::OUT_ARG_DgDp,i,p).none())
         continue;

       activeGArgs |= (!outArgs.get_DgDp(i,p).isEmpty());
     }
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
   for(int i=0;i<Teuchos::as<int>(parameters_.size());i++) {

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
   for(int i=0;i<Teuchos::as<int>(parameters_.size());i++) {

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
    RCP<const LinearObjFactory<Traits> > param_lof = cloneWithNewDomain(*lof_,
                                                                        parameters_[p]->global_indexer);

    // the user wants global sensitivities, hooray! Build and setup the response library
    RCP<ResponseLibrary<Traits> > rLibrary
        = Teuchos::rcp(new ResponseLibrary<Traits>(wc,lof_->getRangeGlobalIndexer(),
                                                   param_lof,true));
    rLibrary->buildResidualResponseEvaluators(physicsBlocks,eqset_factory,bcs,bc_factory,
                                              cm_factory,closure_models,user_data,
                                              write_graphviz_file,graphviz_file_prefix);

    // make sure parameter response library is correct
    parameters_[p]->dfdp_rl = rLibrary;
  }
}

template <typename Scalar>
void panzer::ModelEvaluator<Scalar>::
buildDistroParamDgDp_RL(
       const Teuchos::RCP<panzer::WorksetContainer> & wc,
       const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
       const std::vector<panzer::BC>& /* bcs */,
       const panzer::EquationSetFactory & eqset_factory,
       const panzer::BCStrategyFactory& /* bc_factory */,
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

    // extract the linear object factory that has the correct sizing for
    // the sensitivity vector
    RCP<const LinearObjFactory<Traits> > param_lof = parameters_[p]->dfdp_rl->getLinearObjFactory();
    RCP<const GlobalIndexer > param_ugi = parameters_[p]->global_indexer;

    // the user wants global sensitivities, hooray! Build and setup the response library
    RCP<ResponseLibrary<Traits> > rLibrary
        = Teuchos::rcp(new ResponseLibrary<Traits>(wc,lof_->getRangeGlobalIndexer(), lof_));


    // build evaluators for all flexible responses
    for(std::size_t r=0;r<responses_.size();r++) {
      // only responses with a builder are non null!
      if(responses_[r]->builder==Teuchos::null)
        continue;

      // set the current derivative information in the builder
      // responses_[r]->builder->setDerivativeInformationBase(param_lof,param_ugi);
      responses_[r]->builder->setDerivativeInformation(param_lof);

      // add the response
      rLibrary->addResponse(responses_[r]->name,
                            responses_[r]->wkst_desc,
                            *responses_[r]->builder);
    }

    rLibrary->buildResponseEvaluators(physicsBlocks,eqset_factory,
                                      cm_factory,closure_models,user_data,
                                      write_graphviz_file,graphviz_file_prefix);

    // make sure parameter response library is correct
    parameters_[p]->dgdp_rl = rLibrary;
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
createScalarParameter(const Teuchos::Array<std::string> & in_names,
                      const Teuchos::Array<Scalar> & in_values) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::ptrFromRef;

  TEUCHOS_ASSERT(in_names.size()==in_values.size());

  // Check that the parameters are valid (i.e., they already exist in the parameter library)
  // std::size_t np = in_names.size();
  // for(std::size_t i=0;i<np;i++)
  //   TEUCHOS_TEST_FOR_EXCEPTION(!global_data_->pl->isParameter(in_names[i]),
  //                              std::logic_error,
  //                              "Parameter \"" << in_names[i] << "\" does not exist in parameter library!");

  RCP<ParameterObject> paramObj = rcp(new ParameterObject);

  paramObj->names = rcp(new Teuchos::Array<std::string>(in_names));
  paramObj->is_distributed = false;

  // register all the scalar parameters, setting initial
  for(int i=0;i<in_names.size();i++)
    registerScalarParameter(in_names[i],*global_data_->pl,in_values[i]);

  paramObj->scalar_value = panzer::ParamVec();
  global_data_->pl->fillVector<panzer::Traits::Residual>(*paramObj->names, paramObj->scalar_value);

  // build initial condition vector
  paramObj->space =
    Thyra::locallyReplicatedDefaultSpmdVectorSpace<Scalar>(
      rcp(new Teuchos::MpiComm<long int>(lof_->getComm().getRawMpiComm())),paramObj->names->size());

  // fill vector with parameter values
  Teuchos::ArrayRCP<Scalar> data;
  RCP<Thyra::VectorBase<Scalar> > initial_value = Thyra::createMember(paramObj->space);
  RCP<Thyra::SpmdVectorBase<Scalar> > vec = rcp_dynamic_cast<Thyra::SpmdVectorBase<Scalar> >(initial_value);
  vec->getNonconstLocalData(ptrFromRef(data));
  for (unsigned int i=0; i < paramObj->scalar_value.size(); i++)
    data[i] = in_values[i];

  paramObj->initial_value = initial_value;

  return paramObj;
}

template <typename Scalar>
Teuchos::RCP<typename panzer::ModelEvaluator<Scalar>::ParameterObject>
panzer::ModelEvaluator<Scalar>::
createDistributedParameter(const std::string & key,
                           const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > & vs,
                           const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & initial,
                           const Teuchos::RCP<const GlobalIndexer> & ugi) const
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

template <typename Scalar>
void
panzer::ModelEvaluator<Scalar>::
setParameters(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs) const
{
  for(std::size_t i=0; i < parameters_.size(); i++) {

     // skip non-scalar parameters (for now)
     if(parameters_[i]->is_distributed)
       continue;

     // set parameter values for given parameter vector for all evaluation types
     Teuchos::RCP<const Thyra::VectorBase<Scalar> > p = inArgs.get_p(i);
     if (p != Teuchos::null) {
       for (unsigned int j=0; j < parameters_[i]->scalar_value.size(); j++) {
         parameters_[i]->scalar_value[j].family->setRealValueForAllTypes(Thyra::get_ele(*p,j));
       }
     }

  }
}

template <typename Scalar>
void
panzer::ModelEvaluator<Scalar>::
resetParameters() const
{
  for(std::size_t i=0; i < parameters_.size(); i++) {

     // skip non-scalar parameters (for now)
     if(parameters_[i]->is_distributed)
       continue;

     // Reset each parameter back to its nominal
     for (unsigned int j=0; j < parameters_[i]->scalar_value.size(); j++) {
       parameters_[i]->scalar_value[j].family->setRealValueForAllTypes(Thyra::get_ele(*(parameters_[i]->initial_value),j));
     }

  }
}

#endif // __Panzer_ModelEvaluator_impl_hpp__
