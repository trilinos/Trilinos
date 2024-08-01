// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "PanzerDiscFE_config.hpp"
#include "Panzer_ModelEvaluator_Epetra.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_ParameterList_GlobalEvaluationData.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"

#include "Epetra_CrsMatrix.h"
#include "Epetra_MpiComm.h"
#include "Epetra_LocalMap.h"

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Thyra_EpetraOperatorWrapper.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"

#include <sstream>

namespace {

//! For simple epetra operators
class EpetraLOF_EOpFactory : public Teuchos::AbstractFactory<Epetra_Operator> {
   Teuchos::RCP<Epetra_CrsGraph>  W_graph_;
public:
    EpetraLOF_EOpFactory(const panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> & lof)
       : W_graph_(lof.getGraph(0,0)) {}

    virtual Teuchos::RCP<Epetra_Operator> create() const
    { return Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_)); }
};

}


panzer::ModelEvaluator_Epetra::
ModelEvaluator_Epetra(const Teuchos::RCP<panzer::FieldManagerBuilder>& fmb,
                      const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> >& rLibrary,
                      const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & lof,
                      const std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >& p_names,
                      const std::vector<Teuchos::RCP<Teuchos::Array<double> > >& p_values,
                      const Teuchos::RCP<panzer::GlobalData>& global_data,
                      bool build_transient_support)
  : t_init_(0.0)
  , fmb_(fmb)
  , responseLibrary_(rLibrary)
  , p_names_(p_names)
  , global_data_(global_data)
  , build_transient_support_(build_transient_support)
  , lof_(lof)
  , oneTimeDirichletBeta_on_(false)
  , oneTimeDirichletBeta_(0.0)
{
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  panzer::AssemblyEngine_TemplateBuilder builder(fmb,lof);
  ae_tm_.buildObjects(builder);

  // Setup parameters
  this->initializeParameterVector(p_names_,p_values,global_data->pl);

  // try to determine the runtime linear object factory

  Teuchos::RCP<panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> > ep_lof =
     Teuchos::rcp_dynamic_cast<panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> >(lof);

  // initialize maps, x_dot_init, x0, p_init, g_map, and linear operator factory
  if(ep_lof!=Teuchos::null)
     initializeEpetraObjs(*ep_lof);
  else {
     TEUCHOS_ASSERT(false); // bad news!
  }
}

panzer::ModelEvaluator_Epetra::
ModelEvaluator_Epetra(const Teuchos::RCP<panzer::FieldManagerBuilder>& fmb,
                      const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> >& rLibrary,
                      const Teuchos::RCP<panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> >& lof,
                      const std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >& p_names,
                      const std::vector<Teuchos::RCP<Teuchos::Array<double> > >& p_values,
                      const Teuchos::RCP<panzer::GlobalData>& global_data,
                      bool build_transient_support)
  : t_init_(0.0)
  , fmb_(fmb)
  , responseLibrary_(rLibrary)
  , p_names_(p_names)
  , global_data_(global_data)
  , build_transient_support_(build_transient_support)
  , lof_(lof)
  , oneTimeDirichletBeta_on_(false)
  , oneTimeDirichletBeta_(0.0)
{
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  panzer::AssemblyEngine_TemplateBuilder builder(fmb,lof);
  ae_tm_.buildObjects(builder);

  // Setup parameters
  this->initializeParameterVector(p_names_,p_values,global_data->pl);

  // initailize maps, x_dot_init, x0, p_init, g_map, and linear operator factory
  initializeEpetraObjs(*lof);
}

void panzer::ModelEvaluator_Epetra::initializeEpetraObjs(panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> & lof)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  TEUCHOS_TEST_FOR_EXCEPTION(responseLibrary_==Teuchos::null,std::logic_error,
                     "panzer::ModelEvaluator_Epetra::initializeEpetraObjs: The response library "
                     "was not correctly initialized before calling initializeEpetraObjs.");

  map_x_ = lof.getMap(0);
  x0_ = rcp(new Epetra_Vector(*map_x_));
  x_dot_init_ = rcp(new Epetra_Vector(*map_x_));
  x_dot_init_->PutScalar(0.0);

  // setup parameters (initialize vector from parameter library)
  for (int i=0; i < parameter_vector_.size(); i++) {
    RCP<Epetra_Map> local_map = rcp(new Epetra_LocalMap(static_cast<int>(parameter_vector_[i].size()), 0, map_x_->Comm())) ;
    p_map_.push_back(local_map);

    RCP<Epetra_Vector> ep_vec = rcp(new Epetra_Vector(*local_map));
    ep_vec->PutScalar(0.0);

    for (unsigned int j=0; j < parameter_vector_[i].size(); j++)
      (*ep_vec)[j] = parameter_vector_[i][j].baseValue;

    p_init_.push_back(ep_vec);
  }

  // Initialize the epetra operator factory
  epetraOperatorFactory_ = Teuchos::rcp(new EpetraLOF_EOpFactory(lof));
}

void panzer::ModelEvaluator_Epetra::
initializeParameterVector(const std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >& p_names,
                          const std::vector<Teuchos::RCP<Teuchos::Array<double> > >& p_values,
                          const Teuchos::RCP<panzer::ParamLib>& parameter_library)
{
  parameter_vector_.resize(p_names.size());
  is_distributed_parameter_.resize(p_names.size(),false);
  for (std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >::size_type p = 0;
       p < p_names.size(); ++p) {

    // register all the scalar parameters, setting initial
    for(int i=0;i<p_names[p]->size();i++)
      registerScalarParameter((*p_names[p])[i],*parameter_library,(*p_values[p])[i]);

    parameter_library->fillVector<panzer::Traits::Residual>(*(p_names[p]), parameter_vector_[p]);

    for(int i=0;i<p_names[p]->size();i++) {
      parameter_vector_[p][i].baseValue = (*p_values[p])[i];
      parameter_vector_[p][i].family->setRealValueForAllTypes((*p_values[p])[i]);
    }
  }
}

// Post-Construction methods to add parameters and responses

int panzer::ModelEvaluator_Epetra::
addDistributedParameter(const std::string name,
                        const Teuchos::RCP<Epetra_Map>& global_map,
                        const Teuchos::RCP<Epetra_Import>& importer,
                        const Teuchos::RCP<Epetra_Vector>& ghosted_vector)
{
  // Will push_back a new parameter entry
  int index = static_cast<int>(p_map_.size());

  p_map_.push_back(global_map);
  Teuchos::RCP<Epetra_Vector> ep_vec = Teuchos::rcp(new Epetra_Vector(*global_map));
  ep_vec->PutScalar(0.0);
  p_init_.push_back(ep_vec);

  Teuchos::RCP<Teuchos::Array<std::string> > tmp_names =
    Teuchos::rcp(new Teuchos::Array<std::string>);
  tmp_names->push_back(name);
  p_names_.push_back(tmp_names);

  is_distributed_parameter_.push_back(true);

  // NOTE: we do not add this parameter to the sacado parameter
  // library in the global data object.  That library is for scalars.
  // We will need to do something different if we need sensitivities
  // wrt distributed parameters.

  distributed_parameter_container_.push_back(std::make_tuple(name,index,importer,ghosted_vector));

  return index;
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
panzer::ModelEvaluator_Epetra::get_x_map() const
{
  return map_x_;
}

Teuchos::RCP<const Epetra_Map>
panzer::ModelEvaluator_Epetra::get_f_map() const
{
  return map_x_;
}

Teuchos::RCP<const Epetra_Vector>
panzer::ModelEvaluator_Epetra::get_x_init() const
{
  return x0_;
}

Teuchos::RCP<const Epetra_Vector>
panzer::ModelEvaluator_Epetra::get_x_dot_init() const
{
  return x_dot_init_;
}

double
panzer::ModelEvaluator_Epetra::get_t_init() const
{
  return t_init_;
}

Teuchos::RCP<Epetra_Operator>
panzer::ModelEvaluator_Epetra::create_W() const
{
  return epetraOperatorFactory_->create();
}

Teuchos::RCP<const Epetra_Map>
panzer::ModelEvaluator_Epetra::get_p_map(int l) const
{
  return p_map_[l];
}

Teuchos::RCP<const Teuchos::Array<std::string> >
panzer::ModelEvaluator_Epetra::get_p_names(int l) const
{
  return p_names_[l];
}

Teuchos::RCP<const Epetra_Vector>
panzer::ModelEvaluator_Epetra::get_p_init(int l) const
{
  return p_init_[l];
}

Teuchos::RCP<const Epetra_Map>
panzer::ModelEvaluator_Epetra::get_g_map(int l) const
{
  return g_map_[l];
}

EpetraExt::ModelEvaluator::InArgs
panzer::ModelEvaluator_Epetra::createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(IN_ARG_x,true);
  if (build_transient_support_) {
    inArgs.setSupports(IN_ARG_x_dot,true);
    inArgs.setSupports(IN_ARG_t,true);
    inArgs.setSupports(IN_ARG_alpha,true);
    inArgs.setSupports(IN_ARG_beta,true);
  }
  inArgs.set_Np(p_init_.size());

  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
panzer::ModelEvaluator_Epetra::createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(p_init_.size(), g_map_.size());
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
  outArgs.set_W_properties(
    DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_FULL
      ,true // supportsAdjoint
      )
    );

  // add in df/dp (if appropriate)
  for(std::size_t i=0;i<p_init_.size();i++) {
    if(!is_distributed_parameter_[i])
      outArgs.setSupports(OUT_ARG_DfDp,i,EpetraExt::ModelEvaluator::DerivativeSupport(DERIV_MV_BY_COL));
  }

  // add in dg/dx (if appropriate)
  for(std::size_t i=0;i<g_names_.size();i++) {
    typedef panzer::Traits::Jacobian RespEvalT;

    // check dg/dx and add it in if appropriate
    Teuchos::RCP<panzer::ResponseBase> respJacBase = responseLibrary_->getResponse<RespEvalT>(g_names_[i]);
    if(respJacBase!=Teuchos::null) {
      // cast is guranteed to succeed because of check in addResponse
      Teuchos::RCP<panzer::ResponseMESupportBase<RespEvalT> > resp = Teuchos::rcp_dynamic_cast<panzer::ResponseMESupportBase<RespEvalT> >(respJacBase);

      // class must supoprt a derivative
      if(resp->supportsDerivative())
        outArgs.setSupports(OUT_ARG_DgDx,i,DerivativeSupport(DERIV_TRANS_MV_BY_ROW));
        //outArgs.setSupports(OUT_ARG_DgDx,i,DerivativeSupport(DERIV_LINEAR_OP));
    }
  }

  return outArgs;
}

void panzer::ModelEvaluator_Epetra::
applyDirichletBCs(const Teuchos::RCP<Thyra::VectorBase<double> > & x,
                  const Teuchos::RCP<Thyra::VectorBase<double> > & f) const
{
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::Array;
  //using Teuchos::tuple;
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

  // this is the tempory target
  lof_->initializeContainer(panzer::LinearObjContainer::F,*ae_inargs.container_);

  // here we are building a container, this operation is fast, simply allocating a struct
  const RCP<panzer::ThyraObjContainer<double> > thGlobalContainer =
    Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(ae_inargs.container_,true);

  TEUCHOS_ASSERT(!Teuchos::is_null(thGlobalContainer));

  // Ghosted container objects are zeroed out below only if needed for
  // a particular calculation.  This makes it more efficient than
  // zeroing out all objects in the container here.
  const RCP<panzer::ThyraObjContainer<double> > thGhostedContainer =
    Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(ae_inargs.ghostedContainer_,true);
  TEUCHOS_ASSERT(!Teuchos::is_null(thGhostedContainer));
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
  RCP<panzer::LinearObjContainer> counter = ae_tm_.getAsObject<panzer::Traits::Residual>()->evaluateOnlyDirichletBCs(ae_inargs);

  // allocate the result container
  RCP<panzer::LinearObjContainer> result = lof_->buildLinearObjContainer(); // we use a new global container

  // stuff the evaluate boundary conditions into the f spot of the counter ... the x is already filled
  Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(counter)->set_f_th(thGlobalContainer->get_f_th());

  // stuff the vector that needs applied dirichlet conditions in the the f spot of the result LOC
  Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(result)->set_f_th(f);

  // use the linear object factory to apply the result
  lof_->applyDirichletBCs(*counter,*result);
}

void panzer::ModelEvaluator_Epetra::evalModel( const InArgs& inArgs,
                                               const OutArgs& outArgs ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  evalModel_basic(inArgs,outArgs);
}

void panzer::ModelEvaluator_Epetra::evalModel_basic( const InArgs& inArgs,
                                                      const OutArgs& outArgs ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  // Transient or steady-state evaluation is determined by the x_dot
  // vector.  If this RCP is null, then we are doing a steady-state
  // fill.
  bool has_x_dot = false;
  if (inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x_dot ))
    has_x_dot = nonnull(inArgs.get_x_dot());

  // Make sure construction built in transient support
  TEUCHOS_TEST_FOR_EXCEPTION(has_x_dot && !build_transient_support_, std::runtime_error,
                     "ModelEvaluator was not built with transient support enabled!");

  //
  // Get the output arguments
  //
  const RCP<Epetra_Vector> f_out = outArgs.get_f();
  const RCP<Epetra_Operator> W_out = outArgs.get_W();
  bool requiredResponses = required_basic_g(outArgs);
  bool requiredSensitivities = required_basic_dfdp(outArgs);

  // see if the user wants us to do anything
  if(Teuchos::is_null(f_out) && Teuchos::is_null(W_out) &&
     !requiredResponses && !requiredSensitivities) {
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
  const RCP<const Epetra_Vector> x = inArgs.get_x();
  RCP<const Epetra_Vector> x_dot;

  panzer::AssemblyEngineInArgs ae_inargs;
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

  // handle application of the one time dirichlet beta in the
  // assembly engine. Note that this has to be set explicitly
  // each time because this badly breaks encapsulation. Essentially
  // we must work around the model evaluator abstraction!
  ae_inargs.apply_dirichlet_beta = false;
  if(oneTimeDirichletBeta_on_) {
    ae_inargs.dirichlet_beta = oneTimeDirichletBeta_;
    ae_inargs.apply_dirichlet_beta = true;

    oneTimeDirichletBeta_on_ = false;
  }

  // Set locally replicated scalar input parameters
  for (int i=0; i<inArgs.Np(); i++) {
    Teuchos::RCP<const Epetra_Vector> p = inArgs.get_p(i);
    if ( nonnull(p) && !is_distributed_parameter_[i]) {
      for (unsigned int j=0; j < parameter_vector_[i].size(); j++) {
        parameter_vector_[i][j].baseValue = (*p)[j];
      }
    }
  }

  for (Teuchos::Array<panzer::ParamVec>::size_type i=0; i < parameter_vector_.size(); i++) {
    for (unsigned int j=0; j < parameter_vector_[i].size(); j++) {
      parameter_vector_[i][j].family->setRealValueForAllTypes(parameter_vector_[i][j].baseValue);
    }
  }

  // Perform global to ghost and set distributed parameters
  for (std::vector<std::tuple<std::string,int,Teuchos::RCP<Epetra_Import>,Teuchos::RCP<Epetra_Vector> > >::const_iterator i =
         distributed_parameter_container_.begin(); i != distributed_parameter_container_.end(); ++i) {
    // do export if parameter exists in inArgs
    Teuchos::RCP<const Epetra_Vector> global_vec = inArgs.get_p(std::get<1>(*i));
    if (nonnull(global_vec)) {
      // Only import if the importer is nonnull
      Teuchos::RCP<Epetra_Import> importer = std::get<2>(*i);
      if (nonnull(importer))
	std::get<3>(*i)->Import(*global_vec,*importer,Insert);
    }

    // set in ae_inargs_ string lookup container
    Teuchos::RCP<panzer::EpetraLinearObjContainer> container =
      Teuchos::rcp(new panzer::EpetraLinearObjContainer(p_map_[std::get<1>(*i)],p_map_[std::get<1>(*i)]));
    container->set_x(std::get<3>(*i));
    ae_inargs.addGlobalEvaluationData(std::get<0>(*i),container);
  }

  // here we are building a container, this operation is fast, simply allocating a struct
  const RCP<panzer::EpetraLinearObjContainer> epGlobalContainer =
    Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(ae_inargs.container_);

  TEUCHOS_ASSERT(!Teuchos::is_null(epGlobalContainer));

  // Ghosted container objects are zeroed out below only if needed for
  // a particular calculation.  This makes it more efficient thatn
  // zeroing out all objects in the container here.
  const RCP<panzer::EpetraLinearObjContainer> epGhostedContainer =
    Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(ae_inargs.ghostedContainer_);

  // Set the solution vector (currently all targets require solution).
  // In the future we may move these into the individual cases below.
  // A very subtle (and fragile) point: A non-null pointer in global
  // container triggers export operations during fill.  Also, the
  // introduction of the container is forcing us to cast away const on
  // arguments that should be const.  Another reason to redesign
  // LinearObjContainer layers.
  epGlobalContainer->set_x(Teuchos::rcp_const_cast<Epetra_Vector>(x));
  if (has_x_dot)
    epGlobalContainer->set_dxdt(Teuchos::rcp_const_cast<Epetra_Vector>(x_dot));

  if (!Teuchos::is_null(f_out) && !Teuchos::is_null(W_out)) {

    PANZER_FUNC_TIME_MONITOR("panzer::ModelEvaluator::evalModel(f and J)");

    // Set the targets
    epGlobalContainer->set_f(f_out);
    epGlobalContainer->set_A(Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_out));

    // Zero values in ghosted container objects
    epGhostedContainer->get_f()->PutScalar(0.0);
    epGhostedContainer->get_A()->PutScalar(0.0);

    ae_tm_.getAsObject<panzer::Traits::Jacobian>()->evaluate(ae_inargs);
  }
  else if(!Teuchos::is_null(f_out) && Teuchos::is_null(W_out)) {

    PANZER_FUNC_TIME_MONITOR("panzer::ModelEvaluator::evalModel(f)");

    epGlobalContainer->set_f(f_out);

    // Zero values in ghosted container objects
    epGhostedContainer->get_f()->PutScalar(0.0);

    ae_tm_.getAsObject<panzer::Traits::Residual>()->evaluate(ae_inargs);

    f_out->Update(1.0, *(epGlobalContainer->get_f()), 0.0);
  }
  else if(Teuchos::is_null(f_out) && !Teuchos::is_null(W_out)) {

    PANZER_FUNC_TIME_MONITOR("panzer::ModelEvaluator::evalModel(J)");

    // this dummy nonsense is needed only for scattering dirichlet conditions
    if (Teuchos::is_null(dummy_f_))
      dummy_f_ = Teuchos::rcp(new Epetra_Vector(*(this->get_f_map())));
    epGlobalContainer->set_f(dummy_f_);
    epGlobalContainer->set_A(Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_out));

    // Zero values in ghosted container objects
    epGhostedContainer->get_A()->PutScalar(0.0);

    ae_tm_.getAsObject<panzer::Traits::Jacobian>()->evaluate(ae_inargs);
  }
  // HACK: set A to null before calling responses to avoid touching the
  // the Jacobian after it has been properly assembled.  Should be fixed
  // by using a modified version of ae_inargs instead.
  epGlobalContainer->set_A(Teuchos::null);

  // evaluate responses...uses the stored assembly arguments and containers
  if(requiredResponses) {
     evalModel_basic_g(ae_inargs,inArgs,outArgs);

     // evaluate response derivatives
     if(required_basic_dgdx(outArgs))
       evalModel_basic_dgdx(ae_inargs,inArgs,outArgs);
  }

  if(required_basic_dfdp(outArgs))
     evalModel_basic_dfdp(ae_inargs,inArgs,outArgs);

  // TODO: Clearing all references prevented a seg-fault with Rythmos,
  // which is no longer used. Check if it's still needed.
  epGlobalContainer->set_x(Teuchos::null);
  epGlobalContainer->set_dxdt(Teuchos::null);
  epGlobalContainer->set_f(Teuchos::null);
  epGlobalContainer->set_A(Teuchos::null);

  // forget previous containers
  ae_inargs.container_ = Teuchos::null;
  ae_inargs.ghostedContainer_ = Teuchos::null;
}

void
panzer::ModelEvaluator_Epetra::
evalModel_basic_g(AssemblyEngineInArgs ae_inargs, const InArgs& /* inArgs */, const OutArgs& outArgs) const
{
   // optional sanity check
   // TEUCHOS_ASSERT(required_basic_g(outArgs));

   for(std::size_t i=0;i<g_names_.size();i++) {
      Teuchos::RCP<Epetra_Vector> vec = outArgs.get_g(i);
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

void
panzer::ModelEvaluator_Epetra::
evalModel_basic_dgdx(AssemblyEngineInArgs ae_inargs, const InArgs& /* inArgs */, const OutArgs& outArgs) const
{
   // optional sanity check
   TEUCHOS_ASSERT(required_basic_dgdx(outArgs));

   for(std::size_t i=0;i<g_names_.size();i++) {
      // get "Vector" out of derivative, if its something else, throw an exception
      EpetraExt::ModelEvaluator::Derivative deriv = outArgs.get_DgDx(i);
      if(deriv.isEmpty())
        continue;

      Teuchos::RCP<Epetra_MultiVector> vec = deriv.getMultiVector();

      if(vec!=Teuchos::null) {
        std::string responseName = g_names_[i];
        Teuchos::RCP<panzer::ResponseMESupportBase<panzer::Traits::Jacobian> > resp
            = Teuchos::rcp_dynamic_cast<panzer::ResponseMESupportBase<panzer::Traits::Jacobian> >(responseLibrary_->getResponse<panzer::Traits::Jacobian>(responseName));
        resp->setDerivative(vec);
      }
   }

   // evaluator responses
   responseLibrary_->addResponsesToInArgs<panzer::Traits::Jacobian>(ae_inargs);
   responseLibrary_->evaluate<panzer::Traits::Jacobian>(ae_inargs);
}

void
panzer::ModelEvaluator_Epetra::
evalModel_basic_dfdp(AssemblyEngineInArgs ae_inargs, const InArgs& /* inArgs */, const OutArgs& outArgs) const
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   TEUCHOS_ASSERT(required_basic_dfdp(outArgs));

   // dynamic cast to blocked LOF for now
   RCP<const Thyra::VectorSpaceBase<double> > glblVS = rcp_dynamic_cast<const ThyraObjFactory<double> >(lof_,true)->getThyraRangeSpace();;

   std::vector<std::string> activeParameters;

   // fill parameter vector containers
   int totalParameterCount = 0;
   for(Teuchos::Array<panzer::ParamVec>::size_type i=0; i < parameter_vector_.size(); i++) {
     // have derivatives been requested?
     EpetraExt::ModelEvaluator::Derivative deriv = outArgs.get_DfDp(i);
     if(deriv.isEmpty())
       continue;

     // grab multivector, make sure its the right dimension
     Teuchos::RCP<Epetra_MultiVector> mVec = deriv.getMultiVector();
     TEUCHOS_ASSERT(mVec->NumVectors()==int(parameter_vector_[i].size()));

     for (unsigned int j=0; j < parameter_vector_[i].size(); j++) {

       // build containers for each vector
       RCP<LOCPair_GlobalEvaluationData> loc_pair = Teuchos::rcp(new LOCPair_GlobalEvaluationData(lof_,LinearObjContainer::F));
       RCP<LinearObjContainer> globalContainer = loc_pair->getGlobalLOC();

       // stuff target vector into global container
       RCP<Epetra_Vector> vec = Teuchos::rcpFromRef(*(*mVec)(j));
       RCP<panzer::ThyraObjContainer<double> > thGlobalContainer =
         Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(globalContainer);
       thGlobalContainer->set_f_th(Thyra::create_Vector(vec,glblVS));

       // add container into in args object
       std::string name = "PARAMETER_SENSITIVIES: "+(*p_names_[i])[j];
       ae_inargs.addGlobalEvaluationData(name,loc_pair->getGhostedLOC());
       ae_inargs.addGlobalEvaluationData(name+"_pair",loc_pair);

       activeParameters.push_back(name);
       totalParameterCount++;
     }
   }

   // this communicates to the scatter evaluators so that the appropriate parameters are scattered
   RCP<GlobalEvaluationData> ged_activeParameters = Teuchos::rcp(new ParameterList_GlobalEvaluationData(activeParameters));
   ae_inargs.addGlobalEvaluationData("PARAMETER_NAMES",ged_activeParameters);

   int paramIndex = 0;
   for (Teuchos::Array<panzer::ParamVec>::size_type i=0; i < parameter_vector_.size(); i++) {
     // don't modify the parameter if its not needed
     EpetraExt::ModelEvaluator::Derivative deriv = outArgs.get_DfDp(i);
     if(deriv.isEmpty()) {
       // reinitialize values that should not have sensitivities computed (this is a precaution)
       for (unsigned int j=0; j < parameter_vector_[i].size(); j++) {
         Traits::FadType p = Traits::FadType(totalParameterCount, parameter_vector_[i][j].baseValue);
         parameter_vector_[i][j].family->setValue<panzer::Traits::Tangent>(p);
       }
       continue;
     }
     else {
       // loop over each parameter in the vector, initializing the AD type
       for (unsigned int j=0; j < parameter_vector_[i].size(); j++) {
         Traits::FadType p = Traits::FadType(totalParameterCount, parameter_vector_[i][j].baseValue);
         p.fastAccessDx(paramIndex) = 1.0;
         parameter_vector_[i][j].family->setValue<panzer::Traits::Tangent>(p);
         paramIndex++;
       }
     }
   }

   // make sure that the total parameter count and the total parameter index match up
   TEUCHOS_ASSERT(paramIndex==totalParameterCount);

   if(totalParameterCount>0) {
     ae_tm_.getAsObject<panzer::Traits::Tangent>()->evaluate(ae_inargs);
   }
}

bool panzer::ModelEvaluator_Epetra::required_basic_g(const OutArgs & outArgs) const
{
   // determine if any of the outArgs are not null!
   bool activeGArgs = false;
   for(int i=0;i<outArgs.Ng();i++)
      activeGArgs |= (outArgs.get_g(i)!=Teuchos::null);

   return activeGArgs | required_basic_dgdx(outArgs);
}

bool panzer::ModelEvaluator_Epetra::required_basic_dgdx(const OutArgs & outArgs) const
{
   // determine if any of the outArgs are not null!
   bool activeGArgs = false;
   for(int i=0;i<outArgs.Ng();i++) {
     // no derivatives are supported
     if(outArgs.supports(OUT_ARG_DgDx,i).none())
       continue;

     // this is basically a redundant computation
     activeGArgs |= (!outArgs.get_DgDx(i).isEmpty());
   }

   return activeGArgs;
}

bool panzer::ModelEvaluator_Epetra::required_basic_dfdp(const OutArgs & outArgs) const
{
   // determine if any of the outArgs are not null!
   bool activeFPArgs = false;
   for(int i=0;i<outArgs.Np();i++) {
     // no derivatives are supported
     if(outArgs.supports(OUT_ARG_DfDp,i).none())
       continue;

     // this is basically a redundant computation
     activeFPArgs |= (!outArgs.get_DfDp(i).isEmpty());
   }

   return activeFPArgs;
}

void panzer::ModelEvaluator_Epetra::set_t_init(double t)
{
  t_init_ = t;
}

void panzer::ModelEvaluator_Epetra::
copyThyraIntoEpetra(const Thyra::VectorBase<double>& thyraVec, Epetra_MultiVector& x) const
{

  using Teuchos::rcpFromRef;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::RCP;
  using Teuchos::ArrayView;
  using Teuchos::rcp_dynamic_cast;

  const int numVecs = x.NumVectors();

  TEUCHOS_TEST_FOR_EXCEPTION(numVecs != 1, std::runtime_error,
    "epetraToThyra does not work with MV dimension != 1");

  const RCP<const Thyra::ProductVectorBase<double> > prodThyraVec =
    Thyra::castOrCreateProductVectorBase(rcpFromRef(thyraVec));

  const Teuchos::ArrayView<double> epetraData(x[0], x.Map().NumMyElements());
  // NOTE: See above!

  std::size_t offset = 0;
  const int numBlocks = prodThyraVec->productSpace()->numBlocks();
  for (int b = 0; b < numBlocks; ++b) {
    const RCP<const Thyra::VectorBase<double> > vec_b = prodThyraVec->getVectorBlock(b);
    const RCP<const Thyra::SpmdVectorBase<double> > spmd_b =
           rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(vec_b, true);

    Teuchos::ArrayRCP<const double> thyraData;
    spmd_b->getLocalData(Teuchos::ptrFromRef(thyraData));

    for (Teuchos::ArrayRCP<const double>::size_type i=0; i < thyraData.size(); ++i) {
      epetraData[i+offset] = thyraData[i];
    }
    offset += thyraData.size();
  }

}


void panzer::ModelEvaluator_Epetra::
copyEpetraIntoThyra(const Epetra_MultiVector& x, const Teuchos::Ptr<Thyra::VectorBase<double> > &thyraVec) const
{
  using Teuchos::RCP;
  using Teuchos::ArrayView;
  using Teuchos::rcpFromPtr;
  using Teuchos::rcp_dynamic_cast;

  const int numVecs = x.NumVectors();

  TEUCHOS_TEST_FOR_EXCEPTION(numVecs != 1, std::runtime_error,
    "ModelEvaluator_Epetra::copyEpetraToThyra does not work with MV dimension != 1");

  const RCP<Thyra::ProductVectorBase<double> > prodThyraVec =
    Thyra::castOrCreateNonconstProductVectorBase(rcpFromPtr(thyraVec));

  const Teuchos::ArrayView<const double> epetraData(x[0], x.Map().NumMyElements());
  // NOTE: I tried using Epetra_MultiVector::operator()(int) to return an
  // Epetra_Vector object but it has a defect when Reset(...) is called which
  // results in a memory access error (see bug 4700).

  std::size_t offset = 0;
  const int numBlocks = prodThyraVec->productSpace()->numBlocks();
  for (int b = 0; b < numBlocks; ++b) {
    const RCP<Thyra::VectorBase<double> > vec_b = prodThyraVec->getNonconstVectorBlock(b);
    const RCP<Thyra::SpmdVectorBase<double> > spmd_b =
           rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(vec_b, true);

    Teuchos::ArrayRCP<double> thyraData;
    spmd_b->getNonconstLocalData(Teuchos::ptrFromRef(thyraData));

    for (Teuchos::ArrayRCP<double>::size_type i=0; i < thyraData.size(); ++i) {
      thyraData[i] = epetraData[i+offset];
    }
    offset += thyraData.size();
  }

}


Teuchos::RCP<panzer::ModelEvaluator_Epetra>
panzer::buildEpetraME(const Teuchos::RCP<panzer::FieldManagerBuilder>& fmb,
                      const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> >& rLibrary,
                      const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> >& lof,
                      const std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >& p_names,
	        const std::vector<Teuchos::RCP<Teuchos::Array<double> > >& p_values,
                      const Teuchos::RCP<panzer::GlobalData>& global_data,
                      bool build_transient_support)
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;
   using Teuchos::rcp;

   std::stringstream ss;
   ss << "panzer::buildEpetraME: Linear object factory is incorrect type. Should be one of: ";

   RCP<BlockedEpetraLinearObjFactory<panzer::Traits,int> > ep_lof
       = rcp_dynamic_cast<BlockedEpetraLinearObjFactory<panzer::Traits,int> >(lof);

   // if you can, build from an epetra linear object factory
   if(ep_lof!=Teuchos::null)
     return rcp(new ModelEvaluator_Epetra(fmb,rLibrary,ep_lof,p_names,p_values,global_data,build_transient_support));

   ss << "\"BlockedEpetraLinearObjFactory\", ";

   TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,ss.str());
}

void panzer::ModelEvaluator_Epetra::
setOneTimeDirichletBeta(const double & beta) const
{
  oneTimeDirichletBeta_on_ = true;
  oneTimeDirichletBeta_    = beta;
}