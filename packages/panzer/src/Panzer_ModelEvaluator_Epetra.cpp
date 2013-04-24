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

#include "Panzer_config.hpp"
#include "Panzer_ModelEvaluator_Epetra.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_SGEpetraLinearObjFactory.hpp"
#include "Panzer_SGEpetraLinearObjContainer.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_GlobalData.hpp"

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

#ifdef HAVE_STOKHOS
   #include "Stokhos_EpetraVectorOrthogPoly.hpp"
   #include "Stokhos_EpetraOperatorOrthogPoly.hpp"
#endif

#include <sstream>

namespace {

//! For simple epetra operators
class EpetraLOF_EOpFactory : public Teuchos::AbstractFactory<Epetra_Operator> {
   Teuchos::RCP<Epetra_CrsGraph>  W_graph_;
public:
    EpetraLOF_EOpFactory(const panzer::EpetraLinearObjFactory<panzer::Traits,int> & lof)
       : W_graph_(lof.getGraph()) {}
    
    virtual Teuchos::RCP<Epetra_Operator> create() const
    { return Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_)); }
};

}


panzer::ModelEvaluator_Epetra::
ModelEvaluator_Epetra(const Teuchos::RCP<panzer::FieldManagerBuilder>& fmb,
                      const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> >& rLibrary,
		      const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & lof,
		      const std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >& p_names,
		      const Teuchos::RCP<panzer::GlobalData>& global_data,
		      bool build_transient_support)
  : t_init_(0.0)
  , fmb_(fmb)
  , responseLibrary_(rLibrary)
  , p_names_(p_names)
  , global_data_(global_data)
  , build_transient_support_(build_transient_support)
  , lof_(lof)
  #ifdef HAVE_STOKHOS
  , sg_lof_(Teuchos::null)
  #endif
{
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  panzer::AssemblyEngine_TemplateBuilder builder(fmb,lof);
  ae_tm_.buildObjects(builder);

  // Setup parameters
  this->initializeParameterVector(p_names_,global_data->pl);

  // try to determine the runtime linear object factory

  Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > ep_lof =
     Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjFactory<panzer::Traits,int> >(lof);

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
		      const Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> >& lof,
		      const std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >& p_names,
		      const Teuchos::RCP<panzer::GlobalData>& global_data,
		      bool build_transient_support)
  : t_init_(0.0)
  , fmb_(fmb)
  , responseLibrary_(rLibrary)
  , p_names_(p_names)
  , global_data_(global_data)
  , build_transient_support_(build_transient_support)
  , lof_(lof)
  #ifdef HAVE_STOKHOS
  , sg_lof_(Teuchos::null)
  #endif
{
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  panzer::AssemblyEngine_TemplateBuilder builder(fmb,lof);
  ae_tm_.buildObjects(builder);

  // Setup parameters
  this->initializeParameterVector(p_names_,global_data->pl);

  // initailize maps, x_dot_init, x0, p_init, g_map, and linear operator factory
  initializeEpetraObjs(*lof);
}

#ifdef HAVE_STOKHOS
panzer::ModelEvaluator_Epetra::
ModelEvaluator_Epetra(const Teuchos::RCP<panzer::FieldManagerBuilder>& fmb,
                      const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> >& rLibrary,
		      const Teuchos::RCP<panzer::SGEpetraLinearObjFactory<panzer::Traits,int> >& lof,
		      const std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >& p_names,
		      const Teuchos::RCP<panzer::GlobalData>& global_data,
		      bool build_transient_support) : 
  fmb_(fmb),
  responseLibrary_(rLibrary),
  p_names_(p_names),
  global_data_(global_data),
  build_transient_support_(build_transient_support),
  lof_(lof->getEpetraFactory()),
  sg_lof_(lof)
{
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  panzer::AssemblyEngine_TemplateBuilder builder(fmb,lof);
  ae_tm_.buildObjects(builder);

  // Setup parameters
  this->initializeParameterVector(p_names_,global_data->pl);

  // initailize maps, x_dot_init, x0, p_init, g_map, and linear operator factory
  initializeEpetraObjs(*lof->getEpetraFactory());
}
#endif

void panzer::ModelEvaluator_Epetra::initializeEpetraObjs(panzer::EpetraLinearObjFactory<panzer::Traits,int> & lof)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
 
  TEUCHOS_TEST_FOR_EXCEPTION(responseLibrary_==Teuchos::null,std::logic_error,
                     "panzer::ModelEvaluator_Epetra::initializeEpetraObjs: The response library "
                     "was not correctly initialized before calling initializeEpetraObjs.");

  map_x_ = lof.getMap();
  x0_ = rcp(new Epetra_Vector(*map_x_));
  x_dot_init_ = rcp(new Epetra_Vector(*map_x_));
  x_dot_init_->PutScalar(0.0);

  // setup parameters
  for (std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >::const_iterator p = p_names_.begin(); 
       p != p_names_.end(); ++p) {
    RCP<Epetra_Map> local_map = rcp(new Epetra_LocalMap(static_cast<int>((*p)->size()), 0, map_x_->Comm())) ;
    p_map_.push_back(local_map);
    RCP<Epetra_Vector> ep_vec = rcp(new Epetra_Vector(*local_map));
    ep_vec->PutScalar(0.0);
    p_init_.push_back(ep_vec);
  }
  
  // Initialize the epetra operator factory
  epetraOperatorFactory_ = Teuchos::rcp(new EpetraLOF_EOpFactory(lof));
}

void panzer::ModelEvaluator_Epetra::
initializeParameterVector(const std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >& p_names,
			  const Teuchos::RCP<panzer::ParamLib>& parameter_library)
{
  parameter_vector_.resize(p_names.size());
  is_distributed_parameter_.resize(p_names.size(),false);
  for (std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >::size_type p = 0; 
       p < p_names.size(); ++p) {
    parameter_library->fillVector<panzer::Traits::Residual>(*(p_names[p]), parameter_vector_[p]);
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

  distributed_parameter_container_.push_back(boost::make_tuple(name,index,importer,ghosted_vector));

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

#ifdef HAVE_STOKHOS
  if(!Teuchos::is_null(sg_lof_)) {
     inArgs.setSupports(IN_ARG_x_sg,true);
     // inArgs.setSupports(IN_ARG_x_dot_sg,true); NOT YET!
    
     // no parameter support yet!
     for (std::size_t i=0; i<p_map_.size(); i++)
         inArgs.setSupports(IN_ARG_p_sg, i, true);
   
     inArgs.setSupports(IN_ARG_sg_basis,true);
     inArgs.setSupports(IN_ARG_sg_quadrature,true);
     inArgs.setSupports(IN_ARG_sg_expansion,true);
  }
#endif

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
        outArgs.setSupports(OUT_ARG_DgDx,i,DerivativeSupport(DERIV_LINEAR_OP));
    }
  }

#ifdef HAVE_STOKHOS
  if(!Teuchos::is_null(sg_lof_)) {
     outArgs.setSupports(OUT_ARG_f_sg,true);
     outArgs.setSupports(OUT_ARG_W_sg,true);

     for(std::size_t i=0;i<g_map_.size();i++)
        outArgs.setSupports(OUT_ARG_g_sg,i,true);
  }
#endif

  return outArgs;
}

void panzer::ModelEvaluator_Epetra::evalModel( const InArgs& inArgs, 
					       const OutArgs& outArgs ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  evalModel_basic(inArgs,outArgs); 

  #ifdef HAVE_STOKHOS
    // use x_sg, p_sg to evaluate f_sg, W_sg and associated responses
    if(!Teuchos::is_null(sg_lof_))
       evalModel_sg(inArgs,outArgs);  
  #endif
}

void panzer::ModelEvaluator_Epetra::evalModel_basic( const InArgs& inArgs, 
 					             const OutArgs& outArgs ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  
  // Transient or steady-state evaluation is determined by the x_dot
  // vector.  If this RCP is null, then we are doing a steady-state
  // fill.
  bool is_transient = false;
  if (inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x_dot ))
    is_transient = nonnull(inArgs.get_x_dot());

  // Make sure construction built in transient support
  TEUCHOS_TEST_FOR_EXCEPTION(is_transient && !build_transient_support_, std::runtime_error,
		     "ModelEvaluator was not built with transient support enabled!");

  //
  // Get the output arguments
  //
  const RCP<Epetra_Vector> f_out = outArgs.get_f();
  const RCP<Epetra_Operator> W_out = outArgs.get_W();
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
  const RCP<const Epetra_Vector> x = inArgs.get_x();
  RCP<const Epetra_Vector> x_dot;
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
  
  // Set locally replicated scalar input parameters
  for (int i=0; i<inArgs.Np(); i++) {
    Teuchos::RCP<const Epetra_Vector> p = inArgs.get_p(i);
    if ( nonnull(p) && !is_distributed_parameter_[i]) {
      for (unsigned int j=0; j < parameter_vector_[i].size(); j++)
	parameter_vector_[i][j].baseValue = (*p)[j];
    }
  }

  for (Teuchos::Array<panzer::ParamVec>::size_type i=0; i < parameter_vector_.size(); i++)
    for (unsigned int j=0; j < parameter_vector_[i].size(); j++)
      parameter_vector_[i][j].family->setRealValueForAllTypes(parameter_vector_[i][j].baseValue);

  // Perform global to ghost and set distributed parameters
  for (std::vector<boost::tuple<std::string,int,Teuchos::RCP<Epetra_Import>,Teuchos::RCP<Epetra_Vector> > >::const_iterator i = 
	 distributed_parameter_container_.begin(); i != distributed_parameter_container_.end(); ++i) {
    // do export if parameter exists in inArgs
    Teuchos::RCP<const Epetra_Vector> global_vec = inArgs.get_p(i->get<1>());
    if (nonnull(global_vec)) {
      // Only import if the importer is nonnull
      Teuchos::RCP<Epetra_Import> importer = i->get<2>();
      if (nonnull(importer))
	i->get<3>()->Import(*global_vec,*importer,Insert);
    }

    // set in ae_inargs_ string lookup container
    Teuchos::RCP<panzer::EpetraLinearObjContainer> container = 
      Teuchos::rcp(new panzer::EpetraLinearObjContainer(p_map_[i->get<1>()],p_map_[i->get<1>()]));
    container->set_x(i->get<3>());
    ae_inargs.addGlobalEvaluationData(i->get<0>(),container);
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
  if (is_transient)
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

  // evaluate responses...uses the stored assembly arguments and containers
  if(requiredResponses)
     evalModel_basic_g(ae_inargs,inArgs,outArgs);
  
  // Holding a rcp to f produces a seg fault in Rythmos when the next
  // f comes in and the resulting dtor is called.  Need to discuss
  // with Ross.  Clearing all references here works!

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
evalModel_basic_g(AssemblyEngineInArgs ae_inargs,const InArgs & inArgs,const OutArgs & outArgs) const
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

   // evaluate response derivatives 
   if(required_basic_dgdx(outArgs))
     evalModel_basic_dgdx(ae_inargs,inArgs,outArgs);
}

void 
panzer::ModelEvaluator_Epetra::
evalModel_basic_dgdx(AssemblyEngineInArgs ae_inargs,const InArgs & inArgs,const OutArgs & outArgs) const
{
   // optional sanity check
   TEUCHOS_ASSERT(required_basic_dgdx(outArgs));

   for(std::size_t i=0;i<g_names_.size();i++) {
      // get "Vector" out of derivative, if its something else, throw an exception
      EpetraExt::ModelEvaluator::Derivative deriv = outArgs.get_DgDx(i);
      if(deriv.isEmpty())
        continue;

      Teuchos::RCP<Epetra_Vector> vec = Teuchos::rcp_dynamic_cast<Epetra_Vector>(deriv.getMultiVector(),true);

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

void panzer::ModelEvaluator_Epetra::set_t_init(double t)
{
  t_init_ = t;
}

#ifdef HAVE_STOKHOS
void 
panzer::ModelEvaluator_Epetra::
evalModel_sg(const InArgs & inArgs,const OutArgs & outArgs) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  TEUCHOS_ASSERT(!Teuchos::is_null(sg_lof_));
  
  // NO TRANSIENTS YET!!!!

  //
  // Get the output arguments
  //
  const RCP<Stokhos::EpetraVectorOrthogPoly> f_out = outArgs.get_f_sg();
  const RCP<Stokhos::EpetraOperatorOrthogPoly > W_out = outArgs.get_W_sg();
  bool requiredResponses = required_sg_g(outArgs);

  // see if the user wants us to do anything
  if(Teuchos::is_null(f_out) && Teuchos::is_null(W_out) && !requiredResponses)
     return;

  // the user requested work from this method
  // keep on moving

  // if neccessary build a ghosted container
  if(Teuchos::is_null(sg_ghostedContainer_)) {
     sg_ghostedContainer_ = sg_lof_->buildGhostedLinearObjContainer();
     sg_lof_->initializeGhostedContainer(panzer::LinearObjContainer::X |
                                         panzer::LinearObjContainer::DxDt |
                                         panzer::LinearObjContainer::F |
                                         panzer::LinearObjContainer::Mat, *sg_ghostedContainer_); 
  }

  //
  // Get the input arguments
  //
  const RCP<const Stokhos::EpetraVectorOrthogPoly> x_in = inArgs.get_x_sg();
  panzer::AssemblyEngineInArgs ae_inargs;
  ae_inargs.container_ = sg_lof_->buildLinearObjContainer(); // we use a new global container
  ae_inargs.ghostedContainer_ = sg_ghostedContainer_;        // we can reuse the ghosted container
  ae_inargs.alpha = 0.0;
  ae_inargs.beta = 1.0;
  ae_inargs.evaluate_transient_terms = false;
  
  // here we are building a container, this operation is fast, simply allocating a struct
  const RCP<panzer::SGEpetraLinearObjContainer> sgGlobalContainer = 
    Teuchos::rcp_dynamic_cast<panzer::SGEpetraLinearObjContainer>(ae_inargs.container_);

  // Ghosted container objects are zeroed out below only if needed for
  // a particular calculation.  This makes it more efficient thatn
  // zeroing out all objects in the container here.
  const RCP<panzer::SGEpetraLinearObjContainer> sgGhostedContainer = 
    Teuchos::rcp_dynamic_cast<panzer::SGEpetraLinearObjContainer>(ae_inargs.ghostedContainer_);

  TEUCHOS_ASSERT(!Teuchos::is_null(sgGlobalContainer));
  TEUCHOS_ASSERT(!Teuchos::is_null(sgGhostedContainer));

  
  // Set the solution vector (currently all targets require solution).
  // In the future we may move these into the individual cases below.
  // A very subtle (and fragile) point: A non-null pointer in global
  // container triggers export operations during fill.  Also, the
  // introduction of the container is forcing us to cast away const on
  // arguments that should be const.  Another reason to redesign
  // LinearObjContainer layers.
  //
  // copy sg data structure into linear object container data structure
  {
     TEUCHOS_ASSERT(x_in->size()==(int) sgGlobalContainer->size()); 
     TEUCHOS_ASSERT(x_in->size()==(int) sgGhostedContainer->size()); 
     if(!Teuchos::is_null(W_out)) { TEUCHOS_ASSERT(x_in->size()==W_out->size()); }
     if(!Teuchos::is_null(f_out)) { TEUCHOS_ASSERT(x_in->size()==f_out->size()); }
     
     // loop over each coefficient, setting up in/out arguments for the lo containers
     SGEpetraLinearObjContainer::iterator glbItr = sgGlobalContainer->begin();
     SGEpetraLinearObjContainer::iterator ghsItr = sgGhostedContainer->begin();
     for(int coeff_ind=0;coeff_ind<x_in->size();++coeff_ind,++glbItr,++ghsItr) {
        RCP<EpetraLinearObjContainer> globalContainer = *glbItr;
        RCP<EpetraLinearObjContainer> ghostedContainer = *ghsItr;

        // this is what roger was warning us about!!!!
        globalContainer->set_x(Teuchos::rcp_const_cast<Epetra_Vector>(x_in->getCoeffPtr(coeff_ind)));

        if(!Teuchos::is_null(f_out) && !Teuchos::is_null(W_out)) { // requires residual and jacobian
           globalContainer->set_f(f_out->getCoeffPtr(coeff_ind)); 
           globalContainer->set_A(rcp_dynamic_cast<Epetra_CrsMatrix>(W_out->getCoeffPtr(coeff_ind))); 
 
           ghostedContainer->get_f()->PutScalar(0.0);
           ghostedContainer->get_A()->PutScalar(0.0);
        }
        else if(!Teuchos::is_null(f_out) && Teuchos::is_null(W_out)) {
           globalContainer->set_f(f_out->getCoeffPtr(coeff_ind)); 
 
           // Zero values in ghosted container objects
           ghostedContainer->get_f()->PutScalar(0.0);
        }
        else if(Teuchos::is_null(f_out) && !Teuchos::is_null(W_out)) {

           // this dummy nonsense is needed only for scattering dirichlet conditions
           if(Teuchos::is_null(dummy_f_))
              dummy_f_ = Teuchos::rcp(new Epetra_Vector(*(this->get_f_map())));
           globalContainer->set_f(dummy_f_); 
           globalContainer->set_A(rcp_dynamic_cast<Epetra_CrsMatrix>(W_out->getCoeffPtr(coeff_ind))); 

           // Zero values in ghosted container objects
           ghostedContainer->get_A()->PutScalar(0.0);
        }
     }
  }
  
  if (!Teuchos::is_null(f_out) && !Teuchos::is_null(W_out)) {
    ae_tm_.getAsObject<panzer::Traits::SGJacobian>()->evaluate(ae_inargs);
  }
  else if(!Teuchos::is_null(f_out) && Teuchos::is_null(W_out)) {

    ae_tm_.getAsObject<panzer::Traits::SGResidual>()->evaluate(ae_inargs);

    // f_out->Update(1.0, *(epGlobalContainer->f), 0.0); // WHAT????
  }
  else if(Teuchos::is_null(f_out) && !Teuchos::is_null(W_out)) {
    ae_tm_.getAsObject<panzer::Traits::SGJacobian>()->evaluate(ae_inargs);
  }

  // evaluate responses...uses the stored assembly arguments and containers
  if(requiredResponses)
     evalModel_sg_g(ae_inargs,inArgs,outArgs);

  // forget previous containers
  ae_inargs.container_ = Teuchos::null;
  ae_inargs.ghostedContainer_ = Teuchos::null;
}

bool panzer::ModelEvaluator_Epetra::required_sg_g(const OutArgs & outArgs) const
{
   // determine if any of the outArgs are not null!
   bool activeGArgs = false;
   for(int i=0;i<outArgs.Ng();i++) 
      activeGArgs |= (outArgs.get_g_sg(i)!=Teuchos::null); 

   return activeGArgs;
}

void 
panzer::ModelEvaluator_Epetra::
evalModel_sg_g(AssemblyEngineInArgs ae_inargs,const InArgs & inArgs,const OutArgs & outArgs) const
{
   // build a teuchos comm from an mpi comm
   Teuchos::RCP<Teuchos::Comm<int> > tComm 
      = Teuchos::rcp(new Teuchos::MpiComm<int>(
        Teuchos::opaqueWrapper(dynamic_cast<const Epetra_MpiComm &>(map_x_->Comm()).Comm())));

   // evaluator responses
   responseLibrary_->evaluateVolumeFieldManagers<panzer::Traits::SGResidual>(ae_inargs,*tComm);

   std::vector<Teuchos::RCP<const Response<panzer::Traits> > > responses;
   responseLibrary_->getLabeledVolumeResponses(responses);
   for(std::size_t i=0;i<responses.size();i++) {

      // get destination vector
      Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_vec = outArgs.get_g_sg(i);
      if(sg_vec!=Teuchos::null) {
         // get stochastic value distribute them in vectors
         panzer::Traits::SGType value = responses[i]->getSGValue();
         for(int j=0;j<value.size();j++) {
            Epetra_Vector & vec = (*sg_vec)[j];

            // pull value out of stochastic type
            vec[0] = value.fastAccessCoeff(j);
            for(int l=1;l<vec.MyLength();l++)
               vec[i] = 0.0;
         }

         // zero out uninitialized values: Stokhos didn't require these fields
         for(int j=value.size();j<sg_vec->size();j++) {
            Epetra_Vector & vec = (*sg_vec)[j];
            vec.PutScalar(0.0);
         }
      }
   }
}

#endif

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
		      const Teuchos::RCP<panzer::GlobalData>& global_data,
	              bool build_transient_support)
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;
   using Teuchos::rcp;

   std::stringstream ss;
   ss << "panzer::buildEpetraME: Linear object factory is incorrect type. Should be one of: ";

   RCP<EpetraLinearObjFactory<panzer::Traits,int> > ep_lof 
       = rcp_dynamic_cast<EpetraLinearObjFactory<panzer::Traits,int> >(lof);

   // if you can, build from an epetra linear object factory
   if(ep_lof!=Teuchos::null) 
     return rcp(new ModelEvaluator_Epetra(fmb,rLibrary,ep_lof,p_names,global_data,build_transient_support));

   ss << "\"EpetraLinearObjFactory\", ";

#ifdef HAVE_STOKHOS
   RCP<SGEpetraLinearObjFactory<panzer::Traits,int> > sg_lof 
       = rcp_dynamic_cast<SGEpetraLinearObjFactory<panzer::Traits,int> >(lof);

   // if you can, build from an SG epetra linear object factory
   if(sg_lof!=Teuchos::null) 
     return rcp(new ModelEvaluator_Epetra(fmb,rLibrary,sg_lof,p_names,global_data,build_transient_support));

   ss << "\"SGEpetraLinearObjFactory\", ";
#endif

   TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,ss.str());
}
