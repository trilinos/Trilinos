#include "Panzer_ModelEvaluator_Epetra.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_SGEpetraLinearObjFactory.hpp"
#include "Panzer_SGEpetraLinearObjContainer.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_ResponseLibrary.hpp"

#include "Epetra_CrsMatrix.h"
#include "Epetra_MpiComm.h"
#include "Epetra_LocalMap.h"

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_DefaultMpiComm.hpp"

#ifdef HAVE_STOKHOS
   #include "Stokhos_EpetraVectorOrthogPoly.hpp"
   #include "Stokhos_EpetraOperatorOrthogPoly.hpp"
#endif

panzer::ModelEvaluator_Epetra::
ModelEvaluator_Epetra(const Teuchos::RCP<panzer::FieldManagerBuilder<int,int> >& fmb,
                      const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> >& rLibrary,
		      const Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> >& lof,
		      const std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >& p_names,
		      bool build_transient_support)
  : fmb_(fmb)
  , responseLibrary_(rLibrary)
  , p_names_(p_names)
  , build_transient_support_(build_transient_support)
  , lof_(lof)
  #ifdef HAVE_STOKHOS
  , sg_lof_(Teuchos::null)
  #endif
{
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  panzer::AssemblyEngine_TemplateBuilder<int,int> builder(fmb,lof);
  ae_tm_.buildObjects(builder);

  // initailize maps, x_dot_init, x0, p_init, g_map, and W_graph
  initializeEpetraObjs();
}

#ifdef HAVE_STOKHOS
panzer::ModelEvaluator_Epetra::
ModelEvaluator_Epetra(const Teuchos::RCP<panzer::FieldManagerBuilder<int,int> >& fmb,
                      const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> >& rLibrary,
		      const Teuchos::RCP<panzer::SGEpetraLinearObjFactory<panzer::Traits,int> >& lof,
		      const std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >& p_names,
		      bool build_transient_support) : 
  fmb_(fmb),
  responseLibrary_(rLibrary),
  p_names_(p_names),
  build_transient_support_(build_transient_support),
  lof_(lof->getEpetraFactory()),
  sg_lof_(lof)
{
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  panzer::AssemblyEngine_TemplateBuilder<int,int> builder(fmb,lof);
  ae_tm_.buildObjects(builder);

  // initailize maps, x_dot_init, x0, p_init, g_map, and W_graph
  initializeEpetraObjs();
}
#endif

void panzer::ModelEvaluator_Epetra::initializeEpetraObjs()
{
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
 
  TEST_FOR_EXCEPTION(lof_==Teuchos::null,std::logic_error,
                     "panzer::ModelEvaluator_Epetra::initializeEpetraObjs: The linear object factory "
                     "was not correctly initialized before calling initializeEpetraObjs.");
  TEST_FOR_EXCEPTION(responseLibrary_==Teuchos::null,std::logic_error,
                     "panzer::ModelEvaluator_Epetra::initializeEpetraObjs: The response library "
                     "was not correctly initialized before calling initializeEpetraObjs.");

  map_x_ = lof_->getMap();
  x0_ = rcp(new Epetra_Vector(*map_x_));
  x_dot_init_ = rcp(new Epetra_Vector(*map_x_));
  x_dot_init_->PutScalar(0.0);

  // setup parameters
  for (std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >::const_iterator p = p_names_.begin(); 
       p != p_names_.end(); ++p) {
    RCP<Epetra_Map> local_map = rcp(new Epetra_LocalMap((*p)->size(), 0, map_x_->Comm())) ;
    p_map_.push_back(local_map);
    RCP<Epetra_Vector> ep_vec = rcp(new Epetra_Vector(*local_map));
    ep_vec->PutScalar(0.0);
    p_init_.push_back(ep_vec);
  }
  
  // setup response maps
  for (std::size_t i=0;i<responseLibrary_->getLabeledResponseCount();i++) {
    RCP<Epetra_Map> local_map = rcp(new Epetra_LocalMap(1, 0, map_x_->Comm()));
    g_map_.push_back(local_map);
  }

  // Initialize the graph for W CrsMatrix object
  W_graph_ = lof_->getGraph();
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
  return x0_;
}

Teuchos::RCP<Epetra_Operator>
panzer::ModelEvaluator_Epetra::create_W() const
{
  return Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_));
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
     // for (int i=0; i<number of prameters; i++)
     //    inArgs.setSupports(IN_ARG_p_sg, i, true);
   
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

#ifdef HAVE_STOKHOS
  if(!Teuchos::is_null(sg_lof_)) {
     outArgs.setSupports(OUT_ARG_f_sg,true);
     outArgs.setSupports(OUT_ARG_W_sg,true);
  }
#endif

  return outArgs;
}


void panzer::ModelEvaluator_Epetra::evalModel( const InArgs& inArgs, 
					       const OutArgs& outArgs ) const
{
  // use x,p to evaluate, f, W and associated responses
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
  using Teuchos::rcp_dynamic_cast;
  
  // Transient or steady-state evaluation is determined by the x_dot
  // vector.  If this RCP is null, then we are doing a steady-state
  // fill.
  bool is_transient = false;
  if (inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_x_dot ))
    is_transient = !Teuchos::is_null(inArgs.get_x_dot());

  // Make sure construction built in transient support
  TEST_FOR_EXCEPTION(is_transient && !build_transient_support_, std::runtime_error,
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
     ghostedContainer_ = rcp_dynamic_cast<EpetraLinearObjContainer>(lof_->buildGhostedLinearObjContainer());
     lof_->initializeGhostedContainer(panzer::EpetraLinearObjContainer::X |
                                      panzer::EpetraLinearObjContainer::DxDt |
                                      panzer::EpetraLinearObjContainer::F |
                                      panzer::EpetraLinearObjContainer::Mat, *ghostedContainer_); 
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
  epGlobalContainer->x = Teuchos::rcp_const_cast<Epetra_Vector>(x);
  if (is_transient)
    epGlobalContainer->dxdt = Teuchos::rcp_const_cast<Epetra_Vector>(x_dot);
  
  if (!Teuchos::is_null(f_out) && !Teuchos::is_null(W_out)) {

    // Set the targets
    epGlobalContainer->f = f_out;
    epGlobalContainer->A = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_out);

    // Zero values in ghosted container objects
    epGhostedContainer->f->PutScalar(0.0);
    epGhostedContainer->A->PutScalar(0.0);

    ae_tm_.getAsObject<panzer::Traits::Jacobian>()->evaluate(ae_inargs);
  }
  else if(!Teuchos::is_null(f_out) && Teuchos::is_null(W_out)) {

    epGlobalContainer->f = f_out;

    // Zero values in ghosted container objects
    epGhostedContainer->f->PutScalar(0.0);

    ae_tm_.getAsObject<panzer::Traits::Residual>()->evaluate(ae_inargs);

    f_out->Update(1.0, *(epGlobalContainer->f), 0.0);
  }
  else if(Teuchos::is_null(f_out) && !Teuchos::is_null(W_out)) {

    // this dummy nonsense is needed only for scattering dirichlet conditions
    if (Teuchos::is_null(dummy_f_))
      dummy_f_ = Teuchos::rcp(new Epetra_Vector(*(this->get_f_map())));
    epGlobalContainer->f = dummy_f_; 
    epGlobalContainer->A = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_out);

    // Zero values in ghosted container objects
    epGhostedContainer->A->PutScalar(0.0);

    ae_tm_.getAsObject<panzer::Traits::Jacobian>()->evaluate(ae_inargs);
  }

  // evaluate responses...uses the stored assembly arguments and containers
  if(requiredResponses)
     evalModel_basic_g(ae_inargs,inArgs,outArgs);
  
  // Holding a rcp to f produces a seg fault in Rythmos when the next
  // f comes in and the resulting dtor is called.  Need to discuss
  // with Ross.  Clearing all references here works!

  epGlobalContainer->x = Teuchos::null;
  epGlobalContainer->dxdt = Teuchos::null;
  epGlobalContainer->f = Teuchos::null;
  epGlobalContainer->A = Teuchos::null;

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

   // build a teuchos comm from an mpi comm
   Teuchos::RCP<Teuchos::Comm<int> > tComm 
      = Teuchos::rcp(new Teuchos::MpiComm<int>(
        Teuchos::opaqueWrapper(dynamic_cast<const Epetra_MpiComm &>(map_x_->Comm()).Comm())));

   const std::vector< Teuchos::RCP<std::vector<panzer::Workset> > >& worksetVec = fmb_->getWorksets();
  
   // convert vector of vectors of worksets into a map from block IDs to worksets
   std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > worksets;
   for(std::size_t i=0;i<worksetVec.size();i++) {
      std::string blockId = (*worksetVec[i])[0].block_id;
      worksets[blockId] = worksetVec[i];
   }

   // evaluator responses
   responseLibrary_->evaluateVolumeFieldManagers<panzer::Traits::Residual>(worksets,ae_inargs,*tComm);

   std::vector<Teuchos::RCP<const Response<panzer::Traits> > > responses;
   responseLibrary_->getLabeledVolumeResponses(responses);
   for(std::size_t i=0;i<responses.size();i++) {
      Teuchos::RCP<Epetra_Vector> vec = outArgs.get_g(i);
      if(vec!=Teuchos::null)
         (*vec)[0] = responses[i]->getValue();
   }
}

bool panzer::ModelEvaluator_Epetra::required_basic_g(const OutArgs & outArgs) const
{
   // determine if any of the outArgs are not null!
   bool activeGArgs = false;
   for(int i=0;i<outArgs.Ng();i++) 
      activeGArgs |= (outArgs.get_g(i)!=Teuchos::null); 

   return activeGArgs;
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
     sg_ghostedContainer_ = rcp_dynamic_cast<SGEpetraLinearObjContainer>(sg_lof_->buildGhostedLinearObjContainer());
     sg_lof_->initializeGhostedContainer(panzer::EpetraLinearObjContainer::X |
                                         panzer::EpetraLinearObjContainer::DxDt |
                                         panzer::EpetraLinearObjContainer::F |
                                         panzer::EpetraLinearObjContainer::Mat, *sg_ghostedContainer_); 
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
     std::cout << x_in->size() << " " << W_out->size() << std::endl;
     std::cout << x_in->size() << " " << f_out->size() << std::endl;
     TEUCHOS_ASSERT(x_in->size()==sgGlobalContainer->size()); 
     TEUCHOS_ASSERT(x_in->size()==sgGhostedContainer->size()); 
     if(!Teuchos::is_null(W_out)) { TEUCHOS_ASSERT(x_in->size()==W_out->size()); }
     if(!Teuchos::is_null(f_out)) { TEUCHOS_ASSERT(x_in->size()==f_out->size()); }
     
     // loop over each coefficient, setting up in/out arguments for the lo containers
     SGEpetraLinearObjContainer::iterator glbItr = sgGlobalContainer->begin();
     SGEpetraLinearObjContainer::iterator ghsItr = sgGhostedContainer->begin();
     for(std::size_t coeff_ind=0;coeff_ind<x_in->size();++coeff_ind,++glbItr,++ghsItr) {
        RCP<EpetraLinearObjContainer> globalContainer = *glbItr;
        RCP<EpetraLinearObjContainer> ghostedContainer = *ghsItr;

        // this is what roger was warning us about!!!!
        globalContainer->x = Teuchos::rcp_const_cast<Epetra_Vector>(x_in->getCoeffPtr(coeff_ind));

        if(!Teuchos::is_null(f_out) && !Teuchos::is_null(W_out)) { // requires residual and jacobian
           globalContainer->f = f_out->getCoeffPtr(coeff_ind); 
           globalContainer->A = rcp_dynamic_cast<Epetra_CrsMatrix>(W_out->getCoeffPtr(coeff_ind)); 
 
           ghostedContainer->f->PutScalar(0.0);
           ghostedContainer->A->PutScalar(0.0);
        }
        else if(!Teuchos::is_null(f_out) && Teuchos::is_null(W_out)) {
           globalContainer->f = f_out->getCoeffPtr(coeff_ind); 
 
           // Zero values in ghosted container objects
           ghostedContainer->f->PutScalar(0.0);
        }
        else if(Teuchos::is_null(f_out) && !Teuchos::is_null(W_out)) {

           // this dummy nonsense is needed only for scattering dirichlet conditions
           if(Teuchos::is_null(dummy_f_))
              dummy_f_ = Teuchos::rcp(new Epetra_Vector(*(this->get_f_map())));
           globalContainer->f = dummy_f_; 
           globalContainer->A = rcp_dynamic_cast<Epetra_CrsMatrix>(W_out->getCoeffPtr(coeff_ind)); 

           // Zero values in ghosted container objects
           ghostedContainer->A->PutScalar(0.0);
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
   return false;
}

void 
panzer::ModelEvaluator_Epetra::
evalModel_sg_g(AssemblyEngineInArgs ae_inargs,const InArgs & inArgs,const OutArgs & outArgs) const
{
   // evaluate this
   TEUCHOS_ASSERT(false); // fail until implemented
}

#endif
