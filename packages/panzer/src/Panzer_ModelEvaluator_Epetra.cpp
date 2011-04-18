#include "Panzer_ModelEvaluator_Epetra.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Epetra_CrsMatrix.h"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Epetra_LocalMap.h"

panzer::ModelEvaluator_Epetra::
ModelEvaluator_Epetra(const Teuchos::RCP<panzer::FieldManagerBuilder<int,int> >& fmb,
		      const Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> >& lof,
		      const std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >& p_names,
		      bool build_transient_support) : 
  fmb_(fmb),
  lof_(lof),
  p_names_(p_names),
  build_transient_support_(build_transient_support)
{
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  map_x_ = lof->getMap();
  x0_ = rcp(new Epetra_Vector(*map_x_));
  x_dot_init_ = rcp(new Epetra_Vector(*map_x_));
  x_dot_init_->PutScalar(0.0);
  
  
  for (std::vector<Teuchos::RCP<Teuchos::Array<std::string> > >::const_iterator p = p_names_.begin(); 
       p != p_names_.end(); ++p) {
    RCP<Epetra_Map> local_map = rcp(new Epetra_LocalMap((*p)->size(), 0, map_x_->Comm())) ;
    p_map_.push_back(local_map);
    RCP<Epetra_Vector> ep_vec = rcp(new Epetra_Vector(*local_map));
    ep_vec->PutScalar(0.0);
    p_init_.push_back(ep_vec);
  }

  // Initialize the graph for W CrsMatrix object
  W_graph_ = lof->getGraph();

  panzer::AssemblyEngine_TemplateBuilder<int,int> builder(fmb,lof);
  ae_tm_.buildObjects(builder);

  ae_inargs_ = rcp(new panzer::AssemblyEngineInArgs(lof->buildGhostedLinearObjContainer(),
                                                    lof->buildLinearObjContainer()));
  
  // for convenience save the container objects
  ghostedContainer_ = rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(ae_inargs_->ghostedContainer_);
  container_ = rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(ae_inargs_->container_);

  isInitialized_ = true;
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
  outArgs.set_Np_Ng(p_init_.size(), 0);
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
  outArgs.set_W_properties(
    DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_FULL
      ,true // supportsAdjoint
      )
    );
  return outArgs;
}

void panzer::ModelEvaluator_Epetra::evalModel( const InArgs& inArgs, 
					       const OutArgs& outArgs ) const
{
  using Teuchos::dyn_cast;
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

//   if (is_transient)
//     cout << "ME:evalModel() is transient (x_dot is non-null)!" << endl;
//   else 
//     cout << "ME:evalModel() is steady-state (x_dot is null)!" << endl;

  //
  // Get the input arguments
  //
  const RCP<const Epetra_Vector> x = inArgs.get_x();
  RCP<const Epetra_Vector> x_dot;
  ae_inargs_->alpha = 0.0;
  ae_inargs_->beta = 1.0;
  if (is_transient) {
    x_dot = inArgs.get_x_dot();
    ae_inargs_->alpha = inArgs.get_alpha();
    ae_inargs_->beta = inArgs.get_beta();
    ae_inargs_->time = inArgs.get_t();
  }
  
//   cout << "ME: alpha = " << ae_inargs_->alpha << endl;
//   cout << "ME: beta = " << ae_inargs_->beta << endl;

  //
  // Get the output arguments
  //
  const RCP<Epetra_Vector> f_out = outArgs.get_f();
  const RCP<Epetra_Operator> W_out = outArgs.get_W();
  
  // Zero out global container objects
  const RCP<panzer::EpetraLinearObjContainer> epGlobalContainer = 
    Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(ae_inargs_->container_);
  TEUCHOS_ASSERT(!Teuchos::is_null(epGlobalContainer));
  epGlobalContainer->x = Teuchos::null;
  epGlobalContainer->dxdt = Teuchos::null;
  epGlobalContainer->f = Teuchos::null;
  epGlobalContainer->A = Teuchos::null;

  // Ghosted container objects are zeroed out below only if needed for
  // a particular calculation.  This makes it more efficient thatn
  // zeroing out all objects in the container here.
  const RCP<panzer::EpetraLinearObjContainer> epGhostedContainer = 
    Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(ae_inargs_->ghostedContainer_);
  
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
    //std::cout << "F AND J" << std::endl;
    
    // Set the targets
    epGlobalContainer->f = f_out;
    epGlobalContainer->A = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_out);

    // Zero values in ghosted container objects
    epGhostedContainer->f->PutScalar(0.0);
    epGhostedContainer->A->PutScalar(0.0);

    ae_tm_.getAsObject<panzer::Traits::Jacobian>()->evaluate(*ae_inargs_);
    
    //f_out->Print(std::cout);
    //Epetra_CrsMatrix& J = dynamic_cast<Epetra_CrsMatrix&>(*W_out);
    //J.Print(std::cout);
  }
  else if(!Teuchos::is_null(f_out) && Teuchos::is_null(W_out)) {
    //std::cout << "F only" << std::endl;

    epGlobalContainer->f = f_out;

    // Zero values in ghosted container objects
    epGhostedContainer->f->PutScalar(0.0);

    ae_tm_.getAsObject<panzer::Traits::Residual>()->evaluate(*ae_inargs_);

    f_out->Update(1.0, *(container_->f), 0.0);

    //f_out->Print(std::cout);
  }
  else if(Teuchos::is_null(f_out) && !Teuchos::is_null(W_out)) {
    //std::cout << "J only" << std::endl;

    if (Teuchos::is_null(dummy_f_))
      dummy_f_ = Teuchos::rcp(new Epetra_Vector(*(this->get_f_map())));

    epGlobalContainer->f = dummy_f_;
    epGlobalContainer->A = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_out);

    // Zero values in ghosted container objects
    epGhostedContainer->A->PutScalar(0.0);

    ae_tm_.getAsObject<panzer::Traits::Jacobian>()->evaluate(*ae_inargs_);

    //Epetra_CrsMatrix& J = dynamic_cast<Epetra_CrsMatrix&>(*W_out);
    //J.Print(std::cout);
  }

  // Holding a rcp to f produces a seg fault in Rythmos when the next
  // f comes in and the resulting dtor is called.  Need to discuss
  // with Ross.  Clearing all references here works!
  epGlobalContainer->x = Teuchos::null;
  epGlobalContainer->dxdt = Teuchos::null;
  epGlobalContainer->f = Teuchos::null;
  epGlobalContainer->A = Teuchos::null;

}
