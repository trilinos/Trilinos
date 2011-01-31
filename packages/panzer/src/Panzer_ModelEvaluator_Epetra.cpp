#include "Panzer_ModelEvaluator_Epetra.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Epetra_CrsMatrix.h"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"

panzer::ModelEvaluator_Epetra::
ModelEvaluator_Epetra(Teuchos::RCP<panzer::FieldManagerBuilder<int,int> > fmb,
		      Teuchos::RCP<panzer::EpetraLinearObjFactory<panzer::Traits,int> > lof) : 
  fmb_(fmb),
  lof_(lof)
{
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  map_x_ = lof->getMap();
  x0_ = rcp(new Epetra_Vector(*map_x_));
  p_  = rcp(new Epetra_Vector(*map_x_));

  // Initialize the graph for W CrsMatrix object
  W_graph_ = lof->getGraph();

  panzer::AssemblyEngine_TemplateBuilder<int,int> builder(fmb,lof);
  ae_tm_.buildObjects(builder);

  ae_inargs_ = rcp(new panzer::AssemblyEngineInArgs(lof->buildGhostedLinearObjContainer(),
                                                    lof->buildLinearObjContainer()));
  
  // for convenience save the container objects
  ghostedContainer_ = rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(ae_inargs_->ghostedContainer_);
  container_ = rcp_dynamic_cast<panzer::EpetraLinearObjContainer>(ae_inargs_->container_);

/*
  ae_inargs_.x = rcp(new Epetra_Vector(*(lof->getGhostedMap())));
  ae_inargs_.dxdt = rcp(new Epetra_Vector(*(lof->getGhostedMap())));
  ae_inargs_.f = rcp(new Epetra_Vector(*(lof->getGhostedMap())));
  ae_inargs_.j = rcp(new Epetra_CrsMatrix(View, *(lof->getGhostedGraph())));
*/

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

Teuchos::RCP<Epetra_Operator>
panzer::ModelEvaluator_Epetra::create_W() const
{
  return Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_));
}

EpetraExt::ModelEvaluator::InArgs
panzer::ModelEvaluator_Epetra::createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.setSupports(IN_ARG_x_dot,true);
  inArgs.setSupports(IN_ARG_alpha,true);
  inArgs.setSupports(IN_ARG_beta,true);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
panzer::ModelEvaluator_Epetra::createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
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
  //
  // Get the input arguments
  //
  const Epetra_Vector& x = *inArgs.get_x();
  const Epetra_Vector& x_dot = *inArgs.get_x_dot();
  ae_inargs_->alpha = inArgs.get_alpha();
  ae_inargs_->beta = inArgs.get_beta();

  //
  // Get the output arguments
  //
  Epetra_Vector* f_out = outArgs.get_f().get();
  Epetra_Operator* W_out = outArgs.get_W().get();

  // move global solution and time derivative to ghosted solution
  RCP<Epetra_Import> importer = lof_->getGhostedImport();
  ghostedContainer_->x->Import(x,*importer,Insert);
  ghostedContainer_->dxdt->Import(x_dot,*importer,Insert);

  if (f_out && W_out) {
    ae_tm_.getAsObject<panzer::Traits::Jacobian>()->evaluate(*ae_inargs_);

    // move ghosted residual and Jacobian to global versions
    RCP<Epetra_Export> exporter = lof_->getGhostedExport();

    f_out->PutScalar(0.0);
    f_out->Export(*ghostedContainer_->f, *exporter, Add);

    Epetra_CrsMatrix& J = dynamic_cast<Epetra_CrsMatrix&>(*W_out);
    J.PutScalar(0.0);
    J.Export(*ghostedContainer_->A, *exporter, Add);
  }
  else if(f_out && !W_out) {
    ae_tm_.getAsObject<panzer::Traits::Residual>()->evaluate(*ae_inargs_);

    // move ghosted residual to global residual
    RCP<Epetra_Export> exporter = lof_->getGhostedExport();
    f_out->PutScalar(0.0);
    f_out->Export(*ghostedContainer_->f, *exporter, Add);
  }
  else if(!f_out && W_out) {
    ae_tm_.getAsObject<panzer::Traits::Jacobian>()->evaluate(*ae_inargs_);

    // move ghosted Jacobian to global Jacobian
    RCP<Epetra_Export> exporter = lof_->getGhostedExport();
    Epetra_CrsMatrix& J = dynamic_cast<Epetra_CrsMatrix&>(*W_out);
    J.PutScalar(0.0);
    J.Export(*ghostedContainer_->A, *exporter, Add);
  }

}
