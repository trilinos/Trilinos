
#include <cmath>

#include "Piro_Epetra_MatrixFreeDecorator.hpp"


#include "Teuchos_VerboseObject.hpp"


Piro::Epetra::MatrixFreeDecorator::MatrixFreeDecorator(
                          Teuchos::RCP<EpetraExt::ModelEvaluator>& model_,
                          double lambda_) :
  model(model_),
  lambda(lambda_)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  // Allocate storage for base and perturbed residuals
  fBase = rcp(new Epetra_Vector(*(model->get_f_map())));
  fPert = rcp(new Epetra_Vector(*(model->get_f_map())));
}

Piro::Epetra::MatrixFreeDecorator::~MatrixFreeDecorator()
{
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::MatrixFreeDecorator::get_x_map() const
{
  return model->get_x_map();
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::MatrixFreeDecorator::get_f_map() const
{
  return model->get_f_map();
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::MatrixFreeDecorator::get_p_map(int l) const
{
  return model->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::MatrixFreeDecorator::get_g_map(int j) const
{
  return model->get_g_map(j);
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::MatrixFreeDecorator::get_x_init() const
{
  return model->get_x_init();
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::MatrixFreeDecorator::get_p_init(int l) const
{
  return model->get_p_init(l);
}

Teuchos::RCP<Epetra_Operator> Piro::Epetra::MatrixFreeDecorator::create_W() const
{
  return rcp(this, false);
}

EpetraExt::ModelEvaluator::InArgs Piro::Epetra::MatrixFreeDecorator::createInArgs() const
{
  return model->createInArgs();
}

EpetraExt::ModelEvaluator::OutArgs Piro::Epetra::MatrixFreeDecorator::createOutArgs() const
{
  // Augment outArgs to support W
  EpetraExt::ModelEvaluator::OutArgs outArgs =  model->createOutArgs();
  outArgs.setSupports(OUT_ARG_W,true);
  outArgs.set_W_properties( DerivativeProperties(
      DERIV_LINEARITY_UNKNOWN ,DERIV_RANK_FULL , false));

  return outArgs;
}

void Piro::Epetra::MatrixFreeDecorator::evalModel( const InArgs& inArgs,
                                     const OutArgs& outArgs ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  RCP<Epetra_Operator> W_out = outArgs.get_W();

  if (W_out == Teuchos::null) {
    // Just pass through as is: nothing to Decorate
    model->evalModel(inArgs, outArgs);
  }
  else {
    // Do base case for MatrixFree: set f instead of W
    OutArgs modelOutArgs(outArgs);
    InArgs modelInArgs(inArgs);

    // Store f_out in case it was also requested
    RCP<Epetra_Vector> f_out = outArgs.get_f();

    modelOutArgs.set_f(fBase);
    modelOutArgs.set_W(Teuchos::null);

    //Evaluate the underlying model
    model->evalModel(modelInArgs, modelOutArgs);

    // If f_out was requested, return it.
    if (f_out != Teuchos::null) *f_out = *fBase;
    W_out = rcp(this);

  }
}
/***  Epetra Operator Overloads follow ***/
// Piro::Epetra::MatrixFreeDecorator::

int  Piro::Epetra::MatrixFreeDecorator::SetUseTranspose(bool UseTranspose)
{
   TEST_FOR_EXCEPTION(UseTranspose, std:logic_error, 
     " Piro::Epetra::MatrixFreeDecorator cannot support Transpose operators");
}
