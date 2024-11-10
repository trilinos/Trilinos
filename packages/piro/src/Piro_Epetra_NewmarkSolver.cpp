// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <cmath>

#include "Piro_Epetra_NewmarkSolver.hpp"

#include "EpetraExt_ModelEvaluator.h"

#include "Thyra_VectorBase.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Teuchos_VerboseObject.hpp"

Piro::Epetra::NewmarkSolver::NewmarkSolver(
                          Teuchos::RCP<Teuchos::ParameterList> appParams_,
                          Teuchos::RCP<EpetraExt::ModelEvaluator> origModel_,
                          Teuchos::RCP<NOX::Epetra::Observer> observer_ ) :
  appParams(appParams_),
  observer(observer_)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out = Teuchos::VerboseObjectBase::getDefaultOStream();

  *out << "\nA) Get the base parameter list ...\n";

  RCP<Teuchos::ParameterList> trPL = sublist(appParams, "Newmark", true);
  trPL->validateParameters(*getValidNewmarkParameters(),0);

  {
    const std::string verbosity = trPL->get("Verbosity Level", "VERB_DEFAULT");
    solnVerbLevel = Teuchos::VERB_DEFAULT;
    if      (verbosity == "VERB_NONE")    solnVerbLevel = Teuchos::VERB_NONE;
    else if (verbosity == "VERB_LOW")     solnVerbLevel = Teuchos::VERB_LOW;
    else if (verbosity == "VERB_MEDIUM")  solnVerbLevel = Teuchos::VERB_MEDIUM;
    else if (verbosity == "VERB_HIGH")    solnVerbLevel = Teuchos::VERB_HIGH;
    else if (verbosity == "VERB_EXTREME") solnVerbLevel = Teuchos::VERB_EXTREME;
  }

  numTimeSteps = trPL->get("Num Time Steps", 10);
  t_final = trPL->get("Final Time", 0.1);
  t_init  = trPL->get("Initial Time", 0.0);
  beta    = trPL->get("Beta", 0.25);
  gamma   = trPL->get("Gamma", 0.5);
  delta_t = t_final / numTimeSteps;

  *out << "\nB) Using Newmark Decorator and NOX Solver\n";
  model = Teuchos::rcp(new Piro::Epetra::NewmarkDecorator(origModel_,beta,gamma));

  // Construct NOX solver -- will look for NOX sublist -- this must be set!!
  trPL->sublist("NOX").set("Reset Initial Guess",true);
  noxSolver = Teuchos::rcp(new Piro::Epetra::NOXSolver(trPL, model));

  num_p = model->createInArgs().Np();
  num_g = model->createOutArgs().Ng();

  TEUCHOS_TEST_FOR_EXCEPTION(num_p > 1, Teuchos::Exceptions::InvalidParameter,
                     std::endl << "Error in Piro::Epetra::NewmarkSolver " <<
                     "Not Implemented for Np>1 : " << num_p << std::endl);
  TEUCHOS_TEST_FOR_EXCEPTION(num_g > 1, Teuchos::Exceptions::InvalidParameter,
                     std::endl << "Error in Piro::Epetra::NewmarkSolver " <<
                     "Not Implemented for Ng>1 : " << num_g << std::endl);
}

Piro::Epetra::NewmarkSolver::~NewmarkSolver()
{
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::NewmarkSolver::get_x_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::NewmarkSolver::get_f_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::NewmarkSolver::get_p_map(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l >= num_p || l < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::Epetra::NewmarkSolver::get_p_map():  " <<
                     "Invalid parameter index l = " <<
                     l << std::endl);
  return model->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::NewmarkSolver::get_g_map(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(j > num_g || j < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::Epetra::NewmarkSolver::get_g_map():  " <<
                     "Invalid response index j = " <<
                     j << std::endl);

  if (j < num_g) {
    return model->get_g_map(j);
  } else {
    // j == num_g
    return model->get_x_map();
  }
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::NewmarkSolver::get_x_init() const
{
  Teuchos::RCP<const Epetra_Vector> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::NewmarkSolver::get_p_init(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l >= num_p || l < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::Epetra::NewmarkSolver::get_p_init():  " <<
                     "Invalid parameter index l = " <<
                     l << std::endl);
  return model->get_p_init(l);
}

EpetraExt::ModelEvaluator::InArgs Piro::Epetra::NewmarkSolver::createInArgs() const
{
  //return underlyingME->createInArgs();
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs Piro::Epetra::NewmarkSolver::createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());

  // Ng is 1 bigger then model-Ng so that the solution vector can be an outarg
  outArgs.set_Np_Ng(num_p, num_g+1);

  //EpetraExt::ModelEvaluator::OutArgs model_outargs = model->createOutArgs();
  //for (int i=0; i<num_g; i++)
  //  for (int j=0; j<num_p; j++)
  //    if (!model_outargs.supports(OUT_ARG_DgDp, i, j).none())
  //      outArgs.setSupports(OUT_ARG_DgDp, i, j,
  //                          DerivativeSupport(DERIV_MV_BY_COL));

  return outArgs;
}

void Piro::Epetra::NewmarkSolver::evalModel( const InArgs& inArgs,
                                     const OutArgs& outArgs ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  EpetraExt::ModelEvaluator::InArgs nox_inargs = noxSolver->createInArgs();
  EpetraExt::ModelEvaluator::OutArgs nox_outargs = noxSolver->createOutArgs();

  // Parse InArgs

  RCP<const Epetra_Vector> p_in;
  if (num_p > 0) {
    p_in = inArgs.get_p(0);
    nox_inargs.set_p(0, p_in);
  }

  // Parse OutArgs: always 1 extra
  RCP<Epetra_Vector> g_out;
  if (num_g > 0) {
    g_out = outArgs.get_g(0);
    nox_outargs.set_g(0, g_out);
  }
  RCP<Epetra_Vector> gx_out = outArgs.get_g(num_g);
  if (Teuchos::is_null(gx_out)) {
    // Solution not requested by caller as a response, create local temporary instead
    gx_out = rcp(new Epetra_Vector(*model->get_x_map()));
  }
  nox_outargs.set_g(num_g, gx_out);

  TEUCHOS_TEST_FOR_EXCEPTION(
     model->get_x_init() == Teuchos::null || model->get_x_dot_init() == Teuchos::null
         || model->get_x_dotdot_init() == Teuchos::null,
                     Teuchos::Exceptions::InvalidParameter,
                     std::endl << "Error in Piro::Epetra::NewmarkSolver " <<
                     "Requires x, x_dot, and x_dotdot: " << std::endl);

  RCP<Epetra_Vector> a = rcp(new Epetra_Vector(*model->get_x_dotdot_init()));
  RCP<Epetra_Vector> v = rcp(new Epetra_Vector(*model->get_x_dot_init()));
  RCP<Epetra_Vector> d = rcp(new Epetra_Vector(*model->get_f_map()));
  RCP<Epetra_Vector> v_pred = rcp(new Epetra_Vector(*model->get_f_map()));
  RCP<Epetra_Vector> d_pred = rcp(new Epetra_Vector(*model->get_f_map()));
  RCP<Epetra_Vector> a_old = rcp(new Epetra_Vector(*model->get_f_map()));

   double t = t_init;

   // Observe initial condition
   if (observer != Teuchos::null) observer->observeSolution(*d,t);

   double nrm;
   v->Norm2(&nrm); *out << "Initial Velocity = " << nrm << std::endl;

   //calculate intial acceleration
   {
     *v_pred = *v;
     *d_pred = *d;
     model->injectData(a, v_pred, d_pred, /*delta_t=*/ 0.0, t);
     noxSolver->evalModel(nox_inargs, nox_outargs);
     *a = *gx_out;
     a->Norm2(&nrm); *out << "Calculated a_init = " << nrm << std::endl;
   }

   // Start integration loop
   for (int timeStep = 1; timeStep <= numTimeSteps; timeStep++) {

     t += delta_t;

     *a_old = *a;

     // compute velocity and displacement predictors
     *v_pred = *v;
      v_pred->Update(delta_t*(1.0-gamma),*a_old,1.0);

     *d_pred = *d;
      d_pred->Update(delta_t,*v, delta_t*delta_t/2.0*(1.0-2.0*beta),*a_old,1.0);

     model->injectData(a, v_pred, d_pred, delta_t, t);
     noxSolver->evalModel(nox_inargs, nox_outargs);

     // Copy out final solution from nonlinear solver
     *a =  *gx_out;

     // Compute v and d and new conditions
     v->Update(1.0, *v_pred, delta_t*gamma, *a, 0.0);
     d->Update(1.0, *d_pred, beta*delta_t*delta_t, *a, 0.0);

     // Observe completed time step
     if (observer != Teuchos::null) observer->observeSolution(*d,t);

     if (g_out != Teuchos::null)
       g_out->Print(*out << "Responses at time step(time) = " << timeStep << "("<<t<<")\n");
   }
}

Teuchos::RCP<const Teuchos::ParameterList>
Piro::Epetra::NewmarkSolver::getValidNewmarkParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
     Teuchos::rcp(new Teuchos::ParameterList("ValidNewmarkParams"));;
  validPL->set<int>("Num Time Steps", 0, "");
  validPL->set<double>("Gamma", 0.5, "");
  validPL->set<double>("Beta", 0.25, "");
  validPL->set<double>("Final Time", 1.0, "");
  validPL->set<double>("Initial Time", 0.0, "");
  validPL->set<std::string>("Verbosity Level", "", "");
  validPL->sublist("NOX", false, "");
  return validPL;
}

/****************************************************************************/

Piro::Epetra::NewmarkDecorator::NewmarkDecorator(
                          Teuchos::RCP<EpetraExt::ModelEvaluator>& model_, 
                          double beta_, double gamma_) :
  model(model_),
  beta(beta_),
  gamma(gamma_),
  delta_t(0.0),
  time(0.0)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  displacement = Teuchos::rcp(new Epetra_Vector(*model->get_x_map()));
  velocity = Teuchos::rcp(new Epetra_Vector(*model->get_x_map()));
  velocity_pred = Teuchos::rcp(new Epetra_Vector(*model->get_x_init()));
  displacement_pred = Teuchos::rcp(new Epetra_Vector(*model->get_x_init()));
  acceleration_save = Teuchos::rcp(new Epetra_Vector(*model->get_x_init()));

  Teuchos::RCP<Teuchos::FancyOStream> out
     = Teuchos::VerboseObjectBase::getDefaultOStream();
}

Piro::Epetra::NewmarkDecorator::~NewmarkDecorator()
{
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::NewmarkDecorator::get_x_map() const
{
  return model->get_x_map();
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::NewmarkDecorator::get_f_map() const
{
  return model->get_f_map();
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::NewmarkDecorator::get_p_map(int l) const
{
  return model->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::NewmarkDecorator::get_g_map(int j) const
{
  return model->get_g_map(j);
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::NewmarkDecorator::get_x_init() const
{
  return acceleration_save;
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::NewmarkDecorator::get_x_dot_init() const
{
  return model->get_x_dot_init();
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::NewmarkDecorator::get_x_dotdot_init() const
{
  return model->get_x_dotdot_init();
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::NewmarkDecorator::get_p_init(int l) const
{
  return model->get_p_init(l);
}

Teuchos::RCP<Epetra_Operator> Piro::Epetra::NewmarkDecorator::create_W() const
{
  return model->create_W();
}

EpetraExt::ModelEvaluator::InArgs Piro::Epetra::NewmarkDecorator::createInArgs() const
{
  return model->createInArgs();
}

EpetraExt::ModelEvaluator::OutArgs Piro::Epetra::NewmarkDecorator::createOutArgs() const
{
  return model->createOutArgs();
}

void Piro::Epetra::NewmarkDecorator::injectData(
    const Teuchos::RCP<Epetra_Vector>& a_,
    const Teuchos::RCP<Epetra_Vector>& v_pred_,
    const Teuchos::RCP<Epetra_Vector>& d_pred_,
    double delta_t_,double time_)
{
  *acceleration_save = *a_;
  *velocity_pred     = *v_pred_;
  *displacement_pred = *d_pred_;
  delta_t = delta_t_;
  time = time_;
}

void Piro::Epetra::NewmarkDecorator::evalModel( const InArgs& inArgs,
                                     const OutArgs& outArgs ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Copy outArgs; add time term
  OutArgs modelOutArgs(outArgs);
  InArgs modelInArgs(inArgs);

  // the solution variable in NOX is the acceleration, a_{n+1}
  modelInArgs.set_x_dotdot(inArgs.get_x());

  // compute the velocity, v_{n+1}(a_{n+1}) = velocity_{pred} + \gamma dt a_{n+1}
  velocity->Update(1.0, *velocity_pred, delta_t*gamma, *inArgs.get_x(), 0.0);
  modelInArgs.set_x_dot(velocity);

  // compute the displacement, d_{n+1}(a_{n+1}) = displacement_{pred} + \beta dt^2 a_{n+1}
  displacement->Update(1.0, *displacement_pred, beta*delta_t*delta_t, *inArgs.get_x(), 0.0);
  modelInArgs.set_x(displacement);

  modelInArgs.set_omega(1.0);                 // da/da
  modelInArgs.set_alpha(gamma*delta_t);       // dv/da
  modelInArgs.set_beta(beta*delta_t*delta_t); // dd/da
  modelInArgs.set_t(time);

  // Compute  F(a, v(a), d(a), t) = 0
  model->evalModel(modelInArgs, modelOutArgs);
}
