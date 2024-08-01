// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <cmath>

#include "Piro_Epetra_VelocityVerletSolver.hpp"
#include "Piro_Epetra_InvertMassMatrixDecorator.hpp"

#include "EpetraExt_ModelEvaluator.h"


Piro::Epetra::VelocityVerletSolver::VelocityVerletSolver(Teuchos::RCP<Teuchos::ParameterList> appParams_,
                          Teuchos::RCP<EpetraExt::ModelEvaluator> model_,
                          Teuchos::RCP<NOX::Epetra::Observer> observer_ ) :
  appParams(appParams_),
  model(model_),
  observer(observer_)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out = Teuchos::VerboseObjectBase::getDefaultOStream();

  num_p = model->createInArgs().Np();
  num_g = model->createOutArgs().Ng();

  TEUCHOS_TEST_FOR_EXCEPTION(num_p > 1, Teuchos::Exceptions::InvalidParameter,
                     std::endl << "Error in Piro::Epetra::VelocityVerletSolver " <<
                     "Not Implemented for Np>1 : " << num_p << std::endl);
  TEUCHOS_TEST_FOR_EXCEPTION(num_g > 1, Teuchos::Exceptions::InvalidParameter,
                     std::endl << "Error in Piro::Epetra::VelocityVerletSolver " <<
                     "Not Implemented for Ng>1 : " << num_g << std::endl);

  *out << "\nA) Get the base parameter list ...\n";

  RCP<Teuchos::ParameterList> vvPL = sublist(appParams, "Velocity Verlet", true);
  vvPL->validateParameters(*getValidVelocityVerletParameters(),0);

  {
    const std::string verbosity = vvPL->get("Verbosity Level", "VERB_DEFAULT");
    solnVerbLevel = Teuchos::VERB_DEFAULT;
    if      (verbosity == "VERB_NONE")    solnVerbLevel = Teuchos::VERB_NONE;
    else if (verbosity == "VERB_LOW")     solnVerbLevel = Teuchos::VERB_LOW;
    else if (verbosity == "VERB_MEDIUM")  solnVerbLevel = Teuchos::VERB_MEDIUM;
    else if (verbosity == "VERB_HIGH")    solnVerbLevel = Teuchos::VERB_HIGH;
    else if (verbosity == "VERB_EXTREME") solnVerbLevel = Teuchos::VERB_EXTREME;
  }

  numTimeSteps = vvPL->get("Num Time Steps", 10);
  t_final = vvPL->get("Final Time", 0.1);
  t_init  = vvPL->get("Initial Time", 0.0);
  delta_t = t_final / numTimeSteps;
  
  if (vvPL->get("Invert Mass Matrix", false)) {
    Teuchos::RCP<EpetraExt::ModelEvaluator> origModel = model;
    bool lump = vvPL->get("Lump Mass Matrix", false);
    bool isConstMass = vvPL->get("Constant Mass Matrix", false); 
    *out << "\nB) Using InvertMassMatrix Decorator\n";
    model = Teuchos::rcp(new Piro::Epetra::InvertMassMatrixDecorator(
             sublist(vvPL,"Stratimikos", true), origModel, isConstMass, lump, true));
  }
}

Piro::Epetra::VelocityVerletSolver::~VelocityVerletSolver()
{
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::VelocityVerletSolver::get_x_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::VelocityVerletSolver::get_f_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::VelocityVerletSolver::get_p_map(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l >= num_p || l < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::Epetra::VelocityVerletSolver::get_p_map():  " <<
                     "Invalid parameter index l = " <<
                     l << std::endl);
  return model->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::VelocityVerletSolver::get_g_map(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(j > num_g || j < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::Epetra::VelocityVerletSolver::get_g_map():  " <<
                     "Invalid response index j = " <<
                     j << std::endl);

  if (j < num_g) {
    return model->get_g_map(j);
  } else {
    // j == num_g
    return model->get_x_map();
  }
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::VelocityVerletSolver::get_x_init() const
{
  Teuchos::RCP<const Epetra_Vector> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::VelocityVerletSolver::get_p_init(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l >= num_p || l < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::Epetra::VelocityVerletSolver::get_p_init():  " <<
                     "Invalid parameter index l = " <<
                     l << std::endl);
  return model->get_p_init(l);
}

EpetraExt::ModelEvaluator::InArgs Piro::Epetra::VelocityVerletSolver::createInArgs() const
{
  //return underlyingME->createInArgs();
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p);
//  inArgs.setSupports(IN_ARG_x,true);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs Piro::Epetra::VelocityVerletSolver::createOutArgs() const
{
  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());

  // Ng is 1 bigger then model-Ng so that the solution vector can be an outarg
  outArgs.set_Np_Ng(num_p, num_g+1);

  EpetraExt::ModelEvaluator::OutArgs model_outargs = model->createOutArgs();
  //for (int i=0; i<num_g; i++)
  //  for (int j=0; j<num_p; j++)
  //    if (!model_outargs.supports(OUT_ARG_DgDp, i, j).none())
  //      outArgs.setSupports(OUT_ARG_DgDp, i, j,
  //                          DerivativeSupport(DERIV_MV_BY_COL));

  return outArgs;
}

void Piro::Epetra::VelocityVerletSolver::evalModel( const InArgs& inArgs,
                                     const OutArgs& outArgs ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Parse InArgs

  RCP<const Epetra_Vector> p_in;
  if (num_p > 0) p_in = inArgs.get_p(0);

  // Parse OutArgs: always 1 extra
  RCP<Epetra_Vector> g_out; 
  if (num_g > 0) g_out = outArgs.get_g(0); 
  RCP<Epetra_Vector> gx_out = outArgs.get_g(num_g); 

  TEUCHOS_TEST_FOR_EXCEPTION(
     model->get_x_init() == Teuchos::null || model->get_x_dot_init() == Teuchos::null,
                     Teuchos::Exceptions::InvalidParameter,
                     std::endl << "Error in Piro::Epetra::VelocityVerletSolver " <<
                     "Requires x, and x_dot " << std::endl);

  RCP<Epetra_Vector> x = rcp(new Epetra_Vector(*model->get_x_init()));
  RCP<Epetra_Vector> v = rcp(new Epetra_Vector(*model->get_x_dot_init()));
  RCP<Epetra_Vector> a = rcp(new Epetra_Vector(*model->get_f_map()));
  a->PutScalar(0.0); 

  double t = t_init;

  // Observe initial condition
  if (observer != Teuchos::null) observer->observeSolution(*x, t);

  double vo; v->Norm2(&vo);
  *out << "Initial Velocity = " << vo << std::endl;

   if (Teuchos::VERB_MEDIUM <= solnVerbLevel) *out << std::endl;

   EpetraExt::ModelEvaluator::InArgs model_inargs = model->createInArgs();
   EpetraExt::ModelEvaluator::OutArgs model_outargs = model->createOutArgs();
   model_inargs.set_x(x);
   if (num_p > 0)  model_inargs.set_p(0, p_in);

   model_outargs.set_f(a);
   if (g_out != Teuchos::null) model_outargs.set_g(0, g_out);

   double ddt = 0.5 * delta_t * delta_t;

   // Calculate acceleration at time 0
   model->evalModel(model_inargs, model_outargs);

   for (int timeStep = 1; timeStep <= numTimeSteps; timeStep++) {
 
     x->Update(delta_t, *v, ddt, *a, 1.0);
     t += delta_t; model_inargs.set_t(t);

     v->Update(0.5*delta_t, *a, 1.0);

     //calc a(x,t,p);
     model->evalModel(model_inargs, model_outargs);

     v->Update(0.5*delta_t, *a, 1.0);

     // Observe completed time step
     if (observer != Teuchos::null) observer->observeSolution(*x, t);

     if (g_out != Teuchos::null) 
       g_out->Print(*out << "Responses at time step(time) = " << timeStep << "("<<t<<")");
   }

   // return the final solution as an additional g-vector, if requested
   if (gx_out != Teuchos::null)  *gx_out = *x;
}

Teuchos::RCP<const Teuchos::ParameterList>
Piro::Epetra::VelocityVerletSolver::getValidVelocityVerletParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
     Teuchos::rcp(new Teuchos::ParameterList("ValidVelocityVerletParams"));;
  validPL->set<int>("Num Time Steps", 0, "");
  validPL->set<double>("Final Time", 1.0, "");
  validPL->set<double>("Initial Time", 0.0, "");
  validPL->set<std::string>("Verbosity Level", "", "");
  validPL->set<bool>("Invert Mass Matrix", false, "");
  validPL->set<bool>("Lump Mass Matrix", false, "Boolean telling code whether to lump mass matrix");
  validPL->set<bool>("Constant Mass Matrix", false, "Boolean telling code if mass matrix is constant in time");
  validPL->sublist("Stratimikos", false, "");
  return validPL;
}

