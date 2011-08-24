// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include <cmath>

#include "Piro_Epetra_TrapezoidRuleSolver.hpp"
#include "Piro_ValidPiroParameters.hpp"

#include "EpetraExt_ModelEvaluator.h"

#include "Thyra_VectorBase.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Teuchos_VerboseObject.hpp"

Piro::Epetra::TrapezoidRuleSolver::TrapezoidRuleSolver(
                          Teuchos::RCP<Teuchos::ParameterList> appParams_,
                          Teuchos::RCP<EpetraExt::ModelEvaluator> origModel_,
                          Teuchos::RCP<NOX::Epetra::Observer> observer_ ) :
  appParams(appParams_),
  observer(observer_)
{
  //appParams->validateParameters(*Piro::getValidPiroParameters(),0);

  using Teuchos::RCP;
  using Teuchos::rcp;

  out = Teuchos::VerboseObjectBase::getDefaultOStream();

  *out << "\nA) Get the base parameter list ...\n";

  RCP<Teuchos::ParameterList> trPL = sublist(appParams, "Trapezoid Rule", true);
  trPL->validateParameters(*getValidTrapezoidRuleParameters(),0);

  {
    const string verbosity = trPL->get("Verbosity Level", "VERB_DEFAULT");
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
  delta_t = t_final / numTimeSteps;
  
  *out << "\nB) Using Trapezoid Decorator and NOX Solver\n";
  model = Teuchos::rcp(new Piro::Epetra::TrapezoidDecorator(origModel_));

  // Construct NOX solver -- will look for NOX sublist -- this must be set!!
  trPL->sublist("NOX").set("Reset Initial Guess",true);
  noxSolver = Teuchos::rcp(new Piro::Epetra::NOXSolver(trPL, model));

  num_p = model->createInArgs().Np();
  num_g = model->createOutArgs().Ng();

  TEST_FOR_EXCEPTION(num_p > 1, Teuchos::Exceptions::InvalidParameter,
                     std::endl << "Error in Piro::Epetra::TrapezoidRuleSolver " <<
                     "Not Implemented for Np>1 : " << num_p << std::endl);
  TEST_FOR_EXCEPTION(num_g > 1, Teuchos::Exceptions::InvalidParameter,
                     std::endl << "Error in Piro::Epetra::TrapezoidRuleSolver " <<
                     "Not Implemented for Ng>1 : " << num_g << std::endl);
}

Piro::Epetra::TrapezoidRuleSolver::~TrapezoidRuleSolver()
{
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::TrapezoidRuleSolver::get_x_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::TrapezoidRuleSolver::get_f_map() const
{
  Teuchos::RCP<const Epetra_Map> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::TrapezoidRuleSolver::get_p_map(int l) const
{
  TEST_FOR_EXCEPTION(l >= num_p || l < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::Epetra::TrapezoidRuleSolver::get_p_map():  " <<
                     "Invalid parameter index l = " <<
                     l << std::endl);
  return model->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::TrapezoidRuleSolver::get_g_map(int j) const
{
  TEST_FOR_EXCEPTION(j > num_g || j < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::Epetra::TrapezoidRuleSolver::get_g_map():  " <<
                     "Invalid response index j = " <<
                     j << std::endl);

  if      (j < num_g) return model->get_g_map(j);
  else if (j == num_g) return model->get_x_map();
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::TrapezoidRuleSolver::get_x_init() const
{
  Teuchos::RCP<const Epetra_Vector> neverused;
  return neverused;
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::TrapezoidRuleSolver::get_p_init(int l) const
{
  TEST_FOR_EXCEPTION(l >= num_p || l < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::Epetra::TrapezoidRuleSolver::get_p_init():  " <<
                     "Invalid parameter index l = " <<
                     l << std::endl);
  return model->get_p_init(l);
}

EpetraExt::ModelEvaluator::InArgs Piro::Epetra::TrapezoidRuleSolver::createInArgs() const
{
  //return underlyingME->createInArgs();
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs Piro::Epetra::TrapezoidRuleSolver::createOutArgs() const
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

void Piro::Epetra::TrapezoidRuleSolver::evalModel( const InArgs& inArgs,
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
  nox_outargs.set_g(num_g, gx_out);


  RCP<Epetra_Vector> x = rcp(new Epetra_Vector(*model->get_x_init()));
  RCP<Epetra_Vector> v = rcp(new Epetra_Vector(*model->get_x_dot_init()));
  RCP<Epetra_Vector> a = rcp(new Epetra_Vector(*model->get_f_map()));
  RCP<Epetra_Vector> x_pred = rcp(new Epetra_Vector(*model->get_f_map()));
  RCP<Epetra_Vector> a_old = rcp(new Epetra_Vector(*model->get_f_map()));

  TEST_FOR_EXCEPTION(v == Teuchos::null || x == Teuchos::null, 
                     Teuchos::Exceptions::InvalidParameter,
                     std::endl << "Error in Piro::Epetra::TrapezoidRuleSolver " <<
                     "Requires initial x and x_dot: " << std::endl);
   double nrm;
   v->Norm2(&nrm); *out << "Initial Velocity = " << nrm << endl;

   double t = t_init;

   //calculate intial acceleration using small time step (1.0e-3*delta_t)
   {
     double pert= 1.0e6 * 4.0 / (delta_t * delta_t);
     *x_pred = *x;
     model->injectData(x_pred, x_pred, pert, t);
     noxSolver->evalModel(nox_inargs, nox_outargs);
     a->Update(pert, *gx_out,  -pert, *x_pred,0.0);
     a->Norm2(&nrm); *out << "Calculated a_init = " << nrm << endl;
   }

   // Start integration loop
   double fdt2 = 4.0 / (delta_t * delta_t);
   double dt2f =  delta_t * delta_t / 4.0;
   double hdt  =  delta_t/ 2.0;

   for (int timeStep = 1; timeStep <= numTimeSteps; timeStep++) {
 
     t += delta_t;
 
     *a_old = *a;
     *x_pred = *x;
     x_pred->Update(delta_t, *v, dt2f, *a, 1.0);
     model->injectData(x, x_pred, fdt2, t);

     noxSolver->evalModel(nox_inargs, nox_outargs);
     // Copy out final solution from nonlinear solver
     *x =  *gx_out;
     // Compute a and v and new conditions
     a->Update(fdt2, *x,  -fdt2, *x_pred,0.0);
     v->Update(hdt, *a, hdt, *a_old, 1.0); 

     if (observer != Teuchos::null) observer->observeSolution(*x,t);
     if (g_out != Teuchos::null) 
       g_out->Print(*out << "Responses at time step(time) = " << timeStep << "("<<t<<")\n");
   }
}

Teuchos::RCP<const Teuchos::ParameterList>
Piro::Epetra::TrapezoidRuleSolver::getValidTrapezoidRuleParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
     Teuchos::rcp(new Teuchos::ParameterList("ValidTrapezoidRuleParams"));;
  validPL->set<int>("Num Time Steps", 0, "");
  validPL->set<double>("Final Time", 1.0, "");
  validPL->set<double>("Initial Time", 0.0, "");
  validPL->set<std::string>("Verbosity Level", "", "");
  validPL->sublist("NOX", false, "");
  return validPL;
}

/****************************************************************************/

Piro::Epetra::TrapezoidDecorator::TrapezoidDecorator(
                          Teuchos::RCP<EpetraExt::ModelEvaluator>& model_) :
  model(model_)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  xDotDot = Teuchos::rcp(new Epetra_Vector(*model->get_x_map()));
  x_pred = Teuchos::rcp(new Epetra_Vector(*model->get_x_init()));
  x_save = Teuchos::rcp(new Epetra_Vector(*model->get_x_init()));

  Teuchos::RCP<Teuchos::FancyOStream> out
     = Teuchos::VerboseObjectBase::getDefaultOStream();
}

Piro::Epetra::TrapezoidDecorator::~TrapezoidDecorator()
{
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::TrapezoidDecorator::get_x_map() const
{
  return model->get_x_map();
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::TrapezoidDecorator::get_f_map() const
{
  return model->get_f_map();
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::TrapezoidDecorator::get_p_map(int l) const
{
  return model->get_p_map(l);
}

Teuchos::RCP<const Epetra_Map> Piro::Epetra::TrapezoidDecorator::get_g_map(int j) const
{
  return model->get_g_map(j);
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::TrapezoidDecorator::get_x_init() const
{
  return x_save;
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::TrapezoidDecorator::get_x_dot_init() const
{
  return model->get_x_dot_init();
}

Teuchos::RCP<const Epetra_Vector> Piro::Epetra::TrapezoidDecorator::get_p_init(int l) const
{
  return model->get_p_init(l);
}

Teuchos::RCP<Epetra_Operator> Piro::Epetra::TrapezoidDecorator::create_W() const
{
  return model->create_W();
}

EpetraExt::ModelEvaluator::InArgs Piro::Epetra::TrapezoidDecorator::createInArgs() const
{
  return model->createInArgs();
}

EpetraExt::ModelEvaluator::OutArgs Piro::Epetra::TrapezoidDecorator::createOutArgs() const
{
  return model->createOutArgs();
}

void Piro::Epetra::TrapezoidDecorator::injectData(
    const Teuchos::RCP<Epetra_Vector>& x_, 
    const Teuchos::RCP<Epetra_Vector>& x_pred_, 
    double fdt2_, double time_) 
{
  *x_save = *x_;
  *x_pred = *x_pred_;
  fdt2 = fdt2_;
  time = time_;
}

void Piro::Epetra::TrapezoidDecorator::evalModel( const InArgs& inArgs,
                                     const OutArgs& outArgs ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Copy outArgs; add time term 
  OutArgs modelOutArgs(outArgs);
  InArgs modelInArgs(inArgs);

  xDotDot->Update(fdt2, *inArgs.get_x(), -fdt2, *x_pred, 0.0);
  // WARNING:: Model must use  x_dot  as xDotDot!!!
  modelInArgs.set_x_dot(xDotDot); 

  modelInArgs.set_alpha(fdt2);  // fdt2 = 4/(dt)^2
  modelInArgs.set_beta(1.0);
  modelInArgs.set_t(time);

  //Evaluate the underlying model
  model->evalModel(modelInArgs, modelOutArgs);
}
