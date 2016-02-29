// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#include <cmath>

#include "Thyra_DefaultModelEvaluatorWithSolveFactory.hpp"
#include "Piro_TransientDecorator.hpp"

#define ALBANY_BUILD
#include "Piro_InvertMassMatrixDecorator.hpp"

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Piro::VelocityVerletSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
VelocityVerletSolver(const Teuchos::RCP<Teuchos::ParameterList> &appParams_,
                          const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model_,
                          const Teuchos::RCP<Thyra::AdaptiveSolutionManager> &solMgr_,
                          const Teuchos::RCP<Piro::ObserverBase<Scalar> > &observer_ ) :
  appParams(appParams_),
  model(model_),
  observer(observer_),
  solMgr(solMgr_),
  out(Teuchos::VerboseObjectBase::getDefaultOStream())
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  num_p = model->Np();
  num_g = model->Ng();

  TEUCHOS_TEST_FOR_EXCEPTION(num_p > 1, Teuchos::Exceptions::InvalidParameter,
                     std::endl << "Error in Piro::VelocityVerletSolver " <<
                     "Not Implemented for Np>1 : " << num_p << std::endl);
  TEUCHOS_TEST_FOR_EXCEPTION(num_g > 1, Teuchos::Exceptions::InvalidParameter,
                     std::endl << "Error in Piro::VelocityVerletSolver " <<
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
    Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > origModel = model;
    bool lump=vvPL->get("Lump Mass Matrix", false);
    *out << "\nB) Using InvertMassMatrix Decorator\n";
    model = Teuchos::rcp(new Piro::InvertMassMatrixDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
             sublist(vvPL,"Stratimikos", true), origModel,
             true, lump, true));
  }
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::VelocityVerletSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::get_x() const
{

  return model->get_x();

}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::VelocityVerletSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::get_x_dot() const
{

  return model->get_x_dot();

}


template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::VelocityVerletSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::get_p_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      l >= num_p || l < 0,
      Teuchos::Exceptions::InvalidParameter,
      std::endl <<
      "Error in Piro::VelocityVerletSolver::get_p_map():  " <<
      "Invalid parameter index l = " <<
      l << std::endl);

  return model->get_p_space(l);
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::VelocityVerletSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::get_g_space(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      j > num_g || j < 0,
      Teuchos::Exceptions::InvalidParameter,
      std::endl <<
      "Error in Piro::VelocityVerletSolver::get_g_map():  " <<
      "Invalid response index j = " <<
      j << std::endl);

  if (j < num_g) {
    return model->get_g_space(j);
  } else {
    // j == num_g
    return model->get_x_space();
  }
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::VelocityVerletSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNominalValues() const
{
  Thyra::ModelEvaluatorBase::InArgs<Scalar> result = this->createInArgs();
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> modelNominalValues = model->getNominalValues();
  for (int l = 0; l < num_p; ++l) {
    result.set_p(l, modelNominalValues.get_p(l));
  }
  return result;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::VelocityVerletSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createInArgs() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p);
  return inArgs;
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
Piro::VelocityVerletSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createOutArgsImpl() const
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());

  // One additional response slot for the solution vector
  outArgs.set_Np_Ng(num_p, num_g + 1);

  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = model->createOutArgs();

  if (num_p > 0) {
    // Only one parameter supported
    const int l = 0;

    // Computing the DxDp sensitivity for a transient problem currently requires the evaluation of
    // the mutilivector-based, Jacobian-oriented DfDp derivatives of the underlying transient model.
    const Thyra::ModelEvaluatorBase::DerivativeSupport model_dfdp_support =
      modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, l);
    if (!model_dfdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
      // Ok to return early since only one parameter supported
      return outArgs;
    }

    // Solution sensitivity
    outArgs.setSupports(
        Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,
        num_g,
        l,
        Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);

    if (num_g > 0) {
      // Only one response supported
      const int j = 0;

      const Thyra::ModelEvaluatorBase::DerivativeSupport model_dgdx_support =
        modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, j);
      if (!model_dgdx_support.none()) {
        const Thyra::ModelEvaluatorBase::DerivativeSupport model_dgdp_support =
          modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
        // Response sensitivity
        Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support;
        if (model_dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
          dgdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
        }
        if (model_dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP)) {
          dgdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
        }
        outArgs.setSupports(
            Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,
            j,
            l,
            dgdp_support);
      }
    }
  }

  return outArgs;

}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void
Piro::VelocityVerletSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // TODO: Support more than 1 parameter and 1 response
  const int j = 0;
  const int l = 0;

  // Parse InArgs
  RCP<const Thyra::VectorBase<Scalar> > p_in;
  if (num_p > 0) {
    p_in = inArgs.get_p(l);
  }

  // Parse OutArgs
  RCP<Thyra::VectorBase<Scalar> > g_out;
  if (num_g > 0) {
    g_out = outArgs.get_g(j);
  }
  const RCP<Thyra::VectorBase<Scalar> > gx_out = outArgs.get_g(num_g);

// create a new vector and fill it with the contents of model->get_x()
  // Build a multivector holding x (0th vector), v (1st vector), and a (2nd vector)
//  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > soln = createMembers(model->get_x_space(), 3);

// create a new vector and fill it with the contents of model->get_x()
/*
  Teuchos::RCP<Thyra::VectorBase<Scalar> > x = soln->col(0);
  assign(x.ptr(), *model->getNominalValues().get_x());
  Teuchos::RCP<Thyra::VectorBase<Scalar> > v = soln->col(1);
  assign(v.ptr(), *model->getNominalValues().get_x_dot());
  Teuchos::RCP<Thyra::VectorBase<Scalar> > a = soln->col(2);
  assign(a.ptr(), *model->get_x_dotdot());
*/

  Teuchos::RCP<Thyra::VectorBase<Scalar> > x = model->getNominalValues().get_x()->clone_v();
  Teuchos::RCP<Thyra::VectorBase<Scalar> > v = model->getNominalValues().get_x_dot()->clone_v();

// Note that Thyra doesn't have x_dotdot - go get it from the transient decorator around the Albany model
//  Teuchos::RCP<Thyra::VectorBase<Scalar> > a = model->get_x_dotdot()->clone_v();

  Teuchos::RCP<Thyra::DefaultModelEvaluatorWithSolveFactory<Scalar> >
      DMEWSF(Teuchos::rcp_dynamic_cast<Thyra::DefaultModelEvaluatorWithSolveFactory<Scalar> >(model));

  Teuchos::RCP<const Piro::TransientDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > dec =
       Teuchos::rcp_dynamic_cast<const Piro::TransientDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
           (DMEWSF->getUnderlyingModel());

  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(dec), std::logic_error,
      "Underlying model in VelovityVerletSolver does not cast to a Piro::TransientDecorator<Scalar, LocalOrdinal, GlobalOrdinal, Node>"
      << std::endl);

  Teuchos::RCP<Thyra::VectorBase<Scalar> > a = dec->get_x_dotdot()->clone_v();

  RCP<Thyra::VectorBase<Scalar> > finalSolution;

  // Zero out the acceleration vector
  put_scalar(0.0, a.ptr());

  TEUCHOS_TEST_FOR_EXCEPTION(v == Teuchos::null || x == Teuchos::null,
                     Teuchos::Exceptions::InvalidParameter,
                     std::endl << "Error in Piro::VelocityVerletSolver " <<
                     "Requires initial x and x_dot: " << std::endl);

  Scalar t = t_init;

  // Observe initial condition
//  if (observer != Teuchos::null) observer->observeSolution(*soln, t);

  Scalar vo = norm_2(*v);
  *out << "Initial Velocity = " << vo << std::endl;

   if (Teuchos::VERB_MEDIUM <= solnVerbLevel) *out << std::endl;

   Thyra::ModelEvaluatorBase::InArgs<Scalar> model_inargs = model->createInArgs();
   Thyra::ModelEvaluatorBase::OutArgs<Scalar> model_outargs = model->createOutArgs();
   model_inargs.set_x(x);
   if (num_p > 0)  model_inargs.set_p(0, p_in);

   model_outargs.set_f(a);
   if (g_out != Teuchos::null) model_outargs.set_g(0, g_out);

   Scalar ddt = 0.5 * delta_t * delta_t;

   // Calculate acceleration at time 0
   model->evalModel(model_inargs, model_outargs);

   for (int timeStep = 1; timeStep <= numTimeSteps; timeStep++) {

//     x->Update(delta_t, *v, ddt, *a, 1.0);
     Vp_StV(x.ptr(), delta_t, *v);
     Vp_StV(x.ptr(), ddt, *a);

     t += delta_t;
     model_inargs.set_t(t);

//     v->Update(0.5*delta_t, *a, 1.0);
     Vp_StV(v.ptr(), 0.5 * delta_t, *a);

     //calc a(x,t,p);
     model->evalModel(model_inargs, model_outargs);

//     v->Update(0.5*delta_t, *a, 1.0);
     Vp_StV(v.ptr(), 0.5 * delta_t, *a);

     // Observe completed time step
     if (observer != Teuchos::null) observer->observeSolution(*x, t);

   }

   // return the final solution as an additional g-vector, if requested
   if (finalSolution != Teuchos::null)  finalSolution = x->clone_v();


  // Return the final solution as an additional g-vector, if requested
  if (Teuchos::nonnull(gx_out)) {
    Thyra::copy(*finalSolution, gx_out.ptr());
  }
}

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
Teuchos::RCP<const Teuchos::ParameterList>
Piro::VelocityVerletSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getValidVelocityVerletParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
     Teuchos::rcp(new Teuchos::ParameterList("ValidVelocityVerletParams"));;
  validPL->set<int>("Num Time Steps", 0, "");
  validPL->set<double>("Final Time", 1.0, "");
  validPL->set<double>("Initial Time", 0.0, "");
  validPL->set<std::string>("Verbosity Level", "", "");
  validPL->set<bool>("Invert Mass Matrix", false, "");
  validPL->set<bool>("Lump Mass Matrix", false, "");
  validPL->sublist("Stratimikos", false, "");
  return validPL;
}

