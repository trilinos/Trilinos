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

//#include "Piro_NOXSolver.hpp"


template <typename Scalar>
Piro::NOXSolver<Scalar>::
NOXSolver(Teuchos::RCP<Teuchos::ParameterList> appParams_,
	  Teuchos::RCP< Thyra::ModelEvaluatorDefaultBase<Scalar> > model_) :
  appParams(appParams_),
  model(model_)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out = Teuchos::VerboseObjectBase::getDefaultOStream();

  RCP<Teuchos::ParameterList> noxParams =
	rcp(&(appParams->sublist("NOX")),false);

//   string jacobianSource = appParams->get("Jacobian Operator", "Have Jacobian");

//   if (jacobianSource == "Matrix-Free") {
//     TEUCHOS_TEST_FOR_EXCEPTION(jacobianSource == "Matrix-Free", std::logic_error,
//        "MATRIX_free not yet implemented for Piro Thyra");
//     model = rcp(new Piro::Thyra::MatrixFreeDecorator(model));
//   }

  // Grab some modelEval stuff from underlying model
  num_p = model->createInArgs().Np();
  num_g = model->createOutArgs().Ng();

  solver = rcp(new Thyra::NOXNonlinearSolver);
  solver->setParameterList(noxParams);
  solver->setModel(model);
}

template <typename Scalar>
Piro::NOXSolver<Scalar>::~NOXSolver()
{
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::NOXSolver<Scalar>::get_p_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(l > num_p || l < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::NOXSolver::get_p_map():  " <<
                     "Invalid parameter index l = " <<
                     l << std::endl);
  if      (l < num_p) return model->get_p_space(l);
  else return model->get_x_space(); // l == num_p
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::NOXSolver<Scalar>::get_g_space(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(j > num_g || j < 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl <<
                     "Error in Piro::NOXSolver::get_g_map():  " <<
                     "Invalid response index j = " <<
                     j << std::endl);

  if      (j < num_g) return model->get_g_space(j);
  else return model->get_x_space(); // j == num_g
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::NOXSolver<Scalar>::getNominalValues() const
{
  return model->getNominalValues();
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> Piro::NOXSolver<Scalar>::createInArgs() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p+1);
  return inArgs;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar> Piro::NOXSolver<Scalar>::createOutArgsImpl() const
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());

  // Ng is 1 bigger then model-Ng so that the solution vector can be an outarg
  outArgs.set_Np_Ng(num_p, num_g+1);

  Thyra::ModelEvaluatorBase::OutArgs<Scalar> model_outargs = model->createOutArgs();
  for (int i=0; i<num_g; i++)
    for (int j=0; j<num_p; j++)
      if (!model_outargs.supports( Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, i, j).none())
        outArgs.setSupports( Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, i, j,
             Thyra::ModelEvaluatorBase::DerivativeSupport( Thyra::ModelEvaluatorBase::DERIV_MV_BY_COL));

  return outArgs;
}

template <typename Scalar>
void Piro::NOXSolver<Scalar>::evalModelImpl(
       const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
       const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  //*out << "In eval Modle " << endl;

  // Parse InArgs

  // grab initial guess if user has provided one
  RCP<const Thyra::VectorBase<Scalar> > initial_guess_in;
  initial_guess_in = inArgs.get_p(num_p);

  RCP<const Thyra::VectorBase<Scalar> > p_in;
  if (num_p > 0) p_in = inArgs.get_p(0);

  // Parse OutArgs: always 1 extra
  RCP< Thyra::VectorBase<Scalar> > g_out; 
  if (num_g > 0) g_out = outArgs.get_g(0); 
  RCP< Thyra::VectorBase<Scalar> > gx_out = outArgs.get_g(num_g); 

  // Parse out-args for sensitivity calculation

  ::Thyra::SolveCriteria<double> solve_criteria;
  ::Thyra::SolveStatus<double> solve_status;

  RCP< ::Thyra::VectorBase<double> > initial_guess;
  if(initial_guess_in!=Teuchos::null)
    initial_guess = initial_guess_in->clone_v();
  else
    initial_guess = model->getNominalValues().get_x()->clone_v();

  solve_status = solver->solve(initial_guess.get(), &solve_criteria, NULL);

//   if (solve_status.solveStatus == ::Thyra::SOLVE_STATUS_CONVERGED)
//     std::cout << "Test passed!" << std::endl;

  // return the final solution as an additional g-vector, if requested
  RCP<const Thyra::VectorBase<Scalar> > finalSolution = solver->get_current_x();

  if (gx_out != Teuchos::null)  Thyra::copy(*finalSolution, gx_out.ptr());

  if (g_out != Teuchos::null) {
     // As post-processing step, calc responses at final solution
     Thyra::ModelEvaluatorBase::InArgs<Scalar>  model_inargs = model->createInArgs();
     Thyra::ModelEvaluatorBase::OutArgs<Scalar> model_outargs = model->createOutArgs();
     model_inargs.set_x(finalSolution);
     if (num_p > 0)  model_inargs.set_p(0, p_in);
     if (g_out != Teuchos::null) {
       Thyra::put_scalar(0.0,g_out.ptr());
       model_outargs.set_g(0, g_out);
     }

     model->evalModel(model_inargs, model_outargs);
  }

/*********************  NEED TO CONVERT TO THYRA *******************
  RCP< Thyra::MultiVectorBase<Scalar> > dgdp_out;
  if (num_p>0 && num_g>0)
    dgdp_out = outArgs.get_DgDp(0,0).getMultiVector();

  if (dgdp_out == Teuchos::null) {

     Teuchos::RCP<Epetra_MultiVector> dgdx
          = Teuchos::rcp(new Epetra_MultiVector(finalSolution->Map(),
                                                   dgdp_out->GlobalLength()));
     Teuchos::Array<int> p_indexes =
       outArgs.get_DgDp(0,0).getDerivativeMultiVector().getParamIndexes();
 
     EpetraExt::ModelEvaluator::DerivativeMultiVector dmv_dgdp(dgdp_out,
                                                               DERIV_MV_BY_COL,
                                                               p_indexes);
 
     EpetraExt::ModelEvaluator::InArgs model_inargs = model->createInArgs();
     EpetraExt::ModelEvaluator::OutArgs model_outargs = model->createOutArgs();
     model_inargs.set_x(finalSolution);
     model_inargs.set_p(0, p_in);

     if (g_out != Teuchos::null) {
       g_out->PutScalar(0.0);
       model_outargs.set_g(0, g_out);
     }
     model_outargs.set_DgDp(0,0,dmv_dgdp);
     model_outargs.set_DgDx(0,dgdx);

     model->evalModel(model_inargs, model_outargs);

 
     // (3) Calculate dg/dp = dg/dx*dx/dp + dg/dp
     // This may be the transpose of what we want since we specified
     // we want dg/dp by column in createOutArgs().
     // In this case just interchange the order of dgdx and dxdp
     // We should really probably check what the underlying ME does

     if (Teuchos::VERB_MEDIUM <= solnVerbLevel) cout << " dgdx \n" << *dgdx << endl;
     if (Teuchos::VERB_MEDIUM <= solnVerbLevel) cout << " dxdp \n" << *dxdp << endl;

     dgdp_out->Multiply('T', 'N', 1.0, *dgdx, *dxdp, 1.0);
   }
*********************/
}

