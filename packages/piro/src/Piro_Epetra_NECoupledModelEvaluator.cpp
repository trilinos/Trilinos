/*
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
*/

#include "Piro_Epetra_NECoupledModelEvaluator.hpp"
#include "Piro_Epetra_Factory.hpp"
#include "Piro_Epetra_StokhosSolver.hpp"

#include "Epetra_LocalMap.h"

#include "Teuchos_Assert.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "Stokhos_Epetra.hpp"
#include "Stokhos_StieltjesGramSchmidtBuilder.hpp"
#include "EpetraExt_MultiComm.h"

Piro::Epetra::NECoupledModelEvaluator::
NECoupledModelEvaluator(
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& modelA_, 
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& modelB_,
  const Teuchos::RCP<Teuchos::ParameterList>& piroParamsA_,
  const Teuchos::RCP<Teuchos::ParameterList>& piroParamsB_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_,
  const Teuchos::RCP<const Epetra_Comm>& comm_):
  modelA(modelA_),
  modelB(modelB_),
  piroParamsA(piroParamsA_),
  piroParamsB(piroParamsB_),
  params(params_),
  comm(comm_)
{
  // Setup VerboseObject
  Teuchos::readVerboseObjectSublist(params.get(), this);

  // Create solvers for models A and B
  bool stochastic = params->get("Stochastic", false);
  if (stochastic) {
    sgSolverA = Teuchos::rcp(new Piro::Epetra::StokhosSolver(piroParamsA, 
							     comm));
    sgSolverB = Teuchos::rcp(new Piro::Epetra::StokhosSolver(piroParamsB, 
							     comm));
    sgSolverA->setup(modelA);
    sgSolverB->setup(modelB);
    solverA = sgSolverA;
    solverB = sgSolverB;
  }
  else {
    solverA = Piro::Epetra::Factory::createSolver(piroParamsA, modelA);
    solverB = Piro::Epetra::Factory::createSolver(piroParamsB, modelB);
  }

  // Get connectivity information
  pIndexA = params->get("Model A Parameter Index", 0);
  pIndexB = params->get("Model B Parameter Index", 0);
  gIndexA = params->get("Model A Response Index", 0);
  gIndexB = params->get("Model B Response Index", 0);

  // Get number of parameter and response vectors
  EpetraExt::ModelEvaluator::InArgs solverA_inargs = solverA->createInArgs();
  EpetraExt::ModelEvaluator::InArgs solverB_inargs = solverB->createInArgs();
  EpetraExt::ModelEvaluator::OutArgs solverA_outargs = solverA->createOutArgs();
  EpetraExt::ModelEvaluator::OutArgs solverB_outargs = solverB->createOutArgs();
  num_params_A = solverA_inargs.Np();
  num_params_B = solverB_inargs.Np();
  num_responses_A = solverA_outargs.Ng();
  num_responses_B = solverB_outargs.Ng();

  // Building indexing maps between coupled system parameters/responses and
  // individual components
  num_params_total = num_params_A + num_params_B - 2;
  param_map.resize(num_params_total);
  for (int i=0; i<pIndexA; i++)
    param_map[i] = i;
  for (int i=pIndexA+1; i<num_params_A; i++)
    param_map[i-1] = i;
  for (int i=0; i<pIndexB; i++)
    param_map[num_params_A-1+i] = i;
  for (int i=pIndexB+1; i<num_params_B; i++)
    param_map[num_params_A-2+i] = i;

  num_responses_total = num_responses_A + num_responses_B - 2;
  response_map.resize(num_responses_total);
  for (int i=0; i<pIndexA; i++)
    response_map[i] = i;
  for (int i=pIndexA+1; i<num_responses_A; i++)
    response_map[i-1] = i;
  for (int i=0; i<pIndexB; i++)
    response_map[num_responses_A-1+i] = i;
  for (int i=pIndexB+1; i<num_responses_B; i++)
    response_map[num_responses_A-2+i] = i;

  //
  // The network equations look like:
  //    p_1 - g_2(x_2,p_2) = 0 s.t. f_1(x_1,p_1) = 0
  //    p_2 - g_1(x_1,p_1) = 0 s.t. f_2(x_2,p_2) = 0
  //
  // We define x = [p_1; p_2] and f = [ p_1 - g_2(x_2,p_2); p_2 - g_1(x_1,p_1)]
  //

  // Build x map, which is the product of the pIndexA and pIndexB parameter maps
  // For the time being, we will assume local maps, in the future we need to
  // build proper product maps
  p_map_A = solverA->get_p_map(pIndexA);
  p_map_B = solverB->get_p_map(pIndexB);
  n_p_A = p_map_A->NumGlobalElements();
  n_p_B = p_map_B->NumGlobalElements();
  int nx = n_p_A + n_p_B;
  x_map = Teuchos::rcp(new Epetra_Map(nx, 0, *comm));
  x_overlap_map = Teuchos::rcp(new Epetra_LocalMap(nx, 0, *comm));
  x_importer = Teuchos::rcp(new Epetra_Import(*x_overlap_map, *x_map));
  x_overlap = Teuchos::rcp(new Epetra_Vector(*x_overlap_map));

  // Build f map, which is the product of the gIndexA and gIndexB response maps
  // For the time being, we will assume local maps, in the future we need to
  // build proper product maps
  g_map_A = solverA->get_g_map(gIndexA);
  g_map_B = solverB->get_g_map(gIndexB);
  n_g_A = g_map_A->NumGlobalElements();
  n_g_B = g_map_B->NumGlobalElements();
  int nf = n_g_A + n_g_B;
  f_map = Teuchos::rcp(new Epetra_Map(nf, 0, *comm));
  f_overlap_map = Teuchos::rcp(new Epetra_LocalMap(nf, 0, *comm));
  f_exporter = Teuchos::rcp(new Epetra_Export(*f_overlap_map, *f_map));
  f_overlap = Teuchos::rcp(new Epetra_Vector(*f_overlap_map));

  // Do some consistency checking
  TEUCHOS_TEST_FOR_EXCEPTION(
    !p_map_A->SameAs(*g_map_B), std::logic_error,
    "Model A parameter map for index " << pIndexA << " must be the same " <<
    "map as model B response map for index " << gIndexB << "!");
  TEUCHOS_TEST_FOR_EXCEPTION(
    !p_map_B->SameAs(*g_map_A), std::logic_error,
    "Model B parameter map for index " << pIndexB << " must be the same " <<
    "map as model A response map for index " << gIndexA << "!");

  DerivativeSupport dsA = 
    solverA_outargs.supports(OUT_ARG_DgDp, gIndexA, pIndexA);
  DerivativeSupport dsB = 
    solverB_outargs.supports(OUT_ARG_DgDp, gIndexB, pIndexB);
  supports_W = 
    (dsA.supports(DERIV_MV_BY_COL) || dsA.supports(DERIV_TRANS_MV_BY_ROW)) && 
    (dsB.supports(DERIV_MV_BY_COL) || dsB.supports(DERIV_TRANS_MV_BY_ROW));

  supports_x_sg = 
    solverA_inargs.supports(IN_ARG_p_sg, pIndexA) &&
    solverB_inargs.supports(IN_ARG_p_sg, pIndexB);
  supports_f_sg = 
    solverA_outargs.supports(OUT_ARG_g_sg, gIndexA) &&
    solverB_outargs.supports(OUT_ARG_g_sg, gIndexB);
  DerivativeSupport dsA_sg = 
    solverA_outargs.supports(OUT_ARG_DgDp_sg, gIndexA, pIndexA);
  DerivativeSupport dsB_sg = 
    solverB_outargs.supports(OUT_ARG_DgDp_sg, gIndexB, pIndexB);
  supports_W_sg = 
    (dsA_sg.supports(DERIV_MV_BY_COL) || 
     dsA_sg.supports(DERIV_TRANS_MV_BY_ROW)) && 
    (dsB_sg.supports(DERIV_MV_BY_COL) || 
     dsB_sg.supports(DERIV_TRANS_MV_BY_ROW));

  // Build the Jacobian graph, which looks like
  //   [     I      dg_2/dp_2 ]
  //   [ dg_1/dp_1     I      ]
  if (supports_W || supports_W_sg) {
    W_graph = Teuchos::rcp(new Epetra_CrsGraph(Copy, *f_map, nx));
    for (int i=0; i<f_map->NumMyElements(); i++) {
      int row = f_map->GID(i);
      // Diagonal part
      int index = row;
      W_graph->InsertGlobalIndices(row, 1, &index);
      if (row < n_p_A) {
	// dg_2/dp_2 part
	for (int j=0; j<n_p_B; j++) {
	  index = n_p_A + j;
	  W_graph->InsertGlobalIndices(row, 1, &index);
	}
      }
      else {
	// dg_1/dp_1 part
	for (int j=0; j<n_p_A; j++) {
	  index = j;
	  W_graph->InsertGlobalIndices(row, 1, &index);
	}
      }
    }
    W_graph->FillComplete();

    W_overlap_graph = 
      Teuchos::rcp(new Epetra_CrsGraph(Copy, *f_overlap_map, nx));
    for (int i=0; i<n_p_A; i++) {
      int row = i;

      // Diagonal part
      int index = row;
      W_overlap_graph->InsertGlobalIndices(row, 1, &index);
      
      // dg_2/dp_2 part
      for (int j=0; j<n_p_B; j++) {
	index = n_p_A + j;
	W_overlap_graph->InsertGlobalIndices(row, 1, &index);
      }
    }
    for (int i=0; i<n_p_B; i++) {
      int row = i+n_p_B;

      // Diagonal part
      int index = row;
      W_overlap_graph->InsertGlobalIndices(row, 1, &index);

      // dg_1/dp_1 part
      for (int j=0; j<n_p_A; j++) {
	index = j;
	W_overlap_graph->InsertGlobalIndices(row, 1, &index);
      }
    }
    W_overlap_graph->FillComplete();
    W_overlap = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *W_overlap_graph));
  }

  // Build initial guess
  Epetra_Vector x_init_overlap(*x_overlap_map);
  Teuchos::RCP<const Epetra_Vector> p_init_A = solverA->get_p_init(pIndexA);
  Teuchos::RCP<const Epetra_Vector> p_init_B = solverB->get_p_init(pIndexB);
  for (int i=0; i<n_p_A; i++)
    x_init_overlap[i] = (*p_init_A)[i];
  for (int i=0; i<n_p_B; i++)
    x_init_overlap[n_p_A + i] = (*p_init_B)[i];
  x_init = Teuchos::rcp(new Epetra_Vector(*x_map));
  x_init->Export(x_init_overlap, *x_importer, Insert);

  // Create storage for parameters, responses, and derivatives
  p_A = Teuchos::rcp(new Epetra_Vector(*p_map_A));
  p_B = Teuchos::rcp(new Epetra_Vector(*p_map_B));
  g_A = Teuchos::rcp(new Epetra_Vector(*g_map_A));
  g_B = Teuchos::rcp(new Epetra_Vector(*g_map_B));
  if (supports_W || supports_W_sg) {
    if (dsA.supports(DERIV_MV_BY_COL)) {
      dgdp_A_layout = DERIV_MV_BY_COL;
      dgdp_A = Teuchos::rcp(new Epetra_MultiVector(*g_map_A, n_p_A));
    }
    else {
      dgdp_A_layout = DERIV_TRANS_MV_BY_ROW;
      dgdp_A = Teuchos::rcp(new Epetra_MultiVector(*p_map_A, n_g_A));
    }
    if (dsB.supports(DERIV_MV_BY_COL)) {
      dgdp_B_layout = DERIV_MV_BY_COL;
      dgdp_B = Teuchos::rcp(new Epetra_MultiVector(*g_map_B, n_p_B));
    }
    else {
      dgdp_B_layout = DERIV_TRANS_MV_BY_ROW;
      dgdp_B = Teuchos::rcp(new Epetra_MultiVector(*p_map_B, n_g_B));
    }
  }

  reduce_dimension = params->get("Reduce Dimension", false);
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
Piro::Epetra::NECoupledModelEvaluator::
get_x_map() const
{
  return x_map;
}

Teuchos::RCP<const Epetra_Map>
Piro::Epetra::NECoupledModelEvaluator::
get_f_map() const
{
  return f_map;
}

Teuchos::RCP<const Epetra_Vector>
Piro::Epetra::NECoupledModelEvaluator::
get_x_init() const
{
  return x_init;
}

Teuchos::RCP<const Epetra_Map>
Piro::Epetra::NECoupledModelEvaluator::
get_p_map(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    j >= num_params_total || j < 0, Teuchos::Exceptions::InvalidParameter,
    std::endl <<
    "Error in Piro::Epetra::NECoupledModelEvaluator::get_p_map():  " <<
    "Invalid parameter index j = " << j << std::endl);

  if (j < num_params_A-1)
    return solverA->get_p_map(param_map[j]);
  return solverB->get_p_map(param_map[j]);
}

Teuchos::RCP<const Epetra_Map>
Piro::Epetra::NECoupledModelEvaluator::
get_g_map(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    j >= num_responses_total || j < 0, Teuchos::Exceptions::InvalidParameter,
    std::endl <<
    "Error in Piro::Epetra::NECoupledModelEvaluator::get_g_map():  " <<
    "Invalid response index j = " << j << std::endl);

  if (j < num_responses_A-1)
    return solverA->get_g_map(response_map[j]);
  return solverB->get_g_map(response_map[j]);
}

Teuchos::RCP<const Teuchos::Array<std::string> > 
Piro::Epetra::NECoupledModelEvaluator::
get_p_names(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    j >= num_params_total || j < 0, Teuchos::Exceptions::InvalidParameter,
    std::endl <<
    "Error in Piro::Epetra::NECoupledModelEvaluator::get_p_names():  " <<
    "Invalid parameter index j = " << j << std::endl);

  if (j < num_params_A-1)
    return solverA->get_p_names(param_map[j]);
  return solverB->get_p_names(param_map[j]);
}

Teuchos::RCP<const Epetra_Vector>
Piro::Epetra::NECoupledModelEvaluator::
get_p_init(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    j >= num_params_total || j < 0, Teuchos::Exceptions::InvalidParameter,
    std::endl <<
    "Error in Piro::Epetra::NECoupledModelEvaluator::get_p_init():  " <<
    "Invalid parameter index j = " << j << std::endl);

  if (j < num_params_A-1)
    return solverA->get_p_init(param_map[j]);
  return solverB->get_p_init(param_map[j]);
}

Teuchos::RCP<Epetra_Operator>
Piro::Epetra::NECoupledModelEvaluator::
create_W() const
{
  Teuchos::RCP<Epetra_CrsMatrix> mat =
    Teuchos::rcp(new Epetra_CrsMatrix(Copy, *W_graph));
  mat->FillComplete();
  return mat;
}

EpetraExt::ModelEvaluator::InArgs
Piro::Epetra::NECoupledModelEvaluator::
createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());

  EpetraExt::ModelEvaluator::InArgs solverA_inargs = solverA->createInArgs();
  EpetraExt::ModelEvaluator::InArgs solverB_inargs = solverB->createInArgs();

  // Deterministic InArgs
  inArgs.setSupports(IN_ARG_x, true);
  inArgs.set_Np(num_params_total);
  
  // Stochastic InArgs
  inArgs.setSupports(IN_ARG_x_sg, supports_x_sg);
  for (int i=0; i<num_params_A-1; i++)
    inArgs.setSupports(IN_ARG_p_sg, i, 
		       solverA_inargs.supports(IN_ARG_p_sg, param_map[i])); 
  for (int i=0; i<num_params_B-1; i++)
    inArgs.setSupports(IN_ARG_p_sg, num_params_A-1+i, 
		       solverB_inargs.supports(IN_ARG_p_sg, 
					      param_map[num_params_A-1+i])); 
  inArgs.setSupports(IN_ARG_sg_basis,true);
  inArgs.setSupports(IN_ARG_sg_quadrature,true);
  inArgs.setSupports(IN_ARG_sg_expansion,true);

  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
Piro::Epetra::NECoupledModelEvaluator::
createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());

  EpetraExt::ModelEvaluator::OutArgs solverA_outargs = solverA->createOutArgs();
  EpetraExt::ModelEvaluator::OutArgs solverB_outargs = solverB->createOutArgs();

  // Deterministic OutArgs
  outArgs.setSupports(OUT_ARG_f, true);
  outArgs.setSupports(OUT_ARG_W, supports_W);
  outArgs.set_W_properties(
    DerivativeProperties(DERIV_LINEARITY_NONCONST, DERIV_RANK_FULL, true));
  outArgs.set_Np_Ng(num_params_total, num_responses_total);
  for (int i=0; i<num_params_A-1; i++)
    outArgs.setSupports(
      OUT_ARG_DfDp, i,
      solverA_outargs.supports(OUT_ARG_DgDp, gIndexA, param_map[i])); 
  for (int i=0; i<num_params_B-1; i++)
    outArgs.setSupports(
      OUT_ARG_DfDp, num_params_A-1+i, 
      solverB_outargs.supports(OUT_ARG_DgDp, gIndexB, 
			      param_map[num_params_A-1+i]));
  for (int i=0; i<num_responses_A-1; i++) {
    outArgs.setSupports(OUT_ARG_DgDx, i, 
			solverA_outargs.supports(OUT_ARG_DgDx, 
						 response_map[i]));
    for (int j=0; j<num_params_A-1; j++)
      outArgs.setSupports(OUT_ARG_DgDp, i, j, 
			  solverA_outargs.supports(OUT_ARG_DgDp, 
						   response_map[i],
						   param_map[j])); 
  }
  for (int i=0; i<num_responses_B-1; i++) {
    outArgs.setSupports(
      OUT_ARG_DgDx, num_responses_A-1+i, 
      solverA_outargs.supports(OUT_ARG_DgDx, 
			      response_map[num_responses_A-1+i]));
    for (int j=0; j<num_params_B-1; j++)
      outArgs.setSupports(
	OUT_ARG_DgDp, num_responses_A-1+i, num_params_A-1+j, 
	solverB_outargs.supports(OUT_ARG_DgDp, 
				 response_map[num_responses_A-1+i],
				 param_map[num_params_A-1+j]));
  }

  // Stochastic OutArgs
  outArgs.setSupports(OUT_ARG_f_sg, supports_f_sg);
  outArgs.setSupports(OUT_ARG_W_sg, supports_W_sg);
  for (int i=0; i<num_responses_A-1; i++)
    outArgs.setSupports(
      OUT_ARG_g_sg, i, solverA_outargs.supports(OUT_ARG_g_sg, response_map[i]));
  for (int i=0; i<num_responses_B-1; i++)
    outArgs.setSupports(
      OUT_ARG_g_sg, num_responses_A-1+i, 
      solverB_outargs.supports(OUT_ARG_g_sg, 
			      response_map[num_responses_A-1+i]));
  for (int i=0; i<num_params_A-1; i++)
    outArgs.setSupports(
      OUT_ARG_DfDp_sg, i, 
      solverA_outargs.supports(OUT_ARG_DgDp_sg, gIndexA, param_map[i])); 
  for (int i=0; i<num_params_B-1; i++)
    outArgs.setSupports(
      OUT_ARG_DfDp_sg, num_params_A-1+i, 
      solverB_outargs.supports(OUT_ARG_DgDp_sg, gIndexB, 
			      param_map[num_params_A-1+i]));
  for (int i=0; i<num_responses_A-1; i++) {
    outArgs.setSupports(OUT_ARG_DgDx_sg, i, 
			solverA_outargs.supports(OUT_ARG_DgDx_sg, 
						response_map[i]));
    for (int j=0; j<num_params_A-1; j++)
      outArgs.setSupports(OUT_ARG_DgDp_sg, i, j, 
			  solverA_outargs.supports(OUT_ARG_DgDp_sg, 
						  response_map[i],
						  param_map[j]));
  }
  for (int i=0; i<num_responses_B-1; i++) {
    outArgs.setSupports(
      OUT_ARG_DgDx_sg, num_responses_A-1+i, 
      solverA_outargs.supports(OUT_ARG_DgDx_sg, 
			      response_map[num_responses_A-1+i]));
    for (int j=0; j<num_params_B-1; j++)
      outArgs.setSupports(
	OUT_ARG_DgDp_sg, num_responses_A-1+i, num_params_A-1+j, 
	solverB_outargs.supports(OUT_ARG_DgDp_sg, 
				 response_map[num_responses_A-1+i],
				 param_map[num_params_A-1+j]));
  }
  
  return outArgs;
}

void 
Piro::Epetra::NECoupledModelEvaluator::
evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();

  EpetraExt::ModelEvaluator::InArgs solverA_inargs = solverA->createInArgs();
  EpetraExt::ModelEvaluator::InArgs solverB_inargs = solverB->createInArgs();
  EpetraExt::ModelEvaluator::OutArgs solverA_outargs = solverA->createOutArgs();
  EpetraExt::ModelEvaluator::OutArgs solverB_outargs = solverB->createOutArgs();

  //
  // Deterministic calculation
  //
  bool do_deterministic = false;
  Teuchos::RCP<const Epetra_Vector> x = inArgs.get_x();
  if (x != Teuchos::null) {
    do_deterministic = true;

    // p
    x_overlap->Import(*x, *x_importer, Insert);
    for (int i=0; i<n_p_A; i++)
      (*p_A)[i] = (*x_overlap)[i];
    for (int i=0; i<n_p_B; i++)
      (*p_B)[i] = (*x_overlap)[n_p_A + i];
    solverA_inargs.set_p(pIndexA, p_A);
    solverB_inargs.set_p(pIndexB, p_B);
    for (int i=0; i<num_params_A-1; i++)
      solverA_inargs.set_p(param_map[i], inArgs.get_p(i));
    for (int i=0; i<num_params_B-1; i++)
      solverB_inargs.set_p(param_map[num_params_A-1+i], 
			   inArgs.get_p(num_params_A-1+i));
    
    // f
    Teuchos::RCP<Epetra_Vector> f = outArgs.get_f();
    if (f != Teuchos::null) {
      solverA_outargs.set_g(gIndexA, g_A);
      solverB_outargs.set_g(gIndexB, g_B);
    }

    // W
    Teuchos::RCP<Epetra_Operator> W = outArgs.get_W();
    if (W != Teuchos::null) {
      Derivative dgdp_A_deriv(dgdp_A, dgdp_A_layout);
      Derivative dgdp_B_deriv(dgdp_B, dgdp_B_layout);
      solverA_outargs.set_DgDp(gIndexA, pIndexA, dgdp_A_deriv);
      solverB_outargs.set_DgDp(gIndexB, pIndexB, dgdp_B_deriv);
    }

    // g
    for (int i=0; i<num_responses_A-1; i++)
      solverA_outargs.set_g(response_map[i], outArgs.get_g(i));
    for (int i=0; i<num_responses_B-1; i++)
      solverB_outargs.set_g(response_map[i+num_responses_A-1], 
			    outArgs.get_g(i+num_responses_A-1));

    // dg/dx
    for (int i=0; i<num_responses_A-1; i++)
      if (!solverA_outargs.supports(OUT_ARG_DgDx, response_map[i]).none())
	solverA_outargs.set_DgDx(response_map[i], outArgs.get_DgDx(i));
    for (int i=0; i<num_responses_B-1; i++)
      if (!solverB_outargs.supports(OUT_ARG_DgDx, 
				    response_map[i+num_responses_A-1]).none())
	solverB_outargs.set_DgDx(response_map[i+num_responses_A-1], 
				 outArgs.get_DgDx(i+num_responses_A-1));

    // dg/dp
    for (int i=0; i<num_responses_A-1; i++)
      for (int j=0; j<num_params_A-1; j++)
	if (!solverA_outargs.supports(OUT_ARG_DgDp, response_map[i], 
				      param_map[j]).none())
	  solverA_outargs.set_DgDp(response_map[i], param_map[j], 
				   outArgs.get_DgDp(i,j));
    for (int i=0; i<num_responses_B-1; i++)
      for (int j=0; j<num_params_B-1; j++)
	if (!solverB_outargs.supports(OUT_ARG_DgDp,
				      param_map[j+num_params_A-1],
				      response_map[i+num_responses_A-1]).none())
	  solverB_outargs.set_DgDp(response_map[i+num_responses_A-1], 
				   param_map[j+num_params_A-1],
				   outArgs.get_DgDp(i+num_responses_A-1,
						    j+num_params_A-1));
  }

  //
  // Stochastic calculation
  //
  bool do_stochastic = false;
  InArgs::sg_const_vector_t x_sg;
  if (supports_x_sg)
    x_sg = inArgs.get_x_sg();
  if (x_sg != Teuchos::null) {
    do_stochastic = true;

    // SG data
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis = 
      inArgs.get_sg_basis();
    if (solverA_inargs.supports(IN_ARG_sg_basis))
      solverA_inargs.set_sg_basis(basis);
    if (solverA_inargs.supports(IN_ARG_sg_quadrature))
      solverA_inargs.set_sg_quadrature(inArgs.get_sg_quadrature());
    if (solverA_inargs.supports(IN_ARG_sg_expansion))
      solverA_inargs.set_sg_expansion(inArgs.get_sg_expansion());
    if (solverB_inargs.supports(IN_ARG_sg_basis))
      solverB_inargs.set_sg_basis(basis);
    if (solverB_inargs.supports(IN_ARG_sg_quadrature))
      solverB_inargs.set_sg_quadrature(inArgs.get_sg_quadrature());
    if (solverB_inargs.supports(IN_ARG_sg_expansion))
      solverB_inargs.set_sg_expansion(inArgs.get_sg_expansion());

    if (p_A_sg == Teuchos::null) {
      Teuchos::RCP<const EpetraExt::MultiComm> multiComm;
      multiComm = x_sg->productComm();
      sg_overlap_map =
	Teuchos::rcp(new Epetra_LocalMap(basis->size(), 0, 
					 multiComm->TimeDomainComm()));
      p_A_sg = 
	Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		       basis, sg_overlap_map, p_map_A, multiComm));
      p_B_sg = 
	Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		       basis, sg_overlap_map, p_map_B, multiComm));
      if (supports_f_sg) {
	g_A_sg = 
	  Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
			 basis, sg_overlap_map, g_map_A, multiComm));
	g_B_sg = 
	  Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
			 basis, sg_overlap_map, g_map_B, multiComm));
      }
      if (supports_W_sg) {
	DerivativeSupport dsA_sg = 
	  solverA_outargs.supports(OUT_ARG_DgDp_sg, gIndexA, pIndexA);
	DerivativeSupport dsB_sg = 
	  solverB_outargs.supports(OUT_ARG_DgDp_sg, gIndexB, pIndexB);
	if (dsA_sg.supports(DERIV_MV_BY_COL)) {
	  dgdp_A_sg_layout = DERIV_MV_BY_COL;
	  dgdp_A_sg = 
	    Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			   basis, sg_overlap_map, g_map_A, multiComm, n_p_A));
	}
	else {
	  dgdp_A_sg_layout = DERIV_TRANS_MV_BY_ROW;
	  dgdp_A_sg = 
	    Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			   basis, sg_overlap_map, p_map_A, multiComm, n_g_A));
	}
	if (dsB_sg.supports(DERIV_MV_BY_COL)) {
	  dgdp_B_sg_layout = DERIV_MV_BY_COL;
	  dgdp_B_sg = 
	    Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			   basis, sg_overlap_map, g_map_B, multiComm, n_p_B));
	}
	else {
	  dgdp_B_sg_layout = DERIV_TRANS_MV_BY_ROW;
	  dgdp_B_sg = 
	    Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			   basis, sg_overlap_map, p_map_B, multiComm, n_g_B));
	}
      }
    }
  
    // p_sg
    for (int block=0; block<x_sg->size(); block++) {
      x_overlap->Import((*x_sg)[block], *x_importer, Insert);
      for (int i=0; i<n_p_A; i++)
	(*p_A_sg)[block][i] = (*x_overlap)[i];
      for (int i=0; i<n_p_B; i++)
	(*p_B_sg)[block][i] = (*x_overlap)[n_p_A + i];
    }
    solverA_inargs.set_p_sg(pIndexA, p_A_sg);
    solverB_inargs.set_p_sg(pIndexB, p_B_sg);
    for (int i=0; i<num_params_A-1; i++)
      if (solverA_inargs.supports(IN_ARG_p_sg, param_map[i]))
	solverA_inargs.set_p_sg(param_map[i], inArgs.get_p_sg(i));
    for (int i=0; i<num_params_B-1; i++)
      if (solverB_inargs.supports(IN_ARG_p_sg, param_map[num_params_A-1+i]))
	solverB_inargs.set_p_sg(param_map[num_params_A-1+i], 
				inArgs.get_p_sg(num_params_A-1+i));
      
    // f_sg
    if (supports_f_sg) {
      OutArgs::sg_vector_t f_sg = outArgs.get_f_sg();
      if (f_sg != Teuchos::null) {
	solverA_outargs.set_g_sg(gIndexA, g_A_sg);
	solverB_outargs.set_g_sg(gIndexB, g_B_sg);
      }
    }

    // W_sg
    if (supports_W_sg) {
      OutArgs::sg_operator_t W_sg = outArgs.get_W_sg();
      if (W_sg != Teuchos::null) {
	solverA_outargs.set_DgDp_sg(gIndexA, pIndexA, dgdp_A_sg);
	solverB_outargs.set_DgDp_sg(gIndexB, pIndexB, dgdp_B_sg);
      }
    }

    // g_sg
    for (int i=0; i<num_responses_A-1; i++)
      if (solverA_outargs.supports(OUT_ARG_g_sg, response_map[i]))
	solverA_outargs.set_g_sg(response_map[i], outArgs.get_g_sg(i));
    for (int i=0; i<num_responses_B-1; i++)
      if (solverB_outargs.supports(OUT_ARG_g_sg, 
				   response_map[i+num_responses_A-1]))
	solverB_outargs.set_g_sg(response_map[i+num_responses_A-1], 
				 outArgs.get_g_sg(i+num_responses_A-1));

    // dg/dx_sg
    for (int i=0; i<num_responses_A-1; i++)
      if (!solverA_outargs.supports(OUT_ARG_DgDx_sg, response_map[i]).none())
	solverA_outargs.set_DgDx_sg(response_map[i], outArgs.get_DgDx_sg(i));
    for (int i=0; i<num_responses_B-1; i++)
      if (!solverB_outargs.supports(OUT_ARG_DgDx_sg, 
				    response_map[i+num_responses_A-1]).none())
	solverB_outargs.set_DgDx_sg(response_map[i+num_responses_A-1], 
				    outArgs.get_DgDx_sg(i+num_responses_A-1));

    // dg/dp_sg
    for (int i=0; i<num_responses_A-1; i++)
      for (int j=0; j<num_params_A-1; j++)
	if (!solverA_outargs.supports(OUT_ARG_DgDp_sg, response_map[i], 
				      param_map[j]).none())
	  solverA_outargs.set_DgDp_sg(response_map[i], param_map[j], 
				      outArgs.get_DgDp_sg(i,j));
    for (int i=0; i<num_responses_B-1; i++)
      for (int j=0; j<num_params_B-1; j++)
	if (!solverB_outargs.supports(OUT_ARG_DgDp_sg,
				      param_map[j+num_params_A-1],
				      response_map[i+num_responses_A-1]).none())
	  solverB_outargs.set_DgDp_sg(response_map[i+num_responses_A-1], 
				      param_map[j+num_params_A-1],
				      outArgs.get_DgDp_sg(i+num_responses_A-1,
							  j+num_params_A-1));
  }

  // Stochastic dimension reduction
  EpetraExt::ModelEvaluator::InArgs solverA_inargs_red;
  EpetraExt::ModelEvaluator::InArgs solverB_inargs_red;
  EpetraExt::ModelEvaluator::OutArgs solverA_outargs_red;
  EpetraExt::ModelEvaluator::OutArgs solverB_outargs_red;
  Teuchos::RCP<EpetraExt::ModelEvaluator> solverA_red;
  Teuchos::RCP<EpetraExt::ModelEvaluator> solverB_red;
  Teuchos::RCP<Teuchos::ParameterList> piroParamsA_red;
  Teuchos::RCP<Teuchos::ParameterList> piroParamsB_red;
  do_dimension_reduction(inArgs, 
			 solverA_inargs, solverA_outargs,
			 modelA, solverA, piroParamsA,
			 solverA_inargs_red, solverA_outargs_red,
			 solverA_red, piroParamsA_red,
			 red_basis_vals_A);
  do_dimension_reduction(inArgs,
			 solverB_inargs, solverB_outargs,
			 modelB, solverB, piroParamsB,
			 solverB_inargs_red, solverB_outargs_red,
			 solverB_red, piroParamsB_red,
			 red_basis_vals_B);

  int proc = comm->MyPID();

  // Evaluate models
  {
  TEUCHOS_FUNC_TIME_MONITOR("NECoupledModelEvaluator -- Model A nonlinear elimination");
  if (verbLevel != Teuchos::VERB_NONE)
    *out << "Eliminating model A states...";
  solverA_red->evalModel(solverA_inargs_red, solverA_outargs_red);
  }

  {
  TEUCHOS_FUNC_TIME_MONITOR("NECoupledModelEvaluator -- Model B nonlinear elimination");
  if (verbLevel != Teuchos::VERB_NONE)
    *out << "Eliminating model B states...";
  solverB_red->evalModel(solverB_inargs_red, solverB_outargs_red);
  }

  // Project back to original stochastic bases
  do_dimension_projection(inArgs, solverA_inargs_red, solverA_outargs_red,
			  red_basis_vals_A, solverA_outargs);
  do_dimension_projection(inArgs, solverB_inargs_red, solverB_outargs_red,
			  red_basis_vals_B, solverB_outargs);

  if (do_deterministic) {

    // f
    Teuchos::RCP<Epetra_Vector> f = outArgs.get_f();
    if (f != Teuchos::null) {
      f_overlap->PutScalar(0.0);
      for (int i=0; i<n_p_A; i++)
	(*f_overlap)[i]       = (*p_A)[i] - (*g_B)[i];
      for (int i=0; i<n_p_B; i++)
	(*f_overlap)[i+n_p_A] = (*p_B)[i] - (*g_A)[i];
      f->Export(*f_overlap, *f_exporter, Insert);
    }

    // W
    Teuchos::RCP<Epetra_Operator> W = outArgs.get_W();
    if (W != Teuchos::null) {
      Teuchos::RCP<Epetra_CrsMatrix> W_crs = 
	Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W, true);
      W_overlap->PutScalar(0.0);
      int row, col;
      double val;
      for (int i=0; i<n_p_A; i++) {
	row = i; 

	// Diagonal part
	col = row; 
	val = 1.0;
	W_overlap->ReplaceGlobalValues(row, 1, &val, &col);

	// dg_2/dp_2 part
	for (int j=0; j<n_p_B; j++) {
	  col = n_p_A+j; 
	  if (dgdp_B_layout == DERIV_MV_BY_COL)
	    val = -(*dgdp_B)[j][i];
	  else
	    val = -(*dgdp_B)[i][j];
	  W_overlap->ReplaceGlobalValues(row, 1, &val, &col);
	}
      }
      for (int i=0; i<n_p_B; i++) {
	row = n_p_A + i; 

	// Diagonal part
	col = row; 
	val = 1.0;
	W_overlap->ReplaceGlobalValues(row, 1, &val, &col);

	// dg_1/dp_1 part
	for (int j=0; j<n_p_A; j++) {
	  col = j; 
	  if (dgdp_A_layout == DERIV_MV_BY_COL)
	    val = -(*dgdp_A)[j][i];
	  else
	    val = -(*dgdp_A)[i][j];
	  W_overlap->ReplaceGlobalValues(row, 1, &val, &col);
	}
      }
      W_crs->Export(*W_overlap, *f_exporter, Insert);
    }
  }

  if (do_stochastic) {

    // f_sg
    if (supports_f_sg) {
      OutArgs::sg_vector_t f_sg = outArgs.get_f_sg();
      if (f_sg != Teuchos::null) {
	for (int block=0; block<f_sg->size(); block++) {
	  f_overlap->PutScalar(0.0);
	  for (int i=0; i<n_p_A; i++)
	    (*f_overlap)[i]       = (*p_A_sg)[block][i] - (*g_B_sg)[block][i];
	  for (int i=0; i<n_p_B; i++)
	    (*f_overlap)[i+n_p_A] = (*p_B_sg)[block][i] - (*g_A_sg)[block][i];
	  (*f_sg)[block].Export(*f_overlap, *f_exporter, Insert);
	}
      }
    }

    // W_sg
    if (supports_W_sg) {
      OutArgs::sg_operator_t W_sg = outArgs.get_W_sg();
      if (W_sg != Teuchos::null) {
	for (int block=0; block<W_sg->size(); block++) {
	  Teuchos::RCP<Epetra_CrsMatrix> W_crs = 
	    Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_sg->getCoeffPtr(block), true);
	  W_overlap->PutScalar(0.0);
	  int row, col;
	  double val;
	  for (int i=0; i<n_p_A; i++) {
	    row = i; 

	    // Diagonal part
	    if (block == 0) {
	      col = row; 
	      val = 1.0;
	      W_overlap->ReplaceGlobalValues(row, 1, &val, &col);
	    }

	    // dg_2/dp_2 part
	    for (int j=0; j<n_p_B; j++) {
	      col = n_p_A+j; 
	      if (dgdp_B_layout == DERIV_MV_BY_COL)
		val = -(*dgdp_B_sg)[block][j][i];
	      else
		val = -(*dgdp_B_sg)[block][i][j];
	      W_overlap->ReplaceGlobalValues(row, 1, &val, &col);
	    }
	  }
	  for (int i=0; i<n_p_B; i++) {
	    row = n_p_A + i; 
	    
	    // Diagonal part
	    if (block == 0) {
	      col = row; 
	      val = 1.0;
	      W_overlap->ReplaceGlobalValues(row, 1, &val, &col);
	    }

	    // dg_1/dp_1 part
	    for (int j=0; j<n_p_A; j++) {
	      col = j; 
	      if (dgdp_A_layout == DERIV_MV_BY_COL)
		val = -(*dgdp_A_sg)[block][j][i];
	      else
		val = -(*dgdp_A_sg)[block][i][j];
	      W_overlap->ReplaceGlobalValues(row, 1, &val, &col);
	    }
	  }
	  W_crs->Export(*W_overlap, *f_exporter, Insert);
	}
      }
    }
  }
}

void
Piro::Epetra::NECoupledModelEvaluator:: 
do_dimension_reduction(
  const InArgs& inArgs,
  const InArgs& solver_inargs, 
  const OutArgs& solver_outargs,
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& model,
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& solver,
  const Teuchos::RCP<Teuchos::ParameterList>& solver_params,
  InArgs& reduced_inargs, 
  OutArgs& reduced_outargs,
  Teuchos::RCP<EpetraExt::ModelEvaluator>& reduced_solver,
  Teuchos::RCP<Teuchos::ParameterList>& reduced_params,
  Teuchos::RCP<const Teuchos::Array< Teuchos::Array<double> > >& red_basis_vals) const
{
  TEUCHOS_FUNC_TIME_MONITOR("NECoupledModelEvaluator -- dimension reduction");

  // First copy the in/out args to set everything we don't modify
  reduced_inargs = solver_inargs;
  reduced_outargs = solver_outargs;
  reduced_solver = solver;
  reduced_params = params;

  // Make sure there is something to do
  InArgs::sg_const_vector_t x_sg;
  if (supports_x_sg)
    x_sg = inArgs.get_x_sg();
  if (!reduce_dimension || x_sg == Teuchos::null)
    return;

  Teuchos::RCP<const Stokhos::ProductBasis<int,double> > basis = 
    Teuchos::rcp_dynamic_cast<const Stokhos::ProductBasis<int,double> >(
      inArgs.get_sg_basis(), true);
  Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad =
    inArgs.get_sg_quadrature();
  Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double,StorageType> > expansion
    = inArgs.get_sg_expansion();
 
  // Copy Epetra PCEs into Stokhos PCE objects
  int total_num_p = 0;
  for (int i=0; i<solver_inargs.Np(); i++) {
    if (solver_inargs.supports(IN_ARG_p_sg, i) &&
	solver_inargs.get_p_sg(i) != Teuchos::null) {
      InArgs::sg_const_vector_t p_sg = solver_inargs.get_p_sg(i);
      total_num_p += p_sg->coefficientMap()->NumMyElements();
    }
  }
  int sz = basis->size();
  Teuchos::Array< Stokhos::OrthogPolyApprox<int,double> > p_opa(total_num_p);
  int index = 0;
  for (int i=0; i<solver_inargs.Np(); i++) {
    if (solver_inargs.supports(IN_ARG_p_sg, i) &&
	solver_inargs.get_p_sg(i) != Teuchos::null) {
      InArgs::sg_const_vector_t p_sg = solver_inargs.get_p_sg(i);
      for (int k=0; k<p_sg->coefficientMap()->NumMyElements(); k++)
	p_opa[index+k].reset(basis);
      for (int j=0; j<sz; j++) {
	for (int k=0; k<(*p_sg)[j].MyLength(); k++)
	  p_opa[index+k][j] = (*p_sg)[j][k];
      }
      index += p_sg->coefficientMap()->NumMyElements();
    }
  }

  // Build Stieltjes basis, quadrature, and new PCEs
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > red_basis;
  Teuchos::RCP<const Stokhos::Quadrature<int,double> > red_quad;
  Teuchos::Array<Stokhos::OrthogPolyApprox<int,double> > red_pces;
  bool orthogonalize_bases = params->get("Orthogonalize Bases", false);
  int order = basis->order();
  int new_order = order;
  if (orthogonalize_bases) {
    Stokhos::StieltjesGramSchmidtBuilder<int,double> gs_builder(
      quad, p_opa, new_order, true, false);
    red_basis = gs_builder.getReducedBasis();
    red_quad = gs_builder.getReducedQuadrature();
    red_basis_vals = Teuchos::rcp(&(red_quad->getBasisAtQuadPoints()),false);
    gs_builder.computeReducedPCEs(p_opa, red_pces);
  }
  else {
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double > > >
      coordinate_bases = basis->getCoordinateBases();
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double > > >
      new_coordinate_bases(p_opa.size());
    Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk = 
      expansion->getTripleProduct();
    if (st_quad == Teuchos::null) {
      st_quad = quad;
      // st_quad =
      //   Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(
      // 		 basis, new_order+1));
    }
    for (int i=0; i<p_opa.size(); i++) {
      new_coordinate_bases[i] = Teuchos::rcp(
	new Stokhos::StieltjesPCEBasis<int,double>(
	  new_order, Teuchos::rcp(&(p_opa[i]),false), st_quad, 
	  false, false, true, Cijk));
    }
    Teuchos::RCP<const Stokhos::ProductBasis<int,double> > tensor_basis = 
      Teuchos::rcp(
	new Stokhos::CompletePolynomialBasis<int,double>(new_coordinate_bases)
	);
    red_basis = tensor_basis;
    if (red_basis->dimension() <= 3)
      red_quad = 
	Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(
		       tensor_basis));
    else
#ifdef HAVE_STOKHOS_DAKOTA
      red_quad = 
	Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(
		       tensor_basis, new_order));
#else
    red_quad = 
      Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(
		     tensor_basis));
#endif
    const Teuchos::Array< Teuchos::Array<double> >& points = 
      quad->getQuadPoints();
    const Teuchos::Array< Teuchos::Array<double> >& basis_vals = 
      quad->getBasisAtQuadPoints();
    int nqp = points.size();
    Teuchos::Array<double> p_opa_val(p_opa.size());
    Teuchos::RCP< Teuchos::Array< Teuchos::Array<double> > > ncred_basis_vals
      = Teuchos::rcp(new Teuchos::Array< Teuchos::Array<double> >(nqp));
    for (int i=0; i<nqp; i++) {
      for (int j=0; j<p_opa_val.size(); j++)
	p_opa_val[j] = p_opa[j].evaluate(points[i], basis_vals[i]);
      (*ncred_basis_vals)[i].resize(red_basis->size());
      red_basis->evaluateBases(p_opa_val, (*ncred_basis_vals)[i]);
    }
    red_basis_vals = ncred_basis_vals;
    red_pces.resize(p_opa.size());
    for (int k=0; k<p_opa.size(); k++) {
      red_pces[k].reset(red_basis);
      red_pces[k].term(k, 0) = p_opa[k].mean();
      red_pces[k].term(k, 1) = 1.0; 
    }
  }

  Teuchos::RCP<const EpetraExt::MultiComm> multiComm = x_sg->productComm();
    
  // Copy into Epetra objects
  int red_sz = red_basis->size();
  Teuchos::RCP<const Epetra_BlockMap> red_overlap_map =
    Teuchos::rcp(new Epetra_LocalMap(red_sz, 0, 
				     multiComm->TimeDomainComm()));

  // p_red
  index = 0;
  for (int i=0; i<solver_inargs.Np(); i++) {
    if (solver_inargs.supports(IN_ARG_p_sg, i) &&
	solver_inargs.get_p_sg(i) != Teuchos::null) {
      OutArgs::sg_vector_t p_red = 
	Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		       red_basis, red_overlap_map, solver->get_p_map(i), 
		       multiComm));
      for (int j=0; j<red_sz; j++)
	for (int k=0; k<(*p_red)[j].MyLength(); k++)
	  (*p_red)[j][k] = red_pces[index+k][j];
      index += p_red->coefficientMap()->NumMyElements();
      reduced_inargs.set_p_sg(i, p_red);
    }
  }
  
  for (int i=0; i<solver_outargs.Ng(); i++) {

    // g_red
    if (solver_outargs.supports(OUT_ARG_g_sg, i) &&
	solver_outargs.get_g_sg(i) != Teuchos::null) {
      OutArgs::sg_vector_t g_red = 
	Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		       red_basis, red_overlap_map, solver->get_g_map(i), 
		       multiComm));
      reduced_outargs.set_g_sg(i, g_red);
    }

    // dg/dx_red
    if (!solver_outargs.supports(OUT_ARG_DgDx_sg, i).none()) {
      Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> dgdx_sg = 
	solver_outargs.get_DgDx_sg(i).getMultiVector();
      if (dgdx_sg != Teuchos::null) {
	Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> dgdx_red = 
	  Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			 red_basis, red_overlap_map, 
			 dgdx_sg->coefficientMap(), 
			 multiComm,
			 dgdx_sg->numVectors()));
	reduced_outargs.set_DgDx_sg(
	  i, SGDerivative(dgdx_red,
			  solver_outargs.get_DgDx_sg(i).getMultiVectorOrientation()));
      }
    }

    // dg/dp_red
    for (int j=0; j<solver_outargs.Np(); j++) {
      if (!solver_outargs.supports(OUT_ARG_DgDp_sg, i, j).none()) {
	Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> dgdp_sg = 
	  solver_outargs.get_DgDp_sg(i,j).getMultiVector();
	if (dgdp_sg != Teuchos::null) {
	  Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> dgdp_red = 
	    Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			   red_basis, red_overlap_map, 
			   dgdp_sg->coefficientMap(), 
			   multiComm,
			   dgdp_sg->numVectors()));
	  reduced_outargs.set_DgDp_sg(
	    i, j, SGDerivative(dgdp_red, 
			       solver_outargs.get_DgDp_sg(i,j).getMultiVectorOrientation()));
	}
      }
    }
  }

    
  // Setup new solver
  reduced_params = 
    Teuchos::rcp(new Teuchos::ParameterList(*solver_params));
  Teuchos::ParameterList& red_sg_params = 
    reduced_params->sublist("Stochastic Galerkin");
  red_sg_params.sublist("Basis").set("Stochastic Galerkin Basis",
				     red_basis);
  red_sg_params.sublist("Quadrature").set("Stochastic Galerkin Quadrature",
					  red_quad);
  if (red_sg_params.sublist("Expansion").isParameter("Stochastic Galerkin Expansion"))
    red_sg_params.sublist("Expansion").remove("Stochastic Galerkin Expansion");
  if (red_sg_params.isParameter("Triple Product Tensor"))
    red_sg_params.remove("Triple Product Tensor");
  Teuchos::RCP<Piro::Epetra::StokhosSolver> reduced_piro_solver = 
    Teuchos::rcp(new Piro::Epetra::StokhosSolver(reduced_params, comm));
  reduced_piro_solver->setup(model);
  reduced_solver = reduced_piro_solver;

  if (reduced_inargs.supports(IN_ARG_sg_basis))
    reduced_inargs.set_sg_basis(red_basis);
  if (reduced_inargs.supports(IN_ARG_sg_quadrature))
    reduced_inargs.set_sg_quadrature(red_quad);
  if (reduced_inargs.supports(IN_ARG_sg_expansion))
    reduced_inargs.set_sg_expansion(red_sg_params.sublist("Expansion").get< Teuchos::RCP< Stokhos::OrthogPolyExpansion<int,double,StorageType> > >("Stochastic Galerkin Expansion"));
}

void 
Piro::Epetra::NECoupledModelEvaluator:: 
do_dimension_projection(
  const InArgs& inArgs, 
  const InArgs& reduced_inargs, 
  const OutArgs& reduced_outargs,
  const Teuchos::RCP<const Teuchos::Array< Teuchos::Array<double> > >& red_basis_vals,
  OutArgs& solver_outargs) const
{
  TEUCHOS_FUNC_TIME_MONITOR("NECoupledModelEvaluator -- dimension projection");

  // Make sure there is something to do
  InArgs::sg_const_vector_t x_sg;
  if (supports_x_sg)
    x_sg = inArgs.get_x_sg();
  if (!reduce_dimension || x_sg == Teuchos::null)
    return;

  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis = 
    inArgs.get_sg_basis();
  Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad =
    inArgs.get_sg_quadrature();
  const Teuchos::Array<double>& weights = quad->getQuadWeights();
  const Teuchos::Array< Teuchos::Array<double> >& basis_vals = 
    quad->getBasisAtQuadPoints();
  int nqp = weights.size();
  const Teuchos::Array<double>& norms = basis->norm_squared();

  for (int i=0; i<solver_outargs.Ng(); i++) {

    // g_sg
    if (solver_outargs.supports(OUT_ARG_g_sg, i)) {
      OutArgs::sg_vector_t g_sg = solver_outargs.get_g_sg(i);
      if (g_sg != Teuchos::null) {
	OutArgs::sg_vector_t g_red = reduced_outargs.get_g_sg(i);
	Epetra_Vector g_val(*(g_red->coefficientMap()));
	g_sg->init(0.0);
	for (int qp=0; qp<nqp; qp++) {
	  g_red->evaluate((*red_basis_vals)[qp], g_val);
	  g_sg->sumIntoAllTerms(weights[qp], basis_vals[qp], norms, g_val);
	}
      }
    }

    // dg/dx_sg
    if (!solver_outargs.supports(OUT_ARG_DgDx_sg, i).none()) {
      Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> dgdx_sg = 
	solver_outargs.get_DgDx_sg(i).getMultiVector();
      if (dgdx_sg != Teuchos::null) {
	Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> dgdx_red = 
	  reduced_outargs.get_DgDx_sg(i).getMultiVector();
	Epetra_MultiVector dgdx_val(*(dgdx_red->coefficientMap()), 
				    dgdx_red->numVectors());
	dgdx_sg->init(0.0);
	for (int qp=0; qp<nqp; qp++) {
	  dgdx_red->evaluate((*red_basis_vals)[qp], dgdx_val);
	  dgdx_sg->sumIntoAllTerms(weights[qp], basis_vals[qp], norms, 
				   dgdx_val);
	}
      }
    }

    // dg/dp_sg
    for (int j=0; j<solver_outargs.Np(); j++) {
      if (!solver_outargs.supports(OUT_ARG_DgDp_sg, i, j).none()) {
	Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> dgdp_sg = 
	  solver_outargs.get_DgDp_sg(i,j).getMultiVector();
	if (dgdp_sg != Teuchos::null) {
	  Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> dgdp_red = 
	    reduced_outargs.get_DgDp_sg(i,j).getMultiVector();
	  Epetra_MultiVector dgdp_val(*(dgdp_red->coefficientMap()), 
				      dgdp_red->numVectors());
	  dgdp_sg->init(0.0);
	  for (int qp=0; qp<nqp; qp++) {
	    dgdp_red->evaluate((*red_basis_vals)[qp], dgdp_val);
	    dgdp_sg->sumIntoAllTerms(weights[qp], basis_vals[qp], norms, 
				     dgdp_val);
	  }
	}
      }
    }
  }
}
