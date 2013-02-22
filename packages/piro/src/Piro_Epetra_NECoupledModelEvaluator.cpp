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
#include "Stokhos_ReducedBasisFactory.hpp"
#include "EpetraExt_MultiComm.h"

Piro::Epetra::NECoupledModelEvaluator::
NECoupledModelEvaluator(
  const Teuchos::Array<Teuchos::RCP<EpetraExt::ModelEvaluator> >& models_,
  const Teuchos::Array<Teuchos::RCP<Teuchos::ParameterList> >& piroParams_,
  const Teuchos::RCP<AbstractNetworkModel>& network_model_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_,
  const Teuchos::RCP<const Epetra_Comm>& comm_,
  const Teuchos::Array< Teuchos::RCP<NOX::Epetra::Observer> >& observers_):
  models(models_),
  piroParams(piroParams_),
  network_model(network_model_),
  params(params_),
  comm(comm_),
  observers(observers_)
{
  // Setup VerboseObject
  Teuchos::readVerboseObjectSublist(params.get(), this);

  n_models = models.size();
  solvers.resize(n_models);

  // Create solvers for models A and B
  bool stochastic = params->get("Stochastic", false);
  if (observers.size() < n_models)
    observers.resize(n_models);
  if (stochastic) {
    sgSolvers.resize(n_models);
    for (int i=0; i<n_models; i++) {
      sgSolvers[i] = 
	Teuchos::rcp(new Piro::Epetra::StokhosSolver(piroParams[i], 
						     comm));
      sgSolvers[i]->setup(models[i], observers[i]);
      solvers[i] = sgSolvers[i];
    }
  }
  else {
    for (int i=0; i<n_models; i++)
      solvers[i] = Piro::Epetra::Factory::createSolver(piroParams[i], 
						       models[i]);
  }

  // Get connectivity information
  p_indices = 
    params->get< Teuchos::Array<int> >("Network Coupling Parameter Indices");
  g_indices = 
    params->get< Teuchos::Array<int> >("Network Coupling Response Indices");
  TEUCHOS_ASSERT(p_indices.size() == n_models);
  TEUCHOS_ASSERT(g_indices.size() == n_models);

  // Get number of parameter and response vectors
  solver_inargs.resize(n_models); 
  solver_outargs.resize(n_models);
  num_params.resize(n_models);
  num_responses.resize(n_models);
  num_params_total = 0;
  num_responses_total = 0;
  for (int i=0; i<n_models; i++) {
    solver_inargs[i] = solvers[i]->createInArgs();
    solver_outargs[i] = solvers[i]->createOutArgs();
    num_params[i] = solver_inargs[i].Np();
    num_responses[i] = solver_outargs[i].Ng();
    num_params_total += num_params[i];
    num_responses_total += num_responses[i];
  }
  num_params_total -= n_models;
  num_responses_total -= n_models;
  
  // Building indexing maps between coupled system parameters/responses and
  // individual components
  // Parameter vector i of this model evaluator corresponds to parameter
  // param_map[i].second for model param_map[i].first.  Similarly for the
  // responses
  for (int i=0; i<n_models; i++) {
    for (int j=0; j<num_params[i]; j++)
      if (j != p_indices[i])
	param_map.push_back(std::make_pair(i,j));
    for (int j=0; j<num_responses[i]; j++)
      if (j != g_indices[i])
	response_map.push_back(std::make_pair(i,j));
  }
  TEUCHOS_ASSERT(param_map.size() == num_params_total);
  TEUCHOS_ASSERT(response_map.size() == num_responses_total);

  // Build x map, which is the product of the p_indices parameter maps
  // For the time being, we will assume local maps, in the future we need to
  // build proper product maps
  p_maps.resize(n_models);
  n_p.resize(n_models);
  int nx = 0;
  for (int i=0; i<n_models; i++) {
    p_maps[i] = solvers[i]->get_p_map(p_indices[i]);
    n_p[i] = p_maps[i]->NumGlobalElements();
    nx += n_p[i];
  }
  x_map = Teuchos::rcp(new Epetra_Map(nx, 0, *comm));
  x_overlap_map = Teuchos::rcp(new Epetra_LocalMap(nx, 0, *comm));
  x_importer = Teuchos::rcp(new Epetra_Import(*x_overlap_map, *x_map));
  x_overlap = Teuchos::rcp(new Epetra_Vector(*x_overlap_map));

  // Build f map, which is the product of the g_indices response maps
  // For the time being, we will assume local maps, in the future we need to
  // build proper product maps
  g_maps.resize(n_models);
  n_g.resize(n_models);
  int nf = 0;
  for (int i=0; i<n_models; i++) {
    g_maps[i] = solvers[i]->get_g_map(g_indices[i]);
    n_g[i] = g_maps[i]->NumGlobalElements();
    nf += n_g[i];
  }
  f_map = Teuchos::rcp(new Epetra_Map(nf, 0, *comm));
  f_overlap_map = Teuchos::rcp(new Epetra_LocalMap(nf, 0, *comm));
  f_exporter = Teuchos::rcp(new Epetra_Export(*f_overlap_map, *f_map));
  f_overlap = Teuchos::rcp(new Epetra_Vector(*f_overlap_map));
  
  // Determine what we support
  supports_W = true;
  supports_x_sg = true;
  supports_f_sg = true;
  supports_W_sg = true;
  Teuchos::Array<DerivativeSupport> ds(n_models);
  for (int i=0; i<n_models; i++) {
    ds[i] = 
      solver_outargs[i].supports(OUT_ARG_DgDp, g_indices[i], p_indices[i]);
    DerivativeSupport ds_sg = 
      solver_outargs[i].supports(OUT_ARG_DgDp_sg, g_indices[i], p_indices[i]);
    supports_W = supports_W && 
      (ds[i].supports(DERIV_MV_BY_COL) || 
       ds[i].supports(DERIV_TRANS_MV_BY_ROW));
    supports_x_sg = supports_x_sg &&
      solver_inargs[i].supports(IN_ARG_p_sg, p_indices[i]);
    supports_f_sg = supports_f_sg && 
      solver_outargs[i].supports(OUT_ARG_g_sg, g_indices[i]);
    supports_W_sg = supports_W_sg &&
      (ds_sg.supports(DERIV_MV_BY_COL) || 
       ds_sg.supports(DERIV_TRANS_MV_BY_ROW));
  }

  // Build the Jacobian graph (currently dense)
  if (supports_W || supports_W_sg) {
    W_graph = Teuchos::rcp(new Epetra_CrsGraph(Copy, *f_map, nx));
    int *indices = f_overlap_map->MyGlobalElements();
    for (int i=0; i<f_map->NumMyElements(); i++) {
      int row = f_map->GID(i);
      W_graph->InsertGlobalIndices(row, nx, indices);
    }
    W_graph->FillComplete();

    W_overlap_graph = 
      Teuchos::rcp(new Epetra_CrsGraph(Copy, *f_overlap_map, nx));
    for (int i=0; i<f_overlap_map->NumMyElements(); i++) {
      int row = f_overlap_map->GID(i);
      W_overlap_graph->InsertGlobalIndices(row, nx, indices);
    }
    W_overlap_graph->FillComplete();
    W_overlap = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *W_overlap_graph));
  }

  // Build initial guess
  Epetra_Vector x_init_overlap(*x_overlap_map);
  int offset = 0;
  for (int i=0; i<n_models; i++) {
    Teuchos::RCP<const Epetra_Vector> p_init = 
      solvers[i]->get_p_init(p_indices[i]);
    for (int j=0; j<n_p[i]; j++)
      x_init_overlap[j+offset] = (*p_init)[j];
    offset += n_p[i];
  }
  x_init = Teuchos::rcp(new Epetra_Vector(*x_map));
  x_init->Export(x_init_overlap, *x_importer, Insert);

  // Create storage for parameters, responses, and derivatives
  p.resize(n_models);
  g.resize(n_models);
  dgdp_layout.resize(n_models);
  dgdp.resize(n_models);
  for (int i=0; i<n_models; i++) {
    p[i] = Teuchos::rcp(new Epetra_Vector(*p_maps[i]));
    g[i] = Teuchos::rcp(new Epetra_Vector(*g_maps[i]));
    if (supports_W || supports_W_sg) {
      if (ds[i].supports(DERIV_MV_BY_COL)) {
	dgdp_layout[i] = DERIV_MV_BY_COL;
	dgdp[i] = Teuchos::rcp(new Epetra_MultiVector(*g_maps[i], n_p[i]));
      }
      else {
	dgdp_layout[i] = DERIV_TRANS_MV_BY_ROW;
	dgdp[i] = Teuchos::rcp(new Epetra_MultiVector(*p_maps[i], n_g[i]));
      }
    }
  }

  p_sg.resize(n_models);
  g_sg.resize(n_models);
  dgdp_sg_layout.resize(n_models);
  dgdp_sg.resize(n_models);
  Teuchos::ParameterList& dim_reduct_params = 
    params->sublist("Dimension Reduction");
  if (!dim_reduct_params.isParameter("Reduce Dimension"))
    reduce_dimension.resize(n_models, 0);
  else if (dim_reduct_params.isType<bool>("Reduce Dimension"))
    reduce_dimension.resize(n_models, 
			    dim_reduct_params.get<bool>("Reduce Dimension"));
  else if (dim_reduct_params.isType< Teuchos::Array<int> >("Reduce Dimension"))
    reduce_dimension = 
      dim_reduct_params.get< Teuchos::Array<int> >("Reduce Dimension");
  else
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, 
      "Invalid type for parameter \"Dimension Reduction\"");
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

  return solvers[param_map[j].first]->get_p_map(param_map[j].second);
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

  return solvers[response_map[j].first]->get_g_map(response_map[j].second);
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

  return solvers[param_map[j].first]->get_p_names(param_map[j].second);
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

  return solvers[param_map[j].first]->get_p_init(param_map[j].second);
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

  // Deterministic InArgs
  inArgs.setSupports(IN_ARG_x, true);
  inArgs.set_Np(num_params_total);
  
  // Stochastic InArgs
  if (supports_x_sg) {
    inArgs.setSupports(IN_ARG_x_sg, supports_x_sg);
    inArgs.setSupports(IN_ARG_sg_basis,true);
    inArgs.setSupports(IN_ARG_sg_quadrature,true);
    inArgs.setSupports(IN_ARG_sg_expansion,true);
    for (int i=0; i<num_params_total; i++) {
      inArgs.setSupports(
	IN_ARG_p_sg, i, 
	solver_inargs[param_map[i].first].supports(IN_ARG_p_sg, 
						   param_map[i].second)); 
    }
  }

  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
Piro::Epetra::NECoupledModelEvaluator::
createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());

  // Deterministic OutArgs
  outArgs.setSupports(OUT_ARG_f, true);
  outArgs.setSupports(OUT_ARG_W, supports_W);
  outArgs.set_W_properties(
    DerivativeProperties(DERIV_LINEARITY_NONCONST, DERIV_RANK_FULL, true));
  outArgs.set_Np_Ng(num_params_total, num_responses_total);
  for (int i=0; i<num_params_total; i++)
    outArgs.setSupports(
      OUT_ARG_DfDp, i, solver_outargs[param_map[i].first].supports(
	OUT_ARG_DgDp, g_indices[param_map[i].first], param_map[i].second)); 
  for (int i=0; i<num_responses_total; i++) {
    int model_index = response_map[i].first;
    outArgs.setSupports(
      OUT_ARG_DgDx, i, 
      solver_outargs[model_index].supports(OUT_ARG_DgDx, 
					   response_map[i].second));
    for (int j=0; j<num_params_total; j++) {
      if (param_map[j].first == model_index) {
	outArgs.setSupports(
	  OUT_ARG_DgDp, i, j, 
	  solver_outargs[model_index].supports(OUT_ARG_DgDp, 
					       response_map[i].second,
					       param_map[j].second));
      }
    }
  }

  // Stochastic OutArgs
  if (supports_f_sg)
    outArgs.setSupports(OUT_ARG_f_sg, supports_f_sg);
  if (supports_W_sg)
    outArgs.setSupports(OUT_ARG_W_sg, supports_W_sg);
  for (int i=0; i<num_responses_total; i++)
    outArgs.setSupports(
      OUT_ARG_g_sg, i, solver_outargs[response_map[i].first].supports(
	OUT_ARG_g_sg, response_map[i].second));
  for (int i=0; i<num_params_total; i++)
    outArgs.setSupports(
      OUT_ARG_DfDp_sg, i, solver_outargs[param_map[i].first].supports(
	OUT_ARG_DgDp_sg, g_indices[param_map[i].first], param_map[i].second)); 
  for (int i=0; i<num_responses_total; i++) {
    int model_index = response_map[i].first;
    outArgs.setSupports(
      OUT_ARG_DgDx_sg, i, 
      solver_outargs[model_index].supports(OUT_ARG_DgDx_sg, 
					   response_map[i].second));
    for (int j=0; j<num_params_total; j++) {
      if (param_map[j].first == model_index) {
	outArgs.setSupports(
	  OUT_ARG_DgDp_sg, i, j, 
	  solver_outargs[model_index].supports(OUT_ARG_DgDp_sg, 
					       response_map[i].second,
					       param_map[j].second));
      }
    }
  }
  
  return outArgs;
}

void 
Piro::Epetra::NECoupledModelEvaluator::
evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();

  // Create fresh in/out args for sub-models
  for (int i=0; i<n_models; i++) {
    solver_inargs[i] = solvers[i]->createInArgs();
    solver_outargs[i] = solvers[i]->createOutArgs();
  }

  EpetraExt::ModelEvaluator::InArgs network_inargs = inArgs;
  EpetraExt::ModelEvaluator::OutArgs network_outargs = outArgs;

  //
  // Deterministic calculation
  //
  bool do_deterministic = false;
  Teuchos::RCP<const Epetra_Vector> x = inArgs.get_x();
  if (x != Teuchos::null) {
    do_deterministic = true;

    // p
    x_overlap->Import(*x, *x_importer, Insert);
    int offset = 0;
    for (int i=0; i<n_models; i++) {
      for (int j=0; j<n_p[i]; j++)
	(*p[i])[j] = (*x_overlap)[j+offset];
      offset += n_p[i];
      solver_inargs[i].set_p(p_indices[i], p[i]);
      for (int j=0; j<num_params_total; j++)
	if (param_map[j].first == i)
	  solver_inargs[i].set_p(param_map[j].second, inArgs.get_p(j));
    }
    network_inargs.set_x(x_overlap);
    
    // f
    Teuchos::RCP<Epetra_Vector> f = outArgs.get_f();
    if (f != Teuchos::null) {
      for (int i=0; i<n_models; i++)
	solver_outargs[i].set_g(g_indices[i], g[i]);
      network_outargs.set_f(f_overlap);
    }

    // W
    Teuchos::RCP<Epetra_Operator> W = outArgs.get_W();
    if (W != Teuchos::null) {
      for (int i=0; i<n_models; i++) {
	Derivative dgdp_deriv(dgdp[i], dgdp_layout[i]);
	solver_outargs[i].set_DgDp(g_indices[i], p_indices[i], dgdp_deriv);
      }
      network_outargs.set_W(W_overlap);
    }

    for (int i=0; i<num_responses_total; i++) {
      int model_index = response_map[i].first;
      
      // g
      solver_outargs[model_index].set_g(response_map[i].second, 
					outArgs.get_g(i));
    

      // dg/dx
      if (!solver_outargs[model_index].supports(OUT_ARG_DgDx, 
						response_map[i].second).none())
	solver_outargs[model_index].set_DgDx(response_map[i].second, 
					     outArgs.get_DgDx(i));
    

      // dg/dp
      for (int j=0; j<num_params_total; j++) {
	if (param_map[j].first == model_index) {
	  if (!solver_outargs[model_index].supports(OUT_ARG_DgDp, 
						    response_map[i].second, 
						    param_map[j].second).none())
	    solver_outargs[model_index].set_DgDp(response_map[i].second, 
						 param_map[j].second, 
						 outArgs.get_DgDp(i,j));
	}
      }

    }
    
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
    Teuchos::RCP<const EpetraExt::MultiComm> multiComm = 
      x_sg->productComm();
    if (sg_overlap_map == Teuchos::null)
      sg_overlap_map =
	Teuchos::rcp(new Epetra_LocalMap(basis->size(), 0, 
					 multiComm->TimeDomainComm()));

    if (x_sg_overlap == Teuchos::null)
      x_sg_overlap = 
	Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		       basis, sg_overlap_map, x_overlap_map, multiComm));
    if (supports_f_sg && f_sg_overlap == Teuchos::null)
      f_sg_overlap = 
	Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
		       basis, sg_overlap_map, f_overlap_map, multiComm));
    if (supports_W_sg && W_sg_overlap == Teuchos::null) {
      Teuchos::RCP<const Epetra_Map> domain_base_map = 
	x_overlap_map;
      Teuchos::RCP<const Epetra_Map> range_base_map = 
	f_overlap_map;
      W_sg_overlap =
	Teuchos::rcp(new Stokhos::EpetraOperatorOrthogPoly(
		       basis, sg_overlap_map, domain_base_map, range_base_map, 
		       multiComm));
      for (int block=0; block<W_sg_overlap->size(); block++) {
	Teuchos::RCP<Epetra_Operator> W =
	  Teuchos::rcp(new Epetra_CrsMatrix(Copy, *W_overlap_graph));
	W_sg_overlap->setCoeffPtr(block,W);
      }
    }

    for (int i=0; i<n_models; i++) {
      if (solver_inargs[i].supports(IN_ARG_sg_basis))
	solver_inargs[i].set_sg_basis(basis);
      if (solver_inargs[i].supports(IN_ARG_sg_quadrature))
	solver_inargs[i].set_sg_quadrature(inArgs.get_sg_quadrature());
      if (solver_inargs[i].supports(IN_ARG_sg_expansion))
	solver_inargs[i].set_sg_expansion(inArgs.get_sg_expansion());
      
      if (p_sg[i] == Teuchos::null) {
	p_sg[i] = 
	  Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
			 basis, sg_overlap_map, p_maps[i], multiComm));
	if (supports_f_sg) {
	  g_sg[i] = 
	    Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(
			   basis, sg_overlap_map, g_maps[i], multiComm));
	}
	if (supports_W_sg) {
	  DerivativeSupport ds_sg = 
	    solver_outargs[i].supports(OUT_ARG_DgDp_sg, g_indices[i], 
				       p_indices[i]);
	  if (ds_sg.supports(DERIV_MV_BY_COL)) {
	    dgdp_sg_layout[i] = DERIV_MV_BY_COL;
	    dgdp_sg[i] = 
	      Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			     basis, sg_overlap_map, g_maps[i], multiComm, 
			     n_p[i]));
	  }
	  else {
	    dgdp_sg_layout[i] = DERIV_TRANS_MV_BY_ROW;
	    dgdp_sg[i] = 
	      Teuchos::rcp(new Stokhos::EpetraMultiVectorOrthogPoly(
			     basis, sg_overlap_map, p_maps[i], multiComm, 
			     n_g[i]));
	  }
	}
      }
    }
      
    // p_sg
    for (int block=0; block<x_sg->size(); block++) {
      (*x_sg_overlap)[block].Import((*x_sg)[block], *x_importer, Insert);
      int offset = 0;
      for (int i=0; i<n_models; i++) {
	for (int j=0; j<n_p[i]; j++)
	  (*p_sg[i])[block][j] = (*x_sg_overlap)[block][j+offset];
	offset += n_p[i];
      }
      network_inargs.set_x_sg(x_sg_overlap);
    }
      
    for (int i=0; i<n_models; i++) {
      solver_inargs[i].set_p_sg(p_indices[i], p_sg[i]);
      for (int j=0; j<num_params_total; j++)
	if (param_map[j].first == i)
	  if (solver_inargs[i].supports(IN_ARG_p_sg, param_map[i].second))
	    solver_inargs[i].set_p_sg(param_map[j].second, 
				      inArgs.get_p_sg(j));
    }
    
    // f_sg
    if (supports_f_sg) {
      OutArgs::sg_vector_t f_sg = outArgs.get_f_sg();
      if (f_sg != Teuchos::null) {
	for (int i=0; i<n_models; i++)
	  solver_outargs[i].set_g_sg(g_indices[i], g_sg[i]);
	network_outargs.set_f_sg(f_sg_overlap);
      }
    }
      
    // W_sg
    if (supports_W_sg) {
      OutArgs::sg_operator_t W_sg = outArgs.get_W_sg();
      if (W_sg != Teuchos::null) {
	for (int i=0; i<n_models; i++) {
	  SGDerivative dgdp_deriv(dgdp_sg[i], dgdp_sg_layout[i]);
	  solver_outargs[i].set_DgDp_sg(g_indices[i], p_indices[i], dgdp_deriv);
	}
	network_outargs.set_W_sg(W_sg_overlap);
      }
    }

    for (int i=0; i<num_responses_total; i++) {
      int model_index = response_map[i].first;
      
      // g_sg
      if (solver_outargs[model_index].supports(OUT_ARG_g_sg, 
					       response_map[i].second))
	solver_outargs[model_index].set_g_sg(response_map[i].second, 
					     outArgs.get_g_sg(i));
      
      // dg/dx_sg
      if (!solver_outargs[model_index].supports(OUT_ARG_DgDx_sg, 
						response_map[i].second).none())
	solver_outargs[model_index].set_DgDx_sg(response_map[i].second, 
						outArgs.get_DgDx_sg(i));
      
      
      // dg/dp_sg
      for (int j=0; j<num_params_total; j++) {
	if (param_map[j].first == model_index) {
	  if (!solver_outargs[model_index].supports(OUT_ARG_DgDp_sg, 
						    response_map[i].second, 
						    param_map[j].second).none())
	    solver_outargs[model_index].set_DgDp_sg(response_map[i].second, 
						    param_map[j].second, 
						    outArgs.get_DgDp_sg(i,j));
	}
      }
      
    } 
  }

  // Stochastic dimension reduction
  Teuchos::Array<EpetraExt::ModelEvaluator::InArgs> solver_inargs_red(n_models);
  Teuchos::Array<EpetraExt::ModelEvaluator::OutArgs> solver_outargs_red(n_models);
  Teuchos::Array<Teuchos::RCP<EpetraExt::ModelEvaluator> > solvers_red(n_models);
  
  Teuchos::Array<Teuchos::RCP<Teuchos::ParameterList> > piroParams_red(n_models);
  
  for (int i=0; i<n_models; i++) {
    do_dimension_reduction(i, inArgs, 
			   solver_inargs[i], solver_outargs[i],
			   models[i], solvers[i], piroParams[i],
			   solver_inargs_red[i], solver_outargs_red[i],
			   solvers_red[i], piroParams_red[i]);
  }
  

  // Evaluate models
  if (n_models == 2) {
    {
      TEUCHOS_FUNC_TIME_MONITOR(
	"NECoupledModelEvaluator -- Model 1 nonlinear elimination");
      if (verbLevel != Teuchos::VERB_NONE)
	*out << "Eliminating model " << 1 << " states...";
      solvers_red[0]->evalModel(solver_inargs_red[0], solver_outargs_red[0]);
    }

    {
      TEUCHOS_FUNC_TIME_MONITOR(
	"NECoupledModelEvaluator -- Model 2 nonlinear elimination");
      if (verbLevel != Teuchos::VERB_NONE)
	*out << "Eliminating model " << 2 << " states...";
      solvers_red[1]->evalModel(solver_inargs_red[1], solver_outargs_red[1]);
    }
  }
  else {
    for (int i=0; i<n_models; i++) {
      TEUCHOS_FUNC_TIME_MONITOR(
	"NECoupledModelEvaluator -- Model nonlinear elimination");
      if (verbLevel != Teuchos::VERB_NONE)
	*out << "Eliminating model " << i+1 << " states...";
      solvers_red[i]->evalModel(solver_inargs_red[i], solver_outargs_red[i]);
    }
  }

  // Project back to original stochastic bases
  for (int i=0; i<n_models; i++)
    do_dimension_projection(i, inArgs, solver_inargs_red[i], 
			    solver_outargs_red[i], solver_outargs[i]);

  // Evaluate network model
  network_model->evalModel(solver_inargs, solver_outargs, 
			   network_inargs, network_outargs,
			   n_p, n_g, p, g, dgdp, dgdp_layout,
			   p_sg, g_sg, dgdp_sg, dgdp_sg_layout);

  // Export network residuals, Jacobians, etc...

  // f
  Teuchos::RCP<Epetra_Vector> f = outArgs.get_f();
  if (f != Teuchos::null)
    f->Export(*f_overlap, *f_exporter, Insert);

  // W
  Teuchos::RCP<Epetra_Operator> W = outArgs.get_W();
  if (W != Teuchos::null) {
    Teuchos::RCP<Epetra_CrsMatrix> W_crs = 
      Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W, true);
    W_crs->Export(*W_overlap, *f_exporter, Insert);
  }

  // f_sg
  if (supports_f_sg) {
    OutArgs::sg_vector_t f_sg = outArgs.get_f_sg();
    if (f_sg != Teuchos::null) {
      for (int block=0; block<f_sg->size(); block++)
	(*f_sg)[block].Export((*f_sg_overlap)[block], *f_exporter, Insert);
    }
  }

  // W_sg
  if (supports_W_sg) {
    OutArgs::sg_operator_t W_sg = outArgs.get_W_sg();
    if (W_sg != Teuchos::null) {
      for (int block=0; block<W_sg->size(); block++) {
	Teuchos::RCP<Epetra_CrsMatrix> W_crs = 
	  Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(
	    W_sg->getCoeffPtr(block), true);
	Teuchos::RCP<Epetra_CrsMatrix> W_overlap_crs = 
	  Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(
	    W_sg_overlap->getCoeffPtr(block), true);
	W_crs->Export(*W_overlap_crs, *f_exporter, Insert);
      }
    }
  }
  }


void
Piro::Epetra::NECoupledModelEvaluator:: 
do_dimension_reduction(
  int model_index,
  const InArgs& inArgs,
  const InArgs& solver_inargs, 
  const OutArgs& solver_outargs,
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& model,
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& solver,
  const Teuchos::RCP<Teuchos::ParameterList>& solver_params,
  InArgs& reduced_inargs, 
  OutArgs& reduced_outargs,
  Teuchos::RCP<EpetraExt::ModelEvaluator>& reduced_solver,
  Teuchos::RCP<Teuchos::ParameterList>& reduced_params) const
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
  if (!reduce_dimension[model_index] || x_sg == Teuchos::null)
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
  Teuchos::ParameterList& reduct_params = 
    params->sublist("Dimension Reduction");
  int order = basis->order();
  int new_order = reduct_params.get("Reduced Order", -1);
  if (new_order == -1)
    new_order = order;
  if (st_quad == Teuchos::null) {
    st_quad = quad;
    // st_quad =
    //   Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(
    // 		 basis, new_order+1));
    // st_quad =
    //   Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(
    // 		 basis, 4*new_order+1));
    // std::cout << "st_quad->size() = " << st_quad->size() << std::endl;
  }
  Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk = 
    expansion->getTripleProduct();
  Stokhos::ReducedBasisFactory<int,double> factory(reduct_params);
  Teuchos::RCP< Stokhos::ReducedPCEBasis<int,double> > gs_basis = 
    factory.createReducedBasis(new_order, p_opa, st_quad, Cijk);
  red_basis = gs_basis;
  red_quad = gs_basis->getReducedQuadrature();
  red_pces.resize(p_opa.size());
  for (int i=0; i<p_opa.size(); i++) {
    red_pces[i].reset(red_basis);
    gs_basis->transformFromOriginalBasis(p_opa[i].coeff(), red_pces[i].coeff());
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
  reduced_piro_solver->setup(model, observers[model_index]);
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
  int model_index,
  const InArgs& inArgs, 
  const InArgs& reduced_inargs, 
  const OutArgs& reduced_outargs,
  OutArgs& solver_outargs) const
{
  TEUCHOS_FUNC_TIME_MONITOR("NECoupledModelEvaluator -- dimension projection");

  // Make sure there is something to do
  InArgs::sg_const_vector_t x_sg;
  if (supports_x_sg)
    x_sg = inArgs.get_x_sg();
  if (!reduce_dimension[model_index] || x_sg == Teuchos::null)
    return;

  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis = 
    inArgs.get_sg_basis();
  Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad =
    inArgs.get_sg_quadrature();
  Teuchos::RCP<const Stokhos::ReducedPCEBasis<int,double> > red_basis = 
    Teuchos::rcp_dynamic_cast<const Stokhos::ReducedPCEBasis<int,double> >(reduced_inargs.get_sg_basis());

  for (int i=0; i<solver_outargs.Ng(); i++) {

    // g_sg
    if (solver_outargs.supports(OUT_ARG_g_sg, i)) {
      OutArgs::sg_vector_t g_sg = solver_outargs.get_g_sg(i);
      if (g_sg != Teuchos::null) {
	OutArgs::sg_vector_t g_red = reduced_outargs.get_g_sg(i);
	red_basis->transformToOriginalBasis(
	  (*g_red)[0].Values(), 
	  (*g_sg)[0].Values(), 
	  g_red->coefficientMap()->NumMyElements(), 
	  true);
      }
    }

    // dg/dx_sg
    if (!solver_outargs.supports(OUT_ARG_DgDx_sg, i).none()) {
      Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> dgdx_sg = 
	solver_outargs.get_DgDx_sg(i).getMultiVector();
      if (dgdx_sg != Teuchos::null) {
	Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> dgdx_red = 
	  reduced_outargs.get_DgDx_sg(i).getMultiVector();

	// transformToOriginalBasis() needs the entries for each pce
	  // coefficient stored contiguously.  This isn't the case for the
	  // full multivector (each column along with all of its pce
	  // coefficients is stored in one contiguous chunk).  Thus we need
	  // to transform each column individually
	  int ncol = dgdx_red->numVectors();
	  for (int col=0; col<ncol; col++)
	    red_basis->transformToOriginalBasis(
	      (*dgdx_red)[0](col)->Values(), 
	      (*dgdx_sg)[0](col)->Values(), 
	      dgdx_red->coefficientMap()->NumMyElements(), 
	      true);
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

	  // transformToOriginalBasis() needs the entries for each pce
	  // coefficient stored contiguously.  This isn't the case for the
	  // full multivector (each column along with all of its pce
	  // coefficients is stored in one contiguous chunk).  Thus we need
	  // to transform each column individually
	  int ncol = dgdp_red->numVectors();
	  for (int col=0; col<ncol; col++)
	    red_basis->transformToOriginalBasis(
	      (*dgdp_red)[0](col)->Values(), 
	      (*dgdp_sg)[0](col)->Values(), 
	      dgdp_red->coefficientMap()->NumMyElements(), 
	      true);
	}
      }
    }
  }
}


void 
Piro::Epetra::ParamToResponseNetworkModel::
evalModel(
  const Teuchos::Array<EpetraExt::ModelEvaluator::InArgs>& model_inargs, 
  const Teuchos::Array<EpetraExt::ModelEvaluator::OutArgs>& model_outargs,
  const EpetraExt::ModelEvaluator::InArgs& network_inargs, 
  const EpetraExt::ModelEvaluator::OutArgs& network_outargs,
  const Teuchos::Array<int>& n_p,
  const Teuchos::Array<int>& n_g,
  const Teuchos::Array< Teuchos::RCP<Epetra_Vector> >& p,
  const Teuchos::Array< Teuchos::RCP<Epetra_Vector> >& g,
  const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dgdp,
  const Teuchos::Array<EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation >& dgdp_layout,
  const Teuchos::Array<EpetraExt::ModelEvaluator::OutArgs::sg_vector_t>& p_sg,
  const Teuchos::Array<EpetraExt::ModelEvaluator::OutArgs::sg_vector_t>& g_sg,
  const Teuchos::Array<Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly> >& dgdp_sg,
  const Teuchos::Array<EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation>& dgdp_sg_layout) const
{

  // f
  Teuchos::RCP<Epetra_Vector> f = network_outargs.get_f();
  if (f != Teuchos::null) {
    f->PutScalar(0.0);
    for (int i=0; i<n_p[0]; i++)
      (*f)[i]        = (*p[0])[i] - (*g[1])[i];
    for (int i=0; i<n_p[1]; i++)
      (*f)[i+n_p[0]] = (*p[1])[i] - (*g[0])[i];
  }

  // W
  Teuchos::RCP<Epetra_Operator> W = network_outargs.get_W();
  if (W != Teuchos::null) {
    Teuchos::RCP<Epetra_CrsMatrix> W_crs = 
      Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W, true);
    W_crs->PutScalar(0.0);
    int row, col;
    double val;
    for (int i=0; i<n_p[0]; i++) {
      row = i; 
      
      // Diagonal part
      col = row; 
      val = 1.0;
      W_crs->ReplaceGlobalValues(row, 1, &val, &col);
      
      // dg_2/dp_2 part
      for (int j=0; j<n_p[1]; j++) {
	col = n_p[0]+j; 
	if (dgdp_layout[1] == EpetraExt::ModelEvaluator::DERIV_MV_BY_COL)
	  val = -(*dgdp[1])[j][i];
	else
	  val = -(*dgdp[1])[i][j];
	W_crs->ReplaceGlobalValues(row, 1, &val, &col);
      }
    }
    for (int i=0; i<n_p[1]; i++) {
      row = n_p[0] + i; 
      
      // Diagonal part
      col = row; 
      val = 1.0;
      W_crs->ReplaceGlobalValues(row, 1, &val, &col);
      
      // dg_1/dp_1 part
      for (int j=0; j<n_p[0]; j++) {
	col = j; 
	if (dgdp_layout[0] == EpetraExt::ModelEvaluator::DERIV_MV_BY_COL)
	  val = -(*dgdp[0])[j][i];
	else
	  val = -(*dgdp[0])[i][j];
	W_crs->ReplaceGlobalValues(row, 1, &val, &col);
      }
    }
  }

  // f_sg
  if (network_outargs.supports(EpetraExt::ModelEvaluator::OUT_ARG_f_sg)) {
    EpetraExt::ModelEvaluator::OutArgs::sg_vector_t f_sg = 
      network_outargs.get_f_sg();
    if (f_sg != Teuchos::null) {
      // std::cout << "g_sg[0] = " << *g_sg[0] << std::endl;
      // std::cout << "g_sg[1] = " << *g_sg[1] << std::endl;
      f_sg->init(0.0);
      for (int block=0; block<f_sg->size(); block++) {
	for (int i=0; i<n_p[0]; i++)
	  (*f_sg)[block][i] = 
	    (*p_sg[0])[block][i] - (*g_sg[1])[block][i];
	for (int i=0; i<n_p[1]; i++)
	  (*f_sg)[block][i+n_p[0]] = 
	    (*p_sg[1])[block][i] - (*g_sg[0])[block][i];
      }
      //std::cout << "f_sg = " << *f_sg << std::endl;
    }
  }
  
  // W_sg
  if (network_outargs.supports(EpetraExt::ModelEvaluator::OUT_ARG_W_sg)) {
    EpetraExt::ModelEvaluator::OutArgs::sg_operator_t W_sg = 
      network_outargs.get_W_sg();
    if (W_sg != Teuchos::null) {
      // std::cout << "dgdp_sg[0] = " << *dgdp_sg[0] << std::endl;
      // std::cout << "dgdp_sg[1] = " << *dgdp_sg[1] << std::endl;
      for (int block=0; block<W_sg->size(); block++) {
	Teuchos::RCP<Epetra_CrsMatrix> W_crs = 
	  Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_sg->getCoeffPtr(block), 
						      true);
	W_crs->PutScalar(0.0);
	int row, col;
	double val;
	for (int i=0; i<n_p[0]; i++) {
	  row = i; 
	  
	  // Diagonal part
	  if (block == 0) {
	    col = row; 
	    val = 1.0;
	    W_crs->ReplaceGlobalValues(row, 1, &val, &col);
	  }
	  
	  // dg_2/dp_2 part
	  for (int j=0; j<n_p[1]; j++) {
	    col = n_p[0]+j; 
	    if (dgdp_layout[1] == EpetraExt::ModelEvaluator::DERIV_MV_BY_COL)
	      val = -(*dgdp_sg[1])[block][j][i];
	    else
	      val = -(*dgdp_sg[1])[block][i][j];
	    W_crs->ReplaceGlobalValues(row, 1, &val, &col);
	  }
	}
	for (int i=0; i<n_p[1]; i++) {
	  row = n_p[0] + i; 
	  
	  // Diagonal part
	  if (block == 0) {
	    col = row; 
	    val = 1.0;
	    W_crs->ReplaceGlobalValues(row, 1, &val, &col);
	  }
	  
	  // dg_1/dp_1 part
	  for (int j=0; j<n_p[0]; j++) {
	    col = j; 
	    if (dgdp_layout[0] == EpetraExt::ModelEvaluator::DERIV_MV_BY_COL)
	      val = -(*dgdp_sg[0])[block][j][i];
	    else
	      val = -(*dgdp_sg[0])[block][i][j];
	    W_crs->ReplaceGlobalValues(row, 1, &val, &col);
	  }
	}
      }
    }
  }
}
