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

#include "Piro_Epetra_StokhosSolver.hpp"

#include "Piro_Epetra_SolverFactory.hpp"
#include "Piro_ExtensibleFactory.hpp"

#include "Stokhos.hpp"
#include "Stokhos_Epetra.hpp"

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "NOX_Epetra_ModelEvaluatorInterface.H"
#include "NOX_Epetra_LinearSystem_Stratimikos.H"
#include "NOX_Epetra_LinearSystem_MPBD.hpp"
#include "NOX_Epetra_LinearSystem_SGGS.hpp"
#include "NOX_Epetra_LinearSystem_SGJacobi.hpp"

Piro::Epetra::StokhosSolverFactory::
StokhosSolverFactory(const Teuchos::RCP<Teuchos::ParameterList>& piroParams_,
	      const Teuchos::RCP<const Epetra_Comm>& globalComm) :
  piroParams(piroParams_)
{
  // Setup VerboseObject
  Teuchos::readVerboseObjectSublist(piroParams.get(), this);
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();

  // Validate parameters
  Teuchos::ParameterList& sgParams =
    piroParams->sublist("Stochastic Galerkin");
  sgParams.validateParameters(*getValidSGParameters(),0);

  sgSolverParams = 
    //Teuchos::rcp(&(sgParams.sublist("SG Solver Parameters")),false);
    Teuchos::rcp(new Teuchos::ParameterList(sgParams.sublist("SG Solver Parameters")));

  // Get SG expansion type
  std::string sg_type = sgParams.get("SG Method", "Direct");
  if (sg_type == "Direct" || sg_type == "AD")
    sg_method = SG_AD;
  else if (sg_type == "Global")
    sg_method = SG_GLOBAL;
  else if (sg_type == "Non-intrusive")
    sg_method = SG_NI;
  else if (sg_type == "Multi-point Non-intrusive")
    sg_method = SG_MPNI;
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
		       std::endl << "Error!  ENAT_SGNOXSolver():  " <<
		       "Invalid SG Method  " << sg_type << std::endl);
  
  // Create SG basis
  basis = Stokhos::BasisFactory<int,double>::create(sgParams);
  if (verbLevel != Teuchos::VERB_NONE)
    *out << "Basis size = " << basis->size() << std::endl;

  // Create SG Quadrature
  Teuchos::ParameterList& expParams = sgParams.sublist("Expansion");
  std::string exp_type = expParams.get("Type", "Quadrature");
  if (exp_type == "Quadrature" || 
      sg_method == SG_GLOBAL ||
      sg_method == SG_NI ||
      sg_method == SG_MPNI) {
    quad = Stokhos::QuadratureFactory<int,double>::create(sgParams);
    if (verbLevel != Teuchos::VERB_NONE)
      *out << "Quadrature size = " << quad->size() << std::endl;
  }

  // Create SG expansion & triple-product
  if (sg_method != SG_NI && sg_method != SG_MPNI) {
    expansion = 
      Stokhos::ExpansionFactory<int,double>::create(sgParams);
    Cijk = 
      sgParams.get< Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > >("Triple Product Tensor");
  }

  // Create stochastic parallel distribution
  int num_spatial_procs = 
    sgParams.get("Number of Spatial Processors", -1);
  int num_stoch_blocks;
  if (sg_method == SG_MPNI)
    num_stoch_blocks = quad->size();
  else
    num_stoch_blocks = basis->size();
  sg_comm =
    Stokhos::buildMultiComm(*globalComm, num_stoch_blocks, num_spatial_procs);
  sg_parallel_data =
    Teuchos::rcp(new Stokhos::ParallelData(basis, Cijk, sg_comm, sgParams));

}

void 
Piro::Epetra::StokhosSolverFactory::
resetSolverParameters(const Teuchos::ParameterList& new_solver_params)
{
  *sgSolverParams = new_solver_params;
}

Teuchos::RCP<Stokhos::SGModelEvaluator>
Piro::Epetra::StokhosSolverFactory::
createSGModel(const Teuchos::RCP<EpetraExt::ModelEvaluator>& model_)
{
  Teuchos::ParameterList& sgParams =
    piroParams->sublist("Stochastic Galerkin");
  sgParams.sublist("Basis");

  model = model_;

  // Set up stochastic Galerkin model
  Teuchos::RCP<EpetraExt::ModelEvaluator> sg_model;
  if (sg_method == SG_AD) {
    sg_model = model;
  }
  else if (sg_method == SG_MPNI) {
    int num_mp = quad->size();
    Teuchos::RCP<const Epetra_Comm> mp_comm = 
      Stokhos::getStochasticComm(sg_comm);
    Teuchos::RCP<const Epetra_Map> mp_block_map = 
      Teuchos::rcp(new Epetra_Map(num_mp, 0, *mp_comm));
    Teuchos::RCP<EpetraExt::ModelEvaluator> mp_model = model;

    // Turn mp_model into an MP-nonlinear problem
    Teuchos::RCP<Teuchos::ParameterList> mpParams = 
    Teuchos::rcp(&(sgParams.sublist("MP Solver Parameters")),false);
    Teuchos::RCP<Stokhos::MPModelEvaluator> mp_nonlinear_model =
      Teuchos::rcp(new Stokhos::MPModelEvaluator(mp_model, sg_comm,
						 mp_block_map, mpParams));

    bool use_mpbd_solver = mpParams->get("Use MPBD Solver", false);
    Teuchos::RCP<NOX::Epetra::LinearSystem> linsys;
    Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> nox_interface;
    if (use_mpbd_solver) {
      nox_interface = 
	Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(mp_nonlinear_model));
      Teuchos::RCP<Epetra_Operator> A = 
	mp_nonlinear_model->create_W();
      Teuchos::RCP<Epetra_Operator> M = 
	mp_nonlinear_model->create_WPrec()->PrecOp;
      Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = 
	nox_interface;
      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = 
	nox_interface;
      Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = 
	nox_interface;

      Teuchos::ParameterList& noxParams = piroParams->sublist("NOX");
      Teuchos::ParameterList& printParams = noxParams.sublist("Printing");
      Teuchos::ParameterList& newtonParams = 
	noxParams.sublist("Direction").sublist("Newton");
      Teuchos::ParameterList& noxstratlsParams = 
	newtonParams.sublist("Stratimikos Linear Solver");
      Teuchos::ParameterList& mpbdParams = 
	mpParams->sublist("MPBD Linear Solver");
      mpbdParams.sublist("Deterministic Solver Parameters") = 
	noxstratlsParams;
      Teuchos::RCP<Epetra_Operator> inner_A = model->create_W();
      Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> inner_nox_interface = 
	Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(model));
      Teuchos::RCP<NOX::Epetra::Interface::Required> inner_iReq = 
	inner_nox_interface;
      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> inner_iJac = 
	inner_nox_interface;
      Teuchos::RCP<const Epetra_Vector> inner_u = model->get_x_init();
      Teuchos::RCP<NOX::Epetra::LinearSystem> inner_linsys = 
	Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(
		       printParams, 
		       noxstratlsParams,
		       inner_iJac, inner_A, *inner_u));
      linsys = 
	Teuchos::rcp(new NOX::Epetra::LinearSystemMPBD(printParams, 
						       mpbdParams,
						       inner_linsys,
						       iReq, iJac, A,
						       model->get_x_map()));
    }

    Piro::Epetra::SolverFactory solverFactory;
    {
      const std::string token = "Interface";
      solverFactory.setNOXInterfaceProvider(token, nox_interface);
      Teuchos::sublist(piroParams, token)->set("Type", token);
    }
    {
      const std::string token = "Linear System";
      solverFactory.setNOXLinearSystemProvider(token, linsys);
      Teuchos::sublist(piroParams, token)->set("Type", token);
    }

    // Create solver to map p -> g
    const Teuchos::RCP<EpetraExt::ModelEvaluator>  mp_solver
      = solverFactory.createSolver(piroParams, mp_nonlinear_model);

    // Create MP inverse model evaluator to map p_mp -> g_mp
    Teuchos::Array<int> mp_p_index_map = 
      mp_nonlinear_model->get_p_mp_map_indices();
    Teuchos::Array<int> mp_g_index_map = 
      mp_nonlinear_model->get_g_mp_map_indices();
    Teuchos::Array< Teuchos::RCP<const Epetra_Map> > base_g_maps = 
      mp_nonlinear_model->get_g_mp_base_maps();
    mp_g_index_map.push_back(base_g_maps.size());
    base_g_maps.push_back(model->get_x_map());
    Teuchos::RCP<EpetraExt::ModelEvaluator> mp_inverse_solver =
      Teuchos::rcp(new Stokhos::MPInverseModelEvaluator(mp_solver,
							mp_p_index_map,
							mp_g_index_map,
							base_g_maps));

    // Create MP-based SG Quadrature model evaluator to calculate g_sg
    sg_model =
      Teuchos::rcp(new Stokhos::SGQuadMPModelEvaluator(mp_inverse_solver, 
						       sg_comm, 
						       mp_block_map));
  }
  else {
    Teuchos::RCP<EpetraExt::ModelEvaluator> underlying_model;
    if (sg_method == SG_GLOBAL) {
      underlying_model = model;
    } else {
      Piro::Epetra::SolverFactory solverFactory;
      underlying_model = solverFactory.createSolver(piroParams, model);
    }
    sg_model =
      Teuchos::rcp(new Stokhos::SGQuadModelEvaluator(underlying_model));
  }

  // Set up SG nonlinear model
  sg_nonlin_model =
    Teuchos::rcp(new Stokhos::SGModelEvaluator(sg_model, basis, quad, expansion,
					       sg_parallel_data, 
					       sgSolverParams));

  // Set up stochastic parameters
  // One sublist for each stochastic parameter *vector*, and each parameter
  // vector can provide an initial set of expansion coefficients in the basis.
  // This decouples the stochastic parameters from the SG basis allowing e.g.,
  // more stochastic parameters than fundamental r.v.'s in the basis
  // (for correlation) or fewer.
  Teuchos::ParameterList& sgParameters = sgParams.sublist("SG Parameters");
  bool set_initial_params = sgParameters.get("Set Initial SG Parameters", true);
  if (set_initial_params) {
    int num_param_vectors = 
      sgParameters.get("Number of SG Parameter Vectors", 1);
    Teuchos::Array<double> point(basis->dimension(), 1.0);
    Teuchos::Array<double> basis_vals(basis->size());
    basis->evaluateBases(point, basis_vals);
    int idx=0;
    for (int i=0; i<num_param_vectors; i++) {
      std::stringstream ss;
      ss << "SG Parameter Vector " << i;
      Teuchos::ParameterList& pList = sgParameters.sublist(ss.str());
      int p_vec = pList.get("Parameter Vector Index", i);
      
      // Create sg parameter vector
      Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_p =
	sg_nonlin_model->create_p_sg(p_vec);
      
      // Initalize sg parameter vector
      int num_params = sg_p->coefficientMap()->NumMyElements();
      for (int j=0; j<num_params; j++) {
	std::stringstream ss2;
	ss2 << "Parameter " << j << " Initial Expansion Coefficients";
	Teuchos::Array<double> initial_p_vals;
	initial_p_vals = pList.get(ss2.str(),initial_p_vals);
	if (initial_p_vals.size() == 0) {
	  // Default to mean-zero linear expansion, ie, p_j = \xi_j,
	  // by setting term j+1 to 1 (unnormalized)
	  (*sg_p)[idx+1][j] = 1.0 / basis_vals[idx+1];
	}
	else
	  for (Teuchos::Array<double>::size_type l=0; l<initial_p_vals.size(); 
	       l++)
	    (*sg_p)[l][j] = initial_p_vals[l];
	idx++;
      }

      // Set sg parameter vector
      sg_nonlin_model->set_p_sg_init(p_vec, *sg_p);
    }
  }

  // Setup stochastic initial guess
  if (sg_method != SG_NI && sg_method != SG_MPNI) {
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_x = 
      sg_nonlin_model->create_x_sg();
    sg_x->init(0.0);
    if (sg_x->myGID(0))
      (*sg_x)[0] = *(model->get_x_init());
    sg_nonlin_model->set_x_sg_init(*sg_x);
  }

  return sg_nonlin_model;
}

Teuchos::RCP<NOX::Epetra::Observer>
Piro::Epetra::StokhosSolverFactory::
createSGObserver(const Teuchos::RCP<NOX::Epetra::Observer>& noxObserver)
{
  // Set up Observer to call noxObserver for each vector block
  Teuchos::RCP<NOX::Epetra::Observer> sgnoxObserver;

  Teuchos::ParameterList& sgParams = piroParams->sublist("Stochastic Galerkin");
  if (noxObserver != Teuchos::null && sg_method != SG_NI && sg_method != SG_MPNI) {
    int save_moments = sgParams.get("Save Moments",-1);
    sgnoxObserver = 
      Teuchos::rcp(new Piro::Epetra::StokhosNOXObserver(
        noxObserver, basis, 
        sg_nonlin_model->get_overlap_stochastic_map(),
	model->get_x_map(), 
        sg_nonlin_model->get_x_sg_overlap_map(),
        sg_comm, sg_nonlin_model->get_x_sg_importer(), save_moments));
  }

  return sgnoxObserver;
}

Teuchos::RCP<EpetraExt::ModelEvaluator>
Piro::Epetra::StokhosSolverFactory::
createSGSolver(const Teuchos::RCP<EpetraExt::ModelEvaluator>& sg_model,
               const Teuchos::RCP<NOX::Epetra::Observer>& sg_observer)
{
  // Get SG solver type
  std::string solve_type = sgSolverParams->get("SG Solver Algorithm", "Krylov");
  SG_SOLVER solve_method;
  if (solve_type == "Krylov")
    solve_method = SG_KRYLOV;
  else if (solve_type ==  "Gauss-Seidel")
    solve_method = SG_GS;
  else if (solve_type ==  "Jacobi")
    solve_method = SG_JACOBI; 
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
		       std::endl << "Error!  ENAT_SGNOXSolver():  " <<
		       "Invalid Solver Algorithm  " << solve_type << std::endl);

  Teuchos::RCP<EpetraExt::ModelEvaluator> sg_block_solver;
  if (sg_method != SG_NI && sg_method != SG_MPNI) {
    Piro::Epetra::SolverFactory solverFactory;

    Teuchos::RCP<NOX::Epetra::LinearSystem> sg_linsys = Teuchos::null;
    if (solve_method==SG_GS || solve_method==SG_JACOBI) {
      // Create NOX interface
      Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> det_nox_interface = 
         Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(model));

      // Create NOX linear system object
      Teuchos::RCP<const Epetra_Vector> det_u = model->get_x_init();
      Teuchos::RCP<Epetra_Operator> det_A = model->create_W();
      Teuchos::RCP<NOX::Epetra::Interface::Required> det_iReq = det_nox_interface;
      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> det_iJac = det_nox_interface;
      //Teuchos::ParameterList det_printParams;
      Teuchos::ParameterList& noxParams = piroParams->sublist("NOX");
      Teuchos::ParameterList& det_printParams = noxParams.sublist("Printing");
      Teuchos::ParameterList& printParams = noxParams.sublist("Printing");
      Teuchos::ParameterList& newtonParams = 
	noxParams.sublist("Direction").sublist("Newton");
      Teuchos::ParameterList& det_lsParams = 
	newtonParams.sublist("Stratimikos Linear Solver");
      
      Teuchos::RCP<NOX::Epetra::LinearSystem> det_linsys = 
      	Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(
		     det_printParams, det_lsParams, det_iJac, 
		     det_A, *det_u));

     // Sublist for linear solver for the Newton method
      //Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
      Teuchos::ParameterList& sgjacobiParams = 
                newtonParams.sublist("Linear Solver");
     // Create NOX interface
      Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> nox_interface =
        Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(sg_model));
      Teuchos::RCP<const Epetra_Map> base_map = model->get_x_map();
      Teuchos::RCP<const Epetra_Map> sg_map = sg_model->get_x_map();
      Teuchos::RCP<Epetra_Operator> A = sg_model->create_W();
      Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = nox_interface;
      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = nox_interface; 
    
      if (solve_method==SG_GS) {
       sgjacobiParams.sublist("Deterministic Solver Parameters") = det_lsParams;
       
       sg_linsys =
	 Teuchos::rcp(new NOX::Epetra::LinearSystemSGGS(
	  	       printParams, sgjacobiParams, det_linsys, iReq, iJac, 
	 	       basis, sg_parallel_data, A, base_map, sg_map));
      }

      else if (solve_method==SG_JACOBI) {
       sgjacobiParams.sublist("Deterministic Solver Parameters") = det_lsParams;
       Teuchos::ParameterList& jacobiOpParams =
	 sgjacobiParams.sublist("Jacobi SG Operator");
       jacobiOpParams.set("Only Use Linear Terms", true);
       sg_linsys =
	Teuchos::rcp(new NOX::Epetra::LinearSystemSGJacobi(
		       printParams, sgjacobiParams, det_linsys, iReq, iJac, 
		       basis, sg_parallel_data, A, base_map, sg_map));
      }

      {
        const std::string token = "Linear System";
        solverFactory.setNOXLinearSystemProvider(token, sg_linsys);
        Teuchos::sublist(piroParams, token)->set("Type", token);
      }
    }

    {
      const std::string token = "NOX Observer";
      solverFactory.setNOXObserverProvider(token, sg_observer);
      Teuchos::sublist(piroParams, token)->set("Type", token);
    }

    // Will find preconditioner for Matrix-Free method
    sg_block_solver = solverFactory.createSolver(piroParams, sg_model);
  }
  else 
    sg_block_solver = sg_model;

  return sg_block_solver;
}

Teuchos::RCP<Stokhos::SGInverseModelEvaluator>
Piro::Epetra::StokhosSolverFactory::
createSGSolverAdapter(const Teuchos::RCP<EpetraExt::ModelEvaluator>& sg_solver)
{
  // Create SG Inverse model evaluator
  Teuchos::Array<int> sg_p_index_map = sg_nonlin_model->get_p_sg_map_indices();
  Teuchos::Array<int> sg_g_index_map = sg_nonlin_model->get_g_sg_map_indices();
  Teuchos::Array< Teuchos::RCP<const Epetra_Map> > base_g_maps = 
    sg_nonlin_model->get_g_sg_base_maps();
  // Add sg_u response function supplied by Piro::Epetra::NOXSolver
  if (sg_method != SG_NI && sg_method != SG_MPNI && 
      piroParams->get<std::string>("Solver Type") == "NOX") {
    sg_g_index_map.push_back(base_g_maps.size());
    base_g_maps.push_back(model->get_x_map());
  }
  Teuchos::RCP<Stokhos::SGInverseModelEvaluator> sg_adapter = 
    Teuchos::rcp(new Stokhos::SGInverseModelEvaluator(sg_solver, 
						      sg_p_index_map,
						      sg_g_index_map,
						      base_g_maps));

  return sg_adapter;
}

Teuchos::RCP<EpetraExt::ModelEvaluator>
Piro::Epetra::StokhosSolverFactory::
createRSModel(const Teuchos::RCP<EpetraExt::ModelEvaluator>& sg_model)
{
  // Create ResponseStatistic model evaluator
  Teuchos::Array< Teuchos::RCP<const Epetra_Map> > base_g_maps = 
    sg_nonlin_model->get_g_sg_base_maps();
  // Add sg_u response function supplied by Piro::Epetra::NOXSolver
  if (sg_method != SG_NI && sg_method != SG_MPNI && 
      piroParams->get("Solver Type", "NOX") == "NOX") {
    base_g_maps.push_back(model->get_x_map());
  }
  Teuchos::RCP<const Epetra_BlockMap> block_map =
    sg_nonlin_model->get_overlap_stochastic_map();
  Teuchos::RCP<EpetraExt::ModelEvaluator> rs_model = 
    Teuchos::rcp(new Stokhos::ResponseStatisticModelEvaluator(
		   sg_model, base_g_maps, basis, sg_comm, block_map));

  return rs_model;
}

Teuchos::RCP<const Epetra_Comm>
Piro::Epetra::StokhosSolverFactory::
getSpatialComm() const
{
  return Stokhos::getSpatialComm(sg_comm);
}

Teuchos::RCP<const Epetra_Comm>
Piro::Epetra::StokhosSolverFactory::
getStochasticComm() const
{
  return Stokhos::getStochasticComm(sg_comm);
}

Teuchos::RCP<const EpetraExt::MultiComm>
Piro::Epetra::StokhosSolverFactory::
getGlobalMultiComm() const
{
  return sg_comm;
}

Teuchos::RCP<const Teuchos::ParameterList>
Piro::Epetra::StokhosSolverFactory::getValidSGParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
     rcp(new Teuchos::ParameterList("ValidSGParams"));;
  validPL->sublist("SG Parameters", false, "");
  validPL->sublist("SG Solver Parameters", false, "");
  validPL->sublist("MP Solver Parameters", false, "");
  validPL->sublist("Basis", false, "");
  validPL->sublist("Pseudospectral Operator", false, "");
  validPL->sublist("Expansion", false, "");
  validPL->sublist("Quadrature", false, "");
  validPL->set<std::string>("SG Method", "","");
  validPL->set<std::string>("Triple Product Size", "","");
  validPL->set<bool>("Rebalance Stochastic Graph", false, "");
  validPL->set<int>("Save Moments", -1, "Set to 2 for Mean and Variance. Default writes Coeffs");
  validPL->set<int>("Number of Spatial Processors", -1, "");
  validPL->sublist("Isorropia", false, "");
  validPL->sublist("Response KL", false, "");

  return validPL;
}

