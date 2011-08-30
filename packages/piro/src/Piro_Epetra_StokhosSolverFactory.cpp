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

#include "Piro_Epetra_StokhosSolver.hpp"
#include "Piro_Epetra_Factory.hpp"
#include "Piro_ValidPiroParameters.hpp"

#include "Stokhos.hpp"
#include "Stokhos_Epetra.hpp"

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
  //piroParams->validateParameters(*Piro::getValidPiroParameters(),0);

  // Validate parameters
  Teuchos::ParameterList& sgParams =
    piroParams->sublist("Stochastic Galerkin");
  sgParams.validateParameters(*getValidSGParameters(),0);

  sgSolverParams = 
    Teuchos::rcp(&(sgParams.sublist("SG Solver Parameters")),false);

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
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
		       std::endl << "Error!  ENAT_SGNOXSolver():  " <<
		       "Invalid SG Method  " << sg_type << std::endl);
  
  // Create SG basis
  basis = Stokhos::BasisFactory<int,double>::create(sgParams);
  if (globalComm->MyPID()==0) 
    std::cout << "Basis size = " << basis->size() << std::endl;

  // Create SG Quadrature
  Teuchos::ParameterList& expParams = sgParams.sublist("Expansion");
  std::string exp_type = expParams.get("Type", "Quadrature");
  if (exp_type == "Quadrature" || 
      sg_method == SG_GLOBAL ||
      sg_method == SG_NI ||
      sg_method == SG_MPNI) {
    quad = Stokhos::QuadratureFactory<int,double>::create(sgParams);
    if (globalComm->MyPID()==0) 
      std::cout << "Quadrature size = " << quad->size() << std::endl;
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

Piro::Epetra::StokhosSolverFactory::~StokhosSolverFactory()
{
  // Get rid of circular dependencies
  piroParams->set("Interface", Teuchos::null);
  piroParams->set("Linear System", Teuchos::null);
}

Teuchos::RCP<Stokhos::SGModelEvaluator>
Piro::Epetra::StokhosSolverFactory::
createSGModel(const Teuchos::RCP<EpetraExt::ModelEvaluator>& model_,
	      const Teuchos::RCP<NOX::Epetra::Observer>& noxObserver)
{
  Teuchos::ParameterList& sgParams =
    piroParams->sublist("Stochastic Galerkin");
  Teuchos::ParameterList& sg_basisParams = sgParams.sublist("Basis");
  int dim = sg_basisParams.get<int>("Dimension");

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

      Teuchos::RCP<const Teuchos::ParameterList> ortho_params = 
	Teuchos::rcp(new Teuchos::ParameterList);
      noxstratlsParams.sublist("Stratimikos").sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist("GCRODR").set("Orthogonalization Parameters", ortho_params);


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

      piroParams->set("Interface", nox_interface);
      piroParams->set("Linear System", linsys);
    }

    // Create solver to map p -> g
    Teuchos::RCP<EpetraExt::ModelEvaluator> mp_solver =
      Piro::Epetra::Factory::createSolver(piroParams, mp_nonlinear_model);

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
    if (sg_method == SG_GLOBAL)
      underlying_model = model;
    else 
      underlying_model =
	Piro::Epetra::Factory::createSolver(piroParams, model);
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
  int num_param_vectors = sgParameters.get("Number of SG Parameter Vectors", 1);
  for (int i=0; i<num_param_vectors; i++) {
    std::stringstream ss;
    ss << "SG Parameter Vector " << i;
    Teuchos::ParameterList& pList = sgParameters.sublist(ss.str());
    int p_vec = pList.get("Parameter Vector Index", 0);
    
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
      if (initial_p_vals.size() == 0 && dim == num_params) {
	// Default to mean-zero linear expansion, ie, p_i = \xi_i
	// This only makes sense if the stochastic dimension is equal to the
	// number of parameters
	sg_p->term(j,0)[j] = 0.0;
	sg_p->term(j,1)[j] = 1.0;  // Set order 1 coeff to 1 for this RV
      }
      else
	for (Teuchos::Array<double>::size_type l=0; l<initial_p_vals.size(); 
	     l++)
	  (*sg_p)[l][j] = initial_p_vals[l];
    }

    // Set sg parameter vector
    sg_nonlin_model->set_p_sg_init(p_vec, *sg_p);
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

  // Set up Observer to call noxObserver for each vector block
  Teuchos::RCP<NOX::Epetra::Observer> sgnoxObserver;
  if (noxObserver != Teuchos::null && sg_method != SG_NI && sg_method != SG_MPNI) {
    int save_moments = sgParams.get("Save Moments",-1);
    sgnoxObserver = 
      Teuchos::rcp(new Piro::Epetra::StokhosNOXObserver(
        noxObserver, basis, 
        sg_nonlin_model->get_overlap_stochastic_map(),
	model->get_x_map(), 
        sg_nonlin_model->get_x_sg_overlap_map(),
        sg_comm, sg_nonlin_model->get_x_sg_importer(), save_moments));
    piroParams->set("Observer", sgnoxObserver);
  }

  return sg_nonlin_model;
}

Teuchos::RCP<EpetraExt::ModelEvaluator>
Piro::Epetra::StokhosSolverFactory::
createSGSolver(const Teuchos::RCP<EpetraExt::ModelEvaluator>& sg_model)
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
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
		       std::endl << "Error!  ENAT_SGNOXSolver():  " <<
		       "Invalid Solver Algorithm  " << solve_type << std::endl);

  Teuchos::RCP<EpetraExt::ModelEvaluator> sg_block_solver;
  if (sg_method != SG_NI && sg_method != SG_MPNI) {
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

      piroParams->set("Linear System", sg_linsys);
    }

    // Will find preconditioner for Matrix-Free method
    sg_block_solver = 
      Piro::Epetra::Factory::createSolver(piroParams, sg_model);
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
  validPL->sublist("Expansion", false, "");
  validPL->sublist("Quadrature", false, "");
  validPL->set<std::string>("SG Method", "","");
  validPL->set<std::string>("Triple Product Size", "","");
  validPL->set<bool>("Rebalance Stochastic Graph", false, "");
  validPL->set<int>("Save Moments", -1, "Set to 2 for Mean and Variance. Default writes Coeffs");
  validPL->set<int>("Number of Spatial Processors", -1, "");
  validPL->sublist("Isorropia", false, "");

  return validPL;
}

