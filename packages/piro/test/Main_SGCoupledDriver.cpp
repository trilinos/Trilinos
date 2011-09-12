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

#include <iostream>
#include <string>

#include "MockModelEval_D.hpp"

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

// Piro solvers
#include "Piro_Epetra_StokhosSolver.hpp"
#include "Piro_Epetra_NECoupledModelEvaluator.hpp"

int main(int argc, char *argv[]) {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  int status=0; // 0 = pass, failures are incremented
  int overall_status=0; // 0 = pass, failures are incremented over multiple tests
  bool success=true;

  // Initialize MPI 
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  int Proc=mpiSession.getRank();

  // Create a communicator for Epetra objects
  RCP<const Epetra_Comm> globalComm;
#ifdef HAVE_MPI
  globalComm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  globalComm = rcp(new Epetra_SerialComm);
#endif

  std::string problem1_filename = "input_problem1.xml";
  std::string problem2_filename = "input_problem2.xml";
  std::string coupled_filename = "input_coupled.xml";

  try {
    
    // Setup stochastic coupled problem to get spatial comm's
    RCP<ParameterList> coupledParams = 
      Teuchos::getParametersFromXmlFile(coupled_filename);
    RCP<Piro::Epetra::StokhosSolver> coupledSolver =
      rcp(new Piro::Epetra::StokhosSolver(coupledParams, globalComm));
    RCP<const Epetra_Comm> app_comm = coupledSolver->getSpatialComm();

    // Setup problem 1
    RCP<ParameterList> piroParams1 = 
      Teuchos::getParametersFromXmlFile(problem1_filename);
    RCP<EpetraExt::ModelEvaluator> model1 = rcp(new MockModelEval_D(app_comm));

    // Setup problem 2
    RCP<ParameterList> piroParams2 = 
      Teuchos::getParametersFromXmlFile(problem2_filename);
    RCP<EpetraExt::ModelEvaluator> model2 = rcp(new MockModelEval_D(app_comm));

    // Setup coupled model
    RCP<EpetraExt::ModelEvaluator> coupledModel =
      rcp(new Piro::Epetra::NECoupledModelEvaluator(model1, model2,
						    piroParams1, piroParams2,
						    coupledParams, globalComm));

    // Setup solver
    coupledSolver->setup(coupledModel);

    Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly> x_sg_init =
      coupledSolver->get_x_sg_init();
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> x_sg_init_new =
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(*x_sg_init));
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis =
      coupledSolver->getBasis();
    for (int i=0; i<basis->dimension(); i++)
      (*x_sg_init_new)[i+1].PutScalar(1.0);
    coupledSolver->set_x_sg_init(*x_sg_init_new);
    
    // Solve coupled system
    EpetraExt::ModelEvaluator::InArgs inArgs = coupledSolver->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs outArgs = coupledSolver->createOutArgs();
    for (int i=0; i<inArgs.Np(); i++)
      if (inArgs.supports(EpetraExt::ModelEvaluator::IN_ARG_p_sg, i))
	inArgs.set_p_sg(i, coupledSolver->get_p_sg_init(i));
    for (int i=0; i<outArgs.Ng(); i++) 
      if (outArgs.supports(EpetraExt::ModelEvaluator::OUT_ARG_g_sg, i)) {
	RCP<Stokhos::EpetraVectorOrthogPoly> g_sg = 
	  coupledSolver->create_g_sg(i);
	outArgs.set_g_sg(i, g_sg);
      }
    coupledSolver->evalModel(inArgs, outArgs);

    // Print results
    for (int i=0; i<outArgs.Ng(); i++) {
      if (outArgs.supports(EpetraExt::ModelEvaluator::OUT_ARG_g_sg, i)) {
	Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> g_sg = 
	  outArgs.get_g_sg(i);
	if (g_sg != Teuchos::null) {
	  Epetra_Vector g_mean(*(coupledSolver->get_g_map(i)));
	  Epetra_Vector g_std_dev(*(coupledSolver->get_g_map(i)));
	  g_sg->computeMean(g_mean);
	  g_sg->computeStandardDeviation(g_std_dev);
	  std::cout.precision(12);
	  std::cout << "Response " << i << " Mean =      " << std::endl 
		    << g_mean << std::endl;
	  std::cout << "Response " << i << " Std. Dev. = " << std::endl 
		    << g_std_dev << std::endl;
	  std::cout << "Response vector " << i << ":" << std::endl
		    << *(outArgs.get_g_sg(i)) << std::endl;
	}
      }
    }

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  if (!success) status+=1000;

  overall_status += status;


  if (Proc==0) {
    if (overall_status==0) 
      cout << "\nTEST PASSED\n" << endl;
    else 
      cout << "\nTEST Failed: " << overall_status << "\n" << endl;
  }

  return status;
}
