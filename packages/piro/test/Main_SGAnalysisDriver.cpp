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

#include "MockModelEval_C.hpp"

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

// Piro solvers
#include "Piro_Epetra_StokhosSolverFactory.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Piro_PerformAnalysis.hpp"
#include "Thyra_VectorBase.hpp"

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

  std::string xml_filename = "input_SGAnalysis.xml";

  try {

    // Set up application parameters
    RCP<ParameterList> appParams = 
      Teuchos::getParametersFromXmlFile(xml_filename);
    
    // Create stochastic Galerkin solver factory
    RCP<ParameterList> piroParams = 
      rcp(&(appParams->sublist("Piro")),false);
    Piro::Epetra::StokhosSolverFactory sg_solver_factory(piroParams, 
							 globalComm);

    // Get comm for spatial problem
    RCP<const Epetra_Comm> app_comm = sg_solver_factory.getSpatialComm();
    
    // Create application model evaluator
    RCP<EpetraExt::ModelEvaluator> model = rcp(new MockModelEval_C(app_comm));

    // Setup rest of solver
    RCP<Stokhos::SGModelEvaluator> sg_model = 
      sg_solver_factory.createSGModel(model);
    RCP<EpetraExt::ModelEvaluator> sg_solver =
      sg_solver_factory.createSGSolver(sg_model);
    RCP<EpetraExt::ModelEvaluator> rs_model =
      sg_solver_factory.createRSModel(sg_solver);

    Thyra::EpetraModelEvaluator rs_model_thyra;
    rs_model_thyra.initialize(rs_model, Teuchos::null);

    RCP< ::Thyra::VectorBase<double> > p;
    ParameterList& analysisParams = piroParams->sublist("Analysis");
    int status = Piro::PerformAnalysis(rs_model_thyra, analysisParams, p); 

    std::cout << *p;

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
