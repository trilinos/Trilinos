// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include <string>
#include <fenv.h>

#include "MockModelEval_A_Tpetra.hpp"
#include "MockModelEval_B_Tpetra.hpp"
#include "MockModelEval_B_Tpetra_2_parameters.hpp"
#include "MassSpringDamperModel.hpp"
//#include "ObserveSolution_Epetra.hpp"

#include "Piro_SolverFactory.hpp"
#include "Piro_ProviderHelpers.hpp"

#include "Piro_PerformAnalysis.hpp"

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Piro_StratimikosUtils.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Tpetra_Core.hpp"

#ifdef HAVE_PIRO_IFPACK2
#include "Thyra_Ifpack2PreconditionerFactory.hpp"
#endif

#ifdef HAVE_PIRO_MUELU
#include "Stratimikos_MueLuHelpers.hpp"
#endif

#ifdef HAVE_PIRO_TEMPUS
#include "Piro_TempusSolver.hpp"
#include "Tempus_StepperFactory.hpp"
#endif

const Teuchos::RCP<Piro::TempusSolver<double> > solverNew(
    const Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > &thyraModel,
    const Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > &thyraAdjointModel,
    double finalTime, 
    const std::string sens_method_string)
{
  Teuchos::RCP<Teuchos::ParameterList> analysisPL =
    Teuchos::rcp(new Teuchos::ParameterList("Analysis"));
  auto tempusPL = analysisPL->sublist("Tempus");
  analysisPL->sublist("Tempus").sublist("Tempus", false, "");
  /*
  analysisPL->sublist("Tempus").sublist("Albany Time Step Control Options", false, "");
  analysisPL->sublist("Tempus").sublist("Albany Time Step Control Options", false, "").set<double>("Initial Time", 0.0, "");
  analysisPL->sublist("Tempus").sublist("Albany Time Step Control Options", false, "").set<double>("Initial Time Step", 0.1, "");
  analysisPL->sublist("Tempus").sublist("Albany Time Step Control Options", false, "").set<double>("Minimum Time Step", 0.1, "");
  analysisPL->sublist("Tempus").sublist("Albany Time Step Control Options", false, "").set<double>("Maximum Time Step", 0.1, "");
  analysisPL->sublist("Tempus").sublist("Albany Time Step Control Options", false, "").set<double>("Final Time", 1.0, "");
  analysisPL->sublist("Tempus").sublist("Albany Time Step Control Options", false, "").set<double>("Reduction Factor", 1.0, "");
  analysisPL->sublist("Tempus").sublist("Albany Time Step Control Options", false, "").set<double>("Amplification Factor", 1.0, "");
  */
  analysisPL->sublist("Tempus").sublist("Stratimikos", false, "");
  analysisPL->sublist("Tempus").sublist("NonLinear Solver", false, "");
  //analysisPL->sublist("Tempus").set<std::string>("Verbosity Level", "", "");
  analysisPL->sublist("Tempus").set<bool>("Lump Mass Matrix", false, "Boolean to tell code whether to lump mass matrix");
  analysisPL->sublist("Tempus").set<bool>("Invert Mass Matrix", true, "Boolean to tell code whether or not to invert mass matrix");
  analysisPL->sublist("Tempus").set<bool>("Constant Mass Matrix", false, "Boolean to tell code if mass matrix is constant in time");
  analysisPL->sublist("Tempus").set<bool>("Abort on Failure", true, "");
  analysisPL->sublist("Tempus").set("Integrator Name", "Tempus Integrator");
  analysisPL->sublist("Tempus").sublist("Tempus Integrator").set("Integrator Type", "Integrator Basic");
  analysisPL->sublist("Tempus").sublist("Tempus Integrator").set("Stepper Name", "Tempus Stepper");
  analysisPL->sublist("Tempus").sublist("Tempus Integrator").sublist("Solution History").set("Storage Type", "Unlimited");
  analysisPL->sublist("Tempus").sublist("Tempus Integrator").sublist("Solution History").set("Storage Limit", 20);
  analysisPL->sublist("Tempus").sublist("Tempus Integrator").sublist("Time Step Control").set("Initial Time", 0.0);
  analysisPL->sublist("Tempus").sublist("Tempus Integrator").sublist("Time Step Control").set("Final Time", finalTime);
  analysisPL->sublist("Tempus").sublist("Tempus Stepper").set("Stepper Type", "Backward Euler");
  analysisPL->sublist("Tempus").sublist("Tempus Stepper").set("Zero Initial Guess", false);
  analysisPL->sublist("Tempus").sublist("Tempus Stepper").set("Solver Name", "Demo Solver");
  analysisPL->sublist("Tempus").sublist("Tempus Stepper").sublist("Demo Solver").sublist("NOX").sublist("Direction").set("Method","Newton");
  analysisPL->sublist("Tempus").sublist("Sensitivities", false, "");
  analysisPL->set("Sensitivity Method", "Adjoint");

  return Teuchos::rcp(new Piro::TempusSolver<double>(analysisPL, thyraModel, thyraAdjointModel));
}


#include "Piro_ConfigDefs.hpp"

int main(int argc, char *argv[]) {

  int status=0; // 0 = pass, failures are incremented
  int overall_status=0; // 0 = pass, failures are incremented over multiple tests
  bool success=true;

  // Initialize MPI
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  int Proc=mpiSession.getRank();

  auto appComm = Tpetra::getDefaultComm();

  using Teuchos::RCP;
  using Teuchos::rcp;

  std::string inputFile;
  bool doAll = (argc==1);
  if (argc>1) doAll = !strcmp(argv[1],"-v");

  //Piro::SolverFactory solverFactory;

  for (int iTest=0; iTest<3; iTest++) {

    if (doAll) {
      switch (iTest) {
       case 0: inputFile="input_Analysis_ROL_ReducedSpace_Transient.xml"; break;
       case 1: inputFile="input_Analysis_ROL_ReducedSpace_Transient_2_parameters.xml"; break;
       case 2: inputFile="input_Analysis_ROL_ReducedSpace_Transient_MSD.xml"; break;
       default : std::cout << "iTest logic error " << std::endl; exit(-1);
      }
    }
    else {
      inputFile=argv[1];
      iTest = 999;
    }

    try {

      std::vector<std::string> mockModels = {"MockModelEval_B_Tpetra", "MockModelEval_B_Tpetra_2_parameters", "MassSpringDamperModel"};
      // MockModelEval_A_Tpetra is very slow so it is currently disabled.

      for (auto mockModel : mockModels) {

        if (mockModel=="MockModelEval_B_Tpetra_2_parameters" && inputFile != "input_Analysis_ROL_ReducedSpace_Transient_2_parameters.xml") {
          continue;
        }
        if (mockModel=="MockModelEval_A_Tpetra" && inputFile != "input_Analysis_ROL_ReducedSpace_Transient.xml") {
          continue;
        }
        if (mockModel=="MockModelEval_B_Tpetra" && inputFile != "input_Analysis_ROL_ReducedSpace_Transient.xml") {
          continue;
        }
        if (mockModel=="MassSpringDamperModel" && inputFile != "input_Analysis_ROL_ReducedSpace_Transient_MSD.xml" && inputFile != "input_Analysis_ROL_Tempus_Transient_MSD.xml" && inputFile != "input_Analysis_ROL_ReducedSpace_Transient_Integrated_Response_MSD.xml") {
          continue;
        }
        if (mockModel=="MassSpringDamperModel" && appComm->getSize() > 1) {
          continue;
        }

        // BEGIN Builder
        const RCP<Teuchos::ParameterList> appParams = rcp(new Teuchos::ParameterList("Application Parameters"));
        Teuchos::updateParametersFromXmlFile(inputFile, Teuchos::ptr(appParams.get()));

        const RCP<Teuchos::ParameterList>  probParams = Teuchos::sublist(appParams,"Problem");
        const RCP<Teuchos::ParameterList>  piroParams = Teuchos::sublist(appParams,"Piro");
 
        bool boundConstrained = piroParams->sublist("Analysis").sublist("ROL").get<bool>("Bound Constrained");
     
        // Create (1) a Model Evaluator and (2) a ParameterList
        std::string modelName;
        bool adjoint = (piroParams->get("Sensitivity Method", "Forward") == "Adjoint");
        bool explicitAdjointME = adjoint && piroParams->get("Explicit Adjoint Model Evaluator", false);
        RCP<Thyra::ModelEvaluator<double>> model, adjointModel(Teuchos::null);

        int num_parameters = piroParams->sublist("Analysis").sublist("ROL").get<int>("Number Of Parameters", 1);
        std::vector<int> p_indices(num_parameters);

        for(int i=0; i<num_parameters; ++i) {
          std::ostringstream ss; ss << "Parameter Vector Index " << i;
          p_indices[i] = piroParams->sublist("Analysis").sublist("ROL").get<int>(ss.str(), i);
        }

        if (mockModel=="MockModelEval_A_Tpetra") {
          if(boundConstrained) {
            RCP<Thyra::ModelEvaluator<double>> model_tmp = rcp(new MockModelEval_A_Tpetra(appComm,false,probParams,true));
            model = rcp(new Piro::ProductModelEvaluator<double>(model_tmp,p_indices));
            if(explicitAdjointME) {
              RCP<Thyra::ModelEvaluator<double>> adjointModel_tmp = rcp(new MockModelEval_A_Tpetra(appComm,true));
              adjointModel = rcp(new Piro::ProductModelEvaluator<double>(adjointModel_tmp,p_indices));
            }
            modelName = "A";
          } else   // optimization of problem A often diverges when the parameters are not constrained
            continue;
        }
        else if (mockModel=="MockModelEval_B_Tpetra") { 
          RCP<Thyra::ModelEvaluator<double>> model_tmp = rcp(new MockModelEval_B_Tpetra(appComm,false,probParams,true));
          model = rcp(new Piro::ProductModelEvaluator<double>(model_tmp,p_indices));
          if(explicitAdjointME) {
            RCP<Thyra::ModelEvaluator<double>> adjointModel_tmp = rcp(new MockModelEval_B_Tpetra(appComm,true));
            adjointModel = rcp(new Piro::ProductModelEvaluator<double>(adjointModel_tmp,p_indices));
          }
          modelName = "B";
        }
        else if (mockModel=="MockModelEval_B_Tpetra_2_parameters") {
          RCP<Thyra::ModelEvaluator<double>> model_tmp = rcp(new MockModelEval_B_Tpetra_2_parameters(appComm,false,probParams,true));
          model = rcp(new Piro::ProductModelEvaluator<double>(model_tmp,p_indices));
          if(explicitAdjointME) {
            RCP<Thyra::ModelEvaluator<double>> adjointModel_tmp = rcp(new MockModelEval_B_Tpetra_2_parameters(appComm,true));
            adjointModel = rcp(new Piro::ProductModelEvaluator<double>(adjointModel_tmp,p_indices));
          }
          modelName = "B_2";
        }
        else if (mockModel=="MassSpringDamperModel") {
          auto comm = appComm->createSubcommunicator(std::vector<int>({0}));
          RCP<Thyra::ModelEvaluator<double>> model_tmp = rcp(new MassSpringDamperModel(comm,false,probParams,true,probParams->get<bool>("Response Depends On Solution Time Derivative", false)));
          model = rcp(new Piro::ProductModelEvaluator<double>(model_tmp,p_indices));
          if(explicitAdjointME) {
            RCP<Thyra::ModelEvaluator<double>> adjointModel_tmp = rcp(new MassSpringDamperModel(comm,true,Teuchos::null,false,probParams->get<bool>("Response Depends On Solution Time Derivative", false)));
            adjointModel = rcp(new Piro::ProductModelEvaluator<double>(adjointModel_tmp,p_indices));
          }
          modelName = "MSD";
        }
        if (Proc==0)
          std::cout << "=======================================================================================================\n"
                    << "======  Solving Problem " << modelName << " with input file "<< iTest <<": "<< inputFile <<"\n"
                    << "=======================================================================================================\n"
            << std::endl;
        

        Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

  #ifdef HAVE_PIRO_IFPACK2
        typedef Thyra::PreconditionerFactoryBase<double>              Base;
        typedef Thyra::Ifpack2PreconditionerFactory<Tpetra_CrsMatrix> Impl;
        linearSolverBuilder.setPreconditioningStrategyFactory(
            Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");
  #endif

  #ifdef HAVE_PIRO_MUELU
        using scalar_type = Tpetra::CrsMatrix<>::scalar_type;
        using local_ordinal_type = Tpetra::CrsMatrix<>::local_ordinal_type;
        using global_ordinal_type = Tpetra::CrsMatrix<>::global_ordinal_type;
        using node_type = Tpetra::CrsMatrix<>::node_type;
        Stratimikos::enableMueLu<scalar_type, local_ordinal_type, global_ordinal_type, node_type>(linearSolverBuilder);
  #endif

        const Teuchos::RCP<Teuchos::ParameterList> stratList = Piro::extractStratimikosParams(piroParams);
        linearSolverBuilder.setParameterList(stratList);


        const RCP<Thyra::LinearOpWithSolveFactoryBase<double>> lowsFactory =
            createLinearSolveStrategy(linearSolverBuilder);

        RCP<Thyra::ModelEvaluator<double>> modelWithSolve = rcp(new Thyra::DefaultModelEvaluatorWithSolveFactory<double>(
            model, lowsFactory));
        RCP<Thyra::ModelEvaluator<double>> adjointModelWithSolve(Teuchos::null);
        if(Teuchos::nonnull(adjointModel))
          adjointModelWithSolve= rcp(new Thyra::DefaultModelEvaluatorWithSolveFactory<double>(adjointModel, lowsFactory));

        const RCP<Thyra::ModelEvaluatorDefaultBase<double>> piro = solverNew(Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDefaultBase<double>>(modelWithSolve), Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDefaultBase<double>>(adjointModelWithSolve), 10., "Adjoint");

        // Call the analysis routine
        RCP<Thyra::VectorBase<double>> p;
        status = Piro::PerformAnalysis(*piro, *piroParams, p);

        if (Teuchos::nonnull(p)) { //p might be null if the package ROL is not enabled
          Thyra::DetachedVectorView<double> p_view(p);
          double p_exact[2];
          if (mockModel=="MockModelEval_A_Tpetra") {
            p_exact[0] = 1;
            p_exact[1] = 3;
          }
          if (mockModel=="MockModelEval_B_Tpetra") {
            p_exact[0] = 6;
            p_exact[1] = 4;
          }
          if (mockModel=="MockModelEval_B_Tpetra_2_parameters") {
            p_exact[0] = 4;
            p_exact[1] = 6;
          }
          if (mockModel=="MassSpringDamperModel") {
            p_exact[0] = 1.;
            p_exact[1] = 0.5;
          }
          double tol = 1e-5;

          double l2_diff = std::sqrt(std::pow(p_view(0)-p_exact[0],2) + std::pow(p_view(1)-p_exact[1],2));
          if(l2_diff > tol) {
            status+=100;
            if (Proc==0) {
              std::cout << "\nPiro_AnalysisDrvier:  Optimum parameter values are: {"
                  <<  p_exact[0] << ", " << p_exact[1] << "}, but computed values are: {"
                  <<  p_view(0) << ", " << p_view(1) << "}." <<
                  "\n                      Difference in l2 norm: " << l2_diff << " > tol: " << tol <<   std::endl;
            }
          }
        }
      }
    }
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
    if (!success) status+=1000;

    overall_status += status;
  }  // End loop over tests

  if (Proc==0) {
    if (overall_status==0)
      std::cout << "\nTEST PASSED\n" << std::endl;
    else
      std::cout << "\nTEST Failed: " << overall_status << "\n" << std::endl;
  }

  return status;
}
