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
#include <cmath>

#include "MockModelEval_A_Tpetra.hpp"
#include "MockModelEval_B_Tpetra.hpp"
#include "MockModelEval_B_Tpetra_2_parameters.hpp"
#include "MassSpringDamperModel.hpp"

#include "Piro_SolverFactory.hpp"
#include "Piro_ProviderHelpers.hpp"

#include "Piro_PerformAnalysis.hpp"
#include "Piro_ROL_ObserverBase.hpp"

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

// ROL observer recording the observed solutions, to check that the transient
// optimization calls the observer with the state at the final time.
class RecordingROLObserver : public Piro::ROL_ObserverBase<double> {
public:
  using Piro::ROL_ObserverBase<double>::observeSolution;

  void observeSolution(double stamp, const Thyra::VectorBase<double>& solution,
                       const Teuchos::Ptr<const Thyra::MultiVectorBase<double>>& /*solution_dxdp*/,
                       const Teuchos::Ptr<const Thyra::VectorBase<double>>& /*solutionDot*/,
                       const Teuchos::Ptr<const Thyra::VectorBase<double>>& /*solutionDotDot*/) override {
    ++numSolutions;
    lastStamp = stamp;
    lastSolutionNorm = Thyra::norm_2(solution);
  }

  void observeResponse(int /*iter*/) override { ++numResponses; }

  void parametersChanged() override { ++numParametersChanged; }

  int numSolutions = 0;
  int numResponses = 0;
  int numParametersChanged = 0;
  double lastStamp = -1;
  double lastSolutionNorm = 0;
};

const Teuchos::RCP<Piro::TempusSolver<double> > solverNew(
    const Teuchos::RCP<Teuchos::ParameterList>  piroParams,
    const Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > &thyraModel,
    const Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > &thyraAdjointModel)
{
  return Teuchos::rcp(new Piro::TempusSolver<double>(piroParams, thyraModel, thyraAdjointModel));
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

        if (mockModel=="MockModelEval_B_Tpetra_2_parameters" && inputFile.find("_2_parameters") == std::string::npos) {
          continue;
        }
        if (mockModel=="MockModelEval_A_Tpetra" && inputFile != "input_Analysis_ROL_ReducedSpace_Transient.xml") {
          continue;
        }
        if (mockModel=="MockModelEval_B_Tpetra" && inputFile != "input_Analysis_ROL_ReducedSpace_Transient.xml") {
          continue;
        }
        if (mockModel=="MassSpringDamperModel" && inputFile.find("_MSD") == std::string::npos) {
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
            RCP<Thyra::ModelEvaluator<double>> model_tmp = rcp(new MockModelEval_A_Tpetra(appComm,2,false,probParams,true));
            model = rcp(new Piro::ProductModelEvaluator<double>(model_tmp,p_indices));
            if(explicitAdjointME) {
              RCP<Thyra::ModelEvaluator<double>> adjointModel_tmp = rcp(new MockModelEval_A_Tpetra(appComm,2,true));
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
            RCP<Thyra::ModelEvaluator<double>> adjointModel_tmp = rcp(new MassSpringDamperModel(comm,true,probParams,false,probParams->get<bool>("Response Depends On Solution Time Derivative", false)));
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

        const RCP<Thyra::ModelEvaluatorDefaultBase<double>> piro = 
          solverNew(
            piroParams,
            Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDefaultBase<double>>(modelWithSolve), 
            Teuchos::rcp_dynamic_cast<Thyra::ModelEvaluatorDefaultBase<double>>(adjointModelWithSolve));

        // Call the analysis routine
        RCP<Thyra::VectorBase<double>> p;
        const auto observer = Teuchos::rcp(new RecordingROLObserver());
        status = Piro::PerformAnalysis(*piro, *piroParams, p, observer);

        // The observer is called with the final state for the initial guess and for the accepted iterates
        const bool performOptimization = piroParams->sublist("Analysis").sublist("ROL").get<bool>("Perform Optimization", true);
        if (performOptimization && Teuchos::nonnull(p)) {
          if ((observer->numSolutions < 2) || (observer->numResponses != observer->numSolutions) ||
              !std::isfinite(observer->lastSolutionNorm) || (observer->lastSolutionNorm <= 0.0)) {
            status += 10;
            if (Proc==0)
              std::cout << "\nPiro_AnalysisDriver: unexpected observer calls: " << observer->numSolutions << " solutions, "
                        << observer->numResponses << " responses, last solution norm " << observer->lastSolutionNorm << std::endl;
          }
          else if (Proc==0)
            std::cout << "\nPiro_AnalysisDriver: observer called with " << observer->numSolutions << " solutions, last iteration "
                      << observer->lastStamp << ", parameters changed " << observer->numParametersChanged << " times" << std::endl;
        }

        if (Teuchos::nonnull(p)) { //p might be null if the package ROL is not enabled
          Thyra::DetachedVectorView<double> p_view(p);
          double p_exact[2];
          double tol = 1e-5;
          if (mockModel=="MockModelEval_A_Tpetra") {
            p_exact[0] = 1;
            p_exact[1] = 3;
          }
          else if (mockModel=="MockModelEval_B_Tpetra") {
            p_exact[0] = 6;
            p_exact[1] = 4;
          }
          else if (mockModel=="MockModelEval_B_Tpetra_2_parameters") {
            p_exact[0] = 4;
            p_exact[1] = 6;
          }
          else if (mockModel=="MassSpringDamperModel") {
            p_exact[0] = 1.;
            p_exact[1] = 0.5;
          } else {
           TEUCHOS_ASSERT(false);
          }

          // An input file can provide the expected optimum and the tolerance in a "Test" sublist
          if (appParams->isSublist("Test")) {
            Teuchos::ParameterList& testParams = appParams->sublist("Test");
            if (testParams.isParameter("Expected Optimal Parameters")) {
              const auto expected = testParams.get<Teuchos::Array<double>>("Expected Optimal Parameters");
              TEUCHOS_ASSERT(expected.size() == 2);
              p_exact[0] = expected[0];
              p_exact[1] = expected[1];
            }
            tol = testParams.get<double>("Tolerance", tol);
          }

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
