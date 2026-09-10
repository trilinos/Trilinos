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

#include <iostream>
#include <string>

#include "MockModelEval_H_Tpetra.hpp"
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
#include "Piro_ProductModelEval.hpp"
#include "Piro_HDSA_MD_ROL_Data_Interface.hpp"

#ifdef HAVE_PIRO_MUELU
#include "Stratimikos_MueLuHelpers.hpp"
#endif




#include "Piro_ConfigDefs.hpp"

int main(int argc, char *argv[]) {

  int status=0; // 0 = pass, failures are incremented
  int overall_status=0; // 0 = pass, failures are incremented over multiple tests
  bool success=true;

  // Initialize MPI
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  int Proc=mpiSession.getRank();

  HDSA::ROL_Vector<double>::nr = ROL::Elementwise::NormalRandom<double>(0.0,1.0,42+Proc);


  auto appComm = Tpetra::getDefaultComm();

  using Teuchos::RCP;
  using Teuchos::rcp;

  std::string inputFile;
  bool doAll = (argc==1);
  if (argc>1) doAll = !strcmp(argv[1],"-v");

  Piro::SolverFactory solverFactory;

  for (int iTest=0; iTest<1; iTest++) {

    if (doAll) {
      switch (iTest) {
       case 0: inputFile="input_Analysis_HDSA.xml"; break;
       default : std::cout << "iTest logic error " << std::endl; exit(-1);
      }
    }
    else {
      inputFile=argv[1];
      iTest = 999;
    }

    try {

      std::vector<std::string> mockModels = {"MockModelEval_H_Tpetra"};
      RCP<Thyra::VectorBase<double>> p_opt, true_p_opt;
      std::vector<RCP<Thyra::VectorBase<double>>> p_samples, u_diff_at_samples;
      for (auto mockModel : mockModels) {
        p_samples.clear();
        u_diff_at_samples.clear();

        // BEGIN Builder
        const RCP<Teuchos::ParameterList> appParams = rcp(new Teuchos::ParameterList("Application Parameters"));
        Teuchos::updateParametersFromXmlFile(inputFile, Teuchos::ptr(appParams.get()));

        const RCP<Teuchos::ParameterList>  probParams = Teuchos::sublist(appParams,"Problem");
        const RCP<Teuchos::ParameterList>  piroParams = Teuchos::sublist(appParams,"Piro");
     
        // Create (1) a Model Evaluator and (2) a ParameterList
        std::string modelName;
        bool adjoint = (piroParams->get("Sensitivity Method", "Forward") == "Adjoint");
        bool explicitAdjointME = adjoint && piroParams->get("Explicit Adjoint Model Evaluator", false);
        RCP<Thyra::ModelEvaluator<double>> model, adjointModel(Teuchos::null);

        int num_parameters = piroParams->sublist("Analysis").sublist("HDSA").get<int>("Number Of Parameters", 1);
        std::vector<int> p_indices(num_parameters);

        for(int i=0; i<num_parameters; ++i) {
          std::ostringstream ss; ss << "Parameter Vector Index " << i;
          p_indices[i] = piroParams->sublist("Analysis").sublist("HDSA").get<int>(ss.str(), i);
        }

        if (mockModel=="MockModelEval_H_Tpetra") {
          auto model_tmp = rcp(new MockModelEval_H_Tpetra(appComm,false,probParams,true));
          p_opt = model_tmp->get_p_opt(0);
          true_p_opt = model_tmp->get_true_p_opt(0);
          for (int k=0; k<2; k++) {
            p_samples.push_back(model_tmp->get_param_samples(k));
            u_diff_at_samples.push_back(model_tmp->get_solution_diff_at_samples(k));
          }
          model = rcp(new Piro::ProductModelEvaluator<double>(model_tmp,p_indices));
          if(explicitAdjointME) {
            RCP<Thyra::ModelEvaluator<double>> adjointModel_tmp = rcp(new MockModelEval_H_Tpetra(appComm,true));
            adjointModel = rcp(new Piro::ProductModelEvaluator<double>(adjointModel_tmp,p_indices));
          }
          modelName = "H";
        }
        if (Proc==0)
          std::cout << "=======================================================================================================\n"
                    << "======  Solving Problem " << modelName << " with input file "<< iTest <<": "<< inputFile <<"\n"
                    << "=======================================================================================================\n"
            << std::endl;
        

        Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

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

        const RCP<Thyra::ModelEvaluatorDefaultBase<double>> piro = solverFactory.createSolver(piroParams, modelWithSolve, adjointModelWithSolve);

        // Call the analysis routine
        RCP<Thyra::VectorBase<double>> p;  //the parameter nominal value in modelWithSolve is supposed to be the same as p_opt (low fidelity optimal paramer) 
        status = Piro::PerformAnalysis(*piro, *piroParams, p, Teuchos::null, u_diff_at_samples, p_samples);
        Teuchos::RCP<const Thyra::ProductVectorBase<double> > p_prodvec = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double>>(p);

        auto p_opt_vec = ConverterT::getConstTpetraVector(p_opt);
        auto p_vec = ConverterT::getConstTpetraVector(p_prodvec->getVectorBlock(0));
        auto true_p_vec = ConverterT::getConstTpetraVector(true_p_opt);
        
        Tpetra::MatrixMarket::Writer< Tpetra::MultiVector<> >   vecWriter;

        vecWriter.writeDenseFile("z_initial", p_opt_vec);
        vecWriter.writeDenseFile("z_updated", p_vec);
        vecWriter.writeDenseFile("z_true", true_p_vec);

        Tpetra_Vector diff(p_vec->getMap());

        // diff = p_opt_vec - true_p_vec
        diff.update(1.0, *p_opt_vec, -1.0, *true_p_vec, 0.0 );
        double initRelErr = diff.norm2()/true_p_vec->norm2();

        // diff = p_vec - true_p_vec
        diff.update(1.0, *p_vec, -1.0, *true_p_vec, 0.0 );
        double relErr = diff.norm2()/true_p_vec->norm2();

        double expectedRelErr = 0.00210408306503;
        
        double tol = 1e-6;  
        if(std::abs(relErr-expectedRelErr) > tol) {
          status+=100;
          if (Proc==0) {
            std::cout << "\nPiro_AnalysisDrvier_HDSA:  Expected relative error: "
                  <<  expectedRelErr << ", but computed relative error: " 
                  <<  relErr <<
                  "\n Absolute difference grater than tol: " << tol <<   std::endl;
          }
        }         

        if (Proc==0)
          std::cout << std::setprecision(12) << "Relative error of (initial) optimal parameter from low-fidelity model: " << initRelErr << 
          "\nRelative error of parameter updated with HSDA: " << relErr << std::endl;

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
