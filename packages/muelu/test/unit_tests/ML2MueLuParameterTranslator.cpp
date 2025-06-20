// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <MueLu_config.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <MueLu_ML2MueLuParameterTranslator.hpp>

#ifdef HAVE_MUELU_ML
#include <ml_MultiLevelPreconditioner.h>
#include <ml_RefMaxwell.h>
#endif

namespace MueLuTests {

bool compareLists(Teuchos::ParameterList& L1, Teuchos::ParameterList& L2) {
  return Teuchos::haveSameValuesSorted(L1, L2, true);
}

TEUCHOS_UNIT_TEST(ML2MueLuParameterTranslator, SA) {
  // SetDefaults(SA)
  Teuchos::ParameterList List, goldList;
  MueLu::ML2MueLuParameterTranslator::SetDefaults("SA", List);

#ifdef HAVE_MUELU_ML
  Teuchos::ParameterList mlList;
  ML_Epetra::SetDefaults("SA", mlList);
  mlList.set("parameterlist: syntax", "ml");
  TEST_EQUALITY(compareLists(List, mlList), true);
#endif

  // Gold list
  goldList.set("parameterlist: syntax", "ml");
  goldList.set("default values", "SA");
  goldList.set("max levels", 10);
  goldList.set("prec type", "MGV");
  goldList.set("increasing or decreasing", "increasing");
  goldList.set("aggregation: type", "Uncoupled-MIS");
  goldList.set("aggregation: damping factor", 1.333);
  goldList.set("eigen-analysis: type", "cg");
  goldList.set("eigen-analysis: iterations", 10);
  goldList.set("smoother: sweeps", 2);
  goldList.set("smoother: damping factor", 1.0);
  goldList.set("smoother: pre or post", "both");
  goldList.set("smoother: type", "symmetric Gauss-Seidel");
  goldList.set("coarse: type", "Amesos-KLU");
  goldList.set("coarse: max size", 128);
  goldList.set("coarse: pre or post", "post");
  goldList.set("coarse: sweeps", 1);
  goldList.set("coarse: split communicator", false);

  TEST_EQUALITY(compareLists(List, goldList), true);
}  // SA

TEUCHOS_UNIT_TEST(ML2MueLuParameterTranslator, SA_plus_translate) {
  // SetDefaults(SA)
  Teuchos::ParameterList List, goldList, dummy;
  MueLu::ML2MueLuParameterTranslator::SetDefaults("SA", List);

  std::string syntaxStr = "parameterlist: syntax";
  TEST_EQUALITY(List.isParameter(syntaxStr), true);
  std::string listStr = List.get<std::string>(syntaxStr);
  TEST_EQUALITY(listStr, "ml");

  List.remove(syntaxStr);
  std::string paramXML = MueLu::ML2MueLuParameterTranslator::translate(List, "");
  List                 = *Teuchos::getParametersFromXmlString(paramXML);

  std::cout << "\n-----------------------\n"
            << List << std::endl;

  // gold list
  goldList.set("aggregation: match ML phase1", true);
  goldList.set("aggregation: match ML phase2a", true);
  goldList.set("aggregation: match ML phase2b", true);
  goldList.set("aggregation: use ml scaling of drop tol", true);
  goldList.set("max levels", 10);
  goldList.set("cycle type", "V");
  goldList.set("sa: damping factor", 1.333);
  goldList.set("sa: eigenvalue estimate num iterations", 10);
  goldList.set("smoother: type", "RELAXATION");
  goldList.set("smoother: pre or post", "both");
  goldList.set("repartition: start level", 2);
  goldList.set("repartition: put on single proc", 5000);
  goldList.set("coarse: type", "klu");
  goldList.set("coarse: max size", 128);
  goldList.sublist("smoother: params").set("relaxation: type", "Symmetric Gauss-Seidel");
  goldList.sublist("smoother: params").set("relaxation: sweeps", 2);
  goldList.sublist("smoother: params").set("relaxation: damping factor", 1.0);
  goldList.set("coarse: params", dummy);

  std::cout << "\n-----------------------\n"
            << goldList << std::endl;

  TEST_EQUALITY(compareLists(List, goldList), true);

}  // SA_plus_translate

TEUCHOS_UNIT_TEST(ML2MueLuParameterTranslator, RefMaxwell) {
  // SetDefaults(SA)
  Teuchos::ParameterList List, goldList;
  MueLu::ML2MueLuParameterTranslator::SetDefaultsRefMaxwell(List);

#ifdef HAVE_MUELU_ML
  Teuchos::ParameterList mlList;
  ML_Epetra::SetDefaultsRefMaxwell(mlList);
  mlList.set("parameterlist: syntax", "ml");
  TEST_EQUALITY(compareLists(List, mlList), true);
#endif

  // Gold list
  goldList.set("parameterlist: syntax", "ml");
  goldList.set("default values", "RefMaxwell");
  goldList.set("max levels", 10);
  goldList.set("prec type", "MGV");
  goldList.set("increasing or decreasing", "decreasing");
  goldList.set("aggregation: type", "Uncoupled-MIS");
  goldList.set("aggregation: damping factor", 1.333);
  goldList.set("eigen-analysis: type", "cg");
  goldList.set("eigen-analysis: iterations", 10);
  goldList.set("aggregation: edge prolongator drop threshold", 0.0);
  goldList.set("smoother: sweeps", 2);
  goldList.set("smoother: damping factor", 1.0);
  goldList.set("smoother: pre or post", "both");
  goldList.set("smoother: type", "Chebyshev");
  goldList.set("smoother: Hiptmair efficient symmetric", true);
  goldList.set("subsmoother: type", "Chebyshev");
  goldList.set("subsmoother: Chebyshev alpha", 20.0);
  goldList.set("subsmoother: node sweeps", 4);
  goldList.set("subsmoother: edge sweeps", 4);
  goldList.set("coarse: type", "Amesos-KLU");
  goldList.set("coarse: max size", 128);
  goldList.set("coarse: pre or post", "post");
  goldList.set("coarse: sweeps", 1);
  goldList.set("refmaxwell: 11solver", "edge matrix free");
  goldList.set("refmaxwell: 22solver", "multilevel");
  goldList.set("refmaxwell: mode", "additive");
  goldList.set("zero starting solution", false);
  goldList.sublist("refmaxwell: 11list").set("default values", "SA");
  goldList.sublist("refmaxwell: 11list").set("max levels", 10);
  goldList.sublist("refmaxwell: 11list").set("prec type", "MGV");
  goldList.sublist("refmaxwell: 11list").set("increasing or decreasing", "increasing");
  goldList.sublist("refmaxwell: 11list").set("aggregation: type", "Uncoupled");
  goldList.sublist("refmaxwell: 11list").set("aggregation: damping factor", 0.0);
  goldList.sublist("refmaxwell: 11list").set("eigen-analysis: type", "cg");
  goldList.sublist("refmaxwell: 11list").set("eigen-analysis: iterations", 10);
  goldList.sublist("refmaxwell: 11list").set("smoother: sweeps", 0);
  goldList.sublist("refmaxwell: 11list").set("smoother: damping factor", 1.0);
  goldList.sublist("refmaxwell: 11list").set("smoother: pre or post", "both");
  goldList.sublist("refmaxwell: 11list").set("smoother: type", "symmetric Gauss-Seidel");
  goldList.sublist("refmaxwell: 11list").set("coarse: type", "Amesos-KLU");
  goldList.sublist("refmaxwell: 11list").set("coarse: max size", 128);
  goldList.sublist("refmaxwell: 11list").set("coarse: pre or post", "post");
  goldList.sublist("refmaxwell: 11list").set("coarse: sweeps", 1);
  goldList.sublist("refmaxwell: 11list").set("coarse: split communicator", false);
  goldList.sublist("refmaxwell: 11list").set("cycle applications", 1);
  goldList.sublist("refmaxwell: 11list").set("aggregation: threshold", 0.01);
  goldList.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse").set("default values", "SA");
  goldList.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse").set("max levels", 10);
  goldList.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse").set("prec type", "MGV");
  goldList.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse").set("increasing or decreasing", "increasing");
  goldList.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse").set("aggregation: type", "Uncoupled-MIS");
  goldList.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse").set("aggregation: damping factor", 1.333);
  goldList.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse").set("eigen-analysis: type", "cg");
  goldList.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse").set("eigen-analysis: iterations", 10);
  goldList.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse").set("smoother: sweeps", 2);
  goldList.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse").set("smoother: damping factor", 1.0);
  goldList.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse").set("smoother: pre or post", "both");
  goldList.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse").set("smoother: type", "Chebyshev");
  goldList.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse").set("coarse: type", "Amesos-KLU");
  goldList.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse").set("coarse: max size", 128);
  goldList.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse").set("coarse: pre or post", "post");
  goldList.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse").set("coarse: sweeps", 1);
  goldList.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse").set("coarse: split communicator", false);
  goldList.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse").set("cycle applications", 1);
  goldList.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse").set("aggregation: threshold", 0.01);
  goldList.sublist("refmaxwell: 11list").sublist("edge matrix free: coarse").set("ML label", "coarse (1,1) block");
  goldList.sublist("refmaxwell: 22list").set("default values", "SA");
  goldList.sublist("refmaxwell: 22list").set("max levels", 10);
  goldList.sublist("refmaxwell: 22list").set("prec type", "MGV");
  goldList.sublist("refmaxwell: 22list").set("increasing or decreasing", "increasing");
  goldList.sublist("refmaxwell: 22list").set("aggregation: type", "Uncoupled");
  goldList.sublist("refmaxwell: 22list").set("aggregation: damping factor", 1.333);
  goldList.sublist("refmaxwell: 22list").set("eigen-analysis: type", "cg");
  goldList.sublist("refmaxwell: 22list").set("eigen-analysis: iterations", 10);
  goldList.sublist("refmaxwell: 22list").set("smoother: sweeps", 2);
  goldList.sublist("refmaxwell: 22list").set("smoother: damping factor", 1.0);
  goldList.sublist("refmaxwell: 22list").set("smoother: pre or post", "both");
  goldList.sublist("refmaxwell: 22list").set("smoother: type", "Chebyshev");
  goldList.sublist("refmaxwell: 22list").set("coarse: type", "Amesos-KLU");
  goldList.sublist("refmaxwell: 22list").set("coarse: max size", 128);
  goldList.sublist("refmaxwell: 22list").set("coarse: pre or post", "post");
  goldList.sublist("refmaxwell: 22list").set("coarse: sweeps", 1);
  goldList.sublist("refmaxwell: 22list").set("coarse: split communicator", false);
  goldList.sublist("refmaxwell: 22list").set("cycle applications", 1);
  goldList.sublist("refmaxwell: 22list").set("aggregation: threshold", 0.01);
  goldList.sublist("refmaxwell: 22list").set("ML label", "(2,2) block");

  TEST_EQUALITY(compareLists(List, goldList), true);
}

}  // namespace MueLuTests
