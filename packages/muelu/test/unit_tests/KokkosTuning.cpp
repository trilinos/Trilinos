// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>

#include <MueLu_config.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <MueLu_KokkosTuningInterface.hpp>


namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(KokkosTuningInterface, Basic, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();



  // Driver parameters (like you'd actually use)
  Teuchos::ParameterList baseList;
  Teuchos::ParameterList & pL = baseList.sublist("kokkos tuning: muelu parameter mapping");
  Teuchos::Array<std::string> input_vars(1); input_vars[0] = "Parameter List Item";
  Teuchos::Array<int> i_range{1,6,1};
  Teuchos::Array<double> d_range{5.0,50.0,5.0};
  pL.set("input variables",input_vars);
  pL.sublist("param0").set("muelu parameter","smoother: params||chebyshev: degree");
  pL.sublist("param0").set("discrete range",i_range);
  pL.sublist("param0").set("initial guess",(int) 1);
  pL.sublist("param1").set("muelu parameter","smoother: params||chebyshev: ratio eigenvalue");
  pL.sublist("param1").set("continuous range",d_range);
  pL.sublist("param1").set("initial guess",(double) 10.0);


  // Actually make the interface
  MueLu::KokkosTuningInterface interface(comm);
  interface.SetParameterList(baseList);

  // Call the tuner
  size_t kokkos_context_id;
  Kokkos::Tools::Experimental::begin_context(kokkos_context_id);
  Teuchos::ParameterList outputList;
  interface.SetMueLuParameters(kokkos_context_id,outputList);
  Kokkos::Tools::Experimental::end_context(kokkos_context_id);

  // Check that the output has the varables set to something
  TEST_EQUALITY(outputList.isSublist("smoother: params"),true);
  TEST_EQUALITY(outputList.sublist("smoother: params").isParameter("chebyshev: degree"),true);
  TEST_EQUALITY(outputList.sublist("smoother: params").isParameter("chebyshev: ratio eigenvalue"),true);


}


TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(KokkosTuningInterface, Advanced, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

 RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

 RCP<MueLu::KokkosTuningInterface> interface = rcp(new MueLu::KokkosTuningInterface(comm));

  TEST_INEQUALITY(interface, Teuchos::null);
}


#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(KokkosTuningInterface, Basic, Scalar, LO, GO, Node)   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(KokkosTuningInterface, Advanced, Scalar, LO, GO, Node)


#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
