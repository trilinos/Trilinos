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
#include <Kokkos_Core.hpp>

#include <MueLu_config.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <MueLu_KokkosTuningInterface.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(KokkosTuningInterface, Basic, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  // This test just makes sure the KokkosTuning interface does *something* plausible
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  // Driver parameters (like you'd actually use)
  Teuchos::ParameterList baseList;
  Teuchos::ParameterList& pL = baseList.sublist("kokkos tuning: muelu parameter mapping");
  Teuchos::Array<std::string> input_vars(1);
  input_vars[0] = "Parameter List Item";
  Teuchos::Array<int> i_range{1, 6, 1};
  Teuchos::Array<double> d_range{5.0, 50.0, 5.0};
  pL.set("input variables", input_vars);
  pL.sublist("param0").set("muelu parameter", "smoother: params||chebyshev: degree");
  pL.sublist("param0").set("discrete range", i_range);
  pL.sublist("param0").set("initial guess", (int)1);
  pL.sublist("param1").set("muelu parameter", "smoother: params||chebyshev: ratio eigenvalue");
  pL.sublist("param1").set("continuous range", d_range);
  pL.sublist("param1").set("initial guess", (double)10.0);

  // Actually make the interface
  MueLu::KokkosTuningInterface interface(comm);
  interface.SetParameterList(baseList);

  // Call the tuner
  size_t kokkos_context_id = 0;
  Kokkos::Tools::Experimental::begin_context(kokkos_context_id);
  Teuchos::ParameterList outputList;
  interface.SetMueLuParameters(kokkos_context_id, outputList);
  Kokkos::Tools::Experimental::end_context(kokkos_context_id);

  // Check that the output has the varables set to something
  TEST_EQUALITY(outputList.isSublist("smoother: params"), true);
  TEST_EQUALITY(outputList.sublist("smoother: params").isParameter("chebyshev: degree"), true);
  TEST_EQUALITY(outputList.sublist("smoother: params").isParameter("chebyshev: ratio eigenvalue"), true);
}

namespace TestUtilities {
std::vector<Kokkos::Tools::Experimental::VariableInfo> output_info;
std::vector<Kokkos::Tools::Experimental::VariableInfo> input_info;

const char* str_type(Kokkos_Tools_VariableInfo_ValueType x) {
  if (x == Kokkos_Tools_VariableInfo_ValueType::kokkos_value_double)
    return "double";
  else if (x == Kokkos_Tools_VariableInfo_ValueType::kokkos_value_int64)
    return "int64";
  else if (x == Kokkos_Tools_VariableInfo_ValueType::kokkos_value_string)
    return "string";
  else
    return "unknown";
}

const char* str_category(const enum Kokkos_Tools_VariableInfo_StatisticalCategory x) {
  if (x == Kokkos_Tools_VariableInfo_StatisticalCategory::kokkos_value_categorical)
    return "categorical";
  else if (x == Kokkos_Tools_VariableInfo_StatisticalCategory::kokkos_value_ordinal)
    return "ordinal";
  else if (x == Kokkos_Tools_VariableInfo_StatisticalCategory::kokkos_value_interval)
    return "interval";
  else
    return "unknown";
}

const char* str_candidatevalue(const enum Kokkos_Tools_VariableInfo_CandidateValueType x) {
  if (x == Kokkos_Tools_VariableInfo_CandidateValueType::kokkos_value_set)
    return "set";
  else if (x == Kokkos_Tools_VariableInfo_CandidateValueType::kokkos_value_range)
    return "range";
  else if (x == Kokkos_Tools_VariableInfo_CandidateValueType::kokkos_value_unbounded)
    return "unbounded";
  else
    return "unknown";
}

void print_variable_info(const Kokkos::Tools::Experimental::VariableInfo* info, const int index = 0) {
  std::cout << "[" << index << "] Variable information  = " << str_type(info->type) << "," << str_category(info->category) << "," << str_candidatevalue(info->valueQuantity) << std::endl;
  std::cout << "[" << index << "] Actual values         = " << info->type << "," << info->category << "," << info->valueQuantity << std::endl;
  std::cout << "[" << index << "] Baseline enums        = " << Kokkos_Tools_VariableInfo_ValueType::kokkos_value_double << "," << Kokkos_Tools_VariableInfo_StatisticalCategory::kokkos_value_categorical << "," << Kokkos_Tools_VariableInfo_CandidateValueType::kokkos_value_set << std::endl;
}

void declare_input_type(const char* name, const size_t id,
                        Kokkos::Tools::Experimental::VariableInfo* info) {
  // We copy this data in and assume the default constructor works
  // std::cout<<"DEBUG: calling declare_input_type"<<std::endl;
  // print_variable_info(info,-1);
  input_info.push_back(*info);
}

void declare_output_type(const char* name, const size_t id,
                         Kokkos::Tools::Experimental::VariableInfo* info) {
  // We copy this data in and assume the default constructor works
  // std::cout<<"DEBUG: calling declare_output_type"<<std::endl;
  // print_variable_info(info,-1);
  output_info.push_back(*info);
}

void request_output_values(const size_t context, const size_t num_inputs,
                           const Kokkos::Tools::Experimental::VariableValue* inputs_in,
                           const size_t num_outputs,
                           Kokkos::Tools::Experimental::VariableValue* outputs_in) {
  // This dummy callback will set the output value to one step more than the bottom guy in the range
  // std::cout << "\nDEBUG: request_ouput_values called with " << num_outputs << " outputs" << std::endl;
  for (int i = 0; i < (int)num_outputs; i++) {
    Kokkos::Tools::Experimental::VariableInfo& info = output_info[i];
    // print_variable_info(&info,i);
    if (info.category == Kokkos_Tools_VariableInfo_StatisticalCategory::kokkos_value_interval &&
        info.valueQuantity == Kokkos_Tools_VariableInfo_CandidateValueType::kokkos_value_range) {
      auto range = info.candidates.range;

      // We only suport ranges in this test
      if (info.type == Kokkos_Tools_VariableInfo_ValueType::kokkos_value_int64) {
        outputs_in[i].value.int_value = (int)(range.lower.int_value + range.step.int_value);
        // std::cout << "Setting parameter " << i << " to value = " << outputs_in[i].value.int_value << std::endl;
      } else if (info.type == Kokkos_Tools_VariableInfo_ValueType::kokkos_value_double) {
        outputs_in[i].value.double_value = range.lower.double_value + range.step.double_value;
        // std::cout << "Setting parameter " << i << " to value = " << outputs_in[i].value.double_value << std::endl;
      }
    }
  }
}

}  // namespace TestUtilities

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(KokkosTuningInterface, Advanced, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  // This test makes sure the KokkosTuning interface does something specific
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  // Driver parameters (like you'd actually use)
  Teuchos::ParameterList baseList;
  Teuchos::ParameterList& pL = baseList.sublist("kokkos tuning: muelu parameter mapping");
  Teuchos::Array<std::string> input_vars(1);
  input_vars[0] = "Parameter List Item";
  Teuchos::Array<int> i_range{1, 6, 1};
  Teuchos::Array<double> d_range{5.0, 50.0, 5.0};
  pL.set("input variables", input_vars);
  pL.sublist("param0").set("muelu parameter", "smoother: params||chebyshev: degree");
  pL.sublist("param0").set("discrete range", i_range);
  pL.sublist("param0").set("initial guess", (int)1);
  pL.sublist("param1").set("muelu parameter", "smoother: params||chebyshev: ratio eigenvalue");
  pL.sublist("param1").set("continuous range", d_range);
  pL.sublist("param1").set("initial guess", (double)5.0);

  // Set the callbacks
  Kokkos::Tools::Experimental::set_declare_input_type_callback(TestUtilities::declare_input_type);
  Kokkos::Tools::Experimental::set_declare_output_type_callback(TestUtilities::declare_output_type);
  Kokkos::Tools::Experimental::set_request_output_values_callback(TestUtilities::request_output_values);

  // Actually make the interface
  MueLu::KokkosTuningInterface interface(comm);
  interface.SetParameterList(baseList);

  // Call the tuner
  size_t kokkos_context_id = 0;
  Kokkos::Tools::Experimental::begin_context(kokkos_context_id);
  Teuchos::ParameterList outputList;
  interface.SetMueLuParameters(kokkos_context_id, outputList);
  Kokkos::Tools::Experimental::end_context(kokkos_context_id);

  // Check that the output has the varables set to something
  TEST_EQUALITY(outputList.isSublist("smoother: params"), true);
  TEST_EQUALITY(outputList.sublist("smoother: params").isParameter("chebyshev: degree"), true);
  TEST_EQUALITY(outputList.sublist("smoother: params").isParameter("chebyshev: ratio eigenvalue"), true);

#ifdef KOKKOS_ENABLE_TUNING
  // Check that the variables are set to one step off the bottom
  int degree   = 2;
  double ratio = 10.0;
#else
  // Expect the default values
  int degree   = 1;
  double ratio = 5.0;
#endif
  TEST_EQUALITY(outputList.sublist("smoother: params").get<int>("chebyshev: degree"), degree);
  TEST_FLOATING_EQUALITY(outputList.sublist("smoother: params").get<double>("chebyshev: ratio eigenvalue"), ratio, 1e-10);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(KokkosTuningInterface, Advanced_PLContext, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  // This test makes sure the KokkosTuning interface does something specific
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

  // Driver parameters (like you'd actually use)
  Teuchos::ParameterList baseList;
  Teuchos::ParameterList& pL = baseList.sublist("kokkos tuning: muelu parameter mapping");
  Teuchos::Array<std::string> input_vars(1);
  input_vars[0] = "Parameter List Item";
  Teuchos::Array<int> i_range{1, 6, 1};
  Teuchos::Array<double> d_range{5.0, 50.0, 5.0};
  pL.set("input variables", input_vars);
  pL.sublist("param0").set("muelu parameter", "smoother: params||chebyshev: degree");
  pL.sublist("param0").set("discrete range", i_range);
  pL.sublist("param0").set("initial guess", (int)1);
  pL.sublist("param1").set("muelu parameter", "smoother: params||chebyshev: ratio eigenvalue");
  pL.sublist("param1").set("continuous range", d_range);
  pL.sublist("param1").set("initial guess", (double)5.0);

  // Set the callbacks
  Kokkos::Tools::Experimental::set_declare_input_type_callback(TestUtilities::declare_input_type);
  Kokkos::Tools::Experimental::set_declare_output_type_callback(TestUtilities::declare_output_type);
  Kokkos::Tools::Experimental::set_request_output_values_callback(TestUtilities::request_output_values);

  // Actually make the interface
  MueLu::KokkosTuningInterface interface(comm);

  // Call the tuner, using a PL-based context ID
  size_t kokkos_context_id = 0;
  Kokkos::Tools::Experimental::begin_context(kokkos_context_id);
  interface.SetParameterList(baseList);

  Teuchos::ParameterList outputList;
  outputList.sublist("kokkos tuning: muelu parameter mapping").set("kokkos context id", kokkos_context_id);
  interface.SetMueLuParameters(outputList);

  Kokkos::Tools::Experimental::end_context(kokkos_context_id);

  // Check that the output has the varables set to something
  TEST_EQUALITY(outputList.isSublist("smoother: params"), true);
  TEST_EQUALITY(outputList.sublist("smoother: params").isParameter("chebyshev: degree"), true);
  TEST_EQUALITY(outputList.sublist("smoother: params").isParameter("chebyshev: ratio eigenvalue"), true);

#ifdef KOKKOS_ENABLE_TUNING
  // Check that the variables are set to one step off the bottom
  int degree   = 2;
  double ratio = 10.0;
#else
  // Expect the default values
  int degree   = 1;
  double ratio = 5.0;
#endif
  TEST_EQUALITY(outputList.sublist("smoother: params").get<int>("chebyshev: degree"), degree);
  TEST_FLOATING_EQUALITY(outputList.sublist("smoother: params").get<double>("chebyshev: ratio eigenvalue"), ratio, 1e-10);
}

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                                 \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(KokkosTuningInterface, Basic, Scalar, LO, GO, Node)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(KokkosTuningInterface, Advanced, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(KokkosTuningInterface, Advanced_PLContext, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
