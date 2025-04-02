// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MueLu_KokkosTuningInterface.hpp"

#include <string>
#include <sstream>
#include "Teuchos_Array.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_RawParameterListHelpers.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_Exceptions.hpp"

// ***********************************************************************
/* Notional Parameterlist Structure
   "kokkos tuning: muelu parameter mapping"
     - "input variables"            "{"Chebyshev","parallel_for"}
     - "kokkos context id"          "1"           # Context ID to use, assuming you inject it onto the list
     - "param0"
       - "muelu parameter"          "smoother: params||chebyshev: degree"
       - "discrete range"           "{1,6,1}"    # (low, high, step')
       - "initial guess"            "2"
     - "param1"
       - "muelu parameter"          "smoother: params||chebyshev: eigenvalue ratio"
       - "continuous range"         "{5.0,50.0,5.0}"    # (low, high, step')
       - "initial guess"            "10.0"


  The input variables should be handled by the tuning tool.
  The output variable should be calculated automatically based on the types above.
  Note: We use "||" to indicate sublists.
 */

namespace MueLu {

// FIXME: Will probably need to bump this
namespace KokkosTuningParams {
const int MAX_VALID_PARAMS = 10;
}

// ***********************************************************************
KokkosTuningInterface::KokkosTuningInterface(const Teuchos::RCP<const Teuchos::Comm<int>>& comm)
  : comm_(comm) {
}

// ***********************************************************************
RCP<const Teuchos::ParameterList> KokkosTuningInterface::GetValidParameterList() const {
  RCP<ParameterList> topValidParamList = rcp(new ParameterList());
  ParameterList validParamList;

  ParameterList pl_dummy;
  Teuchos::Array<std::string> ar_dummy;
  std::string s_dummy;

  // Input variables for Kokkos tuning
  validParamList.set<Teuchos::Array<std::string>>("input variables", ar_dummy, "Names of the input variables for Kokkos tuning");

  validParamList.set<size_t>("kokkos context id", Teuchos::OrdinalTraits<size_t>::invalid(), "Context ID for Kokkos tuning");

  for (int i = 0; i < KokkosTuningParams::MAX_VALID_PARAMS; i++) {
    std::ostringstream oss;
    oss << "param" << i;
    const std::string name = oss.str();

    // FIXME: Not validating parameter sublists at present
    validParamList.set<Teuchos::ParameterList>(name, pl_dummy, "Parameter-specific sublist");
  }

  topValidParamList->set<Teuchos::ParameterList>("kokkos tuning: muelu parameter mapping", validParamList, "Sublist for Kokkos tuning of MueLu");

  return topValidParamList;
}

// ***********************************************************************
void KokkosTuningInterface::Setup() {
  // Sanity check
  if (comm_.is_null()) throw std::runtime_error("MueLu::KokkosTuningInterface::Setup(): Communicator cannot be null");

  // Unpack the MueLu Mapping into something actionable
  UnpackMueLuMapping();
}

// ***********************************************************************
void KokkosTuningInterface::UnpackMueLuMapping() {
  const Teuchos::ParameterList& pL = params_.get<Teuchos::ParameterList>("kokkos tuning: muelu parameter mapping");
  namespace KTE                    = Kokkos::Tools::Experimental;

  /********************************/
  /* Process the output variables */
  /********************************/
  out_variables.clear();
  out_names.clear();

  for (int i = 0; i < KokkosTuningParams::MAX_VALID_PARAMS; i++) {
    std::ostringstream oss;
    oss << "param" << i;
    const std::string name = oss.str();

    if (pL.isSublist(name)) {
      const Teuchos::ParameterList& sublist = pL.sublist(name);
      std::string muelu_param               = sublist.get<std::string>("muelu parameter");

      // Infer types from initial guess
      if (sublist.isType<int>("initial guess")) {
        // Discrete range
        int guess                        = sublist.get<int>("initial guess");
        const Teuchos::Array<int>& range = sublist.get<Teuchos::Array<int>>("discrete range");
        TEUCHOS_TEST_FOR_EXCEPTION(range.size() != 3, Exceptions::RuntimeError, "MueLu::KokkosTuningInterface: 'discrete range' needs to be (low, high, step)");

        // Set the VariableInfo
        KTE::VariableInfo out_info;
        out_info.type          = KTE::ValueType::kokkos_value_int64;
        out_info.category      = KTE::StatisticalCategory::kokkos_value_interval;
        out_info.valueQuantity = KTE::CandidateValueType::kokkos_value_range;

        // Unlike the ordinal lists, the ranges get copied into Kokkos
        // TODO: Add support for open/closed ranges
        out_info.candidates = KTE::make_candidate_range((int64_t)range[0], (int64_t)range[1], (int64_t)range[2], false, false);
        size_t var_id       = KTE::declare_output_type(muelu_param, out_info);
        out_variables.push_back(KTE::make_variable_value(var_id, int64_t(guess)));

        out_typenames.push_back("int");
      } else if (sublist.isType<double>("initial guess")) {
        // Continuous range
        double guess                        = sublist.get<double>("initial guess");
        const Teuchos::Array<double>& range = sublist.get<Teuchos::Array<double>>("continuous range");
        TEUCHOS_TEST_FOR_EXCEPTION(range.size() != 3, Exceptions::RuntimeError, "MueLu::KokkosTuningInterface: 'continuous range' needs to be (low, high, step)");

        // Set the VariableInfo
        KTE::VariableInfo out_info;
        out_info.type          = KTE::ValueType::kokkos_value_double;
        out_info.category      = KTE::StatisticalCategory::kokkos_value_interval;
        out_info.valueQuantity = KTE::CandidateValueType::kokkos_value_range;

        // Unlike the ordinal lists, the ranges get copied into Kokkos
        // TODO: Add support for open/closed ranges
        out_info.candidates = KTE::make_candidate_range(range[0], range[1], range[2], false, false);
        size_t var_id       = KTE::declare_output_type(muelu_param, out_info);
        out_variables.push_back(KTE::make_variable_value(var_id, guess));

        out_typenames.push_back("double");
      }
      // TODO: Add support for categorical and set parameters
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::KokkosTuningInterface: We currently only handle int and double ranges.");
      }

      // Stash the parameter name
      out_names.push_back(muelu_param);

    }  // end if pL.isSublist

  }  // end for

  // Sanity check
  TEUCHOS_TEST_FOR_EXCEPTION(out_names.size() != out_variables.size(), Exceptions::RuntimeError, "MueLu::KokkosTuningInterface: Error in option processing");
  TEUCHOS_TEST_FOR_EXCEPTION(out_names.size() != out_typenames.size(), Exceptions::RuntimeError, "MueLu::KokkosTuningInterface: Error in option processing");

  /********************************/
  /* Process the input variables  */
  /********************************/
  in_variables.clear();

  const Teuchos::Array<std::string>& inputs = pL.get<Teuchos::Array<std::string>>("input variables");

  for (int i = 0; i < (int)inputs.size(); i++) {
    // NOTE: The string name is copied in here (unlike double/int) so we don't need to cache.
    KTE::VariableInfo in_info;
    in_info.type  = KTE::ValueType::kokkos_value_string;
    size_t var_id = KTE::declare_input_type(inputs[i].c_str(), in_info);
    in_variables.push_back(KTE::make_variable_value(var_id, inputs[i].c_str()));
  }
}

// ***********************************************************************
std::vector<std::string> KokkosTuningInterface::SplitString(const std::string& base_string, const std::string& delimiter) const {
  std::vector<std::string> tokens;
  size_t start = 0;

  size_t end = base_string.find(delimiter);

  while (end != std::string::npos) {
    tokens.push_back(base_string.substr(start, end - start));
    start = end + delimiter.length();
    end   = base_string.find(delimiter, start);
  }

  // And the final token...
  tokens.push_back(base_string.substr(start, end));

  return tokens;
}

// ***********************************************************************
void KokkosTuningInterface::SetMueLuParameters(size_t kokkos_context_id, Teuchos::ParameterList& mueluParams, bool overwrite) const {
  namespace KTE = Kokkos::Tools::Experimental;
  Teuchos::ParameterList tunedParams;

  if (comm_->getRank() == 0) {
    // Only Rank 0 calls KokkosTuning
    GetOStream(Runtime0) << "MueLu::KokkosTuningInterface: Tuning " << out_variables.size() << " parameters" << std::endl;

    // Set input variable
    if (IsPrint(Runtime1))
      GetOStream(Runtime1) << "Adding " << in_variables.size() << " input variables" << std::endl;

    KTE::set_input_values(kokkos_context_id, in_variables.size(), in_variables.data());

    // Start the tuning
    KTE::request_output_values(kokkos_context_id, out_variables.size(), out_variables.data());

    // Diagnostic output
    if (IsPrint(Runtime1)) {
      GetOStream(Runtime1) << "Tuned Parameters: " << std::endl;
      for (int i = 0; i < (int)out_variables.size(); i++) {
        if (out_typenames[i] == "int")
          GetOStream(Runtime1) << "- " << out_names[i] << ": " << out_variables[i].value.int_value << std::endl;
        else if (out_typenames[i] == "double")
          GetOStream(Runtime1) << "- " << out_names[i] << ": " << out_variables[i].value.double_value << std::endl;
      }
    }

    // Unpack the tuned values
    for (int i = 0; i < (int)out_names.size(); i++) {
      // Because we'll want to set parameters inside sublists we'll allow the "muelu parameter" option to specify sublists with '||' as a nesting delimiter
      // That's really, really unlikely to be part of a parameter list item name, so we'll go with it.
      Teuchos::ParameterList* activeList = &tunedParams;
      std::vector<std::string> treeWalk  = SplitString(out_names[i], "||");

      // Walk down all but the last guy
      for (int j = 0; j < (int)treeWalk.size() - 1; j++) {
        activeList = &(activeList->sublist(treeWalk[j]));
      }

      std::string activeName = treeWalk[treeWalk.size() - 1];

      if (out_typenames[i] == "int")
        activeList->set(activeName, (int)out_variables[i].value.int_value);
      else if (out_typenames[i] == "double")
        activeList->set(activeName, out_variables[i].value.double_value);
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::KokkosTuningInterface: Unknown variable output type");
      }
    }
  }

  Teuchos::updateParametersAndBroadcast(outArg(tunedParams), outArg(mueluParams), *comm_, 0, overwrite);
}

// ***********************************************************************
void KokkosTuningInterface::SetMueLuParameters(Teuchos::ParameterList& mueluParams, bool overwrite) const {
  size_t PL_kokkos_context_id = mueluParams.get<Teuchos::ParameterList>("kokkos tuning: muelu parameter mapping").get<size_t>("kokkos context id");

  SetMueLuParameters(PL_kokkos_context_id, mueluParams, overwrite);
}

}  // namespace MueLu
