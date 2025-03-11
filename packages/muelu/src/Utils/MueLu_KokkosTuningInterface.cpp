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
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_RawParameterListHelpers.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_Exceptions.hpp"



// ***********************************************************************
/* Notional Parameterlist Structure
   "kokkos tuning: muelu parameter mapping"
     - "input variables"            "{"Chebyshev","parallel_for"}
     - "param0"
       - "muelu parameter"          "smoother: sweeps"
       - "discrete range"           "{1,6,1}"    # (low, high, step')
       - "initial guess"            "2"
     - "param1"
       - "muelu parameter"          "chebyshev: eigenvalue ratio"
       - "continuous range"         "{5.0,50.0,5.0}"    # (low, high, step')
       - "initial guess"            "10.0"


  The input variables should be handled by the tuning tool.
  The output variable should be calculated automatically based on the types above
 */


namespace MueLu {

  // FIXME: Will probably need to bump this
  namespace KokkosTuningParams {
    const int MAX_VALID_PARAMS=10;
  };

// ***********************************************************************
RCP<const Teuchos::ParameterList> KokkosTuningInterface::GetValidParameterList() const {
  RCP<ParameterList> topValidParamList = rcp(new ParameterList());
  ParameterList validParamList;

  ParameterList pl_dummy;
  Teuchos::Array<std::string> ar_dummy;
  std::string s_dummy;


  // Input variables for Kokkos tuning
  validParamList.set<Teuchos::Array<std::string>>("input variables", ar_dummy, "Names of the input variables for Kokkos tuning");

  for(int i=0;i<KokkosTuningParams::MAX_VALID_PARAMS; i++) {
    std:: ostringstream oss;
    oss<< "param" <<i;
    const std::string name = oss.str();

    // FIXME: Not validating parameter sublists at present
    validParamList.set<Teuchos::ParameterList>(name,pl_dummy,"Parameter-specific sublist");
  }

  topValidParamList->set<Teuchos::ParameterList>("kokkos tuning: muelu parameter mapping",validParamList,"Sublist for Kokkos tuning of MueLu");

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
  namespace KTE = Kokkos::Tools::Experimental;


  /********************************/
  /* Process the output variables */
  /********************************/
  out_variables.clear();
  out_names.clear();

  for(int i=0;i<KokkosTuningParams::MAX_VALID_PARAMS; i++) {
    std:: ostringstream oss;
    oss<< "param" <<i;
    const std::string name = oss.str();

    if (pL.isSublist(name)) {
      const Teuchos::ParameterList& sublist = pL.sublist(name);
      std::string muelu_param = sublist.get<std::string>("muelu parameter");

      // Infer types from initial guess
      if(sublist.isType<int>("initial guess")) {
        // Discrete range
        int guess = sublist.get<int>("initial guess");
        const Teuchos::Array<int> & range = sublist.get<Teuchos::Array<int> >("discrete range");
        TEUCHOS_TEST_FOR_EXCEPTION(range.size() !=3, Exceptions::RuntimeError, "MueLu::KokkosTuningInterface: 'discrete range' needs to be (low, high, step)");

        // Copy the range to int64_t because Kokkos
        std::vector<int64_t> range64(range.size());
        for(int j=0;j<(int)range.size(); j++)
          range64[j] = range[j];
        int64_ranges.push_back(range64);


        // Set the VariableInfo
        KTE::VariableInfo out_info;
        out_info.type = Kokkos::Tools::Experimental::ValueType::kokkos_value_int64;
        out_info.category = Kokkos::Tools::Experimental::StatisticalCategory::kokkos_value_ordinal;
        out_info.valueQuantity = Kokkos::Tools::Experimental::CandidateValueType::kokkos_value_set;
        out_info.candidates = Kokkos::Tools::Experimental::make_candidate_set(3,int64_ranges[int64_ranges.size()-1].data());
        size_t var_id = Kokkos::Tools::Experimental::declare_output_type(muelu_param,out_info);
        out_variables.push_back(KTE::make_variable_value(var_id,int64_t(guess)));

        out_typenames.push_back("int");
      }
      else if(sublist.isType<double>("initial guess")) {
        // Continuous range
        double guess = sublist.get<double>("initial guess");
        const Teuchos::Array<double> & range = sublist.get<Teuchos::Array<double> >("continuous range");
        TEUCHOS_TEST_FOR_EXCEPTION(range.size() !=3, Exceptions::RuntimeError, "MueLu::KokkosTuningInterface: 'continuous range' needs to be (low, high, step)");


        // Copy the range because Kokkos
        std::vector<double> rangeD(range.size());
        for(int j=0;j<(int)range.size(); j++)
          rangeD[j] = range[j];
        double_ranges.push_back(rangeD);


        // Set the VariableInfo
        KTE::VariableInfo out_info;
        out_info.type = Kokkos::Tools::Experimental::ValueType::kokkos_value_double;
        out_info.category = Kokkos::Tools::Experimental::StatisticalCategory::kokkos_value_ordinal;
        out_info.valueQuantity = Kokkos::Tools::Experimental::CandidateValueType::kokkos_value_set;
        // FIXME: This is a memory error waiting to happen
        out_info.candidates = Kokkos::Tools::Experimental::make_candidate_set(3,double_ranges[double_ranges.size()-1].data());
        size_t var_id = Kokkos::Tools::Experimental::declare_output_type(muelu_param,out_info);
        out_variables.push_back(KTE::make_variable_value(var_id,guess));

        out_typenames.push_back("double");
      }
      /*else if(sublist.isType<std::string>("initial guess")) {
        // Categorical

        }*/
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::KokkosTuningInterface: Only handles int and double parameters");
      }

      // Stash the parameter name
      out_names.push_back(muelu_param);

    }// end if pL.isSublist

  }// end for

  // Sanity check
  TEUCHOS_TEST_FOR_EXCEPTION(out_names.size() != out_variables.size(), Exceptions::RuntimeError, "MueLu::KokkosTuningInterface: Error in option processing");
  TEUCHOS_TEST_FOR_EXCEPTION(out_names.size() != out_typenames.size(), Exceptions::RuntimeError, "MueLu::KokkosTuningInterface: Error in option processing");

  /********************************/
  /* Process the input variables  */
  /********************************/
  in_variables.clear();

  const Teuchos::Array<std::string> & inputs = pL.get<Teuchos::Array<std::string> >("input_variables");
  for (int i=0; i< (int)inputs.size(); i++) {
    in_variables.push_back(KTE::make_variable_value(i,inputs[i].c_str()));
  }


}


// ***********************************************************************
  void KokkosTuningInterface::SetMueLuParameters(size_t kokkos_context_id, Teuchos::ParameterList& mueluParams, bool overwrite) const {
  namespace KTE = Kokkos::Tools::Experimental;
  Teuchos::ParameterList tunedParams;

  if (comm_->getRank() == 0) {
    // Only Rank 0 calls KokkosTuning
    GetOStream(Runtime0) << "MueLu::KokkosTuningInterface: Tuning " << out_variables.size() << " parameters" <<std::endl;

    // Start the tuning
    KTE::request_output_values(kokkos_context_id, out_variables.size(), out_variables.data());

    // Unpack the tuned values
    for(int i=0; i<(int)out_names.size(); i++) {
      // FIXME: Allow for '/' separated sublisted params
      if (out_typenames[i] == "int")
        tunedParams.set(out_names[i], (int) out_variables[i].value.int_value);
      else if (out_typenames[i] == "double")
        tunedParams.set(out_names[i], out_variables[i].value.double_value);
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::KokkosTuningInterface: Unknown variable output type");
      }
    }
  }

  Teuchos::updateParametersAndBroadcast(outArg(tunedParams), outArg(mueluParams), *comm_, 0, overwrite);
}


}  // namespace MueLu
