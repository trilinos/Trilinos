// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MueLu_AvatarInterface.hpp"

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "MueLu_BaseClass.hpp"
#include "Teuchos_RawParameterListHelpers.hpp"

// ***********************************************************************
/* Notional Parameterlist Structure
   "avatar: filestem"              "{'mystem1','mystem2'}"
   "avatar: decision tree files"   "{'mystem1.trees','mystem2.trees'}"
   "avatar: names files"           "{'mystem1.names','mystem2.names'}"
   "avatar: good class"            "1"
   "avatar: heuristic" 		   "1"
   "avatar: bounds file"           "{'bounds.data'}"
   "avatar: muelu parameter mapping"
     - "param0'
       - "muelu parameter"          "aggregation: threshold"
       - "avatar parameter"         "DROP_TOL"
       - "muelu values"             "{0,1e-4,1e-3,1e-2}"
       - "avatar values"            "{1,2,3,4}
     - "param1'
       - "muelu parameter"          "smoother: sweeps"
       - "avatar parameter"         "SWEEPS"
       - "muelu values"             "{1,2,3}"
       - "avatar values"            "{1,2,3}"


   Notional SetMueLuParameters "problemFeatures"  Structure
   "my feature #1"                "246.01"
   "my feature #2"                "123.45"

 */

/*
TODO List:
Modify MueLu
    Parameter name checking (make sure names match between Avatar’s names file and the parameter / feature names that MueLu sees).
*/

#ifdef HAVE_MUELU_AVATAR

extern "C" {
#include "avatar_api.h"
}

namespace MueLu {

// ***********************************************************************
RCP<const ParameterList> AvatarInterface::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  Teuchos::ParameterList pl_dummy;
  Teuchos::Array<std::string> ar_dummy;
  int int_dummy;

  // Files from which to read Avatar trees
  validParamList->set<Teuchos::Array<std::string>>("avatar: decision tree files", ar_dummy, "Names of Avatar decision tree files");

  // Strings from which to read Avatar trees
  validParamList->set<Teuchos::Array<std::string>>("avatar: decision tree strings", ar_dummy, "Avatar decision tree strings");

  // Files from which to read Avatar names
  validParamList->set<Teuchos::Array<std::string>>("avatar: names files", ar_dummy, "Names of Avatar decision names files");

  // Strings from which to read Avatar names
  validParamList->set<Teuchos::Array<std::string>>("avatar: names strings", ar_dummy, "Avatar decision names strings");

  // Filestem arg for Avatar
  validParamList->set<Teuchos::Array<std::string>>("avatar: filestem", ar_dummy, "Filestem for the files Avatar requires");

  // This should be a MueLu parameter-to-Avatar parameter mapping (e.g. if Avatar doesn't like spaces)
  validParamList->set<Teuchos::ParameterList>("avatar: muelu parameter mapping", pl_dummy, "Mapping of MueLu to Avatar Parameters");

  // "Good" Class ID for Avatar
  validParamList->set<int>("avatar: good class", int_dummy, "Numeric code for class Avatar considers to be good");

  // Which drop tol choice heuristic to use
  validParamList->set<int>("avatar: heuristic", int_dummy, "Numeric code for which heuristic we want to use");

  // Bounds file for extrapolation risk
  validParamList->set<Teuchos::Array<std::string>>("avatar: bounds file", ar_dummy, "Bounds file for Avatar extrapolation risk");

  // Add dummy variables at the start
  validParamList->set<int>("avatar: initial dummy variables", int_dummy, "Number of dummy variables to add at the start");

  // Add dummy variables before the class
  validParamList->set<int>("avatar: pre-class dummy variables", int_dummy, "Number of dummy variables to add at the before the class");

  // Value of the dummy variables
  validParamList->set<int>("avatar: dummy value", int_dummy, "Value of the dummy variables to add at the before/after the class");

  return validParamList;
}

// ***********************************************************************
Teuchos::ArrayRCP<std::string> AvatarInterface::ReadFromFiles(const char* paramName) const {
  //  const Teuchos::Array<std::string> & tf = params_.get<const Teuchos::Array<std::string> >(paramName);
  Teuchos::Array<std::string>& tf = params_.get<Teuchos::Array<std::string>>(paramName);
  Teuchos::ArrayRCP<std::string> treelist;
  // Only Proc 0 will read the files and print the strings
  if (comm_->getRank() == 0) {
    treelist.resize(tf.size());
    for (Teuchos_Ordinal i = 0; i < tf.size(); i++) {
      std::fstream file;
      std::stringstream ss;
      file.open(tf[i]);
      ss << file.rdbuf();
      treelist[i] = ss.str();
      file.close();
    }
  }
  return treelist;
}

// ***********************************************************************
void AvatarInterface::Setup() {
  // Sanity check
  if (comm_.is_null()) throw std::runtime_error("MueLu::AvatarInterface::Setup(): Communicator cannot be null");

  // Get the avatar strings (NOTE: Only exist on proc 0)
  if (params_.isParameter("avatar: decision tree strings"))
    avatarStrings_ = Teuchos::arcpFromArray(params_.get<Teuchos::Array<std::string>>("avatar: decision tree strings"));
  else
    avatarStrings_ = ReadFromFiles("avatar: decision tree files");

  if (params_.isParameter("avatar: names strings"))
    namesStrings_ = Teuchos::arcpFromArray(params_.get<Teuchos::Array<std::string>>("avatar: names strings"));
  else
    namesStrings_ = ReadFromFiles("avatar: names files");

  if (params_.isParameter("avatar: bounds file"))
    boundsString_ = ReadFromFiles("avatar: bounds file");

  filestem_ = params_.get<Teuchos::Array<std::string>>("avatar: filestem");

  if (comm_->getRank() == 0) {
    // Now actually set up avatar - Avatar's cleanup routine will free the pointer
    //    Avatar_handle* avatar_train(int argc, char **argv, char* names_file, int names_file_is_a_string, char* train_file, int train_file_is_a_string);
    const int namesfile_is_a_string = 1;
    const int treesfile_is_a_string = 1;
    avatarHandle_                   = avatar_load(const_cast<char*>(filestem_[0].c_str()), const_cast<char*>(namesStrings_[0].c_str()), namesfile_is_a_string, const_cast<char*>(avatarStrings_[0].c_str()), treesfile_is_a_string);
  }

  // Which class does Avatar consider "good"
  avatarGoodClass_ = params_.get<int>("avatar: good class");

  heuristicToUse_ = params_.get<int>("avatar: heuristic");

  // Unpack the MueLu Mapping into something actionable
  UnpackMueLuMapping();
}

// ***********************************************************************
void AvatarInterface::Cleanup() {
  avatar_cleanup(avatarHandle_);
  avatarHandle_ = 0;
}

// ***********************************************************************
void AvatarInterface::GenerateFeatureString(const Teuchos::ParameterList& problemFeatures, std::string& featureString) const {
  // NOTE: Assumes that the features are in the same order Avatar wants them.
  std::stringstream ss;

  // Initial Dummy Variables
  if (params_.isParameter("avatar: initial dummy variables")) {
    int num_dummy   = params_.get<int>("avatar: initial dummy variables");
    int dummy_value = params_.get("avatar: dummy value", 666);

    for (int i = 0; i < num_dummy; i++)
      ss << dummy_value << ",";
  }

  for (Teuchos::ParameterList::ConstIterator i = problemFeatures.begin(); i != problemFeatures.end(); i++) {
    //    const std::string& name = problemFeatures.name(i);
    const Teuchos::ParameterEntry& entry = problemFeatures.entry(i);
    if (i != problemFeatures.begin()) ss << ",";
    entry.leftshift(ss, false);  // Because ss<<entry prints out '[unused]' and we don't want that.
  }
  featureString = ss.str();
}

// ***********************************************************************
void AvatarInterface::UnpackMueLuMapping() {
  const Teuchos::ParameterList& mapping = params_.get<Teuchos::ParameterList>("avatar: muelu parameter mapping");
  // Each MueLu/Avatar parameter pair gets its own sublist.  These must be numerically ordered with no gap

  bool done     = false;
  int idx       = 0;
  int numParams = mapping.numParams();

  mueluParameterName_.resize(numParams);
  avatarParameterName_.resize(numParams);
  mueluParameterValues_.resize(numParams);
  avatarParameterValues_.resize(numParams);

  while (!done) {
    std::stringstream ss;
    ss << "param" << idx;
    if (mapping.isSublist(ss.str())) {
      const Teuchos::ParameterList& sublist = mapping.sublist(ss.str());

      // Get the names
      mueluParameterName_[idx]  = sublist.get<std::string>("muelu parameter");
      avatarParameterName_[idx] = sublist.get<std::string>("avatar parameter");

      // Get the values
      // FIXME: For now we assume that all of these guys are doubles and their Avatar analogues are doubles
      mueluParameterValues_[idx]  = sublist.get<Teuchos::Array<double>>("muelu values");
      avatarParameterValues_[idx] = sublist.get<Teuchos::Array<double>>("avatar values");

      idx++;
    } else {
      done = true;
    }
  }

  if (idx != numParams)
    throw std::runtime_error("MueLu::AvatarInterface::UnpackMueLuMapping(): 'avatar: muelu parameter mapping' has unknown fields");
}
// ***********************************************************************
std::string AvatarInterface::ParamsToString(const std::vector<int>& indices) const {
  std::stringstream ss;
  for (Teuchos_Ordinal i = 0; i < avatarParameterValues_.size(); i++) {
    ss << "," << avatarParameterValues_[i][indices[i]];
  }

  // Pre-Class dummy variables
  if (params_.isParameter("avatar: pre-class dummy variables")) {
    int num_dummy   = params_.get<int>("avatar: pre-class dummy variables");
    int dummy_value = params_.get("avatar: dummy value", 666);
    for (int i = 0; i < num_dummy; i++)
      ss << "," << dummy_value;
  }

  return ss.str();
}

// ***********************************************************************
void AvatarInterface::SetIndices(int id, std::vector<int>& indices) const {
  // The combo numbering here starts with the first guy
  int numParams = (int)avatarParameterValues_.size();
  int curr_id   = id;
  for (int i = 0; i < numParams; i++) {
    int div    = avatarParameterValues_[i].size();
    int mod    = curr_id % div;
    indices[i] = mod;
    curr_id    = (curr_id - mod) / div;
  }
}

// ***********************************************************************
void AvatarInterface::GenerateMueLuParametersFromIndex(int id, Teuchos::ParameterList& pl) const {
  // The combo numbering here starts with the first guy
  int numParams = (int)avatarParameterValues_.size();
  int curr_id   = id;
  for (int i = 0; i < numParams; i++) {
    int div = avatarParameterValues_[i].size();
    int mod = curr_id % div;
    pl.set(mueluParameterName_[i], mueluParameterValues_[i][mod]);
    curr_id = (curr_id - mod) / div;
  }
}

// ***********************************************************************
void AvatarInterface::SetMueLuParameters(const Teuchos::ParameterList& problemFeatures, Teuchos::ParameterList& mueluParams, bool overwrite) const {
  Teuchos::ParameterList avatarParams;
  std::string paramString;

  if (comm_->getRank() == 0) {
    // Only Rank 0 calls Avatar
    if (!avatarHandle_) throw std::runtime_error("MueLu::AvatarInterface::SetMueLuParameters(): Setup has not been run");

    // Turn the problem features into a "trial" string for Avatar
    std::string trialString;
    GenerateFeatureString(problemFeatures, trialString);

    // Compute the number of things we need to test
    int numParams = (int)avatarParameterValues_.size();
    std::vector<int> indices(numParams);
    std::vector<int> sizes(numParams);
    int num_combos = 1;
    for (int i = 0; i < numParams; i++) {
      sizes[i] = avatarParameterValues_[i].size();
      num_combos *= avatarParameterValues_[i].size();
    }
    GetOStream(Runtime0) << "MueLu::AvatarInterface: Testing " << num_combos << " option combinations" << std::endl;

    // For each input parameter to avatar we iterate over its allowable values and then compute the list of options which Avatar
    // views as acceptable
    // FIXME: Find alternative to hard coding malloc size (input deck?)
    int num_classes = avatar_num_classes(avatarHandle_);
    std::vector<int> predictions(num_combos, 0);
    std::vector<float> probabilities(num_classes * num_combos, 0);

    std::string testString;
    for (int i = 0; i < num_combos; i++) {
      SetIndices(i, indices);
      // Now we add the MueLu parameters into one, enormous Avatar trial string and run avatar once
      testString += trialString + ParamsToString(indices) + ",0\n";
    }

    std::cout << "** Avatar TestString ***\n"
              << testString << std::endl;  // DEBUG

    int bound_check = true;
    if (params_.isParameter("avatar: bounds file"))
      bound_check = checkBounds(testString, boundsString_);

    // FIXME: Only send in first tree's string
    // int* avatar_test(Avatar_handle* a, char* test_data_file, int test_data_is_a_string);
    const int test_data_is_a_string = 1;
    avatar_test(avatarHandle_, const_cast<char*>(testString.c_str()), test_data_is_a_string, predictions.data(), probabilities.data());

    // Look at the list of acceptable combinations of options
    std::vector<int> acceptableCombos;
    acceptableCombos.reserve(100);
    for (int i = 0; i < num_combos; i++) {
      if (predictions[i] == avatarGoodClass_) acceptableCombos.push_back(i);
    }
    GetOStream(Runtime0) << "MueLu::AvatarInterface: " << acceptableCombos.size() << " acceptable option combinations found" << std::endl;

    // Did we have any good combos at all?
    int chosen_option_id = 0;
    if (acceptableCombos.size() == 0) {
      GetOStream(Runtime0) << "WARNING: MueLu::AvatarInterface: found *no* combinations of options which it believes will perform well on this problem" << std::endl
                           << "         An arbitrary set of options will be chosen instead" << std::endl;
    } else {
      // If there is only one acceptable combination, use it;
      // otherwise, find the parameter choice with the highest
      // probability of success
      if (acceptableCombos.size() == 1) {
        chosen_option_id = acceptableCombos[0];
      } else {
        switch (heuristicToUse_) {
          case 1:
            chosen_option_id = hybrid(probabilities.data(), acceptableCombos);
            break;
          case 2:
            chosen_option_id = highProb(probabilities.data(), acceptableCombos);
            break;
          case 3:
            // Choose the first option in the list of acceptable
            // combinations; the lowest drop tolerance among the
            // acceptable combinations
            chosen_option_id = acceptableCombos[0];
            break;
          case 4:
            chosen_option_id = lowCrash(probabilities.data(), acceptableCombos);
            break;
          case 5:
            chosen_option_id = weighted(probabilities.data(), acceptableCombos);
            break;
        }
      }
    }

    // If mesh parameters are outside bounding box, set drop tolerance
    // to 0, otherwise use avatar recommended drop tolerance
    if (bound_check == 0) {
      GetOStream(Runtime0) << "WARNING: Extrapolation risk detected, setting drop tolerance to 0" << std::endl;
      GenerateMueLuParametersFromIndex(0, avatarParams);
    } else {
      GenerateMueLuParametersFromIndex(chosen_option_id, avatarParams);
    }
  }

  Teuchos::updateParametersAndBroadcast(outArg(avatarParams), outArg(mueluParams), *comm_, 0, overwrite);
}

int AvatarInterface::checkBounds(std::string trialString, Teuchos::ArrayRCP<std::string> boundsString) const {
  std::stringstream ss(trialString);
  std::vector<double> vect;

  double b;
  while (ss >> b) {
    vect.push_back(b);
    if (ss.peek() == ',') ss.ignore();
  }

  std::stringstream ssBounds(boundsString[0]);
  std::vector<double> boundsVect;

  while (ssBounds >> b) {
    boundsVect.push_back(b);
    if (ssBounds.peek() == ',') ssBounds.ignore();
  }

  int min_idx = (int)std::min(vect.size(), boundsVect.size() / 2);

  bool inbounds = true;
  for (int i = 0; inbounds && i < min_idx; i++)
    inbounds = boundsVect[2 * i] <= vect[i] && vect[i] <= boundsVect[2 * i + 1];

  return (int)inbounds;
}

int AvatarInterface::hybrid(float* probabilities, std::vector<int> acceptableCombos) const {
  float low_crash = probabilities[0];
  float best_prob = probabilities[2];
  float diff;
  int this_combo;
  int chosen_option_id = acceptableCombos[0];
  for (int x = 0; x < (int)acceptableCombos.size(); x++) {
    this_combo = acceptableCombos[x] * 3;
    diff       = probabilities[this_combo] - low_crash;
    // If this parameter combination has a crash
    // probability .2 lower than the current "best", we
    // will use this drop tolerance
    if (diff < -.2) {
      low_crash        = probabilities[this_combo];
      best_prob        = probabilities[this_combo + 2];
      chosen_option_id = acceptableCombos[x];
    }
    // If this parameter combination has the same
    // or slightly lower crash probability than the
    // current best, we compare their "GOOD" probabilities
    else if (diff <= 0 && probabilities[this_combo + 2] > best_prob) {
      low_crash        = probabilities[this_combo];
      best_prob        = probabilities[this_combo + 2];
      chosen_option_id = acceptableCombos[x];
    }
  }
  return chosen_option_id;
}

int AvatarInterface::highProb(float* probabilities, std::vector<int> acceptableCombos) const {
  float high_prob = probabilities[2];
  int this_combo;
  int chosen_option_id = acceptableCombos[0];
  for (int x = 0; x < (int)acceptableCombos.size(); x++) {
    this_combo = acceptableCombos[x] * 3;
    // If this parameter combination has a higher "GOOD"
    // probability, use this combination
    if (probabilities[this_combo + 2] > high_prob) {
      high_prob        = probabilities[this_combo + 2];
      chosen_option_id = acceptableCombos[x];
    }
  }
  return chosen_option_id;
}

int AvatarInterface::lowCrash(float* probabilities, std::vector<int> acceptableCombos) const {
  float low_crash = probabilities[0];
  int this_combo;
  int chosen_option_id = acceptableCombos[0];
  for (int x = 0; x < (int)acceptableCombos.size(); x++) {
    this_combo = acceptableCombos[x] * 3;
    // If this parameter combination has a lower "CRASH"
    // probability, use this combination
    if (probabilities[this_combo] < low_crash) {
      low_crash        = probabilities[this_combo];
      chosen_option_id = acceptableCombos[x];
    }
  }
  return chosen_option_id;
}

int AvatarInterface::weighted(float* probabilities, std::vector<int> acceptableCombos) const {
  float low_crash = probabilities[0];
  float best_prob = probabilities[2];
  float diff;
  int this_combo;
  int chosen_option_id = acceptableCombos[0];
  for (int x = 0; x < (int)acceptableCombos.size(); x++) {
    this_combo = acceptableCombos[x] * 3;
    diff       = probabilities[this_combo] - low_crash;
    // If this parameter combination has a crash
    // probability .2 lower than the current "best", we
    // will use this drop tolerance
    if (diff < -.2) {
      low_crash        = probabilities[this_combo];
      best_prob        = probabilities[this_combo + 2];
      chosen_option_id = acceptableCombos[x];
    }
    // If this parameter combination is within .1
    // or has a slightly lower crash probability than the
    // current best, we compare their "GOOD" probabilities
    else if (diff <= .1 && probabilities[this_combo + 2] > best_prob) {
      low_crash        = probabilities[this_combo];
      best_prob        = probabilities[this_combo + 2];
      chosen_option_id = acceptableCombos[x];
    }
  }
  return chosen_option_id;
}

}  // namespace MueLu

#endif  // HAVE_MUELU_AVATAR
