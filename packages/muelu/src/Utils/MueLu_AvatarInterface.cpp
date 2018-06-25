// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include "MueLu_AvatarInterface.hpp"

#include <string> 
#include <fstream> 
#include <sstream> 
#include <vector> 
#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "MueLu_BaseClass.hpp"


#ifdef HAVE_MUELU_AVATAR
namespace MueLu {

// ***********************************************************************
RCP<const ParameterList> AvatarInterface::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  Teuchos::ParameterList pl_dummy;
  Teuchos::Array<std::string> ar_dummy;
  std::string s_dummy;

  // Files from which to read Avatar trees
  validParamList->set<Teuchos::Array<std::string> >("avatar: decision tree files",ar_dummy,"Names of Avatar decision tree files");

  // This is a general, app-specified parameter list
  validParamList->set<Teuchos::ParameterList>("avatar: problem features",pl_dummy,"Problem features as computed by the application");

  // Avatar command line args
  validParamList->set<std::string>("avatar: args",s_dummy,"Arguments to control the execution of Avatar");

  // This should be a MueLu parameter-to-Avatar parameter mapping (e.g. if Avatar doesn't like spaces)
  validParamList->set<Teuchos::ParameterList>("avatar: muelu parameter mapping",pl_dummy,"Mapping of MueLu to Avatar Parameters");
  
  return validParamList;
}


// ***********************************************************************
Teuchos::ArrayRCP<std::string> AvatarInterface::ReadAvatarStringsFromFiles() const {
  const Teuchos::Array<std::string> & tf = params_.get<const Teuchos::Array<std::string> >("avatar: decision tree files");
  Teuchos::ArrayRCP<std::string> treelist;
  // Only Proc 0 will read the files and print the strings
  if (comm_->getRank() == 0) {
    treelist.resize(tf.size());
    for(Teuchos_Ordinal i=0; i<tf.size(); i++) {
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
  if(comm_.is_null()) throw std::runtime_error("MueLu::AvatarInterface::Setup(): Communicator cannot be null");

  // Get the avatar strings (NOTE: Only exist on proc 0)
  avatarStrings_ = ReadAvatarStringsFromFiles();

  // Now actually set up avatar
  std::string avatarArgs = params_.get("avatar: args",std::string(""));
#if 0
  if(comm_->getRank() == 0)
    avatarHandle_ = Teuchos::rcp(new Avatar(avatarArgs));
#endif

  // Unpack the MueLu Mapping into something actionable
  UnpackMueLuMapping();

  throw std::runtime_error("Not yet implemented!");

}

// ***********************************************************************
void AvatarInterface::GenerateFeatureString(const Teuchos::ParameterList & problemFeatures, std::string & featureString) const {
  // NOTE: Assumes that the features are in the same order Avatar wants them.
  // FIXME: This ordering should be checked against the names file
  std::stringstream ss;
  for(Teuchos::ParameterList::ConstIterator i=problemFeatures.begin(); i != problemFeatures.end(); i++) {
    //    const std::string& name = problemFeatures.name(i);
    const Teuchos::ParameterEntry& entry = problemFeatures.entry(i);
    if(i!=problemFeatures.begin()) ss<<",";
    ss<<entry;
  }
  featureString = ss.str();
}

// ***********************************************************************
void AvatarInterface::UnpackMueLuMapping() {
  const Teuchos::ParameterList & mapping = params_.get<Teuchos::ParameterList>("avatar: muelu parameter mapping");
  // Each MueLu/Avatar parameter pair gets its own sublist.  These must be numerically ordered with no gap

  bool done=false; 
  int idx=0;
  int numParams = mapping.numParams();

  mueluParameterName_.resize(numParams);
  avatarParameterName_.resize(numParams);
  mueluParameterValues_.resize(numParams);
  avatarParameterValues_.resize(numParams);

  while(!done) {
    std::stringstream ss; 
    ss << "param" << idx;
    if(mapping.isSublist(ss.str())) {
      const Teuchos::ParameterList & sublist = mapping.sublist(ss.str());

      // Get the names
      mueluParameterName_[idx] = sublist.get<std::string>("muelu parameter");
      avatarParameterName_[idx] = sublist.get<std::string>("avatar parameter");

      // Get the values
      //FIXME: For now we assume that all of these guys are doubles and their Avatar analogues are ints
      mueluParameterValues_[idx]  = sublist.get<Teuchos::Array<double> >("muelu values");
      avatarParameterValues_[idx] = sublist.get<Teuchos::Array<int> >("avatar values");            
    }
    else {
      done=true;
    }
    idx++;
  }
  
  if(idx!=numParams) 
    throw std::runtime_error("MueLu::AvatarInterface::UnpackMueLuMapping(): 'avatar: muelu parameter mapping' has unknown fields");
}
// ***********************************************************************
std::string AvatarInterface::ParamsToString(const std::vector<int> & indices) const {
  std::stringstream ss;
  for(Teuchos_Ordinal i=0; i<avatarParameterValues_.size(); i++) {
    ss << "," << avatarParameterValues_[i][indices[i]];
  }
  return ss.str();
}

// ***********************************************************************
void AvatarInterface::SetIndices(int id,std::vector<int> & indices) const {
  // The combo numbering here starts with the first guy
  int numParams = (int)avatarParameterValues_.size();
  int curr_id = id;
  for(int i=0; i<numParams; i++) {
    int div = avatarParameterValues_[i].size();
    int mod = curr_id % div;
    indices[i] = mod;
    curr_id = (curr_id - mod)/div;
  }
}



// ***********************************************************************
void AvatarInterface::GenerateMueLuParametersFromIndex(int id,Teuchos::ParameterList & pl) const {
  // The combo numbering here starts with the first guy
  int numParams = (int)avatarParameterValues_.size();
  int curr_id = id;
  for(int i=0; i<numParams; i++) {
    int div = avatarParameterValues_[i].size();
    int mod = curr_id % div;
    pl.set(mueluParameterName_[i],mueluParameterValues_[mod]);
    curr_id = (curr_id - mod)/div;
  }
}


// ***********************************************************************
void AvatarInterface::SetMueLuParameters(const Teuchos::ParameterList & problemFeatures, Teuchos::ParameterList & mueluParams, bool overwrite) const {
  Teuchos::ParameterList avatarParams;
  std::string paramString;

  if (comm_->getRank() == 0) {
    // Only Rank 0 calls Avatar
    if(avatarHandle_.is_null()) throw std::runtime_error("MueLu::AvatarInterface::SetMueLuParameters(): Setup has not been run");

    // Turn the problem features into a "trial" string for Avatar
    std::string trialString;
    GenerateFeatureString(problemFeatures,trialString);
    
    // Compute the number of things we need to test
    int numParams = (int)avatarParameterValues_.size();
    std::vector<int> indices(numParams);
    std::vector<int> sizes(numParams);
    int num_combos = 1;
    for(int i=0; i<numParams; i++) {
      sizes[i]    = avatarParameterValues_[i].size();
      num_combos *= avatarParameterValues_[i].size();
    }
    GetOStream(Runtime0)<< "MueLu::AvatarInterface: Testing "<< num_combos << " option combinations"<<std::endl;

    // For each input parameter to avatar we iterate over its allowable values and then compute the list of options which Avatar
    // views as acceptable
    std::vector<int> avatarOutput;
    {
      std::string avatarString;
      for(int i=0; i<num_combos; i++) {
        SetIndices(i,indices);
        // Now we add the MueLu parameters into one, enormous Avatar trial string and run avatar once
        avatarString += trialString + ParamsToString(indices) + "\n";
      }
      
#if 0
      avatarOutput= avatarHandle_->evaluate(avatarString);
#endif      
    }

    // Look at the list of acceptable combinations of options 
    std::vector<int> acceptableCombos; acceptableCombos.reserve(100);
    for(int i=0; i<num_combos; i++) {    
      if(avatarOutput[i] == 1) acceptableCombos.push_back(i);      
    }
    GetOStream(Runtime0)<< "MueLu::AvatarInterface: "<< acceptableCombos.size() << " acceptable option combinations found"<<std::endl;

    // Did we have any good combos at all?
    int chosen_option_id = 0;
    if(acceptableCombos.size() == 0) { 
      GetOStream(Runtime0) << "WARNING: MueLu::AvatarInterface: found *no* combinations of options which it believes will perform well on this problem" <<std::endl
                           << "         An arbitrary set of options will be chosen instead"<<std::endl;    
    }
    else {
      // As a placeholder, we'll choose the first acceptable combination.  Later we can do something smarter
      chosen_option_id = acceptableCombos[0];
    }

    // Generate the parameterList from the chosen option
    GenerateMueLuParametersFromIndex(chosen_option_id,avatarParams);
  }

# if 0
  Teuchos::updateParametersAndBroadcast(outArg(avatarParams),outArg(mueluParams),*comm_,0,overwrite);
#endif


}






} // namespace MueLu

#endif// HAVE_MUELU_AVATAR

