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
#include <regex>
 
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "MueLu_BaseClass.hpp"


#ifdef HAVE_MUELU_AVATAR
namespace MueLu {

// ***********************************************************************
RCP<const ParameterList> AvatarInterface::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  Teuchos::ParameterList pl_dummy;
  // Parameters for Avatar Input
  validParamList->set<Teuchos::Array<std::string> >("avatar: decision tree files",pl_dummy,"Names of Avatar decision tree files");

  // This is a general, app-specified parameter list
  validParamList->set<Teuchos::ParameterList>("avatar: problem features",pl_dummy,"Problem features as computed by the application");

  // This should be a MueLu parameter-to-Avatar parameter mapping (e.g. if Avatar doesn't like spaces)
  validParamList->set<Teuchos::ParameterList>("avatar: muelu parameter mapping",pl_dummy,"Mapping of MueLu to Avatar Parameters");
  
  return validParamList;
}


// ***********************************************************************
Teuchos::ArrayRCP<std::string> AvatarInterface::ReadAvatarStringsFromFiles() const {
  const Teuchos::ParameterList & pl = params_.get<const Teuchos::ParameterList>("avatar: decision tree files");
  std::regex str_regex("file([0-9]+)");
  Teuchos::ArrayRCP<std::string> treelist(1);

  // Avatar files are named "fileXX" where XX is a non-negative integer
  if (comm_.getRank() == 0) {
    // Read

  }
   
  if (comm_.getSize() > 1) {
    // Broadcast
  }
throw std::runtime_error("Not yet implemented!");
}


// ***********************************************************************
Teuchos::ArrayRCP<std::string> AvatarInterface::ReadAvatarStringsFromParameterList() const {
  const Teuchos::ParameterList & pl = params_.get<const Teuchos::ParameterList>("avatar: decision tree strings");
  std::regex str_regex("string([0-9]+)");

  Teuchos::ArrayRCP<std::string> treelist(1);

  // Avatar strings are named "stringXX" where XX is a non-negative integer
  for(Teuchos::ParameterList::ConstIterator i=pl.begin(); i != pl.end(); i++) {
    const std::string& name = pl.name(i);
    const Teuchos::ParameterEntry& entry = pl.entry(i);
    std::smatch num_match;
 
    if(std::regex_match(name, num_match, str_regex)) {
      if (num_match.size() == 2) {
        int strno = std::stoi(num_match[1].str());
        if (strno < 0) throw std::runtime_error("MueLu::AvatarInterface::SetMueLuParameters(): Tree numbers must be positive");       
        if(treelist.size() < strno) treelist.resize(strno+1);
        treelist[strno] = Teuchos::getValue<std::string>(entry);
      } 
      else 
        throw std::runtime_error("MueLu::AvatarInterface::SetMueLuParameters(): Invalid parameters on avatar string list");
    }
    else 
      throw std::runtime_error("MueLu::AvatarInterface::SetMueLuParameters(): Invalid parameters on avatar string list");
  }

  return treelist;
}


// ***********************************************************************
void AvatarInterface::SetMueLuParameters(Teuchos::ParameterList & pl) const {

  bool have_strings = params_.isParameter("avatar: decision tree strings");
  bool have_files   = params_.isParameter("avatar: decision tree files");

  // Sanity Checks
  if(have_strings && have_files) 
    throw std::runtime_error("MueLu::AvatarInterface::SetMueLuParameters():  Cannot specify both 'avatar: decision tree strings' and 'avatar: decision tree file'");

  // Get the avatar strings
  Teuchos::ArrayRCP<std::string> avatar_strings;
  if(have_files)   avatar_strings = ReadAvatarStringsFromFiles();
  if(have_strings) avatar_strings = ReadAvatarStringsFromParameterList();


  // FIXME: Now actually call Avatar


}






} // namespace MueLu

#endif// HAVE_MUELU_AVATAR

