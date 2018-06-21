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
#include "Teuchos_Array.hpp"
#include "MueLu_BaseClass.hpp"


#ifdef HAVE_MUELU_AVATAR
namespace MueLu {

// ***********************************************************************
RCP<const ParameterList> AvatarInterface::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  Teuchos::ParameterList pl_dummy;
  Teuchos::Array<std::string> ar_dummy;

  // Files from which to read Avatar trees
  validParamList->set<Teuchos::Array<std::string> >("avatar: decision tree files",ar_dummy,"Names of Avatar decision tree files");

  // This is a general, app-specified parameter list
  validParamList->set<Teuchos::ParameterList>("avatar: problem features",pl_dummy,"Problem features as computed by the application");

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

  throw std::runtime_error("Not yet implemented!");

}

// ***********************************************************************
void AvatarInterface::SetMueLuParameters(Teuchos::ParameterList & pl) const {
 

  // FIXME: Now actually call Avatar


}






} // namespace MueLu

#endif// HAVE_MUELU_AVATAR

