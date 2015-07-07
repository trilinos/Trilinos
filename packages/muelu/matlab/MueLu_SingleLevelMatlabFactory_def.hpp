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
#ifndef MUELU_SINGLELEVELMATLABFACTORY_DEF_HPP
#define MUELU_SINGLELEVELMATLABFACTORY_DEF_HPP
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVector.hpp>

#include "MueLu_Aggregates.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_SingleLevelMatlabFactory_decl.hpp"
#include "MueLu_MatlabUtils_decl.hpp"

#ifdef HAVE_MUELU_MATLAB
#include "mex.h"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  SingleLevelMatlabFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SingleLevelMatlabFactory()
    : hasDeclaredInput_(false) { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> SingleLevelMatlabFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
    
    //add all possible input ("needs") data factories as valid factory parameters
    validParamList->set<RCP<const FactoryBase>>("A", Teuchos::null, "Factory for the matrix A.");
    validParamList->set<RCP<const FactoryBase>>("P", Teuchos::null, "Factory for the prolongator.");
    validParamList->set<RCP<const FactoryBase>>("R", Teuchos::null, "Factory for the restrictor.");
    validParamList->set<RCP<const FactoryBase>>("Ptent", Teuchos::null, "Factory for the tentative (unsmoothed) prolongator.");
    validParamList->set<RCP<const FactoryBase>>("Coordinates", Teuchos::null, "Factory for the node coordinates.");
    validParamList->set<RCP<const FactoryBase>>("Nullspace", Teuchos::null, "Factory for the nullspace.");
    validParamList->set<RCP<const FactoryBase>>("Aggregates",  Teuchos::null, "Factory for the aggregates.");
    validParamList->set<RCP<const FactoryBase>>("UnamalgamationInfo", Teuchos::null, "Factory for amalgamation.");   

    validParamList->set<std::string>("Provides"     , "" ,"A comma-separated list of objects provided by the SingleLevelMatlabFactory");
    validParamList->set<std::string>("Needs"        , "", "A comma-separated list of objects needed by the SingleLevelMatlabFactory");
    validParamList->set<std::string>("Function"     , "" , "The name of the Matlab MEX function to call for Build()");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void SingleLevelMatlabFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
    const Teuchos::ParameterList& pL = GetParameterList();

    // Get needs strings
    const std::string str = pL.get<std::string>("Needs");

    needs_.clear();
    // Tokenize the strings
    TokenizeStringAndStripWhiteSpace(str, needs_, " ,;");

    // Declare inputs
    for(size_t i = 0; i < needs_.size(); i++)
      this->Input(currentLevel, needs_[i]);

    hasDeclaredInput_ = true;

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void SingleLevelMatlabFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
    const Teuchos::ParameterList& pL = GetParameterList();
    using Teuchos::rcp;
    using Teuchos::rcp;

    /* NOTE: We need types to call Get functions.  For the moment, I'm just going to infer types from names.
       We'll need to replace this wby modifying the needs list with strings that define types and adding some kind of lookup function instead */

    // NOTE: mexOutput[0] is the "Provides."  Might want to modify to allow for additional outputs
    std::vector<RCP<MuemexArg>> InputArgs;

    // Fine needs
    for(size_t i = 0; i < needs_.size(); i++) {
      if(needs_[i] == "A" || needs_[i] == "P" || needs_[i] == "R" || needs_[i]=="Ptent") {
	RCP<Matrix> mydata = Get<RCP<Matrix>>(currentLevel, needs_[i]);
	InputArgs.push_back(rcp(new MuemexData<RCP<Matrix>>(mydata)));
      }

      if(needs_[i] == "Nullspace" || needs_[i] == "Coordinates") {
	RCP<MultiVector> mydata = Get<RCP<MultiVector>>(currentLevel, needs_[i]);
	InputArgs.push_back(rcp(new MuemexData<RCP<MultiVector> >(mydata)));
      }

      if(needs_[i] == "Aggregates") {
	RCP<Aggregates> mydata = Get<RCP<Aggregates>>(currentLevel, needs_[i]);
	InputArgs.push_back(rcp(new MuemexData<RCP<Aggregates> >(mydata)));
      }

      if(needs_[i] == "UnAmalgamationInfo") {
	RCP<AmalgamationInfo> mydata = Get<RCP<AmalgamationInfo> >(currentLevel,needs_[i]);
	InputArgs.push_back(rcp(new MuemexData<RCP<AmalgamationInfo>>(mydata)));
      }
    }

    // Determine output
    const std::string str_prov = pL.get<std::string>("Provides");
    std::vector<std::string> provides;
    TokenizeStringAndStripWhiteSpace(str_prov, provides, " ,;");

    // Call mex function
    std::string matlabFunction = pL.get<std::string>("Function");
    if(!matlabFunction.length())
      throw std::runtime_error("Invalid matlab function name");
    std::vector<Teuchos::RCP<MuemexArg> > mexOutput = callMatlab(matlabFunction, provides.size(), InputArgs);

    // Set output
    for(size_t i = 0; i < provides.size(); i++)  {
      if(provides[i] == "A" || provides[i] == "P" || provides[i] == "R" || provides[i]=="Ptent") {
	RCP<MuemexData<RCP<Matrix> > > mydata = Teuchos::rcp_static_cast<MuemexData<RCP<Matrix> > >(mexOutput[i]);
	Set(currentLevel, provides[i], mydata->getData());
      }

      if(provides[i] == "Nullspace" || provides[i] == "Coordinates") {
	RCP<MuemexData<RCP<MultiVector> > > mydata = Teuchos::rcp_static_cast<MuemexData<RCP<MultiVector> > >(mexOutput[i]);
	Set(currentLevel, provides[i], mydata->getData());
      }

      if(provides[i] == "Aggregates") {
	RCP<MuemexData<RCP<Aggregates> > > mydata = Teuchos::rcp_static_cast<MuemexData<RCP<Aggregates> > >(mexOutput[i]);
	Set(currentLevel, provides[i], mydata->getData());
      }

      if(provides[i] == "UnAmalgamationInfo") {
	RCP<MuemexData<RCP<AmalgamationInfo> > > mydata = Teuchos::rcp_static_cast<MuemexData<RCP<AmalgamationInfo> > >(mexOutput[i]);
	Set(currentLevel, provides[i], mydata->getData());
      }
    }
  }

} //namespace MueLu

#define MUELU_SINGLELEVELMATLABFACTORY_SHORT
#endif // HAVE_MUELU_MATLAB

#endif // MUELU_SINGLELEVELMATLABFACTORY_DEF_HPP

