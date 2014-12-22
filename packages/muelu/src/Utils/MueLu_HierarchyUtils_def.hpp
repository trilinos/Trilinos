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
#ifndef MUELU_HIERARCHYUTILS_DEF_HPP
#define MUELU_HIERARCHYUTILS_DEF_HPP
#include <MueLu_HierarchyManager.hpp>
#include <MueLu_HierarchyManager.hpp>
#include <MueLu_FactoryManager.hpp>
#include <MueLu_Hierarchy.hpp>

namespace MueLu {

  
  /* Adds the following non-serializable data (A,P,R,Nullspace,Coordinates) from level-specific sublist nonSerialList,
     calling AddNewLevel as appropriate.
  */
  template<class SC, class LO, class GO, class NO>
  void HierarchyUtils<SC, LO, GO, NO>::AddNonSerializableDataToHierarchy(MueLu::HierarchyManager<SC,LO,GO,NO> & HM, Hierarchy & H, const Teuchos::ParameterList & List) {
    using Teuchos::ParameterList;
    ParameterList dummy;

    for(ParameterList::ConstIterator it = List.begin(); it!=List.end(); it++) {
      // Check for mach of the form "levelX" where X is a positive integer
      if(List.isSublist(it->first) && it->first.find("level ")==0) {
	std::string levelstr = it->first.substr(6,std::string::npos);
	int id = (int) strtol(levelstr.c_str(),0,0);
	if(id > 0)  {
	  // Do enough level adding so we can be sure to add the data to the right place
	  for(int i=H.GetNumLevels(); i<=id; i++)
	      H.AddNewLevel();
	  RCP<FactoryManager> Mfact = rcp(new FactoryManager());

	  // Grab the level sublist & loop over parameters
	  const ParameterList & sublist = List.sublist(it->first);
	  for(ParameterList::ConstIterator it2 = sublist.begin(); it2!=sublist.end(); it2++) {	   
	    if(!it2->first.compare("A") || !it2->first.compare("R") || !it2->first.compare("P")) {
	      H.GetLevel(id)->Set(it2->first,Teuchos::getValue<RCP<Matrix > >(it2->second));
	      Mfact->SetFactory(it2->first,MueLu::NoFactory::getRCP());
	    }	    
	    else if (!it2->first.compare("Nullspace") || !it2->first.compare("Coordinates")) {
	      H.GetLevel(id)->Set(it2->first,Teuchos::getValue<RCP<MultiVector > >(it2->second));
	      Mfact->SetFactory(it2->first,MueLu::NoFactory::getRCP());
	    }
	    else
	      throw std::runtime_error("MueLu::Utils::AddNonSerializableDataToHierarchy: List contains unknown data type");
	  }
	  HM.AddFactoryManager(id,1,Mfact);
	}    
      }
    }
  }




}// end namespace
#endif
