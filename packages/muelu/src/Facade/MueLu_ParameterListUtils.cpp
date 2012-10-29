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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#include "MueLu_ParameterListUtils.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {
  
  /* See also: ML_Epetra::UpdateList */
  void MergeParameterList(const Teuchos::ParameterList &source, Teuchos::ParameterList &dest, bool overWrite){
    for(Teuchos::ParameterList::ConstIterator param=source.begin(); param!=source.end(); ++param)
      if (dest.isParameter(source.name(param)) == false || overWrite)
        dest.setEntry(source.name(param),source.entry(param));
  }

  void CreateSublists(const Teuchos::ParameterList &List, Teuchos::ParameterList &newList)
  {
    using Teuchos::ParameterList;
    using std::string;

    newList.setName(List.name());

    // Copy general (= not level-specific) options and sublists to the new list.
    // - Coarse and level-specific parameters are not copied yet. They will be moved to sublists later.
    // - Already existing level-specific lists are copied to the new list but the coarse list is not copied 
    //   yet because it has to be modified before copy (s/coarse/smoother/)
    for (ParameterList::ConstIterator param=List.begin(); param!=List.end(); ++param)
      {
        const string & pname=List.name(param);

        if ((pname.find(" (level",0)  == string::npos || pname.find("smoother: list (level",0) == 0 || pname.find("aggregation: list (level",0) == 0) &&
            (pname.find("coarse: ",0) == string::npos))
          {
            newList.setEntry(pname,List.entry(param));
          }
      } // for

    // Copy of the sublist "coarse: list" to the new list. Change "coarse:" to "smoother:" along the way.
    if (List.isSublist("coarse: list")) {
      const ParameterList &coarseList = List.sublist("coarse: list");
      ParameterList &newCoarseList = newList.sublist("coarse: list");
      for (ParameterList::ConstIterator param=coarseList.begin(); param!=coarseList.end() ; ++param) {
        const string & pname=coarseList.name(param);
      
        if (pname.find("coarse:",0) == 0) {
          // change "coarse: " to "smoother:"
          newCoarseList.setEntry("smoother: "+pname.substr(8),coarseList.entry(param));
        } else {
          newCoarseList.setEntry(pname,coarseList.entry(param));
        }
      }
    } // if

    // Copy of level-specific parameters and coarse parameters to sublist
    for (ParameterList::ConstIterator param=List.begin(); param!=List.end(); ++param)
      {
        const string & pname=List.name(param);
        if (pname.find(" (level",0) != string::npos && pname.find("smoother: list (level",0) != 0 && pname.find("aggregation: list (level",0) != 0)
          {
            // Copy level-specific parameters (smoother and aggregation)
              
            // Scan pname (ex: pname="smoother: type (level 2)")
            string type, option;  
            int levelID=-1;
            {
              typedef Teuchos::ArrayRCP<char>::size_type size_type;    // (!)
              Teuchos::Array<char> ctype  (size_type(pname.size()+1));
              Teuchos::Array<char> coption(size_type(pname.size()+1));
              
              int matched = sscanf(pname.c_str(),"%s %[^(](level %d)", ctype.getRawPtr(), coption.getRawPtr(), &levelID); // use [^(] instead of %s to allow for strings with white-spaces (ex: "ifpack list")
              type = string(ctype.getRawPtr());
              option = string(coption.getRawPtr()); option.resize(option.size () - 1); // remove final white-space
              
              if (matched != 3 || (type != "smoother:" && type != "aggregation:")) {
                TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "MueLu::CreateSublist(), Line " << __LINE__ << ". "
                                           << "Error in creating level-specific sublists" << std::endl
                                           << "Offending parameter: " << pname << std::endl);
              }
            }
            
            // Create/grab the corresponding sublist of newList
            ParameterList &newSubList = newList.sublist(type + " list (level " + Teuchos::toString(levelID) + ")");
            // Shove option w/o level number into sublist
            newSubList.setEntry(type + " " + option,List.entry(param));
            
          } else if (pname.find("coarse:",0) == 0 && pname != "coarse: list") {
          // Copy coarse parameters
          ParameterList &newCoarseList = newList.sublist("coarse: list"); // the coarse sublist is created only if there is at least one "coarse:" parameter
          newCoarseList.setEntry("smoother: "+pname.substr(8),List.entry(param)); // change "coarse: " to "smoother:"
        } // end if
        
      } // for
    
  } //MueLu::CreateSublist()
  
  // Usage: GetMLSubList(paramList, "smoother", 2);
  const Teuchos::ParameterList & GetMLSubList(const Teuchos::ParameterList & paramList, const std::string & type, int levelID) {
    static const Teuchos::ParameterList emptyParamList;

    char levelChar[11];
    sprintf(levelChar, "(level %d)", levelID);
    std::string levelStr(levelChar);
    
    if (paramList.isSublist(type + ": list " + levelStr)) {
      return paramList.sublist(type + ": list " + levelStr);
    } else {
      return emptyParamList;
    }
  }

  // Extract all the parameters that begin with "str:" (but skip sublist)
  Teuchos::RCP<Teuchos::ParameterList> ExtractSetOfParameters(const Teuchos::ParameterList & paramList, const std::string & str) {
    Teuchos::RCP<Teuchos::ParameterList> subList = rcp(new Teuchos::ParameterList());
  
    for (Teuchos::ParameterList::ConstIterator param = paramList.begin(); param != paramList.end(); ++param) {
      const Teuchos::ParameterEntry & entry = paramList.entry(param);
      const std::string             & pname = paramList.name(param);
      if (pname.find(str+":",0) == 0 && !entry.isList()) {
        subList->setEntry(pname,entry);
      }
    }

    return subList;
  }

} // namespace MueLu
