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

#include <Xpetra_Matrix.hpp>
#include <Xpetra_Operator.hpp>

#include "MueLu_HierarchyUtils_decl.hpp"
#include "MueLu_HierarchyManager.hpp"
#include "MueLu_SmootherBase.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_FactoryManager.hpp"

//TODO/FIXME: DeclareInput(, **this**) cannot be used here
#ifdef HAVE_MUELU_INTREPID2
#include "Kokkos_DynRankView.hpp"
#endif

namespace MueLu {


  // Adds the following non-serializable data (A,P,R,Nullspace,Coordinates) from level-specific sublist nonSerialList,
  // calling AddNewLevel as appropriate.
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void HierarchyUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddNonSerializableDataToHierarchy(HierarchyManager& HM, Hierarchy& H, const ParameterList& paramList) {

    for (ParameterList::ConstIterator it = paramList.begin(); it != paramList.end(); it++) {
      const std::string& levelName = it->first;
      // Check for mach of the form "level X" where X is a positive integer
      if (paramList.isSublist(levelName) && levelName.find("level ") == 0 && levelName.size() > 6) {
        int levelID = strtol(levelName.substr(6).c_str(), 0, 0);
        if (levelID > 0)
        {
          // Do enough level adding so we can be sure to add the data to the right place
          for (int i = H.GetNumLevels(); i <= levelID; i++)
            H.AddNewLevel();
        }
        RCP<Level> level = H.GetLevel(levelID);

        RCP<FactoryManager> M = Teuchos::rcp_dynamic_cast<FactoryManager>(HM.GetFactoryManager(levelID));
        TEUCHOS_TEST_FOR_EXCEPTION(M.is_null(), Exceptions::InvalidArgument, "MueLu::Utils::AddNonSerializableDataToHierarchy: cannot get FactoryManager");

        // Grab the level sublist & loop over parameters
        const ParameterList& levelList = paramList.sublist(levelName);
        for (ParameterList::ConstIterator it2 = levelList.begin(); it2 != levelList.end(); it2++) {
          const std::string& name = it2->first;
          TEUCHOS_TEST_FOR_EXCEPTION(name != "A" && name != "P" && name != "R" && name != "K"  && name != "M" &&
                                     name != "Nullspace" && name != "Coordinates" && name != "pcoarsen: element to node map" &&
                                     !IsParamMuemexVariable(name), Exceptions::InvalidArgument,
                                     "MueLu::Utils::AddNonSerializableDataToHierarchy: parameter list contains unknown data type");

          if (name == "A") {
            level->Set(name, Teuchos::getValue<RCP<Matrix > > (it2->second),NoFactory::get());
            M->SetFactory(name, NoFactory::getRCP()); // TAW: not sure about this: be aware that this affects all levels
                                                      //      However, A is accessible through NoFactory anyway, so it should
                                                      //      be fine here.
          }
          else if(name == "P" || name == "R" || name == "K" || name == "M") {
            level->AddKeepFlag(name,NoFactory::get(),MueLu::UserData);
            level->Set(name, Teuchos::getValue<RCP<Matrix > >     (it2->second), NoFactory::get());
          }
          else if (name == "Nullspace")
          {
            level->AddKeepFlag(name,NoFactory::get(),MueLu::UserData);
            level->Set(name, Teuchos::getValue<RCP<MultiVector > >(it2->second), NoFactory::get());
            //M->SetFactory(name, NoFactory::getRCP()); // TAW: generally it is a bad idea to overwrite the factory manager data here
                                                        // One should do this only in very special cases
          }
          else if(name == "Coordinates") //Scalar of Coordinates MV is always double
          {
            level->AddKeepFlag(name,NoFactory::get(),MueLu::UserData);
            level->Set(name, Teuchos::getValue<RCP<RealValuedMultiVector> >(it2->second), NoFactory::get());
            //M->SetFactory(name, NoFactory::getRCP()); // TAW: generally it is a bad idea to overwrite the factory manager data here
          }
#ifdef HAVE_MUELU_INTREPID2
          else if (name == "pcoarsen: element to node map")
          {
            level->AddKeepFlag(name,NoFactory::get(),MueLu::UserData);
            level->Set(name, Teuchos::getValue<RCP<Kokkos::DynRankView<LocalOrdinal,typename Node::device_type> > >(it2->second), NoFactory::get());
          }
#endif
          else
#ifdef HAVE_MUELU_MATLAB
          {
            //Custom variable for Muemex
            size_t typeNameStart = name.find_first_not_of(' ');
            size_t typeNameEnd = name.find(' ', typeNameStart);
            std::string typeName = name.substr(typeNameStart, typeNameEnd - typeNameStart);
            std::transform(typeName.begin(), typeName.end(), typeName.begin(), ::tolower);
            level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
            if(typeName == "matrix")
              level->Set(name, Teuchos::getValue<RCP<Matrix> >(it2->second), NoFactory::get());
            else if(typeName == "multivector")
              level->Set(name, Teuchos::getValue<RCP<MultiVector> >(it2->second), NoFactory::get());
            else if(typeName == "map")
              level->Set(name, Teuchos::getValue<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >(it2->second), NoFactory::get());
            else if(typeName == "ordinalvector")
              level->Set(name, Teuchos::getValue<RCP<Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > >(it2->second), NoFactory::get());
            else if(typeName == "scalar")
              level->Set(name, Teuchos::getValue<Scalar>(it2->second), NoFactory::get());
            else if(typeName == "double")
              level->Set(name, Teuchos::getValue<double>(it2->second), NoFactory::get());
            else if(typeName == "complex")
              level->Set(name, Teuchos::getValue<std::complex<double> >(it2->second), NoFactory::get());
            else if(typeName == "int")
              level->Set(name, Teuchos::getValue<int>(it2->second), NoFactory::get());
            else if(typeName == "string")
              level->Set(name, Teuchos::getValue<std::string>(it2->second), NoFactory::get());
          }
#else
          {
            throw std::runtime_error("Invalid non-serializable data on list");
          }
#endif
        }
      } else if (paramList.isSublist(levelName) && levelName.find("user data") != std::string::npos) {
        // So far only put data on level 0
        int levelID = 0;
        RCP<Level> level = H.GetLevel(levelID);

        RCP<FactoryManager> M = Teuchos::rcp_dynamic_cast<FactoryManager>(HM.GetFactoryManager(levelID));
        TEUCHOS_TEST_FOR_EXCEPTION(M.is_null(), Exceptions::InvalidArgument, "MueLu::Utils::AddNonSerializableDataToHierarchy: cannot get FactoryManager");

        // Grab the user data sublist & loop over parameters
        const ParameterList& userList = paramList.sublist(levelName);
        for (ParameterList::ConstIterator it2 = userList.begin(); it2 != userList.end(); it2++) {
          const std::string& name = it2->first;
          TEUCHOS_TEST_FOR_EXCEPTION(name != "P" && name != "R" &&
                                     name != "Nullspace" && name != "Coordinates" &&
                                     !IsParamValidVariable(name), Exceptions::InvalidArgument,
                                     "MueLu::Utils::AddNonSerializableDataToHierarchy: user data parameter list contains unknown data type!");
          if( name == "P" || name == "R") {
            level->AddKeepFlag(name,NoFactory::get(),MueLu::UserData);
            level->Set(name, Teuchos::getValue<RCP<Matrix > >     (it2->second), NoFactory::get());
          } else if (name == "Nullspace") {
            level->AddKeepFlag(name,NoFactory::get(),MueLu::UserData);
            level->Set(name, Teuchos::getValue<RCP<MultiVector > >(it2->second), NoFactory::get());
            //M->SetFactory(name, NoFactory::getRCP()); // TAW: generally it is a bad idea to overwrite the factory manager data here
                                                        // One should do this only in very special cases
          } else if(name == "Coordinates") {//Scalar of Coordinates MV is always double
            level->AddKeepFlag(name,NoFactory::get(),MueLu::UserData);
            level->Set(name, Teuchos::getValue<RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > >(it2->second), NoFactory::get());
            //M->SetFactory(name, NoFactory::getRCP()); // TAW: generally it is a bad idea to overwrite the factory manager data here
          } else {
            //Custom variable
            size_t typeNameStart = name.find_first_not_of(' ');
            size_t typeNameEnd   = name.find(' ', typeNameStart);
            std::string typeName = name.substr(typeNameStart, typeNameEnd - typeNameStart);
            size_t varNameStart  = name.find_first_not_of(' ', typeNameEnd);
            std::string varName  = name.substr(varNameStart, name.size());
            std::transform(typeName.begin(), typeName.end(), typeName.begin(), ::tolower);
            level->AddKeepFlag(varName, NoFactory::get(), MueLu::UserData);
            if(typeName == "matrix")
              level->Set(varName, Teuchos::getValue<RCP<Matrix> >(it2->second), NoFactory::get());
            else if(typeName == "multivector")
              level->Set(varName, Teuchos::getValue<RCP<MultiVector> >(it2->second), NoFactory::get());
            else if(typeName == "map")
              level->Set(varName, Teuchos::getValue<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >(it2->second), NoFactory::get());
            else if(typeName == "ordinalvector")
              level->Set(varName, Teuchos::getValue<RCP<Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > >(it2->second), NoFactory::get());
            else if(typeName == "scalar")
              level->Set(varName, Teuchos::getValue<Scalar>(it2->second), NoFactory::get());
            else if(typeName == "double")
              level->Set(varName, Teuchos::getValue<double>(it2->second), NoFactory::get());
            else if(typeName == "complex")
              level->Set(varName, Teuchos::getValue<std::complex<double> >(it2->second), NoFactory::get());
            else if(typeName == "int")
              level->Set(varName, Teuchos::getValue<int>(it2->second), NoFactory::get());
            else if(typeName == "string")
              level->Set(varName, Teuchos::getValue<std::string>(it2->second), NoFactory::get());
            else if(typeName == "array<go>")
              level->Set(varName, Teuchos::getValue<Array<GlobalOrdinal> > (it2->second), NoFactory::get());
            else if(typeName == "array<lo>")
              level->Set(varName, Teuchos::getValue<Array<LocalOrdinal> >(it2->second), NoFactory::get());
            else if(typeName == "arrayrcp<lo>")
              level->Set(varName, Teuchos::getValue<ArrayRCP<LocalOrdinal> >(it2->second), NoFactory::get());
            else if(typeName == "arrayrcp<go>")
              level->Set(varName, Teuchos::getValue<ArrayRCP<GlobalOrdinal> >(it2->second), NoFactory::get());
            else
              throw std::runtime_error("Invalid non-serializable data on list");
          }
        }
      }
    }
  }
} // namespace MueLu

#define MUELU_HIERARCHY_UTILS_SHORT
#endif // MUELU_HIERARCHYHELPERS_DEF_HPP
