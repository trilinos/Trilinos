// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_HIERARCHYUTILS_DEF_HPP
#define MUELU_HIERARCHYUTILS_DEF_HPP

#include "Teuchos_ScalarTraits.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_Operator.hpp>

#include "MueLu_HierarchyUtils_decl.hpp"
#include "MueLu_HierarchyManager.hpp"
#include "MueLu_FactoryManager.hpp"

// TODO/FIXME: DeclareInput(, **this**) cannot be used here
#ifdef HAVE_MUELU_INTREPID2
#include "Kokkos_DynRankView.hpp"
#endif

namespace MueLu {

// Copy object from one hierarchy to another calling AddNewLevel as appropriate.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void HierarchyUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CopyBetweenHierarchies(Hierarchy& fromHierarchy, Hierarchy& toHierarchy, const std::string fromLabel, const std::string toLabel, const std::string dataType) {
  // add any necessary levels
  for (int i = toHierarchy.GetNumLevels(); i < fromHierarchy.GetNumLevels(); i++)
    toHierarchy.AddNewLevel();

  for (int i = 0; i < fromHierarchy.GetNumLevels(); i++) {
    RCP<Level> fromLevel = fromHierarchy.GetLevel(i);
    RCP<Level> toLevel   = toHierarchy.GetLevel(i);

    TEUCHOS_TEST_FOR_EXCEPTION(dataType != "RCP<Matrix>" && dataType != "RCP<const Import>", Exceptions::InvalidArgument,
                               std::string("MueLu::Utils::CopyBetweenHierarchies: unknown data type(") + dataType + ")");
    if (fromLevel->IsAvailable(fromLabel)) {
      if (dataType == "RCP<Matrix>") {
        // Normally, we should only do
        //      toLevel->Set(toLabel,fromLevel->Get<RCP<Matrix> >(fromLabel));
        // The logic below is meant to handle a special case when we
        // repartition a processor away, leaving behind a RCP<Operator> on
        // on the level instead of an RCP<Matrix>

        auto tempOp     = fromLevel->Get<RCP<Operator>>(fromLabel);
        auto tempMatrix = rcp_dynamic_cast<Matrix>(tempOp);
        if (!tempMatrix.is_null())
          toLevel->Set(toLabel, tempMatrix);
        else
          toLevel->Set(toLabel, tempOp);
      }
      if (dataType == "RCP<const Import>") {
        toLevel->Set(toLabel, fromLevel->Get<RCP<const Import>>(fromLabel));
      }
    }
  }
}

// Adds the following non-serializable data (A,P,R,Nullspace,Coordinates) from level-specific sublist nonSerialList,
// calling AddNewLevel as appropriate.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void HierarchyUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddNonSerializableDataToHierarchy(HierarchyManager& HM, Hierarchy& H, const ParameterList& nonSerialList) {
  typedef typename Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType,
                                       LocalOrdinal, GlobalOrdinal, Node>
      realvaluedmultivector_type;

  for (ParameterList::ConstIterator nonSerialEntry = nonSerialList.begin(); nonSerialEntry != nonSerialList.end(); nonSerialEntry++) {
    const std::string& levelName = nonSerialEntry->first;
    // Check for match of the form "level X" where X is a positive integer
    if (nonSerialList.isSublist(levelName) && levelName.find("level ") == 0 && levelName.size() > 6) {
      int levelID = strtol(levelName.substr(6).c_str(), 0, 0);
      if (levelID > 0) {
        // Do enough level adding so we can be sure to add the data to the right place
        for (int i = H.GetNumLevels(); i <= levelID; i++)
          H.AddNewLevel();
      }
      RCP<Level> level = H.GetLevel(levelID);

      RCP<FactoryManager> M = Teuchos::rcp_dynamic_cast<FactoryManager>(HM.GetFactoryManager(levelID));
      TEUCHOS_TEST_FOR_EXCEPTION(M.is_null(), Exceptions::InvalidArgument, "MueLu::Utils::AddNonSerializableDataToHierarchy: cannot get FactoryManager");

      // Grab the level sublist & loop over parameters
      const ParameterList& levelList = nonSerialList.sublist(levelName);
      for (ParameterList::ConstIterator levelListEntry = levelList.begin(); levelListEntry != levelList.end(); levelListEntry++) {
        const std::string& name = levelListEntry->first;
        TEUCHOS_TEST_FOR_EXCEPTION(name != "A" && name != "P" && name != "R" && name != "K" && name != "M" && name != "Mdiag" &&
                                       name != "D0" && name != "Dk_1" && name != "Dk_2" &&
                                       name != "Mk_one" && name != "Mk_1_one" && name != "M1_beta" && name != "M1_alpha" &&
                                       name != "invMk_1_invBeta" && name != "invMk_2_invAlpha" &&
                                       name != "M1" && name != "Ms" && name != "M0inv" &&
                                       name != "Pnodal" && name != "NodeMatrix" && name != "NodeAggMatrix" &&
                                       name != "Nullspace" && name != "Coordinates" && name != "pcoarsen: element to node map" &&
                                       name != "Node Comm" && name != "DualNodeID2PrimalNodeID" && name != "Primal interface DOF map" &&
                                       name != "dropMap1" && name != "dropMap2" &&
                                       !IsParamMuemexVariable(name),
                                   Exceptions::InvalidArgument,
                                   std::string("MueLu::Utils::AddNonSerializableDataToHierarchy: parameter list contains unknown data type(") + name + ")");

        // Get a valid communicator and lib
        RCP<const Teuchos::Comm<int>> comm;
        if (!level->GetComm().is_null())
          comm = level->GetComm();
        else if (level->IsAvailable("A")) {
          RCP<Matrix> mat;
          level->Get("A", mat);
          comm = mat->getMap()->getComm();
        } else {
          RCP<Level> level0 = H.GetLevel(0);
          if (!level0->GetComm().is_null())
            comm = level0->GetComm();
          else {
            RCP<Matrix> mat;
            level0->Get("A", mat);
            comm = mat->getMap()->getComm();
          }
        }
        Xpetra::UnderlyingLib lib = level->lib();

        if (name == "A") {
          RCP<Matrix> mat;
          if (levelListEntry->second.isType<std::string>())
            // We might also want to read maps here.
            mat = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read(Teuchos::getValue<std::string>(levelListEntry->second), lib, comm);
          else
            mat = Teuchos::getValue<RCP<Matrix>>(levelListEntry->second);
          level->Set(name, mat, NoFactory::get());
          M->SetFactory(name, NoFactory::getRCP());  // TAW: not sure about this: be aware that this affects all levels
                                                     //      However, A is accessible through NoFactory anyway, so it should
                                                     //      be fine here.
        } else if (name == "P" || name == "R" || name == "K" || name == "M") {
          if (levelListEntry->second.isType<RCP<Operator>>()) {
            RCP<Operator> mat;
            mat = Teuchos::getValue<RCP<Operator>>(levelListEntry->second);

            RCP<const FactoryBase> fact = M->GetFactory(name);
            level->AddKeepFlag(name, fact.get(), MueLu::UserData);
            level->Set(name, mat, fact.get());

            level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
            level->Set(name, mat, NoFactory::get());
          } else {
            RCP<Matrix> mat;
            if (levelListEntry->second.isType<std::string>())
              // We might also want to read maps here.
              mat = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read(Teuchos::getValue<std::string>(levelListEntry->second), lib, comm);
            else
              mat = Teuchos::getValue<RCP<Matrix>>(levelListEntry->second);

            RCP<const FactoryBase> fact = M->GetFactory(name);
            level->AddKeepFlag(name, fact.get(), MueLu::UserData);
            level->Set(name, mat, fact.get());

            level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
            level->Set(name, mat, NoFactory::get());
          }
        } else if (name == "D0" || name == "Dk_1" || name == "Dk_2" ||
                   name == "Mk_one" || name == "Mk_1_one" || name == "M1_beta" || name == "M1_alpha" ||
                   name == "invMk_1_invBeta" || name == "invMk_2_invAlpha" ||
                   name == "M1" || name == "Ms" || name == "M0inv" ||
                   name == "Pnodal" || name == "NodeMatrix" || name == "NodeAggMatrix") {
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          if (levelListEntry->second.isType<RCP<Operator>>())
            level->Set(name, Teuchos::getValue<RCP<Operator>>(levelListEntry->second), NoFactory::get());
          else
            level->Set(name, Teuchos::getValue<RCP<Matrix>>(levelListEntry->second), NoFactory::get());
        } else if (name == "Mdiag") {
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          level->Set(name, Teuchos::getValue<RCP<Vector>>(levelListEntry->second), NoFactory::get());
        } else if (name == "Nullspace") {
          RCP<MultiVector> vec;
          if (levelListEntry->second.isType<std::string>()) {
            TEUCHOS_ASSERT(level->IsAvailable("A"));
            RCP<Matrix> mat;
            level->Get("A", mat);
            auto map = mat->getMap();
            vec      = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ReadMultiVector(Teuchos::getValue<std::string>(levelListEntry->second), map);
          } else
            vec = Teuchos::getValue<RCP<MultiVector>>(levelListEntry->second);
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          level->Set(name, vec, NoFactory::get());
          // M->SetFactory(name, NoFactory::getRCP()); // TAW: generally it is a bad idea to overwrite the factory manager data here
          //  One should do this only in very special cases
        } else if (name == "Coordinates")  // Scalar of Coordinates MV is always double
        {
          RCP<realvaluedmultivector_type> vec;
          if (levelListEntry->second.isType<std::string>()) {
            TEUCHOS_ASSERT(level->IsAvailable("A"));
            RCP<Matrix> mat;
            level->Get("A", mat);
            size_t blkSize         = mat->GetFixedBlockSize();
            RCP<const Map> nodeMap = mat->getRowMap();
            if (blkSize > 1) {
              // Create a nodal map, as coordinates have not been expanded to a DOF map yet.
              RCP<const Map> dofMap = mat->getRowMap();
              GO indexBase          = dofMap->getIndexBase();
              size_t numLocalDOFs   = dofMap->getLocalNumElements();
              TEUCHOS_TEST_FOR_EXCEPTION(numLocalDOFs % blkSize, Exceptions::RuntimeError,
                                         "HierarchyUtils: block size (" << blkSize << ") is incompatible with the number of local dofs in a row map (" << numLocalDOFs);
              ArrayView<const GO> GIDs = dofMap->getLocalElementList();

              Array<GO> nodeGIDs(numLocalDOFs / blkSize);
              for (size_t i = 0; i < numLocalDOFs; i += blkSize)
                nodeGIDs[i / blkSize] = (GIDs[i] - indexBase) / blkSize + indexBase;

              Xpetra::global_size_t INVALID = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
              nodeMap                       = MapFactory::Build(dofMap->lib(), INVALID, nodeGIDs(), indexBase, dofMap->getComm());
            }
            vec = Xpetra::IO<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node>::ReadMultiVector(Teuchos::getValue<std::string>(levelListEntry->second), nodeMap);
          } else
            vec = Teuchos::getValue<RCP<realvaluedmultivector_type>>(levelListEntry->second);
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          level->Set(name, vec, NoFactory::get());
          // M->SetFactory(name, NoFactory::getRCP()); // TAW: generally it is a bad idea to overwrite the factory manager data here
        } else if (name == "Node Comm") {
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          level->Set(name, Teuchos::getValue<RCP<const Teuchos::Comm<int>>>(levelListEntry->second), NoFactory::get());
        } else if (name == "DualNodeID2PrimalNodeID") {
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          level->Set(name, Teuchos::getValue<RCP<std::map<LO, LO>>>(levelListEntry->second), NoFactory::get());
        } else if (name == "Primal interface DOF map") {
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          level->Set(name, Teuchos::getValue<RCP<const Map>>(levelListEntry->second), NoFactory::get());
        } else if (name == "dropMap1") {
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          level->Set(name, Teuchos::getValue<RCP<const Map>>(levelListEntry->second), NoFactory::get());
        } else if (name == "dropMap2") {
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          level->Set(name, Teuchos::getValue<RCP<const Map>>(levelListEntry->second), NoFactory::get());
        }
#ifdef HAVE_MUELU_INTREPID2
        else if (name == "pcoarsen: element to node map") {
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          level->Set(name, Teuchos::getValue<RCP<Kokkos::DynRankView<LocalOrdinal, typename Node::device_type>>>(levelListEntry->second), NoFactory::get());
        }
#endif
        else
#ifdef HAVE_MUELU_MATLAB
        {
          // Custom variable for Muemex
          size_t typeNameStart = name.find_first_not_of(' ');
          size_t typeNameEnd   = name.find(' ', typeNameStart);
          std::string typeName = name.substr(typeNameStart, typeNameEnd - typeNameStart);
          std::transform(typeName.begin(), typeName.end(), typeName.begin(), ::tolower);
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          if (typeName == "matrix")
            level->Set(name, Teuchos::getValue<RCP<Matrix>>(levelListEntry->second), NoFactory::get());
          else if (typeName == "multivector")
            level->Set(name, Teuchos::getValue<RCP<MultiVector>>(levelListEntry->second), NoFactory::get());
          else if (typeName == "map")
            level->Set(name, Teuchos::getValue<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>>(levelListEntry->second), NoFactory::get());
          else if (typeName == "ordinalvector")
            level->Set(name, Teuchos::getValue<RCP<Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>>(levelListEntry->second), NoFactory::get());
          else if (typeName == "scalar")
            level->Set(name, Teuchos::getValue<Scalar>(levelListEntry->second), NoFactory::get());
          else if (typeName == "double")
            level->Set(name, Teuchos::getValue<double>(levelListEntry->second), NoFactory::get());
          else if (typeName == "complex")
            level->Set(name, Teuchos::getValue<std::complex<double>>(levelListEntry->second), NoFactory::get());
          else if (typeName == "int")
            level->Set(name, Teuchos::getValue<int>(levelListEntry->second), NoFactory::get());
          else if (typeName == "string")
            level->Set(name, Teuchos::getValue<std::string>(levelListEntry->second), NoFactory::get());
        }
#else
        {
          throw std::runtime_error("Invalid non-serializable data on list");
        }
#endif
      }
    } else if (nonSerialList.isSublist(levelName) && levelName.find("user data") != std::string::npos) {
      // So far only put data on level 0
      int levelID      = 0;
      RCP<Level> level = H.GetLevel(levelID);

      RCP<FactoryManager> M = Teuchos::rcp_dynamic_cast<FactoryManager>(HM.GetFactoryManager(levelID));
      TEUCHOS_TEST_FOR_EXCEPTION(M.is_null(), Exceptions::InvalidArgument, "MueLu::Utils::AddNonSerializableDataToHierarchy: cannot get FactoryManager");

      // Grab the user data sublist & loop over parameters
      const ParameterList& userList = nonSerialList.sublist(levelName);
      for (ParameterList::ConstIterator userListEntry = userList.begin(); userListEntry != userList.end(); userListEntry++) {
        const std::string& name = userListEntry->first;
        TEUCHOS_TEST_FOR_EXCEPTION(name != "P" && name != "R" && name != "K" && name != "M" && name != "Mdiag" &&
                                       name != "D0" && name != "Dk_1" && name != "Dk_2" &&
                                       name != "Mk_one" && name != "Mk_1_one" && name != "M1_beta" && name != "M1_alpha" &&
                                       name != "invMk_1_invBeta" && name != "invMk_2_invAlpha" &&
                                       name != "M1" && name != "Ms" && name != "M0inv" &&
                                       name != "NodeMatrix" &&
                                       name != "Nullspace" && name != "Coordinates" && name != "pcoarsen: element to node map" &&
                                       name != "Node Comm" && name != "DualNodeID2PrimalNodeID" && name != "Primal interface DOF map" &&
                                       name != "dropMap1" && name != "dropMap2" &&
                                       name != "output stream" &&
                                       !IsParamValidVariable(name),
                                   Exceptions::InvalidArgument,
                                   std::string("MueLu::Utils::AddNonSerializableDataToHierarchy: user data parameter list contains unknown data type (") + name + ")");
        if (name == "P" || name == "R" || name == "K" || name == "M" ||
            name == "D0" || name == "Dk_1" || name == "Dk_2" ||
            name == "Mk_one" || name == "Mk_1_one" || name == "M1_beta" || name == "M1_alpha" ||
            name == "invMk_1_invBeta" || name == "invMk_2_invAlpha" ||
            name == "M1" || name == "Ms" || name == "M0inv" ||
            name == "NodeMatrix") {
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          level->Set(name, Teuchos::getValue<RCP<Matrix>>(userListEntry->second), NoFactory::get());
        } else if (name == "Mdiag") {
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          level->Set(name, Teuchos::getValue<RCP<Vector>>(userListEntry->second), NoFactory::get());
        } else if (name == "Nullspace") {
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          level->Set(name, Teuchos::getValue<RCP<MultiVector>>(userListEntry->second), NoFactory::get());
          // M->SetFactory(name, NoFactory::getRCP()); // TAW: generally it is a bad idea to overwrite the factory manager data here
          //  One should do this only in very special cases
        } else if (name == "Coordinates") {  // Scalar of Coordinates MV is always double
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          level->Set(name, Teuchos::getValue<RCP<realvaluedmultivector_type>>(userListEntry->second), NoFactory::get());
        } else if (name == "Node Comm") {
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          level->Set(name, Teuchos::getValue<RCP<const Teuchos::Comm<int>>>(userListEntry->second), NoFactory::get());
        } else if (name == "DualNodeID2PrimalNodeID") {
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          level->Set(name, Teuchos::getValue<RCP<std::map<LO, LO>>>(userListEntry->second), NoFactory::get());
        } else if (name == "Primal interface DOF map") {
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          level->Set(name, Teuchos::getValue<RCP<const Map>>(userListEntry->second), NoFactory::get());
        } else if (name == "dropMap1") {
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          level->Set(name, Teuchos::getValue<RCP<const Map>>(userListEntry->second), NoFactory::get());
        } else if (name == "dropMap2") {
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          level->Set(name, Teuchos::getValue<RCP<const Map>>(userListEntry->second), NoFactory::get());
        }
#ifdef HAVE_MUELU_INTREPID2
        else if (name == "pcoarsen: element to node map") {
          level->AddKeepFlag(name, NoFactory::get(), MueLu::UserData);
          level->Set(name, Teuchos::getValue<RCP<Kokkos::DynRankView<LocalOrdinal, typename Node::device_type>>>(userListEntry->second), NoFactory::get());
        }
#endif
        else if (name == "output stream") {
          H.SetMueLuOStream(Teuchos::getValue<RCP<Teuchos::FancyOStream>>(userListEntry->second));
        } else {
          // Custom variable
          size_t typeNameStart = name.find_first_not_of(' ');
          size_t typeNameEnd   = name.find(' ', typeNameStart);
          std::string typeName = name.substr(typeNameStart, typeNameEnd - typeNameStart);
          size_t varNameStart  = name.find_first_not_of(' ', typeNameEnd);
          std::string varName  = name.substr(varNameStart, name.size());
          std::transform(typeName.begin(), typeName.end(), typeName.begin(), ::tolower);
          level->AddKeepFlag(varName, NoFactory::get(), MueLu::UserData);
          if (typeName == "matrix")
            level->Set(varName, Teuchos::getValue<RCP<Matrix>>(userListEntry->second), NoFactory::get());
          else if (typeName == "multivector")
            level->Set(varName, Teuchos::getValue<RCP<MultiVector>>(userListEntry->second), NoFactory::get());
          else if (typeName == "vector")
            level->Set(varName, Teuchos::getValue<RCP<Vector>>(userListEntry->second), NoFactory::get());
          else if (typeName == "map")
            level->Set(varName, Teuchos::getValue<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>>(userListEntry->second), NoFactory::get());
          else if (typeName == "ordinalvector")
            level->Set(varName, Teuchos::getValue<RCP<Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>>(userListEntry->second), NoFactory::get());
          else if (typeName == "scalar")
            level->Set(varName, Teuchos::getValue<Scalar>(userListEntry->second), NoFactory::get());
          else if (typeName == "double")
            level->Set(varName, Teuchos::getValue<double>(userListEntry->second), NoFactory::get());
          else if (typeName == "complex")
            level->Set(varName, Teuchos::getValue<std::complex<double>>(userListEntry->second), NoFactory::get());
          else if (typeName == "int")
            level->Set(varName, Teuchos::getValue<int>(userListEntry->second), NoFactory::get());
          else if (typeName == "string")
            level->Set(varName, Teuchos::getValue<std::string>(userListEntry->second), NoFactory::get());
          else if (typeName == "array<go>")
            level->Set(varName, Teuchos::getValue<Array<GlobalOrdinal>>(userListEntry->second), NoFactory::get());
          else if (typeName == "array<lo>")
            level->Set(varName, Teuchos::getValue<Array<LocalOrdinal>>(userListEntry->second), NoFactory::get());
          else if (typeName == "arrayrcp<lo>")
            level->Set(varName, Teuchos::getValue<ArrayRCP<LocalOrdinal>>(userListEntry->second), NoFactory::get());
          else if (typeName == "arrayrcp<go>")
            level->Set(varName, Teuchos::getValue<ArrayRCP<GlobalOrdinal>>(userListEntry->second), NoFactory::get());
          else
            throw std::runtime_error("Invalid non-serializable data on list");
        }
      }
      // level->print(std::cout, MueLu::Debug);
    }
  }
}
}  // namespace MueLu

#define MUELU_HIERARCHY_UTILS_SHORT
#endif  // MUELU_HIERARCHYHELPERS_DEF_HPP
