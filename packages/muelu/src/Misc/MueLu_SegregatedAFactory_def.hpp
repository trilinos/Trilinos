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
//                    Ray Tuminaro      (rstumin@sandia.gov)
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_SEGREGATEDAFACTORY_DEF_HPP
#define MUELU_SEGREGATEDAFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixFactory.hpp>

#include "MueLu_SegregatedAFactory_decl.hpp"

#include "MueLu_FactoryManager.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> SegregatedAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
#undef SET_VALID_ENTRY

  validParamList->set<RCP<const FactoryBase>> ("A", Teuchos::null, "Generating factory of the matrix A to be filtered");
  validParamList->set<std::string>            ("droppingScheme", "vague", "Strategy to drop entries from matrix A based on the input of some map(s) [blockmap, map-pair]");

  validParamList->set<RCP<const FactoryBase>> ("dropMap1", Teuchos::null, "Generating factory for dropMap1");
  validParamList->set<bool> ("Call ReduceAll on dropMap1", false, "Boolean for calling reduceAll on dropMap1");

  validParamList->set<RCP<const FactoryBase>> ("dropMap2", Teuchos::null, "Generating factory for dropMap2'");
  validParamList->set<bool> ("Call ReduceAll on dropMap2", false, "Boolean for calling reduceAll on dropMap2");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SegregatedAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {

  const ParameterList& pL = GetParameterList();

  TEUCHOS_TEST_FOR_EXCEPTION(pL.get<RCP<const FactoryBase>>("A")==Teuchos::null, Exceptions::InvalidArgument,
                             "Please specify a generating factory for the matrix \"A\" to be filtered.")
  TEUCHOS_TEST_FOR_EXCEPTION(pL.get<std::string>("droppingScheme")=="vague", Exceptions::InvalidArgument,
                             "Input map type not selected. Please select one of the available strategies.")
  TEUCHOS_TEST_FOR_EXCEPTION(
          (pL.get<std::string>("droppingScheme") != "blockmap" && pL.get<std::string>("droppingScheme") != "map-pair"),
          Exceptions::InvalidArgument,
          "Unknown User Input: map type (=" << pL.get<std::string>("droppingScheme") << ")")

  Input(currentLevel, "A");

  if (currentLevel.GetLevelID() == 0){
    // Not needed, as the map is provided as user data
    currentLevel.DeclareInput("dropMap1", NoFactory::get(), this);
    currentLevel.DeclareInput("dropMap2", NoFactory::get(), this);
  }
  else {
    // check whether user has provided a specific name for "map: factory"
    Input(currentLevel, "dropMap1");
    Input(currentLevel, "dropMap2");
  }

} // DeclareInput

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SegregatedAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const
{
  // Call a specialized build routine based on the format of user-given input
  const ParameterList &pL = GetParameterList();
  const std::string parameterName = "droppingScheme";
  if (pL.get<std::string>(parameterName) == "blockmap"){
    BuildBasedOnBlockmap(currentLevel);
  }
  else if (pL.get<std::string>(parameterName) == "map-pair") {
    BuildBasedOnMapPair(currentLevel);
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::InvalidArgument,
                               "MueLu::SegregatedAFactory::Build(): Unknown map type of user input. "
                               "Set a valid value for the parameter \"" << parameterName << "\".")
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SegregatedAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildBasedOnBlockmap(Level &currentLevel) const{
  FactoryMonitor m(*this, "Matrix filtering (segregation, blockmap)", currentLevel);

  RCP<Matrix> Ain = Get< RCP<Matrix> >(currentLevel, "A");
  RCP<const Map> dropMap1         = Teuchos::null;
  const std::string dropMap1Name  = "dropMap1";

  // fetch maps from level
  if (currentLevel.GetLevelID() == 0) {
    dropMap1 = currentLevel.Get<RCP<const Map>>(dropMap1Name, NoFactory::get());
    GetOStream(Statistics0) << "User provided dropMap1 \"" << dropMap1Name << "\": length dimension=" << dropMap1->getGlobalNumElements() << std::endl;
  }
  else {
    dropMap1 = Get<RCP<const Map>>(currentLevel, dropMap1Name);
  }
  TEUCHOS_ASSERT(!dropMap1.is_null());

  // create new empty Operator
  Teuchos::RCP<Matrix> Aout = MatrixFactory::Build(Ain->getRowMap(), Ain->getGlobalMaxNumRowEntries());

  size_t numLocalRows = Ain->getLocalNumRows();
  for(size_t row=0; row<numLocalRows; row++) {
    GlobalOrdinal grid = Ain->getRowMap()->getGlobalElement(row); // global row id
    bool isInMap = dropMap1->isNodeGlobalElement(grid);

    // extract row information from input matrix
    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> vals;
    Ain->getLocalRowView(row, indices, vals);

    // just copy all values in output
    Teuchos::ArrayRCP<GlobalOrdinal> indout(indices.size(), Teuchos::ScalarTraits<GlobalOrdinal>::zero());
    Teuchos::ArrayRCP<Scalar> valout(indices.size(), Teuchos::ScalarTraits<Scalar>::zero());

    size_t nNonzeros = 0;
    for(size_t i=0; i<(size_t)indices.size(); i++) {
      GlobalOrdinal gcid = Ain->getColMap()->getGlobalElement(indices[i]); // global column id
      bool isInMap2 = dropMap1->isNodeGlobalElement(gcid);

      if (isInMap == isInMap2) {
        indout [nNonzeros] = gcid;
        valout [nNonzeros] = vals[i];
        nNonzeros++;
      }
    }
    indout.resize(nNonzeros);
    valout.resize(nNonzeros);

    Aout->insertGlobalValues(Ain->getRowMap()->getGlobalElement(row), indout.view(0,indout.size()), valout.view(0,valout.size()));
  }

  Aout->fillComplete(Ain->getDomainMap(), Ain->getRangeMap());

  // copy block size information
  Aout->SetFixedBlockSize(Ain->GetFixedBlockSize());

  GetOStream(Statistics0, 0) << "Nonzeros in A (input): " << Ain->getGlobalNumEntries() << ", Nonzeros after filtering A: " << Aout->getGlobalNumEntries() << std::endl;

  Set(currentLevel, "A", Aout);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SegregatedAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildBasedOnMapPair(Level &currentLevel) const
{
  FactoryMonitor m(*this, "Matrix filtering (segregation, map-pair)", currentLevel);

  RCP<Matrix> Ain = Get< RCP<Matrix> >(currentLevel, "A");

  const ParameterList &pL = GetParameterList();
  auto comm = Ain->getRowMap()->getComm();

  // fetch maps from level
  RCP<const Map> dropMap1 = Teuchos::null;
  RCP<const Map> dropMap2 = Teuchos::null;

  const std::string dropMap1Name = "dropMap1";
  const std::string dropMap2Name = "dropMap2";

  if (currentLevel.GetLevelID() == 0) {
    dropMap1 = currentLevel.Get<RCP<const Map>>(dropMap1Name, NoFactory::get());
    dropMap2 = currentLevel.Get<RCP<const Map>>(dropMap2Name, NoFactory::get());
    GetOStream(Statistics0) << "User provided dropMap1 \"" << dropMap1Name << "\": length dimension=" << dropMap1->getGlobalNumElements() << std::endl;
    GetOStream(Statistics0) << "User provided dropMap2 \"" << dropMap2Name << "\": length dimension=" << dropMap2->getGlobalNumElements() << std::endl;
  }
  else {
    dropMap1 = Get<RCP<const Map>>(currentLevel, dropMap1Name);
    dropMap2 = Get<RCP<const Map>>(currentLevel, dropMap2Name);
  }

  TEUCHOS_ASSERT(!dropMap1.is_null());
  TEUCHOS_ASSERT(!dropMap2.is_null());

  // create new empty Operator
  Teuchos::RCP<Matrix> Aout = MatrixFactory::Build(Ain->getRowMap(), Ain->getGlobalMaxNumRowEntries());

  Teuchos::RCP<const Map> finalDropMap1 = Teuchos::null;
  Teuchos::RCP<const Map> finalDropMap2 = Teuchos::null;

  if (pL.get<bool>("Call ReduceAll on dropMap1")){
    finalDropMap1 = CreateRedundantMaps(dropMap1);
  }
  else{
    // if reduceAll is not called, we simply work with the local dropMap
    finalDropMap1 = dropMap1;
  }

  if (pL.get<bool>("Call ReduceAll on dropMap2")){
    finalDropMap2 = CreateRedundantMaps(dropMap2);
  }
  else{
    // if reduceAll is not called, we simply work with the local dropMap
    finalDropMap2 = dropMap2;
  }

  // Start copying the matrix row by row and dropping any entries that are contained as a combination of entries of
  // dropMap1 and dropMap2
  size_t numLocalMatrixRows = Ain->getLocalNumRows();

  for(size_t row=0; row<numLocalMatrixRows; row++) {
    GlobalOrdinal grid = Ain->getRowMap()->getGlobalElement(row); // global row id
    bool rowIsInMap1 = finalDropMap1->isNodeGlobalElement(grid);
    bool rowIsInMap2 = finalDropMap2->isNodeGlobalElement(grid);

    // extract row information from input matrix
    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> vals;
    Ain->getLocalRowView(row, indices, vals);

    // just copy all values in output
    Teuchos::ArrayRCP<GlobalOrdinal> indout(indices.size(),Teuchos::ScalarTraits<GlobalOrdinal>::zero());
    Teuchos::ArrayRCP<Scalar>        valout(indices.size(),Teuchos::ScalarTraits<Scalar>::zero());

    size_t nNonzeros = 0;
    for(size_t i=0; i<(size_t)indices.size(); i++)
    {
      GlobalOrdinal gcid = Ain->getColMap()->getGlobalElement(indices[i]); // global column id
      bool colIsInMap1 = finalDropMap1->isNodeGlobalElement(gcid);
      bool colIsInMap2 = finalDropMap2->isNodeGlobalElement(gcid);

      if ((rowIsInMap1 && colIsInMap2) || (rowIsInMap2 && colIsInMap1)){
        // do nothing == drop this entry
      }
      else {
        indout [nNonzeros] = gcid;
        valout [nNonzeros] = vals[i];
        nNonzeros++;
      }
    }
    indout.resize(nNonzeros);
    valout.resize(nNonzeros);
    Aout->insertGlobalValues(Ain->getRowMap()->getGlobalElement(row), indout.view(0,indout.size()), valout.view(0,valout.size()));
  }

  Aout->fillComplete(Ain->getDomainMap(), Ain->getRangeMap());
  // copy block size information
  Aout->SetFixedBlockSize(Ain->GetFixedBlockSize());
  Aout->fillComplete(Ain->getDomainMap(), Ain->getRangeMap());

  // copy block size information
  Aout->SetFixedBlockSize(Ain->GetFixedBlockSize());

  GetOStream(Statistics0, 0) << "Nonzeros in A (input): " << Ain->getGlobalNumEntries() << ", Nonzeros after filtering A: " << Aout->getGlobalNumEntries() << std::endl;

  Set(currentLevel, "A", Aout);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
        SegregatedAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateRedundantMaps(
                Teuchos::RCP<const Map> localDropMap) const {

  Teuchos::RCP<const Teuchos::Comm<int>> comm = localDropMap->getComm();

  Array<GO> localDropMapGIDList = localDropMap->getLocalElementList();
  const int GIDListSize = localDropMap->getMaxAllGlobalIndex()+1;

//  Create a list of GID with only an incomplete/partial set of GIDs, which can then be completed by reduceAll
  Array<GO> partialDropMapGIDList(GIDListSize, -Teuchos::ScalarTraits<GlobalOrdinal>::one());
  Array<GO> redundantDropMapGIDList(GIDListSize, -Teuchos::ScalarTraits<GlobalOrdinal>::one());

  for(GO gid : localDropMapGIDList){
    partialDropMapGIDList[gid] = gid;
  }

  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, GIDListSize, &partialDropMapGIDList[0],
                     &redundantDropMapGIDList[0]);
  redundantDropMapGIDList.erase(std::remove(redundantDropMapGIDList.begin(), redundantDropMapGIDList.end(), -1),
                                 redundantDropMapGIDList.end());
  Teuchos::RCP<const Map> redundantDropMap = MapFactory::Build(
          localDropMap->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), redundantDropMapGIDList, 0, comm);

  return redundantDropMap;
}

} //namespace MueLu

#endif  // MUELU_SEGREGATEDAFACTORY_DEF_HPP
