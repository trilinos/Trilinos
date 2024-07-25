// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SEGREGATEDAFACTORY_DEF_HPP
#define MUELU_SEGREGATEDAFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixFactory.hpp>

#include "MueLu_SegregatedAFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities_decl.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> SegregatedAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase>>("A", Teuchos::null, "Generating factory of the matrix A to be filtered");
  validParamList->set<std::string>("droppingScheme", "vague", "Strategy to drop entries from matrix A based on the input of some map(s) [blockmap, map-pair]");

  validParamList->set<RCP<const FactoryBase>>("dropMap1", Teuchos::null, "Generating factory for dropMap1");  ////
  validParamList->set<RCP<const FactoryBase>>("dropMap2", Teuchos::null, "Generating factory for dropMap2'");

  return validParamList;
}  // GetValidParameterList

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SegregatedAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
  const ParameterList &pL = GetParameterList();

  TEUCHOS_TEST_FOR_EXCEPTION(pL.get<RCP<const FactoryBase>>("A") == Teuchos::null, Exceptions::InvalidArgument,
                             "Please specify a generating factory for the matrix \"A\" to be filtered.")
  TEUCHOS_TEST_FOR_EXCEPTION(pL.get<std::string>("droppingScheme") == "vague", Exceptions::InvalidArgument,
                             "Input map type not selected. Please select one of the available strategies.")
  TEUCHOS_TEST_FOR_EXCEPTION(
      (pL.get<std::string>("droppingScheme") != "blockmap" && pL.get<std::string>("droppingScheme") != "map-pair"),
      Exceptions::InvalidArgument,
      "Unknown User Input: droppingScheme (=" << pL.get<std::string>("droppingScheme") << ")")

  Input(currentLevel, "A");

  if (pL.get<std::string>("droppingScheme") == "blockmap") {
    if (currentLevel.GetLevelID() == 0) {
      currentLevel.DeclareInput("dropMap1", NoFactory::get(), this);
    } else {
      Input(currentLevel, "dropMap1");
    }
  } else if (pL.get<std::string>("droppingScheme") == "map-pair") {
    if (currentLevel.GetLevelID() == 0) {
      currentLevel.DeclareInput("dropMap1", NoFactory::get(), this);
      currentLevel.DeclareInput("dropMap2", NoFactory::get(), this);
    } else {
      Input(currentLevel, "dropMap1");
      Input(currentLevel, "dropMap2");
    }
  } else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::InvalidArgument, "Unknown droppingScheme.")

}  // DeclareInput

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SegregatedAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const {
  // Call a specialized build routine based on the format of user-given input
  const ParameterList &pL         = GetParameterList();
  const std::string parameterName = "droppingScheme";
  if (pL.get<std::string>(parameterName) == "blockmap") {
    BuildBasedOnBlockmap(currentLevel);
  } else if (pL.get<std::string>(parameterName) == "map-pair") {
    BuildBasedOnMapPair(currentLevel);
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::InvalidArgument,
                               "MueLu::SegregatedAFactory::Build(): Unknown map type of user input. "
                               "Set a valid value for the parameter \""
                                   << parameterName << "\".")
  }
}  // Build

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SegregatedAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildBasedOnBlockmap(Level &currentLevel) const {
  FactoryMonitor m(*this, "Matrix filtering (segregation, blockmap)", currentLevel);

  RCP<Matrix> Ain                = Get<RCP<Matrix>>(currentLevel, "A");
  RCP<const Map> dropMap1        = Teuchos::null;
  const std::string dropMap1Name = "dropMap1";

  // fetch maps from level
  if (currentLevel.GetLevelID() == 0) {
    dropMap1 = currentLevel.Get<RCP<const Map>>(dropMap1Name, NoFactory::get());
    GetOStream(Statistics0) << "User provided dropMap1 \"" << dropMap1Name << "\": length dimension=" << dropMap1->getGlobalNumElements() << std::endl;
  } else {
    dropMap1 = Get<RCP<const Map>>(currentLevel, dropMap1Name);
  }
  TEUCHOS_ASSERT(!dropMap1.is_null());

  // create new empty Operator
  Teuchos::RCP<Matrix> Aout = MatrixFactory::Build(Ain->getRowMap(), Ain->getGlobalMaxNumRowEntries());

  size_t numLocalRows = Ain->getLocalNumRows();
  for (size_t row = 0; row < numLocalRows; row++) {
    GlobalOrdinal grid = Ain->getRowMap()->getGlobalElement(row);  // global row id
    bool isInMap       = dropMap1->isNodeGlobalElement(grid);

    // extract row information from input matrix
    auto lclMat  = Ain->getLocalMatrixHost();
    auto rowView = lclMat.row(row);

    // just copy all values in output
    Teuchos::ArrayRCP<GO> indout(rowView.length, Teuchos::ScalarTraits<GO>::zero());
    Teuchos::ArrayRCP<SC> valout(rowView.length, Teuchos::ScalarTraits<SC>::zero());

    size_t nNonzeros = 0;
    for (LO jj = 0; jj < rowView.length; ++jj) {
      LO lcid       = rowView.colidx(jj);
      GO gcid       = Ain->getColMap()->getGlobalElement(lcid);
      auto val      = rowView.value(jj);
      bool isInMap2 = dropMap1->isNodeGlobalElement(gcid);

      if (isInMap == isInMap2) {
        indout[nNonzeros] = gcid;
        valout[nNonzeros] = val;
        nNonzeros++;
      }
    }
    indout.resize(nNonzeros);
    valout.resize(nNonzeros);

    Aout->insertGlobalValues(Ain->getRowMap()->getGlobalElement(row), indout.view(0, indout.size()), valout.view(0, valout.size()));
  }

  Aout->fillComplete(Ain->getDomainMap(), Ain->getRangeMap());

  // copy block size information
  Aout->SetFixedBlockSize(Ain->GetFixedBlockSize());

  GetOStream(Statistics0, 0) << "Nonzeros in A (input): " << Ain->getGlobalNumEntries() << ", Nonzeros after filtering A: " << Aout->getGlobalNumEntries() << std::endl;

  Set(currentLevel, "A", Aout);
}  // BuildBasedOnBlockmap

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SegregatedAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildBasedOnMapPair(Level &currentLevel) const {
  FactoryMonitor m(*this, "Matrix filtering (segregation, map-pair)", currentLevel);

  RCP<Matrix> Ain = Get<RCP<Matrix>>(currentLevel, "A");

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
  } else {
    dropMap1 = Get<RCP<const Map>>(currentLevel, dropMap1Name);
    dropMap2 = Get<RCP<const Map>>(currentLevel, dropMap2Name);
  }

  TEUCHOS_ASSERT(!dropMap1.is_null());
  TEUCHOS_ASSERT(!dropMap2.is_null());

  // create new empty Operator
  Teuchos::RCP<Matrix> Aout = MatrixFactory::Build(Ain->getRowMap(), Ain->getGlobalMaxNumRowEntries());

  // import the dropping information from other procs for off Rank entries
  Teuchos::RCP<const Map> finalDropMap1 = Teuchos::null;
  Teuchos::RCP<const Map> finalDropMap2 = Teuchos::null;

  finalDropMap1 = MueLu::importOffRankDroppingInfo(dropMap1, Ain);
  finalDropMap2 = MueLu::importOffRankDroppingInfo(dropMap2, Ain);

  // Start copying the matrix row by row and dropping any entries that are contained as a combination of entries of
  // dropMap1 and dropMap2
  size_t numLocalMatrixRows = Ain->getLocalNumRows();

  for (size_t row = 0; row < numLocalMatrixRows; row++) {
    GO grid          = Ain->getRowMap()->getGlobalElement(row);  // global row id
    bool rowIsInMap1 = finalDropMap1->isNodeGlobalElement(grid);
    bool rowIsInMap2 = finalDropMap2->isNodeGlobalElement(grid);

    // extract row information from input matrix
    auto lclMat  = Ain->getLocalMatrixHost();
    auto rowView = lclMat.row(row);

    // just copy all values in output
    Teuchos::ArrayRCP<GO> indout(rowView.length, Teuchos::ScalarTraits<GO>::zero());
    Teuchos::ArrayRCP<SC> valout(rowView.length, Teuchos::ScalarTraits<SC>::zero());

    size_t nNonzeros = 0;
    for (LO jj = 0; jj < rowView.length; ++jj) {
      LO lcid          = rowView.colidx(jj);
      GO gcid          = Ain->getColMap()->getGlobalElement(lcid);  // global column id
      auto val         = rowView.value(jj);
      bool colIsInMap1 = finalDropMap1->isNodeGlobalElement(gcid);
      bool colIsInMap2 = finalDropMap2->isNodeGlobalElement(gcid);

      if ((rowIsInMap1 && colIsInMap2) || (rowIsInMap2 && colIsInMap1)) {
        // do nothing == drop this entry
      } else {
        indout[nNonzeros] = gcid;
        valout[nNonzeros] = val;
        nNonzeros++;
      }
    }
    indout.resize(nNonzeros);
    valout.resize(nNonzeros);
    Aout->insertGlobalValues(Ain->getRowMap()->getGlobalElement(row), indout.view(0, indout.size()), valout.view(0, valout.size()));
  }

  Aout->fillComplete(Ain->getDomainMap(), Ain->getRangeMap());

  // copy block size information
  Aout->SetFixedBlockSize(Ain->GetFixedBlockSize());

  GetOStream(Statistics0, 0) << "Nonzeros in A (input): " << Ain->getGlobalNumEntries() << ", Nonzeros after filtering A: " << Aout->getGlobalNumEntries() << std::endl;

  currentLevel.Set("A", Aout, this);
}  // BuildBasedOnMapPair

}  // namespace MueLu

#endif  // MUELU_SEGREGATEDAFACTORY_DEF_HPP
