// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_DROPNEGATIVEENTRIESFACTORY_DEF_HPP
#define MUELU_DROPNEGATIVEENTRIESFACTORY_DEF_HPP

#include <Xpetra_CrsGraph.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_CrsGraphFactory.hpp>

#include "MueLu_DropNegativeEntriesFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> DropNegativeEntriesFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
#undef SET_VALID_ENTRY

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A used for filtering");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void DropNegativeEntriesFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "A");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void DropNegativeEntriesFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Matrix filtering (springs)", currentLevel);

  RCP<Matrix> Ain = Get<RCP<Matrix> >(currentLevel, "A");

  LocalOrdinal nDofsPerNode = Ain->GetFixedBlockSize();

  // create new empty Operator
  Teuchos::RCP<Matrix> Aout = MatrixFactory::Build(Ain->getRowMap(), Ain->getGlobalMaxNumRowEntries());

  size_t numLocalRows = Ain->getLocalNumRows();
  for (size_t row = 0; row < numLocalRows; row++) {
    GlobalOrdinal grid = Ain->getRowMap()->getGlobalElement(row);

    int rDofID = Teuchos::as<int>(grid % nDofsPerNode);

    // extract row information from input matrix
    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> vals;
    Ain->getLocalRowView(row, indices, vals);

    // just copy all values in output
    Teuchos::ArrayRCP<GlobalOrdinal> indout(indices.size(), Teuchos::ScalarTraits<GlobalOrdinal>::zero());
    Teuchos::ArrayRCP<Scalar> valout(indices.size(), Teuchos::ScalarTraits<Scalar>::zero());

    size_t nNonzeros = 0;
    for (size_t i = 0; i < (size_t)indices.size(); i++) {
      GlobalOrdinal gcid = Ain->getColMap()->getGlobalElement(indices[i]);  // global column id

      int cDofID = Teuchos::as<int>(gcid % nDofsPerNode);
      if (rDofID == cDofID && Teuchos::ScalarTraits<Scalar>::magnitude(vals[i]) >= Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero())) {
        indout[nNonzeros] = gcid;
        valout[nNonzeros] = vals[i];
        nNonzeros++;
      }
    }
    indout.resize(nNonzeros);
    valout.resize(nNonzeros);

    Aout->insertGlobalValues(Ain->getRowMap()->getGlobalElement(row), indout.view(0, indout.size()), valout.view(0, valout.size()));
  }

  Aout->fillComplete(Ain->getDomainMap(), Ain->getRangeMap());

  // copy block size information
  Aout->SetFixedBlockSize(nDofsPerNode);

  GetOStream(Statistics0, 0) << "Nonzeros in A (input): " << Ain->getGlobalNumEntries() << ", Nonzeros after filtering A: " << Aout->getGlobalNumEntries() << std::endl;

  Set(currentLevel, "A", Aout);
}

}  // namespace MueLu

#endif  // MUELU_DROPNEGATIVEENTRIESFACTORY_DEF_HPP
