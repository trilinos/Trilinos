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
//#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> SegregatedAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
#undef  SET_VALID_ENTRY

    validParamList->set< RCP<const FactoryBase> >("A",              Teuchos::null, "Generating factory of the matrix A used for filtering");

    validParamList->set< std::string >           ("map: name", "", "Name of map (Xpetra::Map) provided by user containing the special DOFs.");
    validParamList->set< std::string >           ("map: factory", "null", "Name of generating factory for 'map: name'");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void SegregatedAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "A");

    const ParameterList& pL = GetParameterList();
    std::string mapName                       = pL.get<std::string> ("map: name");
    std::string mapFact                       = pL.get<std::string> ("map: factory");
    if(mapFact == "NoFactory") {
      currentLevel.DeclareInput(mapName,MueLu::NoFactory::get(),this);
    } else if (mapFact != "null") {
      Teuchos::RCP<const FactoryBase> fact = currentLevel.GetFactoryManager()->GetFactory(mapFact);
      currentLevel.DeclareInput(mapName,fact.get(),this);
    } else {
      currentLevel.DeclareInput(mapName,NULL,this);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void SegregatedAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
    FactoryMonitor m(*this, "Matrix filtering (segregation)", currentLevel);


    RCP<Matrix> Ain = Get< RCP<Matrix> >(currentLevel, "A");

    const ParameterList& pL = GetParameterList();
    std::string mapName  = pL.get<std::string> ("map: name");
    std::string mapFact  = pL.get<std::string> ("map: factory");

    Teuchos::RCP<const FactoryBase> mapFactory = Teuchos::null;
    if(mapFact == "NoFactory") {
      mapFactory = MueLu::NoFactory::getRCP();
    } else if (mapFact != "null") {
      mapFactory = currentLevel.GetFactoryManager()->GetFactory(mapFact);
    }

    Teuchos::RCP<Map> map = currentLevel.Get< Teuchos::RCP<Map> >(mapName,mapFactory.get());

    // create new empty Operator
    Teuchos::RCP<Matrix> Aout = MatrixFactory::Build(Ain->getRowMap(), Ain->getGlobalMaxNumRowEntries());

    size_t numLocalRows = Ain->getLocalNumRows();
    for(size_t row=0; row<numLocalRows; row++) { // how can i replace this by a parallel for?
      GlobalOrdinal grid = Ain->getRowMap()->getGlobalElement(row); // global row id
      bool isInMap = map->isNodeGlobalElement(grid);

      // extract row information from input matrix
      Teuchos::ArrayView<const LocalOrdinal> indices;
      Teuchos::ArrayView<const Scalar> vals;
      Ain->getLocalRowView(row, indices, vals);

      // just copy all values in output
      Teuchos::ArrayRCP<GlobalOrdinal> indout(indices.size(),Teuchos::ScalarTraits<GlobalOrdinal>::zero());
      Teuchos::ArrayRCP<Scalar>        valout(indices.size(),Teuchos::ScalarTraits<Scalar>::zero());

      size_t nNonzeros = 0;
      for(size_t i=0; i<(size_t)indices.size(); i++) { // or can be parallelize this loop?
        GlobalOrdinal gcid = Ain->getColMap()->getGlobalElement(indices[i]); // global column id
        bool isInMap2 = map->isNodeGlobalElement(gcid);

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

} //namespace MueLu

#endif // MUELU_SEGREGATEDAFACTORY_DEF_HPP
