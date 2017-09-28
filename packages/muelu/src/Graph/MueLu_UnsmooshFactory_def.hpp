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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef PACKAGES_MUELU_SRC_GRAPH_MUELU_UNSMOOSHFACTORY_DEF_HPP_
#define PACKAGES_MUELU_SRC_GRAPH_MUELU_UNSMOOSHFACTORY_DEF_HPP_

#include "MueLu_Monitor.hpp"

#include "MueLu_UnsmooshFactory_decl.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  UnsmooshFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::UnsmooshFactory() { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> UnsmooshFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
    validParamList->set< RCP<const FactoryBase> >("A",                  Teuchos::null, "Generating factory for unamalgamated matrix. Row map of (unamalgamted) output prolongation operator should match row map of this A.");
    validParamList->set< RCP<const FactoryBase> >("P",                  Teuchos::null, "Generating factory of the (amalgamated) prolongator P");
    validParamList->set< RCP<const FactoryBase> >("DofStatus",          Teuchos::null, "Generating factory for dofStatus array (usually the VariableDofLaplacdianFactory)");

    validParamList->set< int  >                  ("maxDofPerNode", 1,     "Maximum number of DOFs per node");
    validParamList->set< bool >                  ("fineIsPadded" , false, "true if finest level input matrix is padded");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void UnsmooshFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    //const ParameterList& pL = GetParameterList();
    Input(fineLevel, "A");
    Input(coarseLevel, "P");

    // DofStatus only provided on the finest level (by user)
    // On the coarser levels it is auto-generated using the DBC information from the unamalgamated matrix A
    if(fineLevel.GetLevelID() == 0)
      Input(fineLevel, "DofStatus");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void UnsmooshFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);
    typedef Teuchos::ScalarTraits<SC> STS;

    const ParameterList  & pL = GetParameterList();

    // extract matrices (unamalgamated A and amalgamated P)
    RCP<Matrix> unamalgA = Get< RCP<Matrix> >(fineLevel,   "A");
    RCP<Matrix> amalgP   = Get< RCP<Matrix> >(coarseLevel, "P");

    // extract user parameters
    int maxDofPerNode = pL.get<int> ("maxDofPerNode");
    bool fineIsPadded = pL.get<bool>("fineIsPadded");

    // get dofStatus information
    // On the finest level it is provided by the user. On the coarser levels it is constructed
    // using the DBC information of the matrix A
    Teuchos::Array<char> dofStatus;
    if(fineLevel.GetLevelID() == 0) {
      dofStatus = Get<Teuchos::Array<char> >(fineLevel, "DofStatus");
      TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(dofStatus.size()) == Teuchos::as<size_t>(unamalgA->getRowMap()->getNodeNumElements()), MueLu::Exceptions::RuntimeError,"MueLu::UnsmooshFactory::Build: User provided dofStatus on level 0 does not fit to size of unamalgamted A");
    } else {
      // dof status is the dirichlet information of unsmooshed/unamalgamated A (fine level)
      dofStatus = Teuchos::Array<char>(unamalgA->getRowMap()->getNodeNumElements() /*amalgP->getRowMap()->getNodeNumElements() * maxDofPerNode*/,'s');

      bool bHasZeroDiagonal = false;
      Teuchos::ArrayRCP<const bool> dirOrNot = MueLu::Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DetectDirichletRowsExt(*unamalgA,bHasZeroDiagonal,STS::magnitude(0.5));

      TEUCHOS_TEST_FOR_EXCEPTION(dirOrNot.size() != dofStatus.size(), MueLu::Exceptions::RuntimeError,"MueLu::UnsmooshFactory::Build: inconsistent number of coarse DBC array and dofStatus array. dirOrNot.size() = " << dirOrNot.size() << " dofStatus.size() = " << dofStatus.size());
      for(decltype(dirOrNot.size()) i = 0; i < dirOrNot.size(); ++i) {
        if(dirOrNot[i] == true) dofStatus[i] = 'p';
      }
    }

    // TODO: TAW the following check is invalid for SA-AMG based input prolongators
    //TEUCHOS_TEST_FOR_EXCEPTION(amalgP->getDomainMap()->isSameAs(*amalgP->getColMap()) == false, MueLu::Exceptions::RuntimeError,"MueLu::UnsmooshFactory::Build: only support for non-overlapping aggregates. (column map of Ptent must be the same as domain map of Ptent)");

    // extract CRS information from amalgamated prolongation operator
    Teuchos::ArrayRCP<const size_t> amalgRowPtr(amalgP->getNodeNumRows());
    Teuchos::ArrayRCP<const LocalOrdinal> amalgCols(amalgP->getNodeNumEntries());
    Teuchos::ArrayRCP<const Scalar> amalgVals(amalgP->getNodeNumEntries());
    Teuchos::RCP<CrsMatrixWrap> amalgPwrap = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(amalgP);
    Teuchos::RCP<CrsMatrix> amalgPcrs = amalgPwrap->getCrsMatrix();
    amalgPcrs->getAllValues(amalgRowPtr, amalgCols, amalgVals);

    // calculate number of dof rows for new prolongator
    size_t paddedNrows = amalgP->getRowMap()->getNodeNumElements() * Teuchos::as<size_t>(maxDofPerNode);

    // reserve CSR arrays for new prolongation operator
    Teuchos::ArrayRCP<size_t> newPRowPtr(paddedNrows+1);
    Teuchos::ArrayRCP<LocalOrdinal> newPCols(amalgP->getNodeNumEntries() * maxDofPerNode);
    Teuchos::ArrayRCP<Scalar> newPVals(amalgP->getNodeNumEntries() * maxDofPerNode);

    size_t rowCount = 0; // actual number of (local) in unamalgamated prolongator
    if(fineIsPadded == true || fineLevel.GetLevelID() > 0) {

      // build prolongation operator for padded fine level matrices.
      // Note: padded fine level dofs are transfered by injection.
      // That is, these interpolation stencils do not take averages of
      // coarse level variables. Further, fine level Dirichlet points
      // also use injection.

      size_t cnt = 0; // local id counter
      for (decltype(amalgRowPtr.size()) i = 0; i < amalgRowPtr.size() - 1; i++) {
        // determine number of entries in amalgamated dof row i
        size_t rowLength = amalgRowPtr[i+1] - amalgRowPtr[i];

        // loop over dofs per node (unamalgamation)
        for(int j = 0; j < maxDofPerNode; j++) {
          newPRowPtr[i*maxDofPerNode+j] = cnt;
          if (dofStatus[i*maxDofPerNode+j] == 's') { // add only "standard" dofs to unamalgamated prolongator
            // loop over column entries in amalgamated P
            for (size_t k = 0; k < rowLength; k++) {
              newPCols[cnt  ] = amalgCols[k+amalgRowPtr[i]] * maxDofPerNode + j;
              newPVals[cnt++] = amalgVals[k+amalgRowPtr[i]];
            }

          }
        }
      }

      newPRowPtr[paddedNrows] = cnt; // close row CSR array
      rowCount = paddedNrows;
    } else {
      // Build prolongation operator for non-padded fine level matrices.
      // Need to map from non-padded dofs to padded dofs. For this, look
      // at the status array and skip padded dofs.

      size_t cnt = 0; // local id counter

      for (decltype(amalgRowPtr.size()) i = 0; i < amalgRowPtr.size() - 1; i++) {
        // determine number of entries in amalgamated dof row i
        size_t rowLength = amalgRowPtr[i+1] - amalgRowPtr[i];

        // loop over dofs per node (unamalgamation)
        for(int j = 0; j < maxDofPerNode; j++) {
          // no interpolation for padded fine dofs as they do not exist

          if (dofStatus[i*maxDofPerNode+j] == 's') { // add only "standard" dofs to unamalgamated prolongator
            newPRowPtr[rowCount++] = cnt;
            // loop over column entries in amalgamated P
            for (size_t k = 0; k < rowLength; k++) {
              newPCols[cnt  ] = amalgCols[k+amalgRowPtr[i]] * maxDofPerNode + j;
              newPVals[cnt++] = amalgVals[k+amalgRowPtr[i]];
            }

          }
          if (dofStatus[i*maxDofPerNode+j] == 'd') { // Dirichlet handling
            newPRowPtr[rowCount++] = cnt;
          }
        }
      }
      newPRowPtr[rowCount] = cnt; // close row CSR array
    } // fineIsPadded == false

    // generate coarse domain map
    // So far no support for gid offset or strided maps. This information
    // could be gathered easily from the unamalgamated fine level operator A.
    std::vector<size_t> stridingInfo(1,maxDofPerNode);

    GlobalOrdinal nCoarseDofs = amalgP->getDomainMap()->getNodeNumElements() * maxDofPerNode;
    GlobalOrdinal indexBase   = amalgP->getDomainMap()->getIndexBase();
    RCP<const Map> coarseDomainMap = StridedMapFactory::Build(amalgP->getDomainMap()->lib(),
        Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
        nCoarseDofs,
        indexBase,
        stridingInfo,
        amalgP->getDomainMap()->getComm(),
        -1 /* stridedBlockId */,
        0  /*domainGidOffset */);

    size_t nColCoarseDofs = Teuchos::as<size_t>(amalgP->getColMap()->getNodeNumElements() * maxDofPerNode);
    Teuchos::Array<GlobalOrdinal> unsmooshColMapGIDs(nColCoarseDofs);
    for(size_t c = 0; c < amalgP->getColMap()->getNodeNumElements(); ++c) {
      GlobalOrdinal gid = amalgP->getColMap()->getGlobalElement(c) * maxDofPerNode;
      for(int i = 0; i < maxDofPerNode; ++i) {
        unsmooshColMapGIDs[c * maxDofPerNode + i] = gid + i;
      }
    }
    Teuchos::RCP<Map> coarseColMap = MapFactory::Build(amalgP->getDomainMap()->lib(),
        Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
        unsmooshColMapGIDs(), //View,
        indexBase,
        amalgP->getDomainMap()->getComm());

    // Assemble unamalgamated P
    Teuchos::RCP<CrsMatrix> unamalgPCrs = CrsMatrixFactory::Build(unamalgA->getRowMap(),coarseColMap, 3);
    for (decltype(rowCount) i = 0; i < rowCount; i++) {
      unamalgPCrs->insertLocalValues(i, newPCols.view(newPRowPtr[i],newPRowPtr[i+1]-newPRowPtr[i]),
          newPVals.view(newPRowPtr[i],newPRowPtr[i+1]-newPRowPtr[i]));
    }
    unamalgPCrs->fillComplete(coarseDomainMap,unamalgA->getRowMap());

    Teuchos::RCP<Matrix> unamalgP = Teuchos::rcp(new CrsMatrixWrap(unamalgPCrs));

    Set(coarseLevel,"P",unamalgP);
  }


} /* MueLu */


#endif /* PACKAGES_MUELU_SRC_GRAPH_MUELU_UNSMOOSHFACTORY_DEF_HPP_ */
