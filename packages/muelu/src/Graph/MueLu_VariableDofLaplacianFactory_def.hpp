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
#ifndef PACKAGES_MUELU_SRC_GRAPH_MUELU_VARIABLEDOFLAPLACIANFACTORY_DEF_HPP_
#define PACKAGES_MUELU_SRC_GRAPH_MUELU_VARIABLEDOFLAPLACIANFACTORY_DEF_HPP_


#include "MueLu_Monitor.hpp"

#include "MueLu_VariableDofLaplacianFactory_decl.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> VariableDofLaplacianFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
/*
#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("aggregation: drop tol");
    SET_VALID_ENTRY("aggregation: Dirichlet threshold");
    SET_VALID_ENTRY("aggregation: drop scheme");
    {
      typedef Teuchos::StringToIntegralParameterEntryValidator<int> validatorType;
      validParamList->getEntry("aggregation: drop scheme").setValidator(
        rcp(new validatorType(Teuchos::tuple<std::string>("classical", "distance laplacian"), "aggregation: drop scheme")));
    }
#undef  SET_VALID_ENTRY
    validParamList->set< bool >                  ("lightweight wrap",           true, "Experimental option for lightweight graph access");
*/

    validParamList->set< RCP<const FactoryBase> >("A",                  Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("Coordinates",        Teuchos::null, "Generating factory for Coordinates");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  VariableDofLaplacianFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::VariableDofLaplacianFactory() { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void VariableDofLaplacianFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
    //const ParameterList& pL = GetParameterList();
    Input(currentLevel, "A");
    Input(currentLevel, "Coordinates");

    //if (currentLevel.GetLevelID() == 0) // TODO check for finest level (special treatment)
    if (currentLevel.IsAvailable("DofPresent", NoFactory::get())) {
      currentLevel.DeclareInput("DofPresent", NoFactory::get(), this);
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void VariableDofLaplacianFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);
    typedef Teuchos::ScalarTraits<SC> STS;

    const ParameterList  & pL = GetParameterList();

    RCP<Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");

    typedef Xpetra::MultiVector<double,LO,GO,NO> dxMV;
    RCP<dxMV> Coords = Get< RCP<Xpetra::MultiVector<double,LO,GO,NO> > >(currentLevel, "Coordinates");

    int maxDofPerNode = 3;    // TODO add parameter from parameter list for this...
    Scalar dirDropTol = 1e-5; // TODO parameter from parameter list ("ML advnaced Dirichlet: threshold"), should be magnitude type?
    Scalar amalgDropTol = 1.8e-9; // TODO parameter from parameter list ("variable DOF amalgamation: threshold")

    bool bHasZeroDiagonal = false;
    Teuchos::ArrayRCP<const bool> dirOrNot = MueLu::Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DetectDirichletRowsExt(*A,bHasZeroDiagonal,dirDropTol);

    // TODO check availability + content (length)
    Teuchos::ArrayRCP<const bool> dofPresent = currentLevel.Get< Teuchos::ArrayRCP<const bool> >("DofPresent", NoFactory::get());

    // map[k] indicates that the kth dof in the variable dof matrix A would
    // correspond to the map[k]th dof in the padded system. If, i.e., it is
    // map[35] = 39 then dof no 35 in the variable dof matrix A corresponds to
    // row map id 39 in an imaginary padded matrix Apadded.
    // The padded system is never built but would be the associated matrix if
    // every node had maxDofPerNode dofs.
    std::vector<LocalOrdinal> map(A->getNodeNumRows());
    this->buildPaddedMap(dofPresent, map, A->getNodeNumRows());

    //GetOStream(Parameters0) << "lightweight wrap = " << doExperimentalWrap << std::endl;


    // filled with local node ids (associated with each dof)
    std::vector<LocalOrdinal> myLocalNodeIds(A->getColMap()->getNodeNumElements()); // possible maximum (we need the ghost nodes, too)

    // assign the local node ids for the ghosted nodes
    size_t nLocalNodes, nLocalPlusGhostNodes;
    this->assignGhostLocalNodeIds(A->getRowMap(), A->getColMap(), myLocalNodeIds, map, maxDofPerNode, nLocalNodes, nLocalPlusGhostNodes, A->getRowMap()->getComm());

    // fill nodal maps

    Teuchos::ArrayView< const GlobalOrdinal > myGids = A->getColMap()->getNodeElementList();

    // vector containing row/col gids of amalgamated matrix (with holes)

    size_t nLocalDofs = A->getRowMap()->getNodeNumElements();
    size_t nLocalPlusGhostDofs = A->getColMap()->getNodeNumElements();

    std::vector<GlobalOrdinal> amalgRowMapGIDs(nLocalNodes);
    std::vector<GlobalOrdinal> amalgColMapGIDs(nLocalPlusGhostNodes);

    // initialize
    size_t count = 0;
    if (nLocalDofs > 0) {
      amalgRowMapGIDs[count] = myGids[0];
      amalgColMapGIDs[count] = myGids[0];
      count++;
    }

    // dofs belonging to the same node must be consecutively only for the local part of myLocalNodeIds[]
    // fill amalgRowMapGIDs + local part of amalgColMapGIDs
    for(size_t i = 1; i < nLocalDofs; i++) {
      if(myLocalNodeIds[i] != myLocalNodeIds[i-1]) {
        amalgRowMapGIDs[count] = myGids[i];
        amalgColMapGIDs[count] = myGids[i];
        count++;
      }
    }

    ArrayView<GlobalOrdinal> amalgRowMapGIDsView(amalgRowMapGIDs.size() ? &amalgRowMapGIDs[0] : 0, amalgRowMapGIDs.size());
    Teuchos::RCP<Map> amalgRowMap = MapFactory::Build(A->getRowMap()->lib(),
               Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
               amalgRowMapGIDsView,
               A->getRowMap()->getIndexBase(),
               A->getRowMap()->getComm());

    // TODO check me!
    // copy GIDs from nodal vector to dof vector
    Teuchos::RCP<Vector> localGIDsSrc = VectorFactory::Build(A->getRowMap(0),true);
    Teuchos::ArrayRCP< Scalar > localGIDsSrcData = localGIDsSrc->getDataNonConst(0);

    for (int i = 0; i < myLocalNodeIds.size(); i++)
      localGIDsSrcData[i] = amalgColMapGIDs[ myLocalNodeIds[i]];

    Teuchos::RCP<Import> importer = ImportFactory::Build(A->getRowMap(), A->getColMap());
    Teuchos::RCP<Vector> localGIDsTarget = VectorFactory::Build(A->getColMap(0),true);

    localGIDsTarget->doImport(*localGIDsSrc, *importer, Xpetra::INSERT);
    Teuchos::ArrayRCP< const Scalar > localGIDsTargetData = localGIDsTarget->getData(0);

    // copy from dof vector to nodal vector
    for (int i = 0; i < myLocalNodeIds.size(); i++)
      amalgColMapGIDs[ myLocalNodeIds[i]] = localGIDsTargetData[i];

    ArrayView<GlobalOrdinal> amalgColMapGIDsView(amalgColMapGIDs.size() ? &amalgRowMapGIDs[0] : 0, amalgColMapGIDs.size());
    Teuchos::RCP<Map> amalgColMap = MapFactory::Build(A->getColMap()->lib(),
               Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
               amalgColMapGIDsView,
               A->getColMap()->getIndexBase(),
               A->getColMap()->getComm());

    // end fill nodal maps

    // start variable dof amalgamation

    Teuchos::RCP<CrsMatrixWrap> Awrap = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A);
    Teuchos::RCP<CrsMatrix> Acrs = Awrap->getCrsMatrix();

    Teuchos::ArrayRCP<const size_t> rowptr(Acrs->getNodeNumRows());
    Teuchos::ArrayRCP<const LocalOrdinal> colind(Acrs->getNodeNumEntries());
    Teuchos::ArrayRCP<const Scalar> values(Acrs->getNodeNumEntries());
    Acrs->getAllValues(rowptr, colind, values);

    // create arrays for amalgamated matrix
    Teuchos::ArrayRCP<size_t> amalgRowPtr(nLocalNodes+1);
    Teuchos::ArrayRCP<LocalOrdinal> amalgCols(rowptr[rowptr.size()-1]);
    //Teuchos::ArrayRCP<const Scalar> values(Acrs->getNodeNumEntries());

    size_t nNonZeros = 0;
    std::vector<bool> isNonZero(nLocalPlusGhostDofs,false);
    std::vector<size_t> nonZeroList(nLocalPlusGhostDofs);  // ???

    // also used in DetectDirichletExt
    Teuchos::RCP<Vector> diagVec = VectorFactory::Build(A->getRowMap());
    A->getLocalDiagCopy(*diagVec);
    Teuchos::ArrayRCP< const Scalar > diagVecData = diagVec->getData(0);


    LocalOrdinal oldBlockRow = 0;
    LocalOrdinal blockRow, blockColumn;
    size_t newNzs = 0;
    amalgRowPtr[0] = newNzs;

    bool doNotDrop = false;
    if (amalgDropTol == Teuchos::ScalarTraits<Scalar>::zero()) doNotDrop = true;
    if (values.size() == 0) doNotDrop = true;

    for(size_t i = 0; i < rowptr.size()-1; i++) {
      blockRow = std::floor<LocalOrdinal>( map[i] / maxDofPerNode);
      if (blockRow != oldBlockRow) {
        // zero out info recording nonzeros in oldBlockRow
        for(size_t j = 0; j < nNonZeros; j++) isNonZero[nonZeroList[j]] = false;
        nNonZeros = 0;
        amalgRowPtr[blockRow] = newNzs; // record start of next row
      }
      for (size_t j = rowptr[i]; j < rowptr[i+1]; j++) {
        if(doNotDrop == true ||
            ( STS::magnitude(values[j] / sqrt(STS::magnitude(diagVecData[i]) * STS::magnitude(diagVecData[colind[j]]))   ) >= amalgDropTol )) {
          blockColumn = myLocalNodeIds[colind[j]];
          if(isNonZero[blockColumn] == false) {
            isNonZero[blockColumn] = true;
            nonZeroList[nNonZeros++] = blockColumn;
            amalgCols[newNzs++] = blockColumn;
          }
        }
      }
      oldBlockRow = blockRow;
    }
    amalgRowPtr[blockRow+1] = newNzs;

    TEUCHOS_TEST_FOR_EXCEPTION((blockRow+1 != nLocalNodes) && (nLocalNodes !=0), MueLu::Exceptions::RuntimeError, "VariableDofsPerNodeAmalgamation: error, computed # block rows (" << blockRow+1 <<") != nLocalNodes (" << nLocalNodes <<")");

    amalgCols.resize(amalgRowPtr[nLocalNodes]);

    // end variableDofAmalg
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void VariableDofLaplacianFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::buildPaddedMap(const Teuchos::ArrayRCP<const bool> & dofPresent, std::vector<LocalOrdinal> & map, size_t nDofs) const {
    size_t count = 0;
    for (size_t i = 0; i < dofPresent.size(); i++)
      if(dofPresent[i] == true) map[count++] = Teuchos::as<LocalOrdinal>(i);
    TEUCHOS_TEST_FOR_EXCEPTION(nDofs != count, MueLu::Exceptions::RuntimeError, "VariableDofLaplacianFactory::buildPaddedMap: #dofs in dofPresent does not match the expected value (number of rows of A): " << nDofs << " vs. " << count);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void VariableDofLaplacianFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::assignGhostLocalNodeIds(const Teuchos::RCP<const Map> & rowDofMap, const Teuchos::RCP<const Map> & colDofMap, std::vector<LocalOrdinal> & myLocalNodeIds, const std::vector<LocalOrdinal> & dofMap, size_t maxDofPerNode, size_t& nLocalNodes, size_t& nLocalPlusGhostNodes, Teuchos::RCP< const Teuchos::Comm< int > > comm) const {

    size_t nLocalDofs = rowDofMap->getNodeNumElements();
    size_t nLocalPlusGhostDofs = colDofMap->getNodeNumElements(); // TODO remove parameters

    // create importer for dof-based information
    Teuchos::RCP<Import> importer = ImportFactory::Build(rowDofMap, colDofMap);

    // create a vector living on column map of A (dof based)
    Teuchos::RCP<Vector> localNodeIdsTemp = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(rowDofMap,true);
    Teuchos::RCP<Vector> localNodeIds = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(colDofMap,true);

    // fill local dofs (padded local ids)
    Teuchos::ArrayRCP< Scalar > localNodeIdsTempData = localNodeIdsTemp->getDataNonConst(0);
    for(size_t i = 0; i < localNodeIdsTemp->getLocalLength(); i++)
      localNodeIdsTempData[i] = std::floor<LocalOrdinal>( dofMap[i] / maxDofPerNode );

    localNodeIds->doImport(*localNodeIdsTemp, *importer, Xpetra::INSERT);
    Teuchos::ArrayRCP< const Scalar > localNodeIdsData = localNodeIds->getData(0);

    // Note: localNodeIds contains local ids for the padded version as vector values


    // we use Scalar instead of int as type
    Teuchos::RCP<Vector> myProcTemp = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(rowDofMap,true);
    Teuchos::RCP<Vector> myProc = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(colDofMap,true);

    // fill local dofs (padded local ids)
    Teuchos::ArrayRCP< Scalar > myProcTempData = myProcTemp->getDataNonConst(0);
    for(size_t i = 0; i < myProcTemp->getLocalLength(); i++)
      myProcTempData[i] = Teuchos::as<Scalar>(comm->getRank());
    myProc->doImport(*myProcTemp, *importer, Xpetra::INSERT);
    Teuchos::ArrayRCP<Scalar > myProcData = myProc->getDataNonConst(0); // we have to modify the data (therefore the non-const version)

    // At this point, the ghost part of localNodeIds corresponds to the local ids
    // associated with the current owning processor. We want to convert these to
    // local ids associated with the processor on which these are ghosts.
    // Thus we have to re-number them. In doing this re-numbering we must make sure
    // that we find all ghosts with the same id & proc and assign a unique local
    // id to this group (id&proc). To do this find, we sort all ghost entries in
    // localNodeIds that are owned by the same processor. Then we can look for
    // duplicates (i.e., several ghost entries corresponding to dofs with the same
    // node id) easily and make sure these are all assigned to the same local id.
    // To do the sorting we'll make a temporary copy of the ghosts via tempId and
    // tempProc and sort this multiple times for each group owned by the same proc.


    std::vector<size_t> location(nLocalPlusGhostDofs - nLocalDofs + 1);
    std::vector<size_t> tempId  (nLocalPlusGhostDofs - nLocalDofs + 1);
    std::vector<size_t> tempProc(nLocalPlusGhostDofs - nLocalDofs + 1);

    size_t notProcessed = nLocalDofs; // iteration index over all ghosted dofs
    size_t tempIndex = 0;
    size_t first = tempIndex;
    int neighbor;

    while (notProcessed < nLocalPlusGhostDofs) {
      neighbor = Teuchos::as<int>(myProcData[notProcessed]); // get processor id of not-processed element
      first    = tempIndex;
      location[tempIndex] = notProcessed;
      tempId[tempIndex++] = localNodeIdsData[notProcessed];
      myProcData[notProcessed] = Teuchos::as<Scalar>(-1 - neighbor);

      for(size_t i = notProcessed + 1; i < nLocalPlusGhostDofs; i++) {
        if(myProcData[i] == neighbor) {
          location[tempIndex] = i;
          tempId[tempIndex++] = localNodeIdsData[i];
          myProcData[i] = -1; // mark as visited
        }
      }
      this->MueLu_az_sort(&(tempId[first]), tempIndex - first, &(location[first]), NULL);
      for(size_t i = first; i < tempIndex; i++) tempProc[i] = neighbor;

      // increment index. Find next notProcessed dof index corresponding to first non-visited element
      notProcessed++;
      while ( (notProcessed < nLocalPlusGhostDofs) && (myProcData[notProcessed] < 0))
        notProcessed++;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(tempIndex != nLocalPlusGhostDofs-nLocalDofs, MueLu::Exceptions::RuntimeError,"Number of nonzero ghosts is inconsistent.");

    // Now assign ids to all ghost nodes (giving the same id to those with the
    // same myProc[] and the same local id on the proc that actually owns the
    // variable associated with the ghost

    nLocalNodes = 0; // initialize return value
    if(nLocalDofs > 0) nLocalNodes = localNodeIdsData[nLocalDofs-1] + 1;

    nLocalPlusGhostNodes = nLocalNodes; // initialize return value
    if(nLocalDofs < nLocalPlusGhostDofs) nLocalPlusGhostNodes++; // 1st ghost node is unique (not accounted for). number will be increased later, if there are more ghost nodes

    // check if two adjacent ghost dofs correspond to different nodes. To do this,
    // check if they are from different processors or whether they have different
    // local node ids

    // loop over all (remaining) ghost dofs
    size_t lagged = -1;
    for (size_t i = nLocalDofs+1; i < nLocalPlusGhostDofs; i++) {
      lagged = nLocalPlusGhostNodes-1;

      // i is a new unique ghost node (not already accounted for)
      if ((tempId[i-nLocalDofs] != tempId[i-1-nLocalDofs]) ||
          (tempProc[i-nLocalDofs] != tempProc[i-1-nLocalDofs]))
        nLocalPlusGhostNodes++; // update number of ghost nodes
      tempId[i-1-nLocalDofs] = lagged;
    }
    if (nLocalPlusGhostDofs > nLocalDofs)
      tempId[nLocalPlusGhostDofs-1-nLocalDofs] = nLocalPlusGhostNodes - 1;

    // fill myLocalNodeIds array. Start with local part (not ghosted)
    for(size_t i = 0; i < nLocalDofs; i++)
      myLocalNodeIds[i] = std::floor<LocalOrdinal>( dofMap[i] / maxDofPerNode );

    // copy ghosted nodal ids into myLocalNodeIds
    for(size_t i = nLocalDofs; i < nLocalPlusGhostDofs; i++)
      myLocalNodeIds[location[i-nLocalDofs]] = tempId[i-nLocalDofs];

  }


  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void VariableDofLaplacianFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MueLu_az_sort(size_t list[], size_t N, size_t list2[], Scalar list3[]) const {
    /* local variables */

    size_t RR, K, l, r, j, flag, i;
    size_t RR2;
    Scalar RR3;

    /*********************** execution begins ******************************/

    if (N <= 1) return;

    l   = N / 2 + 1;
    r   = N - 1;
    l   = l - 1;
    RR  = list[l - 1];
    K   = list[l - 1];

    if ((list2 != NULL) && (list3 != NULL)) {
      RR2 = list2[l - 1];
      RR3 = list3[l - 1];
      while (r != 0) {
        j = l;
        flag = 1;

        while (flag == 1) {
          i = j;
          j = j + j;

          if (j > r + 1)
            flag = 0;
          else {
            if (j < r + 1)
              if (list[j] > list[j - 1]) j = j + 1;

            if (list[j - 1] > K) {
              list[ i - 1] = list[ j - 1];
              list2[i - 1] = list2[j - 1];
              list3[i - 1] = list3[j - 1];
            }
            else {
              flag = 0;
            }
          }
        }

        list[ i - 1] = RR;
        list2[i - 1] = RR2;
        list3[i - 1] = RR3;

        if (l == 1) {
          RR  = list [r];
          RR2 = list2[r];
          RR3 = list3[r];

          K = list[r];
          list[r ] = list[0];
          list2[r] = list2[0];
          list3[r] = list3[0];
          r = r - 1;
        }
        else {
          l   = l - 1;
          RR  = list[ l - 1];
          RR2 = list2[l - 1];
          RR3 = list3[l - 1];
          K   = list[l - 1];
        }
      }

      list[ 0] = RR;
      list2[0] = RR2;
      list3[0] = RR3;
    }
    else if (list2 != NULL) {
      RR2 = list2[l - 1];
      while (r != 0) {
        j = l;
        flag = 1;

        while (flag == 1) {
          i = j;
          j = j + j;

          if (j > r + 1)
            flag = 0;
          else {
            if (j < r + 1)
              if (list[j] > list[j - 1]) j = j + 1;

            if (list[j - 1] > K) {
              list[ i - 1] = list[ j - 1];
              list2[i - 1] = list2[j - 1];
            }
            else {
              flag = 0;
            }
          }
        }

        list[ i - 1] = RR;
        list2[i - 1] = RR2;

        if (l == 1) {
          RR  = list [r];
          RR2 = list2[r];

          K = list[r];
          list[r ] = list[0];
          list2[r] = list2[0];
          r = r - 1;
        }
        else {
          l   = l - 1;
          RR  = list[ l - 1];
          RR2 = list2[l - 1];
          K   = list[l - 1];
        }
      }

      list[ 0] = RR;
      list2[0] = RR2;
    }
    else if (list3 != NULL) {
      RR3 = list3[l - 1];
      while (r != 0) {
        j = l;
        flag = 1;

        while (flag == 1) {
          i = j;
          j = j + j;

          if (j > r + 1)
            flag = 0;
          else {
            if (j < r + 1)
              if (list[j] > list[j - 1]) j = j + 1;

            if (list[j - 1] > K) {
              list[ i - 1] = list[ j - 1];
              list3[i - 1] = list3[j - 1];
            }
            else {
              flag = 0;
            }
          }
        }

        list[ i - 1] = RR;
        list3[i - 1] = RR3;

        if (l == 1) {
          RR  = list [r];
          RR3 = list3[r];

          K = list[r];
          list[r ] = list[0];
          list3[r] = list3[0];
          r = r - 1;
        }
        else {
          l   = l - 1;
          RR  = list[ l - 1];
          RR3 = list3[l - 1];
          K   = list[l - 1];
        }
      }

      list[ 0] = RR;
      list3[0] = RR3;

    }
    else {
      while (r != 0) {
        j = l;
        flag = 1;

        while (flag == 1) {
          i = j;
          j = j + j;

          if (j > r + 1)
            flag = 0;
          else {
            if (j < r + 1)
              if (list[j] > list[j - 1]) j = j + 1;

            if (list[j - 1] > K) {
              list[ i - 1] = list[ j - 1];
            }
            else {
              flag = 0;
            }
          }
        }

        list[ i - 1] = RR;

        if (l == 1) {
          RR  = list [r];

          K = list[r];
          list[r ] = list[0];
          r = r - 1;
        }
        else {
          l   = l - 1;
          RR  = list[ l - 1];
          K   = list[l - 1];
        }
      }

      list[ 0] = RR;
    }

  }
} /* MueLu */


#endif /* PACKAGES_MUELU_SRC_GRAPH_MUELU_VARIABLEDOFLAPLACIANFACTORY_DEF_HPP_ */
