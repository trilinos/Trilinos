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
#ifndef MUELU_COALESCEDROPFACTORY_DEF_HPP
#define MUELU_COALESCEDROPFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsGraph.hpp>
#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_StridedMap.hpp>
#include <Xpetra_Map.hpp>

#include "MueLu_CoalesceDropFactory_decl.hpp"

#include "MueLu_Exceptions.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_LWGraph.hpp"
#include "MueLu_PreDropFunctionBaseClass.hpp"
#include "MueLu_PreDropFunctionConstVal.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_AmalgamationFactory.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {

    typedef  Teuchos::ScalarTraits<Scalar> TST;
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("A",                  Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Generating factory for UnAmalgamationInfo");
    validParamList->set< RCP<const FactoryBase> >("Coordinates",        Teuchos::null, "Generating factory for Coordinates");
    validParamList->set< bool >                  ("lightweight wrap",   false,         "Experimental option for lightweight graph access");
    validParamList->set< Scalar >                ("Dirichlet detection threshold", TST::zero(), "Threshold for determining whether entries are zero during Dirichlet row detection");
    validParamList->set< Scalar >                ("aggregation threshold", TST::zero(), "Aggregation dropping threshold");
    validParamList->set< std::string >           ("algorithm",          "original",    "Dropping algorithm");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CoalesceDropFactory() : predrop_(Teuchos::null) { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    Input(currentLevel, "A");
    Input(currentLevel, "UnAmalgamationInfo");

    const ParameterList  & pL = GetParameterList();
    if (pL.get<bool>("lightweight wrap") == true && pL.get<std::string>("algorithm") == "laplacian")
      Input(currentLevel, "Coordinates");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    typedef Teuchos::ScalarTraits<Scalar> STS;

    if (predrop_ != Teuchos::null)
      GetOStream(Parameters0, 0) << predrop_->description();

    RCP<Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");

    const ParameterList  & pL = GetParameterList();
    bool doExperimentalWrap = pL.get<bool>("lightweight wrap");

    GetOStream(Parameters0, 0) << "CoalesceDropFactory::Build : lightweight wrap = " << doExperimentalWrap << std::endl;

    if (doExperimentalWrap) {
      std::string algo = pL.get<std::string>("algorithm");

      TEUCHOS_TEST_FOR_EXCEPTION(predrop_ != null    && algo != "original", Exceptions::RuntimeError, "Dropping function must not be provided for \"" << algo << "\" algorithm");
      TEUCHOS_TEST_FOR_EXCEPTION(algo != "original" && algo != "laplacian", Exceptions::RuntimeError, "\"algorithm\" must be one of (original|laplacian)");

      Scalar threshold = Teuchos::as<Scalar>(pL.get<double>("aggregation threshold"));
      GetOStream(Runtime0, 0) << "algorithm = \"" << algo << "\": threshold = " << threshold << std::endl;
      Set(currentLevel, "Filtering", (threshold != STS::zero()));

      const typename STS::magnitudeType dirichletThreshold = STS::magnitude(pL.get<Scalar>("Dirichlet detection threshold"));

      GlobalOrdinal numDropped = 0, numTotal = 0;
      if (algo == "original") {
        if (predrop_ == null) {
          // ap: this is a hack: had to declare predrop_ as mutable
          predrop_ = rcp(new PreDropFunctionConstVal(threshold));
        }

        if (predrop_ != null) {
          RCP<PreDropFunctionConstVal> predropConstVal = rcp_dynamic_cast<PreDropFunctionConstVal>(predrop_);
          TEUCHOS_TEST_FOR_EXCEPTION(predropConstVal == Teuchos::null, Exceptions::BadCast,
                                     "MueLu::CoalesceFactory::Build: cast to PreDropFunctionConstVal failed.");
          // If a user provided a predrop function, it overwrites the XML threshold parameter
          Scalar newt = predropConstVal->GetThreshold();
          if (newt != threshold) {
            GetOStream(Warnings0,0) << "switching threshold parameter from " << threshold << " (list) to " << newt << " (user function" << std::endl;
            threshold = newt;
          }
        }

        // At this points we either have
        //     (predrop_ != null)
        // Therefore, it is sufficient to check only threshold

        // Detect and record rows that correspond to Dirichlet boundary conditions
        const ArrayRCP<const bool > boundaryNodes = MueLu::Utils<SC,LO,GO,NO>::DetectDirichletRows(*A, dirichletThreshold);

        if ( (A->GetFixedBlockSize() == 1) && (threshold == STS::zero()) ) {
          // Case 1:  scalar problem, no dropping => just use matrix graph
          RCP<GraphBase> graph = rcp(new Graph(A->getCrsGraph(), "graph of A"));
          graph->SetBoundaryNodeMap(boundaryNodes);

          if (GetVerbLevel() & Statistics0) {
            GO numLocalBoundaryNodes  = 0;
            GO numGlobalBoundaryNodes = 0;
            for (LO i = 0; i < boundaryNodes.size(); ++i)
              if (boundaryNodes[i])
                numLocalBoundaryNodes++;
            RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
            sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
            GetOStream(Statistics0, 0) << "Detected " << numGlobalBoundaryNodes << " Dirichlet nodes" << std::endl;
          }

          Set(currentLevel, "DofsPerNode", 1);
          Set(currentLevel, "Graph", graph);

        } else if ( (A->GetFixedBlockSize() == 1) && threshold != STS::zero() ) {
          // Case 2:  scalar problem with dropping => record the column indices of undropped entries, but still use original
          //                                          graph's map information, e.g., whether index is local

          // allocate space for the local graph
          ArrayRCP<LocalOrdinal> rows    = ArrayRCP<LO>(A->getNodeNumRows()+1);
          ArrayRCP<LocalOrdinal> columns = ArrayRCP<LO>(A->getNodeNumEntries());

          RCP<Vector> ghostedDiag = MueLu::Utils<SC,LO,GO,NO>::GetMatrixOverlappedDiagonal(*A);
          const ArrayRCP<const SC> ghostedDiagVals = ghostedDiag->getData(0);
          const ArrayRCP<bool> amalgBoundaryNodes(A->getNodeNumRows(),false);

          LocalOrdinal realnnz = 0;

          rows[0] = 0;
          for (LocalOrdinal row = 0; row < Teuchos::as<LocalOrdinal>(A->getRowMap()->getNodeNumElements()); ++row) {
            size_t nnz = A->getNumEntriesInLocalRow(row);
            ArrayView<const LocalOrdinal> indices;
            ArrayView<const Scalar>       vals;
            A->getLocalRowView(row, indices, vals);

            //FIXME the current predrop function uses the following
            //FIXME    if(std::abs(vals[k]) > std::abs(threshold_) || grow == gcid )
            //FIXME but the threshold doesn't take into account the rows' diagonal entries
            //FIXME For now, hardwiring the dropping in here

            LocalOrdinal rownnz = 0;
            for (LocalOrdinal colID = 0; colID < Teuchos::as<LocalOrdinal>(nnz); colID++) {
              LocalOrdinal col = indices[colID];

              // we avoid a square root by using squared values
              typename STS::magnitudeType aiiajj = STS::magnitude(threshold*threshold * ghostedDiagVals[col]*ghostedDiagVals[row]);  // eps^2*|a_ii|*|a_jj|
              typename STS::magnitudeType aij    = STS::magnitude(vals[colID]*vals[colID]);                                          // |a_ij|^2

              if (aij > aiiajj || row == col) {
                columns[realnnz++] = col;
                rownnz++;
              } else
                numDropped++;
            }
            if (rownnz == 1) {
              // If the only element remaining after filtering is diagonal, mark node as bounday
              // FIXME: this should really be replaced by the following
              //    if (indices.size() == 1 && indices[0] == row)
              //        amalgBoundaryNodes[row] = true;
              // We do not do it this way now because there is no framework for distinguishing isolated
              // and boundary nodes in the aggregation algorithms
              amalgBoundaryNodes[row] = true;
            }
            rows[row+1] = realnnz;
          }
          columns.resize(realnnz);

          numTotal = A->getNodeNumEntries();

          RCP<GraphBase> graph = rcp(new LWGraph(rows, columns, A->getRowMap(), A->getColMap(), "amalgamated graph of A"));
          graph->SetBoundaryNodeMap(boundaryNodes);
          Set(currentLevel, "Graph",       graph);
          Set(currentLevel, "DofsPerNode", 1);

        } else if ( (A->GetFixedBlockSize() > 1) && (threshold == STS::zero()) ) {
          // Case 3:  Multiple DOF/node problem without dropping
          // TODO
          throw Exceptions::NotImplemented("Fast CoalesceDrop with multiple DOFs is not yet implemented.");

        } else if ( (A->GetFixedBlockSize() > 1) && (threshold != STS::zero()) ) {
          // Case 4:  Multiple DOF/node problem with dropping
          // TODO
          throw Exceptions::NotImplemented("Fast CoalesceDrop with multiple DOFs and dropping is not yet implemented.");
        }

      } else if (algo == "laplacian") {
        LO blkSize   = A->GetFixedBlockSize();
        GO indexBase = A->getRowMap()->getIndexBase();

        // [*0*] : FIXME
        // ap: somehow, if I move this line to [*1*], Belos throws an error
        // I'm not sure what's going on. Do we always have to Get data, if we did
        // DeclareInput for it?
        RCP<MultiVector> Coords = Get< RCP<MultiVector> >(currentLevel, "Coordinates");

        // Detect and record rows that correspond to Dirichlet boundary conditions
        // TODO If we use ArrayRCP<LO>, then we can record boundary nodes as usual.  Size
        // TODO the array one bigger than the number of local rows, and the last entry can
        // TODO hold the actual number of boundary nodes.  Clever, huh?
        const ArrayRCP<const bool > pointBoundaryNodes = MueLu::Utils<SC,LO,GO,NO>::DetectDirichletRows(*A, dirichletThreshold);

        if ( (blkSize == 1) && (threshold == STS::zero()) ) {
          // Trivial case: scalar problem, no dropping. Can return original graph
          RCP<GraphBase> graph = rcp(new Graph(A->getCrsGraph(), "graph of A"));
          graph->SetBoundaryNodeMap(pointBoundaryNodes);

          if (GetVerbLevel() & Statistics0) {
            GO numLocalBoundaryNodes  = 0;
            GO numGlobalBoundaryNodes = 0;
            for (LO i = 0; i < pointBoundaryNodes.size(); ++i)
              if (pointBoundaryNodes[i])
                numLocalBoundaryNodes++;
            RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
            sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
            GetOStream(Statistics0, 0) << "Detected " << numGlobalBoundaryNodes << " Dirichlet nodes" << std::endl;
          }

          Set(currentLevel, "DofsPerNode", blkSize);
          Set(currentLevel, "Graph",       graph);

        } else {
          // ap: We make quite a few assumptions here; general case may be a lot different,
          // but much much harder to implement. We assume that:
          //  1) all maps are standard maps, not strided maps
          //  2) global indices of dofs in A are related to dofs in coordinates in a simple arithmetic
          //     way: rows i*blkSize, i*blkSize+1, ..., i*blkSize + (blkSize-1) correspond to node i
          //
          // NOTE: Potentially, some of the code below could be simplified with UnAmalgamationInfo,
          // but as I totally don't understand that code, here is my solution

          // [*1*]: see [*0*]

          // Check that the number of local coordinates is consistent with the #rows in A
          std::string msg = "MueLu::CoalesceDropFactory::Build : coordinate vector length is incompatible with number of rows in A.  The vector length should be the same as the number of mesh points.";
          TEUCHOS_TEST_FOR_EXCEPTION(A->getRowMap()->getNodeNumElements()/blkSize != Coords->getLocalLength(), Exceptions::Incompatible, msg);

          const RCP<const Map> colMap = A->getColMap();
          RCP<const Map> uniqueMap, nonUniqueMap;
          if (blkSize == 1) {
            uniqueMap    = A->getRowMap();
            nonUniqueMap = A->getColMap();

          } else {
            uniqueMap    = Coords->getMap();
            TEUCHOS_TEST_FOR_EXCEPTION(uniqueMap->getIndexBase() != indexBase, Exceptions::Incompatible, "Different index bases for matrix and coordinates");

            // Amalgamate column map
            ArrayView<const GO> elementAList = colMap->getNodeElementList();
            size_t              numElements  = elementAList.size();
            Array<GO>           elementList(numElements);
            std::set<GO>        filter;        // TODO:  replace std::set with an object having faster lookup/insert, hashtable for instance

            LO numRows = 0;
            for (LO id = 0; id < static_cast<LO>(numElements); id++) {
              GO amalgID = (elementAList[id] - indexBase)/blkSize + indexBase;
              if (filter.find(amalgID) == filter.end()) {
                elementList[numRows++] = amalgID;
                filter.insert(amalgID);
              }
            }
            elementList.resize(numRows);

#ifdef HAVE_MUELU_DEBUG
            // At the moment we assume that first GIDs in nonUniqueMap coincide with those in uniqueMap
            for (LO row = 0; row < uniqueMap->getNodeNumElements(); row++)
              TEUCHOS_TEST_FOR_EXCEPTION(uniqueMap->getGlobalElement(row) != elementList[row], Exceptions::RuntimeError, "row = " << row << ", uniqueMap GID = " << uniqueMap->getGlobalElement(row) << ", nonUniqueMap GID = " << elementList[row] << std::endl);
#endif

            nonUniqueMap = MapFactory::Build(uniqueMap->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), elementList, indexBase, uniqueMap->getComm());
          }
          LocalOrdinal numRows = Teuchos::as<LocalOrdinal>(uniqueMap->getNodeNumElements());

          RCP<MultiVector>          ghostedCoords;
          RCP<Vector>               ghostedLaplDiag;
          Teuchos::ArrayRCP<Scalar> ghostedLaplDiagData;
          if (threshold != STS::zero()) {
            // Get ghost coordinates
            RCP<const Import> importer = ImportFactory::Build(uniqueMap, nonUniqueMap);
            ghostedCoords = MultiVectorFactory::Build(nonUniqueMap, Coords->getNumVectors());
            ghostedCoords->doImport(*Coords, *importer, Xpetra::INSERT);

            // Construct Distance Laplacian diagonal
            RCP<Vector>      localLaplDiag     = VectorFactory::Build(uniqueMap);
            ArrayRCP<Scalar> localLaplDiagData = localLaplDiag->getDataNonConst(0);
            for (LocalOrdinal row = 0; row < numRows; row++) {
              ArrayView<const LocalOrdinal> indices;
              Array<LocalOrdinal>           indicesExtra;

              if (blkSize == 1) {
                ArrayView<const Scalar> vals;
                A->getLocalRowView(row, indices, vals);

              } else {
                // Merge rows of A
                std::set<LO> cols;
                for (LO j = 0; j < blkSize; ++j) {
                  ArrayView<const LocalOrdinal> inds;
                  ArrayView<const Scalar>       vals;
                  A->getLocalRowView(row*blkSize+j, inds, vals);
                  for (LO k = 0; k < inds.size(); k++) {
                    // TODO: speed this up by using something like map for translation
                    LO  dofLID = inds[k];
                    GO  dofGID = colMap->getGlobalElement(dofLID);
                    GO nodeGID = (dofGID-indexBase)/blkSize + indexBase;
                    LO nodeLID = nonUniqueMap->getLocalElement(nodeGID);
                    cols.insert(nodeLID);
                  }
                }
                indicesExtra.resize(cols.size());
                size_t pos = 0;
                for (typename std::set<LocalOrdinal>::const_iterator it = cols.begin(); it != cols.end(); it++)
                  indicesExtra[pos++] = *it;
                indices = indicesExtra;
              }

              LocalOrdinal nnz = indices.size();
              for (LocalOrdinal colID = 0; colID < nnz; colID++) {
                LocalOrdinal col = indices[colID];

                if (row != col)
                  localLaplDiagData[row] += STS::one()/MueLu::Utils<SC,LO,GO,NO>::Distance2(*ghostedCoords, row, col);
              }
            }
            ghostedLaplDiag = VectorFactory::Build(nonUniqueMap);
            ghostedLaplDiag->doImport(*localLaplDiag, *importer, Xpetra::INSERT);
            ghostedLaplDiagData = ghostedLaplDiag->getDataNonConst(0);

          } else {
            GetOStream(Runtime0,0) << "Skipping distance laplacian construction due to 0 threshold" << std::endl;
          }

          // NOTE: ghostedLaplDiagData might be zero if we don't actually calculate the laplacian

          // allocate space for the local graph
          ArrayRCP<LocalOrdinal> rows    = ArrayRCP<LO>(numRows+1);
          ArrayRCP<LocalOrdinal> columns = ArrayRCP<LO>(A->getNodeNumEntries());

          const ArrayRCP<bool> amalgBoundaryNodes(numRows, false);

          LocalOrdinal realnnz = 0;
          rows[0] = 0;
          for (LocalOrdinal row = 0; row < numRows; row++) {
            ArrayView<const LocalOrdinal> indices;
            Array<LocalOrdinal>           indicesExtra;

            if (blkSize == 1) {
              ArrayView<const Scalar>     vals;
              A->getLocalRowView(row, indices, vals);

            } else {
              // FIXME: for now, if any of the rows in the row block is Dirichlet, we
              // assume that all are
              // This may not be true in general, for instance we might have mixed b.c.
              // where pressure is Dirichlet and velocities are not
              bool isBoundary = false;
              for (LO j = 0; j < blkSize; j++)
                if (pointBoundaryNodes[row*blkSize+j])
                  isBoundary = true;

              // Merge rows of A
              std::set<LO> cols;
              if (!isBoundary) {
                for (LO j = 0; j < blkSize; j++) {
                  ArrayView<const LocalOrdinal> inds;
                  ArrayView<const Scalar>       vals;
                  A->getLocalRowView(row*blkSize+j, inds, vals);
                  for (LO k = 0; k < inds.size(); k++) {
                    // TODO: speed this up by using something like map for translation
                    LO  dofLID = inds[k];
                    GO  dofGID = colMap->getGlobalElement(dofLID);
                    GO nodeGID = (dofGID-indexBase)/blkSize + indexBase;
                    LO nodeLID = nonUniqueMap->getLocalElement(nodeGID);
                    cols.insert(nodeLID);
                  }
                }
              } else {
                cols.insert(row);
              }
              indicesExtra.resize(cols.size());
              size_t pos = 0;
              for (typename std::set<LocalOrdinal>::const_iterator it = cols.begin(); it != cols.end(); it++)
                indicesExtra[pos++] = *it;
              indices = indicesExtra;
            }
            numTotal += indices.size();

            LocalOrdinal nnz = indices.size(), rownnz = 0;
            if (threshold != STS::zero()) {
              for (LocalOrdinal colID = 0; colID < nnz; colID++) {
                LocalOrdinal col = indices[colID];

                if (row == col) {
                  columns[realnnz++] = col;
                  rownnz++;
                  continue;
                }

                Scalar laplVal = STS::one() / MueLu::Utils<SC,LO,GO,NO>::Distance2(*ghostedCoords, row, col);
                typename STS::magnitudeType aiiajj = STS::magnitude(threshold*threshold * ghostedLaplDiagData[row]*ghostedLaplDiagData[col]);
                typename STS::magnitudeType aij    = STS::magnitude(laplVal*laplVal);

                if (aij > aiiajj) {
                  columns[realnnz++] = col;
                  rownnz++;
                } else {
                  numDropped++;
                }
              }

            } else {
              // Skip laplace calculation and threshold comparison for zero threshold
              for (LocalOrdinal colID = 0; colID < nnz; colID++) {
                LocalOrdinal col = indices[colID];
                columns[realnnz++] = col;
                rownnz++;
              }
            }

            if (rownnz == 1) {
              // If the only element remaining after filtering is diagonal, mark node as bounday
              // FIXME: this should really be replaced by the following
              //    if (indices.size() == 1 && indices[0] == row)
              //        boundaryNodes[row] = true;
              // We do not do it this way now because there is no framework for distinguishing isolated
              // and boundary nodes in the aggregation algorithms
              amalgBoundaryNodes[row] = true;
            }
            rows[row+1] = realnnz;
          } //for (LocalOrdinal row = 0; row < numRows; row++)
          columns.resize(realnnz);

          RCP<GraphBase> graph = rcp(new LWGraph(rows, columns, uniqueMap, nonUniqueMap, "amalgamated graph of A"));
          graph->SetBoundaryNodeMap(amalgBoundaryNodes);

          if (GetVerbLevel() & Statistics0) {
            GO numLocalBoundaryNodes  = 0;
            GO numGlobalBoundaryNodes = 0;

            for (LO i = 0; i < amalgBoundaryNodes.size(); ++i)
              if (amalgBoundaryNodes[i])
                numLocalBoundaryNodes++;

            RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
            sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
            GetOStream(Statistics0, 0) << "Detected " << numGlobalBoundaryNodes << " agglomerated Dirichlet nodes"
                                       << " using threshold " << dirichletThreshold << std::endl;
          }

          Set(currentLevel, "Graph",       graph);
          Set(currentLevel, "DofsPerNode", blkSize);
        } //if ( (blkSize == 1) && (threshold == STS::zero()) )
      }

      if (GetVerbLevel() & Statistics0) {
          RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
          GlobalOrdinal numGlobalTotal, numGlobalDropped;
          sumAll(comm, numTotal,   numGlobalTotal);
          sumAll(comm, numDropped, numGlobalDropped);
          GetOStream(Statistics0, -1) << "Number of dropped entries in amalgamated graph: " << numGlobalDropped << "/" << numGlobalTotal
              << " (" << 100*Teuchos::as<double>(numGlobalDropped)/Teuchos::as<double>(numGlobalTotal) << "%)" << std::endl;
      }

    } else {

      Set(currentLevel, "Filtering", (predrop_ != Teuchos::null));
      //what Tobias has implemented

      RCP<const Map> rowMap = A->getRowMap();
      RCP<const Map> colMap = A->getColMap();

      LO blockdim = 1;                          // block dim for fixed size blocks
      GO indexBase = rowMap->getIndexBase();    // index base of maps
      GO offset    = indexBase;                 // global offset of dof gids

      // 1) check for blocking/striding information
      if(A->IsView("stridedMaps") &&
         Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap("stridedMaps")) != Teuchos::null) {
        Xpetra::viewLabel_t oldView = A->SwitchToView("stridedMaps"); // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
        RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap());
        TEUCHOS_TEST_FOR_EXCEPTION(strMap == Teuchos::null,Exceptions::BadCast,"MueLu::CoalesceFactory::Build: cast to strided row map failed.");
        blockdim = strMap->getFixedBlockSize(); // TODO shorten code
        offset   = strMap->getOffset();
        oldView = A->SwitchToView(oldView);
        GetOStream(Statistics0, -1) << "CoalesceDropFactory::Build():" << " found blockdim=" << blockdim << " from strided maps. offset=" << offset << std::endl;
      } else GetOStream(Statistics0, -1) << "CoalesceDropFactory::Build(): no striding information available. Use blockdim=1 with offset=0" << std::endl;

      // 2) build (un)amalgamation information
      //    prepare generation of nodeRowMap (of amalgamated matrix)
      // TODO: special handling for blockdim=1
      RCP<AmalgamationInfo> amalInfo = Get< RCP<AmalgamationInfo> >(currentLevel, "UnAmalgamationInfo");
      RCP<std::map<GlobalOrdinal,std::vector<GlobalOrdinal> > > nodegid2dofgids = amalInfo->GetGlobalAmalgamationParams();
      RCP<std::vector<GlobalOrdinal> > gNodeIds = amalInfo->GetNodeGIDVector();
      GlobalOrdinal cnt_amalRows = amalInfo->GetNumberOfNodes();

      // inter processor communication: sum up number of block ids
      GlobalOrdinal num_blockids = 0;
      Teuchos::reduceAll<int,GlobalOrdinal>(*(A->getRowMap()->getComm()),Teuchos::REDUCE_SUM, cnt_amalRows, Teuchos::ptr(&num_blockids) );
      GetOStream(Statistics0, -1) << "CoalesceDropFactory::SetupAmalgamationData()" << " # of amalgamated blocks=" << num_blockids << std::endl;

      // 3) generate row map for amalgamated matrix (graph of A)
      //    with same distribution over all procs as row map of A
      Teuchos::ArrayRCP<GlobalOrdinal> arr_gNodeIds = Teuchos::arcp( gNodeIds );
      Teuchos::RCP<Map> nodeMap = MapFactory::Build(A->getRowMap()->lib(), num_blockids, arr_gNodeIds(), A->getRowMap()->getIndexBase(), A->getRowMap()->getComm());
      GetOStream(Statistics0, -1) << "CoalesceDropFactory: nodeMap " << nodeMap->getNodeNumElements() << "/" << nodeMap->getGlobalNumElements() << " elements" << std::endl;

      /////////////////////// experimental
      // vector of boundary node GIDs on current proc
      //RCP<std::map<GlobalOrdinal,bool> > gBoundaryNodes = Teuchos::rcp(new std::map<GlobalOrdinal,bool>);
      ////////////////////////////////////

      // 4) create graph of amalgamated matrix
      RCP<CrsGraph> crsGraph = CrsGraphFactory::Build(nodeMap, 10, Xpetra::DynamicProfile);

      // 5) do amalgamation. generate graph of amalgamated matrix
      for(LocalOrdinal row=0; row<Teuchos::as<LocalOrdinal>(A->getRowMap()->getNodeNumElements()); row++) {
        // get global DOF id
        GlobalOrdinal grid = rowMap->getGlobalElement(row);

        // translate grid to nodeid
        GlobalOrdinal nodeId = AmalgamationFactory::DOFGid2NodeId(grid, A, blockdim, offset);

        size_t nnz = A->getNumEntriesInLocalRow(row);
        Teuchos::ArrayView<const LocalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> vals;
        A->getLocalRowView(row, indices, vals);
        //TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Exceptions::RuntimeError, "MueLu::CoalesceFactory::Amalgamate: number of nonzeros not equal to number of indices? Error.");

        RCP<std::vector<GlobalOrdinal> > cnodeIds = Teuchos::rcp(new std::vector<GlobalOrdinal>);  // global column block ids
        LocalOrdinal realnnz = 0;
        for(LocalOrdinal col=0; col<Teuchos::as<LocalOrdinal>(nnz); col++) {
          //TEUCHOS_TEST_FOR_EXCEPTION(A->getColMap()->isNodeLocalElement(indices[col])==false,Exceptions::RuntimeError, "MueLu::CoalesceFactory::Amalgamate: Problem with columns. Error.");
          GlobalOrdinal gcid = colMap->getGlobalElement(indices[col]); // global column id

          if((predrop_ == Teuchos::null && vals[col]!=0.0) ||
             (predrop_ != Teuchos::null && predrop_->Drop(row,grid, col,indices[col],gcid,indices,vals) == false)) {
            GlobalOrdinal cnodeId = AmalgamationFactory::DOFGid2NodeId(gcid, A, blockdim, offset);
            cnodeIds->push_back(cnodeId);
            realnnz++; // increment number of nnz in matrix row
          }
        }

        // todo avoid duplicate entries in cnodeIds

        ////////////////// experimental
        //if(gBoundaryNodes->count(nodeId) == 0)
        //  (*gBoundaryNodes)[nodeId] = false;  // new node GID (probably no Dirichlet bdry node)
        //if(realnnz == 1)
        //  (*gBoundaryNodes)[nodeId] = true; // if there's only one nnz entry the node has some Dirichlet bdry dofs
        ///////////////////////////////

        Teuchos::ArrayRCP<GlobalOrdinal> arr_cnodeIds = Teuchos::arcp( cnodeIds );

        //TEUCHOS_TEST_FOR_EXCEPTION(crsGraph->getRowMap()->isNodeGlobalElement(nodeId)==false,Exceptions::RuntimeError, "MueLu::CoalesceFactory::Amalgamate: global row id does not belong to current proc. Error.");
        if(arr_cnodeIds.size() > 0 )
          crsGraph->insertGlobalIndices(nodeId, arr_cnodeIds());
      }
      // fill matrix graph
      crsGraph->fillComplete(nodeMap,nodeMap);

      ///////////////// experimental
      //LocalOrdinal nLocalBdryNodes = 0;
      //GlobalOrdinal nGlobalBdryNodes = 0;
      //Array<GlobalOrdinal> bdryNodeIds;
      //for(typename std::map<GlobalOrdinal,bool>::iterator it = gBoundaryNodes->begin(); it!=gBoundaryNodes->end(); it++) {
      //  if ((*it).second == true) {
      //    nLocalBdryNodes++;
      //    bdryNodeIds.push_back((*it).first);
      //  }
      //}
      //Teuchos::reduceAll<int,GlobalOrdinal>(*(A->getRowMap()->getComm()),Teuchos::REDUCE_SUM, nLocalBdryNodes, Teuchos::ptr(&nGlobalBdryNodes) );
      //GetOStream(Debug, 0) << "CoalesceDropFactory::SetupAmalgamationData()" << " # detected Dirichlet boundary nodes = " << nGlobalBdryNodes << std::endl;

      //RCP<const Map> gBoundaryNodeMap = MapFactory::Build(nodeMap->lib(),
      //                                                    Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
      //                                                    bdryNodeIds,
      //                                                    nodeMap->getIndexBase(), nodeMap->getComm());
      //////////////////////////////

      // 6) create MueLu Graph object
      RCP<GraphBase> graph = rcp(new Graph(crsGraph, "amalgamated graph of A"));

      // 7) store results in Level
      //graph->SetBoundaryNodeMap(gBoundaryNodeMap);
      Set(currentLevel, "DofsPerNode", blockdim);
      Set(currentLevel, "Graph", graph);

    } //if (doExperimentalWrap) ... else ...

  } //Build

} //namespace MueLu

#endif // MUELU_COALESCEDROPFACTORY_DEF_HPP
