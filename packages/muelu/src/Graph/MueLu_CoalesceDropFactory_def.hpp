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
#ifndef MUELU_COALESCEDROPFACTORY_DEF_HPP
#define MUELU_COALESCEDROPFACTORY_DEF_HPP

#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_CrsGraph.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_StridedMap.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>

#include "MueLu_CoalesceDropFactory_decl.hpp"

#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_LWGraph.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PreDropFunctionBaseClass.hpp"
#include "MueLu_PreDropFunctionConstVal.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

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

    validParamList->set< RCP<const FactoryBase> >("A",                  Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Generating factory for UnAmalgamationInfo");
    validParamList->set< RCP<const FactoryBase> >("Coordinates",        Teuchos::null, "Generating factory for Coordinates");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CoalesceDropFactory() : predrop_(Teuchos::null) { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
    Input(currentLevel, "A");
    Input(currentLevel, "UnAmalgamationInfo");

    const ParameterList& pL = GetParameterList();
    if (pL.get<bool>("lightweight wrap") == true) {
      if (pL.get<std::string>("aggregation: drop scheme") == "distance laplacian")
        Input(currentLevel, "Coordinates");

    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    typedef Xpetra::MultiVector<double,LO,GO,NO> dxMV;

    typedef Teuchos::ScalarTraits<SC> STS;

    if (predrop_ != Teuchos::null)
      GetOStream(Parameters0) << predrop_->description();

    RCP<Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");
    RCP<AmalgamationInfo> amalInfo = Get< RCP<AmalgamationInfo> >(currentLevel, "UnAmalgamationInfo");

    const ParameterList  & pL = GetParameterList();
    bool doExperimentalWrap = pL.get<bool>("lightweight wrap");

    GetOStream(Parameters0) << "lightweight wrap = " << doExperimentalWrap << std::endl;

    // decide wether to use the fast-track code path for standard maps or the somewhat slower
    // code path for non-standard maps
    /*bool bNonStandardMaps = false;
    if (A->IsView("stridedMaps") == true) {
      Teuchos::RCP<const Map> myMap = A->getRowMap("stridedMaps");
      Teuchos::RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(myMap);
      TEUCHOS_TEST_FOR_EXCEPTION(strMap == null, Exceptions::RuntimeError, "Map is not of type StridedMap");
      if (strMap->getStridedBlockId() != -1 || strMap->getOffset() > 0)
        bNonStandardMaps = true;
    }*/

    if (doExperimentalWrap) {
      std::string algo = pL.get<std::string>("aggregation: drop scheme");

      TEUCHOS_TEST_FOR_EXCEPTION(predrop_ != null   && algo != "classical", Exceptions::RuntimeError, "Dropping function must not be provided for \"" << algo << "\" algorithm");
      TEUCHOS_TEST_FOR_EXCEPTION(algo != "classical" && algo != "distance laplacian", Exceptions::RuntimeError, "\"algorithm\" must be one of (classical|distance laplacian)");

      SC threshold = as<SC>(pL.get<double>("aggregation: drop tol"));
      GetOStream(Runtime0) << "algorithm = \"" << algo << "\": threshold = " << threshold << ", blocksize = " << A->GetFixedBlockSize() << std::endl;
      Set(currentLevel, "Filtering", (threshold != STS::zero()));

      const typename STS::magnitudeType dirichletThreshold = STS::magnitude(as<SC>(pL.get<double>("aggregation: Dirichlet threshold")));

      GO numDropped = 0, numTotal = 0;
      std::string graphType = "unamalgamated"; //for description purposes only
      if (algo == "classical") {
        if (predrop_ == null) {
          // ap: this is a hack: had to declare predrop_ as mutable
          predrop_ = rcp(new PreDropFunctionConstVal(threshold));
        }

        if (predrop_ != null) {
          RCP<PreDropFunctionConstVal> predropConstVal = rcp_dynamic_cast<PreDropFunctionConstVal>(predrop_);
          TEUCHOS_TEST_FOR_EXCEPTION(predropConstVal == Teuchos::null, Exceptions::BadCast,
                                     "MueLu::CoalesceFactory::Build: cast to PreDropFunctionConstVal failed.");
          // If a user provided a predrop function, it overwrites the XML threshold parameter
          SC newt = predropConstVal->GetThreshold();
          if (newt != threshold) {
            GetOStream(Warnings0) << "switching threshold parameter from " << threshold << " (list) to " << newt << " (user function" << std::endl;
            threshold = newt;
          }
        }
        // At this points we either have
        //     (predrop_ != null)
        // Therefore, it is sufficient to check only threshold
        if (A->GetFixedBlockSize() == 1 && threshold == STS::zero() && A->hasCrsGraph()) {
          // Case 1:  scalar problem, no dropping => just use matrix graph
          RCP<GraphBase> graph = rcp(new Graph(A->getCrsGraph(), "graph of A"));
          // Detect and record rows that correspond to Dirichlet boundary conditions
          ArrayRCP<const bool > boundaryNodes;
          boundaryNodes = MueLu::Utilities<SC,LO,GO,NO>::DetectDirichletRows(*A, dirichletThreshold);
          graph->SetBoundaryNodeMap(boundaryNodes);
          numTotal = A->getNodeNumEntries();

          if (GetVerbLevel() & Statistics1) {
            GO numLocalBoundaryNodes  = 0;
            GO numGlobalBoundaryNodes = 0;
            for (LO i = 0; i < boundaryNodes.size(); ++i)
              if (boundaryNodes[i])
                numLocalBoundaryNodes++;
            RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
            MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
            GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes << " Dirichlet nodes" << std::endl;
          }

          Set(currentLevel, "DofsPerNode", 1);
          Set(currentLevel, "Graph", graph);

        } else if ( (A->GetFixedBlockSize() == 1 && threshold != STS::zero()) || 
                    (A->GetFixedBlockSize() == 1 && threshold == STS::zero() && !A->hasCrsGraph())) {
          // Case 2:  scalar problem with dropping => record the column indices of undropped entries, but still use original
          //                                          graph's map information, e.g., whether index is local
          // OR a matrix without a CrsGraph

          // allocate space for the local graph
          ArrayRCP<LO> rows   (A->getNodeNumRows()+1);
          ArrayRCP<LO> columns(A->getNodeNumEntries());

          RCP<Vector> ghostedDiag = MueLu::Utilities<SC,LO,GO,NO>::GetMatrixOverlappedDiagonal(*A);
          const ArrayRCP<const SC> ghostedDiagVals = ghostedDiag->getData(0);
          const ArrayRCP<bool>     boundaryNodes(A->getNodeNumRows(), false);

          LO realnnz = 0;

          rows[0] = 0;
          for (LO row = 0; row < Teuchos::as<LO>(A->getRowMap()->getNodeNumElements()); ++row) {
            size_t nnz = A->getNumEntriesInLocalRow(row);
            ArrayView<const LO> indices;
            ArrayView<const SC> vals;
            A->getLocalRowView(row, indices, vals);

            //FIXME the current predrop function uses the following
            //FIXME    if(std::abs(vals[k]) > std::abs(threshold_) || grow == gcid )
            //FIXME but the threshold doesn't take into account the rows' diagonal entries
            //FIXME For now, hardwiring the dropping in here

            LO rownnz = 0;
            for (LO colID = 0; colID < Teuchos::as<LO>(nnz); colID++) {
              LO col = indices[colID];

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
              // If the only element remaining after filtering is diagonal, mark node as boundary
              // FIXME: this should really be replaced by the following
              //    if (indices.size() == 1 && indices[0] == row)
              //        boundaryNodes[row] = true;
              // We do not do it this way now because there is no framework for distinguishing isolated
              // and boundary nodes in the aggregation algorithms
              boundaryNodes[row] = true;
            }
            rows[row+1] = realnnz;
          }
          columns.resize(realnnz);
          numTotal = A->getNodeNumEntries();

          RCP<GraphBase> graph = rcp(new LWGraph(rows, columns, A->getRowMap(), A->getColMap(), "thresholded graph of A"));
          graph->SetBoundaryNodeMap(boundaryNodes);
          if (GetVerbLevel() & Statistics1) {
            GO numLocalBoundaryNodes  = 0;
            GO numGlobalBoundaryNodes = 0;
            for (LO i = 0; i < boundaryNodes.size(); ++i)
              if (boundaryNodes[i])
                numLocalBoundaryNodes++;
            RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
            MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
            GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes << " Dirichlet nodes" << std::endl;
          }
          Set(currentLevel, "Graph",       graph);
          Set(currentLevel, "DofsPerNode", 1);

        } else if (A->GetFixedBlockSize() > 1 && threshold == STS::zero()) {
          // Case 3:  Multiple DOF/node problem without dropping
          const RCP<const Map> rowMap = A->getRowMap();
          const RCP<const Map> colMap = A->getColMap();

          graphType = "amalgamated";

          // build node row map (uniqueMap) and node column map (nonUniqueMap)
          // the arrays rowTranslation and colTranslation contain the local node id
          // given a local dof id. The data is calculated by the AmalgamationFactory and
          // stored in the variable container "UnAmalgamationInfo"
          RCP<const Map> uniqueMap = amalInfo->getNodeRowMap();
          RCP<const Map> nonUniqueMap = amalInfo->getNodeColMap();
          Array<LO> rowTranslation = *(amalInfo->getRowTranslation());
          Array<LO> colTranslation = *(amalInfo->getColTranslation());

          // get number of local nodes
          LO numRows = Teuchos::as<LocalOrdinal>(uniqueMap->getNodeNumElements());

          // Allocate space for the local graph
          ArrayRCP<LO> rows    = ArrayRCP<LO>(numRows+1);
          ArrayRCP<LO> columns = ArrayRCP<LO>(A->getNodeNumEntries());

          const ArrayRCP<bool> amalgBoundaryNodes(numRows, false);

          // Detect and record rows that correspond to Dirichlet boundary conditions
          // TODO If we use ArrayRCP<LO>, then we can record boundary nodes as usual.  Size
          // TODO the array one bigger than the number of local rows, and the last entry can
          // TODO hold the actual number of boundary nodes.  Clever, huh?
          ArrayRCP<const bool > pointBoundaryNodes;
          pointBoundaryNodes = MueLu::Utilities<SC,LO,GO,NO>::DetectDirichletRows(*A, dirichletThreshold);

          // extract striding information
          LO blkSize = A->GetFixedBlockSize();     //< the full block size (number of dofs per node in strided map)
          LO blkId   = -1;                         //< the block id within the strided map (or -1 if it is a full block map)
          LO blkPartSize = A->GetFixedBlockSize(); //< stores the size of the block within the strided map
          if (A->IsView("stridedMaps") == true) {
            Teuchos::RCP<const Map> myMap = A->getRowMap("stridedMaps");
            Teuchos::RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(myMap);
            TEUCHOS_TEST_FOR_EXCEPTION(strMap == null, Exceptions::RuntimeError, "Map is not of type StridedMap");
            blkSize = Teuchos::as<const LO>(strMap->getFixedBlockSize());
            blkId   = strMap->getStridedBlockId();
            if (blkId > -1)
              blkPartSize = Teuchos::as<LO>(strMap->getStridingData()[blkId]);
          }

          // loop over all local nodes
          LO realnnz = 0;
          rows[0] = 0;
          Array<LO> indicesExtra;
          for (LO row = 0; row < numRows; row++) {
            ArrayView<const LO> indices;
            indicesExtra.resize(0);

            // The amalgamated row is marked as Dirichlet iff all point rows are Dirichlet
            // Note, that pointBoundaryNodes lives on the dofmap (and not the node map).
            // Therefore, looping over all dofs is fine here. We use blkPartSize as we work
            // with local ids.
            // TODO: Here we have different options of how to define a node to be a boundary (or Dirichlet)
            // node.
            bool isBoundary = false;
            isBoundary = true;
            for (LO j = 0; j < blkPartSize; j++) {
              if (!pointBoundaryNodes[row*blkPartSize+j]) {
                isBoundary = false;
                break;
              }
            }

            // Merge rows of A
            // The array indicesExtra contains local column node ids for the current local node "row"
            if (!isBoundary)
              MergeRows(*A, row, indicesExtra, colTranslation);
            else
              indicesExtra.push_back(row);
            indices   = indicesExtra;
            numTotal += indices.size();

            // add the local column node ids to the full columns array which
            // contains the local column node ids for all local node rows
            LO nnz = indices.size(), rownnz = 0;
            for (LO colID = 0; colID < nnz; colID++) {
              LO col = indices[colID];
              columns[realnnz++] = col;
              rownnz++;
            }

            if (rownnz == 1) {
              // If the only element remaining after filtering is diagonal, mark node as boundary
              // FIXME: this should really be replaced by the following
              //    if (indices.size() == 1 && indices[0] == row)
              //        boundaryNodes[row] = true;
              // We do not do it this way now because there is no framework for distinguishing isolated
              // and boundary nodes in the aggregation algorithms
              amalgBoundaryNodes[row] = true;
            }
            rows[row+1] = realnnz;
          } //for (LO row = 0; row < numRows; row++)
          columns.resize(realnnz);

          RCP<GraphBase> graph = rcp(new LWGraph(rows, columns, uniqueMap, nonUniqueMap, "amalgamated graph of A"));
          graph->SetBoundaryNodeMap(amalgBoundaryNodes);

          if (GetVerbLevel() & Statistics1) {
            GO numLocalBoundaryNodes  = 0;
            GO numGlobalBoundaryNodes = 0;

            for (LO i = 0; i < amalgBoundaryNodes.size(); ++i)
              if (amalgBoundaryNodes[i])
                numLocalBoundaryNodes++;

            RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
            MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
            GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes
                                       << " agglomerated Dirichlet nodes" << std::endl;
          }

          Set(currentLevel, "Graph",       graph);
          Set(currentLevel, "DofsPerNode", blkSize); // full block size

        } else if (A->GetFixedBlockSize() > 1 && threshold != STS::zero()) {
          // Case 4:  Multiple DOF/node problem with dropping
          const RCP<const Map> rowMap = A->getRowMap();
          const RCP<const Map> colMap = A->getColMap();

          graphType = "amalgamated";

          // build node row map (uniqueMap) and node column map (nonUniqueMap)
          // the arrays rowTranslation and colTranslation contain the local node id
          // given a local dof id. The data is calculated by the AmalgamationFactory and
          // stored in the variable container "UnAmalgamationInfo"
          RCP<const Map> uniqueMap = amalInfo->getNodeRowMap();
          RCP<const Map> nonUniqueMap = amalInfo->getNodeColMap();
          Array<LO> rowTranslation = *(amalInfo->getRowTranslation());
          Array<LO> colTranslation = *(amalInfo->getColTranslation());

          // get number of local nodes
          LO numRows = Teuchos::as<LocalOrdinal>(uniqueMap->getNodeNumElements());

          // Allocate space for the local graph
          ArrayRCP<LO> rows    = ArrayRCP<LO>(numRows+1);
          ArrayRCP<LO> columns = ArrayRCP<LO>(A->getNodeNumEntries());

          const ArrayRCP<bool> amalgBoundaryNodes(numRows, false);

          // Detect and record rows that correspond to Dirichlet boundary conditions
          // TODO If we use ArrayRCP<LO>, then we can record boundary nodes as usual.  Size
          // TODO the array one bigger than the number of local rows, and the last entry can
          // TODO hold the actual number of boundary nodes.  Clever, huh?
          ArrayRCP<const bool > pointBoundaryNodes;
          pointBoundaryNodes = MueLu::Utilities<SC,LO,GO,NO>::DetectDirichletRows(*A, dirichletThreshold);

          // extract striding information
          LO blkSize = A->GetFixedBlockSize();     //< the full block size (number of dofs per node in strided map)
          LO blkId   = -1;                         //< the block id within the strided map (or -1 if it is a full block map)
          LO blkPartSize = A->GetFixedBlockSize(); //< stores the size of the block within the strided map
          if (A->IsView("stridedMaps") == true) {
            Teuchos::RCP<const Map> myMap = A->getRowMap("stridedMaps");
            Teuchos::RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(myMap);
            TEUCHOS_TEST_FOR_EXCEPTION(strMap == null, Exceptions::RuntimeError, "Map is not of type StridedMap");
            blkSize = Teuchos::as<const LO>(strMap->getFixedBlockSize());
            blkId   = strMap->getStridedBlockId();
            if (blkId > -1)
              blkPartSize = Teuchos::as<LO>(strMap->getStridingData()[blkId]);
          }

          // extract diagonal data for dropping strategy
          RCP<Vector> ghostedDiag = MueLu::Utilities<SC,LO,GO,NO>::GetMatrixOverlappedDiagonal(*A);
          const ArrayRCP<const SC> ghostedDiagVals = ghostedDiag->getData(0);

          // loop over all local nodes
          LO realnnz = 0;
          rows[0] = 0;
          Array<LO> indicesExtra;
          for (LO row = 0; row < numRows; row++) {
            ArrayView<const LO> indices;
            indicesExtra.resize(0);

            // The amalgamated row is marked as Dirichlet iff all point rows are Dirichlet
            // Note, that pointBoundaryNodes lives on the dofmap (and not the node map).
            // Therefore, looping over all dofs is fine here. We use blkPartSize as we work
            // with local ids.
            // TODO: Here we have different options of how to define a node to be a boundary (or Dirichlet)
            // node.
            bool isBoundary = false;
            isBoundary = true;
            for (LO j = 0; j < blkPartSize; j++) {
              if (!pointBoundaryNodes[row*blkPartSize+j]) {
                isBoundary = false;
                break;
              }
            }

            // Merge rows of A
            // The array indicesExtra contains local column node ids for the current local node "row"
            if (!isBoundary)
              MergeRowsWithDropping(*A, row, ghostedDiagVals, threshold, indicesExtra, colTranslation);
            else
              indicesExtra.push_back(row);
            indices   = indicesExtra;
            numTotal += indices.size();

            // add the local column node ids to the full columns array which
            // contains the local column node ids for all local node rows
            LO nnz = indices.size(), rownnz = 0;
            for (LO colID = 0; colID < nnz; colID++) {
              LO col = indices[colID];
              columns[realnnz++] = col;
              rownnz++;
            }

            if (rownnz == 1) {
              // If the only element remaining after filtering is diagonal, mark node as boundary
              // FIXME: this should really be replaced by the following
              //    if (indices.size() == 1 && indices[0] == row)
              //        boundaryNodes[row] = true;
              // We do not do it this way now because there is no framework for distinguishing isolated
              // and boundary nodes in the aggregation algorithms
              amalgBoundaryNodes[row] = true;
            }
            rows[row+1] = realnnz;
          } //for (LO row = 0; row < numRows; row++)
          columns.resize(realnnz);

          RCP<GraphBase> graph = rcp(new LWGraph(rows, columns, uniqueMap, nonUniqueMap, "amalgamated graph of A"));
          graph->SetBoundaryNodeMap(amalgBoundaryNodes);

          if (GetVerbLevel() & Statistics1) {
            GO numLocalBoundaryNodes  = 0;
            GO numGlobalBoundaryNodes = 0;

            for (LO i = 0; i < amalgBoundaryNodes.size(); ++i)
              if (amalgBoundaryNodes[i])
                numLocalBoundaryNodes++;

            RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
            MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
            GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes
                                       << " agglomerated Dirichlet nodes" << std::endl;
          }

          Set(currentLevel, "Graph",       graph);
          Set(currentLevel, "DofsPerNode", blkSize); // full block size
        }

      } else if (algo == "distance laplacian") {
        LO blkSize   = A->GetFixedBlockSize();
        GO indexBase = A->getRowMap()->getIndexBase();

        // [*0*] : FIXME
        // ap: somehow, if I move this line to [*1*], Belos throws an error
        // I'm not sure what's going on. Do we always have to Get data, if we did
        // DeclareInput for it?
        RCP<dxMV> Coords = Get< RCP<Xpetra::MultiVector<double,LO,GO,NO> > >(currentLevel, "Coordinates");

        // Detect and record rows that correspond to Dirichlet boundary conditions
        // TODO If we use ArrayRCP<LO>, then we can record boundary nodes as usual.  Size
        // TODO the array one bigger than the number of local rows, and the last entry can
        // TODO hold the actual number of boundary nodes.  Clever, huh?
        ArrayRCP<const bool > pointBoundaryNodes;
        pointBoundaryNodes = MueLu::Utilities<SC,LO,GO,NO>::DetectDirichletRows(*A, dirichletThreshold);

        if ( (blkSize == 1) && (threshold == STS::zero()) ) {
          // Trivial case: scalar problem, no dropping. Can return original graph
          RCP<GraphBase> graph = rcp(new Graph(A->getCrsGraph(), "graph of A"));
          graph->SetBoundaryNodeMap(pointBoundaryNodes);
          graphType="unamalgamated";
          numTotal = A->getNodeNumEntries();

          if (GetVerbLevel() & Statistics1) {
            GO numLocalBoundaryNodes  = 0;
            GO numGlobalBoundaryNodes = 0;
            for (LO i = 0; i < pointBoundaryNodes.size(); ++i)
              if (pointBoundaryNodes[i])
                numLocalBoundaryNodes++;
            RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
            MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
            GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes << " Dirichlet nodes" << std::endl;
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
          TEUCHOS_TEST_FOR_EXCEPTION(A->getRowMap()->getNodeNumElements()/blkSize != Coords->getLocalLength(), Exceptions::Incompatible,
                                     "Coordinate vector length (" << Coords->getLocalLength() << ") is incompatible with number of rows in A (" << A->getRowMap()->getNodeNumElements() << ") by modulo block size ("<< blkSize <<").");

          const RCP<const Map> colMap = A->getColMap();
          RCP<const Map> uniqueMap, nonUniqueMap;
          Array<LO>      colTranslation;
          if (blkSize == 1) {
            uniqueMap    = A->getRowMap();
            nonUniqueMap = A->getColMap();
            graphType="unamalgamated";

          } else {
            uniqueMap    = Coords->getMap();
            TEUCHOS_TEST_FOR_EXCEPTION(uniqueMap->getIndexBase() != indexBase, Exceptions::Incompatible,
                                       "Different index bases for matrix and coordinates");

            AmalgamationFactory::AmalgamateMap(*(A->getColMap()), *A, nonUniqueMap, colTranslation);

            graphType = "amalgamated";
          }
          LO numRows = Teuchos::as<LocalOrdinal>(uniqueMap->getNodeNumElements());

          RCP<dxMV>             ghostedCoords;
          RCP<Vector>           ghostedLaplDiag;
          Teuchos::ArrayRCP<SC> ghostedLaplDiagData;
          if (threshold != STS::zero()) {
            // Get ghost coordinates
            RCP<const Import> importer;
            {
              SubFactoryMonitor m1(*this, "Import construction", currentLevel);
              if (blkSize == 1 && A->getCrsGraph()->getImporter() != Teuchos::null) {
                GetOStream(Warnings1) << "Using existing importer from matrix graph" << std::endl;
                importer = A->getCrsGraph()->getImporter();
              } else {
                GetOStream(Warnings0) << "Constructing new importer instance" << std::endl;
                importer = ImportFactory::Build(uniqueMap, nonUniqueMap);
              }
            } //subtimer
            ghostedCoords = Xpetra::MultiVectorFactory<double,LO,GO,NO>::Build(nonUniqueMap, Coords->getNumVectors());
            {
            SubFactoryMonitor m1(*this, "Coordinate import", currentLevel);
            ghostedCoords->doImport(*Coords, *importer, Xpetra::INSERT);
            } //subtimer

            // Construct Distance Laplacian diagonal
            RCP<Vector>  localLaplDiag     = VectorFactory::Build(uniqueMap);
            ArrayRCP<SC> localLaplDiagData = localLaplDiag->getDataNonConst(0);
            Array<LO> indicesExtra;
            Teuchos::Array<Teuchos::ArrayRCP<const double>> coordData;
            if (threshold != STS::zero()) {
              const size_t numVectors = ghostedCoords->getNumVectors();
              coordData.reserve(numVectors);
              for (size_t j = 0; j < numVectors; j++) {
                Teuchos::ArrayRCP<const double> tmpData=ghostedCoords->getData(j);
                coordData.push_back(tmpData);
              }
            }
            {
            SubFactoryMonitor m1(*this, "Laplacian local diagonal", currentLevel);
            for (LO row = 0; row < numRows; row++) {
              ArrayView<const LO> indices;

              if (blkSize == 1) {
                ArrayView<const SC> vals;
                A->getLocalRowView(row, indices, vals);

              } else {
                // Merge rows of A
                indicesExtra.resize(0);
                MergeRows(*A, row, indicesExtra, colTranslation);
                indices = indicesExtra;
              }

              LO nnz = indices.size();
              for (LO colID = 0; colID < nnz; colID++) {
                const LO col = indices[colID];

                if (row != col) {
                  localLaplDiagData[row] += STS::one()/MueLu::Utilities<double,LO,GO,NO>::Distance2(coordData, row, col);
                }
              }
            }
            } //subtimer
            {
            SubFactoryMonitor m1(*this, "Laplacian distributed diagonal", currentLevel);
            ghostedLaplDiag = VectorFactory::Build(nonUniqueMap);
            ghostedLaplDiag->doImport(*localLaplDiag, *importer, Xpetra::INSERT);
            ghostedLaplDiagData = ghostedLaplDiag->getDataNonConst(0);
            } //subtimer

          } else {
            GetOStream(Runtime0) << "Skipping distance laplacian construction due to 0 threshold" << std::endl;
          }

          // NOTE: ghostedLaplDiagData might be zero if we don't actually calculate the laplacian

          // allocate space for the local graph
          ArrayRCP<LO> rows    = ArrayRCP<LO>(numRows+1);
          ArrayRCP<LO> columns = ArrayRCP<LO>(A->getNodeNumEntries());

          const ArrayRCP<bool> amalgBoundaryNodes(numRows, false);

          LO realnnz = 0;
          rows[0] = 0;
          Array<LO> indicesExtra;
          {
          SubFactoryMonitor m1(*this, "Laplacian dropping", currentLevel);
          Teuchos::Array<Teuchos::ArrayRCP<const double>> coordData;
          if (threshold != STS::zero()) {
            const size_t numVectors = ghostedCoords->getNumVectors();
            coordData.reserve(numVectors);
            for (size_t j = 0; j < numVectors; j++) {
              Teuchos::ArrayRCP<const double> tmpData=ghostedCoords->getData(j);
              coordData.push_back(tmpData);
            }
          }
          for (LO row = 0; row < numRows; row++) {
            ArrayView<const LO> indices;
            indicesExtra.resize(0);

            if (blkSize == 1) {
              ArrayView<const SC>     vals;
              A->getLocalRowView(row, indices, vals);

            } else {
              // The amalgamated row is marked as Dirichlet iff all point rows are Dirichlet
              bool isBoundary = false;
              isBoundary = true;
              for (LO j = 0; j < blkSize; j++) {
                if (!pointBoundaryNodes[row*blkSize+j]) {
                  isBoundary = false;
                  break;
                }
              }

              // Merge rows of A
              if (!isBoundary)
                MergeRows(*A, row, indicesExtra, colTranslation);
              else
                indicesExtra.push_back(row);
              indices = indicesExtra;
            }
            numTotal += indices.size();

            LO nnz = indices.size(), rownnz = 0;
            if (threshold != STS::zero()) {
              for (LO colID = 0; colID < nnz; colID++) {
                LO col = indices[colID];

                if (row == col) {
                  columns[realnnz++] = col;
                  rownnz++;
                  continue;
                }

                SC laplVal = STS::one() / MueLu::Utilities<double,LO,GO,NO>::Distance2(coordData, row, col);
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
              for (LO colID = 0; colID < nnz; colID++) {
                LO col = indices[colID];
                columns[realnnz++] = col;
                rownnz++;
              }
            }

            if (rownnz == 1) {
              // If the only element remaining after filtering is diagonal, mark node as boundary
              // FIXME: this should really be replaced by the following
              //    if (indices.size() == 1 && indices[0] == row)
              //        boundaryNodes[row] = true;
              // We do not do it this way now because there is no framework for distinguishing isolated
              // and boundary nodes in the aggregation algorithms
              amalgBoundaryNodes[row] = true;
            }
            rows[row+1] = realnnz;
          } //for (LO row = 0; row < numRows; row++)
          } //subtimer
          columns.resize(realnnz);

          RCP<GraphBase> graph;
          {
          SubFactoryMonitor m1(*this, "Build amalgamated graph", currentLevel);
          graph = rcp(new LWGraph(rows, columns, uniqueMap, nonUniqueMap, "amalgamated graph of A"));
          graph->SetBoundaryNodeMap(amalgBoundaryNodes);
          } //subtimer

          if (GetVerbLevel() & Statistics1) {
            GO numLocalBoundaryNodes  = 0;
            GO numGlobalBoundaryNodes = 0;

            for (LO i = 0; i < amalgBoundaryNodes.size(); ++i)
              if (amalgBoundaryNodes[i])
                numLocalBoundaryNodes++;

            RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
            MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
            GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes << " agglomerated Dirichlet nodes"
                                       << " using threshold " << dirichletThreshold << std::endl;
          }

          Set(currentLevel, "Graph",       graph);
          Set(currentLevel, "DofsPerNode", blkSize);
        }
      }

      if ((GetVerbLevel() & Statistics1) && !(A->GetFixedBlockSize() > 1 && threshold != STS::zero())) {
          RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
          GO numGlobalTotal, numGlobalDropped;
          MueLu_sumAll(comm, numTotal,   numGlobalTotal);
          MueLu_sumAll(comm, numDropped, numGlobalDropped);
          GetOStream(Statistics1) << "Number of dropped entries in " << graphType << " matrix graph: " << numGlobalDropped << "/" << numGlobalTotal;
          if (numGlobalTotal != 0)
            GetOStream(Statistics1) << " (" << 100*Teuchos::as<double>(numGlobalDropped)/Teuchos::as<double>(numGlobalTotal) << "%)";
          GetOStream(Statistics1) << std::endl;
      }

    } else {
      //what Tobias has implemented

      SC threshold = as<SC>(pL.get<double>("aggregation: drop tol"));
      //GetOStream(Runtime0) << "algorithm = \"" << algo << "\": threshold = " << threshold << ", blocksize = " << A->GetFixedBlockSize() << std::endl;
      GetOStream(Runtime0) << "algorithm = \"" << "failsafe" << "\": threshold = " << threshold << ", blocksize = " << A->GetFixedBlockSize() << std::endl;
      Set(currentLevel, "Filtering", (threshold != STS::zero()));

      RCP<const Map> rowMap = A->getRowMap();
      RCP<const Map> colMap = A->getColMap();

      LO blockdim = 1;                          // block dim for fixed size blocks
      GO indexBase = rowMap->getIndexBase();    // index base of maps
      GO offset    = 0;

      // 1) check for blocking/striding information
      if(A->IsView("stridedMaps") &&
         Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap("stridedMaps")) != Teuchos::null) {
        Xpetra::viewLabel_t oldView = A->SwitchToView("stridedMaps"); // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
        RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap());
        TEUCHOS_TEST_FOR_EXCEPTION(strMap == Teuchos::null,Exceptions::BadCast,"MueLu::CoalesceFactory::Build: cast to strided row map failed.");
        blockdim = strMap->getFixedBlockSize();
        offset   = strMap->getOffset();
        oldView = A->SwitchToView(oldView);
        GetOStream(Statistics1) << "CoalesceDropFactory::Build():" << " found blockdim=" << blockdim << " from strided maps. offset=" << offset << std::endl;
      } else GetOStream(Statistics1) << "CoalesceDropFactory::Build(): no striding information available. Use blockdim=1 with offset=0" << std::endl;

      // 2) get row map for amalgamated matrix (graph of A)
      //    with same distribution over all procs as row map of A
      RCP<const Map> nodeMap = amalInfo->getNodeRowMap();
      GetOStream(Statistics1) << "CoalesceDropFactory: nodeMap " << nodeMap->getNodeNumElements() << "/" << nodeMap->getGlobalNumElements() << " elements" << std::endl;

      // 3) create graph of amalgamated matrix
      RCP<CrsGraph> crsGraph = CrsGraphFactory::Build(nodeMap, 10, Xpetra::DynamicProfile);

      LO numRows = A->getRowMap()->getNodeNumElements();
      LO numNodes = nodeMap->getNodeNumElements();
      const ArrayRCP<bool> amalgBoundaryNodes(numNodes, false);
      const ArrayRCP<int>  numberDirichletRowsPerNode(numNodes, 0); // helper array counting the number of Dirichlet nodes associated with node
      bool bIsDiagonalEntry = false;       // boolean flag stating that grid==gcid

      // 4) do amalgamation. generate graph of amalgamated matrix
      //    Note, this code is much more inefficient than the leightwight implementation
      //    Most of the work has already been done in the AmalgamationFactory
      for(LO row=0; row<numRows; row++) {
        // get global DOF id
        GO grid = rowMap->getGlobalElement(row);

        // reinitialize boolean helper variable
        bIsDiagonalEntry = false;

        // translate grid to nodeid
        GO nodeId = AmalgamationFactory::DOFGid2NodeId(grid, blockdim, offset, indexBase);

        size_t nnz = A->getNumEntriesInLocalRow(row);
        Teuchos::ArrayView<const LO> indices;
        Teuchos::ArrayView<const SC> vals;
        A->getLocalRowView(row, indices, vals);

        RCP<std::vector<GO> > cnodeIds = Teuchos::rcp(new std::vector<GO>);  // global column block ids
        LO realnnz = 0;
        for(LO col=0; col<Teuchos::as<LO>(nnz); col++) {
          GO gcid = colMap->getGlobalElement(indices[col]); // global column id

          if(vals[col]!=0.0) {
            GO cnodeId = AmalgamationFactory::DOFGid2NodeId(gcid, blockdim, offset, indexBase);
            cnodeIds->push_back(cnodeId);
            realnnz++; // increment number of nnz in matrix row
            if (grid == gcid) bIsDiagonalEntry = true;
          }
        }

        if(realnnz == 1 && bIsDiagonalEntry == true) {
          LO lNodeId = nodeMap->getLocalElement(nodeId);
          numberDirichletRowsPerNode[lNodeId] += 1;      // increment Dirichlet row counter associated with lNodeId
          if (numberDirichletRowsPerNode[lNodeId] == blockdim) // mark full Dirichlet nodes
            amalgBoundaryNodes[lNodeId] = true;
        }

        Teuchos::ArrayRCP<GO> arr_cnodeIds = Teuchos::arcp( cnodeIds );

        if(arr_cnodeIds.size() > 0 )
          crsGraph->insertGlobalIndices(nodeId, arr_cnodeIds());
      }
      // fill matrix graph
      crsGraph->fillComplete(nodeMap,nodeMap);

      // 5) create MueLu Graph object
      RCP<GraphBase> graph = rcp(new Graph(crsGraph, "amalgamated graph of A"));

      // Detect and record rows that correspond to Dirichlet boundary conditions
      graph->SetBoundaryNodeMap(amalgBoundaryNodes);

      if (GetVerbLevel() & Statistics1) {
        GO numLocalBoundaryNodes  = 0;
        GO numGlobalBoundaryNodes = 0;
        for (LO i = 0; i < amalgBoundaryNodes.size(); ++i)
          if (amalgBoundaryNodes[i])
            numLocalBoundaryNodes++;
        RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
        MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
        GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes << " Dirichlet nodes" << std::endl;
      }

      // 6) store results in Level
      //graph->SetBoundaryNodeMap(gBoundaryNodeMap);
      Set(currentLevel, "DofsPerNode", blockdim);
      Set(currentLevel, "Graph", graph);

    } //if (doExperimentalWrap) ... else ...

  } //Build

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MergeRows(const Matrix& A, const LO row, Array<LO>& cols, const Array<LO>& translation) const {
    typedef typename ArrayView<const LO>::size_type size_type;

    // extract striding information
    LO blkSize = A.GetFixedBlockSize();  //< stores the size of the block within the strided map
    if (A.IsView("stridedMaps") == true) {
      Teuchos::RCP<const Map> myMap = A.getRowMap("stridedMaps");
      Teuchos::RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(myMap);
      TEUCHOS_TEST_FOR_EXCEPTION(strMap == null, Exceptions::RuntimeError, "Map is not of type StridedMap");
      if (strMap->getStridedBlockId() > -1)
        blkSize = Teuchos::as<LO>(strMap->getStridingData()[strMap->getStridedBlockId()]);
    }

    // count nonzero entries in all dof rows associated with node row
    size_t nnz = 0, pos = 0;
    for (LO j = 0; j < blkSize; j++)
      nnz += A.getNumEntriesInLocalRow(row*blkSize+j);

    if (nnz == 0) {
      cols.resize(0);
      return;
    }

    cols.resize(nnz);

    // loop over all local dof rows associated with local node "row"
    ArrayView<const LO> inds;
    ArrayView<const SC> vals;
    for (LO j = 0; j < blkSize; j++) {
      A.getLocalRowView(row*blkSize+j, inds, vals);
      size_type numIndices = inds.size();

      if (numIndices == 0) // skip empty dof rows
        continue;

      // cols: stores all local node ids for current local node id "row"
      cols[pos++] = translation[inds[0]];
      for (size_type k = 1; k < numIndices; k++) {
        LO nodeID = translation[inds[k]];
        // Here we try to speed up the process by reducing the size of an array
        // to sort. This works if the column nonzeros belonging to the same
        // node are stored consequently.
        if (nodeID != cols[pos-1])
          cols[pos++] = nodeID;
      }
    }
    cols.resize(pos);
    nnz = pos;

    // Sort and remove duplicates
    std::sort(cols.begin(), cols.end());
    pos = 0;
    for (size_t j = 1; j < nnz; j++)
      if (cols[j] != cols[pos])
        cols[++pos] = cols[j];
    cols.resize(pos+1);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MergeRowsWithDropping(const Matrix& A, const LO row, const ArrayRCP<const SC>& ghostedDiagVals, SC threshold, Array<LO>& cols, const Array<LO>& translation) const {
    typedef typename ArrayView<const LO>::size_type size_type;
    typedef Teuchos::ScalarTraits<SC> STS;

    // extract striding information
    LO blkSize = A.GetFixedBlockSize();  //< stores the size of the block within the strided map
    if (A.IsView("stridedMaps") == true) {
      Teuchos::RCP<const Map> myMap = A.getRowMap("stridedMaps");
      Teuchos::RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(myMap);
      TEUCHOS_TEST_FOR_EXCEPTION(strMap == null, Exceptions::RuntimeError, "Map is not of type StridedMap");
      if (strMap->getStridedBlockId() > -1)
        blkSize = Teuchos::as<LO>(strMap->getStridingData()[strMap->getStridedBlockId()]);
    }

    // count nonzero entries in all dof rows associated with node row
    size_t nnz = 0, pos = 0;
    for (LO j = 0; j < blkSize; j++)
      nnz += A.getNumEntriesInLocalRow(row*blkSize+j);

    if (nnz == 0) {
      cols.resize(0);
      return;
    }

    cols.resize(nnz);

    // loop over all local dof rows associated with local node "row"
    ArrayView<const LO> inds;
    ArrayView<const SC> vals;
    for (LO j = 0; j < blkSize; j++) {
      A.getLocalRowView(row*blkSize+j, inds, vals);
      size_type numIndices = inds.size();

      if (numIndices == 0) // skip empty dof rows
        continue;

      // cols: stores all local node ids for current local node id "row"
      LO prevNodeID = -1;
      for (size_type k = 0; k < numIndices; k++) {
        LO dofID = inds[k];
        LO nodeID = translation[inds[k]];

        // we avoid a square root by using squared values
        typename STS::magnitudeType aiiajj = STS::magnitude(threshold*threshold*ghostedDiagVals[dofID]*ghostedDiagVals[row*blkSize+j]); // eps^2 * |a_ii| * |a_jj|
        typename STS::magnitudeType aij = STS::magnitude(vals[k]*vals[k]);

        // check dropping criterion
        if (aij > aiiajj || (row*blkSize+j == dofID)) {
          // accept entry in graph

          // Here we try to speed up the process by reducing the size of an array
          // to sort. This works if the column nonzeros belonging to the same
          // node are stored consequently.
          if (nodeID != prevNodeID) {
            cols[pos++] = nodeID;
            prevNodeID = nodeID;
          }
        }
      }
    }
    cols.resize(pos);
    nnz = pos;

    // Sort and remove duplicates
    std::sort(cols.begin(), cols.end());
    pos = 0;
    for (size_t j = 1; j < nnz; j++)
      if (cols[j] != cols[pos])
        cols[++pos] = cols[j];
    cols.resize(pos+1);

    return;
  }
} //namespace MueLu

#endif // MUELU_COALESCEDROPFACTORY_DEF_HPP
