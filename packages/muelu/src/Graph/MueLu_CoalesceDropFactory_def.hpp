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
    if (pL.get<std::string>("algorithm") == "laplacian")
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
      TEUCHOS_TEST_FOR_EXCEPTION(algo != "original" && algo != "laplacian", Exceptions::RuntimeError, "\"algorithm\" must be one of (threshold|laplacian)");

      Scalar threshold = Teuchos::as<Scalar>(pL.get<double>("aggregation threshold"));
      GetOStream(Runtime0, 0) << "algorithm = \"" << algo << "\": threshold = " << threshold << std::endl;

      LocalOrdinal numDropped = 0;
      if (algo == "original") {
        if (predrop_ == null) {
          // ap: this is a hack: had to declare predrop_ as mutable
          predrop_ = rcp(new PreDropFunctionConstVal(threshold));
          GetOStream(Runtime0, 0) << "predrop_ constructed" << std::endl;
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

        const typename STS::magnitudeType dirichletThreshold = STS::magnitude(pL.get<Scalar>("Dirichlet detection threshold"));

        if ( (A->GetFixedBlockSize() == 1) && (threshold == STS::zero()) ) {
          // Case 1:  scalar problem, no dropping => just use matrix graph
          RCP<GraphBase> graph = rcp(new Graph(A->getCrsGraph(), "graph of A"));
          //Detect and record rows that correspond to Dirichlet boundary conditions
          const ArrayRCP<const bool > boundaryNodes = MueLu::Utils<SC,LO,GO,NO>::DetectDirichletRows(*A,dirichletThreshold);
          graph->SetBoundaryNodeMap(boundaryNodes);
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
          const ArrayRCP<bool> boundaryNodes(A->getNodeNumRows(),false);

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

            for (LocalOrdinal colID = 0; colID < Teuchos::as<LocalOrdinal>(nnz); colID++) {
              LocalOrdinal col = indices[colID];

              // we avoid a square root by using squared values
              typename STS::magnitudeType aiiajj = STS::magnitude(threshold*threshold * ghostedDiagVals[col]*ghostedDiagVals[row]);  // eps^2*|a_ii|*|a_jj|
              typename STS::magnitudeType aij    = STS::magnitude(vals[colID]*vals[colID]);                                          // |a_ij|^2

              if (aij > aiiajj || row == col)
                columns[realnnz++] = col;
              else
                numDropped++;
            }
            if (realnnz == 1) boundaryNodes[row] = true;
            rows[row+1] = realnnz;
          }

          RCP<GraphBase> graph = rcp(new LWGraph(rows, columns, A->getCrsGraph(), "amalgamated graph of A"));
          graph->SetBoundaryNodeMap(boundaryNodes);
          Set(currentLevel, "Graph", graph);
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
        // Trivial case: scalar problem, no dropping. Can return original graph
        if ( (A->GetFixedBlockSize() == 1) && (threshold == STS::zero()) ) {
          RCP<GraphBase> graph = rcp(new Graph(A->getCrsGraph(), "graph of A"));
          Set(currentLevel, "DofsPerNode", 1);
          Set(currentLevel, "Graph", graph);
        }

        // TODO: need to be very careful with maps for this case
        if (A->GetFixedBlockSize() > 1)
          throw Exceptions::NotImplemented("Fast CoalesceDrop with multiple DOFS is not yet implemented.");

        // allocate space for the local graph
        ArrayRCP<LocalOrdinal> rows    = ArrayRCP<LO>(A->getNodeNumRows()+1);
        ArrayRCP<LocalOrdinal> columns = ArrayRCP<LO>(A->getNodeNumEntries());

        const RCP<const Map> uniqueMap    = A->getDomainMap();
        const RCP<const Map> nonUniqueMap = A->getColMap();
        RCP<const Import>        importer = ImportFactory::Build(uniqueMap, nonUniqueMap);

        // Get ghost coordinates
        RCP<MultiVector>        Coords = Get< RCP<MultiVector> >(currentLevel, "Coordinates");
        RCP<MultiVector> ghostedCoords = MultiVectorFactory::Build(nonUniqueMap, Coords->getNumVectors());
        ghostedCoords->doImport(*Coords, *importer, Xpetra::INSERT);

        LocalOrdinal numRows = Teuchos::as<LocalOrdinal>(A->getRowMap()->getNodeNumElements());

        // Construct Distance Laplacian diagonal
        RCP<Vector> localLaplDiag = VectorFactory::Build(uniqueMap);
        Teuchos::ArrayRCP<Scalar> localLaplDiagData = localLaplDiag->getDataNonConst(0);
        for (LocalOrdinal row = 0; row < numRows; row++) {
          LocalOrdinal nnz = Teuchos::as<LocalOrdinal>(A->getNumEntriesInLocalRow(row));
          ArrayView<const LocalOrdinal> indices;
          ArrayView<const Scalar>       vals;
          A->getLocalRowView(row, indices, vals);

          for (LocalOrdinal colID = 0; colID < nnz; colID++) {
            LocalOrdinal col = indices[colID];

            if (row != col)
              localLaplDiagData[row] += STS::one()/MueLu::Utils<SC,LO,GO,NO>::Distance2(*ghostedCoords, row, col);
          }
        }
        RCP<Vector> ghostedLaplDiag = VectorFactory::Build(nonUniqueMap);
        ghostedLaplDiag->doImport(*localLaplDiag, *importer, Xpetra::INSERT);
        Teuchos::ArrayRCP<Scalar> ghostedLaplDiagData = ghostedLaplDiag->getDataNonConst(0);

        LocalOrdinal realnnz = 0;

        rows[0] = 0;
        for (LocalOrdinal row = 0; row < numRows; row++) {
          LocalOrdinal nnz = Teuchos::as<LocalOrdinal>(A->getNumEntriesInLocalRow(row));
          ArrayView<const LocalOrdinal> indices;
          ArrayView<const Scalar>       vals;
          A->getLocalRowView(row, indices, vals);

          for (LocalOrdinal colID = 0; colID < nnz; colID++) {
            LocalOrdinal col = indices[colID];

            if (row == col) {
              columns[realnnz++] = col;
              continue;
            }

            Scalar laplVal = STS::one() / MueLu::Utils<SC,LO,GO,NO>::Distance2(*ghostedCoords, row, col);
            typename STS::magnitudeType aiiajj = STS::magnitude(threshold*threshold * ghostedLaplDiagData[row]*ghostedLaplDiagData[col]);
            typename STS::magnitudeType aij    = STS::magnitude(laplVal*laplVal);

            if (aij > aiiajj)
              columns[realnnz++] = col;
            else
              numDropped++;
    
          }
          rows[row+1] = realnnz;
        }

        RCP<GraphBase> graph = rcp(new LWGraph(rows, columns, A->getCrsGraph(), "amalgamated graph of A"));
        Set(currentLevel, "Graph", graph);
        Set(currentLevel, "DofsPerNode", 1);
      }

      GetOStream(Statistics0, 0) << "number of dropped " << numDropped << " (" << 100*Teuchos::as<double>(numDropped)/A->getNodeNumEntries() << "%)" << std::endl;

    } else {

      //what Tobias has implemented

      LocalOrdinal blockdim = 1;         // block dim for fixed size blocks
      GlobalOrdinal offset = 0;          // global offset of dof gids

      // 1) check for blocking/striding information
      if(A->IsView("stridedMaps") &&
         Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap("stridedMaps")) != Teuchos::null) {
        Xpetra::viewLabel_t oldView = A->SwitchToView("stridedMaps"); // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
        RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap());
        TEUCHOS_TEST_FOR_EXCEPTION(strMap == Teuchos::null,Exceptions::BadCast,"MueLu::CoalesceFactory::Build: cast to strided row map failed.");
        blockdim = strMap->getFixedBlockSize(); // TODO shorten code
        offset   = strMap->getOffset();
        oldView = A->SwitchToView(oldView);
        GetOStream(Debug, 0) << "CoalesceDropFactory::Build():" << " found blockdim=" << blockdim << " from strided maps. offset=" << offset << std::endl;
      } else GetOStream(Debug, 0) << "CoalesceDropFactory::Build(): no striding information available. Use blockdim=1 with offset=0" << std::endl;

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
      GetOStream(Debug, 0) << "CoalesceDropFactory::SetupAmalgamationData()" << " # of amalgamated blocks=" << num_blockids << std::endl;

      // 3) generate row map for amalgamated matrix (graph of A)
      //    with same distribution over all procs as row map of A
      Teuchos::ArrayRCP<GlobalOrdinal> arr_gNodeIds = Teuchos::arcp( gNodeIds );
      Teuchos::RCP<Map> nodeMap = MapFactory::Build(A->getRowMap()->lib(), num_blockids, arr_gNodeIds(), A->getRowMap()->getIndexBase(), A->getRowMap()->getComm());
      GetOStream(Debug, 0) << "CoalesceDropFactory: nodeMap " << nodeMap->getNodeNumElements() << "/" << nodeMap->getGlobalNumElements() << " elements" << std::endl;

      /////////////////////// experimental
      // vector of boundary node GIDs on current proc
      //RCP<std::map<GlobalOrdinal,bool> > gBoundaryNodes = Teuchos::rcp(new std::map<GlobalOrdinal,bool>);
      ////////////////////////////////////

      // 4) create graph of amalgamated matrix
      RCP<CrsGraph> crsGraph = CrsGraphFactory::Build(nodeMap, 10, Xpetra::DynamicProfile);

      // 5) do amalgamation. generate graph of amalgamated matrix
      for(LocalOrdinal row=0; row<Teuchos::as<LocalOrdinal>(A->getRowMap()->getNodeNumElements()); row++) {
        // get global DOF id
        GlobalOrdinal grid = A->getRowMap()->getGlobalElement(row);

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
          GlobalOrdinal gcid = A->getColMap()->getGlobalElement(indices[col]); // global column id

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
