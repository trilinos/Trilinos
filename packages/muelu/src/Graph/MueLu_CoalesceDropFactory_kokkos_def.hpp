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
#ifndef MUELU_COALESCEDROPFACTORY_KOKKOS_DEF_HPP
#define MUELU_COALESCEDROPFACTORY_KOKKOS_DEF_HPP

#ifdef HAVE_MUELU_KOKKOS_REFACTOR
#include <Kokkos_CrsMatrix.hpp>

#include "MueLu_CoalesceDropFactory_kokkos_decl.hpp"

#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities_kokkos.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> CoalesceDropFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
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
    validParamList->set< bool >                  ("lightweight wrap",            true, "Experimental option for lightweight graph access");

    validParamList->set< RCP<const FactoryBase> >("A",                  Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Generating factory for UnAmalgamationInfo");
    validParamList->set< RCP<const FactoryBase> >("Coordinates",        Teuchos::null, "Generating factory for Coordinates");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoalesceDropFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
    Input(currentLevel, "A");
    Input(currentLevel, "UnAmalgamationInfo");

    const ParameterList& pL = GetParameterList();
    if (pL.get<bool>("lightweight wrap") == true) {
      if (pL.get<std::string>("aggregation: drop scheme") == "distance laplacian")
        Input(currentLevel, "Coordinates");
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoalesceDropFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    typedef Teuchos::ScalarTraits<SC> STS;
    SC zero = STS::zero(), one = STS::one();

    RCP<Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");
    RCP<AmalgamationInfo> amalInfo = Get< RCP<AmalgamationInfo> >(currentLevel, "UnAmalgamationInfo");

    const ParameterList& pL = GetParameterList();

    // bool doLightWeightWrap = pL.get<bool>("lightweight wrap");
    GetOStream(Warnings0) << "lightweight wrap is deprecated" << std::endl;

    std::string algo = pL.get<std::string>("aggregation: drop scheme");

    SC threshold = as<SC>(pL.get<double>("aggregation: drop tol"));
    GetOStream(Runtime0) << "algorithm = \"" << algo << "\": threshold = " << threshold << ", blocksize = " << A->GetFixedBlockSize() << std::endl;

    Set(currentLevel, "Filtering", (threshold != STS::zero()));

    const typename STS::magnitudeType dirichletThreshold = STS::magnitude(as<SC>(pL.get<double>("aggregation: Dirichlet threshold")));

    GO numDropped = 0, numTotal = 0;
    std::string graphType = "unamalgamated"; //for description purposes only

    RCP<LWGraph_kokkos> graph;
    if (algo == "classical") {
      if (A->GetFixedBlockSize() == 1 && threshold == STS::zero()) {
        //
        // Case 1:  scalar problem without dropping
        //
        graph = rcp(new LWGraph_kokkos(A->getLocalMatrix().graph, A->getDomainMap(), A->getRangeMap(), "graph of A"));

        // Detect and record rows that correspond to Dirichlet boundary conditions
        auto boundaryNodes = MueLu::Utilities_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Node>::DetectDirichletRows(*A, dirichletThreshold);
        graph->SetBoundaryNodeMap(boundaryNodes);

        numTotal = A->getNodeNumEntries();

        if (GetVerbLevel() & Statistics0) {
          GO numLocalBoundaryNodes  = 0;
          GO numGlobalBoundaryNodes = 0;
          Kokkos::parallel_reduce("CoalesceDropF:Build:case1_bnd", boundaryNodes.dimension_0(), KOKKOS_LAMBDA(const LO i, GO& n) {
            if (boundaryNodes[i])
              n++;
          }, numLocalBoundaryNodes);

          RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
          MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
          GetOStream(Statistics0) << "Detected " << numGlobalBoundaryNodes << " Dirichlet nodes" << std::endl;
        }

        Set(currentLevel, "DofsPerNode",  1);
        Set(currentLevel, "Graph",        graph);

      } else if (A->GetFixedBlockSize() == 1 && threshold != STS::zero()) {

        Set(currentLevel, "DofsPerNode",  1);
        Set(currentLevel, "Graph",        graph);

        //
        // Case 2:  scalar problem with dropping
        //
#if 0
        Kokkos::View rows;
        Kokkos::View cols;

        RCP<Vector>  ghostedDiagVec = Utils_kokkos:GetMatrixOverlappedDiagonal(*A);
        Kokkos::View ghostedDiag    = ghostedDiag->getData(0);

          // allocate space for the local graph
          ArrayRCP<LO> rows   (A->getNodeNumRows()+1);
          ArrayRCP<LO> columns(A->getNodeNumEntries());

          RCP<Vector> ghostedDiag = MueLu::Utils<SC,LO,GO,NO>::GetMatrixOverlappedDiagonal(*A);
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
          if (GetVerbLevel() & Statistics0) {
            GO numLocalBoundaryNodes  = 0;
            GO numGlobalBoundaryNodes = 0;
            for (LO i = 0; i < boundaryNodes.size(); ++i)
              if (boundaryNodes[i])
                numLocalBoundaryNodes++;
            RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
            MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
            GetOStream(Statistics0) << "Detected " << numGlobalBoundaryNodes << " Dirichlet nodes" << std::endl;
          }
          Set(currentLevel, "Graph",       graph);
          Set(currentLevel, "DofsPerNode", 1);

#endif
      }
    }

  }
}
#endif // HAVE_MUELU_KOKKOS_REFACTOR
#endif // MUELU_COALESCEDROPFACTORY_KOKKOS_DEF_HPP
