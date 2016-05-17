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

  namespace { // anonymous

    template<class LocalOrdinal, class RowType>
    class ScanFunctor {
    public:
      ScanFunctor(RowType rows) : rows_(rows) { }

      KOKKOS_INLINE_FUNCTION
      void operator()(const LocalOrdinal i, LocalOrdinal& upd, const bool& final) const {
        upd += rows_(i);
        if (final)
          rows_(i) = upd;
      }

    private:
      RowType rows_;
    };

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  RCP<const ParameterList> CoalesceDropFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::GetValidParameterList() const {
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

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CoalesceDropFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::DeclareInput(Level &currentLevel) const {
    Input(currentLevel, "A");
    Input(currentLevel, "UnAmalgamationInfo");

    const ParameterList& pL = GetParameterList();
    if (pL.get<bool>("lightweight wrap") == true) {
      if (pL.get<std::string>("aggregation: drop scheme") == "distance laplacian")
        Input(currentLevel, "Coordinates");
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CoalesceDropFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::Build(Level& currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    typedef Teuchos::ScalarTraits<SC> STS;
    SC zero = STS::zero(), one = STS::one();

    auto A         = Get< RCP<Matrix> >(currentLevel, "A");
    LO   blkSize   = A->GetFixedBlockSize();
    GO   indexBase = A->getRowMap()->getIndexBase();

    auto amalInfo = Get< RCP<AmalgamationInfo> >(currentLevel, "UnAmalgamationInfo");

    const ParameterList& pL = GetParameterList();

    // bool doLightWeightWrap = pL.get<bool>("lightweight wrap");
    GetOStream(Warnings0) << "lightweight wrap is deprecated" << std::endl;

    std::string algo = pL.get<std::string>("aggregation: drop scheme");

    double threshold = pL.get<double>("aggregation: drop tol");
    GetOStream(Runtime0) << "algorithm = \"" << algo << "\": threshold = " << threshold << ", blocksize = " << A->GetFixedBlockSize() << std::endl;

    Set(currentLevel, "Filtering", (threshold != STS::zero()));

    const typename STS::magnitudeType dirichletThreshold = STS::magnitude(as<SC>(pL.get<double>("aggregation: Dirichlet threshold")));

    GO numDropped = 0, numTotal = 0;

    RCP<LWGraph_kokkos> graph;
    LO                  dofsPerNode = -1;

    typedef typename LWGraph_kokkos::boundary_nodes_type boundary_nodes_type;
    boundary_nodes_type boundaryNodes;

    if (algo == "classical") {

      if (blkSize == 1 && threshold == STS::zero()) {
        //
        // Case 1:  scalar problem without dropping
        //
        graph = rcp(new LWGraph_kokkos(A->getLocalMatrix().graph, A->getRowMap(), A->getColMap(), "graph of A"));

        // Detect and record rows that correspond to Dirichlet boundary conditions
        boundaryNodes = Utilities_kokkos::DetectDirichletRows(*A, dirichletThreshold);
        graph->SetBoundaryNodeMap(boundaryNodes);

        numTotal = A->getNodeNumEntries();

        dofsPerNode = 1;

      } else if (blkSize == 1 && threshold != STS::zero()) {
        //
        // Case 2:  scalar problem with filtering
        //
        RCP<Vector> ghostedDiagVec = Utilities_kokkos::GetMatrixOverlappedDiagonal(*A);
        auto        ghostedDiag    = ghostedDiagVec->template getLocalView<typename NO::device_type>();

        typedef typename Matrix::local_matrix_type          local_matrix_type;
        typedef typename LWGraph_kokkos::local_graph_type   kokkos_graph_type;
        typedef typename kokkos_graph_type::row_map_type    row_map_type;
        typedef typename kokkos_graph_type::entries_type    entries_type;

        LO   numRows      = A->getNodeNumRows();
        auto kokkosMatrix = A->getLocalMatrix();

        typedef Kokkos::ArithTraits<SC>     ATS;
        typedef typename ATS::magnitudeType magnitudeType;

        // Stage 1: calculate the number of remaining entries per row
        typename row_map_type::non_const_type rows("row_map", numRows+1);       // rows(0) = 0 automatically

        LO realnnz = 0;
        Kokkos::parallel_reduce("MueLu:CoalesceDropF:Build:scalar_filter:stage1_reduce", numRows, KOKKOS_LAMBDA(const LO row, LO& nnz) {
          auto rowView = kokkosMatrix.template row<LO>(row);
          auto length  = rowView.length;

          LO rownnz = 0;
          for (decltype(length) colID = 0; colID < length; colID++) {
            LO col = rowView.colidx(colID);

            // Avoid square root by using squared values
            magnitudeType aiiajj = threshold*threshold * ATS::magnitude(ghostedDiag(row, 0))*ATS::magnitude(ghostedDiag(col, 0));   // eps^2*|a_ii|*|a_jj|
            magnitudeType aij2   = ATS::magnitude(rowView.value(colID)) * ATS::magnitude(rowView.value(colID));                     // |a_ij|^2

            if (aij2 > aiiajj || row == col)
              rownnz++;
          }
          rows(row+1) = rownnz;
          nnz += rownnz;
        }, realnnz);

        // parallel_scan (exclusive)
        // NOTE: parallel_scan with KOKKOS_LAMBDA does not work with CUDA for now
#if 0
        Kokkos::parallel_scan("MueLu:CoalesceDropF:Build:scalar_filter:stage1_scan", numRows+1, KOKKOS_LAMBDA(const LO i, LO& upd, const bool& final) {
          upd += rows(i);
          if (final)
            rows(i) = upd;
        });
#else
        ScanFunctor<LO,decltype(rows)> scanFunctor(rows);
        Kokkos::parallel_scan("MueLu:CoalesceDropF:Build:scalar_filter:stage1_scan", numRows+1, scanFunctor);
#endif


        // Stage 2: fill in the column indices
        typename boundary_nodes_type::non_const_type bndNodes("boundaryNodes", numRows);
        typename entries_type::non_const_type        cols    ("entries",       realnnz);
        Kokkos::parallel_reduce("MueLu:CoalesceDropF:Build:scalar_filter:stage2_reduce", numRows, KOKKOS_LAMBDA(const LO row, GO& dropped) {
          auto rowView = kokkosMatrix.template row<LO>(row);
          auto length = rowView.length;

          LO rownnz = 0;
          for (decltype(length) colID = 0; colID < length; colID++) {
            LO col = rowView.colidx(colID);

            // Avoid square root by using squared values
            magnitudeType aiiajj = threshold*threshold * ATS::magnitude(ghostedDiag(row, 0))*ATS::magnitude(ghostedDiag(col, 0));   // eps^2*|a_ii|*|a_jj|
            magnitudeType aij2   = ATS::magnitude(rowView.value(colID)) * ATS::magnitude(rowView.value(colID));                     // |a_ij|^2

            if (aij2 > aiiajj || row == col) {
              cols(rows(row) + rownnz) = col;
              rownnz++;
            } else {
              dropped++;
            }
            if (rownnz == 1) {
              // If the only element remaining after filtering is diagonal, mark node as boundary
              // FIXME: this should really be replaced by the following
              //    if (indices.size() == 1 && indices[0] == row)
              //        boundaryNodes[row] = true;
              // We do not do it this way now because there is no framework for distinguishing isolated
              // and boundary nodes in the aggregation algorithms
              bndNodes(row) = true;
            }
          }
        }, numDropped);
        boundaryNodes = bndNodes;

        kokkos_graph_type kokkosGraph(cols, rows);

        graph = rcp(new LWGraph_kokkos(kokkosGraph, A->getRowMap(), A->getColMap(), "filtered graph of A"));
        graph->SetBoundaryNodeMap(boundaryNodes);

        numTotal = A->getNodeNumEntries();

        dofsPerNode = 1;

      } else if (blkSize > 1 && threshold == STS::zero()) {
        //
        // Case 3:  block problem without filtering
        //
        throw Exceptions::RuntimeError("Block systems without filtering are not implemented");

        // Detect and record rows that correspond to Dirichlet boundary conditions
        // boundary_nodes_type pointBoundaryNodes = Utilities_kokkos::DetectDirichletRows(*A, dirichletThreshold);

      } else if (blkSize > 1 && threshold != STS::zero()) {
        //
        // Case 4:  block problem with filtering
        //
        throw Exceptions::RuntimeError("Block systems with filtering are not implemented");

        // Detect and record rows that correspond to Dirichlet boundary conditions
        // boundary_nodes_type pointBoundaryNodes = Utilities_kokkos::DetectDirichletRows(*A, dirichletThreshold);
      }

    } else if (algo == "distance laplacian") {
      typedef Xpetra::MultiVector<double,LO,GO,NO> doubleMultiVector;

      auto coords = Get<RCP<doubleMultiVector> >(currentLevel, "Coordinates");

      if (A->GetFixedBlockSize() == 1 && threshold == STS::zero()) {
        //
        // Case 1:  scalar problem without dropping
        //
        graph = rcp(new LWGraph_kokkos(A->getLocalMatrix().graph, A->getRowMap(), A->getColMap(), "graph of A"));

        // Detect and record rows that correspond to Dirichlet boundary conditions
        boundaryNodes = Utilities_kokkos::DetectDirichletRows(*A, dirichletThreshold);
        graph->SetBoundaryNodeMap(boundaryNodes);

        numTotal = A->getNodeNumEntries();

        dofsPerNode = 1;

      } else if (blkSize == 1 && threshold != STS::zero()) {
        //
        // Case 2:  scalar problem with filtering
        //
        throw Exceptions::RuntimeError("not implemented");

      } else if (blkSize > 1 && threshold == STS::zero()) {
        //
        // Case 3:  block problem without filtering
        //
        throw Exceptions::RuntimeError("Block systems without filtering are not implemented");

        // Detect and record rows that correspond to Dirichlet boundary conditions
        boundary_nodes_type pointBoundaryNodes = Utilities_kokkos::DetectDirichletRows(*A, dirichletThreshold);

      } else if (blkSize > 1 && threshold != STS::zero()) {
        //
        // Case 4:  block problem with filtering
        //
        throw Exceptions::RuntimeError("Block systems with filtering are not implemented");

        // Detect and record rows that correspond to Dirichlet boundary conditions
        boundary_nodes_type pointBoundaryNodes = Utilities_kokkos::DetectDirichletRows(*A, dirichletThreshold);
      }

    } else if (algo == "distance laplacian") {
      typedef Xpetra::MultiVector<double,LO,GO,NO> doubleMultiVector;

      auto coords = Get<RCP<doubleMultiVector> >(currentLevel, "Coordinates");

      if (A->GetFixedBlockSize() == 1 && threshold == STS::zero()) {
        //
        // Case 1:  scalar problem without dropping
        //
        graph = rcp(new LWGraph_kokkos(A->getLocalMatrix().graph, A->getRowMap(), A->getColMap(), "graph of A"));

        // Detect and record rows that correspond to Dirichlet boundary conditions
        boundaryNodes = Utilities_kokkos::DetectDirichletRows(*A, dirichletThreshold);
        graph->SetBoundaryNodeMap(boundaryNodes);

        numTotal = A->getNodeNumEntries();

        dofsPerNode = 1;

      } else if (blkSize == 1 && threshold != STS::zero()) {
        //
        // Case 2:  scalar problem with filtering
        //
        throw Exceptions::RuntimeError("not implemented");

      } else if (blkSize > 1 && threshold == STS::zero()) {
        //
        // Case 3:  block problem without filtering
        //
        throw Exceptions::RuntimeError("Block systems without filtering are not implemented");

        // Detect and record rows that correspond to Dirichlet boundary conditions
        // boundary_nodes_type pointBoundaryNodes = Utilities_kokkos::DetectDirichletRows(*A, dirichletThreshold);

      } else if (blkSize > 1 && threshold != STS::zero()) {
        //
        // Case 4:  block problem with filtering
        //
        throw Exceptions::RuntimeError("Block systems with filtering are not implemented");

        // Detect and record rows that correspond to Dirichlet boundary conditions
        // boundary_nodes_type pointBoundaryNodes = Utilities_kokkos::DetectDirichletRows(*A, dirichletThreshold);
      }


    }

    if (GetVerbLevel() & Statistics1) {
      GO numLocalBoundaryNodes  = 0;
      GO numGlobalBoundaryNodes = 0;
      Kokkos::parallel_reduce("MueLu:CoalesceDropF:Build:bnd", boundaryNodes.dimension_0(), KOKKOS_LAMBDA(const LO i, GO& n) {
        if (boundaryNodes(i))
          n++;
      }, numLocalBoundaryNodes);

      RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
      MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
      GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes << " Dirichlet nodes" << std::endl;
    }

    if ((GetVerbLevel() & Statistics1) && threshold != STS::zero()) {
      RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
      GO numGlobalTotal, numGlobalDropped;
      MueLu_sumAll(comm, numTotal,   numGlobalTotal);
      MueLu_sumAll(comm, numDropped, numGlobalDropped);
      if (numGlobalTotal != 0) {
        GetOStream(Statistics1) << "Number of dropped entries: "
            << numGlobalDropped << "/" << numGlobalTotal
            << " (" << 100*Teuchos::as<double>(numGlobalDropped)/Teuchos::as<double>(numGlobalTotal) << "%)" << std::endl;
      }
    }

    Set(currentLevel, "DofsPerNode",  dofsPerNode);
    Set(currentLevel, "Graph",        graph);
  }
}
#endif // HAVE_MUELU_KOKKOS_REFACTOR
#endif // MUELU_COALESCEDROPFACTORY_KOKKOS_DEF_HPP
