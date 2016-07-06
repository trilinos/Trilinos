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
#ifndef MUELU_FILTEREDAFACTORY_KOKKOS_DEF_HPP
#define MUELU_FILTEREDAFACTORY_KOKKOS_DEF_HPP

#ifdef HAVE_MUELU_KOKKOS_REFACTOR

#include "MueLu_FilteredAFactory_kokkos_decl.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixFactory.hpp>

#include "MueLu_FactoryManager.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  RCP<const ParameterList> FilteredAFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("filtered matrix: use lumping");
    SET_VALID_ENTRY("filtered matrix: reuse graph");
    SET_VALID_ENTRY("filtered matrix: reuse eigenvalue");
#undef  SET_VALID_ENTRY

    validParamList->set< RCP<const FactoryBase> >("A",              Teuchos::null, "Generating factory of the matrix A used for filtering");
    validParamList->set< RCP<const FactoryBase> >("Graph",          Teuchos::null, "Generating factory for coalesced filtered graph");
    validParamList->set< RCP<const FactoryBase> >("Filtering",      Teuchos::null, "Generating factory for filtering boolean");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void FilteredAFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "A");
    Input(currentLevel, "Filtering");
    Input(currentLevel, "Graph");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void FilteredAFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  Build(Level& currentLevel) const {
    FactoryMonitor m(*this, "Matrix filtering", currentLevel);

    RCP<Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");
    if (Get<bool>(currentLevel, "Filtering") == false) {
      GetOStream(Runtime0) << "Filtered matrix is not being constructed as no filtering is being done" << std::endl;
      Set(currentLevel, "A", A);
      return;
    }

    const ParameterList& pL = GetParameterList();
    bool lumping = pL.get<bool>("filtered matrix: use lumping");
    if (lumping)
      GetOStream(Runtime0) << "Lumping dropped entries" << std::endl;

    RCP<LWGraph_kokkos> graph = Get< RCP<LWGraph_kokkos> >(currentLevel, "Graph");

    RCP<ParameterList> fillCompleteParams(new ParameterList);
    fillCompleteParams->set("No Nonlocal Changes", true);

    RCP<Matrix> filteredA;
    if (pL.get<bool>("filtered matrix: reuse graph")) {
      filteredA = MatrixFactory::Build(A->getCrsGraph());
      filteredA->fillComplete(fillCompleteParams);

      BuildReuse(*A, *graph, lumping, *filteredA);

    } else {
      filteredA = MatrixFactory::Build(A->getRowMap(), A->getColMap(), A->getNodeMaxNumRowEntries(), Xpetra::StaticProfile);

      BuildNew(*A, *graph, lumping, *filteredA);

      filteredA->fillComplete(A->getDomainMap(), A->getRangeMap(), fillCompleteParams);
    }

    filteredA->SetFixedBlockSize(A->GetFixedBlockSize());

    if (pL.get<bool>("filtered matrix: reuse eigenvalue")) {
      // Reuse max eigenvalue from A
      // It is unclear what eigenvalue is the best for the smoothing, but we already may have
      // the D^{-1}A estimate in A, may as well use it.
      // NOTE: ML does that too
      filteredA->SetMaxEigenvalueEstimate(A->GetMaxEigenvalueEstimate());
    }

    Set(currentLevel, "A", filteredA);
  }

  // Both Epetra and Tpetra matrix-matrix multiply use the following trick:
  // if an entry of the left matrix is zero, it does not compute or store the
  // zero value.
  //
  // This trick allows us to bypass constructing a new matrix. Instead, we
  // make a deep copy of the original one, and fill it in with zeros, which
  // are ignored during the prolongator smoothing.
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void FilteredAFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  BuildReuse(const Matrix& A, const LWGraph_kokkos& graph, const bool lumping, Matrix& filteredA) const {
    SC zero = Teuchos::ScalarTraits<SC>::zero();

    size_t blkSize = A.GetFixedBlockSize();

    auto localA  = A        .getLocalMatrix();
    auto localFA = filteredA.getLocalMatrix();

    Kokkos::View<char*> filter("filter", blkSize * graph.GetImportMap()->getNodeNumElements(), 0);

    size_t numGRows = graph.GetNodeNumVertices();
    Kokkos::parallel_for("MueLu:FilteredAF:BuildReuse:for", numGRows, KOKKOS_LAMBDA(const size_t i) {
      // Set up filtering array
      typename LWGraph_kokkos::row_type indsG = graph.getNeighborVertices(i);
      for (size_t j = 0; j < indsG.size(); j++)
        for (size_t k = 0; k < blkSize; k++)
          filter(indsG(j)*blkSize + k) = 1;

      for (size_t k = 0; k < blkSize; k++) {
        LO row = i*blkSize + k;

        auto rowA = localA.row (row);
        auto nnz = rowA.length;

        if (nnz == 0)
          continue;

        auto rowFA = localFA.row (row);

        for (decltype(nnz) j = 0; j < nnz; j++)
          rowFA.value(j) = rowA.value(j);

        if (lumping == false) {
          for (decltype(nnz) j = 0; j < nnz; j++)
            if (!filter(rowA.colidx(j)))
              rowFA.value(j) = zero;

        } else {
          LO diagIndex = -1;
          SC diagExtra = zero;

          for (decltype(nnz) j = 0; j < nnz; j++) {
            if (filter(rowA.colidx(j))) {
              if (rowA.colidx(j) == row) {
                // Remember diagonal position
                diagIndex = j;
              }
              continue;
            }

            // TAW: 3/14/2016: not sure whether this works with CUDA
            diagExtra += Teuchos::as<SC>(rowFA.value(j));

            rowFA.value(j) = zero;
          }

          // Lump dropped entries
          // NOTE
          //  * Does it make sense to lump for elasticity?
          //  * Is it different for diffusion and elasticity?
          if (diagIndex != -1)
            rowFA.value(diagIndex) += diagExtra;
        }
      }

      // Reset filtering array
      for (size_t j = 0; j < indsG.size(); j++)
        for (size_t k = 0; k < blkSize; k++)
          filter(indsG(j)*blkSize + k) = 0;
    });
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void FilteredAFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::
  BuildNew(const Matrix& A, const LWGraph_kokkos& graph, const bool lumping, Matrix& filteredA) const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Not implemented");
  }

} //namespace MueLu

#endif // HAVE_MUELU_KOKKOS_REFACTOR
#endif // MUELU_FILTEREDAFACTORY_KOKKOS_DEF_HPP
