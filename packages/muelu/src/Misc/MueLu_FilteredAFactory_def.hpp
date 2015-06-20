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
#ifndef MUELU_FILTEREDAFACTORY_DEF_HPP
#define MUELU_FILTEREDAFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixFactory.hpp>

#include "MueLu_FilteredAFactory_decl.hpp"

#include "MueLu_FactoryManager.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("filtered matrix: use lumping");
    SET_VALID_ENTRY("filtered matrix: reuse graph");
    SET_VALID_ENTRY("filtered matrix: reuse eigenvalue");
#undef  SET_VALID_ENTRY

    validParamList->set< RCP<const FactoryBase> >("A",              Teuchos::null, "Generating factory of the matrix A used for filtering");
    validParamList->set< RCP<const FactoryBase> >("Graph",          Teuchos::null, "Generating fatory for coalesced filtered graph");
    validParamList->set< RCP<const FactoryBase> >("Filtering",      Teuchos::null, "Generating factory for filtering boolean");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "A");
    Input(currentLevel, "Filtering");
    Input(currentLevel, "Graph");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
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

    RCP<GraphBase> G = Get< RCP<GraphBase> >(currentLevel, "Graph");

    RCP<ParameterList> fillCompleteParams(new ParameterList);
    fillCompleteParams->set("No Nonlocal Changes", true);

    RCP<Matrix> filteredA;
    if (pL.get<bool>("filtered matrix: reuse graph")) {
      filteredA = MatrixFactory::Build(A->getCrsGraph());
      filteredA->resumeFill();

      BuildReuse(*A, *G, lumping, *filteredA);

      filteredA->fillComplete(fillCompleteParams);

    } else {
      filteredA = MatrixFactory::Build(A->getRowMap(), A->getColMap(), A->getNodeMaxNumRowEntries(), Xpetra::StaticProfile);

      BuildNew(*A, *G, lumping, *filteredA);

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

// Epetra's API allows direct access to row array.
// Tpetra's API does not, providing only ArrayView<const .>
// But in most situations we are currently interested in, it is safe to assume
// that the view is to the actual data. So this macro directs the code to do
// const_cast, and modify the entries directly. This allows us to avoid
// replaceLocalValues() call which is quite expensive due to all the searches.
#define ASSUME_DIRECT_ACCESS_TO_ROW

  // Both Epetra and Tpetra matrix-matrix multiply use the following trick:
  // if an entry of the left matrix is zero, it does not compute or store the
  // zero value.
  //
  // This trick allows us to bypass constructing a new matrix. Instead, we
  // make a deep copy of the original one, and fill it in with zeros, which
  // are ignored during the prolongator smoothing.
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  BuildReuse(const Matrix& A, const GraphBase& G, const bool lumping, Matrix& filteredA) const {
    SC zero = Teuchos::ScalarTraits<SC>::zero();

    size_t blkSize = A.GetFixedBlockSize();

    ArrayView<const LO> inds;
    ArrayView<const SC> valsA;
#ifdef ASSUME_DIRECT_ACCESS_TO_ROW
    ArrayView<SC>       vals;
#else
    Array<SC>           vals;
#endif

    Array<char> filter(blkSize * G.GetImportMap()->getNodeNumElements(), 0);

    size_t numGRows = G.GetNodeNumVertices();
    for (size_t i = 0; i < numGRows; i++) {
      // Set up filtering array
      ArrayView<const LO> indsG = G.getNeighborVertices(i);
      for (size_t j = 0; j < as<size_t>(indsG.size()); j++)
        for (size_t k = 0; k < blkSize; k++)
          filter[indsG[j]*blkSize+k] = 1;

      for (size_t k = 0; k < blkSize; k++) {
        LO row = i*blkSize + k;

        A.getLocalRowView(row, inds, valsA);

        size_t nnz = inds.size();
        if (nnz == 0)
          continue;

#ifdef ASSUME_DIRECT_ACCESS_TO_ROW
        // Transform ArrayView<const SC> into ArrayView<SC>
        ArrayView<const SC> vals1;
        filteredA.getLocalRowView(row, inds, vals1);
        vals = ArrayView<SC>(const_cast<SC*>(vals1.getRawPtr()), nnz);

        memcpy(vals.getRawPtr(), valsA.getRawPtr(), nnz*sizeof(SC));
#else
        vals = Array<SC>(valsA);
#endif

        if (lumping == false) {
          for (size_t j = 0; j < nnz; j++)
            if (!filter[inds[j]])
              vals[j] = zero;

        } else {
          LO diagIndex = -1;
          SC diagExtra = zero;

          for (size_t j = 0; j < nnz; j++) {
            if (filter[inds[j]]) {
              if (inds[j] == row) {
                // Remember diagonal position
                diagIndex = j;
              }
              continue;
            }

            diagExtra += vals[j];

            vals[j] = zero;
          }

          // Lump dropped entries
          // NOTE
          //  * Does it make sense to lump for elasticity?
          //  * Is it different for diffusion and elasticity?
          if (diagIndex != -1)
            vals[diagIndex] += diagExtra;
        }

#ifndef ASSUME_DIRECT_ACCESS_TO_ROW
        // Because we used a column map in the construction of the matrix
        // we can just use insertLocalValues here instead of insertGlobalValues
        filteredA.replaceLocalValues(row, inds, vals);
#endif
      }

      // Reset filtering array
      for (size_t j = 0; j < as<size_t> (indsG.size()); j++)
        for (size_t k = 0; k < blkSize; k++)
          filter[indsG[j]*blkSize+k] = 0;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  BuildNew(const Matrix& A, const GraphBase& G, const bool lumping, Matrix& filteredA) const {
    SC zero = Teuchos::ScalarTraits<SC>::zero();

    size_t blkSize = A.GetFixedBlockSize();

    ArrayView<const LO> indsA;
    ArrayView<const SC> valsA;
    Array<LO>           inds;
    Array<SC>           vals;

    Array<char> filter(blkSize * G.GetImportMap()->getNodeNumElements(), 0);

    size_t numGRows = G.GetNodeNumVertices();
    for (size_t i = 0; i < numGRows; i++) {
      // Set up filtering array
      ArrayView<const LO> indsG = G.getNeighborVertices(i);
      for (size_t j = 0; j < as<size_t>(indsG.size()); j++)
        for (size_t k = 0; k < blkSize; k++)
          filter[indsG[j]*blkSize+k] = 1;

      for (size_t k = 0; k < blkSize; k++) {
        LO row = i*blkSize + k;

        A.getLocalRowView(row, indsA, valsA);

        size_t nnz = indsA.size();
        if (nnz == 0)
          continue;

        inds.resize(indsA.size());
        vals.resize(valsA.size());

        size_t numInds = 0;
        if (lumping == false) {
          for (size_t j = 0; j < nnz; j++)
            if (filter[indsA[j]]) {
              inds[numInds] = indsA[j];
              vals[numInds] = valsA[j];
              numInds++;
            }

        } else {
          LO diagIndex = -1;
          SC diagExtra = zero;

          for (size_t j = 0; j < nnz; j++) {
            if (filter[indsA[j]]) {
              inds[numInds] = indsA[j];
              vals[numInds] = valsA[j];

              // Remember diagonal position
              if (inds[numInds] == row)
                diagIndex = numInds;

              numInds++;

            } else {
              diagExtra += valsA[j];
            }
          }

          // Lump dropped entries
          // NOTE
          //  * Does it make sense to lump for elasticity?
          //  * Is it different for diffusion and elasticity?
          if (diagIndex != -1)
            vals[diagIndex] += diagExtra;
        }
        inds.resize(numInds);
        vals.resize(numInds);

        // Because we used a column map in the construction of the matrix
        // we can just use insertLocalValues here instead of insertGlobalValues
        filteredA.insertLocalValues(row, inds, vals);
      }

      // Reset filtering array
      for (size_t j = 0; j < as<size_t> (indsG.size()); j++)
        for (size_t k = 0; k < blkSize; k++)
          filter[indsG[j]*blkSize+k] = 0;
    }
  }

} //namespace MueLu

#endif // MUELU_FILTEREDAFACTORY_DEF_HPP
