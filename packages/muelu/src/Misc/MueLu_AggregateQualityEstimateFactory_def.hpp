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
#ifndef MUELU_AGGREGATEQUALITYESTIMATEFACTORY_DEF_HPP
#define MUELU_AGGREGATEQUALITYESTIMATEFACTORY_DEF_HPP

#include "MueLu_AggregateQualityEstimateFactory_decl.hpp"

#include "MueLu_Level.hpp"

#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_LAPACK.hpp>

#include "MueLu_Aggregates_decl.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  AggregateQualityEstimateFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AggregateQualityEstimateFactory()
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  AggregateQualityEstimateFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~AggregateQualityEstimateFactory() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregateQualityEstimateFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {

      Input(currentLevel, "A");
      Input(currentLevel, "Aggregates");
      Input(currentLevel, "CoarseMap");

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregateQualityEstimateFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level & currentLevel) const {

      return BuildP(currentLevel);

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregateQualityEstimateFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level & currentLevel) const {

      RCP<Matrix> A = Get<RCP<Matrix>>(currentLevel, "A");
      RCP<Aggregates> aggregates = Get<RCP<Aggregates>>(currentLevel, "Aggregates");

      RCP<const Map> map = Get< RCP<const Map> >(currentLevel, "CoarseMap");
      RCP<MultiVector> aggregate_qualities = MultiVectorFactory::Build(map, 1);

      assert(!aggregates->AggregatesCrossProcessors());
      ComputeAggregateQualities(A, aggregates, aggregate_qualities);

      Set(currentLevel, "AggregateQualities", aggregate_qualities);

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregateQualityEstimateFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ComputeAggregateQualities(RCP<Matrix> A, RCP<Aggregates> aggs, RCP<MultiVector> agg_qualities) const {

      const double SCALAR_ZERO = Teuchos::ScalarTraits<SC>::zero();
      const double SCALAR_ONE = Teuchos::ScalarTraits<SC>::one();
      const double SCALAR_TWO = SCALAR_ONE + SCALAR_ONE;

      const double LO_ZERO = Teuchos::OrdinalTraits<LO>::zero();
      const double LO_ONE = Teuchos::OrdinalTraits<LO>::one();
      
      // Reorder local aggregate information into a format amenable to computing
      // per-aggregate quantities. Specifically, we compute a format
      // similar to compressed sparse row format for sparse matrices in which
      // we store all the local vertices in a single array in blocks corresponding
      // to aggregates. (This array is aggSortedVertices.) We then store a second
      // array (aggsToIndices) whose k-th element stores the index of the first
      // vertex in aggregate k in the array aggSortedVertices.

      const RCP<LOMultiVector> vertex2AggId = aggs->GetVertex2AggId();
      const ArrayRCP<const LO> vertex2AggIdData = vertex2AggId->getData(0);

      LO numAggs = aggs->GetNumAggregates();
      Teuchos::ArrayRCP<LO> aggSizes = aggs->ComputeAggregateSizes();

      Teuchos::ArrayRCP<LO> aggsToIndices(numAggs+LO_ONE,LO_ZERO);

      for (LO i=0;i<numAggs;++i) {
	  aggsToIndices[i+LO_ONE] = aggsToIndices[i] + aggSizes[i];
      }

      LO numNodes = vertex2AggId->getLocalLength();
      Teuchos::ArrayRCP<LO> aggSortedVertices(numNodes,-LO_ONE);
      Teuchos::ArrayRCP<LO> vertexInsertionIndexByAgg(numNodes,LO_ZERO);
      
      for (LO i=0;i<numNodes;++i) {

	  LO aggId = vertex2AggIdData[i];
	  aggSortedVertices[aggsToIndices[aggId]+vertexInsertionIndexByAgg[aggId]] = i;
	  vertexInsertionIndexByAgg[aggId]++;

      }

      // Compute the per-aggregate quality estimate

      typedef Teuchos::SerialDenseMatrix<LO,SC> DenseMatrix;
      typedef Teuchos::SerialDenseVector<LO,SC> DenseVector;

      ArrayView<const LO> rowIndices;
      ArrayView<const SC> rowValues;
      Teuchos::LAPACK<LO,SC> myLapack;

#include <iostream>

      // Iterate over each aggregate to compute the quality estimate
      for (LO aggId=LO_ZERO; aggId<numAggs; ++aggId) {

          LO aggSize = aggSizes[aggId];
	  DenseMatrix A_aggPart(aggSize, aggSize, true);
	  DenseVector offDiagonalAbsoluteSums(aggSize, true);

	  // Iterate over each node in the aggregate
	  for (LO idx=LO_ZERO; idx<aggSize; ++idx) {
              
              LO nodeId = aggSortedVertices[idx+aggsToIndices[aggId]];
	      A->getLocalRowView(nodeId, rowIndices, rowValues);

	      // Iterate over each element in the row corresponding to the current node
	      for (LO elem=LO_ZERO; elem<rowIndices.size();++elem) {

		  LO nodeId2 = rowIndices[elem];
		  SC val = rowValues[elem];

		  LO idxInAgg = -LO_ONE; // -1 if element is not in aggregate

		  // Check whether the element belongs in the aggregate. If it does
		  // find, its index. Otherwise, add it's value to the off diagonal
		  // sums
		  for (LO idx2=LO_ZERO; idx2<aggSize; ++idx2) {

		      if (aggSortedVertices[idx2+aggsToIndices[aggId]] == nodeId2) {

			  // Element does belong to aggregate
			  idxInAgg = idx2;
			  break;

		      }

		  }

		  if (idxInAgg == -LO_ONE) { // Element does not belong to aggregate

		      offDiagonalAbsoluteSums[idx] += fabs(val);

		  } else { // Element does belong to aggregate

		      A_aggPart(idx,idxInAgg) = val;
		      
		  }

	      }

	  }

	  // Construct a diagonal matrix consisting of the diagonal
	  // of A_aggPart
	  DenseMatrix A_aggPartDiagonal(aggSize, aggSize, true);
	  SC diag_sum = SCALAR_ZERO;
	  for (int i=0;i<aggSize;++i) {
	      A_aggPartDiagonal(i,i) = A_aggPart(i,i);
	      diag_sum += A_aggPart(i,i);
	  }

	  DenseMatrix ones(aggSize, aggSize, false);
	  ones.putScalar(SCALAR_ONE);

	  // Compute matrix on top of generalized Rayleigh quotient
	  // topMatrix = A_aggPartDiagonal - A_aggPartDiagonal*ones*A_aggPartDiagonal/diag_sum;
	  DenseMatrix tmp(aggSize, aggSize, false);
	  DenseMatrix topMatrix(A_aggPartDiagonal);

	  tmp.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, SCALAR_ONE, ones, A_aggPartDiagonal, SCALAR_ZERO);
	  topMatrix.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, -SCALAR_ONE/diag_sum, A_aggPartDiagonal, tmp, SCALAR_ONE);

	  // Compute matrix on bottom of generalized Rayleigh quotient
	  DenseMatrix bottomMatrix(A_aggPart);
	  SC matrixNorm = A_aggPart.normInf();
	  for (int i=0;i<aggSize;++i){
	      // Include a small perturbation to the bottom matrix to make it nonsingular
	      bottomMatrix(i,i) -= offDiagonalAbsoluteSums(i) - 1e-12*matrixNorm;
	  }

	  // Compute generalized eigenvalues

	  // GGES only works on real matrices
	  assert(!Teuchos::ScalarTraits<SC>::isComplex);

	  LO sdim, info;
	  DenseVector alpha_real(aggSize, false);
	  DenseVector alpha_imag(aggSize, false);
	  DenseVector beta(aggSize, false);

	  DenseVector workArray(14*aggSize, false);

	  LO (*ptr2func)(SC*, SC*, SC*);
	  ptr2func = NULL;
	  LO* bwork = NULL;
	  SC* vl = NULL;
	  SC* vr = NULL;
	  
	  const char NO='N';
	  myLapack.GGES(NO,NO,NO,ptr2func,aggSize,topMatrix.values(),aggSize,
			bottomMatrix.values(), aggSize, &sdim, alpha_real.values(),
			alpha_imag.values(), beta.values(), vl, aggSize,
			vr, aggSize, workArray.values(), workArray.length(),
			bwork, &info);

	  assert(info == LO_ZERO);

	  SC maxEigenVal = SCALAR_ZERO;

	  for (int i=LO_ZERO;i<aggSize;++i) {
	      
	      assert(fabs(alpha_imag[i]) <= 1e-8*fabs(alpha_real[i])); // Eigenvalues should be nearly real
	      maxEigenVal = std::max(maxEigenVal, alpha_real[i]/beta[i]);

	  }

	  (agg_qualities->getDataNonConst(0))[aggId] = SCALAR_TWO*maxEigenVal;

      }

			  
  }
  

} // namespace MueLu

#define MUELU_AGGREGATEQUALITYESTIMATE_SHORT
#endif // MUELU_AGGREGATEQUALITYESTIMATE_DEF_HPP
