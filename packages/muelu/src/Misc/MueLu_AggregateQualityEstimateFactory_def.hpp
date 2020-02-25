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
#include <Xpetra_IO.hpp>
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_Utilities.hpp"

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
  RCP<const ParameterList> AggregateQualityEstimateFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("aggregate qualities: good aggregate threshold");
    SET_VALID_ENTRY("aggregate qualities: file output");
    SET_VALID_ENTRY("aggregate qualities: file base");
    SET_VALID_ENTRY("aggregate qualities: check symmetry");
#undef  SET_VALID_ENTRY

    validParamList->set< RCP<const FactoryBase> >("A",                  Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("Aggregates",         Teuchos::null, "Generating factory of the aggregates");
    validParamList->set< RCP<const FactoryBase> >("CoarseMap",          Teuchos::null, "Generating factory of the coarse map");

    return validParamList;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregateQualityEstimateFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level & currentLevel) const {

    FactoryMonitor m(*this, "Build", currentLevel);

    RCP<Matrix> A = Get<RCP<Matrix>>(currentLevel, "A");
    RCP<Aggregates> aggregates = Get<RCP<Aggregates>>(currentLevel, "Aggregates");

    RCP<const Map> map = Get< RCP<const Map> >(currentLevel, "CoarseMap");
    RCP<Xpetra::MultiVector<magnitudeType,LO,GO,NO>> aggregate_qualities = Xpetra::MultiVectorFactory<magnitudeType,LO,GO,NO>::Build(map, 1);

    assert(!aggregates->AggregatesCrossProcessors());
    ComputeAggregateQualities(A, aggregates, aggregate_qualities);

    Set(currentLevel, "AggregateQualities", aggregate_qualities);

    OutputAggQualities(currentLevel, aggregate_qualities);

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregateQualityEstimateFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ConvertAggregatesData(RCP<const Aggregates> aggs, ArrayRCP<LO>& aggSortedVertices, ArrayRCP<LO>& aggsToIndices, ArrayRCP<LO>& aggSizes) {

    // Reorder local aggregate information into a format amenable to computing
    // per-aggregate quantities. Specifically, we compute a format
    // similar to compressed sparse row format for sparse matrices in which
    // we store all the local vertices in a single array in blocks corresponding
    // to aggregates. (This array is aggSortedVertices.) We then store a second
    // array (aggsToIndices) whose k-th element stores the index of the first
    // vertex in aggregate k in the array aggSortedVertices.

    const LO LO_ZERO = Teuchos::OrdinalTraits<LO>::zero();
    const LO LO_ONE = Teuchos::OrdinalTraits<LO>::one();

    LO numAggs = aggs->GetNumAggregates();
    aggSizes = aggs->ComputeAggregateSizes();

    aggsToIndices = ArrayRCP<LO>(numAggs+LO_ONE,LO_ZERO);

    for (LO i=0;i<numAggs;++i) {
      aggsToIndices[i+LO_ONE] = aggsToIndices[i] + aggSizes[i];
    }

    const RCP<LOMultiVector> vertex2AggId = aggs->GetVertex2AggId();
    const ArrayRCP<const LO> vertex2AggIdData = vertex2AggId->getData(0);

    LO numNodes = vertex2AggId->getLocalLength();
    aggSortedVertices = ArrayRCP<LO>(numNodes,-LO_ONE);
    std::vector<LO> vertexInsertionIndexByAgg(numNodes,LO_ZERO);

    for (LO i=0;i<numNodes;++i) {

      LO aggId = vertex2AggIdData[i];
      if (aggId<0 || aggId>=numAggs) continue;

      aggSortedVertices[aggsToIndices[aggId]+vertexInsertionIndexByAgg[aggId]] = i;
      vertexInsertionIndexByAgg[aggId]++;

    }


  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregateQualityEstimateFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ComputeAggregateQualities(RCP<const Matrix> A, RCP<const Aggregates> aggs, RCP<Xpetra::MultiVector<magnitudeType,LO,GO,Node>> agg_qualities) const {

    const SC SCALAR_ONE = Teuchos::ScalarTraits<SC>::one();
    const SC SCALAR_TWO = SCALAR_ONE + SCALAR_ONE;

    const LO LO_ZERO = Teuchos::OrdinalTraits<LO>::zero();
    const LO LO_ONE = Teuchos::OrdinalTraits<LO>::one();

    using MT = magnitudeType;
    const MT MT_ZERO = Teuchos::ScalarTraits<MT>::zero();
    const MT MT_ONE = Teuchos::ScalarTraits<MT>::one();
    ParameterList pL = GetParameterList();

    RCP<const Matrix> AT = A;

    bool check_symmetry = pL.get<bool>("aggregate qualities: check symmetry");
    if (check_symmetry) {

      RCP<MultiVector> x = MultiVectorFactory::Build(A->getMap(), 1, false);
      x->Xpetra_randomize();

      RCP<MultiVector> tmp = MultiVectorFactory::Build(A->getMap(), 1, false);

      A->apply(*x, *tmp, Teuchos::NO_TRANS); // tmp now stores A*x
      A->apply(*x, *tmp, Teuchos::TRANS, -SCALAR_ONE, SCALAR_ONE); // tmp now stores A*x - A^T*x

      Array<magnitudeType> tmp_norm(1);
      tmp->norm2(tmp_norm());

      Array<magnitudeType> x_norm(1);
      tmp->norm2(x_norm());

      if (tmp_norm[0] > 1e-10*x_norm[0]) {
        std::string transpose_string = "transpose";
        RCP<ParameterList> whatever;
        AT = Utilities::Transpose(*rcp_const_cast<Matrix>(A), true, transpose_string, whatever);

        assert(A->getMap()->isSameAs( *(AT->getMap()) ));
      }

    }

    // Reorder local aggregate information into a format amenable to computing
    // per-aggregate quantities. Specifically, we compute a format
    // similar to compressed sparse row format for sparse matrices in which
    // we store all the local vertices in a single array in blocks corresponding
    // to aggregates. (This array is aggSortedVertices.) We then store a second
    // array (aggsToIndices) whose k-th element stores the index of the first
    // vertex in aggregate k in the array aggSortedVertices.      

    ArrayRCP<LO> aggSortedVertices, aggsToIndices, aggSizes;
    ConvertAggregatesData(aggs, aggSortedVertices, aggsToIndices, aggSizes);

    LO numAggs = aggs->GetNumAggregates();

    // Compute the per-aggregate quality estimate

    typedef Teuchos::SerialDenseMatrix<LO,MT> DenseMatrix;
    typedef Teuchos::SerialDenseVector<LO,MT> DenseVector;

    ArrayView<const LO> rowIndices;
    ArrayView<const SC> rowValues;
    ArrayView<const SC> colValues;
    Teuchos::LAPACK<LO,MT> myLapack;

    // Iterate over each aggregate to compute the quality estimate
    for (LO aggId=LO_ZERO; aggId<numAggs; ++aggId) {

      LO aggSize = aggSizes[aggId];
      DenseMatrix A_aggPart(aggSize, aggSize, true);
      DenseVector offDiagonalAbsoluteSums(aggSize, true);

      // Iterate over each node in the aggregate
      for (LO idx=LO_ZERO; idx<aggSize; ++idx) {

        LO nodeId = aggSortedVertices[idx+aggsToIndices[aggId]];
        A->getLocalRowView(nodeId, rowIndices, rowValues);
        AT->getLocalRowView(nodeId, rowIndices, colValues);

        // Iterate over each element in the row corresponding to the current node
        for (LO elem=LO_ZERO; elem<rowIndices.size();++elem) {

          LO nodeId2 = rowIndices[elem];
          SC val = (rowValues[elem]+colValues[elem])/SCALAR_TWO;

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

            offDiagonalAbsoluteSums[idx] += Teuchos::ScalarTraits<SC>::magnitude(val);

          } else { // Element does belong to aggregate

            A_aggPart(idx,idxInAgg) = Teuchos::ScalarTraits<SC>::real(val);

          }

        }

      }

      // Construct a diagonal matrix consisting of the diagonal
      // of A_aggPart
      DenseMatrix A_aggPartDiagonal(aggSize, aggSize, true);
      MT diag_sum = MT_ZERO;
      for (int i=0;i<aggSize;++i) {
        A_aggPartDiagonal(i,i) = Teuchos::ScalarTraits<SC>::real(A_aggPart(i,i));
        diag_sum += Teuchos::ScalarTraits<SC>::real(A_aggPart(i,i));
      }

      DenseMatrix ones(aggSize, aggSize, false);
      ones.putScalar(MT_ONE);

      // Compute matrix on top of generalized Rayleigh quotient
      // topMatrix = A_aggPartDiagonal - A_aggPartDiagonal*ones*A_aggPartDiagonal/diag_sum;
      DenseMatrix tmp(aggSize, aggSize, false);
      DenseMatrix topMatrix(A_aggPartDiagonal);

      tmp.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, MT_ONE, ones, A_aggPartDiagonal, MT_ZERO);
      topMatrix.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, -MT_ONE/diag_sum, A_aggPartDiagonal, tmp, MT_ONE);

      // Compute matrix on bottom of generalized Rayleigh quotient
      DenseMatrix bottomMatrix(A_aggPart);
      MT matrixNorm = A_aggPart.normInf();
      for (int i=0;i<aggSize;++i){
        // Include a small perturbation to the bottom matrix to make it nonsingular
        bottomMatrix(i,i) -= offDiagonalAbsoluteSums(i) - 1e-12*matrixNorm;
      }

      // Compute generalized eigenvalues
      LO sdim, info;
      DenseVector alpha_real(aggSize, false);
      DenseVector alpha_imag(aggSize, false);
      DenseVector beta(aggSize, false);

      DenseVector workArray(14*(aggSize+1), false);

      LO (*ptr2func)(MT*, MT*, MT*);
      ptr2func = NULL;
      LO* bwork = NULL;
      MT* vl = NULL;
      MT* vr = NULL;

      const char compute_flag ='N';
      myLapack.GGES(compute_flag,compute_flag,compute_flag,ptr2func,aggSize,
                    topMatrix.values(),aggSize,bottomMatrix.values(),aggSize,&sdim,
                    alpha_real.values(),alpha_imag.values(),beta.values(),vl,aggSize,
                    vr,aggSize,workArray.values(),workArray.length(),bwork,
                    &info);

      TEUCHOS_ASSERT(info == LO_ZERO);

      MT maxEigenVal = MT_ZERO;

      for (int i=LO_ZERO;i<aggSize;++i) {
        // NOTE: In theory, the eigenvalues should be nearly real
        //TEUCHOS_ASSERT(fabs(alpha_imag[i]) <= 1e-8*fabs(alpha_real[i])); // Eigenvalues should be nearly real
        maxEigenVal = std::max(maxEigenVal, alpha_real[i]/beta[i]);

      }

      (agg_qualities->getDataNonConst(0))[aggId] = (MT_ONE+MT_ONE)*maxEigenVal;

    }


  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregateQualityEstimateFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::OutputAggQualities(const Level& level, RCP<const Xpetra::MultiVector<magnitudeType,LO,GO,Node>> agg_qualities) const {

    ParameterList pL = GetParameterList();

    magnitudeType good_agg_thresh = Teuchos::as<magnitudeType>(pL.get<double>("aggregate qualities: good aggregate threshold"));
    using MT = magnitudeType;

    ArrayRCP<const MT> data = agg_qualities->getData(0);

    LO num_bad_aggs = 0;
    MT worst_agg = 0.0;

    for (size_t i=0;i<agg_qualities->getLocalLength();++i) {

      if (data[i] > good_agg_thresh) num_bad_aggs++;
      worst_agg = std::max(worst_agg, data[i]);

    }

    if (num_bad_aggs == 0) {
      GetOStream(Statistics1) << "All aggregates passed the quality measure. Worst aggregate had quality " << worst_agg << std::endl;
    } else {
      GetOStream(Statistics1) << num_bad_aggs << " out of " << agg_qualities->getLocalLength() << " did not pass the quality measure. Worst aggregate had quality " << worst_agg << std::endl;
    }

    if (pL.get<bool>("aggregate qualities: file output")) {
      std::string filename = pL.get<std::string>("aggregate qualities: file base")+"."+std::to_string(level.GetLevelID());
      Xpetra::IO<magnitudeType,LO,GO,Node>::Write(filename, *agg_qualities);
    }

  }



} // namespace MueLu

#endif // MUELU_AGGREGATEQUALITYESTIMATE_DEF_HPP
