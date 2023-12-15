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
#ifndef MUELU_NOTAYAGGREGATIONFACTORY_DECL_HPP_
#define MUELU_NOTAYAGGREGATIONFACTORY_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>

#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_GraphBase_fwd.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"

#include "MueLu_NotayAggregationFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class NotayAggregationFactory : public SingleLevelFactoryBase {
#undef MUELU_NOTAYAGGREGATIONFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name typedefs
  //@{
  using local_matrix_type = typename Matrix::local_matrix_type;
  using device_type       = typename local_matrix_type::device_type;
  using execution_space   = typename device_type::execution_space;
  using magnitude_type    = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using impl_scalar_type  = typename Kokkos::ArithTraits<Scalar>::val_type;
  using row_sum_type      = typename Kokkos::View<impl_scalar_type*, Kokkos::LayoutLeft, device_type>;
  //@}

  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  NotayAggregationFactory(){};

  //! Destructor.
  virtual ~NotayAggregationFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! @name Set/get methods.
  //@{

  // Options shared by all aggregation algorithms

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  /*! @brief Build aggregates. */
  void Build(Level& currentLevel) const;

  /*! @brief Initial aggregation phase. */
  void BuildInitialAggregates(const Teuchos::ParameterList& params,
                              const RCP<const Matrix>& A,
                              const ArrayView<const LO>& orderingVector,
                              const magnitude_type kappa,
                              Aggregates& aggregates,
                              std::vector<unsigned>& aggStat,
                              LO& numNonAggregatedNodes,
                              LO& numDirichletNodes) const;

  /*! @brief Further aggregation phase increases coarsening rate by a factor of ~2 per iteration. */
  void BuildFurtherAggregates(const Teuchos::ParameterList& params,
                              const RCP<const Matrix>& A,
                              const Teuchos::ArrayView<const LO>& orderingVector,
                              const local_matrix_type& coarseA,
                              const magnitude_type kappa,
                              const row_sum_type& rowSum,
                              std::vector<LO>& localAggStat,
                              Array<LO>& localVertex2AggID,
                              LO& numLocalAggregates,
                              LO& numNonAggregatedNodes) const;

  void BuildOnRankLocalMatrix(const local_matrix_type& localA,
                              local_matrix_type& onRankA) const;

  /*! @brief Construction of a local prolongator with values equal to 1.0.  */
  void BuildIntermediateProlongator(const LO numRows,
                                    const LO numDirichletNodes,
                                    const LO numLocalAggregates,
                                    const ArrayView<const LO>& localVertex2AggID,
                                    local_matrix_type& intermediateP) const;

  /*! @brief Implementation of a local Galerkin projection called inside BuildFurtherAggregates. */
  void BuildCoarseLocalMatrix(const local_matrix_type& intermediateP,
                              local_matrix_type& coarseA) const;

  /*! @brief Wrapper for kokkos-kernels' spgemm that takes in CrsMatrix. */
  void localSpGEMM(const local_matrix_type& A,
                   const local_matrix_type& B,
                   const std::string matrixLabel,
                   local_matrix_type& C) const;

  //@}

 private:
};  // class NotayAggregationFactory

}  // namespace MueLu

#define MUELU_NOTAYAGGREGATIONFACTORY_SHORT
#endif /* MUELU_NOTAYAGGREGATIONFACTORY_DECL_HPP_ */
