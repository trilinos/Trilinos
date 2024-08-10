// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_NOTAYAGGREGATIONFACTORY_DECL_HPP_
#define MUELU_NOTAYAGGREGATIONFACTORY_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>

#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_LWGraph_fwd.hpp"
#include "MueLu_LWGraph_kokkos_fwd.hpp"
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
