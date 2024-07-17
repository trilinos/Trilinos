// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_AGGREGATEQUALITYESTIMATEFACTORY_DECL_HPP
#define MUELU_AGGREGATEQUALITYESTIMATEFACTORY_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_AggregateQualityEstimateFactory_fwd.hpp"

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>

#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_Level_fwd.hpp"

namespace MueLu {

/*!
  @class AggregateQualityEstimateFactory class.
  @brief An factory which assigns each aggregate a quality
  estimate. Originally developed by Napov and Notay in the
  context of plain aggregation, while this quality estimate
  does not correspond to a robust convergence guarentee (as
  it does for plain aggregation), we find empirically that
  it is a good way of discovering poorly constructed aggregates
  even in the smoothed aggregation context.

  Napov, A., & Notay, Y. (2012). An algebraic multigrid method
  with guaranteed convergence rate. SIAM journal on scientific
  computing, 34(2), A1079-A1109.
*/

template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class AggregateQualityEstimateFactory : public SingleLevelFactoryBase {
#undef MUELU_AGGREGATEQUALITYESTIMATEFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;

  //! Constructor.
  AggregateQualityEstimateFactory();

  //! Destructor.
  virtual ~AggregateQualityEstimateFactory();

  //@}

  RCP<const ParameterList> GetValidParameterList() const;

  //! @name Input
  //@{

  /*! @brief Specifies the data that this class needs, and the factories that generate that data.

    If the Build method of this class requires some data, but the generating factory is not specified in DeclareInput, then this class
    will fall back to the settings in FactoryManager.
  */
  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  //! Build aggregate quality esimates with this factory.
  void Build(Level& currentLevel) const;

  //@}

  //! @name Utility method to convert aggregate data to a convenient format.
  //@{

  //! Build aggregate quality esimates with this factory.
  static void ConvertAggregatesData(RCP<const Aggregates> aggs, ArrayRCP<LO>& aggSortedVertices, ArrayRCP<LO>& aggsToIndices, ArrayRCP<LO>& aggSizes);

  //@}

 private:
  //! @name Internal method for computing aggregate quality.
  //@{

  void ComputeAggregateQualities(RCP<const Matrix> A, RCP<const Aggregates> aggs, RCP<Xpetra::MultiVector<magnitudeType, LO, GO, Node>> agg_qualities) const;

  void ComputeAggregateSizes(RCP<const Matrix> A, RCP<const Aggregates> aggs, RCP<LocalOrdinalVector> agg_sizes) const;

  //@}

  //! @name Internal method for outputting aggregate quality
  //@{

  void OutputAggQualities(const Level& level, RCP<const Xpetra::MultiVector<magnitudeType, LO, GO, Node>> agg_qualities) const;

  void OutputAggSizes(const Level& level, RCP<const LocalOrdinalVector> agg_sizes) const;

  //@}

};  // class AggregateQualityEsimateFactory();

}  // namespace MueLu

#define MUELU_AGGREGATEQUALITYESTIMATEFACTORY_SHORT
#endif  // MUELU_DEMOFACTORY_DECL_HPP
