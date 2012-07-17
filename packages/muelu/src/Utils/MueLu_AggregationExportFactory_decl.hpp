/*
 * MueLu_AggregationExportFactory_decl.hpp
 *
 *  Created on: Feb 10, 2012
 *      Author: wiesner
 */

#ifndef MUELU_AGGREGATIONEXPORTFACTORY_DECL_HPP_
#define MUELU_AGGREGATIONEXPORTFACTORY_DECL_HPP_

#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_CrsOperator_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_AggregationExportFactory_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"

namespace MueLu {

  class Level;

  /*!
    @class ThresholdAFilterFactory class.
    @brief Factory for building a thresholded operator.

  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class AggregationExportFactory : public TwoLevelFactoryBase {
#undef MUELU_AGGREGATIONEXPORTFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    AggregationExportFactory(const std::string outputFileName = "aggs_level%LEVELID_proc%PROCID.out", const FactoryBase* AggFact = NULL, const FactoryBase* CoalesceDropFact = NULL);

    //! Destructor.
    virtual ~AggregationExportFactory();
    //@}

    //! Input
    //@{

    void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    void Build(Level &fineLevel, Level &coarseLevel) const;

    //@}
   

  private:

    /*! @brief ComputeAggregateSizes
     * computes the size of the aggregates (in DOFs). This routine should be the same as in TentativePFactory
     */
    void ComputeAggregateSizes(const Aggregates& aggregates, const AmalgamationInfo& amalgInfo, Teuchos::ArrayRCP<LocalOrdinal> & aggSizes) const;

    /*! @brief ComputeAggregateToRowMap
     * computes the map: aggregate id -> row map. This routine should be the same as in TentativePFactory
     */
    void ComputeAggregateToRowMap(const Aggregates& aggregates, const AmalgamationInfo& amalgInfo, const Teuchos::ArrayRCP<LocalOrdinal> & aggSizes, Teuchos::ArrayRCP<Teuchos::ArrayRCP<GlobalOrdinal> > & aggToRowMap) const;

    
    std::string replaceAll(std::string result, const std::string& replaceWhat, const std::string& replaceWithWhat) const;

  private:
    std::string outputFileName_;            ///< filename template for output
    const FactoryBase* AggFact_;            ///< factory which created aggregates
    const FactoryBase* CoalesceDropFact_;   ///< CoalesceAndDropFactory (needed for DofsPerNode variable)


  }; // class AggregationExportFactory

} // namespace MueLu

#define MUELU_AGGREGATIONEXPORTFACTORY_SHORT

#endif /* MUELU_AGGREGATIONEXPORTFACTORY_DECL_HPP_ */
