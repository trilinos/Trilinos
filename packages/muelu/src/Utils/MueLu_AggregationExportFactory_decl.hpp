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
#include "MueLu_AmalgamationFactory_fwd.hpp"
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
    AggregationExportFactory(const std::string outputFileName = "aggs_level%LEVELID_proc%PROCID.out", const FactoryBase* AggFact = NULL, const FactoryBase* CoalesceDropFact = NULL, const FactoryBase* AmalgFact = NULL);

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

    std::string replaceAll(std::string result, const std::string& replaceWhat, const std::string& replaceWithWhat) const;

  private:
    std::string outputFileName_;            ///< filename template for output
    const FactoryBase* AggFact_;            ///< factory which created aggregates
    const FactoryBase* CoalesceDropFact_;   ///< CoalesceAndDropFactory (needed for DofsPerNode variable)
    const FactoryBase* AmalgFact_;          ///< AmalgamationFactory (needed for UnAmalgamationInfo variable)

  }; // class AggregationExportFactory

} // namespace MueLu

#define MUELU_AGGREGATIONEXPORTFACTORY_SHORT

#endif /* MUELU_AGGREGATIONEXPORTFACTORY_DECL_HPP_ */
