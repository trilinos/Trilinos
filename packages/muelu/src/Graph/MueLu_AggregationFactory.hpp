#ifndef MUELU_AGGREGATIONFACTORY_HPP
#define MUELU_AGGREGATIONFACTORY_HPP

#include "MueLu_Aggregates.hpp"
#include "MueLu_AggregationOptions.hpp"
#include "MueLu_AggAlgorithm.hpp"
#include "MueLu_AggAlgorithm2.hpp"

#include <iostream>

namespace MueLu {

/*!
  @class AggregationFactory class.
  @brief Factory for coarsening a graph.
*/

template<class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
class AggregationFactory {
#include "MueLu_UseShortNames.hpp"

  private:
      //! aggregation algorithm type
      std::string Algorithm_;

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    AggregationFactory() {}

    //! Destructor.
    virtual ~AggregationFactory() {}
    //@}

    //! @name Build methods.
    //@{

    //! Build aggregates.
    Teuchos::RCP<Aggregates> Build(Graph const &graph, AggregationOptions const &options) const
    {
      Teuchos::OSTab tab(this->out_);

      Teuchos::RCP<Aggregates> aggregates = MueLu_Aggregate_CoarsenUncoupled(options,*graph);
      std::string name = "UC_CleanUp";
      MueLu_AggregateLeftOvers(options, *aggregates, name, *graph);
      return aggregates;
    }
    //@}

    //! @name Set/Get methods.
    //@{
    void SetAlgorithm(std::string const &type) {
      Algorithm_ = type;
    }

    std::string GetAlgorithm() {
      return Algorithm_;
    }
    //@}

}; //class AggregationFactory

} //namespace MueLu

#define MUELU_AGGREGATIONFACTORY_SHORT

#endif //ifndef MUELU_AGGREGATIONFACTORY_HPP
