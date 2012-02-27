#ifndef PANZER_RESPONSE_AGGREGATOR_FACTORY_HPP
#define PANZER_RESPONSE_AGGREGATOR_FACTORY_HPP

#include "Panzer_ResponseAggregator_Manager.hpp"

namespace panzer {

  /** \brief Adds user defined aggregator types to the response library

      User defined aggregators are added to response library by adding
      aggregator builders to the ResponseAggregatorManager via the
      method defineAggregatorTypeFromBuilder(...).  This class is a
      pure virtual interface for users to add builders to the
      aggregator manager.
  */
  template <typename Traits>
  class ResponseAggregatorFactory {

  public:

    virtual ~ResponseAggregatorFactory() {}

    /** \brief Registers user defined response aggregator types with a response aggregator manager */
    virtual void addResponseTypes(panzer::ResponseAggregator_Manager<Traits>& ram) const = 0;

  };

}

#endif
