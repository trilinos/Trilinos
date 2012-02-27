#ifndef USER_APP_RESPONSE_AGGREGATOR_FACTORY_HPP
#define USER_APP_RESPONSE_AGGREGATOR_FACTORY_HPP

#include "Panzer_ResponseAggregator_Factory.hpp"

// User Defined Builders (in this case we reuse a panzer builder)
#include "Panzer_ResponseAggregator_IPCoordinates.hpp"

namespace user_app {
  
  /** Simple example of a user defined response aggregator factory. */
  template <typename Traits>
  class MyResponseAggregatorFactory : public panzer::ResponseAggregatorFactory<Traits> 
  {
    
    void addResponseTypes(panzer::ResponseAggregator_Manager<Traits>& ram) const
    {


      // Add integration point coordinates
      {
	panzer::ResponseAggregator_IPCoordinates_Builder ipCoordBuilder;
	ipCoordBuilder.setGlobalIndexer(ram.getGlobalIndexer());
	ipCoordBuilder.setLinearObjFactory(ram.getLinearObjFactory());
	ram.defineAggregatorTypeFromBuilder("IP Coordinates",ipCoordBuilder);
      }

      // Add others...

    }
    
  };

}

#endif
