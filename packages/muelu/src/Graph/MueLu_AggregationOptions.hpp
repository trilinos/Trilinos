#ifndef MUELU_AGGREGATIONOPTIONS_HPP
#define MUELU_AGGREGATIONOPTIONS_HPP

/*******************************************************************************
   User-requested options affecting aggregation. 
*******************************************************************************/

namespace MueLu {


  namespace AggOptions {

    /* Options defining how to pick-up the next root node in the local aggregation procedure */
    enum Ordering {
      NATURAL = 0, /* node ordering   */
      RANDOM  = 1, /* random ordering */
      GRAPH   = 2  /* graph ordering  */
    };

  } // namespace AggOptions

  using namespace AggOptions;

  class AggregationOptions : public Teuchos::Describable {
    
  public:
    
    AggregationOptions() {} //TODO: init
    ~AggregationOptions() {}

    inline void SetPrintFlag(double printFlag) { printFlag_ = printFlag; } //TODO: to be removed
    inline void SetOrdering(Ordering ordering)      { ordering_ = ordering;   }
    inline void SetMinNodesPerAggregate(int minNodesPerAggregate) { minNodesPerAggregate_ = minNodesPerAggregate; }
    inline void SetMaxNeighAlreadySelected(int maxNeighAlreadySelected) { maxNeighAlreadySelected_ = maxNeighAlreadySelected; }
    inline void SetPhase3AggCreation(double phase3AggCreation) { phase3AggCreation_ = phase3AggCreation; }
    
    inline double GetPrintFlag() const { return printFlag_; } //TODO: to be removed
    inline Ordering GetOrdering() const     { return ordering_;   }
    inline int GetMinNodesPerAggregate() const { return minNodesPerAggregate_; }
    inline int GetMaxNeighAlreadySelected() const { return maxNeighAlreadySelected_; }
    inline double GetPhase3AggCreation() const { return phase3AggCreation_; }
        
  private:
    double printFlag_;
    Ordering ordering_;  /**<  natural, random, graph           */
    int    minNodesPerAggregate_;      /**<  aggregate size control           */
    int    maxNeighAlreadySelected_;   /**<  complexity control               */
    double phase3AggCreation_;         /**<  Steers how the MIS  and Uncoupled 
                                          handle phase 3 of aggregation.     
                                          Values near 0 create few additional
                                          aggregates.Large values create many
                                          additional aggregates. Convergence 
                                          can be improve by new  
                                          aggregates but nonzero fill-in     
                                          increases on coarse meshes.        
                                          Default: .5                      */
    
  };
  
}

#endif

// TOOD: Default parameters
// TODO: Use a ParameterList instead ?
