#ifndef MUELU_AGGOPTIONS_HPP
#define MUELU_AGGOPTIONS_HPP

/*******************************************************************************
   User-requested options affecting aggregation. 
*******************************************************************************/
typedef struct MueLu_AggOptions_Struct
{
   double printFlag;
   int    ordering;                  /**<  natural, random, graph           */
   int    minNodesPerAggregate;      /**<  aggregate size control           */
   int    maxNeighAlreadySelected;   /**<  complexity control               */
   double phase3AggCreation;         /**<  Steers how the MIS  and Uncoupled 
                                           handle phase 3 of aggregation.     
                                           Values near 0 create few additional
                                           aggregates.Large values create many
                                           additional aggregates. Convergence 
                                           can be improve by new  
                                           aggregates but nonzero fill-in     
                                           increases on coarse meshes.        
                                           Default: .5                      */
} MueLu_AggOptions;

#endif
