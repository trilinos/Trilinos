#ifndef MUELU_AGGOPTIONS_HPP
#define MUELU_AGGOPTIONS_HPP

/*******************************************************************************
   User-requested options affecting aggregation. 
*******************************************************************************/
typedef struct MueLoo_AggOptions_Struct
{
   double print_flag;
   int    ordering;                  /**<  natural, random, graph           */
   int    min_nodes_per_aggregate;   /**<  aggregate size control           */
   int    max_neigh_already_selected;/**<  complexity control               */
   double phase3_agg_creation;       /**<  Steers how the MIS  and Uncoupled 
                                           handle phase 3 of aggregation.     
                                           Values near 0 create few additional
                                           aggregates.Large values create many
                                           additional aggregates. Convergence 
                                           can be improve convergence by new  
                                           aggregates but nonzero fill-in     
                                           increases on coarse meshes.        
                                           Default: .5                      */
} MueLoo_AggOptions;

#endif
