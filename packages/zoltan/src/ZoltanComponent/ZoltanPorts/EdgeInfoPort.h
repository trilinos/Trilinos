/*
  An Port that allows the querying for type of edges, their numbers etc.
  Hasn't been debated on.
  Jaideep Ray, SNL, Livermore, 08/20/02
*/

#ifndef EdgeInfoPortHSeen
#define EdgeInfoPortHSeen

#include "cca.h"

// namespace conversion
#include "CONV_NS.h"

namespace LoadPartitionerSpace
{
  class EdgeInfo : public virtual CONV_NS(Port)
  {
    public  :

    EdgeInfo() : CONV_NS(Port)() {} 

    virtual ~EdgeInfo() {} 

    /// Querying of edges in a mesh entity
    //@{
    /** Set the type of entity you are querying - could be mesh points,
	mesh cells, particles etc.
	@param name : identifies what entity you want.
    */
    virtual int SetTypeOfEdge(char *name) = 0 ;

    /// Returns number of edges for a specified entity
    /** An entity, specified by a GID and LID, will need to communicate with
	other neighboring entities - these connections are edges. Given a 
	GID, return me the number of connections it has.
	@param num_gid_entries : number of elements in a 1D array to specify a GID
	@param num_lid_entries : number of elements in a 1D array to specify a LID
	@param gid : an array with the GID in it; an array num_gid_entries long
	@param lid : an array with the LID in it; an array num_lid_entries long.
	@return number of edges from the conectivity graph. (N)
    */
    virtual int GetNumberOfEdges(int num_gid_entries, int num_lid_entries,
				 unsigned int *gids, unsigned int *lids) = 0 ;

    /// Given an object, find its neighbors with whom it has connections
    /** Given an entity, specified by a GID and LID, determine the neighbors
	it is connected to, the processors on which these neighbors live and 
	the edge weights for the edges terminating at this entity.
	@param num_gid_entries : number of elements in a 1D array to specify a GID
	@param num_lid_entries : number of elements in a 1D array to specify a LID
	@param gid : an array with the GID in it.
	@param lid : an array with the LID in it.
	@param neighor_gids : an empty array, N long, to be filled with the neighboring
	entities' GIDs.
	@param neighor_procs : an empty array, N long, to be filled with the 
	neighboring entities' processors.
	@param wgt_reqd : does the caller require you to fill in edge weights ?
	@param weights : an empty array, N long, to be filled up by all the edgess' 
	weights.
	@return -1 in case of error.
    */
    virtual int EdgesListForThisEntity(int num_gid_entries, int num_lid_entries,
				       unsigned int *gid, unsigned int *lid, 
				       unsigned int *neighbor_gids, 
				       int *neighbor_procs, 
				       bool wgt_reqd, float *edge_weights) = 0 ;
    
    //@}
  };
};
#endif
