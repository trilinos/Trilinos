/*
  An Port that allows the querying for type of entities, their numbers etc.
  Hasn't been debated on.
  Jaideep Ray, SNL, Livermore, 08/20/02
*/

#ifndef EntityInfoPortHSeen
#define EntityInfoPortHSeen

#include "cca.h"

// namespace conversion
#include "CONV_NS.h"

namespace LoadPartitionerSpace
{
  class EntityInfo : public virtual CONV_NS(Port)
  {
    public  :

    EntityInfo() : CONV_NS(Port)() {} 

    virtual ~EntityInfo() {} 

    /// Querying of entities in a mesh
    //@{
    /** Set the type of entity you are querying - could be mesh points,
	mesh cells, particles etc.
	@param name : identifies what entity you want.
    */
    virtual int SetTypeOfEntity(char *name) = 0 ;

    /// Returns number of entities that you specified above that exist on this 
    /// processor. Let this be L.
    virtual int GetNumberOfEntities() = 0 ;

    /// List of Entities on this proc
    /** Someone queries GetNumberOfEntities and allocates space for the
	GIDs of all the entities on this proc. On passing in the # of
	elements in a 1D array that specify a GID or LID, a call on this
	method with those empty array will fill them up with the GIDs and
	LIDs of all entities of the type specified.
	@param num_gid_entries : number of elements in a 1D array to specify a GID
	@param num_lid_entries : number of elements in a 1D array to specify a LID
	@param gid : a empty array, L*num_gid_entrie long, to be filled up by all 
	GIDs on this proc.
	@param lid : a empty array, L*num_lid_entries long, to be filled up by all 
	LIDs on this proc.
	@param wgt_reqd : is the component required to fill in weights of the entities
	too.
	@param weights : an empty array, L long,  to be filled up by all the 
	entities' weights. This is done only if wgt_reqd is true.
	@return -1 in case of error.
    */
    virtual int EntitiesOnThisProc(int num_gid_entries, int num_lid_entries,
				   unsigned int *gids, unsigned int *lids, 
				   bool wgt_reqd, float *weights) = 0 ;
    
    /// List of Entities on this proc, created by iterating
    /** 
	This is another means of creating an entities' list, as done by
	the method above. The caller calls the next 2 methods, over and over again,
	which return an entity specification and its weight.
	@param num_gid_entries : number of elements in a 1D array to specify a GID
	@param num_lid_entries : number of elements in a 1D array to specify a LID
	@param first_gid : a empty array, num_gid_entries long, to be filled up by 
	the GID of the first entity  on this proc.
	@param first_lid : a empty array, num_lig_entries long, to be filled up by 
	the LIDs of the 1st entity on this proc.
	@param wgt_reqd : true or false - does the caller require the weight of the
	entity too ?
	@param first_weights : if the previous is true, a number to be filled up in 
	by  the 1st entity's weight.
	@param prev_gid : an array, num_gid_entries long, filled with the GID of an
	entity. The callee is expected to furnish the GID of the entity after this
	one in next_gid
	@param prev_lid : Same as above, but for an LID.
	@param next_gid : an empty array, num_gid_entries long, to be filled by the
	callee with the GID of the entity occuring after the entity specified by
	prev_gid in the entity list.
	@param next_lid : same as above, but for LID.
	@return -1 in case of error. -2 if there are no more entities to be returned.
    */
    virtual int GetFirstEntity(int num_gid_entries, int num_lid_entries,
			       unsigned int *first_gid, unsigned int *first_lid,
			       bool wgt_reqd, float *first_weight) = 0 ;

    /// Same as above but for the next entity.
    virtual int GetNextEntity(int num_gid_entries, int num_lid_entries,
			       unsigned int *prev_gid, unsigned int *prev_lid,
			       unsigned int *next_gid, unsigned int *next_lid,
			       bool wgt_reqd, float *next_weight) = 0 ;
    //@}
  };
};
#endif
