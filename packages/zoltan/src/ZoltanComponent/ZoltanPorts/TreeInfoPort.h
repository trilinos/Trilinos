/*
  An Port that allows the querying for type of entities, their numbers etc,
  when stored in a tree format.
  Hasn't been debated on.
  Jaideep Ray, SNL, Livermore, 08/20/02
*/

#ifndef TreeInfoPortHSeen
#define TreeInfoPortHSeen

#include "cca.h"

namespace LoadPartitionerSpace
{
  class TreeInfo : public virtual CONV_NS(Port)
  {
    public  :

    TreeInfo() : CONV_NS(Port)() {} 

    virtual ~TreeInfo() {} 

    /// Querying of entities in a mesh
    //@{
    /** Set the type of entity you are querying - could be mesh points,
	mesh cells, particles etc.
	@param name : identifies what entity you want.
    */
    virtual int SetTypeOfElement(char *name) = 0 ;

    /// Returns number of coarse elements that you specified above. This number
    /// is the sum across all processors. Let the answer be M.
    virtual int GetNumberOfCoarseElements() = 0 ;
    
    /// List of coarse elements  on this proc
    /**
       Someone queries GetNumberOfCoarseElements and allocates space for the
	GIDs of all the entities in the domain. This method is then called, and
	returns : 
	<ol>
	    <li> GIDs and LIDs of all entities in the domain (not just this proc).</li>
	    <li> whether these entities are on-proc or not. </li>
	    <li> the number of vertices associated with each entity. On a given proc
	         if 2 elements share the same vertex, the number of that common vertex 
		 should be unique on that proc.</li>
            <li> the vertex list. If entities 0 to (i-1) are stored in 0 to N-1 of the
	         vertex list, the i^th entity's vertices are in N:N+num_vertices[i].
	    <li> whether these vertices need to be traversed in a user defined order
	         or not </li>
	    <li> if these are to be traversed in a user-defined manner, which vertex
	          does one start with, for each entity </li>
	    <li> ditto as above, but which vertex does one exit from.</li>
        </ol>
	@param num_gid_entries : number of elements in a 1D array to specify a GID
	@param num_lid_entries : number of elements in a 1D array to specify a LID
	@param gid : a empty array, M*num_gid_entries long, to be filled up with 
	all GIDs on all procs.
	@param lid : a empty array, M*num_lid_entries long, to be filled up with
	all LIDs on all procs.
	@param assigned : a empty array, M long, to be filled up with whether 
	the entities are on this processor or not.
	@param num_vert : an empty array, M long, to be filled up with the number of 
	vertices associated with each entity.
	@param vertex_list : an empty array, Q = sum_{i=0}^{M} num_vertices[i] long,
	listing all the vertices of all the entities on the coarse mesh on all procs
	combined.
	@param in_order : a single value, denoting if the vertices in the entities
	need to be traversed in a particular user-defined manner.
	@param in_vertex : an empty array, M long, to be filled up with the starting
	vertex (for traversal purposes) for all entities across all procs. Applicable
	only if in_order is true.
	@param out_vertex : an empty array, M long, to be filled up with the ending
	vertex (for traversal purposes) for all entities across all procs. Applicable
	only if in_order is true.
	@return -1 in case of error.
    */
    virtual int CoarseMeshEntitiesList(int num_gid_entries, int num_lid_entries,
				       unsigned int *gids, unsigned int *lids, 
				       bool *assigned, int *num_vertices, 
				       unsigned int *vertex_list, bool *in_order, 
				       unsigned int *in_vertex, 
				       unsigned int *out_vertex) = 0 ;
    

    /// List of coarse mesh entities, one by one
    /** The next two methods help in the creation of a CoarseMeshEntitiesList,
	by being called one-at-a-time, for each entity.
	@param num_gid_entries : number of elements in a 1D array to specify a GID
	@param num_lid_entries : number of elements in a 1D array to specify a LID
	@param prev_gid : A filled array, num_gid_entries long, containing the GID
	of an entity. The callee is expected to furinish the GID of the entity after
	this one. This is applicable only in GetNextCoarseMeshEntity().
	@param prev_lid : same as above, but for an LID.
	@param next_gid : a empty array, num_gid_entries long, to be filled up by the 
	GID of the 1st coarse mesh entity.
	@param next_lid : same as above, but for LID.
	@param assigned : a boolen to be filled up with whether the first coarse mesh
	entity is on this processor or not.
	@param num_vert : a integerto be filled up with the number of vertices 
	associated with the 1st coarse mesh entity.
	@param vertex_list : an empty array, num_vertices long,	listing all the 
	vertices of the 1st coarse mesh entity on the coarse mesh on all procs
	combined.
	@param in_order : a single value, denoting if the vertices in need to be 
	traversed in a particular user-defined manner.
	@param in_vertex : an integer to be filled up with the starting
	vertex (for traversal purposes) for the first coarse mesh entity. 
	Applicable only if in_order is true.
	@param out_vertex : an integer, to be filled up with the ending
	vertex (for traversal purposes) for the 1st coarse mesh entity. Applicable 
	only if in_order is true.
	@return -1  if error. -2 if there are no more entites to be returned
    */
    virtual int GetFirstCoarseMeshEntity(int num_gid_entries, int num_lid_entries,
					 unsigned int *gid, unsigned int *lid, 
					 bool *assigned, int *num_vertices, 
					 unsigned int *vertex_list, bool *in_order, 
					 unsigned int *in_vertex, 
					 unsigned int *out_vertex) = 0 ;
    
    /// The same as above but for the next and the next and the next entities.
    virtual int GetNextCoarseMeshEntity(int num_gid_entries, int num_lid_entries,
					unsigned int *prev_gid, unsigned int *prev_lid, 
					unsigned int *next_gid, unsigned int *next_lid,
					bool *assigned, int *num_vertices, 
					unsigned int *vertex_list, 
					unsigned int *in_vertex, 
					unsigned int *out_vertex) = 0 ;
    //@}

    /// Querying of children of an entity.
    //@{
    /// Get number of children of an entity.
    /** 
	Given a GID and a LID return the number of children an entity has
	@param num_gid_entries : number of entries in a 1D array that specifies a 
	GID.
	@param num_lid_entries : ditto, but for an LID.
	@param gid : the GID of the entity we are interested in; passed in by callee.
	@param lid : ditto, but for an LID.
	@return number of children of this entity. Let this number be L
    */
    virtual int GetNumberOfChildren(int num_gid_entries, int num_lid_entries,
				    unsigned int *gid, unsigned int *lid) = 0 ;

    
    /// Get the children of an entity.
    /**
       Given an entity, get a list of all its children
       @param num_gid_entries : number of entries in a 1D array that specifies a 
       GID.
       @param num_lid_entries : ditto, but for an LID.
       @param gid : GID of the entity whose children we need.
       @param lid : its LID
       @param child_gids : an empty array, L*num_gid_entries long, which will be 
       filled up with the GIDs of the children on return
       @param child_lids : ditto, but lids
       @param assigned : a boolean array, L long, indicating if the children are 
       on-proc or not, on return.
       @param num_vertices : an empty array, L long, containing the number of vertics
       @param child_vertices : upon return, a pointer into an array, allocated
       by the callee, of the vertices of the children.
       @param refinement_type : an integer denoting how the refinement was achieved.
       <ol>
          <li> 1 = refinement by bisecting triangles </li>
	  <li> 2 = refinement by quadrasection of quadrilaterals. </li>
	  <li> 3 = some other refinement technique. </li>
	  <li> 4 = traverse the children in the order they are provided. </li>
       </ol>
       @param in_vertex : an array, L long, filled up with the starting vertex of
       each of children for the entity specified. This is used only if refinement_type
       is 4.
       @param out_vertex : same as above but for exiting vertex.
       @return -1 if an error, 0 otherwise.
    */
    virtual int GetChildrenListOfEntity(int num_gid_entries, int num_lid_entries,
					unsigned int *gid, unsigned int *lid,
					unsigned int *child_gids, 
					unsigned int *child_lids, bool *assigned,
					int *num_vertices, int **child_vertices,
					int *refinement_type, 
					unsigned int *in_vertex, 
					unsigned int *out_vertex) = 0 ;
    
    /// Weight of children of a given entity.
    /**
       @param num_gid_entries : number of entries in a 1D array that specifies a 
       GID.
       @param num_lid_entries : ditto, but for an LID.
       @param gid : GID of the entity whose weight is desired
       @param lid : its LID
       @param wgt_reqd : does the calle need to fill in the weight of this entity
       @param weight : a number, to be filled in by the weight of the entity specified
       by the GID.
       @return -1 on error, 0 otherwise.
    */
    virtual int GetChildrenWeight(int num_gid_entries, int num_lid_entries,
				  unsigned int *gid, unsigned int *lid,
				  bool wgt_reqd, float *weight) = 0 ;

    /// Cleanup ; remove all memory allocated to process children info for this
    /// coarse mesh element.
    virtual int cleanup() = 0 ;
    //@}
  };
};
#endif
