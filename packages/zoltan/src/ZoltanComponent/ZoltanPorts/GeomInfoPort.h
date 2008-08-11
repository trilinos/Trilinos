/*
  An Port that allows the querying for type of geometry of an entity
  Hasn't been debated on.
  Jaideep Ray, SNL, Livermore, 08/20/02
*/

#ifndef GeomInfoPortHSeen
#define GeomInfoPortHSeen

#include "cca.h"

// namespace conversion
#include "CONV_NS.h"

namespace LoadPartitionerSpace
{
  class GeomInfo : public virtual CONV_NS(Port)
  {
    public  :

    GeomInfo() : CONV_NS(Port)() {} 

    virtual ~GeomInfo() {} 

    /// Querying of geometries in a mesh
    //@{
    /** Set the type of entity you are querying - could be mesh points,
	mesh cells, particles etc.
	@param name : identifies what entity you want.
    */
    virtual int SetTypeOfGeom(char *name) = 0 ;

    /// Returns number of dimensions this particular geometry is . 
    /**
       Essentially telling the caller that the geometry is 1D/2D/3D and so
       it will need a 1D array 1/2/3 long to hold its "description", usually 
       position. Let the answer be Q.
    */
    virtual int GetNumberOfDims() = 0 ;

    /// Get geometry info for an entity. 
    /** 
	This method is a means of getting the "description", usually
	position of an entity of geometry as set by SetTypeOfGeom().
	@param num_gid_entries : number of elements in a 1D array to specify a GID
	@param num_lid_entries : number of elements in a 1D array to specify a LID
	@param gid : GID of the entity the caller is interested in; an array 
	num_gid_entries long.
	@param lid : ditto, but its LID.
	@param geom_vec : a double array, Q long, to be filled up by the 
	callee with the "description" of the entity, usually its position.
	@return -1 in case of error.
    */
    virtual int GetGeomDescription(int num_gid_entries, int num_lid_entries,
				   unsigned int *gids, unsigned int *lids, 
				   double *geom_vec) = 0 ;
    
    //@}
  };
};
#endif
