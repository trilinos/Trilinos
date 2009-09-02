/*
  A port for a distributed data directory, of help to people trying
  to migrate data after a new partition has been proposed.

  Has not been debated on.
  Jaideep Ray, SNL, Livermore, CA 08/21/02
*/

#ifndef DistDataDirHSeen
#define  DistDataDirHSeen

#include "cca.h"
#include "mpi.h"
#include "EntityList.h"

// namespace conversion
#include "CONV_NS.h"

namespace LoadPartitionerSpace 
{
  class DistDataDirectoryPort : public virtual CONV_NS(Port)
  {
    
    public :

    DistDataDirectoryPort() : CONV_NS(Port)() {}

    virtual ~DistDataDirectoryPort() {}

    //@{
    /// Create/destroy a directory
    /** Create or destroy a distributed directory of elements identified
	by their global ids.
	@param no_entries_in_gid : no of ints needed to represent a gid
	@param no_entried_in_lid : ditto for local ids
	@param user_def_field_length : field length per entity 
	@param hash_table_length : 
	@return -1 : something went wrong and the directory was not made.
    */
    virtual int CreateDistDir(MPI_Comm *comm, int no_entries_in_gid,
			      int no_entries_in_lid, int user_def_field_length, 
			      int hash_table_length) = 0 ;

    /// Destroy a directory
    virtual void DestroyDir() = 0 ;
    //@}

    //@{
    /// Put/get these entities
    /** Put the entities in the EntityList in the directory 
	@param NewEntities : List of entities to be shoved into the directory
	@return -1 if there's an error, else 0.
    */
    virtual int UpdateDirWithNewEntityList( EntityList *NewEntities) = 0 ;

    /** Find the entities in the EntityList from the directory.
	You need to allocate the space for the local IDs and user defined info;
	I'll fill them up.
	@param FindTheseEntities : List of entities to be found in the directory
	@return -1 if there's an error, else 0.
    */
    virtual int FindDirEntries( EntityList *FindTheseEntities) = 0 ;

    /** Remove the entities in the list from the directory
	@param RemoveTheseEntities : list of entities to be removed from the dir.
	@return -1 in case of error, else 0
    */
    virtual int RemoveDirEntries( EntityList *RemoveTheseEntities ) = 0 ;
    //@}

    /// Print all the directory entries
    virtual void PrintDirEntries() = 0 ;

    /// Print statistics re directory
    virtual void PrintStats() = 0 ;
  } ;
};
#endif
