/*
  A port that packs up objects to be migrated, unpacks them, and preprocesses/
  postproceses them

  as not been debated on

  J. Ray, SNL, Livermore, 09/21/02
*/

#ifndef DataMigrationPortHSeen
#define DataMigrationPortHSeen

#include "cca.h"
#include "EntityList.h"
#include "LoadBalancer.h"

// namespace conversion
#include "CONV_NS.h"

namespace LoadPartitionerSpace
{
  class DataMigrationPort  : public virtual CONV_NS(Port)
  {

    public :

    DataMigrationPort()  : CONV_NS(Port)() {} 

    virtual ~DataMigrationPort() {}

    /// Entity specific methods
    /** These methods perform operations on an entity-by-entity basis,
	usually to pack/unpack them.
    */
    //@{
    /// Size of a given object
    /** Given an entity, identified by its local and global id, how much 
	memory does it take ? Answer in bytes.
	@param num_gid_entries : number of entries in a 1D array making up a GID.
	@param num_lid_entries : number of entries in a 1D array making up a LID.
	@param gid : the element we are talking about
	@param lid : its local id
	@return : the memory this enetity needs to be stored, in bytes
    */
    virtual int DataBufferForThisEntity(int num_gid_entries, 
					int num_lid_entries,
					 unsigned int *gid,
					 unsigned int *lid) = 0 ;

    /// Pack the given object into the given buffer
    /** @param buffer_size : the size of the buffer supplied
	@param buffer : place to put the packed object
	@param dest_proc : where this will be sent.
	return -1 : if the buffer returned is too small.
    */
    virtual int PackThisEntity(int num_gid_entries, int num_lid_entries,
			       unsigned int *gid, unsigned int *lid, 
			       int dest_proc, int buffer_size, char *buffer) = 0 ;

    /// Unpack object
    /** @param buffer : memory containing the packed object received from
	some other proc. Unpack from here, and shove into your data structure
	as per gid and lid
	@return -1 if the buffer is too small to be holding the data.
    */
    virtual int UnpackThisEntity(int num_gid_entries, unsigned int *gid, 
				 int buffer_size, char *buffer) = 0 ;
    //@}

    ///Pre, mid- and post-migration processing
    /** The sequence is
	Preprocess - pack entities - PostPacking - transfer and unpack entities -
	post-unpacking-processing
    */
    //@{
    /** This method is invoked, by the component orchestrating the moving
	of data first. At this point the orchestrator on proc P knows which 
	entities are expected to come to P and which of P's current entities
	will move out. Using this info, it calls on a component providing this
	port to pre-process the domain - perhaps collate all field together, set up
	MPI_Structs etc. This is also the point where one traverses the grid to
	see whose neighbors are being shunted off and to mark the (to-be) off-processors
	neighbor's proc_id. No packing of entities has occured yet. After this method,
	the packing of entities start.
	@param MyEntitiesGoingOut : my entities leaving
	@param OthersEntitiesComingIn : entites coming from other procs
	@return -1 if there is an error, else 0
    */
    virtual int PrePackingProcessing( EntityList *MyEntitiesGoingOut,
				      EntityList *OthersEntitiesComingIn) = 0 ;

    /** This is a method called after the entities to be exported have been packed.
	Do any book-keeping here. After this method has been called, the communications
	start and the data gets moved elsewhere and unpacking starts.
	@return -1 : Error.
    */
    virtual int PostPackingProcessing( EntityList *MyEntitiesGoingOut,
				       EntityList *OthersEntitiesComingIn) = 0 ;
				 
    /** This is called after unpacking. Might be a good place to shove unpacked data
	into your data-structure.
	@return -1 : Error.
    */
    virtual int PostUnpackingProcessing( EntityList *MyEntitiesGoingOut,
					 EntityList *OthersEntitiesComingIn) = 0 ;
  };
};
#endif
						 
