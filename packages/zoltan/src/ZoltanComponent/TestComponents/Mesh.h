/*
  The header for the Mesh component. This will have the ability to read
  Chaco & Nemesis mesh files in ../ZoltanTests. It will make extensive use
  of the code in ../../Zoltan/driver and ../../Zoltan/ch
*/

#ifndef MeshHSeen
#define MeshHSeen

// CCA standard definition of Port and Component
#include "cca.h"
#include "stdPorts.h"

//Ports I present
#include "EntityInfoPort.h"
#include "GeomInfoPort.h"
#include "EdgeInfoPort.h"
#include "IOPort.h"
#include "DataMigrationPort.h"

// C++ 
#include <iostream>
#include <string>

// Specific to the driver that Zoltan has in ../../Zoltan/zdrive/
// Contains the definition of MESH_INFO, PARIO_INFO, PROB_INFO
#include "dr_const.h"
#include "dr_input_const.h"

// namespace conversion
#include "CONV_NS.h"

namespace ZoltanTestSpace
{
  class Mesh : public virtual CONV_NS(Component), 
               public virtual ::LoadPartitionerSpace::EntityInfo,
               public virtual ::LoadPartitionerSpace::GeomInfo,
               public virtual ::LoadPartitionerSpace::EdgeInfo,
               public virtual ::LoadPartitionerSpace::DataMigrationPort,
               public virtual ::ZoltanTestSpace::IOPort
  {
    public :

    Mesh() ;

    virtual ~Mesh() ;

    virtual void setServices( CONV_NS(Services) *svc) ;

    /*****************************************************************/

    /* Methods in ::LoadPartitionerSpace::EntityInfo */

    /*
      What kind of entity do you want ? Vertex, nodes etc. This particular
      implementation only supports vertices by default, and nothing else
    */
    virtual int SetTypeOfEntity(char *name) ;

    // How many entities exist on this proc ?
    virtual int GetNumberOfEntities() ;

    // Now that you know how many entities exist on this proc, get their
    // GIDs and LIDs.
    virtual int EntitiesOnThisProc(int num_gid_entries, int num_lid_entries,
				   unsigned int *gids, unsigned int *lids,
				   bool wgt_reqd, float *weights) ;

    // Accessing these entities via iterator
    virtual int GetFirstEntity(int num_gid_entries, int num_lid_entries,
			       unsigned int *first_gid, 
			       unsigned int *first_lid, bool wgt_reqd, 
			       float *first_weight) ;

    virtual int GetNextEntity(int num_gid_entries, int num_lid_entries,
			      unsigned int *prev_gid, unsigned int *prev_lid,
			      unsigned int *next_gid, unsigned int *next_lid,
			      bool wgt_reqd, float *next_weight) ;

    /******************************************************************/

    /* Methods in ::LoadPartitionerSpace::GeomInfoPort */

    // Set the type of Geometry you want (rectangles, triangles etc)
    // We only support the default
    virtual int SetTypeOfGeom(char *name) ;

    // Returns whether 1D/2D/3D
    virtual int GetNumberOfDims() ;

    // Get the geometry description
    virtual int GetGeomDescription(int num_gid_entries, int num_lid_entries,
				   unsigned int *gids, unsigned int *lids,
				   double *geom_vec) ;

    /*******************************************************************/

    /* Methods in ::LoadPartitionerSpace::EdgeInfoPort */
    
    // What type of Edge ? We only support the default.
    virtual int SetTypeOfEdge(char *name) ;

    // Get number of edges for the Entity specified by this GID or LID.
    virtual int GetNumberOfEdges(int num_gid_entries, int num_lid_entries,
				 unsigned int *gids, unsigned int *lids) ;

    // Get the edge list for the entity specified by the GID and/or LID
    virtual int EdgesListForThisEntity(int num_gid_entries, int num_lid_entries,
				       unsigned int *gid, unsigned int *lid,
				       unsigned int *neighbor_gids, 
				       int *neighbor_procs, bool wgt_reqd,
				       float *edge_weights) ;

    /********************************************************************/
    // From ::LoadPartionerSpace::DataMigrationPort

    virtual int DataBufferForThisEntity(int num_gid_entries, int num_lid_entries,
					unsigned int *gid, unsigned int *lid) ;

    virtual int PrePackingProcessing( 
                       ::LoadPartitionerSpace::EntityList *MyEntitiesGoingOut,
                       ::LoadPartitionerSpace::EntityList *OthersEntitiesComingIn) ;

    virtual int PostPackingProcessing( 
                       ::LoadPartitionerSpace::EntityList *MyEntitiesGoingOut,
                       ::LoadPartitionerSpace::EntityList *OthersEntitiesComingIn) ;
 
    virtual int PostUnpackingProcessing( 
                        ::LoadPartitionerSpace::EntityList *MyEntitiesGoingOut,
                        ::LoadPartitionerSpace::EntityList *OthersEntitiesComingIn) ;

    virtual int PackThisEntity(int num_gid_entries, int num_lid_entries,
                               unsigned int *gid, unsigned int *lid, 
                               int dest_proc, int buffer_size, char *buffer) ;
    
    virtual int UnpackThisEntity(int num_gid_entries, unsigned int *gid, 
                                 int buffer_size, char *buffer) ;
 
    /********************************************************************/
    // From ZoltanTestSpace::IOPort
    virtual int ProduceGnuplotOutput(char *tag) ;

    virtual int ProduceFileOutput(char *tag) ;

    virtual int DeleteAllAndRestart() ;

    virtual int SetCmdFilename(char *) ;
    /*******************************************************************/
    private :

    void init() ;
    void read_mesh() ;
    void initialize_mesh() ;

    MESH_INFO mesh ;
    PARIO_INFO pio_info ;
    PROB_INFO prob ;

    bool is_init, cmd_filename_set ;
    int Proc, Num_Procs ;
    std::string cmd_filename ;
    CONV_NS(Services) *pSvc ;
  } ;
} ;
#endif
      
    
