/*
  This is the partitioner factory. On being called with the right inputs,
  it creates the required partitioner.

  Jaideep Ray, SNL, 08/23/02
*/

#ifndef PartitionerFactoryHSeen
#define  PartitionerFactoryHSeen

#include "cca.h"
#include "stdPorts.h"
#include "mpi.h"

// Ports I provide ;
#include "LoadBalancer.h"
#include "LoadBalancerFactoryPort.h"
#include "DataLocationPort.h"
#include "EntityList.h"

// Ports I use
#include "EntityInfoPort.h"
#include "GeomInfoPort.h"
#include "EdgeInfoPort.h"
#include "TreeInfoPort.h"
#include "DataMigrationPort.h"


// An object I use
#include "EntityListImpl.h"

//C++ STL
#include <vector>

// Zoltan
#include "include/zoltan.h"

//Sandia's own definition of the (user-) configurable parameters port
#include "jc++/jc++.h"
#include "jc++/util/jc++util.h"
#include "parameters/parametersStar.h"
#include "port/supportInterfaces.h"
#include "dc/port/DefaultParameterPort.h"

// namespace change
#include "CONV_NS.h"

namespace ZoltanSpace
{
  class PartitionerFactory_JR : 
    public virtual ::LoadPartitionerSpace::LoadBalancerFactory,
    public virtual ::LoadPartitionerSpace::DataLocationPort,
    public virtual CONV_NS(Component)
  {
    public :

    PartitionerFactory_JR() ;

    virtual ~PartitionerFactory_JR() ;

    virtual void setServices( CONV_NS(Services) *svc) ;

    /* 
       create a partitioner over the parocs in the communicator A. The good
       values of name are (case-sensitive!)
       "RCB" : Recursive coordinate bisection
       "RIB" : Recursive inertial bisection
       "HSFC" : Hilbert space-filling curve
       "REFTREE" : Refinement-tree
       "PARMETIS" : One of ParMetis's partitioners
          - ParMetis has a few partitioners. You specify which particular
	    one once the Parmetis uber-partitioner exists.
       "JOSTLE" : One of Jostle's partitioners.
          - Not currently supported; needs a license.
       "OCTPART" : Octree-based partitioning.
    */
    ::LoadPartitionerSpace::LoadBalancer *CreateLB_VM_pair(MPI_Comm *A, char *name) ;

    // Destrou this load-balancer.
    virtual void Destroy_LB_VM_pair( ::LoadPartitionerSpace::LoadBalancer *p)  ;

    // -----------------------------------------------------------------------
    
    // Implementation of the DataLocation Port
    
    // Given a list of entitities coming to me, I ask this partitioner to
    // compute which of my entities are needed by other procs. The partitioner
    // factory will allocate an EntitityList for the entities going off elsewhere and
    // will set a pointer (that the caller supplies) to this EntityList
    int ComputeDestinations( ::LoadPartitionerSpace::LoadBalancer *A_lb, 
		      ::LoadPartitionerSpace::EntityList *OffProcEntitiesComingToMe,
		      ::LoadPartitionerSpace::EntityList **MyEntitiesMigratingElsewhere) ;

    // Given a list of my entities going elsewhere, i ask the partitioner
    // which off-proc entities are coming to me. The partitioner factory will
    // allocate an list of entities coming to the caller and will set a pointer
    // supplied by the caller to this list.
    int ComputeSources( ::LoadPartitionerSpace::LoadBalancer *A_lb, 
		        ::LoadPartitionerSpace::EntityList *MyEntitiesMigratingElsewhere,
			::LoadPartitionerSpace::EntityList **OffProcEntitiesComingToMe) ;


    // Given a list of entities I am exporting off and a list that I am
    // importing, go do this migration.
    int MigrateEntities( ::LoadPartitionerSpace::LoadBalancer *A, 
			 ::LoadPartitionerSpace::EntityList *MyEntitiesMigratingElsewhere,
			 ::LoadPartitionerSpace::EntityList *OffProcEntitiesComingToMe) ;

    // Needs to be called after ComputeSources() and Compute Destination are
    // called to de-allocate the EntityLists that are allocated.
    void Cleanup( ::LoadPartitionerSpace::LoadBalancer *lb) ;

    // ------------------------------------------------------------------------
    // Remove this LB from my list
    void removeLB(int i) { my_lbs[i] = 0 ; }

    int SetZoltanQueriesEntityInfoPort(struct Zoltan_Struct *zz)  ;
    int SetZoltanQueriesGeomInfoPort(struct Zoltan_Struct *zz) ;
    int SetZoltanQueriesEdgeInfoPort(struct Zoltan_Struct *zz) ;
    int SetZoltanQueriesTreeInfoPort(struct Zoltan_Struct *zz) ;
    int SetZoltanQueriesDataMigrationPort(struct Zoltan_Struct *zz) ;

    private :

    void init() ;
    bool is_init ;
    CONV_NS(Services) *pSvc;

    ::LoadPartitionerSpace::EntityInfo *pEntityInfo ;
    ::LoadPartitionerSpace::GeomInfo *pGeomInfo ;
    ::LoadPartitionerSpace::EdgeInfo *pEdgeInfo ;
    ::LoadPartitionerSpace::TreeInfo *pTreeInfo ;
    ::LoadPartitionerSpace::DataMigrationPort *pDataMig ;

    ::std::vector< ::LoadPartitionerSpace::LoadBalancer * > my_lbs ;
    
    BoolParameter *pIterate ;
    ConfigurableParameterPort *pp;

    // Temporary EntityLists
    EntityListImpl *temp_incoming, *temp_outgoing ;
  } ;
} ;
#endif
