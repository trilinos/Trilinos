#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Tuple.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_Utilities.hpp"

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/FieldData.hpp>

namespace panzer {

  /** This returns all the global IDs and coordinates for 
    * a particular side. By "all" that means across all processors.
    */
  std::pair<Teuchos::RCP<std::vector<std::size_t> >,
            Teuchos::RCP<std::vector<Teuchos::Tuple<double,3> > > >
  getSideIdsAndCoords(const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
                      const std::string & sideName);

  /** This returns the locally owned global IDs and coordinates for 
    * a particular side. 
    */
  std::pair<Teuchos::RCP<std::vector<std::size_t> >,
            Teuchos::RCP<std::vector<Teuchos::Tuple<double,3> > > >
  getLocalSideIdsAndCoords(const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
                           const std::string & sideName);

  /** This returns the locally resident (includes ghosted) global IDs
    * for a particular side. 
    */
  Teuchos::RCP<std::vector<std::size_t> >
  getLocalSideIds(const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
                  const std::string & sideName);

  /** Determine a map from the specified side to the set of coordinates
    * and Ids passed in. A vector of pairs that maps from (passed in gids)->(locally owned gids)
    * is returned.
    */
  template <typename MatchObject>
  Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >
  getLocallyMatchedSideIds(const std::vector<std::size_t> & side_ids,
                           const std::vector<Teuchos::Tuple<double,3> > & side_coords,
                           const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
                           const std::string & sideName,const MatchObject & matcher);

  template <typename MatchObject>
  Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >
  matchPeriodicSides(const std::string & left,const std::string & right,
                     const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
                     const MatchObject & matcher);

   
  /** Builds a vector of local ids and their matching global indices.
    * This requires a previously discovered vector of pairs of locally matched
    * ids to distribute. This vector comes from the getLocallyMatchedSideIds.
    */
  Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >
  getGlobalPairing(const std::vector<std::size_t> & locallyRequiredIds,
                   const std::vector<std::pair<std::size_t,std::size_t> > & locallyMatchedIds,
                   const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,bool failure);

  TEUCHOS_UNIT_TEST(periodic_bcs, sorted_permutation)
  {
     std::vector<double> vec(5.0);
     std::vector<std::size_t> permute;
     vec[0] = 0.0; 
     vec[1] = 4.0; 
     vec[2] = 2.0; 
     vec[3] = 3.0; 
     vec[4] = 1.0; 

     panzer_stk::sorted_permutation(vec,permute);

     TEST_EQUALITY(permute.size(),5);
     for(std::size_t i=0;i<permute.size();i++)  
        TEST_EQUALITY(vec[permute[i]],(double) i);
  }

  TEUCHOS_UNIT_TEST(periodic_bcs, getSideIdsAndCoords)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",6);
       pl->set("Y Elements",4);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    // setup DOF manager
    /////////////////////////////////////////////
    const Teuchos::RCP<panzer::ConnManager<int,int> > conn_manager 
           = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));
 
    // run tests
    /////////////////////////////////////////////
    std::pair<RCP<std::vector<std::size_t> >,
              RCP<std::vector<Tuple<double,3> > > > idsAndCoords = getSideIdsAndCoords(mesh,"left");
    std::vector<std::size_t> & sideIds = *idsAndCoords.first;
    std::vector<Tuple<double,3> > & sideCoords = *idsAndCoords.second;

    TEST_EQUALITY(sideIds.size(),5);

    TEST_EQUALITY(sideCoords.size(),5);

    std::vector<std::size_t> permute;
    panzer_stk::sorted_permutation(sideIds,permute);
    for(std::size_t i=0;i<permute.size();i++) {
       std::size_t p = permute[i];

       TEST_EQUALITY(sideIds[p]-1,i*13);

       TEST_FLOATING_EQUALITY(sideCoords[p][0],0.0,1e-14);
       TEST_FLOATING_EQUALITY(sideCoords[p][1],i*1.0/4.0,1e-14);
    }

  }
 
  struct CoordMatcher {
     double error_;
     int index_;
     CoordMatcher(int index) : error_(1e-8),index_(index) {}
     CoordMatcher(int index,double error) : error_(error),index_(index) {}

     bool operator()(const Teuchos::Tuple<double,3> & a,
                     const Teuchos::Tuple<double,3> & b) const
     { return std::fabs(a[index_]-b[index_])<error_; /* I'm being lazy here! */ }
  };

  TEUCHOS_UNIT_TEST(periodic_bcs, getLocallyMatchedSideIds)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    int rank = Comm.MyPID(); 
    int procCnt = Comm.NumProc(); 

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",6);
       pl->set("Y Elements",4);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    // setup DOF manager
    /////////////////////////////////////////////
    const Teuchos::RCP<panzer::ConnManager<int,int> > conn_manager 
           = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));
 
    // run tests
    /////////////////////////////////////////////
    std::pair<RCP<std::vector<std::size_t> >,
              RCP<std::vector<Tuple<double,3> > > > idsAndCoords = getSideIdsAndCoords(mesh,"left");
    std::vector<std::size_t> & sideIds = *idsAndCoords.first;
    std::vector<Tuple<double,3> > & sideCoords = *idsAndCoords.second;

    CoordMatcher matcher(1);
    Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > matchedIds 
          = getLocallyMatchedSideIds(sideIds,sideCoords,mesh,"right",matcher);

    if(rank==procCnt-1) {
       for(std::size_t i=0;i<matchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*matchedIds)[i];
          TEST_EQUALITY(pair.first,pair.second-12);
       }
    }
    else {
       TEST_EQUALITY(matchedIds->size(),0)
    }
  }

  TEUCHOS_UNIT_TEST(periodic_bcs, getGlobalPairing)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    int rank = Comm.MyPID(); 

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",6);
       pl->set("Y Elements",4);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    // setup DOF manager
    /////////////////////////////////////////////
    const Teuchos::RCP<panzer::ConnManager<int,int> > conn_manager 
           = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));
 
    // run tests
    /////////////////////////////////////////////
    Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > locallyMatchedIds;
    {
       // next line requires global communication
       std::pair<RCP<std::vector<std::size_t> >,
                 RCP<std::vector<Tuple<double,3> > > > idsAndCoords = getSideIdsAndCoords(mesh,"left");
       std::vector<std::size_t> & sideIds = *idsAndCoords.first;
       std::vector<Tuple<double,3> > & sideCoords = *idsAndCoords.second;

       CoordMatcher matcher(1);
       locallyMatchedIds = getLocallyMatchedSideIds(sideIds,sideCoords,mesh,"right",matcher);
    }

    Teuchos::RCP<std::vector<std::size_t> > locallyRequiredIds = getLocalSideIds(mesh,"left");

    // next line requires communication
    Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds
          = getGlobalPairing(*locallyRequiredIds,*locallyMatchedIds,mesh,false);

    if(rank==0) {
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
          TEST_EQUALITY(pair.first,pair.second-12);
       }
    }
    else {
       TEST_EQUALITY(globallyMatchedIds->size(),0);
    }
  }

  TEUCHOS_UNIT_TEST(periodic_bcs, matchPeriodicSides)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    int rank = Comm.MyPID(); 

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",6);
       pl->set("Y Elements",4);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    // run tests
    {
       CoordMatcher matcher(1);
       Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds
             = matchPeriodicSides("left","right",mesh,matcher);
   
       // match left & right sides
       if(rank==0) {
          for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
             std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
             TEST_EQUALITY(pair.first,pair.second-12);
          }
       }
       else
          TEST_EQUALITY(globallyMatchedIds->size(),0);
    }

    {
       CoordMatcher matcher(0);
       Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds
             = matchPeriodicSides("top","bottom",mesh,matcher);

       Teuchos::RCP<std::vector<std::size_t> > locallyRequiredIds = getLocalSideIds(mesh,"top");

       TEST_EQUALITY(globallyMatchedIds->size(),locallyRequiredIds->size()); 
   
       // match top & bottom sides
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
          TEST_EQUALITY(pair.first,pair.second+52);
       }
    }

    // test a failure case!
    {
       CoordMatcher matcherX(0);
       CoordMatcher matcherY(1);
       
    //    TEST_THROW(matchPeriodicSides("left","bottom",mesh,matcherX),std::logic_error);
    //    TEST_THROW(matchPeriodicSides("top","right",mesh,matcherY),std::logic_error);
       TEST_THROW(matchPeriodicSides("top","right",mesh,matcherX),std::logic_error);
    //    TEST_THROW(matchPeriodicSides("bottom","left",mesh,matcherY),std::logic_error);
    }
  }

  template <typename MatchObject>
  Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >
  matchPeriodicSides(const std::string & left,const std::string & right,
                     const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
                     const MatchObject & matcher)
  {
    //
    // Overview:
    // A three step algorithm
    //    1. Figure out all nodes and their coordinates that live on the "left"
    //       - Distribute information globally
    //    2. Match the global nodes on the "left" to locally owned nodes on the "right"
    //       - only local work required
    //    3. If a processor requires a node on the left (if it is owned or ghosted)
    //       communicate matching conditions from the right boundary 
    //
    // Note: The matching check could definitely be spead up with a sorting operation
    // Note: The communication could be done in a way that requires less global communication
    //       Essentially doing step one in a Many-2-Many way as opposed to an All-2-All
    //

    using Teuchos::Tuple;
    using Teuchos::RCP;
    using Teuchos::rcp;

    // First step is to globally distribute all the node ids and coordinates
    // on the left hand side: requires All-2-All!
    /////////////////////////////////////////////////////////////////////////
    std::pair<RCP<std::vector<std::size_t> >,
              RCP<std::vector<Tuple<double,3> > > > idsAndCoords = getSideIdsAndCoords(mesh,left);
    std::vector<std::size_t> & sideIds = *idsAndCoords.first;
    std::vector<Tuple<double,3> > & sideCoords = *idsAndCoords.second;

    // Now using only local operations, find the right hand side nodes owned
    // by this processor and the matching ones on the left that were previously calculated
    /////////////////////////////////////////////////////////////////////////
    Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > locallyMatchedIds;

    bool failure = false;
    try {
       locallyMatchedIds = getLocallyMatchedSideIds(sideIds,sideCoords,mesh,right,matcher);
    } catch(std::logic_error & e) {
       locallyMatchedIds = Teuchos::rcp(new std::vector<std::pair<std::size_t,std::size_t> >);
       failure = true;
    } 

    // Get the ids on the left required by this processor (they maybe ghosted), 
    // and using the matched ids computed above over all processors, find the
    // corrsponding node on the right boundary.
    /////////////////////////////////////////////////////////////////////////

    // next line requires communication
    Teuchos::RCP<std::vector<std::size_t> > locallyRequiredIds = getLocalSideIds(mesh,left);
    Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds
          = getGlobalPairing(*locallyRequiredIds,*locallyMatchedIds,mesh,failure);

    // now you have a pair of ids that maps ids on the left required by this processor 
    // to ids on the right
    
    return globallyMatchedIds;
  }

  Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >
  getGlobalPairing(const std::vector<std::size_t> & locallyRequiredIds,
                   const std::vector<std::pair<std::size_t,std::size_t> > & locallyMatchedIds,
                   const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,bool failure)
  {
     Epetra_MpiComm Comm(mesh->getBulkData()->parallel());

     // this is needed to prevent hanging: it is unfortunately expensive
     // need a better way!
     int myVal = failure ? 1 : 0;
     int sumVal = 0;
     Comm.SumAll(&myVal,&sumVal,1);
     TEUCHOS_ASSERT(sumVal==0);

     std::vector<int> requiredInts(locallyRequiredIds.size());
     for(std::size_t i=0;i<requiredInts.size();i++) 
        requiredInts[i] = locallyRequiredIds[i];

     std::vector<int> providedInts(locallyMatchedIds.size());
     for(std::size_t i=0;i<locallyMatchedIds.size();i++) 
        providedInts[i] = locallyMatchedIds[i].first;

     // maps and communciation all set up
     Epetra_Map requiredMap(-1,requiredInts.size(),&requiredInts[0],0,Comm);
     Epetra_Map providedMap(-1,providedInts.size(),&providedInts[0],0,Comm);
     Epetra_Import importer(requiredMap,providedMap); 
     
     // this is what to distribute
     Epetra_IntVector providedVector(providedMap);
     for(std::size_t i=0;i<locallyMatchedIds.size();i++) 
        providedVector[i] = locallyMatchedIds[i].second;

     // vector to fill
     Epetra_IntVector requiredVector(requiredMap);
     TEUCHOS_ASSERT(requiredVector.Import(providedVector,importer,Insert)==0);
     int * myMappedIds = requiredVector.Values();

     Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > result
           = Teuchos::rcp(new std::vector<std::pair<std::size_t,std::size_t> >(requiredInts.size()));
     for(std::size_t i=0;i<result->size();i++) {
        (*result)[i].first = requiredInts[i];
        (*result)[i].second = myMappedIds[i];
     } 
    
     return result;
  }


  template <typename MatchObject>
  Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >
  getLocallyMatchedSideIds(const std::vector<std::size_t> & side_ids,
                           const std::vector<Teuchos::Tuple<double,3> > & side_coords,
                           const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
                           const std::string & sideName,const MatchObject & matcher)
  {
     using Teuchos::RCP;
     using Teuchos::Tuple;

     RCP<std::vector<std::pair<std::size_t,std::size_t> > > result
           = Teuchos::rcp(new std::vector<std::pair<std::size_t,std::size_t> >); 

     // grab local IDs and coordinates on this side
     //////////////////////////////////////////////////////////////////

     std::pair<Teuchos::RCP<std::vector<std::size_t> >,
               Teuchos::RCP<std::vector<Teuchos::Tuple<double,3> > > > sidePair =
            getLocalSideIdsAndCoords(mesh,sideName);

     std::vector<std::size_t> & local_side_ids = *sidePair.first;
     std::vector<Teuchos::Tuple<double,3> > & local_side_coords = *sidePair.second;

     bool checkProb = false;
     std::vector<bool> side_flags(side_ids.size(),false);

     // do a slow search for the coordinates: this _can_ be sped
     // up! (considered searches in sorted ranges using the sorted_permutation function)
     ////////////////////////////////////////////////////////
     for(std::size_t localNode=0;localNode<local_side_ids.size();localNode++) { 
        std::size_t local_gid = local_side_ids[localNode];
        const Tuple<double,3> & local_coord = local_side_coords[localNode];

        // loop over globally distributed coordinates and fine a match
        for(std::size_t globalNode=0;globalNode<side_ids.size();globalNode++) { 
           std::size_t global_gid = side_ids[globalNode];
           const Tuple<double,3> & global_coord = side_coords[globalNode];

           if(matcher(global_coord,local_coord)) {
              if(side_flags[globalNode]) // has this node been matched by this
                 checkProb = true;       // processor?

              result->push_back(std::make_pair(global_gid,local_gid));
              side_flags[globalNode] = true; 
              continue;
           }
        }
     }

     // make sure you matched everything you can: If this throws...it can 
     // cause the process to hang!
     TEUCHOS_ASSERT(not checkProb);
     TEUCHOS_ASSERT(local_side_ids.size()==result->size());

     return result;
  }


  /** This returns the locally resident (includes ghosted) global IDs
    * for a particular side. 
    */
  Teuchos::RCP<std::vector<std::size_t> >
  getLocalSideIds(const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
                  const std::string & sideName)
  {
     Teuchos::RCP<stk::mesh::MetaData> metaData = mesh->getMetaData();
     Teuchos::RCP<stk::mesh::BulkData> bulkData = mesh->getBulkData();

     // grab nodes owned by requested side
     /////////////////////////////////////////////
     stk::mesh::Part * side = metaData->get_part(sideName,"Some error message");
     stk::mesh::Selector mySides = *side;
 
     std::vector<stk::mesh::Bucket*> nodeBuckets;
     stk::mesh::get_buckets(mySides,bulkData->buckets(mesh->getNodeRank()),nodeBuckets);

     // build id vector
     ////////////////////////////////////////////
     std::size_t nodeCount = 0;
     for(std::size_t b=0;b<nodeBuckets.size();b++)
        nodeCount += nodeBuckets[b]->size();

     Teuchos::RCP<std::vector<std::size_t> > sideIds
        = Teuchos::rcp(new std::vector<std::size_t>(nodeCount));

     // loop over node buckets
     for(std::size_t b=0,index=0;b<nodeBuckets.size();b++) {
        stk::mesh::Bucket & bucket = *nodeBuckets[b]; 
           
        for(std::size_t n=0;n<bucket.size();n++,index++)
           (*sideIds)[index] = bucket[n].identifier();
     }

     return sideIds;
  }

  std::pair<Teuchos::RCP<std::vector<std::size_t> >,
            Teuchos::RCP<std::vector<Teuchos::Tuple<double,3> > > >
  getLocalSideIdsAndCoords(const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
                           const std::string & sideName)
  {
     unsigned physicalDim = mesh->getDimension();
     
     Teuchos::RCP<stk::mesh::MetaData> metaData = mesh->getMetaData();
     Teuchos::RCP<stk::mesh::BulkData> bulkData = mesh->getBulkData();

     // grab nodes owned by requested side
     /////////////////////////////////////////////
     stk::mesh::Part * side = metaData->get_part(sideName,"Some error message");
     stk::mesh::Selector mySides = (*side) & metaData->locally_owned_part();
 
     std::vector<stk::mesh::Bucket*> nodeBuckets;
     stk::mesh::get_buckets(mySides,bulkData->buckets(mesh->getNodeRank()),nodeBuckets);

     // build id vector
     ////////////////////////////////////////////
     std::size_t nodeCount = 0;
     for(std::size_t b=0;b<nodeBuckets.size();b++)
        nodeCount += nodeBuckets[b]->size();

     Teuchos::RCP<std::vector<std::size_t> > sideIds
        = Teuchos::rcp(new std::vector<std::size_t>(nodeCount));
     Teuchos::RCP<std::vector<Teuchos::Tuple<double,3> > > sideCoords
        = Teuchos::rcp(new std::vector<Teuchos::Tuple<double,3> >(nodeCount));

     // loop over node buckets
     for(std::size_t b=0,index=0;b<nodeBuckets.size();b++) {
        stk::mesh::Bucket & bucket = *nodeBuckets[b]; 
        stk::mesh::BucketArray<panzer_stk::STK_Interface::VectorFieldType> array(mesh->getCoordinatesField(),bucket);
           
        for(std::size_t n=0;n<bucket.size();n++,index++) {
           (*sideIds)[index] = bucket[n].identifier();
           Teuchos::Tuple<double,3> & coord = (*sideCoords)[index];
           
           // copy coordinates into multi vector
           for(std::size_t d=0;d<physicalDim;d++)
              coord[d] = array(d,n);
        }
     }

     return std::make_pair(sideIds,sideCoords);
  }

  std::pair<Teuchos::RCP<std::vector<std::size_t> >,
            Teuchos::RCP<std::vector<Teuchos::Tuple<double,3> > > >
  getSideIdsAndCoords(const Teuchos::RCP<panzer_stk::STK_Interface> & mesh,
                const std::string & sideName)
  {
     Epetra_MpiComm Comm(mesh->getBulkData()->parallel());

     unsigned physicalDim = mesh->getDimension();
 
     // grab local IDs and coordinates on this side
     // and build local epetra vector
     //////////////////////////////////////////////////////////////////

     std::pair<Teuchos::RCP<std::vector<std::size_t> >,
               Teuchos::RCP<std::vector<Teuchos::Tuple<double,3> > > > sidePair =
            getLocalSideIdsAndCoords(mesh,sideName);

     std::vector<std::size_t> & local_side_ids = *sidePair.first;
     std::vector<Teuchos::Tuple<double,3> > & local_side_coords = *sidePair.second;
     int nodeCount = local_side_ids.size();

     // build local Epetra objects
     Epetra_Map idMap(-1,nodeCount,0,Comm);
     Teuchos::RCP<Epetra_IntVector> localIdVec = Teuchos::rcp(new Epetra_IntVector(idMap));
     Teuchos::RCP<Epetra_MultiVector> localCoordVec = Teuchos::rcp(new Epetra_MultiVector(idMap,physicalDim));

     // copy local Ids into Epetra vector
     for(std::size_t n=0;n<local_side_ids.size();n++) {
        std::size_t nodeId = local_side_ids[n];
        Teuchos::Tuple<double,3> & coords = local_side_coords[n];

        (*localIdVec)[n] = nodeId;
        for(unsigned d=0;d<physicalDim;d++)
           (*(*localCoordVec)(d))[n] = coords[d];
     }

     // fully distribute epetra vector across all processors 
     // (these are "distributed" or "dist" objects)
     //////////////////////////////////////////////////////////////

     int dist_nodeCount = idMap.NumGlobalElements();

     // build global epetra objects
     Epetra_LocalMap distMap(dist_nodeCount,0,Comm);
     Teuchos::RCP<Epetra_IntVector> distIdVec = Teuchos::rcp(new Epetra_IntVector(distMap));
     Teuchos::RCP<Epetra_MultiVector> distCoordVec = Teuchos::rcp(new Epetra_MultiVector(distMap,physicalDim));

     // export to the localVec object from the "vector" object
     Epetra_Import importer(distMap,idMap);
     TEUCHOS_ASSERT(distIdVec->Import(*localIdVec,importer,Insert)==0);
     TEUCHOS_ASSERT(distCoordVec->Import(*localCoordVec,importer,Insert)==0);

     // convert back to generic stl vector objects
     ///////////////////////////////////////////////////////////

     Teuchos::RCP<std::vector<std::size_t> > dist_side_ids
        = Teuchos::rcp(new std::vector<std::size_t>(dist_nodeCount));
     Teuchos::RCP<std::vector<Teuchos::Tuple<double,3> > > dist_side_coords
        = Teuchos::rcp(new std::vector<Teuchos::Tuple<double,3> >(dist_nodeCount));

     // copy local Ids into Epetra vector
     for(std::size_t n=0;n<dist_side_ids->size();n++) {
        (*dist_side_ids)[n] = (*distIdVec)[n];

        Teuchos::Tuple<double,3> & coords = (*dist_side_coords)[n];
        for(unsigned d=0;d<physicalDim;d++)
           coords[d] = (*(*distCoordVec)(d))[n];
     }

     return std::make_pair(dist_side_ids,dist_side_coords);
  }

}
