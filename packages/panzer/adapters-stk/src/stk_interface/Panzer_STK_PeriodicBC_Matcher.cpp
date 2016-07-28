// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include "Panzer_STK_PeriodicBC_Matcher.hpp"

#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/FieldBase.hpp>

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Import.h"

#include "Teuchos_FancyOStream.hpp"

namespace panzer_stk {

namespace periodic_helpers {

Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >
getGlobalPairing(const std::vector<std::size_t> & locallyRequiredIds,
                 const std::vector<std::pair<std::size_t,std::size_t> > & locallyMatchedIds,
                 const STK_Interface & mesh,bool failure)
{
   Epetra_MpiComm Comm(mesh.getBulkData()->parallel());

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
   int* requiredIntsPtr = NULL;
   if (requiredInts.size() > 0)
     requiredIntsPtr = &requiredInts[0];
   int* providedIntsPtr = NULL;
   if (providedInts.size() > 0)
     providedIntsPtr = &providedInts[0];
   Epetra_Map requiredMap(-1,requiredInts.size(),requiredIntsPtr,0,Comm);
   Epetra_Map providedMap(-1,providedInts.size(),providedIntsPtr,0,Comm);
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



/** This returns the locally resident (includes ghosted) global IDs
  * for a particular side. 
  */
Teuchos::RCP<std::vector<std::size_t> >
getLocalSideIds(const STK_Interface & mesh,
                const std::string & sideName, const std::string type_)
{
   Teuchos::RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();
   Teuchos::RCP<stk::mesh::BulkData> bulkData = mesh.getBulkData();

   // grab nodes owned by requested side
   /////////////////////////////////////////////
   std::stringstream ss;
   ss << "Can't find part=\"" << sideName << "\"" << std::endl;
   stk::mesh::Part * side = metaData->get_part(sideName,ss.str().c_str());
   stk::mesh::Selector mySides = *side & (metaData->locally_owned_part() | metaData->globally_shared_part());

   stk::mesh::EntityRank rank;
   unsigned int offset = 0; // offset to avoid giving nodes, edges, faces the same sideId 
   if(type_ == "coord"){
     rank = mesh.getNodeRank();
   } else if(type_ == "edge"){
     rank = mesh.getEdgeRank();
     offset = mesh.getMaxEntityId(mesh.getNodeRank());
   } else if(type_ == "face"){
     rank = mesh.getFaceRank();
     offset = mesh.getMaxEntityId(mesh.getNodeRank())+mesh.getMaxEntityId(mesh.getEdgeRank());
   } else {
     ss << "Can't do BCs of type " << type_  << std::endl;
     TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error, ss.str())
   }

   std::vector<stk::mesh::Bucket*> const& nodeBuckets =
     bulkData->get_buckets(rank, mySides);

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
        (*sideIds)[index] = bulkData->identifier(bucket[n]) + offset;
   }

   return sideIds;
}

std::pair<Teuchos::RCP<std::vector<std::size_t> >,
          Teuchos::RCP<std::vector<Teuchos::Tuple<double,3> > > >
getLocalSideIdsAndCoords(const STK_Interface & mesh,
                         const std::string & sideName, const std::string type_)
{
   unsigned physicalDim = mesh.getDimension();
   
   Teuchos::RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();
   Teuchos::RCP<stk::mesh::BulkData> bulkData = mesh.getBulkData();

   // grab nodes owned by requested side
   /////////////////////////////////////////////
   std::stringstream ss;
   ss << "Can't find part=\"" << sideName << "\"" << std::endl;
   stk::mesh::Part * side = metaData->get_part(sideName,ss.str().c_str());
   stk::mesh::Selector mySides = (*side) & metaData->locally_owned_part();

   stk::mesh::EntityRank rank;
   const STK_Interface::VectorFieldType * field = 0;
   unsigned int offset = 0;
   if(type_ == "coord"){
     rank = mesh.getNodeRank();
     field = & mesh.getCoordinatesField();
   } else if(type_ == "edge"){
     rank = mesh.getEdgeRank();
     field = & mesh.getEdgesField();
     offset = mesh.getMaxEntityId(mesh.getNodeRank());
   } else if(type_ == "face"){
     rank = mesh.getFaceRank();
     field = & mesh.getFacesField();
     offset = mesh.getMaxEntityId(mesh.getNodeRank())+mesh.getMaxEntityId(mesh.getEdgeRank());
   } else {
     ss << "Can't do BCs of type " << type_  << std::endl;
     TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error, ss.str())
   }

   std::vector<stk::mesh::Bucket*> const& nodeBuckets =
     bulkData->get_buckets(rank, mySides);

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
      double const* array = stk::mesh::field_data(*field, bucket);

      for(std::size_t n=0;n<bucket.size();n++,index++) {
         (*sideIds)[index] = bulkData->identifier(bucket[n]) + offset;
         Teuchos::Tuple<double,3> & coord = (*sideCoords)[index];

         // copy coordinates into multi vector
         for(std::size_t d=0;d<physicalDim;d++)
            coord[d] = array[physicalDim*n + d];
      }
   }

   return std::make_pair(sideIds,sideCoords);
}

std::pair<Teuchos::RCP<std::vector<std::size_t> >,
          Teuchos::RCP<std::vector<Teuchos::Tuple<double,3> > > >
getSideIdsAndCoords(const STK_Interface & mesh,
              const std::string & sideName, const std::string type_)
{
   Epetra_MpiComm Comm(mesh.getBulkData()->parallel());

   unsigned physicalDim = mesh.getDimension();
 
   // grab local IDs and coordinates on this side
   // and build local epetra vector
   //////////////////////////////////////////////////////////////////

   std::pair<Teuchos::RCP<std::vector<std::size_t> >,
             Teuchos::RCP<std::vector<Teuchos::Tuple<double,3> > > > sidePair =
          getLocalSideIdsAndCoords(mesh,sideName,type_);

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

} // end periodic_helpers

} // end panzer_stk
