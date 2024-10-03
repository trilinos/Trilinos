// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_STK_PeriodicBC_Matcher.hpp"

#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/FieldBase.hpp>

#include "Panzer_NodeType.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Vector.hpp"

#include "Teuchos_FancyOStream.hpp"

namespace panzer_stk {

namespace periodic_helpers {

Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >
getGlobalPairing(const std::vector<std::size_t> & locallyRequiredIds,
                 const std::vector<std::pair<std::size_t,std::size_t> > & locallyMatchedIds,
                 const STK_Interface & mesh,bool failure)
{
   using LO = panzer::LocalOrdinal;
   using GO = panzer::GlobalOrdinal;
   using NODE = panzer::TpetraNodeType;
   using Map = Tpetra::Map<LO,GO,NODE>;
   using Importer = Tpetra::Import<LO,GO,NODE>;

   auto comm = mesh.getComm();

   // this is needed to prevent hanging: it is unfortunately expensive
   // need a better way!
   int myVal = failure ? 1 : 0;
   int sumVal = 0;
   Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_SUM,1,&myVal,&sumVal);
   TEUCHOS_ASSERT(sumVal==0);

   std::vector<GO> requiredInts(locallyRequiredIds.size());
   for(std::size_t i=0;i<requiredInts.size();i++) 
      requiredInts[i] = locallyRequiredIds[i];

   std::vector<GO> providedInts(locallyMatchedIds.size());
   for(std::size_t i=0;i<locallyMatchedIds.size();i++) 
      providedInts[i] = locallyMatchedIds[i].first;

   // maps and communciation all set up
   auto computeInternally = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
   Teuchos::RCP<Map> requiredMap = Teuchos::rcp(new Map(computeInternally,Teuchos::ArrayView<const GO>(requiredInts),0,comm));
   Teuchos::RCP<Map> providedMap = Teuchos::rcp(new Map(computeInternally,Teuchos::ArrayView<const GO>(providedInts),0,comm));
   Importer importer(providedMap,requiredMap);

   // this is what to distribute
   Tpetra::Vector<GO,LO,GO,NODE> providedVector(providedMap);
   {
     auto pvHost = providedVector.getLocalViewHost(Tpetra::Access::OverwriteAll);
     for(std::size_t i=0;i<locallyMatchedIds.size();i++)
       pvHost(i,0) = locallyMatchedIds[i].second;
   }

   Tpetra::Vector<GO,LO,GO,NODE> requiredVector(requiredMap);
   requiredVector.doImport(providedVector,importer,Tpetra::INSERT);
   

   Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > result
         = Teuchos::rcp(new std::vector<std::pair<std::size_t,std::size_t> >(requiredInts.size()));

   auto rvHost = requiredVector.getLocalViewHost(Tpetra::Access::ReadOnly);
   for(std::size_t i=0;i<result->size();i++) {
      (*result)[i].first = requiredInts[i];
      (*result)[i].second = rvHost(i,0);
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
   panzer::GlobalOrdinal offset = 0; // offset to avoid giving nodes, edges, faces the same sideId
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
   stk::mesh::EntityId offset = 0;
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

         // need to ensure that higher dimensions are properly zeroed
         // required for 1D periodic boundary conditions
         for(std::size_t d=physicalDim;d<3;d++)
           coord[d] = 0;
      }
   }

   return std::make_pair(sideIds,sideCoords);
}

std::pair<Teuchos::RCP<std::vector<std::size_t> >,
          Teuchos::RCP<std::vector<Teuchos::Tuple<double,3> > > >
getSideIdsAndCoords(const STK_Interface & mesh,
              const std::string & sideName, const std::string type_)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using LO = panzer::LocalOrdinal;
   using GO = panzer::GlobalOrdinal;
   using NODE = panzer::TpetraNodeType;
   using Map = Tpetra::Map<LO,GO,NODE>;
   using Importer = Tpetra::Import<LO,GO,NODE>;

   // Epetra_MpiComm Comm(mesh.getBulkData()->parallel());
   auto comm = mesh.getComm();

   unsigned physicalDim = mesh.getDimension();
 
   // grab local IDs and coordinates on this side
   // and build local epetra vector
   //////////////////////////////////////////////////////////////////

   std::pair<Teuchos::RCP<std::vector<std::size_t> >,
             Teuchos::RCP<std::vector<Teuchos::Tuple<double,3> > > > sidePair =
          getLocalSideIdsAndCoords(mesh,sideName,type_);

   std::vector<std::size_t> & local_side_ids = *sidePair.first;
   std::vector<Teuchos::Tuple<double,3> > & local_side_coords = *sidePair.second;
   std::size_t nodeCount = local_side_ids.size();

   // build local Tpetra objects
   auto computeInternally = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
   RCP<Map> idMap_ = rcp(new Map(computeInternally,nodeCount,0,comm));
   RCP<Tpetra::Vector<GO,LO,GO,NODE>> localIdVec_ = rcp(new Tpetra::Vector<GO,LO,GO,NODE>(idMap_));
   RCP<Tpetra::MultiVector<double,LO,GO,NODE>> localCoordVec_ = rcp(new Tpetra::MultiVector<double,LO,GO,NODE>(idMap_,physicalDim));

   // copy local Ids and coords into Tpetra vectors
   {
     auto lidHost = localIdVec_->getLocalViewHost(Tpetra::Access::OverwriteAll);
     auto lcoordHost = localCoordVec_->getLocalViewHost(Tpetra::Access::OverwriteAll);
     for(std::size_t n=0;n<local_side_ids.size();n++) {
       std::size_t nodeId = local_side_ids[n];
       Teuchos::Tuple<double,3> & coords = local_side_coords[n];

       lidHost(n,0) = static_cast<GO>(nodeId);
       for(unsigned d=0;d<physicalDim;d++)
         lcoordHost(n,d) = coords[d];
     }
   }

   // fully distribute epetra vector across all processors 
   // (these are "distributed" or "dist" objects)
   //////////////////////////////////////////////////////////////

   std::size_t dist_nodeCount = idMap_->getGlobalNumElements();

   // build global Tpetra objects
   RCP<Map> distMap_ = rcp(new Map(dist_nodeCount,0,comm,Tpetra::LocallyReplicated));
   RCP<Tpetra::Vector<GO,LO,GO,NODE>> distIdVec_ = rcp(new Tpetra::Vector<GO,LO,GO,NODE>(distMap_));
   RCP<Tpetra::MultiVector<double,LO,GO,NODE>> distCoordVec_ = rcp(new Tpetra::MultiVector<double,LO,GO,NODE>(distMap_,physicalDim));

   // export to the localVec object from the "vector" object
   Importer importer_(idMap_,distMap_);
   distIdVec_->doImport(*localIdVec_,importer_,Tpetra::INSERT);
   distCoordVec_->doImport(*localCoordVec_,importer_,Tpetra::INSERT);

   // convert back to generic stl vector objects
   ///////////////////////////////////////////////////////////

   Teuchos::RCP<std::vector<std::size_t> > dist_side_ids
      = Teuchos::rcp(new std::vector<std::size_t>(dist_nodeCount));
   Teuchos::RCP<std::vector<Teuchos::Tuple<double,3> > > dist_side_coords
      = Teuchos::rcp(new std::vector<Teuchos::Tuple<double,3> >(dist_nodeCount));

   // copy local Ids from Tpetra vector
   const auto didHost = distIdVec_->getLocalViewHost(Tpetra::Access::ReadOnly);
   const auto dcoordHost = distCoordVec_->getLocalViewHost(Tpetra::Access::ReadOnly);
   for(std::size_t n=0;n<dist_side_ids->size();++n) {
     (*dist_side_ids)[n] = didHost(n,0);

      Teuchos::Tuple<double,3> & coords = (*dist_side_coords)[n];
      for(unsigned d=0;d<physicalDim;++d) {
        coords[d] = dcoordHost(n,d);
      }
      // ensure that higher dimensions are zero
      for(unsigned d=physicalDim;d<3;++d)
         coords[d] = 0;
   }

   return std::make_pair(dist_side_ids,dist_side_coords);
}

} // end periodic_helpers

} // end panzer_stk
