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

#ifdef PANZER_HAVE_STKSEARCH
namespace panzer_stk {

namespace periodic_helpers {

void fillLocalSearchVector(const STK_Interface & mesh, SphereIdVector & searchVector, const double & error,
                           const std::string & sideName, const std::string & type_, const bool & getGhostedIDs) 
{

   // empty optional arguments for real call below
   // this is partially for backwards compatability but also recognizes that for the first
   // pairing, these vectors are not needed
   // they are also not used for side B in every case

   std::vector<std::string> matchedSides;
   std::vector<SearchId> potentialIDsToRemap;

   fillLocalSearchVector(mesh,searchVector,error,sideName,type_,getGhostedIDs,matchedSides,potentialIDsToRemap);

   return;

}

void fillLocalSearchVector(const STK_Interface & mesh, SphereIdVector & searchVector, const double & error,
                           const std::string & sideName, const std::string & type_, const bool & getGhostedIDs,
                           const std::vector<std::string> & matchedSides, std::vector<SearchId> & potentialIDsToRemap)
{
   unsigned physicalDim = mesh.getDimension();
   
   Teuchos::RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();
   Teuchos::RCP<stk::mesh::BulkData> bulkData = mesh.getBulkData();

   const unsigned parallelRank = bulkData->parallel_rank();

   // grab entities owned by requested side
   /////////////////////////////////////////////
   std::stringstream ss;
   ss << "Can't find a sideset named \"" << sideName << "\" in the mesh" << std::endl;
   stk::mesh::Part * side = metaData->get_part(sideName,ss.str().c_str());

   // if ghosted IDs are requested, add in the shared portion
   stk::mesh::Selector mySides = getGhostedIDs ? 
      (*side) & (metaData->locally_owned_part() | metaData->globally_shared_part()) : 
      (*side) & metaData->locally_owned_part();

   stk::mesh::EntityRank rank;
   const STK_Interface::VectorFieldType * field = 0;
   // the sorting further downstream only uses ids which are offset
   // based on the entity type
   stk::mesh::EntityId offset = 0;
   if(type_ == "coord"){
     rank = mesh.getNodeRank();
     field = & mesh.getCoordinatesField();
     // no offset
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

   // remove previously matched nodes

   stk::mesh::Selector intersection;  // empty set to start

   if (matchedSides.size()>0) {
      for (size_t j=0; j<matchedSides.size(); ++j) {
         auto previouslyMatched = matchedSides[j];
         // add in the overlap between the requested side and the previously matched side
         intersection = intersection | (mySides & *(metaData->get_part(previouslyMatched)));
      }
      // remove all overlaps
      mySides = mySides - intersection;
   }

   // get buckets
   std::vector<stk::mesh::Bucket*> const& entityBuckets =
     bulkData->get_buckets(rank, mySides);

   // loop over entity buckets

   for(size_t b=0;b<entityBuckets.size();b++) {
      stk::mesh::Bucket & bucket = *entityBuckets[b];
      double const* array = stk::mesh::field_data(*field, bucket);

      // loop over entities
      for(size_t n=0;n<bucket.size();n++) {

         double coord[3]; // coordinates
         // copy coordinates into multi vector
         for(size_t d=0;d<physicalDim;d++)
            coord[d] = array[physicalDim*n + d];

         // need to ensure that higher dimensions are properly zeroed
         // required for 1D periodic boundary conditions
         for(size_t d=physicalDim;d<3;d++)
           coord[d] = 0;

         // add to the coordinate and id to the search vector
         // a tolerance can be specified
         // TODO allow for relative tolerances...
         stk::search::Point<double> center(coord[0],coord[1],coord[2]);
         stk::mesh::EntityKey trueKey = bulkData->entity_key(bucket[n]);
         stk::mesh::EntityKey shiftedKey(trueKey.rank(), trueKey.id()+offset);
         SearchId search_id(shiftedKey, parallelRank);

         searchVector.emplace_back( Sphere(center, error), search_id);
      }
   }

   // for multiperiodic case, populate the potentialIDsToRemap vector with the IDs that have
   // already been matched and fall on this side

   if (matchedSides.size()>0) {
      TEUCHOS_ASSERT(potentialIDsToRemap.size()==0);
      // reset mySides
      mySides = getGhostedIDs ? 
         (*side) & (metaData->locally_owned_part() | metaData->globally_shared_part()) : 
         (*side) & metaData->locally_owned_part();

      std::vector<stk::mesh::Bucket*> const & intersectionEntityBuckets = 
        bulkData->get_buckets(rank, mySides & intersection);

      // loop over entity buckets

      for(size_t b=0;b<intersectionEntityBuckets.size();b++) {
         stk::mesh::Bucket & bucket = *intersectionEntityBuckets[b];
         // loop over entities
         for(size_t n=0;n<bucket.size();n++) {
            stk::mesh::EntityKey trueKey = bulkData->entity_key(bucket[n]);
            stk::mesh::EntityKey shiftedKey(trueKey.rank(), trueKey.id()+offset);
            potentialIDsToRemap.emplace_back(shiftedKey, parallelRank);
         }
      }
   }

   return;
}

const std::vector<double> computeGlobalCentroid(const STK_Interface & mesh, const std::string & sideName)
{
   // TODO too much replicated code here 
   unsigned physicalDim = mesh.getDimension();
   
   Teuchos::RCP<stk::mesh::MetaData> metaData = mesh.getMetaData();
   Teuchos::RCP<stk::mesh::BulkData> bulkData = mesh.getBulkData();

   // grab requested side
   /////////////////////////////////////////////
   std::stringstream ss;
   ss << "Can't find part=\"" << sideName << "\"" << std::endl;
   stk::mesh::Part * side = metaData->get_part(sideName,ss.str().c_str());
   stk::mesh::Selector mySides = (*side) & metaData->locally_owned_part();

   // get node buckets
   std::vector<stk::mesh::Bucket*> const& entityBuckets =
     bulkData->get_buckets(mesh.getNodeRank(), mySides);

   // loop over node buckets
   double localCentroid[3] = {0.,0.,0.};
   int localNodeCount = 0;
   for(size_t b=0;b<entityBuckets.size();b++) {
      stk::mesh::Bucket & bucket = *entityBuckets[b];
      double const* array = stk::mesh::field_data(mesh.getCoordinatesField(), bucket);

      // loop over nodes 
      for(size_t n=0;n<bucket.size();n++) {

         ++localNodeCount;
         // sum (note that unused dimensions are skipped)
         for(size_t d=0;d<physicalDim;d++)
            localCentroid[d] += array[physicalDim*n + d];

      }
   }
   int globalNodeCount = 0;

   auto comm = mesh.getComm();

   double globalCentroid[3] = { };

   Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM,1,&localNodeCount,&globalNodeCount);
   Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM,3,&localCentroid[0],&globalCentroid[0]);

   std::vector<double> result = {0.,0.,0.};
   for (size_t d=0;d<physicalDim;d++)
      result[d] = globalCentroid[d]/globalNodeCount;

   return result;
}

void updateMapping(Teuchos::RCP<std::vector<std::pair<size_t,size_t> > > & currentMatches,
                   const std::vector<std::pair<size_t,size_t> > & previousMatches,
                   const std::vector<SearchId> & IDsToRemap, const STK_Interface & mesh)
{

   using LO = panzer::LocalOrdinal;
   using GO = panzer::GlobalOrdinal;
   using NODE = panzer::TpetraNodeType;
   using Map = Tpetra::Map<LO,GO,NODE>;
   using Importer = Tpetra::Import<LO,GO,NODE>;

   auto comm = mesh.getComm();

   // store maps
   // this is necessary because of the uniqueness requirements
   // and convenient to update the map

   std::map<size_t,size_t> myPreviousAtoB,myCurrentAtoB;
   std::map<size_t,std::vector<size_t> > myPreviousBtoA;
   for (size_t i=0;i<previousMatches.size();++i) {
      myPreviousAtoB[previousMatches[i].first] = previousMatches[i].second;
      // may not be one-to-one
      myPreviousBtoA[previousMatches[i].second].push_back(previousMatches[i].first);
   }
   for (size_t i=0;i<currentMatches->size();++i)
      myCurrentAtoB[(*currentMatches)[i].first] = (*currentMatches)[i].second;

   // find which IDs we need to query to get THEIR A to B map
   // this means we need to the the B id of our previous match for the IDsToRemap
   std::vector<GO> requestedAIDs;

   for (auto & id : IDsToRemap)
      requestedAIDs.push_back(myPreviousAtoB[id.id().id()]);

   // quick and dirty way to get rid of repeated entries on each process
   std::set<GO> uniqueAIDs(requestedAIDs.begin(),requestedAIDs.end());
   requestedAIDs = std::vector<GO>(uniqueAIDs.begin(),uniqueAIDs.end());

   // find which A ids we have in our new mapping
   std::vector<GO> newAIDs(currentMatches->size());

   for (size_t i=0;i<currentMatches->size();++i) {
      newAIDs[i] = (*currentMatches)[i].first;
   }

   // create teuchos maps
   auto computeInternally = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
   Teuchos::RCP<const Map> testMap = Teuchos::rcp(new Map(computeInternally,&newAIDs[0],newAIDs.size(),0,comm));
   Teuchos::RCP<const Map> newAIDsMap; 
   // source must be unique across communicator
   if (!testMap->isOneToOne()){
      newAIDsMap = Tpetra::createOneToOne<LO,GO,NODE>(testMap);
   } else {
      newAIDsMap = testMap;
   }

   Teuchos::RCP<Map> requestedAIDsMap = Teuchos::rcp(new Map(computeInternally,&requestedAIDs[0],requestedAIDs.size(),0,comm));

   Importer importer(newAIDsMap,requestedAIDsMap);

   // send out my A to B map
   Tpetra::Vector<GO,LO,GO,NODE> newBIDs(newAIDsMap);
   auto newBIDsHost = newBIDs.getLocalViewHost(Tpetra::Access::OverwriteAll);
   auto myGIDs = newAIDsMap->getMyGlobalIndices();
   for (size_t i=0;i<myGIDs.size();++i)
      newBIDsHost(i,0) = myCurrentAtoB[myGIDs[i]];

   Tpetra::Vector<GO,LO,GO,NODE> requestedBIDs(requestedAIDsMap);
   requestedBIDs.doImport(newBIDs,importer,Tpetra::INSERT);
   auto requestedBIDsHost = requestedBIDs.getLocalViewHost(Tpetra::Access::ReadOnly);

   // now overwrite previous map where necessary...
   // what about error checking? what is the default is something is requested but not there?
   // TODO this assumes that teuchos maps and communication does not
   // alter the ordering in anyway so that AIDs and IDsToRemap correspond appropriately
   size_t ind = 0;
   for (const auto & id : requestedAIDs) {
      // get the corresponding ids to update in the previous map
      for (const auto & idToUpdate : myPreviousBtoA[id])
         // update with the final B id
         myPreviousAtoB[idToUpdate] = requestedBIDsHost(ind,0);
      ++ind;
   }
   
   // and add to new map...
   // needs to respect the previous ordering or type_vec in getPeriodicNodePairing will be wrong
   for (const auto & AB : previousMatches) {
      // so we get the A ids in order
      auto id = AB.first;
      // and use the updated previous A to B map
      (*currentMatches).emplace_back(std::pair<size_t,size_t>(id,myPreviousAtoB[id]));
   }

   return;
}

void appendMapping(Teuchos::RCP<std::vector<std::pair<size_t,size_t> > > & currentMatches,
                   const std::vector<std::pair<size_t,size_t> > & previousMatches)
{
   // add previous mapping to new map
   for (const auto & AB : previousMatches)
      (*currentMatches).push_back(AB);

   return;
}

} // end periodic_helpers

} // end panzer_stk
#endif
