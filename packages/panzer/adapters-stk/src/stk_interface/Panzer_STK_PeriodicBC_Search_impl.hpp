// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Tuple.hpp"
#include "Teuchos_RCP.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"

#include "Teuchos_FancyOStream.hpp"

#include <set>

namespace panzer_stk {
namespace periodic_helpers {

template <typename Matcher>
Teuchos::RCP<std::vector<std::pair<size_t,size_t> > >
matchPeriodicSidesSearch(const std::string & sideA,const std::string & sideB,
                         const STK_Interface & mesh,
                         const Matcher & matcher, const std::string type_)
{

  // this function should be called the first time a match is made 
  // we create a few empty objects for the previous mapping

  std::vector<std::string> matchedSides;
  std::vector<std::pair<size_t,size_t> > previousMatches;

  // then pass along

  return matchPeriodicSidesSearch(sideA,sideB,mesh,matcher,matchedSides,previousMatches,type_);

}


template <typename Matcher>
Teuchos::RCP<std::vector<std::pair<size_t,size_t> > >
matchPeriodicSidesSearch(const std::string & sideA,const std::string & sideB,
                         const STK_Interface & mesh,
                         const Matcher & matcher, const std::vector<std::string> & matchedSides, 
                         const std::vector<std::pair<size_t,size_t> > & previousMatches,
                         const std::string type_)
{
  using Teuchos::Tuple;
  using Teuchos::RCP;
  using Teuchos::rcp;

  auto myRank = mesh.getBulkData()->parallel_rank();

  SphereIdVector coordsIdsA, coordsIdsB;
  std::vector<SearchId> IDsToRemap;

  // populate the search vectors with the node coordinates and ids
  // we will always search for ghosted IDs and repeats only on side A

  auto error = matcher.getAbsoluteTolerance();

  fillLocalSearchVector(mesh,coordsIdsA,error,sideA,type_,true,matchedSides,IDsToRemap);
  fillLocalSearchVector(mesh,coordsIdsB,error,sideB,type_,false);

  // apply the matcher transform to side B to effectively align the periodic entities
  // to do so we need the centroid of the other side

  // requires communication
  std::vector<double> centroidA = computeGlobalCentroid(mesh,sideA);

  // now transform
  transformLocalSearchVector(coordsIdsB,matcher,centroidA);

  // now we find the matches
  SearchPairVector results;
  stk::search::coarse_search(coordsIdsA,coordsIdsB,stk::search::KDTREE,mesh.getBulkData()->parallel(),results);

  // the results are pairs of matched A and B entity keys
  // if my process has a match, it will be stored
  // hence, we only keep the results if the A key proc matches our rank
  // so each process has a map myAIDs --> BIDs

  // we store this A to B map, adding it to the pre-existing
  // map of local periodic nodes to their matches, if necessary
  // note the ids have been adjusted for entity type already
  Teuchos::RCP<std::vector<std::pair<size_t,size_t> > > myMap
    = Teuchos::rcp(new std::vector<std::pair<size_t,size_t>>());

  for (size_t i=0; i<results.size(); ++i) {
     if (results[i].first.proc() == myRank) {
        // first id grabs the entity key which has another id and the entity rank
        (*myMap).emplace_back( 
           std::pair<size_t,size_t>(results[i].first.id().id(),results[i].second.id().id()) );
     }
  }

  TEUCHOS_TEST_FOR_EXCEPTION((*myMap).size()!=coordsIdsA.size(),std::logic_error,
                             "matchPeriodicSidesSearch: error in local match. "
                             "Number of matched IDs not equal to number of requested matches!");

  if (matchedSides.size()>0) {
    // guaranteed to have previous matches and they are of the same entity type
    // in this case we need to handle multiperiodicity
    updateMapping(myMap,previousMatches,IDsToRemap,mesh);
  } else if (previousMatches.size()>0) {
    // we have previous matches, but they are of a different entity type
    // in this case we just append the previous matches unaltered
    appendMapping(myMap,previousMatches);
  }

  return myMap;

}

template<typename Matcher> void
transformLocalSearchVector( SphereIdVector & searchVectorSideA, const Matcher & matcher, const std::vector<double> & centroidSideB)
{

   // loop over sphereIds objects and shift center according to the matcher's periodic transform

   for (auto && sphereIdSideA : searchVectorSideA ) 
      matcher.transform(&sphereIdSideA.first.center()[0],centroidSideB); 

   return;
}

} // end periodic_helpers
} // end panzer_stk
