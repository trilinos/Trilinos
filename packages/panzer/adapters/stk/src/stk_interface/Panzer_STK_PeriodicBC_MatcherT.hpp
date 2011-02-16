#include "Teuchos_Tuple.hpp"
#include "Teuchos_RCP.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_Utilities.hpp"

namespace panzer_stk {
namespace periodic_helpers {

template <typename Matcher>
Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >
matchPeriodicSides(const std::string & left,const std::string & right,
                     const STK_Interface & mesh,
                     const Matcher & matcher)
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
            RCP<std::vector<Tuple<double,3> > > > idsAndCoords = panzer_stk::periodic_helpers::getSideIdsAndCoords(mesh,left);
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

template <typename Matcher>
Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >
getLocallyMatchedSideIds(const std::vector<std::size_t> & side_ids,
                         const std::vector<Teuchos::Tuple<double,3> > & side_coords,
                         const panzer_stk::STK_Interface & mesh,
                         const std::string & sideName,const Matcher & matcher)
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

} // end periodic_helpers
} // end panzer_stk
