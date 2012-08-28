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

#include "Teuchos_Tuple.hpp"
#include "Teuchos_RCP.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_Utilities.hpp"

#include "Teuchos_FancyOStream.hpp"

#include <set>

namespace panzer_stk {
namespace periodic_helpers {

template <typename Matcher>
Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >
matchPeriodicSides(const std::string & left,const std::string & right,
                     const STK_Interface & mesh,
                     const Matcher & matcher,
                     const std::vector<std::pair<std::size_t,std::size_t> >  & ownedToMapped)
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
  // Note: The matching check could definitely be sped up with a sorting operation
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

     Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
     out.setShowProcRank(true);
     out.setOutputToRootOnly(-1);

     out << "Not all sides matched expect failure: \n" << e.what() << std::endl;
  } 

  // Get the ids on the left required by this processor (they maybe ghosted), 
  // and using the matched ids computed above over all processors, find the
  // corrsponding node on the right boundary.
  /////////////////////////////////////////////////////////////////////////

  // THIS COMPLEX ALGORITHM COMPENSATES FOR THE PREVIOUS PAIRING BY UPDATING
  // THE REQUEST LIST. THERE ALSO IS A REQUIREMENT OF PER PROC UNIQUENESS IN THE EPETRA
  // IMPORT.  SO THERE IS SOME ADDITIONAL WORK DONE TO REMOVE REPEATED ENTRIES

  // build reverse map 
  std::map<std::size_t,std::vector<std::size_t> > reverseMap;
  for(std::size_t i=0;i<ownedToMapped.size();i++)
     reverseMap[ownedToMapped[i].second].push_back(ownedToMapped[i].first);
  
  Teuchos::RCP<std::vector<std::size_t> > locallyRequiredIds = getLocalSideIds(mesh,left);
  std::vector<std::size_t> saved_locallyRequiredIds = *locallyRequiredIds; // will be required
                                                                           // to check owner/ghostship
                                                                           // of IDs
  std::vector<std::pair<std::size_t,std::size_t> > unusedOwnedToMapped;

  // apply previous mappings to this set of local IDs
  for(std::size_t i=0;i<ownedToMapped.size();i++) {
     std::size_t owned = ownedToMapped[i].first;
     std::size_t mapped = ownedToMapped[i].second;

     std::vector<std::size_t>::iterator itr 
        = std::find(locallyRequiredIds->begin(),locallyRequiredIds->end(),owned);

     if(itr!=locallyRequiredIds->end())
       *itr = mapped;
     else 
       unusedOwnedToMapped.push_back(ownedToMapped[i]);
  }

  // build  a unique vector of locally matched IDs
  std::vector<std::size_t> unique_locallyRequiredIds;
  {
     std::set<std::size_t> s;
     s.insert(locallyRequiredIds->begin(),locallyRequiredIds->end()); 
     unique_locallyRequiredIds.insert(unique_locallyRequiredIds.begin(),s.begin(),s.end());
  }

  // next line requires communication
  Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds
        = getGlobalPairing(unique_locallyRequiredIds,*locallyMatchedIds,mesh,failure);

  // this nasty bit of code gurantees that only owned/ghosted IDs will
  // end up in the globallyMatchedIds vector.  It uses a search on
  // only locally owned IDs to make sure that they are locally owned.
  // this is a result (and fix) for some of the complexity incurred by
  // the reverseMap above.
  {
     // fill up set with current globally matched ids (not neccessarily owned/ghosted)
     std::set<std::pair<std::size_t,std::size_t> > gmi_set;
     gmi_set.insert(globallyMatchedIds->begin(),globallyMatchedIds->end()); 
   
     // for each globally matched ID, update IDs from the previous 
     // run (i.e. from ownedToMapped) using the reverseMap
     for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
        std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
        const std::vector<std::size_t> & others = reverseMap[pair.first];
        
        // add in reverse map (note other[j] is guranteed to be local to this processor
        // if it was when ownedToMapped was passed in)
        for(std::size_t j=0;j<others.size();j++) 
           gmi_set.insert(std::make_pair(others[j],pair.second)); 
   
        // remove ids that are not ghosted/owned by this processor
        if(std::find(saved_locallyRequiredIds.begin(),
                     saved_locallyRequiredIds.end(),
                     pair.first)==saved_locallyRequiredIds.end()) {
           gmi_set.erase(pair);
        }
     }

     // clear old data, and populate with new matched ids
     globallyMatchedIds->clear();
     globallyMatchedIds->insert(globallyMatchedIds->begin(),gmi_set.begin(),gmi_set.end());
  }

  // now you have a pair of ids that maps ids on the left required by this processor 
  // to ids on the right
  globallyMatchedIds->insert(globallyMatchedIds->end(),unusedOwnedToMapped.begin(),unusedOwnedToMapped.end());

  return globallyMatchedIds;
}

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

     Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
     out.setShowProcRank(true);
     out.setOutputToRootOnly(-1);

     out << "Not all sides matched expect failure: \n" << e.what() << std::endl;
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
   TEUCHOS_TEST_FOR_EXCEPTION(checkProb,std::logic_error,
                      "getLocallyMatchedSideIds: checkProb failed");
   TEUCHOS_TEST_FOR_EXCEPTION(local_side_ids.size()!=result->size(),std::logic_error,
                      "getLocallyMatchedSideIds: not all local side ids are satisfied!");

   return result;
}

} // end periodic_helpers
} // end panzer_stk
