// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_STK_PeriodicBC_Matcher_hpp__
#define __Panzer_STK_PeriodicBC_Matcher_hpp__

#include "Teuchos_Tuple.hpp"
#include "Teuchos_RCP.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"

#ifdef PANZER_HAVE_STKSEARCH
#include "stk_search/CoarseSearch.hpp"
#endif

namespace panzer_stk {

/** These functions are utilities to support the implementation of
  * periodic boundary conditions.  They should not be used by externally
  * as their interface is likely to change.
  */
namespace periodic_helpers {

#ifdef PANZER_HAVE_STKSEARCH
  // Copied from PeriodicBoundarySearch
  typedef stk::search::IdentProc<stk::mesh::EntityKey> SearchId;
  typedef stk::search::Sphere<double> Sphere;
  typedef std::vector< std::pair<Sphere,SearchId> > SphereIdVector;
  typedef std::vector<std::pair<SearchId,SearchId> > SearchPairVector;
  typedef std::vector<std::pair<stk::mesh::EntityKey,stk::mesh::EntityKey> > SearchPairSet;
#endif

   /** Construct the vector pair (similar to <code>getLocallyMatchedPair</code>)
     * usign specified side sets, mesh object, and matcher object. This
     * is primarily a utility function.
     */
   template <typename Matcher>
   Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >
   matchPeriodicSides(const std::string & left,const std::string & right,
                     const STK_Interface & mesh,
                     const Matcher & matcher, const std::string type_ = "coord");

   template <typename Matcher>
   Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >
   matchPeriodicSides(const std::string & left,const std::string & right,
                     const STK_Interface & mesh,
                     const Matcher & matcher,
                     const std::vector<std::pair<std::size_t,std::size_t> > & current, const std::string type_ = "coord");

   /** This returns all the global IDs and coordinates for 
     * a particular side. By "all" that means across all processors.
     */
   std::pair<Teuchos::RCP<std::vector<std::size_t> >,
             Teuchos::RCP<std::vector<Teuchos::Tuple<double,3> > > >
   getSideIdsAndCoords(const STK_Interface & mesh,
                       const std::string & sideName, const std::string type_ = "coord");
 
   /** This returns the locally owned global IDs and coordinates for 
     * a particular side. 
     */
   std::pair<Teuchos::RCP<std::vector<std::size_t> >,
             Teuchos::RCP<std::vector<Teuchos::Tuple<double,3> > > >
   getLocalSideIdsAndCoords(const STK_Interface & mesh,
                            const std::string & sideName, const std::string type_ = "coord");
 
   /** This returns the locally resident (includes ghosted) global IDs
     * for a particular side. 
     */
   Teuchos::RCP<std::vector<std::size_t> >
   getLocalSideIds(const STK_Interface & mesh,
                   const std::string & sideName, const std::string type_ = "coord");
 
   /** Determine a map from the specified side to the set of coordinates
     * and Ids passed in. A vector of pairs that maps from (passed in gids)->(locally owned gids)
     * is returned.
     */
   template <typename Matcher>
   Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >
   getLocallyMatchedSideIds(const std::vector<std::size_t> & side_ids,
                            const std::vector<Teuchos::Tuple<double,3> > & side_coords,
                            const STK_Interface & mesh,
                            const std::string & sideName,const Matcher & matcher, const std::string type_ = "coord");
 
   /** Builds a vector of local ids and their matching global indices.
     * This requires a previously discovered vector of pairs of locally matched
     * ids to distribute. This vector comes from the getLocallyMatchedSideIds.
     *
     * \param[in] locallyRequiredIds Global IDs required by this processor to find
     *                               a matched pair. This condition was not satisfied
     *                               locally.
     * \param[in] locallyMatchedIds Global IDs that this processor has matched.
     */
   Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >
   getGlobalPairing(const std::vector<std::size_t> & locallyRequiredIds,
                    const std::vector<std::pair<std::size_t,std::size_t> > & locallyMatchedIds,
                    const STK_Interface & mesh,bool failure);

#ifdef PANZER_HAVE_STKSEARCH
   /** Construct the local search vector containing the coordinates 
     * and entity keys of the requested entity (and potentially ghosted entities). 
     * Optionally, a list of previously matched sides can be included to avoid matching the same node multiple times.
     * In this case, a list of these nodes can be returned for remapping in the multiperiodic case.
     * 
     * \param[in] mesh STK mesh interface
     * \param[out] searchVector Search vector (to be filled)
     * \param[in] error Tolerance in the coordinate position for matching
     * \param[in] sideName Side name to search for
     * \param[in] type_ Entity type
     * \param[in] getGhostedIDs Flag for including ghosted nodes in the search vector
     * \param[in] matchedSides Vector of the sides that have been previously matched. Should be specific to the \p type_ requested.
     * \param[out] potentialIDsToRemap IDs of the nodes which have been previously matched and are part of the requested side
     */
  void fillLocalSearchVector(const STK_Interface & mesh, SphereIdVector & searchVector, const double & error,
                             const std::string & sideName, const std::string & type_, const bool & getGhostedIDs, 
                             const std::vector<std::string> & matchedSides, std::vector<SearchId> & potentialIDsToRemap);

  void fillLocalSearchVector(const STK_Interface & mesh, SphereIdVector & searchVector, const double & error,
                             const std::string & sideName, const std::string & type_, const bool & getGhostedIDs = false);

  /** Compute the centroid of a given side
   * \param[in] mesh STK mesh interface
   * \param[in] sideName Side name to search
   * \returns The vector (3-D) of the computed centroid
   */

  const std::vector<double> computeGlobalCentroid(const STK_Interface & mesh, const std::string & sideName);

  /** Apply the given matcher's transform to the sideA sphere vector.
   * Side A's and side B's sphere vectors can then be compared directly for coordinate intersection.
   * \param[inout] searchVectorSideA Seach vector to transform
   * \param[in] matcher Periodic matcher with transform rule
   * \param[in] centroidSideB Centroid of side B
   */
  template<typename Matcher> void
  transformLocalSearchVector(SphereIdVector & searchVectorSideA, const Matcher & matcher, const std::vector<double> & centroidSideB );

  /** Match the nodes on side A with nodes on side B.
   * If a node has been previously matched, it does not get matched a second time.
   * Rather, its mapping is updated to respect the mapping between its match and a third node.
   * \param[in] sideA Side name for side A
   * \param[in] sideB Side name for side B 
   * \param[in] mesh STK mesh interface
   * \param[in] matcher The periodic matcher object (contains matcher rule)
   * \param[in] matchedSides List of the previously matched sides (side A's)
   * \param[in] previousMatches Previously matched nodes and their matches
   * \param[in] type_ Entity type
   * 
   * \returns An (updated, if necessary) mapping of local periodic nodes and their matches
   */

  template <typename Matcher>
  Teuchos::RCP<std::vector<std::pair<size_t,size_t> > > 
  matchPeriodicSidesSearch(const std::string & sideA,const std::string & sideB,
                           const STK_Interface & mesh,
                           const Matcher & matcher, const std::vector<std::string> & matchedSides,
                           const std::vector<std::pair<size_t,size_t> > & previousMatches,
                           const std::string type_ = "coord");
 
  template <typename Matcher>
  Teuchos::RCP<std::vector<std::pair<size_t,size_t> > > 
  matchPeriodicSidesSearch(const std::string & sideA, const std::string & sideB,
                           const STK_Interface & mesh,
                           const Matcher & matcher, const std::string type_ = "coord");

  /** Update the current A to B map to include previous matches and remap IDs as appropriate.
   * This should be called when an entity type has already been matched.
   * Otherwise, \ref appendMapping may be appropriate.
   * \param[inout] currentMatches The current A to B map
   * \param[in] previousMatches A to B map from earlier periodic matches
   * \param[in] IDsToRemap Multiply periodic IDs that need to be remapped
   * \param[in] mesh STK mesh interface
   */

  void updateMapping(Teuchos::RCP<std::vector<std::pair<size_t,size_t> > > & currentMatches,
                     const std::vector<std::pair<size_t,size_t> > & previousMatches,
                     const std::vector<SearchId> & IDsToRemap, const STK_Interface & mesh);

  /** Update the current A to B map to append previous matches.
   * This should be called when different entity type has already been matched.
   * Otherwise, \ref updateMapping may be appropriate.
   * \param[inout] currentMatches The current A to B map
   * \param[in] previousMatches A to B map from earlier periodic matches
   */

  void appendMapping(Teuchos::RCP<std::vector<std::pair<size_t,size_t> > > & currentMatches,
                     const std::vector<std::pair<size_t,size_t> > & previousMatches);
#endif
} // end periodic_helpers

/** Simply returns a vector of pairs that match
  * the IDs owned by this processor to their
  * matching IDs on the periodic boundary.
  *
  * Notice that the matched boundaries are not 
  * specified by this object. This is done in the 
  * inherited class.
  */
class PeriodicBC_MatcherBase {
public:
   virtual ~PeriodicBC_MatcherBase() {}

   /** Simply returns a vector of pairs that match
     * the IDs owned by this processor to their
     * matching IDs on the periodic boundary.
     *
     * \returns A vector of pairs. The first entry in the
     *          pair is the global node ID of a node used
     *          on this processor. The second is the global
     *          node ID (not necessarily on this processor)
     *          that replaces it.
     */
   virtual 
   Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >
   getMatchedPair(const STK_Interface & mesh,
                  const Teuchos::RCP<const std::vector<std::pair<std::size_t,std::size_t> > > & currentState = Teuchos::null
                  ) const = 0;

#ifdef PANZER_HAVE_STKSEARCH
  /** Simply returns a vector of pairs that match
   * the IDs owned by this processor to their
   * matching IDs on the periodic boundary.
   * 
   * \param[in] mesh STK mesh interface
   * \param[in] matchedSides List of the previously matched sides
   * \param[in] currentState RCP to the existing map (may be null)
   * 
   * \note This version of getMatchedPair uses STKSearch to accelerate the matching process.
   *
   * \returns A vector of pairs. The first entry in the
   *          pair is the global node ID of a node used
   *          on this processor. The second is the global
   *          node ID (not necessarily on this processor)
   *          that replaces it. 
   */
   virtual 
   Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >
   getMatchedPair(const STK_Interface & mesh, const std::vector<std::string> & matchedSides,
                  const Teuchos::RCP<const std::vector<std::pair<std::size_t,std::size_t> > > & currentState = Teuchos::null
                  ) const = 0;
#endif 

   /** Return a one line string that describes this periodic
     * boundary condition.
     */
   virtual std::string getString() const = 0;

   /** Return a one line string that describes the type of periodic
     * boundary condition.
     */
   virtual std::string getType() const = 0;

   /// Returns the sideset name for the left side.
   virtual std::string getLeftSidesetName() const = 0;

   /// Returns the sideset name for the right side.
   virtual std::string getRightSidesetName() const = 0;

  /// Attempts to cast the underlying matcher to type T. Returns nullptr if cast fails.
  template<typename T>
  const T* getAs() const {return dynamic_cast<const T*>(this);}
};

/** Default implementation class for the periodic boundary conditions.
  * This simply takes a <code>Matcher</code> object as the coordinate matching
  * operation. For details on the <code>Matcher</code> see the
  * <code>CoordMatcher</code> struct.
  */
template <typename Matcher>
class PeriodicBC_Matcher : public PeriodicBC_MatcherBase {
public:
   PeriodicBC_Matcher(const std::string & left, const std::string & right,const Matcher & matcher, const std::string type = "coord")
      : left_(left), right_(right), matcher_(matcher), type_(type) {}
   PeriodicBC_Matcher(const PeriodicBC_Matcher & src)
      : left_(src.left_), right_(src.right_), matcher_(src.matcher_), type_(src.type_) {}

   /** Simply returns a vector of pairs that match
     * the IDs owned by this processor to their
     * matching IDs on the periodic boundary.
     *
     * \returns A vector of pairs. The first entry in the
     *          pair is the global node ID of a node used
     *          on this processor. The second is the global
     *          node ID (not necessarily on this processor)
     *          that replaces it.
     */
   Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >
   getMatchedPair(const STK_Interface & mesh,
                  const Teuchos::RCP<const std::vector<std::pair<std::size_t,std::size_t> > > & currentState = Teuchos::null
                  ) const
   { 
      if(currentState==Teuchos::null) 
         return periodic_helpers::matchPeriodicSides(left_,right_,mesh,matcher_,type_); 
      else
         return periodic_helpers::matchPeriodicSides(left_,right_,mesh,matcher_,*currentState,type_); 
   }

#ifdef PANZER_HAVE_STKSEARCH
  /** Simply returns a vector of pairs that match
   * the IDs owned by this processor to their
   * matching IDs on the periodic boundary.
   * 
   * \param[in] mesh STK mesh interface
   * \param[in] matchedSides List of the previously matched sides
   * \param[in] currentState RCP to the existing map (may be null)
   * 
   * \note This version of getMatchedPair uses STKSearch to accelerate the matching process.
   *
   * \returns A vector of pairs. The first entry in the
   *          pair is the global node ID of a node used
   *          on this processor. The second is the global
   *          node ID (not necessarily on this processor)
   *          that replaces it. 
   */
   Teuchos::RCP<std::vector<std::pair<std::size_t,std::size_t> > >
   getMatchedPair(const STK_Interface & mesh, const std::vector<std::string> & matchedSides,
                  const Teuchos::RCP<const std::vector<std::pair<std::size_t,std::size_t> > > & currentState = Teuchos::null
                  ) const
   { 
      if(currentState==Teuchos::null) 
         return periodic_helpers::matchPeriodicSidesSearch(left_,right_,mesh,matcher_,type_); 
      else
         return periodic_helpers::matchPeriodicSidesSearch(left_,right_,mesh,matcher_,matchedSides,*currentState,type_); 
   }
#endif

   std::string getString() const 
   { 
      std::stringstream ss;
      ss << "condition: " << matcher_.getString() << ", sides = [ "
         << "\"" << left_ << "\", "
         << "\"" << right_ << "\" ]";
      return ss.str();
   }

   std::string getType() const 
   {return type_;}

   std::string getLeftSidesetName() const
   {return left_;}

   std::string getRightSidesetName() const
   {return right_;}

   const Matcher& getMatcher() const
   {return matcher_;}

private:
   PeriodicBC_Matcher(); // hidden!

   std::string left_;  // here left & right are stand in names just so
   std::string right_; // that we realize that these boundaries are
                               // opposite of each other.
   Matcher matcher_;

   std::string type_; // type of periodic BC: coord, edge, face

};

/** A simple constructor function for building a matcher object.
  * This prevents the need to directly instantiate the templated 
  * derived class. It is a convenience.
  */
template <typename Matcher>
Teuchos::RCP<PeriodicBC_MatcherBase> 
buildPeriodicBC_Matcher(const std::string & left, const std::string & right, const Matcher & matcher, const std::string type = "coord")
{ return Teuchos::rcp(new PeriodicBC_Matcher<Matcher>(left,right,matcher,type)); }

} // end panzer_stk

#include "Panzer_STK_PeriodicBC_Matcher_impl.hpp"

#ifdef PANZER_HAVE_STKSEARCH
#include "Panzer_STK_PeriodicBC_Search_impl.hpp"
#endif

#endif
