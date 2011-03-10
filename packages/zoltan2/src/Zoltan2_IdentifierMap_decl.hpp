// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// Questions? Contact Lee Ann Riesen (lriesen@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_IDENTIFIERMAP_DECL_HPP_
#define _ZOLTAN2_IDENTIFIERMAP_DECL_HPP_

/*! \file Zoltan2_IdentifierMap.hpp
*/

#include <stdexcept>
#include <vector>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Map.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_TestForException.hpp>

/*! Z2
    \brief A namespace for objects that are not part of the Zoltan2 user interface.
*/

namespace Z2
{

/*! Z2::IdentifierMap
    \brief An IdentifierMap manages a global space of unique object identifiers.

    The Zoltan2 caller may use arbitrary data types for global ID (AppGID) and 
    local ID (AppLID).  For example, if the objects being partitioned are matrix 
    non-zeros, the AppGID may be an (i,j) pair.  The AppLID is optional and exists
    for the convenience of the caller.  It may be for example an index into an array, 
    or a pointer to memory.  The application may obtain a solution in terms of the 
    AppLID (if it was provided) or the AppGID.

    The Identifier map is templated on the AppGID, AppLID, GNO and LNO.  The GNO
    and LNO are the global identifier and local index types used internally by Zoltan2.  

    The GNO must be a signed or unsigned char, int or long.  Where long long is
    available it can be a signed long long.  (It must be type valid for use as a 
    Tpetra GlobalOrdinal.) The GNO type defaults to the AppGID type.  If the 
    AppGID is not a GlobalOrdinal type then the Caller must 
    specify a valid GNO type that is large enough to enumerate all of the application's 
    global IDs.  It is more efficient to use AppGIDs that are GlobalOrdinal types because
    if they are not, Zoltan2 must translate back and forth between the application
    global ID and its internal global identifier.

    Before including this header file in a compilation, the macro APPGID_IS_NOT_GNO
    must be defined if the application's global ID type is not the same as the GNO type. 

    The LNO defaults to an int.  It is used to index and count local objects.  If there
    is some reason that Zoltan2 could save memory by using a smaller data type, or
    requires a larger data type to count local objects, then the caller should
    define a different LNO type.

    TODO - trim down comments and code once we get it all straight
           replace new/delete with memory wrappers
           exception handling
           use Kokkos node 
*/

template<typename AppGID, typename AppLID, typename GNO=AppGID, typename LNO=int>
  class IdentifierMap{

private:

  Teuchos::RCP<const Teuchos::Comm<int> > comm;        // Get this from InputAdapter

  Teuchos::RCP<std::vector<AppGID> > myGids;  // This proc's global IDs from InputAdapter
  Teuchos::RCP<std::vector<AppLID> > myLids;  // This proc's local IDs from InputAdapter

  Teuchos::RCP<Tpetra::Map<LNO, GNO> > globalMap;

  Tpetra::Vector<AppGID, LNO, GNO> *gidList;   // needed if the GNOs are not the AppGIDs
  Tpetra::Vector<AppLID, LNO, GNO> *lidList;   // needed if local IDs were supplied

  bool consecutive;
  GNO base;

public:

  /*! Constructor - Must be called by all processes 
   *
   *  We are not taking a Kokkos node here.  In the future we should.
   *
   *  Problem: we have to template the communicator on int, because Tpetra_Map
   *  only takes Teuchos::Comm<int>. 
   */

  explicit IdentifierMap(Teuchos::RCP<const Teuchos::Comm<int> > &in_comm, 
                         Teuchos::RCP<std::vector<AppGID> > &gids, 
                         Teuchos::RCP<std::vector<AppLID> > &lids);


  /*! Constructor */
  IdentifierMap();

  /*! Destructor */
  ~IdentifierMap();

  /*! Copy Constructor */
  IdentifierMap(const IdentifierMap &id);

  /*! Assignment operator */
  IdentifierMap &operator=(const IdentifierMap &id);

  /*! Set or reset the communicator*/
  void setComm(Teuchos::RCP<const Teuchos::Comm<int> > &in_comm);

  /*! Set or reset application global IDs for this process*/
  void setGlobalIds(Teuchos::RCP<std::vector<AppGID> > &ids);

  /*! Set or reset application local IDs for this process*/
  void setLocalIds(Teuchos::RCP<std::vector<AppLID> > &ids);
};

}   // end namespace Z2

#endif /* _ZOLTAN2_IDENTIFIERMAP_DECL_HPP_ */
