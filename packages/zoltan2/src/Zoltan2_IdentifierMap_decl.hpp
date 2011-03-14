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
#include <Teuchos_Hashtable.hpp>

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

  // These are RCP'd so the Tpetra::Vectors will be automatically deallocated for us.

  Teuchos::RCP<Tpetra::Map<LNO, GNO> > globalMap;
 
  Teuchos::RCP<Tpetra::Vector<AppGID, LNO, GNO> > gidList;   // if the GNOs are not AppGIDs
  Teuchos::RCP<Tpetra::Vector<AppLID, LNO, GNO> > lidList;   // if local IDs were supplied

  Teuchos::RCP<Teuchos::Hashtable<int, GNO> >  gidHash;   // if GNOs are not AppGIDs

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


  /*! Constructor 
      This constructor does not need to be called by all processes.
   */
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

  /*! Return true if we are using the application global IDs 
   *  for our internal global numbers 
   */
  bool gnosAreGids();

  /*! Map application global IDs to internal global numbers or vice versa.

      This is a local call.  If gid is a vector of application global IDs, then
      gno will be set to the corresponding internal global numbers.  If the
      gno vector contains internal global numbers, they will be translated
      to application global IDs.  The application global IDs must be from
      those supplied by this process.
   */
  void gidTranslate(std::vector<AppGID> &gid, std::vector<GNO> &gno);

  /*! Map application local IDs to internal global numbers or vice versa.

      This is a local call.  If lid is a vector of application local IDs, then
      gno will be set to the corresponding internal global numbers.  If the
      gno vector contains internal global numbers, they will be translated
      to application local IDs.  The application local IDs must be from
      those supplied by this process.
   */
  void lidTranslate(std::vector<AppLID> &lid, std::vector<GNO> &gno);

  /*! Returns a smart pointer to the global number map.  

   *  May be needed by the Model.  Should not be required by the InputAdapter.
   */
  Teuchos::RCP<Tpetra::Map<LNO, GNO> > &getGlobalMap();

  /*! Map application global IDs to internal global numbers or vice versa.

      All processes must call this.  The global IDs or global numbers
      supplied may belong to another process.  This method will fill
      in the empty vector with the corresponding id, and will fill the
      proc vector with the owner of the global ID.
   */
  void gidGlobalTranslate(std::vector<AppGID> &gid, std::vector<GNO> &gno, std::vector<int> &proc);
};

}   // end namespace Z2

#endif /* _ZOLTAN2_IDENTIFIERMAP_DECL_HPP_ */
