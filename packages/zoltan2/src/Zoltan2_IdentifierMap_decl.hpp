// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_IDENTIFIERMAP_DECL_HPP_
#define _ZOLTAN2_IDENTIFIERMAP_DECL_HPP_

/*! \file Zoltan2_IdentifierMap.hpp
*/

#include <vector>
#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Hashtable.hpp>

namespace Z2
{

/*! Z2::IdentifierMap
    \brief An IdentifierMap manages a global space of unique object identifiers.

    TODO - trim down comments and code once we get it all straight
           replace new/delete with memory wrappers
           exception handling
           use Kokkos node 
*/

template<typename AppLID, typename AppGID, typename LNO=int, 
  typename GNO=AppGID> 
    class IdentifierMap{

private:

  // Input communicator

  Teuchos::RCP<const Teuchos::Comm<int> > _comm;

  // Application global and local IDs

  Teuchos::RCP<std::vector<AppGID> > _myGids; 
  Teuchos::RCP<std::vector<AppLID> > _myLids;

  // In the case of consecutive ordinal application global IDs,
  // gnoDist[p] is the first global number on process p, and
  // we don't need the _gidHash.

  Teuchos::ArrayRCP<GNO> _gnoDist;

  // A hash table from application global ID key to our local index.

  Teuchos::RCP<Teuchos::Hashtable<std::string, LNO> >  _gidHash;

  // A hash table from application local ID key to our local index.

  Teuchos::RCP<Teuchos::Hashtable<std::string, LNO> >  _lidHash;

  GNO _globalNumberOfIds;
  LNO _localNumberOfIds;
  bool _haveLocalIds;
  int _myRank;
  int _numProcs;

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

  /*! Initialize object if not done in the constructor */
  void initialize(Teuchos::RCP<const Teuchos::Comm<int> > &in_comm,
                  Teuchos::RCP<std::vector<AppGID> > &gids,
                  Teuchos::RCP<std::vector<AppLID> > &lids);

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

  /*! Map application global IDs to internal global numbers or vice versa.

      All processes must call this.  The global IDs or global numbers
      supplied may belong to another process.  This method will fill
      in the empty vector with the corresponding id, and will fill the
      proc vector with the owner of the global ID.
   */
  void gidGlobalTranslate(const std::vector<AppGID> &gid, 
    std::vector<GNO> &gno, std::vector<int> &proc);
};

}   // end namespace Z2

#endif /* _ZOLTAN2_IDENTIFIERMAP_DECL_HPP_ */
