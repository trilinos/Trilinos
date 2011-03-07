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

#ifndef _ZOLTAN2_IDENTIFIERMAP_HPP_
#define _ZOLTAN2_IDENTIFIERMAP_HPP_

/*! \file Zoltan2_IdentifierMap.hpp
*/

#include <stdexcept>
#include <vector>
#include <Tpetra_Vector.hpp>
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

  Teuchos::Comm<GNO> &comm;

  Tpetra::Vector<AppGID, LNO, GNO> *gidList;
  Tpetra::Vector<AppLID, LNO, GNO> *lidList;

#ifdef APPGID_IS_NOT_GNO
  bool GidGnoMap(true);
  bool consecutive(true);
  GNO base(0);
  std::vector<GNO> gnoDist;
#else
  bool GidGnoMap(false);
  bool consecutive;
  GNO base;
#endif

public:

  /*! Constructor */
  IdentifierMap(Teuchos::Comm<GNO> &in_comm, std::vector<AppGID> &gids, std::vector<AppLID> &lids) 
         : comm(in_comm), gidList(NULL), lidList(NULL)){ 

    Teuchos::TEST_FOR_EXCEPTION((gids.size() > 0), std::logic_error, "empty gid list")

    LNO numGids = static_cast<LNO>(gids.size());
    GNO globalNumGids(0);

    Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, GNO(1), &numGids, &globalNumGids);

#ifdef APPGID_IS_NOT_GNO

    int nprocs = comm.getSize();
    GNO *num = new [nprocs] GNO;
    
    Teuchos::scan(comm, Teuchos::REDUCE_SUM, GNO(1), &numGids, num);

    gnoDist.reserve(nprocs+1);
    gnoDist[0] = 0;

    for (int i=1; i < nprocs+1; i++)
      gnoDist[i] = gnoDist[i-1] + num[i-1];

    delete [] num;

#else
    GNO min(0), max(0), globalMin(0), globalMax(0);
    min = max = static_cast<GNO>(gids[0]);

    for (LNO i=1; i < numGids; i++){
      if (gids[i] < min)
        min = static_cast<GNO>(gids[i]);
      else if (gids[i] > max)
        max = static_cast<GNO>(gids[i]);
    }

    Teuchos::reduceAll(comm, Teuchos::REDUCE_MIN, GNO(1), &min, &globalMin);
    Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, GNO(1), &max, &globalMax);

    if (globalMax - globalMin + 1 == globalNumGids){
      consecutive=true;
      base = globalMin;
    }
#endif

  }

  /*! Constructor */
  IdentifierMap() : comm(), gidList(NULL), lidList(NULL)){ }

  /*! Destructor */
  ~IdentifierMap(){}

  /*! Copy Constructor */
  IdentifierMap(const IdentifierMap &id){
  }

  /*! Assignment operator */
  IdentifierMap &operator=(const IdentifierMap &id){
  }

#if 0

  /*! Methods to duplicate functionality of Zoltan_DD_*   ?? Does this all belong here??
   */

  void update();
  void find();
  void gid2gno();
  void gno2lno();
  void etc();

  /*! Assign global numbers randomly across processes.
   */
  bool randomize_ownership();
#endif
};

}   // end namespace Z2

#endif /* _ZOLTAN2_IDENTIFIERMAP_HPP_ */
