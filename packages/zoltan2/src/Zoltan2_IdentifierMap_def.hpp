// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_IDENTIFIERMAP_DEF_HPP_
#define _ZOLTAN2_IDENTIFIERMAP_DEF_HPP_

/*! \file Zoltan2_IdentifierMap.hpp
*/

#include <stdexcept>
#include <vector>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Map.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_ArrayView.hpp>

namespace Z2
{

template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  IdentifierMap<AppGID,AppLID,GNO,LNO>::IdentifierMap(
    Teuchos::RCP<const Teuchos::Comm<int> > &in_comm, 
    Teuchos::RCP<std::vector<AppGID> > &gids, 
    Teuchos::RCP<std::vector<AppLID> > &lids) 
         : comm(in_comm), myGids(gids), myLids(lids), 
           globalMap(NULL), gidList(NULL), lidList(NULL) { 


    GNO numIds[2] = {static_cast<GNO>(gids->size()), static_cast<GNO>(lids->size())};
    GNO globalNumIds[2];

    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, int(2), numIds, globalNumIds);

    TEST_FOR_EXCEPTION((globalNumIds[1] > 0) && (globalNumIds[0] != globalNumIds[1]), 
                        std::runtime_error,
                       "local IDs are provided but number of global IDs"
                       " does not equal number of local IDs");

#ifdef APPGID_IS_NOT_GNO

    // We can't use the application global ID as our global number.  We'll assign
    //   unique global numbers to each global ID.  
    //
    //   1. We need to store the global numbers in a distributed map because in most 
    //      problems the Model will need this.
    //
    //   2. We need a distributed vector where the value for a given global number
    //      is the application global ID.
    //
    //   3. We need a local hash table to look up the global number for a global ID.
    //
    //   4. If local IDs were supplied by the application, we need a vector analygous to #2
    //        and hash table analygous to #3.


    consecutive = true;   // Tpetra::Map gives processes consecutive GNOs
    base = 0;

    globalMap = Teuchos::rcp(new Tpetra::Map<LNO, GNO>(globalNumIds[0], numIds[0], comm));

    Teuchos::ArrayView<AppGID> gidArray(*gids);

    gidList = new Tpetra::Vector<AppGID, LNO, GNO>(globalMap, gidArray);

    if (globalNumIds[1] > 0){

      Teuchos::ArrayView<AppLID> lidArray(*lids);

      lidList = new Tpetra::Vector<AppLID, LNO, GNO>(globalMap, lidArray);
    }

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
template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  IdentifierMap<AppGID,AppLID,GNO,LNO>::IdentifierMap() : gidList(NULL), lidList(NULL){ }

  /*! Destructor */
template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  IdentifierMap<AppGID,AppLID,GNO,LNO>::~IdentifierMap() 
  {
    if (gidList) delete gidList;
    if (lidList) delete lidList;
  }

  /*! Copy Constructor */
template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  IdentifierMap<AppGID,AppLID,GNO,LNO>::IdentifierMap(const IdentifierMap &id)
  {
  }

  /*! Assignment operator */
template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  IdentifierMap<AppGID,AppLID,GNO,LNO> &
    IdentifierMap<AppGID,AppLID,GNO,LNO>::operator=(const IdentifierMap<AppGID,AppLID,GNO,LNO> &id)
  {
  }

  /*! Set or reset the communicator*/
template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  void 
    IdentifierMap<AppGID,AppLID,GNO,LNO>::setComm(Teuchos::RCP<const Teuchos::Comm<int> > &in_comm) 
  {
    comm = in_comm;
  }

  /*! Set or reset application global IDs for this process*/
template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  void 
    IdentifierMap<AppGID,AppLID,GNO,LNO>::setGlobalIds(Teuchos::RCP<std::vector<AppGID> > &ids)
  {
    myGids = ids;
  }

  /*! Set or reset application local IDs for this process*/
template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  void 
    IdentifierMap<AppGID,AppLID,GNO,LNO>::setLocalIds(Teuchos::RCP<std::vector<AppLID> > &ids)
  {
    myLids = ids;
  }
}   // end namespace Z2

#endif /* _ZOLTAN2_IDENTIFIERMAP_DEF_HPP_ */
