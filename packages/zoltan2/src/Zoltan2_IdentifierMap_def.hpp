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

#ifndef _ZOLTAN2_IDENTIFIERMAP_DEF_HPP_
#define _ZOLTAN2_IDENTIFIERMAP_DEF_HPP_

/*! \file Zoltan2_IdentifierMap.hpp
*/

#include <stdexcept>
#include <vector>
#include <Tpetra_Vector.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_TestForException.hpp>

namespace Z2
{

template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  IdentifierMap<AppGID,AppLID,GNO,LNO>::IdentifierMap(Teuchos::Comm<GNO> &in_comm, std::vector<AppGID> &gids, std::vector<AppLID> &lids) 
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
template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  IdentifierMap<AppGID,AppLID,GNO,LNO>::IdentifierMap() : comm(), gidList(NULL), lidList(NULL)){ }

  /*! Destructor */
template<typename AppGID, typename AppLID, typename GNO, typename LNO>
  IdentifierMap<AppGID,AppLID,GNO,LNO>::~IdentifierMap() 
  {
    if (gidList) free [] gidList;
    if (lidList) free [] lidList;
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

}   // end namespace Z2

#endif /* _ZOLTAN2_IDENTIFIERMAP_DEF_HPP_ */
