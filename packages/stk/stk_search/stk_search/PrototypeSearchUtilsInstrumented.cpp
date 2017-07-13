
#include <stk_util/parallel/Parallel.hpp>    // for parallel_machine_size, etc
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommSparse.hpp>  // for comm_recv_sizes
#include <stk_util/parallel/MPI.hpp>

#include <iomanip>

#include <stk_search/CommonSearchUtilsInstrumented.hpp>
#include <stk_search/PrototypeSearchUtilsInstrumented.hpp>


namespace stk {
namespace search {
namespace experimental {


std::ostream &operator<<(std::ostream &os, const stk::search::experimental::BoxedTilingIndices &bti) {
  os << "{(" << bti.xIdx << " " << bti.yIdx << " " << bti.zIdx << ") " << bti.box << "}";
  return os;
}


// Greedy brute-force algorithm, just the first thing that came to mind.  Don't know
// if the output is actually optimal, which needs to be better defined anyway.
void optimizeTile3D(SuperTile3D &tile, int maxTiles)
{
  double inputTileLen[3] = {tile.xLen, tile.yLen, tile.zLen};
  double tileLen[3]      = {tile.xLen, tile.yLen, tile.zLen};
  int    numTiles[3]     = {tile.numXTiles, tile.numYTiles, tile.numZTiles};

  bool refinable[3] = {true, true, true};
  int mostRefinable = 0;
  while (0 <= mostRefinable) {
    for (int i = 0; i < 3; ++i) {
      int prod = (  (numTiles[0] + (i == 0 ? 1 : 0))
                  * (numTiles[1] + (i == 1 ? 1 : 0))
                  * (numTiles[2] + (i == 2 ? 1 : 0)) );
      if (prod > maxTiles) {
        refinable[i] = false;
      }
    }
    mostRefinable = -1;
    double biggest = 0;
    for (int i = 0; i < 3; ++i) {
      if (!refinable[i]) {
        continue;
      }
      if (tileLen[i] > biggest) {
        biggest = tileLen[i];
        mostRefinable = i;
      }
    }
    if (mostRefinable >= 0) {
      numTiles[mostRefinable] += 1;
      tileLen[mostRefinable] = inputTileLen[mostRefinable] / numTiles[mostRefinable];
    }
  }

  tile.xLen = tileLen[0];
  tile.yLen = tileLen[1];
  tile.zLen = tileLen[2];
  tile.numXTiles = numTiles[0];
  tile.numYTiles = numTiles[1];
  tile.numZTiles = numTiles[2];
}

struct TilingIndicesMsg {
  int xIdx, yIdx, zIdx;
  int pRank;
};

void findGhostingRegionResidents(const SuperTile3D &tilingPattern, MPI_Comm comm, int mpiRank, int mpiSize,
                                 const std::vector<TilingIndices> &tilesOccupied,
                                 ResidencyMapT &residencyMap)
{
  // std::cout << "p_" << mpiRank << ": " << "fGRR -- tilesOccupied.size() returns " << tilesOccupied.size() << std::endl;
  residencyMap.clear();

  stk::CommSparse commsparse(comm);

  for(int phase = 0; phase < 2; ++phase) {
    for (const TilingIndices &tileId : tilesOccupied) {
      int proc = getGhostingRegionRank(tileId, tilingPattern, mpiSize);
      if (mpiRank != proc) {
        stk::CommBuffer &procBuff = commsparse.send_buffer(proc);
        TilingIndicesMsg msg{tileId.xIdx, tileId.yIdx, tileId.zIdx, mpiRank};
        procBuff.pack<TilingIndicesMsg>(msg);
      }
      else if (phase == 1) {
        residencyMap[tileId].push_back(mpiRank);
      }
    }
    if (phase == 0) {
      commsparse.allocate_buffers();
    }
    else {
      commsparse.communicate();
    }
  }

  for(int proc=0; proc < mpiSize; ++proc) {
    if (proc == mpiRank) {
      continue;
    }
    stk::CommBuffer& procBuff = commsparse.recv_buffer(proc);
    int numItems = procBuff.remaining()/sizeof(TilingIndicesMsg);
    for (int i = 0; i < numItems; ++i) {
      TilingIndicesMsg msg;
      procBuff.unpack<TilingIndicesMsg>(msg);
      TilingIndices tileId{msg.xIdx, msg.yIdx, msg.zIdx};
      residencyMap[tileId].push_back(msg.pRank);
    }
  }
}

struct BoxedTilingIndicesMsg : public TilingIndicesMsg {
  Box3D box;
};

void findGhostingRegionResidentsDDRobust(const SuperTile3D &tilingPattern, MPI_Comm comm, int mpiRank, int mpiSize,
                                         const std::vector<BoxedTilingIndices> &tilesOccupied,
                                         ResidencyBoxMapT &residencyMap)
{
  // std::cout << "p_" << mpiRank << ": " << "fGRR -- tilesOccupied.size() returns " << tilesOccupied.size() << std::endl;

  residencyMap.clear();

  stk::CommSparse commsparse(comm);

  int sendCount = 0;
  for(int phase = 0; phase < 2; ++phase) {
    for (const BoxedTilingIndices &tileId : tilesOccupied) {
      int proc = getGhostingRegionRank(tileId, tilingPattern, mpiSize);
      if (mpiRank != proc) {
        stk::CommBuffer &procBuff = commsparse.send_buffer(proc);
        BoxedTilingIndicesMsg msg;
        msg.xIdx = tileId.xIdx;
        msg.yIdx = tileId.yIdx;
        msg.zIdx = tileId.zIdx;
        msg.pRank = mpiRank;
        msg.box = tileId.box;
        procBuff.pack<BoxedTilingIndicesMsg>(msg);

        if (phase != 0) {
          ++sendCount;
        }
      }
      else if (phase == 1){
        residencyMap[tileId].push_back(OwnedBox3D{mpiRank,tileId.box});
      }
    }
    if (phase == 0) {
      commsparse.allocate_buffers();
    }
    else {
      commsparse.communicate();
    }
  }
  // std::cout  << "p_" << mpiRank << ": (domain) sent " << sendCount << " boxes" << std::endl;

  int recvCount = 0;
  for(int proc=0; proc < mpiSize; ++proc) {
    if (proc == mpiRank) {
      continue;
    }
    stk::CommBuffer& procBuff = commsparse.recv_buffer(proc);
    int numItems = procBuff.remaining()/sizeof(BoxedTilingIndicesMsg);
    for (int i = 0; i < numItems; ++i) {
      BoxedTilingIndicesMsg msg;
      procBuff.unpack<BoxedTilingIndicesMsg>(msg);
      TilingIndices tileId{msg.xIdx, msg.yIdx, msg.zIdx};
      residencyMap[tileId].push_back(OwnedBox3D{msg.pRank, msg.box});

      ++recvCount;
    }
  }
  // std::cout  << "p_" << mpiRank << ": (ghosting region) received " << recvCount << " boxes" << std::endl;
}

enum NeighborMsgTypeEnum {
  DOMAIN_NEIGHBORING_MSG,
  RANGE_NEIGHBORING_MSG
};

struct NeighboringRankMsg {
  NeighborMsgTypeEnum msgType;
  int neighborRank;
};

void witnessAndComputeNeighborRanks(MPI_Comm comm, int mpiSize,
                                    const ResidencyMapT &myDomainResidents,
                                    const ResidencyMapT &myRangeResidents,
                                    std::vector<int> &neighborDomainRanks,
                                    std::vector<int> &neighborRangeRanks)
{
  int myRank = -1;
  MPI_Comm_rank(comm, &myRank);

  ResidencyMapT::const_iterator domResIter = myDomainResidents.begin();
  ResidencyMapT::const_iterator domResEnd  = myDomainResidents.end();
  ResidencyMapT::const_iterator rngResEnd  = myRangeResidents.end();

  std::vector<std::pair<int, int> > neighborPairs;

  //
  // Witness, and per-Ghosting-Region neighbor computation
  //
  for (; domResIter != domResEnd; ++domResIter) {
    const TilingIndices &tileId = domResIter->first;
    const std::vector<int> &domainResidents = domResIter->second;
    ResidencyMapT::const_iterator rngResProbe = myRangeResidents.find(tileId);
    if (rngResProbe != rngResEnd) {
      const std::vector<int> &rangeResidents = rngResProbe->second;
      for (int domainRank : domainResidents) {
        for (int rangeRank : rangeResidents) {
          if (domainRank != rangeRank) {
            neighborPairs.push_back(std::pair<int,int>(domainRank, rangeRank));
          }
        }
      }
    }
  }

  //
  // Communication and finish global neighbor computation.
  //

  stk::CommSparse commsparse(comm);
  for(int phase = 0; phase < 2; ++phase) {
    for (const std::pair<int,int> &neighborPair : neighborPairs) {
      int domainProc = neighborPair.first;
      int rangeProc  = neighborPair.second;
      stk::CommBuffer &domainProcBuff = commsparse.send_buffer(domainProc);
      domainProcBuff.pack<NeighboringRankMsg>(NeighboringRankMsg{RANGE_NEIGHBORING_MSG, rangeProc});
      stk::CommBuffer &rangeProcBuff = commsparse.send_buffer(rangeProc);
      rangeProcBuff.pack<NeighboringRankMsg>(NeighboringRankMsg{DOMAIN_NEIGHBORING_MSG, domainProc});
    }
    if (phase == 0) {
      commsparse.allocate_buffers();
    }
    else {
      commsparse.communicate();
    }
  }
  std::set<int> nbrDomRanks, nbrRngRanks;
  for(int proc=0; proc < mpiSize; ++proc) {
    stk::CommBuffer& procBuff = commsparse.recv_buffer(proc);
    int numItems = procBuff.remaining()/sizeof(NeighboringRankMsg);
    for (int i = 0; i < numItems; ++i) {
      NeighboringRankMsg msg;
      procBuff.unpack<NeighboringRankMsg>(msg);
      if (msg.msgType == DOMAIN_NEIGHBORING_MSG) {
        if (nbrDomRanks.find(msg.neighborRank) == nbrDomRanks.end()) {
          nbrDomRanks.insert(msg.neighborRank);
        }
      }
      else if (nbrRngRanks.find(msg.neighborRank) == nbrRngRanks.end()) {
        nbrRngRanks.insert(msg.neighborRank);
      }
    }
  }
  neighborDomainRanks.clear();
  neighborRangeRanks.clear();
  for (int p : nbrDomRanks) { neighborDomainRanks.push_back(p); }
  for (int p : nbrRngRanks) { neighborRangeRanks.push_back(p); }
}


struct NeighboringProcBoxMsg {
  NeighborMsgTypeEnum msgType;
  OwnedBox3D          rankBox;
};

void witnessAndComputeNeighborRanksDDRobust(MPI_Comm mpiComm, int mpiSize,
                                            const ResidencyBoxMapT &myDomainResidents,
                                            const ResidencyBoxMapT &myRangeResidents,
                                            std::vector<OwnedBox3D> &neighborDomainRanks,
                                            std::vector<OwnedBox3D> &neighborRangeRanks)
{
  int myRank = -1;
  MPI_Comm_rank(mpiComm, &myRank);

  ResidencyBoxMapT::const_iterator domResIter = myDomainResidents.begin();
  ResidencyBoxMapT::const_iterator domResEnd  = myDomainResidents.end();
  ResidencyBoxMapT::const_iterator rngResEnd  = myRangeResidents.end();

  std::vector<std::pair<OwnedBox3D, OwnedBox3D> > neighborPairs;

  //
  // Witness, and per-Ghosting-Region neighbor computation
  //
  for (; domResIter != domResEnd; ++domResIter) {
    const TilingIndices &tileId = domResIter->first;
    const std::vector<OwnedBox3D> &domainResidents = domResIter->second;
    ResidencyBoxMapT::const_iterator rngResProbe = myRangeResidents.find(tileId);
    if (rngResProbe != rngResEnd) {
      const std::vector<OwnedBox3D> &rangeResidents = rngResProbe->second;
      for (const OwnedBox3D &domainBox : domainResidents) {
        for (const OwnedBox3D &rangeBox : rangeResidents) {
          if ((domainBox.owner != rangeBox.owner) && intersects(domainBox.box, rangeBox.box)) {
            neighborPairs.push_back(std::pair<OwnedBox3D, OwnedBox3D>(domainBox, rangeBox));
          }
        }
      }
    }
  }

  // std::cout << "p_" << myRank << ": I witness " << neighborPairs.size() << " neighbor pairs" << std::endl;

  //
  // Communication and finish global neighbor computation.
  //

  stk::CommSparse commsparse(mpiComm);
  for(int phase = 0; phase < 2; ++phase) {
    for (const std::pair<OwnedBox3D,OwnedBox3D> &neighborPair : neighborPairs) {
      OwnedBox3D domainBoxProc = neighborPair.first;
      OwnedBox3D rangeBoxProc  = neighborPair.second;
      int domainProc = domainBoxProc.owner;
      int rangeProc  = rangeBoxProc.owner;
      stk::CommBuffer &domainProcBuff = commsparse.send_buffer(domainProc);
      domainProcBuff.pack<NeighboringProcBoxMsg>(NeighboringProcBoxMsg{RANGE_NEIGHBORING_MSG, rangeBoxProc});
      stk::CommBuffer &rangeProcBuff = commsparse.send_buffer(rangeProc);
      rangeProcBuff.pack<NeighboringProcBoxMsg>(NeighboringProcBoxMsg{DOMAIN_NEIGHBORING_MSG, domainBoxProc});
    }
    if (phase == 0) {
      commsparse.allocate_buffers();
    }
    else {
      commsparse.communicate();
    }
  }
  std::set<OwnedBox3D> nbrDomRanks, nbrRngRanks;
  for(int proc=0; proc < mpiSize; ++proc) {
    stk::CommBuffer& procBuff = commsparse.recv_buffer(proc);
    int numItems = procBuff.remaining()/sizeof(NeighboringProcBoxMsg);
    for (int i = 0; i < numItems; ++i) {
      NeighboringProcBoxMsg msg;
      procBuff.unpack<NeighboringProcBoxMsg>(msg);
      if (msg.msgType == DOMAIN_NEIGHBORING_MSG) {
        if (nbrDomRanks.find(msg.rankBox) == nbrDomRanks.end()) {
          nbrDomRanks.insert(msg.rankBox);
        }
      }
      else if (nbrRngRanks.find(msg.rankBox) == nbrRngRanks.end()) {
        nbrRngRanks.insert(msg.rankBox);
      }
    }
  }
  neighborDomainRanks.clear();
  neighborRangeRanks.clear();
  for (const OwnedBox3D &procBox : nbrDomRanks) { neighborDomainRanks.push_back(procBox); }
  for (const OwnedBox3D &procBox : nbrRngRanks) { neighborRangeRanks.push_back(procBox); }
}

} } }
