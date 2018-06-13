
// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stk_search/PrototypeSearchUtilsInstrumented.hpp>
#include <stk_util/parallel/CommSparse.hpp>    // for CommSparse
#include <stk_util/parallel/ParallelComm.hpp>  // for CommBuffer
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################




namespace stk {
namespace search {
namespace experimental {


// Greedy brute-force algorithm, just the first thing that came to mind.  Don't know
// if the output is actually optimal, which needs to be better defined anyway.
void optimizeTiling3D(SuperTile3D &tile, int maxTiles)
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
            neighborPairs.emplace_back(domainRank, rangeRank);
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


} } }
