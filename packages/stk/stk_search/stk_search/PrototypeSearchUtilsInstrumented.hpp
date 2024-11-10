// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef PROTOTYPE_SEARCH_UTILS_INSTRUMENTED_H_
#define PROTOTYPE_SEARCH_UTILS_INSTRUMENTED_H_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <math.h>                                        // for floor
#include <stddef.h>                                      // for size_t
#include <algorithm>                                     // for max
#include <iostream>                                      // for operator<<, etc
#include <map>                                           // for map, etc
#include <set>                                           // for set
#include <stk_search/CommonSearchUtilsInstrumented.hpp>  // for SplitTimer
#include <stk_util/parallel/CommNeighbors.hpp>
#include <tuple>                                         // for get, tuple
#include <unordered_map>                                 // for hash
#include <utility>                                       // for pair
#include <vector>                                        // for vector
#include "mpi.h"                                         // for MPI_Comm, etc
#include "stk_search/BoundingBox.hpp"                    // for add_to_box, etc
#include "stk_search/Box.hpp"                            // for Box
#include "stk_search/CommonSearchUtil.hpp"
#include "stk_search/kdtree/KDTree.hpp"
#include "stk_search/kdtree/KDTree_BoundingBox.hpp"
#include "stk_util/parallel/CommSparse.hpp"              // for CommSparse
#include "stk_util/parallel/ParallelComm.hpp"            // for CommBuffer, etc
namespace stk { class CommBufferV; }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################




namespace stk {
namespace search {
namespace experimental {


struct GhostingSearchTimeBreakdown {
  double computeProcBoundingVolume;
  double findNeighbors_rasterize;
  double findNeighbors_ghostingRegionResidents;
  double findNeighbors_witnessAndComputeNeighbors;
  double findNeighbors_communicateNeighborObjectBBs;
  double intersectLocalRangeBoxesWithNeighborBVs;
  double communicateRangeBoxesToNeighbors;
  double writeRangeBoxesAndGhostIdentifiers;

  void reset() {
    computeProcBoundingVolume
        = findNeighbors_rasterize
        = findNeighbors_ghostingRegionResidents
        = findNeighbors_witnessAndComputeNeighbors
        = findNeighbors_communicateNeighborObjectBBs
        = intersectLocalRangeBoxesWithNeighborBVs
        = communicateRangeBoxesToNeighbors
        = writeRangeBoxesAndGhostIdentifiers = 0.0;
  }

  std::ostream &streamit(std::ostream &os) const {
    std::cout << "    computeProcBoundingVolume = " << computeProcBoundingVolume << std::endl
              << "     findNeighborsFromScratch = " << (findNeighbors_rasterize
                                                        + findNeighbors_ghostingRegionResidents
                                                        + findNeighbors_witnessAndComputeNeighbors
                                                        + findNeighbors_communicateNeighborObjectBBs) << std::endl
              << "                           rasterize = " << findNeighbors_rasterize << std::endl
              << "             ghostingRegionResidents = " << findNeighbors_ghostingRegionResidents << std::endl
              << "          witnessAndComputeNeighbors = " << findNeighbors_witnessAndComputeNeighbors << std::endl
              << "        communicateNeighborObjectBBs = " << findNeighbors_communicateNeighborObjectBBs << std::endl
              << "    isectLocalRngBsWithNbrBVs = " << intersectLocalRangeBoxesWithNeighborBVs << std::endl
              << "    commRangeBoxesToNeighbors = " << communicateRangeBoxesToNeighbors << std::endl
              << "  writeRngBoxesAndGhostIdents = " << writeRangeBoxesAndGhostIdentifiers << std::endl;
    return os;
  }
};


struct SuperTile3D {

  double xLen;
  double yLen;
  double zLen;

  int numXTiles;
  int numYTiles;
  int numZTiles;
};


struct TilingIndices {
  int xIdx, yIdx, zIdx;
};

typedef std::map<TilingIndices, std::vector<int> > ResidencyMapT;

enum NeighborMsgTypeEnum {
  DOMAIN_NEIGHBORING_MSG,
  RANGE_NEIGHBORING_MSG
};

struct NeighboringRankMsg {
  NeighborMsgTypeEnum msgType;
  int neighborRank;
};

struct TilingIndicesMsg {
  int xIdx, yIdx, zIdx;
  int pRank;
};


template <typename T>
struct GS_Types {

  typedef stk::search::Box<T> Box3D;

  struct OwnedBox3D
  {
    int   owner;
    Box3D box;

    typedef T value_type;
  };

  typedef std::map<TilingIndices, std::vector<OwnedBox3D> > ResidencyBoxMapT;

  struct NeighboringProcBoxMsg {
    NeighborMsgTypeEnum msgType;
    OwnedBox3D          rankBox;
  };

  struct BoxedTilingIndices : public TilingIndices {
    Box3D box;
  };

  struct BoxedTilingIndicesMsg : public TilingIndicesMsg {
    Box3D box;
  };
};


inline bool operator<(const TilingIndices &a, const TilingIndices &b) {
  return ((a.xIdx < b.xIdx)
          || ((a.xIdx == b.xIdx)
              && ((a.yIdx < b.yIdx)
                  || ((a.yIdx == b.yIdx) && (a.zIdx < b.zIdx)))));
}


inline bool operator<(const typename GS_Types<float>::OwnedBox3D &a, const typename GS_Types<float>::OwnedBox3D &b)
{
  if (a.owner < b.owner) {
    return true;
  }
  if (a.owner > b.owner) {
    return false;
  }
  const typename GS_Types<float>::Box3D &boxA = a.box;
  const typename GS_Types<float>::Box3D &boxB = b.box;

  return ((boxA.get_x_min() < boxB.get_x_min())
          || ((boxA.get_x_min() == boxB.get_x_min())
              && ((boxA.get_y_min() < boxB.get_y_min())
                  || ((boxA.get_y_min() == boxB.get_y_min())
                      && ((boxA.get_z_min() < boxB.get_z_min())
                          || ((boxA.get_z_min() == boxB.get_z_min())
                              && ((boxA.get_x_max() < boxB.get_x_max())
                                  || ((boxA.get_x_max() == boxB.get_x_max())
                                      && ((boxA.get_y_max() < boxB.get_y_max())
                                          || ((boxA.get_y_max() == boxB.get_y_max())
                                              && ((boxA.get_z_max() < boxB.get_z_max()))))))))))));
}


inline bool operator<(const typename GS_Types<double>::OwnedBox3D &a, const typename GS_Types<double>::OwnedBox3D &b)
{
  if (a.owner < b.owner) {
    return true;
  }
  if (a.owner > b.owner) {
    return false;
  }
  const typename GS_Types<double>::Box3D &boxA = a.box;
  const typename GS_Types<double>::Box3D &boxB = b.box;

  return ((boxA.get_x_min() < boxB.get_x_min())
          || ((boxA.get_x_min() == boxB.get_x_min())
              && ((boxA.get_y_min() < boxB.get_y_min())
                  || ((boxA.get_y_min() == boxB.get_y_min())
                      && ((boxA.get_z_min() < boxB.get_z_min())
                          || ((boxA.get_z_min() == boxB.get_z_min())
                              && ((boxA.get_x_max() < boxB.get_x_max())
                                  || ((boxA.get_x_max() == boxB.get_x_max())
                                      && ((boxA.get_y_max() < boxB.get_y_max())
                                          || ((boxA.get_y_max() == boxB.get_y_max())
                                              && ((boxA.get_z_max() < boxB.get_z_max()))))))))))));
}


void findGhostingRegionResidents(const SuperTile3D &tilingPattern, MPI_Comm comm, int mpiRank, int mpiSize,
                                 const std::vector<TilingIndices> &tilesOccupied,
                                 ResidencyMapT &residencyMapOfGhostingRegion);

template <typename T>
void findGhostingRegionResidentsDDEfficient(const SuperTile3D &tilingPattern, MPI_Comm comm, int mpiRank, int mpiSize,
                                         const std::vector<typename GS_Types<T>::BoxedTilingIndices> &tilesOccupied,
                                         typename GS_Types<T>::ResidencyBoxMapT &residencyMap)
{
  typedef typename GS_Types<T>::BoxedTilingIndices    BoxedTilingIndices;
  typedef typename GS_Types<T>::BoxedTilingIndicesMsg BoxedTilingIndicesMsg;
  typedef typename GS_Types<T>::OwnedBox3D            OwnedBox3D;

  // std::cout << "p_" << mpiRank << ": " << "fGRR -- tilesOccupied.size() returns " << tilesOccupied.size() << std::endl;

  residencyMap.clear();

  stk::CommSparse commsparse(comm);

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

    }
  }
}


inline int getGhostingRegionRank(const TilingIndices &tileId, const SuperTile3D &tilingPattern, int mpiSize) {
  int xIdx = tileId.xIdx % tilingPattern.numXTiles;
  int yIdx = tileId.yIdx % tilingPattern.numYTiles;
  int zIdx = tileId.zIdx % tilingPattern.numZTiles;

  xIdx = (xIdx < 0 ? xIdx + tilingPattern.numXTiles : xIdx);
  yIdx = (yIdx < 0 ? yIdx + tilingPattern.numYTiles : yIdx);
  zIdx = (zIdx < 0 ? zIdx + tilingPattern.numZTiles : zIdx);

  int proc = 0;
  proc += xIdx * tilingPattern.numYTiles * tilingPattern.numZTiles;
  proc += yIdx * tilingPattern.numZTiles;
  proc += zIdx;
  return proc;
}


template<typename BoxType>
BoxType computeGlobalBoundingBox(MPI_Comm mpiComm, const BoxType &box)
{
  double minCornerIn[3], maxCornerIn[3];
  minCornerIn[0] = box.get_x_min();
  minCornerIn[1] = box.get_y_min();
  minCornerIn[2] = box.get_z_min();
  maxCornerIn[0] = box.get_x_max();
  maxCornerIn[1] = box.get_y_max();
  maxCornerIn[2] = box.get_z_max();

  double minCornerOut[3], maxCornerOut[3];
  MPI_Allreduce(minCornerIn, minCornerOut, 3, MPI_DOUBLE, MPI_MIN, mpiComm);
  MPI_Allreduce(maxCornerIn, maxCornerOut, 3, MPI_DOUBLE, MPI_MAX, mpiComm);

  BoxType retval(minCornerOut[0], minCornerOut[1], minCornerOut[2],
                   maxCornerOut[0], maxCornerOut[1], maxCornerOut[2]);
  return retval;
}


void optimizeTiling3D(SuperTile3D &tile, int maxTiles);

template<typename BoxedType>
void rasterizeToTiling(const SuperTile3D &tilingPattern,
                       const stk::search::ObjectBoundingBox_T<BoxedType> &box,
                       std::vector<TilingIndices> &tilesOccupied)
{
  double xTileLen = tilingPattern.xLen;
  double yTileLen = tilingPattern.yLen;
  double zTileLen = tilingPattern.zLen;

  double boxXMin = box.GetBox().get_x_min();
  double boxYMin = box.GetBox().get_y_min();
  double boxZMin = box.GetBox().get_z_min();
  double boxXMax = box.GetBox().get_x_max();
  double boxYMax = box.GetBox().get_y_max();
  double boxZMax = box.GetBox().get_z_max();

  int firstXInd = floor(boxXMin / xTileLen);
  int firstYInd = floor(boxYMin / yTileLen);
  int firstZInd = floor(boxZMin / zTileLen);
  int lastXInd  = floor(boxXMax / xTileLen);
  int lastYInd  = floor(boxYMax / yTileLen);
  int lastZInd  = floor(boxZMax / zTileLen);

  tilesOccupied.clear();
  for (int i = firstXInd; i <= lastXInd; ++i) {
    for (int j = firstYInd; j <= lastYInd; ++j) {
      for (int k = firstZInd; k <= lastZInd; ++k) {
        tilesOccupied.push_back(TilingIndices{i, j, k});
      }
    }
  }
}

struct HashRegionId {
  inline std::size_t operator()(const std::tuple<int,int,int>& id) const
  {
    // Try to avoid collisions within 3D patches up to 1000 * cell_diameter.
    return std::hash<int>()((std::get<0>(id) << 20) ^ (std::get<1>(id) << 10) ^ std::get<2>(id));
  }
};

template<typename ObjType>
void rasterizeToTilingDDEfficient(const SuperTile3D &tilingPattern,
                               const std::vector<ObjType>& objs,
                               std::vector<typename GS_Types<typename ObjType::value_type>::BoxedTilingIndices> &tilesOccupied)
{
  typedef typename ObjType::value_type  valueType;
  typedef stk::search::Box<valueType>   Box;
  typedef std::tuple<int,int,int>       RegionId;

  //  Need to do experiments to see what works best as a function of the
  //  number and spatial spread (wrt ghosting region cell diameter).
#ifdef RASTERIZE_FOR_GHOSTING_DDE_USE_UNORDERED_MAP
  typedef std::unordered_map<RegionId, Box, HashRegionId> OccupancyMapT;
#else
  typedef std::map<RegionId, Box> OccupancyMapT;
#endif
  OccupancyMapT occupancyMap;

  double xTileLen = tilingPattern.xLen;
  double yTileLen = tilingPattern.yLen;
  double zTileLen = tilingPattern.zLen;

  for (const ObjType &obj : objs) {
    Box box;
    stk::search::add_to_box(box, obj);

    double boxXMin = box.get_x_min();
    double boxYMin = box.get_y_min();
    double boxZMin = box.get_z_min();
    double boxXMax = box.get_x_max();
    double boxYMax = box.get_y_max();
    double boxZMax = box.get_z_max();

    int firstXInd = floor(boxXMin / xTileLen);
    int firstYInd = floor(boxYMin / yTileLen);
    int firstZInd = floor(boxZMin / zTileLen);
    int lastXInd  = floor(boxXMax / xTileLen);
    int lastYInd  = floor(boxYMax / yTileLen);
    int lastZInd  = floor(boxZMax / zTileLen);

    for (int i = firstXInd; i <= lastXInd; ++i) {
      for (int j = firstYInd; j <= lastYInd; ++j) {
        for (int k = firstZInd; k <= lastZInd; ++k) {
          RegionId regionId(i,j,k);
          auto probe = occupancyMap.find(regionId);
          if (probe == occupancyMap.end()) {
            // TRS: Clang does not like emplace_hint,
            // use insert instead it should be overloaded.
            //occupancyMap.emplace_hint(probe, regionId, box);
            occupancyMap.insert(probe, std::pair<RegionId,Box>(regionId, box));
          }
          else {
            stk::search::add_to_box(probe->second, box);
          }
        }
      }
    }
  }

  // Now write the output vector.
  tilesOccupied.clear();
  tilesOccupied.resize(occupancyMap.size());
  int i = 0;
  for (auto mappedBox : occupancyMap) {
    const RegionId &tileId = mappedBox.first;
    Box &bbox = mappedBox.second;
    typename GS_Types<typename ObjType::value_type>::BoxedTilingIndices idxBox;
    idxBox.xIdx = std::get<0>(tileId);
    idxBox.yIdx = std::get<1>(tileId);
    idxBox.zIdx = std::get<2>(tileId);
    idxBox.box = bbox;
    tilesOccupied[i] = idxBox;
    ++i;
  }
}


void witnessAndComputeNeighborRanks(MPI_Comm mpiComm, int pSize,
                                    const ResidencyMapT &myDomainResidents,
                                    const ResidencyMapT &myRangeResidents,
                                    std::vector<int> &neighborDomainRanks, std::vector<int> &neighborRangeRanks);

template <typename T>
void witnessAndComputeNeighborRanksDDEfficient(MPI_Comm mpiComm, int pSize,
                                               const typename GS_Types<T>::ResidencyBoxMapT &myDomainResidents,
                                               const typename GS_Types<T>::ResidencyBoxMapT &myRangeResidents,
                                               std::vector<std::pair<typename GS_Types<T>::OwnedBox3D, typename GS_Types<T>::OwnedBox3D> > &neighborPairs)
{
  typedef typename GS_Types<T>::ResidencyBoxMapT ResidencyBoxMapT;
  typedef typename GS_Types<T>::OwnedBox3D       OwnedBox3D;

  int myRank = -1;
  MPI_Comm_rank(mpiComm, &myRank);

  typename ResidencyBoxMapT::const_iterator domResIter = myDomainResidents.begin();
  typename ResidencyBoxMapT::const_iterator domResEnd  = myDomainResidents.end();
  typename ResidencyBoxMapT::const_iterator rngResEnd  = myRangeResidents.end();

  //
  // Witness, and per-Ghosting-Region neighbor computation
  //
  for (; domResIter != domResEnd; ++domResIter) {
    const TilingIndices &tileId = domResIter->first;
    const std::vector<OwnedBox3D> &domainResidents = domResIter->second;
    typename ResidencyBoxMapT::const_iterator rngResProbe = myRangeResidents.find(tileId);
    if (rngResProbe != rngResEnd) {
      const std::vector<OwnedBox3D> &rangeResidents = rngResProbe->second;
      for (const OwnedBox3D &domainBox : domainResidents) {
        for (const OwnedBox3D &rangeBox : rangeResidents) {
          if ((domainBox.owner != rangeBox.owner) && intersects(domainBox.box, rangeBox.box)) {
            neighborPairs.emplace_back(domainBox, rangeBox);
          }
        }
      }
    }
  }

  // std::cout << "p_" << myRank << ": I witness " << neighborPairs.size() << " neighbor pairs" << std::endl;

}



template<typename DomainBox, typename RangeBox>
void communicateNeighborObjectBBs(
    MPI_Comm mpiComm,
    const std::vector<int> &neighborDomainRanks, const std::vector<int> &neighborRangeRanks,
    const stk::search::ObjectBoundingBox_T<DomainBox> domainObjBBox,
    const stk::search::ObjectBoundingBox_T<RangeBox> rangeObjBBox,
    std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > &boxA_proc_box_array,
    std::vector<stk::search::ObjectBoundingBox_T<RangeBox> > &boxB_proc_box_array)
{
  typedef stk::search::ObjectBoundingBox_T<DomainBox> DomObjBBoxT;

  stk::CommNeighbors commneighborsR2D(mpiComm, neighborRangeRanks, neighborDomainRanks);

  for (int nbrRangeProc : neighborRangeRanks) {
    stk::CommBufferV& procBuff = commneighborsR2D.send_buffer(nbrRangeProc);
    procBuff.pack<DomObjBBoxT>(domainObjBBox);
  }
  commneighborsR2D.communicate();
  for(int nbrDomProc : neighborDomainRanks) {
    stk::CommBufferV& procBuff = commneighborsR2D.recv_buffer(nbrDomProc  );
    DomObjBBoxT box;
    procBuff.unpack(box);
    boxA_proc_box_array[box.get_object_number()] = box;
  }

  typedef stk::search::ObjectBoundingBox_T<RangeBox> RngObjBBoxT;

  stk::CommNeighbors commneighborsD2R(mpiComm, neighborDomainRanks, neighborRangeRanks);

  for (int nbrDomainProc : neighborDomainRanks) {
    stk::CommBufferV& procBuff = commneighborsD2R.send_buffer(nbrDomainProc);
    procBuff.pack<RngObjBBoxT>(rangeObjBBox);
  }
  commneighborsD2R.communicate();
  for(int nbrRngProc : neighborRangeRanks) {
    stk::CommBufferV& procBuff = commneighborsD2R.recv_buffer(nbrRngProc);
    RngObjBBoxT box;
    procBuff.unpack(box);
    boxB_proc_box_array[box.get_object_number()] = box;
  }
}

template <typename T>
void communicateNeighborObjectBoundingBoxesDDEfficient(
    MPI_Comm comm, int mpiSize, const std::vector<std::pair<typename GS_Types<T>::OwnedBox3D, typename GS_Types<T>::OwnedBox3D> > &neighborPairs,
    std::vector<stk::search::ObjectBoundingBox_T<typename GS_Types<T>::Box3D> > &boxesA_proc_box_array,
    std::vector<stk::search::ObjectBoundingBox_T<typename GS_Types<T>::Box3D> > &boxesB_proc_box_array)
{
  typedef typename GS_Types<T>::NeighboringProcBoxMsg NeighboringProcBoxMsg;
  typedef typename GS_Types<T>::OwnedBox3D            OwnedBox3D;
  //
  // Communication and finish global neighbor computation.
  //

  stk::CommSparse commsparse(comm);
  for(int phase = 0; phase < 2; ++phase) {
    for (const std::pair<OwnedBox3D, OwnedBox3D> &neighborPair : neighborPairs) {
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

  typedef stk::search::ObjectBoundingBox_T<typename GS_Types<T>::Box3D> ObjBoxT;

  boxesA_proc_box_array.clear();
  for (const typename GS_Types<T>::OwnedBox3D &nbrBox : nbrDomRanks) {
    boxesA_proc_box_array.push_back(ObjBoxT(nbrBox.box, nbrBox.owner));
  }
  boxesB_proc_box_array.clear();
  for (const typename GS_Types<T>::OwnedBox3D &nbrBox : nbrRngRanks) {
    boxesB_proc_box_array.push_back(ObjBoxT(nbrBox.box, nbrBox.owner));
  }
}



template<typename FilterBoxType, typename BaseBoxType>
void filterProcBoundingBoxes(int skipRank, const FilterBoxType &filterBox,
                             std::vector<stk::search::ObjectBoundingBox_T<BaseBoxType> > &objectBoxes)
{
  std::vector<stk::search::ObjectBoundingBox_T<BaseBoxType> > collected;
  int numProcs = objectBoxes.size();
  for (int iproc = 0; iproc < numProcs; ++iproc) {
    if ((iproc != skipRank) && intersects(objectBoxes[iproc].GetBox(), filterBox)) {
      collected.push_back(objectBoxes[iproc]);
    }
  }

  objectBoxes = collected;
}


template<typename DomainBox, typename RangeBox>
void findGhostingNeighborsFromScratch(MPI_Comm comm, int mpiRank, int mpiSize,
                                      const stk::search::ObjectBoundingBox_T<DomainBox> boxA_proc,
                                      const stk::search::ObjectBoundingBox_T<RangeBox> boxB_proc,
                                      std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > &boxA_proc_box_array,
                                      std::vector<stk::search::ObjectBoundingBox_T<RangeBox> > &boxB_proc_box_array,
                                      experimental::GhostingSearchTimeBreakdown &timeBreakdown,
                                      stk::search::SplitTimer &stopwatch)
{
  DomainBox globalDomBox = computeGlobalBoundingBox(comm, boxA_proc.GetBox());
  timeBreakdown.computeProcBoundingVolume += stopwatch.split();

  SuperTile3D tiling = {globalDomBox.get_x_max() - globalDomBox.get_x_min(),
                        globalDomBox.get_y_max() - globalDomBox.get_y_min(),
                        globalDomBox.get_z_max() - globalDomBox.get_z_min(),
                        1, 1, 1};

  // Guess at good trade-off between AABBs overlapping more ghosting regions
  // versus having more neighbors in a region.
  const int cellDiamFactor = 2;
  int maxCells =  mpiSize / (cellDiamFactor * cellDiamFactor * cellDiamFactor);
  // optimizeTiling3D(tile, pSize);
  optimizeTiling3D(tiling, std::max(1, maxCells));

  // Ghosting Regions are identified by their TilingIndices.
  // We need find which regions our domain and our range occupies (3D raster).
  std::vector<TilingIndices> myDomainOccupancies;
  std::vector<TilingIndices> myRangeOccupancies;
  rasterizeToTiling(tiling, boxA_proc, myDomainOccupancies);
  rasterizeToTiling(tiling, boxB_proc, myRangeOccupancies);
  timeBreakdown.findNeighbors_rasterize += stopwatch.split();

  // In role of Decomp Region: Send Ghosting Regions corresponding myDomainOccupancies info.
  // In role of Decomp Region: Send Ghosting Regions corresponding myRangeOccupancies info.
  // In role of Ghosting Region: Receive Domain occupancy info.
  // In role of Ghosting Region: Receive Range occupancy info.
  ResidencyMapT myDomainResidents;
  ResidencyMapT myRangeResidents;
  findGhostingRegionResidents(tiling, comm, mpiRank, mpiSize, myDomainOccupancies, myDomainResidents);
  findGhostingRegionResidents(tiling, comm, mpiRank, mpiSize, myRangeOccupancies, myRangeResidents);
  timeBreakdown.findNeighbors_ghostingRegionResidents += stopwatch.split();

  // Now can witness all pairs of (domain_rank, range_rank) neighbor relationships that exist
  // on my Ghosting Region....
  // In role of Ghosting Region: send each range_rank the (domain_aabb, domain_rank) information
  // about its neighbor relationships I witness.
  // In role of Decomp Region: Receive my neighbor information.
  std::vector<int> domainRanksNeighboringMyRange, rangeRanksNeighboringMyDomain;
  witnessAndComputeNeighborRanks(comm, mpiSize, myDomainResidents, myRangeResidents,
                                 domainRanksNeighboringMyRange, rangeRanksNeighboringMyDomain);
  timeBreakdown.findNeighbors_witnessAndComputeNeighbors += stopwatch.split();

  communicateNeighborObjectBBs(comm, domainRanksNeighboringMyRange, rangeRanksNeighboringMyDomain,
                               boxA_proc, boxB_proc, boxA_proc_box_array, boxB_proc_box_array);

  filterProcBoundingBoxes(mpiRank, boxB_proc.GetBox(), boxA_proc_box_array);
  filterProcBoundingBoxes(mpiRank, boxA_proc.GetBox(), boxB_proc_box_array);

  timeBreakdown.findNeighbors_communicateNeighborObjectBBs += stopwatch.split();
}

template<typename DomainBox, typename RangeBox>
void findGhostingNeighborsFromScratchDDEfficient(MPI_Comm comm, int mpiRank, int mpiSize,
                                                 const stk::search::ObjectBoundingBox_T<DomainBox> boxA_proc,
                                                 const stk::search::ObjectBoundingBox_T<RangeBox> boxB_proc,
                                                 const std::vector<DomainBox> &localDomain,
                                                 const std::vector<RangeBox>  &localRange,
                                                 std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > &boxA_proc_box_array,
                                                 std::vector<stk::search::ObjectBoundingBox_T<RangeBox> > &boxB_proc_box_array,
                                                 experimental::GhostingSearchTimeBreakdown &timeBreakdown,
                                                 stk::search::SplitTimer &stopwatch)
{
  typedef typename DomainBox::value_type ValT;

  DomainBox globalDomBox = computeGlobalBoundingBox(comm, boxA_proc.GetBox());

  timeBreakdown.computeProcBoundingVolume += stopwatch.split();

  SuperTile3D tiling = {globalDomBox.get_x_max() - globalDomBox.get_x_min(),
                        globalDomBox.get_y_max() - globalDomBox.get_y_min(),
                        globalDomBox.get_z_max() - globalDomBox.get_z_min(),
                        1, 1, 1};

  // Guess at good trade-off between AABBs overlapping more ghosting regions
  // versus having more neighbors in a region.
  const int cellDiamFactor = 2;
  int maxCells =  mpiSize / (cellDiamFactor * cellDiamFactor * cellDiamFactor);
  // optimizeTiling3D(tile, pSize);
  optimizeTiling3D(tiling, std::max(1, maxCells));
  // if (pRank == 0) {
  //   std::cout << "SuperTile3D: dims = (" << tile.xLen << "," << tile.yLen << "," << tile.zLen << ")"
  //             << "  subtiles = ("  << tile.numXTiles << "," << tile.numYTiles << "," << tile.numZTiles << ")" << std::endl;
  // }

  // Ghosting Regions are identified by their TilingIndices.
  // We need find which regions our domain and our range occupies (3D raster).
  std::vector<typename GS_Types<ValT>::BoxedTilingIndices> myDomainOccupanciesDDE;
  std::vector<typename GS_Types<ValT>::BoxedTilingIndices> myRangeOccupanciesDDE;
  rasterizeToTilingDDEfficient(tiling, localDomain, myDomainOccupanciesDDE);
  rasterizeToTilingDDEfficient(tiling, localRange, myRangeOccupanciesDDE);

  // if ((pRank == 0) || (pRank == 63)) {
  //   std::cout << "p_" << pRank << ": myDomainOccupancies = {";
  //   for (auto bti : myDomainOccupanciesDDE) {
  //     std::cout << "{(" << bti.xIdx << "," << bti.yIdx << "," <<  bti.zIdx << ") " << bti.box << "} ";
  //   }
  //   std::cout <<"}" << std::endl;
  // }

  timeBreakdown.findNeighbors_rasterize += stopwatch.split();

  // In role of Decomp Region: Send Ghosting Regions corresponding myDomainOccupancies info.
  // In role of Decomp Region: Send Ghosting Regions corresponding myRangeOccupancies info.
  // In role of Ghosting Region: Receive Domain occupancy info.
  // In role of Ghosting Region: Receive Range occupancy info.
  typename GS_Types<ValT>::ResidencyBoxMapT myDomainResidentsDDE;
  typename GS_Types<ValT>::ResidencyBoxMapT myRangeResidentsDDE;
  findGhostingRegionResidentsDDEfficient<ValT>(tiling, comm, mpiRank, mpiSize, myDomainOccupanciesDDE, myDomainResidentsDDE);
  findGhostingRegionResidentsDDEfficient<ValT>(tiling, comm, mpiRank, mpiSize, myRangeOccupanciesDDE, myRangeResidentsDDE);

  timeBreakdown.findNeighbors_ghostingRegionResidents += stopwatch.split();

  // Now can witness all pairs of (domain_rank, range_rank) neighbor relationships that exist
  // on my Ghosting Region....
  // In role of Ghosting Region: send each range_rank the (domain_aabb, domain_rank) information
  // about its neighbor relationships I witness.
  // In role of Decomp Region: Receive my neighbor information.
  std::vector<std::pair<typename GS_Types<ValT>::OwnedBox3D, typename GS_Types<ValT>::OwnedBox3D> > neighborPairs;
  witnessAndComputeNeighborRanksDDEfficient<ValT>(comm, mpiSize, myDomainResidentsDDE, myRangeResidentsDDE,
                                                  neighborPairs);

  timeBreakdown.findNeighbors_witnessAndComputeNeighbors += stopwatch.split();

  communicateNeighborObjectBoundingBoxesDDEfficient<ValT>(comm, mpiSize, neighborPairs,
                                                          boxA_proc_box_array, boxB_proc_box_array);

  timeBreakdown.findNeighbors_communicateNeighborObjectBBs += stopwatch.split();
}

enum FindNeighborsAlgorithmChoice {
  ORIGINAL_FIND_NEIGHBORS,
  SCALABLE_FIND_NEIGHBORS,
  SCALABLE_AND_DISCONTINUOUS_DECOMP_EFFICIENT_FIND_NEIGHBORS
};

template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
class GhostingSearcher {

  typedef typename DomainObjType::value_type domainValueType;
  typedef stk::search::Box<domainValueType>  DomainBox;
  typedef typename RangeBoxType::value_type  rangeValueType;
  typedef stk::search::Box<rangeValueType>   RangeBox;
  typedef typename RangeIdentifier::ident_type GlobalIdType;
  typedef std::pair<RangeBoxType, GlobalIdType> BoxIdPair;

 public:
  GhostingSearcher(
      const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
      const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
      std::vector<RangeBoxType >& rangeBoxes,
      std::vector<RangeIdentifier>& rangeGhostIdentifiers, MPI_Comm comm,
      experimental::GhostingSearchTimeBreakdown &timeBreakdown);

  void searchFromScratch(
      const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
      const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
      std::vector<RangeBoxType >& rangeBoxes,
      std::vector<RangeIdentifier>& rangeGhostIdentifiers,
      experimental::GhostingSearchTimeBreakdown &timeBreakdown,
      FindNeighborsAlgorithmChoice findNeighborsAlgorithm);

  void updateAndSearch(
      const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
      const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
      std::vector<RangeBoxType >& rangeBoxes,
      std::vector<RangeIdentifier>& rangeGhostIdentifiers,
      experimental::GhostingSearchTimeBreakdown &timeBreakdown,
      FindNeighborsAlgorithmChoice findNeighborsAlgorithm);

 private:
  MPI_Comm mpiComm;
  int      pSize;
  int      pRank;

  void findNeighborsFromScratch(
      const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
      const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
      std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > &boxA_proc_box_array,
      std::vector<stk::search::ObjectBoundingBox_T<RangeBox> > &boxB_proc_box_array,
      experimental::GhostingSearchTimeBreakdown &timeBreakdown,
      SplitTimer &stopwatch);

  void findNeighborsFromScratchScalable(
      const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
      const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
      std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > &boxA_proc_box_array,
      std::vector<stk::search::ObjectBoundingBox_T<RangeBox> > &boxB_proc_box_array,
      experimental::GhostingSearchTimeBreakdown &timeBreakdown,
      SplitTimer &stopwatch);

  void findNeighborsFromScratchDDEfficient(
      const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
      const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
      std::vector<stk::search::ObjectBoundingBox_T<typename GS_Types<typename DomainObjType::value_type>::Box3D> > &boxesA_proc_box_array,
      std::vector<stk::search::ObjectBoundingBox_T<typename GS_Types<typename RangeBoxType::value_type>::Box3D> > &boxesB_proc_box_array,
      experimental::GhostingSearchTimeBreakdown &timeBreakdown,
      SplitTimer &stopwatch);

  void intersectLocalRangeBoxesWithNeighborBVs(
      std::vector<stk::search::ObjectBoundingBox_T<DomainBox> >& boxA_proc_box_array,
      const unsigned numBoxRange,
      const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
      std::vector<std::vector<BoxIdPair> >& send_list,
      std::vector<RangeBoxType>& rangeBoxes,
      experimental::GhostingSearchTimeBreakdown& timeBreakdown,
      SplitTimer& stopwatch);

  void intersectLocalRangeBoxesWithNeighborBVsDDEfficient(
      std::vector<stk::search::ObjectBoundingBox_T<typename GS_Types<typename DomainObjType::value_type>::Box3D> >& boxA_proc_box_array,
      const unsigned numBoxRange,
      const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
      std::vector<std::vector<BoxIdPair> >& send_list,
      std::vector<RangeBoxType>& rangeBoxes,
      experimental::GhostingSearchTimeBreakdown& timeBreakdown,
      SplitTimer& stopwatch);

  void commIntersectingRangeBoxesToNeighbors(
      std::vector<std::vector<BoxIdPair> > &send_list,
      std::vector<std::vector<BoxIdPair> > &recv_list,
      experimental::GhostingSearchTimeBreakdown& timeBreakdown,
      SplitTimer& stopwatch);

  void writeRangeBoxesAndGhostIdentifiers(
      const std::vector<std::vector<BoxIdPair> >& recv_list,
      SplitTimer stopwatch, std::vector<RangeIdentifier>& rangeGhostIdentifiers,
      std::vector<RangeBoxType>& rangeBoxes,
      experimental::GhostingSearchTimeBreakdown& timeBreakdown);
};


template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType, RangeBoxType>::GhostingSearcher(
      const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
      const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
      std::vector<RangeBoxType >& rangeBoxes,
      std::vector<RangeIdentifier>& rangeGhostIdentifiers, MPI_Comm comm,
      experimental::GhostingSearchTimeBreakdown &timeBreakdown)
      : mpiComm(comm), pSize(-1), pRank(-1)
{
  MPI_Comm_rank(mpiComm, &pRank);
  MPI_Comm_size(mpiComm, &pSize);
}


template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
void GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType, RangeBoxType>::findNeighborsFromScratch(
    const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
    const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
    std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > &boxA_proc_box_array,
    std::vector<stk::search::ObjectBoundingBox_T<RangeBox> > &boxB_proc_box_array,
    experimental::GhostingSearchTimeBreakdown &timeBreakdown,
    SplitTimer &stopwatch)
{
  //
  //  Compute the processor local bounding boxes for the box sets
  //  Store the boxes in unique entries in a global processor bounding box array.
  //
  stk::search::ObjectBoundingBox_T<DomainBox> boxA_proc;
  boxA_proc.set_object_number(pRank);
  const unsigned numBoxDomain = local_domain.size();
  for(unsigned iboxA = 0; iboxA < numBoxDomain; ++iboxA) {
    stk::search::add_to_box(boxA_proc.GetBox(), local_domain[iboxA].first);
  }

  stk::search::ObjectBoundingBox_T<RangeBox> boxB_proc;
  boxB_proc.set_object_number(pRank);
  const unsigned numBoxRange = local_range.size();
  for(unsigned iboxB = 0; iboxB < numBoxRange; ++iboxB) {
    stk::search::add_to_box(boxB_proc.GetBox(), local_range[iboxB].first);
  }

  timeBreakdown.computeProcBoundingVolume += stopwatch.split();

  //
  //  Do a global communication to communicate all processor boxA bounding boxes
  //  to all processors in the group
  //
  stk::search::all_gather_helper(boxA_proc, boxA_proc_box_array, mpiComm);
  for(int iproc = 0; iproc < pSize; ++iproc) {
    boxA_proc_box_array[iproc].set_object_number(iproc);
  }
  stk::search::all_gather_helper(boxB_proc, boxB_proc_box_array, mpiComm);
  for(int iproc = 0; iproc < pSize; ++iproc) {
    boxB_proc_box_array[iproc].set_object_number(iproc);
  }

  filterProcBoundingBoxes(pRank, boxB_proc.GetBox(), boxA_proc_box_array);
  filterProcBoundingBoxes(pRank, boxA_proc.GetBox(), boxB_proc_box_array);

  timeBreakdown.findNeighbors_communicateNeighborObjectBBs += stopwatch.split();

}


template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
void GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType, RangeBoxType>::findNeighborsFromScratchScalable(
    const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
    const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
    std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > &boxA_proc_box_array,
    std::vector<stk::search::ObjectBoundingBox_T<RangeBox> > &boxB_proc_box_array,
    experimental::GhostingSearchTimeBreakdown &timeBreakdown,
    SplitTimer &stopwatch)
{
  //
  //  Compute the processor local bounding boxes for the box sets
  //  Store the boxes in unique entries in a global processor bounding box array.
  //
  stk::search::ObjectBoundingBox_T<DomainBox> boxA_proc;
  boxA_proc.set_object_number(pRank);
  const unsigned numBoxDomain = local_domain.size();
  for(unsigned iboxA = 0; iboxA < numBoxDomain; ++iboxA) {
    stk::search::add_to_box(boxA_proc.GetBox(), local_domain[iboxA].first);
  }
  boxA_proc_box_array[pRank] = boxA_proc;

  stk::search::ObjectBoundingBox_T<RangeBox> boxB_proc;
  boxB_proc.set_object_number(pRank);
  const unsigned numBoxRange = local_range.size();
  for(unsigned iboxB = 0; iboxB < numBoxRange; ++iboxB) {
    stk::search::add_to_box(boxB_proc.GetBox(), local_range[iboxB].first);
  }
  boxB_proc_box_array[pRank] = boxB_proc;

  findGhostingNeighborsFromScratch(mpiComm, pRank, pSize, boxA_proc, boxB_proc,
                                   boxA_proc_box_array, boxB_proc_box_array,
                                   timeBreakdown, stopwatch);
}


template<typename BoxableObjType, typename BoxType, typename IdentType>
void ComputeBoxVector(const std::vector<std::pair<BoxableObjType, IdentType> >& input,
                     std::vector<BoxType> &output)
{
  const unsigned numObj = input.size();
  output.resize(input.size());
  for (unsigned i = 0; i < numObj; ++i) {
    BoxType box;
    stk::search::add_to_box(box, input[i].first);
    output[i] = box;
  }
}


template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
void GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType, RangeBoxType>::findNeighborsFromScratchDDEfficient(
    const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
    const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
    std::vector<stk::search::ObjectBoundingBox_T<typename GS_Types<typename DomainObjType::value_type>::Box3D> > &boxesA_proc_box_array,
    std::vector<stk::search::ObjectBoundingBox_T<typename GS_Types<typename RangeBoxType::value_type>::Box3D> > &boxesB_proc_box_array,
    experimental::GhostingSearchTimeBreakdown &timeBreakdown,
    SplitTimer &stopwatch)
{
  std::vector<DomainBox> localDomain;
  std::vector<RangeBox>  localRange;
  ComputeBoxVector(local_domain, localDomain);
  ComputeBoxVector(local_range, localRange);

  //
  //  Compute the processor local bounding boxes for the box sets
  //  Store the boxes in unique entries in a global processor bounding box array.
  //
  stk::search::ObjectBoundingBox_T<DomainBox> boxA_proc;
  boxA_proc.set_object_number(pRank);
  const unsigned numBoxDomain = localDomain.size();
  for(unsigned iboxA = 0; iboxA < numBoxDomain; ++iboxA) {
    stk::search::add_to_box(boxA_proc.GetBox(), localDomain[iboxA]);
  }

  stk::search::ObjectBoundingBox_T<RangeBox> boxB_proc;
  boxB_proc.set_object_number(pRank);
  const unsigned numBoxRange = localRange.size();
  for(unsigned iboxB = 0; iboxB < numBoxRange; ++iboxB) {
    stk::search::add_to_box(boxB_proc.GetBox(), localRange[iboxB]);
  }

  findGhostingNeighborsFromScratchDDEfficient(mpiComm, pRank, pSize, boxA_proc, boxB_proc,
                                              localDomain, localRange,
                                              boxesA_proc_box_array, boxesB_proc_box_array,
                                              timeBreakdown, stopwatch);

}



template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
void GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType,
  RangeBoxType>::intersectLocalRangeBoxesWithNeighborBVs(
    std::vector<stk::search::ObjectBoundingBox_T<DomainBox> >& boxA_proc_box_array,
    const unsigned numBoxRange,
    const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
    std::vector<std::vector<BoxIdPair> >& send_list,
    std::vector<RangeBoxType>& rangeBoxes,
    experimental::GhostingSearchTimeBreakdown& timeBreakdown,
    SplitTimer& stopwatch)
{
  //  Create a hierarchy of boxA processor bounding boxes.
  //  This hierarchy will be used to search for overlaps between processors and
  //  objects.
  stk::search::ProximitySearchTree_T<DomainBox> boxA_box_hierarchy(
      boxA_proc_box_array);

  //  Determine what to ghost.  If a boxB box from this processor overlaps another processor's
  //  processor all-boxA box, then we need to ghost the boxB's data to that other processor.
  //  (The stricter criteria used by the full BoxA_BoxB_Ghost function would make sure that
  //  the individual boxB box overlaps some individual boxA box from the other processor.)
  std::vector<int> proc_list(pSize);
  for (unsigned int iboxB = 0; iboxB < numBoxRange; ++iboxB) {
    boxA_box_hierarchy.SearchForOverlap(local_range[iboxB].first, proc_list);
    for (unsigned i = 0; i < proc_list.size(); ++i) {
      int overlapping_proc = proc_list[i];
      GlobalIdType id = local_range[iboxB].second.id();
      send_list[overlapping_proc].push_back(BoxIdPair(rangeBoxes[iboxB], id));
    }
  }
  timeBreakdown.intersectLocalRangeBoxesWithNeighborBVs += stopwatch.split();
}


template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
void GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType,
  RangeBoxType>::intersectLocalRangeBoxesWithNeighborBVsDDEfficient(
    std::vector<stk::search::ObjectBoundingBox_T<typename GS_Types<typename DomainObjType::value_type>::Box3D> >& boxA_proc_box_array,
    const unsigned numBoxRange,
    const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
    std::vector<std::vector<BoxIdPair> >& send_list,
    std::vector<RangeBoxType>& rangeBoxes,
    experimental::GhostingSearchTimeBreakdown& timeBreakdown,
    SplitTimer& stopwatch)
{
  //  Create a hierarchy of boxA processor bounding boxes.
  //  This hierarchy will be used to search for overlaps between processors and
  //  objects.
  stk::search::ProximitySearchTree_T<typename GS_Types<typename DomainObjType::value_type>::Box3D> boxA_box_hierarchy(
      boxA_proc_box_array);

  //  Determine what to ghost.  If a boxB box from this processor overlaps another processor's
  //  processor all-boxA box, then we need to ghost the boxB's data to that other processor.
  //  (The stricter criteria used by the full BoxA_BoxB_Ghost function would make sure that
  //  the individual boxB box overlaps some individual boxA box from the other processor.)
  std::vector<int> proc_list(pSize);
  for (unsigned int iboxB = 0; iboxB < numBoxRange; ++iboxB) {
    std::set<std::pair<int,int>> found;
    boxA_box_hierarchy.SearchForOverlap(local_range[iboxB].first, proc_list);
    for (unsigned i = 0; i < proc_list.size(); ++i) {
      int overlapping_proc = proc_list[i];
      GlobalIdType id = local_range[iboxB].second.id();
      std::pair<int,int> isectItem(overlapping_proc,id);
      auto probe = found.find(isectItem);
      if (probe == found.end()) {
        send_list[overlapping_proc].push_back(BoxIdPair(rangeBoxes[iboxB], id));
        found.insert(probe, isectItem);
      }
    }
  }
  timeBreakdown.intersectLocalRangeBoxesWithNeighborBVs += stopwatch.split();
}


template<typename DomainIdentifier, typename RangeIdentifier,
    typename DomainObjType, typename RangeBoxType>
void GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType,
    RangeBoxType>::commIntersectingRangeBoxesToNeighbors(
    std::vector<std::vector<BoxIdPair> > &send_list,
    std::vector<std::vector<BoxIdPair> > &recv_list,
    experimental::GhostingSearchTimeBreakdown& timeBreakdown,
    SplitTimer& stopwatch)
{
  // Deep inside the implementation, there is an MPI_Allreduce, but at 8K ranks,
  // this is outweighed by each of the current versions of findNeighbors(..).
  stk::parallel_data_exchange_t(send_list, recv_list, mpiComm);
  timeBreakdown.communicateRangeBoxesToNeighbors += stopwatch.split();
}


template<typename DomainIdentifier, typename RangeIdentifier,
    typename DomainObjType, typename RangeBoxType>
void GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType,
    RangeBoxType>::writeRangeBoxesAndGhostIdentifiers(
    const std::vector<std::vector<BoxIdPair> >& recv_list, SplitTimer stopwatch,
    std::vector<RangeIdentifier>& rangeGhostIdentifiers,
    std::vector<RangeBoxType>& rangeBoxes,
    experimental::GhostingSearchTimeBreakdown& timeBreakdown) {
  rangeGhostIdentifiers.clear();
  for (size_t i = 0; i < recv_list.size(); ++i) {
    for (size_t j = 0; j < recv_list[i].size(); ++j) {
      const BoxIdPair& recvd_boxIdPair = recv_list[i][j];
      rangeBoxes.push_back(recvd_boxIdPair.first);
      rangeGhostIdentifiers.push_back(
          RangeIdentifier(recvd_boxIdPair.second, i));
    }
  }
  timeBreakdown.writeRangeBoxesAndGhostIdentifiers += stopwatch.split();
}


template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
void
GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType, RangeBoxType>::searchFromScratch(
    const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
    const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
    std::vector<RangeBoxType >& rangeBoxes,
    std::vector<RangeIdentifier>& rangeGhostIdentifiers,
    experimental::GhostingSearchTimeBreakdown &timeBreakdown,
    FindNeighborsAlgorithmChoice findNeighborsAlgorithm)
{
  typedef typename DomainObjType::value_type ValT;
  // Pretend-o-type for TDD.

  const unsigned numBoxRange  = local_range.size();

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (size_t i = 0; i < numBoxRange; ++i) {
    rangeBoxes[i] = local_range[i].first;
  }
  if(pSize == 0) {
    return;
  }

  SplitTimer stopwatch;

  std::vector<std::vector<BoxIdPair> > recv_list(pSize);
  std::vector<std::vector<BoxIdPair> > send_list(pSize);

  if (findNeighborsAlgorithm == SCALABLE_AND_DISCONTINUOUS_DECOMP_EFFICIENT_FIND_NEIGHBORS) {
    std::vector<stk::search::ObjectBoundingBox_T<typename GS_Types<ValT>::Box3D> > boxA_proc_box_array(pSize);
    std::vector<stk::search::ObjectBoundingBox_T<typename GS_Types<ValT>::Box3D> > boxB_proc_box_array(pSize);
    findNeighborsFromScratchDDEfficient(local_domain, local_range,
                                        boxA_proc_box_array, boxB_proc_box_array,
                                        timeBreakdown, stopwatch);
    intersectLocalRangeBoxesWithNeighborBVsDDEfficient(boxA_proc_box_array, numBoxRange,
                                                    local_range, send_list, rangeBoxes,
                                                    timeBreakdown, stopwatch);
  }
  else {
    std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > boxA_proc_box_array(pSize);
    std::vector<stk::search::ObjectBoundingBox_T<RangeBox> > boxB_proc_box_array(pSize);
    if (findNeighborsAlgorithm == SCALABLE_FIND_NEIGHBORS){
      findNeighborsFromScratchScalable(local_domain, local_range, boxA_proc_box_array, boxB_proc_box_array,
                                       timeBreakdown, stopwatch);
    }
    else { // ORIGINAL_FIND_NEIGHBORS
      findNeighborsFromScratch(local_domain, local_range, boxA_proc_box_array, boxB_proc_box_array,
                               timeBreakdown, stopwatch);
    }
    intersectLocalRangeBoxesWithNeighborBVs(boxA_proc_box_array, numBoxRange,
                                            local_range, send_list, rangeBoxes,
                                            timeBreakdown, stopwatch);
  }

  commIntersectingRangeBoxesToNeighbors(send_list, recv_list, timeBreakdown,
                                        stopwatch);
  writeRangeBoxesAndGhostIdentifiers(recv_list, stopwatch,
                                     rangeGhostIdentifiers, rangeBoxes,
                                     timeBreakdown);
}


template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
void
GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType, RangeBoxType>::updateAndSearch(
    const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
    const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
    std::vector<RangeBoxType >& rangeBoxes,
    std::vector<RangeIdentifier>& rangeGhostIdentifiers,
    experimental::GhostingSearchTimeBreakdown &timeBreakdown,
    FindNeighborsAlgorithmChoice findNeighborsAlg)
{
  // Pretend-o-type!
  searchFromScratch(local_domain, local_range, rangeBoxes, rangeGhostIdentifiers, timeBreakdown,
                    findNeighborsAlg);
}


}}}


template <typename T>
std::ostream &operator<<(std::ostream &os, const typename stk::search::experimental::GS_Types<T>::BoxedTilingIndices &bti) {
  os << "{(" << bti.xIdx << " " << bti.yIdx << " " << bti.zIdx << ") " << bti.box << "}";
  return os;
}

#endif // PROTOTYPE_SEARCH_UTILS_INSTRUMENTED_H_
