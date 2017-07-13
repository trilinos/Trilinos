#ifndef PROTOTYPE_SEARCH_UTILS_INSTRUMENTED_H_
#define PROTOTYPE_SEARCH_UTILS_INSTRUMENTED_H_

#include <tuple>
#include <unordered_set>

#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_size, etc
#include <stk_util/parallel/CommNeighbors.hpp>
#include <stk_util/parallel/MPI.hpp>

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

  std::ostream &streamit(std::ostream &os) {
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


//
// For robustness wrt discontinuous decompositions, we'll must be able to associate a
// decomp region with multiple bounding volumes (AABBs for now).
//

typedef stk::search::Box<double> Box3D;

struct OwnedBox3D
{
  int   owner;
  Box3D box;

  typedef double value_type;
};

typedef std::map<TilingIndices, std::vector<OwnedBox3D> > ResidencyBoxMapT;

struct BoxedTilingIndices : public TilingIndices {
  Box3D box;
};


inline bool operator<(const TilingIndices &a, const TilingIndices &b) {
  return ((a.xIdx < b.xIdx)
          || ((a.xIdx == b.xIdx)
              && ((a.yIdx < b.yIdx)
                  || ((a.yIdx == b.yIdx) && (a.zIdx < b.zIdx)))));
}

inline bool operator<(const OwnedBox3D &a, const OwnedBox3D &b)
{
  if (a.owner < b.owner) {
    return true;
  }
  if (a.owner > b.owner) {
    return false;
  }
  const Box3D &boxA = a.box;
  const Box3D &boxB = b.box;

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

void findGhostingRegionResidentsDDRobust(const SuperTile3D &tilingPattern, MPI_Comm comm, int mpiRank, int mpiSize,
                                         const std::vector<BoxedTilingIndices> &tilesOccupied,
                                         ResidencyBoxMapT &residencyMapOfGhostingRegion);

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

void optimizeTile3D(SuperTile3D &tile, int maxTiles);

template<typename BoxedType>
void rasterizeToTiling(const SuperTile3D &tilingPattern, int decompRank,
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


template<typename ObjType, typename IdentType>
void rasterizeToTilingDDRobust(const SuperTile3D &tilingPattern, int decompRank,
                               const std::vector<std::pair<ObjType, IdentType> >& objs,
                               std::vector<BoxedTilingIndices> &tilesOccupied)
{
  typedef typename ObjType::value_type  valueType;
  typedef stk::search::Box<valueType>   Box;
  typedef std::tuple<int,int,int>       RegionId;

  std::map<RegionId, Box> occupancyMap;

  double xTileLen = tilingPattern.xLen;
  double yTileLen = tilingPattern.yLen;
  double zTileLen = tilingPattern.zLen;

  for (const std::pair<ObjType, IdentType> &obj : objs) {
    Box box;
    stk::search::add_to_box(box, obj.first);

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

    tilesOccupied.clear();
    for (int i = firstXInd; i <= lastXInd; ++i) {
      for (int j = firstYInd; j <= lastYInd; ++j) {
        for (int k = firstZInd; k <= lastZInd; ++k) {
          RegionId regionId(i,j,k);
          if (occupancyMap.find(regionId) == occupancyMap.end()) {
            occupancyMap[regionId] = box;
          }
          else {
            stk::search::add_to_box(occupancyMap[regionId], box);
          }
        }
      }
    }
  }
  // Now need to write the output vector.
  for (auto mappedBox : occupancyMap) {
    const RegionId &tileId = mappedBox.first;
    const Box &bbox        = mappedBox.second;
    BoxedTilingIndices idxBox;
    idxBox.xIdx = std::get<0>(tileId);
    idxBox.yIdx = std::get<1>(tileId);
    idxBox.zIdx = std::get<2>(tileId);
    idxBox.box = bbox;
    tilesOccupied.push_back(idxBox);
  }
}

void witnessAndComputeNeighborRanks(MPI_Comm mpiComm, int pSize,
                                    const ResidencyMapT &myDomainResidents,
                                    const ResidencyMapT &myRangeResidents,
                                    std::vector<int> &neighborDomainRanks, std::vector<int> &neighborRangeRanks);

void witnessAndComputeNeighborRanksDDRobust(MPI_Comm mpiComm, int pSize,
                                            const ResidencyBoxMapT &myDomainResidents,
                                            const ResidencyBoxMapT &myRangeResidents,
                                            std::vector<OwnedBox3D> &neighborDomainRanks,
                                            std::vector<OwnedBox3D> &neighborRangeRanks);


enum FindNeighborsAlgorithmChoice {
  ORIGINAL_FIND_NEIGHBORS,
  SCALABLE_FIND_NEIGHBORS,
  SCALABLE_AND_DISCONTINUOUS_DECOMP_ROBUST_FIND_NEIGHBORS
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

  DomainBox computeGlobalBoundingBox(const DomainBox &box);

  void findNeighborsFromScratch(
      const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
      const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
      std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > &boxA_proc_box_array,
      experimental::GhostingSearchTimeBreakdown &timeBreakdown,
      SplitTimer &stopwatch);

  void findNeighborsFromScratchScalable(
      const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
      const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
      std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > &boxA_proc_box_array,
      experimental::GhostingSearchTimeBreakdown &timeBreakdown,
      SplitTimer &stopwatch);

  void findNeighborsFromScratchDDRobust(
      const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
      const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
      std::vector<stk::search::ObjectBoundingBox_T<Box3D> > &boxesA_proc_box_array,
      experimental::GhostingSearchTimeBreakdown &timeBreakdown,
      SplitTimer &stopwatch);

  void communicateDomainNeighborObjectBoundingBoxes(
      const std::vector<int> &neighborDomainRanks, const std::vector<int> &neighborRangeRanks,
      const stk::search::ObjectBoundingBox_T<DomainBox> &domainObjBBox,
      std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > &boxA_proc_box_array);

  void communicateDomainNeighborObjectBoundingBoxesDDRobust(
      const std::vector<OwnedBox3D> &neighborDomainRanks, const std::vector<OwnedBox3D> &neighborRangeRanks,
      const std::vector<BoxedTilingIndices> &domainOccupancies,
      std::vector<stk::search::ObjectBoundingBox_T<Box3D> > &boxesA_proc_box_array);

  void intersectLocalRangeBoxesWithNeighborBVs(
      std::vector<stk::search::ObjectBoundingBox_T<DomainBox> >& boxA_proc_box_array,
      const unsigned numBoxRange,
      const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
      std::vector<std::vector<BoxIdPair> >& send_list,
      std::vector<RangeBoxType>& rangeBoxes,
      experimental::GhostingSearchTimeBreakdown& timeBreakdown,
      SplitTimer& stopwatch);

  void intersectLocalRangeBoxesWithNeighborBVsDDRobust(
      std::vector<stk::search::ObjectBoundingBox_T<Box3D> >& boxA_proc_box_array,
      const unsigned numBoxRange,
      const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
      std::vector<std::vector<BoxIdPair> >& send_list,
      std::vector<RangeBoxType>& rangeBoxes,
      experimental::GhostingSearchTimeBreakdown& timeBreakdown,
      SplitTimer& stopwatch);

  void commIntersectingRangeBoxesToNeighbors(
      std::vector<std::vector<BoxIdPair> > send_list,
      std::vector<std::vector<BoxIdPair> >& recv_list,
      experimental::GhostingSearchTimeBreakdown& timeBreakdown,
      SplitTimer& stopwatch);

  void writeRangeBoxesAndGhostIdentifiers(
      const std::vector<std::vector<BoxIdPair> >& recv_list,
      SplitTimer stopwatch, std::vector<RangeIdentifier>& rangeGhostIdentifiers,
      std::vector<RangeBoxType>& rangeBoxes,
      experimental::GhostingSearchTimeBreakdown& timeBreakdown);

  void filterBoxAProcArrayByLocalRange(
      const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
      std::vector<stk::search::ObjectBoundingBox_T<DomainBox> >& boxA_proc_box_array);
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
typename GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType, RangeBoxType>::DomainBox
GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType, RangeBoxType>::computeGlobalBoundingBox(const DomainBox &box)
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

  DomainBox retval(minCornerOut[0], minCornerOut[1], minCornerOut[2],
                   maxCornerOut[0], maxCornerOut[1], maxCornerOut[2]);

  return retval;
}

template<typename DomainIdentifier, typename RangeIdentifier,
    typename DomainObjType, typename RangeBoxType>
void GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType,
    RangeBoxType>::filterBoxAProcArrayByLocalRange(
    const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
    std::vector<stk::search::ObjectBoundingBox_T<DomainBox> >& boxA_proc_box_array)
{
  stk::search::Box<typename RangeBoxType::value_type> localRangeAABB;
  int numRangeBoxes = local_range.size();
  for (int i = 0; i < numRangeBoxes; ++i) {
    stk::search::add_to_box(localRangeAABB, local_range[i].first);
  }

  std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > neighbors;
  for (int iproc = 0; iproc < pSize; ++iproc) {
    if ((iproc != pRank) && intersects(boxA_proc_box_array[iproc].GetBox(), localRangeAABB)) {
      neighbors.push_back(boxA_proc_box_array[iproc]);
    }
  }

  boxA_proc_box_array = neighbors;
}


template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
void GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType, RangeBoxType>::findNeighborsFromScratch(
    const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
    const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
    std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > &boxA_proc_box_array,
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

  timeBreakdown.computeProcBoundingVolume += stopwatch.split();

  //
  //  Do a global communication to communicate all processor boxA bounding boxes
  //  to all processors in the group
  //
  instrumented::GlobalBoxCombine(boxA_proc_box_array, mpiComm);
  for(int iproc = 0; iproc < pSize; ++iproc) {
    boxA_proc_box_array[iproc].set_object_number(iproc);
  }
  filterBoxAProcArrayByLocalRange(local_range, boxA_proc_box_array);
  timeBreakdown.findNeighbors_communicateNeighborObjectBBs += stopwatch.split();

}


template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
void GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType, RangeBoxType>::communicateDomainNeighborObjectBoundingBoxes(
    const std::vector<int> &neighborDomainRanks, const std::vector<int> &neighborRangeRanks,
    const stk::search::ObjectBoundingBox_T<DomainBox> &domainObjBBox,
    std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > &boxA_proc_box_array)
{
  typedef stk::search::ObjectBoundingBox_T<DomainBox> ObjBBoxT;

  stk::CommNeighbors commneighbors(mpiComm, neighborRangeRanks, neighborDomainRanks);

  for (int nbrRangeProc : neighborRangeRanks) {
    stk::CommBufferV& procBuff = commneighbors.send_buffer(nbrRangeProc);
    procBuff.pack<ObjBBoxT>(domainObjBBox);
  }
  commneighbors.communicate();
  for(int nbrDomProc : neighborDomainRanks) {
    stk::CommBufferV& procBuff = commneighbors.recv_buffer(nbrDomProc  );
    ObjBBoxT box;
    procBuff.unpack(box);
    boxA_proc_box_array[box.get_object_number()] = box;
  }
}

template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
void GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType, RangeBoxType>::communicateDomainNeighborObjectBoundingBoxesDDRobust(
    const std::vector<OwnedBox3D> &neighborDomainRanks, const std::vector<OwnedBox3D> &neighborRangeRanks,
    const std::vector<BoxedTilingIndices> &domainOccupancies,
    std::vector<stk::search::ObjectBoundingBox_T<Box3D> > &boxesA_proc_box_array)
{
  typedef stk::search::ObjectBoundingBox_T<Box3D> ObjBoxT;

  for(auto nbrBox : neighborDomainRanks) {
    boxesA_proc_box_array.push_back(ObjBoxT(nbrBox.box, nbrBox.owner));
  }
}


template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
void GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType, RangeBoxType>::findNeighborsFromScratchScalable(
    const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
    const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
    std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > &boxA_proc_box_array,
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

  DomainBox globalDomBox = computeGlobalBoundingBox(boxA_proc.GetBox());

  timeBreakdown.computeProcBoundingVolume += stopwatch.split();

  SuperTile3D tile = {globalDomBox.get_x_max() - globalDomBox.get_x_min(),
                      globalDomBox.get_y_max() - globalDomBox.get_y_min(),
                      globalDomBox.get_z_max() - globalDomBox.get_z_min(),
                      1, 1, 1};

  // optimizeTile3D(tile, pSize);
  optimizeTile3D(tile, std::max(1, pSize/8));

  // Ghosting Regions are identified by their TilingIndices.
  // We need find which regions our domain and our range occupies (3D raster).
  std::vector<TilingIndices> myDomainOccupancies;
  std::vector<TilingIndices> myRangeOccupancies;
  rasterizeToTiling(tile, pRank, boxA_proc, myDomainOccupancies);
  rasterizeToTiling(tile, pRank, boxB_proc, myRangeOccupancies);
  timeBreakdown.findNeighbors_rasterize += stopwatch.split();

  // In role of Decomp Region: Send Ghosting Regions corresponding myDomainOccupancies info.
  // In role of Decomp Region: Send Ghosting Regions corresponding myRangeOccupancies info.
  // In role of Ghosting Region: Receive Domain occupancy info.
  // In role of Ghosting Region: Receive Range occupancy info.
  ResidencyMapT myDomainResidents;
  ResidencyMapT myRangeResidents;
  findGhostingRegionResidents(tile, mpiComm, pRank, pSize, myDomainOccupancies, myDomainResidents);
  findGhostingRegionResidents(tile, mpiComm, pRank, pSize, myRangeOccupancies, myRangeResidents);
  timeBreakdown.findNeighbors_ghostingRegionResidents += stopwatch.split();

  // Now can witness all pairs of (domain_rank, range_rank) neighbor relationships that exist
  // on my Ghosting Region....
  // In role of Ghosting Region: send each range_rank the (domain_aabb, domain_rank) information
  // about its neighbor relationships I witness.
  // In role of Decomp Region: Receive my neighbor information.
  std::vector<int> domainRanksNeighboringMyRange, rangeRanksNeighboringMyDomain;
  witnessAndComputeNeighborRanks(mpiComm, pSize, myDomainResidents, myRangeResidents,
                                 domainRanksNeighboringMyRange, rangeRanksNeighboringMyDomain);
  timeBreakdown.findNeighbors_witnessAndComputeNeighbors += stopwatch.split();

  communicateDomainNeighborObjectBoundingBoxes(domainRanksNeighboringMyRange, rangeRanksNeighboringMyDomain,
                                               boxA_proc, boxA_proc_box_array);
  filterBoxAProcArrayByLocalRange(local_range, boxA_proc_box_array);
  timeBreakdown.findNeighbors_communicateNeighborObjectBBs += stopwatch.split();

}



template<typename DomainIdentifier, typename RangeIdentifier, typename DomainObjType, typename RangeBoxType>
void GhostingSearcher<DomainIdentifier, RangeIdentifier, DomainObjType, RangeBoxType>::findNeighborsFromScratchDDRobust(
    const std::vector<std::pair<DomainObjType, DomainIdentifier> >& local_domain,
    const std::vector<std::pair<RangeBoxType, RangeIdentifier> >& local_range,
    std::vector<stk::search::ObjectBoundingBox_T<Box3D> > &boxesA_proc_box_array,
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

  DomainBox globalDomBox = computeGlobalBoundingBox(boxA_proc.GetBox());

  timeBreakdown.computeProcBoundingVolume += stopwatch.split();

  SuperTile3D tile = {globalDomBox.get_x_max() - globalDomBox.get_x_min(),
                      globalDomBox.get_y_max() - globalDomBox.get_y_min(),
                      globalDomBox.get_z_max() - globalDomBox.get_z_min(),
                      1, 1, 1};

  // optimizeTile3D(tile, pSize);
  optimizeTile3D(tile, std::max(1, pSize/8));
  // if (pRank == 0) {
  //   std::cout << "SuperTile3D: dims = (" << tile.xLen << "," << tile.yLen << "," << tile.zLen << ")"
  //             << "  subtiles = ("  << tile.numXTiles << "," << tile.numYTiles << "," << tile.numZTiles << ")" << std::endl;
  // }

  // Ghosting Regions are identified by their TilingIndices.
  // We need find which regions our domain and our range occupies (3D raster).
  std::vector<BoxedTilingIndices> myDomainOccupanciesDDR;
  std::vector<BoxedTilingIndices> myRangeOccupanciesDDR;
  rasterizeToTilingDDRobust(tile, pRank, local_domain, myDomainOccupanciesDDR);
  rasterizeToTilingDDRobust(tile, pRank, local_range, myRangeOccupanciesDDR);

  // if ((pRank == 0) || (pRank == 63)) {
  //   std::cout << "p_" << pRank << ": myDomainOccupancies = {";
  //   for (auto bti : myDomainOccupanciesDDR) {
  //     std::cout << "{(" << bti.xIdx << "," << bti.yIdx << "," <<  bti.zIdx << ") " << bti.box << "} ";
  //   }
  //   std::cout <<"}" << std::endl;
  // }

  timeBreakdown.findNeighbors_rasterize += stopwatch.split();

  // In role of Decomp Region: Send Ghosting Regions corresponding myDomainOccupancies info.
  // In role of Decomp Region: Send Ghosting Regions corresponding myRangeOccupancies info.
  // In role of Ghosting Region: Receive Domain occupancy info.
  // In role of Ghosting Region: Receive Range occupancy info.
  ResidencyBoxMapT myDomainResidentsDDR;
  ResidencyBoxMapT myRangeResidentsDDR;
  findGhostingRegionResidentsDDRobust(tile, mpiComm, pRank, pSize, myDomainOccupanciesDDR, myDomainResidentsDDR);
  findGhostingRegionResidentsDDRobust(tile, mpiComm, pRank, pSize, myRangeOccupanciesDDR, myRangeResidentsDDR);

  timeBreakdown.findNeighbors_ghostingRegionResidents += stopwatch.split();

  // Now can witness all pairs of (domain_rank, range_rank) neighbor relationships that exist
  // on my Ghosting Region....
  // In role of Ghosting Region: send each range_rank the (domain_aabb, domain_rank) information
  // about its neighbor relationships I witness.
  // In role of Decomp Region: Receive my neighbor information.
  std::vector<OwnedBox3D> domainRanksNeighboringMyRangeDDR, rangeRanksNeighboringMyDomainDDR;
  witnessAndComputeNeighborRanksDDRobust(mpiComm, pSize, myDomainResidentsDDR, myRangeResidentsDDR,
                                         domainRanksNeighboringMyRangeDDR, rangeRanksNeighboringMyDomainDDR);

  timeBreakdown.findNeighbors_witnessAndComputeNeighbors += stopwatch.split();

  communicateDomainNeighborObjectBoundingBoxesDDRobust(domainRanksNeighboringMyRangeDDR, rangeRanksNeighboringMyDomainDDR,
                                                       myDomainOccupanciesDDR, boxesA_proc_box_array);

  timeBreakdown.findNeighbors_communicateNeighborObjectBBs += stopwatch.split();
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
  RangeBoxType>::intersectLocalRangeBoxesWithNeighborBVsDDRobust(
    std::vector<stk::search::ObjectBoundingBox_T<Box3D> >& boxA_proc_box_array,
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
  stk::search::ProximitySearchTree_T<Box3D> boxA_box_hierarchy(
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
    std::vector<std::vector<BoxIdPair> > send_list,
    std::vector<std::vector<BoxIdPair> >& recv_list,
    experimental::GhostingSearchTimeBreakdown& timeBreakdown,
    SplitTimer& stopwatch) {
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

  if (findNeighborsAlgorithm == SCALABLE_AND_DISCONTINUOUS_DECOMP_ROBUST_FIND_NEIGHBORS) {
    std::vector<stk::search::ObjectBoundingBox_T<Box3D> > boxA_proc_box_array;
    findNeighborsFromScratchDDRobust(local_domain, local_range, boxA_proc_box_array,
                                     timeBreakdown, stopwatch);
    std::vector<std::vector<BoxIdPair> > send_list(pSize);
    intersectLocalRangeBoxesWithNeighborBVsDDRobust(boxA_proc_box_array, numBoxRange,
                                                    local_range, send_list, rangeBoxes,
                                                    timeBreakdown, stopwatch);
    commIntersectingRangeBoxesToNeighbors(send_list, recv_list, timeBreakdown,
                                          stopwatch);
  }
  else if (findNeighborsAlgorithm == SCALABLE_FIND_NEIGHBORS){
    std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > boxA_proc_box_array(pSize);
    findNeighborsFromScratchScalable(local_domain, local_range, boxA_proc_box_array,
                                    timeBreakdown, stopwatch);
    std::vector<std::vector<BoxIdPair> > send_list(pSize);
    intersectLocalRangeBoxesWithNeighborBVs(boxA_proc_box_array, numBoxRange,
                                            local_range, send_list, rangeBoxes,
                                            timeBreakdown, stopwatch);
    commIntersectingRangeBoxesToNeighbors(send_list, recv_list, timeBreakdown,
                                          stopwatch);
  }
  else { // ORIGINAL_FIND_NEIGHBORS
    std::vector<stk::search::ObjectBoundingBox_T<DomainBox> > boxA_proc_box_array(pSize);
    findNeighborsFromScratchScalable(local_domain, local_range, boxA_proc_box_array,
                                    timeBreakdown, stopwatch);
    std::vector<std::vector<BoxIdPair> > send_list(pSize);
    intersectLocalRangeBoxesWithNeighborBVs(boxA_proc_box_array, numBoxRange,
                                            local_range, send_list, rangeBoxes,
                                            timeBreakdown, stopwatch);
    commIntersectingRangeBoxesToNeighbors(send_list, recv_list, timeBreakdown,
                                          stopwatch);
  }

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

#endif // PROTOTYPE_SEARCH_UTILS_INSTRUMENTED_H_
