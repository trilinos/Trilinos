#ifndef KRINO_KRINO_REFINEMENT_AKRI_REFINERUTILS_HPP_
#define KRINO_KRINO_REFINEMENT_AKRI_REFINERUTILS_HPP_
#include <algorithm>
#include <array>
#include <numeric>
#include <vector>
#include <Akri_QualityMetric.hpp>

namespace krino {

template<class CHILDDESCRIPTION, std::size_t NUMREFINEDPARENTNODES, std::size_t NUMNODES, std::size_t NUMSIDES, std::size_t NUMCHILDELEMENTS>
void append_child_elements(const std::array<unsigned,NUMREFINEDPARENTNODES> & permutedParentNodeOrdinals,
    const std::array<unsigned,NUMSIDES> & permutedParentSideOrdinals,
    const std::array<std::array<int,NUMNODES>,NUMCHILDELEMENTS> & childElementNodeIndices,
    const std::array<std::array<int,NUMSIDES>,NUMCHILDELEMENTS> & childElementSideIndices,
    std::vector<CHILDDESCRIPTION> & childElemDescs)
{
  const std::size_t oldSize = childElemDescs.size();
  childElemDescs.resize(oldSize + NUMCHILDELEMENTS);
  for (std::size_t i=0; i<NUMCHILDELEMENTS; ++i)
  {
    for (std::size_t iNode=0; iNode<NUMNODES; ++iNode)
      childElemDescs[oldSize + i].nodeIds[iNode] = permutedParentNodeOrdinals[childElementNodeIndices[i][iNode]];
    for (std::size_t iSide=0; iSide<NUMSIDES; ++iSide)
    {
      const int childElementSideId = childElementSideIndices[i][iSide];
      //Buggy nvidia compiler flags this as "integer conversion resulted in a change of sign": childElemDescs[oldSize + i].sideIds[iSide] = (childElementSideId<0) ? -1 : permutedParentSideOrdinals[childElementSideId];
      if (childElementSideId >= 0)
        childElemDescs[oldSize + i].sideIds[iSide] = permutedParentSideOrdinals[childElementSideId];
      else
        childElemDescs[oldSize + i].sideIds[iSide] = -1;
    }
  }
}

template <size_t SIZE, typename CONTAINER>
std::array<int,SIZE> get_rank_of_nodes_based_on_coordinates(const CONTAINER & nodeCoords)
{
  // initialize original index locations
  std::array<size_t,SIZE> index;
  std::iota(index.begin(), index.end(), 0);

  std::stable_sort(index.begin(), index.end(),
       [&nodeCoords](size_t i1, size_t i2) {return is_less_than_in_x_then_y_then_z(nodeCoords[i1], nodeCoords[i2]);});

  // invert indices to get node rank by coordinates
  std::array<int,SIZE> rank;
  for (size_t i=0; i < SIZE; ++i)
    rank[index[i]] = i;
  return rank;
}

}

#endif /* KRINO_KRINO_REFINEMENT_AKRI_REFINERUTILS_HPP_ */
