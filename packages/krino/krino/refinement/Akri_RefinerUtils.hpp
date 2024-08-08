#ifndef KRINO_KRINO_REFINEMENT_AKRI_REFINERUTILS_HPP_
#define KRINO_KRINO_REFINEMENT_AKRI_REFINERUTILS_HPP_
#include <array>
#include <vector>

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
      childElemDescs[oldSize + i].sideIds[iSide] = (childElementSideId<0) ? -1 : permutedParentSideOrdinals[childElementSideId];
    }
  }
}

}

#endif /* KRINO_KRINO_REFINEMENT_AKRI_REFINERUTILS_HPP_ */
