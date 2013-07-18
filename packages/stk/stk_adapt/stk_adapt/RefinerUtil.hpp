#ifndef stk_adapt_RefinerUtil_hpp
#define stk_adapt_RefinerUtil_hpp

#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <cmath>
#include <utility>
#include <math.h>
#include <map>
#include <set>
#include <vector>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include <stk_percept/stk_mesh.hpp>
#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/ProgressMeter.hpp>
#include <stk_adapt/UniformRefinerPattern.hpp>
#include <stk_adapt/Refiner.hpp>
#include <stk_adapt/NodeRegistry.hpp>

#include <stk_adapt/SubDimCell.hpp>

#include <stk_adapt/RefinementInfoByType.hpp>


namespace stk {
namespace adapt {

class RefinerUtil
{
public:

  static BlockNamesType 
  getBlockNames(std::string& block_name, unsigned proc_rank, percept::PerceptMesh& eMesh);


  static BlockNamesType 
  correctBlockNamesForPartPartConsistency(percept::PerceptMesh& eMesh, BlockNamesType& blocks);

  /// For use with the Encore/Percept interface.
  /// for each element in elements_to_unref, add parents to the list, repeatedly, 
  ///   for num_levels_to_add times, thus adding grand-parents, great-grandparents, ...
  /// If num_levels_to_add is < 0, add all up to root (ie, infinite number of levels)
  /// Throw an error if the parents are already in the list.

  static void
  addAncestorsToUnrefineList(percept::PerceptMesh& eMesh, int num_levels_to_add, ElementUnrefineCollection& elements_to_unref);

};

}
}

#endif
