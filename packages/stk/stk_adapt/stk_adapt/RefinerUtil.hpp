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
#include <stk_adapt/Colorer.hpp>

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
};

}
}

#endif
