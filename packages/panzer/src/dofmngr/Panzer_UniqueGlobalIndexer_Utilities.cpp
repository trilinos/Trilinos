#include "Panzer_config.hpp"

#include "Panzer_UniqueGlobalIndexer_Utilities.hpp"

using Teuchos::RCP;

namespace panzer {
namespace orientation_helpers {

void computePatternEdgeIndices(const FieldPattern & pattern,std::vector<std::pair<int,int> > & edgeIndices)
{
   unsigned dim = 1;
   shards::CellTopology cellTopo = pattern.getCellTopology();
   for(unsigned e=0;e<cellTopo.getEdgeCount();e++) {
      // get local vertex ids for a this edge
      unsigned local_v0 = cellTopo.getNodeMap(dim,e,0);
      unsigned local_v1 = cellTopo.getNodeMap(dim,e,1);

      // get sub cell indices for geometric pattern
      const std::vector<int> & v0_indices = pattern.getSubcellIndices(0,local_v0);
      const std::vector<int> & v1_indices = pattern.getSubcellIndices(0,local_v1);

      TEUCHOS_ASSERT(v0_indices.size()>0); // there must be a node
      TEUCHOS_ASSERT(v1_indices.size()>0); // there must be a node

      // take the first index on each vertex and make a edge lookup
      edgeIndices.push_back(std::make_pair(v0_indices[0],v1_indices[0]));
   }
}

} // end orientation_helpers
} // end panzer
