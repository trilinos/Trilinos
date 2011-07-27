#ifndef stk_adapt_RefinementInfoByType_hpp
#define stk_adapt_RefinementInfoByType_hpp

#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>
#include <stdint.h>
#include <limits>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include <stk_percept/PerceptMesh.hpp>
#include <stk_adapt/UniformRefinerPattern.hpp>



namespace stk {
  namespace adapt {

    struct RefinementInfoByType
    {
      //typedef uint64_t RefinementInfoCount ;
      typedef unsigned RefinementInfoCount ;

      RefinementInfoCount m_numOrigElems;
      RefinementInfoCount m_numNewElems;
      shards::CellTopology m_topology;
      RefinementInfoCount m_numOrigNodes;
      RefinementInfoCount m_numNewNodes;

      static void printTable(std::ostream& os, std::vector< RefinementInfoByType >& refinementInfoByType, int iRefinePass, bool printAll = false);
      static void countCurrentNodes(stk::percept::PerceptMesh& eMesh, std::vector< RefinementInfoByType >& refinementInfoByType);

    };

  }
}

#endif
