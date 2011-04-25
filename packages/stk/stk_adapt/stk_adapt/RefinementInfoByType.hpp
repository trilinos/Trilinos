#ifndef stk_adapt_RefinementInfoByType_hpp
#define stk_adapt_RefinementInfoByType_hpp

#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include <stk_percept/PerceptMesh.hpp>
#include <stk_adapt/UniformRefinerPattern.hpp>



namespace stk {
  namespace adapt {

    struct RefinementInfoByType
    {
      unsigned m_numOrigElems;
      unsigned m_numNewElems;
      shards::CellTopology m_topology;

      static void printTable(std::ostream& os, std::vector< RefinementInfoByType >& refinementInfoByType, int iRefinePass, bool printAll = false);

    };

  }
}

#endif
