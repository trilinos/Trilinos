#ifndef BasisTable_hpp
#define BasisTable_hpp

#include "Teuchos_RCP.hpp"
#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>
#include <Shards_CellTopology.hpp>
#include <stk_percept/function/MDArray.hpp>


namespace Intrepid {
	template<class Scalar, class ArrayScalar>
	class Basis;
}

namespace stk {
  namespace percept {
    class BasisTable
    {
    public:
      typedef Intrepid::Basis<double, MDArray > BasisType;
      typedef Teuchos::RCP<BasisType>           BasisTypeRCP;
      typedef std::map<unsigned, BasisTypeRCP > BasisTableMap;
      static BasisTypeRCP getBasis(shards::CellTopology& topo);
      static void setupBasisTable();

    private:
      static BasisTableMap m_basisTable;


    };

  }
}

#endif
