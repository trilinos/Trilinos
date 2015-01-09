#ifndef GeometryVerifier_hpp
#define GeometryVerifier_hpp

#include <iostream>


#include "Shards_CellTopology.hpp"
//#include "Teuchos_GlobalMPISession.hpp"

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

//using namespace Intrepid;


namespace stk_classic
{
  namespace percept
  {


    class GeometryVerifier
    {
      //typedef std::set<std::pair<char > invalid_edge_set_type;
      bool m_dump;
      double m_badJacobian;  // default and settable value for bad jacobian
      double getEquiVol(shards::CellTopology& cell_topo);

    public:
      GeometryVerifier(bool dump=false, double badJac=1.e-10);
      bool isGeometryBad(stk_classic::mesh::BulkData& bulk, bool printTable=false);
    };

  }//namespace percept
}//namespace stk_classic


#endif
