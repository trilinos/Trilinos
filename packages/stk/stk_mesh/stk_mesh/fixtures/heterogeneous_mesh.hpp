#ifndef STK_MESH_HETEROGENEOUS_MESH_HPP
#define STK_MESH_HETEROGENEOUS_MESH_HPP

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

namespace stk {
  namespace mesh {
    class MetaData;
    class BulkData;

    namespace fixtures {

      typedef mesh::Field<double,mesh::Cartesian> VectorFieldType ;

      void heterogeneous_mesh_meta_data(stk::mesh::MetaData & meta_data, VectorFieldType & node_coord);

      void heterogeneous_mesh_bulk_data(stk::mesh::BulkData & bulk_data, const VectorFieldType & node_coord);
    }
  }
}
#endif
