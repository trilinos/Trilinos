#ifndef KRINO_KRINO_KRINO_LIB_AKRI_MESHINTERFACE_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_MESHINTERFACE_HPP_

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino {

class MeshInterface {
public:
  virtual ~MeshInterface() {}
  virtual void populate_mesh(const stk::mesh::BulkData::AutomaticAuraOption auto_aura_option = stk::mesh::BulkData::AUTO_AURA) = 0;
  virtual stk::mesh::MetaData & meta_data() = 0;
  virtual const stk::mesh::MetaData & meta_data() const = 0;
  virtual stk::mesh::BulkData & bulk_data() = 0;
  virtual const stk::mesh::BulkData & bulk_data() const = 0;
};

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_MESHINTERFACE_HPP_ */
