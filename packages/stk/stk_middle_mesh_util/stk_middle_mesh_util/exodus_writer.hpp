#ifdef STK_BUILT_IN_SIERRA

#ifndef STK_INTERFACE_EXODUS_WRITER
#define STK_INTERFACE_EXODUS_WRITER

#include "field.hpp"
#include "mesh.hpp"

#include "stk_io/DatabasePurpose.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_topology/topology.hpp"

namespace stk {
namespace middle_mesh {
namespace stk_interface {
namespace impl {

class ExodusWriter
{
  public:
    ExodusWriter(std::shared_ptr<mesh::Mesh> mesh)
      : m_bulkDataOutPtr(::stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create())
      , m_metaDataOutPtr(m_bulkDataOutPtr->mesh_meta_data_ptr())
      , m_mesh(mesh)
    {
      m_metaDataOutPtr->use_simple_fields();
      initialize_output_mesh();
      create_part("block_1");
    }

    void write(const std::string& fname)
    {
      auto stkVerts = create_part_verts();
      create_part_elements(stkVerts);
      write_output_stk(fname);
    }

  private:
    void initialize_output_mesh();

    void create_part(const std::string& name);

    mesh::FieldPtr<::stk::mesh::Entity> create_part_verts();

    void create_part_elements(mesh::FieldPtr<::stk::mesh::Entity> stkVerts);

    void write_output_stk(const std::string& fname);

    std::shared_ptr<::stk::mesh::BulkData> m_bulkDataOutPtr;
    std::shared_ptr<::stk::mesh::MetaData> m_metaDataOutPtr;
    ::stk::mesh::Part* m_part = nullptr;
    std::shared_ptr<mesh::Mesh> m_mesh;
};

} // namespace impl

} // namespace stk_interface
} // namespace middle_mesh
} // namespace stk
#endif

#endif
