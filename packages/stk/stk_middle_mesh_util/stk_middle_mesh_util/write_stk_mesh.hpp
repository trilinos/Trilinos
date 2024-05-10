
#ifndef STK_MIDDLE_MESH_UTILS_WRITE_STK_MESH_H
#define STK_MIDDLE_MESH_UTILS_WRITE_STK_MESH_H

#include "stk_io/DatabasePurpose.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_topology/topology.hpp"

#include "stk_middle_mesh/mesh.hpp"
#include "stk_middle_mesh/nonconformal_abstract.hpp"
#include "create_stk_mesh.hpp"

namespace stk {
namespace middle_mesh {
namespace stk_interface {
namespace impl {

class StkMeshWriter
{
    using MeshFieldPtr = StkMeshCreator::MeshFieldPtr;

  public:
    using VertIdType = StkMeshCreator::VertIdType;
    using FieldType       = StkMeshCreator::FieldType;

    StkMeshWriter(stk::mesh::MetaData& metaDataIn, stk::mesh::BulkData& bulkDataIn, stk::mesh::MetaData& metaDataOut,
                  stk::mesh::BulkData& bulkDataOut, FieldType* gidField)
      : // m_meta_data_in(meta_data_in),
      m_bulkDataIn(bulkDataIn)
      , m_metaDataOut(metaDataOut)
      , m_bulkDataOut(bulkDataOut)
      , m_gidField(gidField)
    {}

    // writes meshes defined by classes to stk Part
    void populate_part(std::shared_ptr<nonconformal::impl::NonconformalAbstract> nonconformalAbstract, MeshPart& pmesh1,
                       MeshPart& pmesh2, stk::mesh::Part* part);

  private:
    using MNodeFieldPtr = mesh::FieldPtr<stk::mesh::Entity>;

    MNodeFieldPtr create_part_verts(std::shared_ptr<nonconformal::impl::NonconformalAbstract> nonconformalAbstract,
                                    stk::mesh::Part* part);

    void create_part_elements(std::shared_ptr<nonconformal::impl::NonconformalAbstract> nonconformalAbstract,
                              MNodeFieldPtr stkVerts, stk::mesh::Part* part, std::vector<stk::mesh::Entity>& entities);

    void write_gid_field(std::shared_ptr<nonconformal::impl::NonconformalAbstract> nonconformalAbstract,
                         std::vector<stk::mesh::Entity> entities, MeshPart& pmesh1, MeshPart& pmesh2);

    // Stk mesh that the GID field describes
    // stk::mesh::MetaData& m_meta_data_in;
    stk::mesh::BulkData& m_bulkDataIn;

    // Stk mesh writing to
    stk::mesh::MetaData& m_metaDataOut;
    stk::mesh::BulkData& m_bulkDataOut;
    FieldType* m_gidField;
};

} // namespace impl

} // namespace stk_interface
} // namespace middle_mesh
} // namespace stk
#endif

