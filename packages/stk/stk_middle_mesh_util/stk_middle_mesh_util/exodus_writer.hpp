
#ifndef STK_MIDDLE_MESH_UTILS_STK_INTERFACE_EXODUS_WRITER
#define STK_MIDDLE_MESH_UTILS_STK_INTERFACE_EXODUS_WRITER

#include "stk_middle_mesh/field.hpp"
#include "stk_middle_mesh/mesh.hpp"
#include "field_output_adaptor.hpp"

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
    ExodusWriter(std::shared_ptr<mesh::Mesh> mesh, const std::vector<FieldOutputAdaptorPtr>& fields = {})
      : m_bulkDataOutPtr(::stk::mesh::MeshBuilder(mesh->get_comm()).set_spatial_dimension(3).create())
      , m_metaDataOutPtr(m_bulkDataOutPtr->mesh_meta_data_ptr())
      , m_mesh(mesh)
      , m_meshFields(fields)
    {
      initialize_output_mesh();
      create_part("block_1");
      declare_fields();
    }

    void write(const std::string& fname)
    {
      auto stkVerts = create_part_verts();
      create_part_elements(stkVerts);
      copy_fields_to_stk_mesh();
      write_output_stk(fname);
    }

  private:
    using StkFieldType    = ::stk::mesh::Field<double>;

    struct StkFields
    {
      StkFieldType* vertField = nullptr;
      StkFieldType* elemField = nullptr;
    };

    enum class PartToUse
    {
      Vert,
      Triangle,
      Quad
    };

    void initialize_output_mesh();

    void create_part(const std::string& name);

    void declare_fields();

    mesh::FieldPtr<::stk::mesh::Entity> create_part_verts();

    void create_part_elements(mesh::FieldPtr<::stk::mesh::Entity> stkVerts);

    void copy_fields_to_stk_mesh();

    void copy_field_to_stk_mesh(PartToUse partEnum, FieldOutputAdaptorPtr meshField, StkFieldType& stkField);

    void write_output_stk(const std::string& fname);

    std::shared_ptr<::stk::mesh::BulkData> m_bulkDataOutPtr;
    std::shared_ptr<::stk::mesh::MetaData> m_metaDataOutPtr;
    stk::mesh::Part* m_vertPart = nullptr;
    stk::mesh::Part* m_triPart  = nullptr;
    stk::mesh::Part* m_quadPart = nullptr;
    std::shared_ptr<mesh::Mesh> m_mesh;
    std::vector<FieldOutputAdaptorPtr> m_meshFields;
    std::vector<StkFields> m_stkFields;
};

} // namespace impl

} // namespace stk_interface
} // namespace middle_mesh
} // namespace stk
#endif

