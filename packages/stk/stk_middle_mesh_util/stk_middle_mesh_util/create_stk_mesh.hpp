#ifndef STK_MIDDLE_MESH_UTILS_CREATE_STK_MESH
#define STK_MIDDLE_MESH_UTILS_CREATE_STK_MESH

#include "stk_io/DatabasePurpose.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/MetaData.hpp"

#include "stk_middle_mesh/field.hpp"
#include "stk_middle_mesh/mesh.hpp"

namespace stk {
namespace middle_mesh {
namespace stk_interface {

struct MeshPart
{
    using VertIdType = double; // Note: this should be an integer
                                    //       type, but IOSS doesn't
                                    //       support integer valued
                                    //       Attribute fields
                                    //       integers up to 2^53 are
                                    //       exactly representable by
                                    //       doubles, so use this for now  
    using MeshFieldPtr = mesh::FieldPtr<::stk::mesh::SideSetEntry>;
    std::shared_ptr<mesh::Mesh> mesh;
    MeshFieldPtr stkEls; // TODO: should be stk_faces
    stk::mesh::Field<VertIdType>* stkVertField;
    stk::mesh::Part* stkPart;
};

class StkMeshCreator
{
  public:
    using VertIdType = MeshPart::VertIdType;
    using FieldType    = stk::mesh::Field<VertIdType>;
    using MeshFieldPtr = MeshPart::MeshFieldPtr;
    static std::string vertex_field_name();

    explicit StkMeshCreator(const std::string& fname, const std::string& autodecompMethod="NONE",
                            MPI_Comm comm = MPI_COMM_WORLD) // CHECK: ALLOW MPI_COMM_WORLD
        : m_bulkDataPtr(::stk::mesh::MeshBuilder(comm).set_spatial_dimension(3).create()),
          m_metaDataPtr(m_bulkDataPtr->mesh_meta_data_ptr()),
          m_autodecompMethod(autodecompMethod)
    {
      declare_stk_vert_field();
      load_mesh(fname);
    }

    explicit StkMeshCreator(::stk::mesh::BulkData * bulkDataPtr)
      : m_bulkDataPtr(std::shared_ptr<::stk::mesh::BulkData>(bulkDataPtr, [](auto *){})) //no delete antipattern
      , m_metaDataPtr(m_bulkDataPtr->mesh_meta_data_ptr())
    {
      declare_stk_vert_field();
    }

    explicit StkMeshCreator(std::shared_ptr<::stk::mesh::BulkData> bulkDataPtr)
      : m_bulkDataPtr(bulkDataPtr) //no delete antipattern
      , m_metaDataPtr(m_bulkDataPtr->mesh_meta_data_ptr())
    {
      declare_stk_vert_field();
    }

    MeshPart create_mesh_from_part(const std::string& name);

    // copies the coordinates from the Mesh to the STK mesh
    void write_back_coords(std::shared_ptr<mesh::Mesh> mesh, const std::string& name);

    stk::mesh::MetaData& get_meta_data() { return *m_metaDataPtr; }

    stk::mesh::BulkData& get_bulk_data() { return *m_bulkDataPtr; }

    std::shared_ptr<stk::mesh::MetaData> get_meta_data_ptr() { return m_metaDataPtr; }

    std::shared_ptr<stk::mesh::BulkData> get_bulk_data_ptr() { return m_bulkDataPtr; }

  private:
    void declare_stk_vert_field();

    void load_mesh(const std::string& fname);

    void create_nodes(std::shared_ptr<mesh::Mesh> mesh);

    void create_faces_from_sideset(std::shared_ptr<mesh::Mesh> mesh, MeshFieldPtr stkEls);

    void create_faces_from_shells(std::shared_ptr<mesh::Mesh> mesh, MeshFieldPtr stkEls);

    void setup_vert_sharing(std::shared_ptr<mesh::Mesh> mesh);

    void create_edges(std::shared_ptr<mesh::Mesh> mesh, stk::mesh::EntityRank rank);

    void setup_edge_sharing(std::shared_ptr<mesh::Mesh> mesh, MeshFieldPtr stkEls);

    std::shared_ptr<::stk::mesh::BulkData> m_bulkDataPtr;
    std::shared_ptr<::stk::mesh::MetaData> m_metaDataPtr;
    stk::mesh::Part* m_part;
    FieldType* m_stkNodeField;
    std::string m_autodecompMethod;
};


} // namespace stk_interface
} // namespace middle_mesh
} // namespace stk
#endif

