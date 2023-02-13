#ifdef STK_BUILT_IN_SIERRA

#ifndef CREATE_STK_MESH
#define CREATE_STK_MESH

#include "stk_io/DatabasePurpose.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_topology/topology.hpp"

#include "field.h"
#include "mesh.h"

namespace stk {
namespace middle_mesh {
namespace stk_interface {
namespace impl {

struct MeshPart
{
    using MeshFieldPtr = mesh::FieldPtr<::stk::mesh::SideSetEntry>;
    std::shared_ptr<mesh::Mesh> mesh;
    MeshFieldPtr stkEls; // TODO: should be stk_faces
};

class StkMeshCreator
{
  public:
    using FieldScalarType = double; // Note: this should be an integer
                                    //       type, but IOSS doesn't
                                    //       support integer valued
                                    //       Attribute fields
                                    //       integers up to 2^53 are
                                    //       exactly representable by
                                    //       doubles, so use this for now
    using FieldType    = ::stk::mesh::Field<FieldScalarType>;
    using MeshFieldPtr = MeshPart::MeshFieldPtr;

    explicit StkMeshCreator(const std::string& fname)
      : m_bulkDataPtr(::stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create())
      , m_metaDataPtr(m_bulkDataPtr->mesh_meta_data_ptr())
    {
      declare_stk_vert_field();
      load_mesh(fname);
    }

    MeshPart create_mesh_from_part(const std::string& name);

    // copies the coordinates from the Mesh to the STK mesh
    void write_back_coords(std::shared_ptr<mesh::Mesh> mesh, const std::string& name);

    ::stk::mesh::MetaData& get_meta_data() { return *m_metaDataPtr; }

    ::stk::mesh::BulkData& get_bulk_data() { return *m_bulkDataPtr; }

  private:
    void declare_stk_vert_field();

    void load_mesh(const std::string& fname);

    void create_nodes(std::shared_ptr<mesh::Mesh> mesh);

    void create_faces_from_sideset(std::shared_ptr<mesh::Mesh> mesh, MeshFieldPtr stkEls);

    void create_faces_from_shells(std::shared_ptr<mesh::Mesh> mesh, MeshFieldPtr stkEls);

    std::shared_ptr<::stk::mesh::BulkData> m_bulkDataPtr;
    std::shared_ptr<::stk::mesh::MetaData> m_metaDataPtr;
    ::stk::mesh::Part* m_part;
    FieldType* m_stkNodeField;
};

} // namespace impl
} // namespace stk_interface
} // namespace middle_mesh
} // namespace stk
#endif

#endif
