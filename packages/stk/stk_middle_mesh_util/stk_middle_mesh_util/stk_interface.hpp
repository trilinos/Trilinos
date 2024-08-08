
#ifndef STK_MIDDLE_MESH_UTILS_STK_INTERFACE
#define STK_MIDDLE_MESH_UTILS_STK_INTERFACE

#include "stk_middle_mesh/mesh_input.hpp"
#include <set>
#include <string>
#include <vector>

#include "stk_io/DatabasePurpose.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_topology/topology.hpp"

#include "stk_middle_mesh/nonconformal_standard.hpp"
#include "create_stk_mesh.hpp"
#include "write_stk_mesh.hpp"

namespace stk {
namespace middle_mesh {
namespace stk_interface {
namespace impl {

class StkInterface
{
    using FieldType    = StkMeshCreator::FieldType;
    using MeshFieldPtr = StkMeshCreator::MeshFieldPtr;

  public:
    explicit StkInterface(const mesh::impl::MeshInput& input)
      : m_input(input)
      , m_stkMeshCreator(input.fnameIn)
      , m_bulkDataOutPtr(stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create())
      , m_metaDataOutPtr(m_bulkDataOutPtr->mesh_meta_data_ptr())
    {
      check_sideset_size();
      check_sideset_uniqueness();
      initialize_output_mesh();

      check_sidesets_exist();
      create_parts();
      create_field();
      check_node_counts();
    }

    void compute_interface(const middle_mesh::impl::NonconformalOpts& opts);

    void write_output();

  private:
    void check_sideset_size();

    void check_sidesets_exist();

    // for each input pair, checks that the number of nodes on the
    // first part plus the number of nodes on the second part equals
    // the number of nodes on the union of both parts.
    // This is useful because it verifies that there are no nodes shared
    // between the two parts (which happens sometimes depending on how
    // the user generated the mesh)
    void check_node_counts();

    void check_name_exists(const std::set<std::string>& names, const std::string& name);

    void check_sideset_uniqueness();

    void check_name_uniqueness(std::set<std::string>& names, const std::string& name);

    void initialize_output_mesh();

    // create one part for each nonconformal interface
    // must be called while MetaData is still modifiable
    void create_parts();

    // creates a field on all the newly-created parts, storing two integer
    // values per element
    // must be called while MetaData is still modifiable
    void create_field();

    // write_names: if true, write names of sidesets that make up the
    //              nonconformal interfaces to the file
    void write_output_stk(const std::string& fname, stk::mesh::BulkData& bulkData, bool writeNames = false);

    std::string get_part_name(mesh::impl::MeshInput::NamePair& p);

    stk::mesh::MetaData& get_input_meta_data() { return m_stkMeshCreator.get_meta_data(); }

    stk::mesh::BulkData& get_input_bulk_data() { return m_stkMeshCreator.get_bulk_data(); }

    mesh::impl::MeshInput m_input;
    StkMeshCreator m_stkMeshCreator;
    std::shared_ptr<stk::mesh::BulkData> m_bulkDataOutPtr;
    std::shared_ptr<stk::mesh::MetaData> m_metaDataOutPtr;
    std::vector<stk::mesh::Part*> m_parts; // newly created parts on
                                           // mesh_out
    FieldType* m_gidField;                 // classification information relating
                                           // mesh_out to mesh_in
};

} // namespace impl
} // namespace stk_interface
} // namespace middle_mesh
} // namespace stk
#endif

