#ifndef stk_io_util_Gmesh_STKmesh_Fixture_hpp
#define stk_io_util_Gmesh_STKmesh_Fixture_hpp

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_io/util/UseCase_mesh.hpp>

#include <vector>
#include <string>

namespace stk {
namespace io {
namespace util {

/**
 * This class implements a Stk-mesh based fixture that uses a generated
 * mesh as the basis of the fixture.
 */
class Gmesh_STKmesh_Fixture
{
 public:

  /**
   * Construct a fixture. Note that the fixture won't be completed until commit
   * is called; the intent is to give the client a chance to make additional
   * changes to the meta-data.
   *
   * @param comm The comm object for all processors using the fixture
   * @param gmesh_spec The specification for the mesh. See Iogn::GeneratedMesh
   * for documentation on how to specify meshes.
   */
  Gmesh_STKmesh_Fixture(stk::ParallelMachine comm, const std::string& gmesh_spec);

  /**
   * Commits the meta-data of the mesh and populates the bulk-data. Don't call
   * this until you are done modifying the meta-data.
   */
  void commit();

  /**
   * For a given surface, return the number of elements in the surface
   *
   * @param surf_id The surface we are interested in.
   */
  size_t getSurfElemCount(size_t surf_id) const;

  /**
   * For a given surface, return the relevant dimension and expected value
   * of that dimension. For example, for surface PY, (1, m_num_y) would be
   * returned; 1 refers to the Y dimension and m_num_y is the expected
   * Y-coordinate value for all the nodes on the PY surface.
   *
   * @surf_id The surface we are interested in.
   */
  std::pair<int, double> getSurfCoordInfo(size_t surf_id) const;

  /**
   * Get the total number of side entities in this mesh.
   */
  size_t getSideCount() const;

  /**
   * Get the total number of elements in this mesh.
   */
  size_t getElemCount() const;

  /**
   * Get the total number of nodes in this mesh.
   */
  size_t getNodeCount() const;

  /**
   * Get the names of all the sideset parts.
   */
  const std::vector<std::string> & getSidesetNames() const
  { return m_sideset_names; }

  /**
   * Get all the sideset parts.
   */
  const stk::mesh::PartVector & getSideParts() const
  { return m_sideset_parts; }

  /**
   * Get a reference to the meta data for the stk-mesh.
  const stk::mesh::MetaData & getMetaData() const
  { return m_meta_data.get_meta_data(m_meta_data); }
   */

  stk::mesh::MetaData & getMetaData()
  { return m_meta_data.get_meta_data(m_meta_data); }

  const stk::mesh::fem::FEMMetaData & getFEMMetaData() const
  { return m_meta_data; }

  stk::mesh::fem::FEMMetaData & getFEMMetaData()
  { return m_meta_data; }

  /**
   * Get a reference to the bulk data for the stk-mesh.
   */
  const stk::mesh::BulkData & getBulkData() const
  { return m_bulk_data; }

  stk::mesh::BulkData & getBulkData()
  { return m_bulk_data; }

 private:
  ///> The meta data for the stk-mesh
  stk::mesh::fem::FEMMetaData m_meta_data;

  ///> The bulk data for the stk-mesh
  stk::mesh::BulkData m_bulk_data;

  /**
   * The mesh-data for the mesh. This is a special object that maintains some
   * state between the meta data and bulk data portions of the mesh generation
   * process for use cases.
   */
  stk::io::util::MeshData m_mesh_data;

  ///> The names of all the side parts
  std::vector<std::string> m_sideset_names;

  ///> Collection of all the side parts
  stk::mesh::PartVector m_sideset_parts;

  ///> The number of elements in the X dimension for the mesh
  int m_num_x;

  ///> The number of elements in the Y dimension for the mesh
  int m_num_y;

  ///> The number of elements in the Z dimension for the mesh
  int m_num_z;
};

}//namespace util
}//namespace io
}//namespace stk

#endif
