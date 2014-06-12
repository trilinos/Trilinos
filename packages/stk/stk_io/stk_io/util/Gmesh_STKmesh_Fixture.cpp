
#include "Gmesh_STKmesh_Fixture.hpp"
#include <Ioss_Property.h>              // for Property
#include <Ioss_Region.h>                // for Region
#include <generated/Iogn_DatabaseIO.h>  // for DatabaseIO
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include "Ioss_DatabaseIO.h"            // for DatabaseIO
#include "Teuchos_RCP.hpp"              // for RCP::operator->
#include "Teuchos_RCPDecl.hpp"          // for RCP
#include "mpi.h"                        // for ompi_communicator_t
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
#include "stk_io/StkMeshIoBroker.hpp"   // for StkMeshIoBroker
#include "stk_mesh/base/Types.hpp"      // for PartVector
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
namespace stk { namespace mesh { struct ConnectivityMap; } }



namespace stk {
namespace io {
namespace util {
///////////////////////////////////////////////////////////////////////////////
Gmesh_STKmesh_Fixture::Gmesh_STKmesh_Fixture(   stk::ParallelMachine comm
						, const std::string& gmesh_spec
						, bool use_64bit_int_IO_api
						, stk::mesh::ConnectivityMap * connectivity_map
						)
  : m_mesh_data(comm, connectivity_map)

///////////////////////////////////////////////////////////////////////////////
{
  if (use_64bit_int_IO_api) {
    m_mesh_data.property_add(Ioss::Property("INTEGER_SIZE_API", 8));
  }
  size_t ifh = m_mesh_data.add_mesh_database(gmesh_spec, "generated", stk::io::READ_MESH);
  m_mesh_data.set_active_mesh(ifh);
  m_mesh_data.create_input_mesh();

  const Iogn::DatabaseIO* database =
    dynamic_cast<const Iogn::DatabaseIO*>(m_mesh_data.get_input_io_region()->get_database());
//  database->set_int_byte_size_api(Ioss::USE_INT64_API);

  // get face parts names; need to convert these to strings
  const std::vector<std::string> sideset_names = database->get_sideset_names();

  for (std::vector<std::string>::const_iterator itr = sideset_names.begin();
       itr != sideset_names.end(); ++itr) {
    m_sideset_names.push_back(*itr);
    m_sideset_parts.push_back(m_mesh_data.meta_data().get_part(*itr));
  }
}

///////////////////////////////////////////////////////////////////////////////
void Gmesh_STKmesh_Fixture::commit()
///////////////////////////////////////////////////////////////////////////////
{
  m_mesh_data.meta_data().commit();
  m_mesh_data.populate_bulk_data();
}

}}}
