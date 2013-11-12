#include "Gmesh_STKmesh_Fixture.hpp"

#include <stk_mesh/base/Part.hpp>

#include <stk_mesh/base/MetaData.hpp>

#include <stk_io/IossBridge.hpp>

#include <generated/Iogn_DatabaseIO.h>
#include <generated/Iogn_GeneratedMesh.h>

#include <init/Ionit_Initializer.h>
#include <Ioss_Region.h>

#include <Shards_BasicTopologies.hpp>

#include <stdexcept>
#include <sstream>

using namespace stk::io::util;

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
  m_mesh_data.open_mesh_database(gmesh_spec, "generated");
  m_mesh_data.create_input_mesh();

  const Iogn::DatabaseIO* database =
    dynamic_cast<const Iogn::DatabaseIO*>(m_mesh_data.input_io_region()->get_database());
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

