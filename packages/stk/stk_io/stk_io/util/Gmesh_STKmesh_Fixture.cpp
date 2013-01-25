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

static const size_t spatial_dimension = 3;

///////////////////////////////////////////////////////////////////////////////
Gmesh_STKmesh_Fixture::Gmesh_STKmesh_Fixture(stk::ParallelMachine comm,
                                             const std::string& gmesh_spec)
///////////////////////////////////////////////////////////////////////////////
{
  m_mesh_data.create_input_mesh("generated", gmesh_spec, comm);

  const Iogn::DatabaseIO* database =
    dynamic_cast<const Iogn::DatabaseIO*>(m_mesh_data.input_io_region()->get_database());

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

