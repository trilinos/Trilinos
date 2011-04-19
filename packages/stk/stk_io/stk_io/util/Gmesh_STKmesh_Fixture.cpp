#include "Gmesh_STKmesh_Fixture.hpp"

#include <stk_mesh/base/Part.hpp>

#include <stk_mesh/base/MetaData.hpp>

#include <stk_io/IossBridge.hpp>

#include <generated/Iogn_DatabaseIO.h>
#include <generated/Iogn_GeneratedMesh.h>

#include <init/Ionit_Initializer.h>

#include <Shards_BasicTopologies.hpp>

#include <stdexcept>
#include <sstream>

using namespace stk::io::util;

static const size_t spatial_dimension = 3;

///////////////////////////////////////////////////////////////////////////////
Gmesh_STKmesh_Fixture::Gmesh_STKmesh_Fixture(stk::ParallelMachine comm,
                                             const std::string& gmesh_spec)
///////////////////////////////////////////////////////////////////////////////
  : 
    m_meta_data(spatial_dimension),
    m_bulk_data(m_meta_data.get_meta_data(m_meta_data), comm),
    m_num_x(0),
    m_num_y(0),
    m_num_z(0)
{
  // Initialize IO system.  Registers all element types and storage
  // types and the exodusII default database type.
  Ioss::Init::Initializer init_db;

  stk::io::util::MeshData mesh_data;
  stk::io::util::create_input_mesh("generated", gmesh_spec,
                                   ".", comm,
                                   m_meta_data,
                                   m_mesh_data,
                                   false);

  const Iogn::DatabaseIO* database =
    dynamic_cast<const Iogn::DatabaseIO*>(m_mesh_data.m_region->get_database());

  // compute m_num_{x|y|z}
  m_num_x = database->get_generated_mesh()->get_num_x();
  m_num_y = database->get_generated_mesh()->get_num_y();
  m_num_z = database->get_generated_mesh()->get_num_z();

  // get face parts names; need to convert these to strings
  const std::vector<std::string> sideset_names = database->get_sideset_names();

  for (std::vector<std::string>::const_iterator itr = sideset_names.begin();
       itr != sideset_names.end(); ++itr) {
    m_sideset_names.push_back(*itr);
    m_sideset_parts.push_back(m_meta_data.get_part(*itr));
  }
}

///////////////////////////////////////////////////////////////////////////////
void Gmesh_STKmesh_Fixture::commit()
///////////////////////////////////////////////////////////////////////////////
{
  m_meta_data.commit();

  stk::io::util::populate_bulk_data(m_bulk_data, m_mesh_data, "generated");
}

///////////////////////////////////////////////////////////////////////////////
size_t Gmesh_STKmesh_Fixture::getSurfElemCount(size_t surf_id) const
///////////////////////////////////////////////////////////////////////////////
{
  if (Iogn::GeneratedMesh::PZ < surf_id) {
    throw std::runtime_error("Invalid surface id");
  }

  switch( surf_id )
  {
    case Iogn::GeneratedMesh::MX:
    case Iogn::GeneratedMesh::PX:
      return m_num_y*m_num_z;
      break;
    case Iogn::GeneratedMesh::MY:
    case Iogn::GeneratedMesh::PY:
      return m_num_x*m_num_z;
      break;
    case Iogn::GeneratedMesh::MZ:
    case Iogn::GeneratedMesh::PZ:
      return m_num_x*m_num_y;
      break;
  }
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
std::pair<int, double>
Gmesh_STKmesh_Fixture::getSurfCoordInfo(size_t surf_id) const
///////////////////////////////////////////////////////////////////////////////
{
  if (Iogn::GeneratedMesh::PZ < surf_id) {
    throw std::runtime_error("Invalid surface id");
  }

  switch( surf_id )
  {
    case Iogn::GeneratedMesh::MX:
      return std::make_pair<int, double>(0, 0.0);
      break;
    case Iogn::GeneratedMesh::PX:
      return std::make_pair<int, double>(0, (double)m_num_x);
      break;
    case Iogn::GeneratedMesh::MY:
      return std::make_pair<int, double>(1, 0.0);
      break;
    case Iogn::GeneratedMesh::PY:
      return std::make_pair<int, double>(1, (double)m_num_y);
      break;
    case Iogn::GeneratedMesh::MZ:
      return std::make_pair<int, double>(2, 0.0);
      break;
    case Iogn::GeneratedMesh::PZ:
      return std::make_pair<int, double>(2, (double)m_num_z);
      break;
  }

  return std::make_pair<int, double>(-1, -1.0);
}

///////////////////////////////////////////////////////////////////////////////
size_t Gmesh_STKmesh_Fixture::getSideCount() const
///////////////////////////////////////////////////////////////////////////////
{
  return 2 * (m_num_x*m_num_y + m_num_x*m_num_z + m_num_y*m_num_z);
}

///////////////////////////////////////////////////////////////////////////////
size_t Gmesh_STKmesh_Fixture::getElemCount() const
///////////////////////////////////////////////////////////////////////////////
{
  return m_num_x * m_num_y * m_num_z;
}

///////////////////////////////////////////////////////////////////////////////
size_t Gmesh_STKmesh_Fixture::getNodeCount() const
///////////////////////////////////////////////////////////////////////////////
{
  return (m_num_x+1)*(m_num_y+1)*(m_num_z+1);
}
