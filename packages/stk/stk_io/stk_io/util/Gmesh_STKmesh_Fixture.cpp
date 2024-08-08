// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "Gmesh_STKmesh_Fixture.hpp"
#include <Ioss_Property.h>                 // for Property
#include <Ioss_Region.h>                   // for Region
#include <generated/Iogn_DatabaseIO.h>     // for DatabaseIO
#include <cstddef>                         // for size_t
#include <stk_mesh/base/MetaData.hpp>      // for MetaData
#include "Ioss_DatabaseIO.h"               // for DatabaseIO
#include "mpi.h"                           // for ompi_communicator_t
#include "stk_io/DatabasePurpose.hpp"      // for READ_MESH
#include "stk_io/StkMeshIoBroker.hpp"      // for StkMeshIoBroker
#include "stk_mesh/base/Types.hpp"         // for PartVector
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace io {
namespace util {

Gmesh_STKmesh_Fixture::Gmesh_STKmesh_Fixture(stk::ParallelMachine comm,
                                             const std::string& gmesh_spec,
                                             bool use_64bit_int_IO_api)
  : m_mesh_data(comm)
{
  if (use_64bit_int_IO_api) {
    m_mesh_data.property_add(Ioss::Property("INTEGER_SIZE_API", 8));
  }
  size_t ifh = m_mesh_data.add_mesh_database(gmesh_spec, "generated", stk::io::READ_MESH);
  m_mesh_data.set_active_mesh(ifh);
  m_mesh_data.create_input_mesh();

  auto iossRegion = m_mesh_data.get_input_ioss_region();
  const Iogn::DatabaseIO* database = dynamic_cast<const Iogn::DatabaseIO*>(iossRegion->get_database());
  //  database->set_int_byte_size_api(Ioss::USE_INT64_API);

  // get face parts names; need to convert these to strings
  const std::vector<std::string> sideset_names = database->get_sideset_names();

  for (std::vector<std::string>::const_iterator itr = sideset_names.begin(); itr != sideset_names.end(); ++itr) {
    m_sideset_names.push_back(*itr);
    m_sideset_parts.push_back(m_mesh_data.meta_data().get_part(*itr));
  }
}

void Gmesh_STKmesh_Fixture::commit()
{
  m_mesh_data.meta_data().commit();
  m_mesh_data.populate_bulk_data();
}

namespace simple_fields {

Gmesh_STKmesh_Fixture::Gmesh_STKmesh_Fixture(stk::ParallelMachine comm,
                                             const std::string& gmesh_spec,
                                             bool use_64bit_int_IO_api)
  : m_mesh_data(comm)
{
  if (use_64bit_int_IO_api) {
    m_mesh_data.property_add(Ioss::Property("INTEGER_SIZE_API", 8));
  }
  size_t ifh = m_mesh_data.add_mesh_database(gmesh_spec, "generated", stk::io::READ_MESH);
  m_mesh_data.set_active_mesh(ifh);
  m_mesh_data.create_input_mesh();

  auto iossRegion = m_mesh_data.get_input_ioss_region();
  const Iogn::DatabaseIO* database = dynamic_cast<const Iogn::DatabaseIO*>(iossRegion->get_database());
  //  database->set_int_byte_size_api(Ioss::USE_INT64_API);

  // get face parts names; need to convert these to strings
  const std::vector<std::string> sideset_names = database->get_sideset_names();

  for (std::vector<std::string>::const_iterator itr = sideset_names.begin(); itr != sideset_names.end(); ++itr) {
    m_sideset_names.push_back(*itr);
    m_sideset_parts.push_back(m_mesh_data.meta_data().get_part(*itr));
  }
}

void Gmesh_STKmesh_Fixture::commit()
{
  m_mesh_data.meta_data().commit();
  m_mesh_data.populate_bulk_data();
}

}

}}}
