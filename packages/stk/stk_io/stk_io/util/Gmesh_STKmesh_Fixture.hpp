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

#ifndef stk_io_util_Gmesh_STKmesh_Fixture_hpp
#define stk_io_util_Gmesh_STKmesh_Fixture_hpp

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stk_io/StkMeshIoBroker.hpp>      // for StkMeshIoBroker
#include <stk_mesh/base/Types.hpp>         // for PartVector
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <string>                          // for string
#include <vector>                          // for vector
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class MetaData; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

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
  Gmesh_STKmesh_Fixture(stk::ParallelMachine comm,
                        const std::string& gmesh_spec,
                        bool use_64bit_int_IO_api=false);

  /**
   * Commits the meta-data of the mesh and populates the bulk-data. Don't call
   * this until you are done modifying the meta-data.
   */
  void commit();

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
   */
  const stk::mesh::MetaData & getMetaData() const
  { return m_mesh_data.meta_data(); }

  stk::mesh::MetaData & getMetaData()
  { return m_mesh_data.meta_data(); }

  /**
   * Get a reference to the bulk data for the stk-mesh.
   */
  const stk::mesh::BulkData & getBulkData() const { return m_mesh_data.bulk_data(); }
        stk::mesh::BulkData & getBulkData()       { return m_mesh_data.bulk_data(); }

 private:
  /**
   * The mesh-data for the mesh. This is a special object that maintains some
   * state between the meta data and bulk data portions of the mesh generation
   * process for use cases.
   */
  stk::io::StkMeshIoBroker m_mesh_data;

  ///> The names of all the side parts
  std::vector<std::string> m_sideset_names;

  ///> Collection of all the side parts
  stk::mesh::PartVector m_sideset_parts;
};

namespace simple_fields {

/**
 * This class implements a Stk-mesh based fixture that uses a generated
 * mesh as the basis of the fixture.
 */
class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
Gmesh_STKmesh_Fixture
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
  Gmesh_STKmesh_Fixture(stk::ParallelMachine comm,
                        const std::string& gmesh_spec,
                        bool use_64bit_int_IO_api=false);

  /**
   * Commits the meta-data of the mesh and populates the bulk-data. Don't call
   * this until you are done modifying the meta-data.
   */
  void commit();

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
   */
  const stk::mesh::MetaData & getMetaData() const
  { return m_mesh_data.meta_data(); }

  stk::mesh::MetaData & getMetaData()
  { return m_mesh_data.meta_data(); }

  /**
   * Get a reference to the bulk data for the stk-mesh.
   */
  const stk::mesh::BulkData & getBulkData() const { return m_mesh_data.bulk_data(); }
        stk::mesh::BulkData & getBulkData()       { return m_mesh_data.bulk_data(); }

 private:
  /**
   * The mesh-data for the mesh. This is a special object that maintains some
   * state between the meta data and bulk data portions of the mesh generation
   * process for use cases.
   */
  stk::io::StkMeshIoBroker m_mesh_data;

  ///> The names of all the side parts
  std::vector<std::string> m_sideset_names;

  ///> Collection of all the side parts
  stk::mesh::PartVector m_sideset_parts;
};

} // namespace simple_fields

}//namespace util
}//namespace io
}//namespace stk

#endif
