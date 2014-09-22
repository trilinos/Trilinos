// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <iostream>
#include <set>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/util/string_case_compare.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/TopologyDimensions.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>

namespace stk_examples {

typedef stk::mesh::Field<double>                                ScalarField ;
typedef stk::mesh::Field<int>                                   ScalarIntField ;
typedef stk::mesh::Field<double, stk::mesh::Cartesian>          CartesianField ;
typedef stk::mesh::Field<double, stk::mesh::FullTensor>         FullTensorField ;
typedef stk::mesh::Field<double, stk::mesh::SymmetricTensor>    SymmetricTensorField ;


CartesianField &
declare_vector_field_on_all_nodes(
  stk::mesh::MetaData & meta_data , const std::string & s , unsigned n1 )
{
  return stk::mesh::put_field( meta_data.declare_field<CartesianField>(stk::topology::NODE_RANK, s), meta_data.universal_part() , n1 );
}


CartesianField &
declare_vector_field_on_all_elements(
  stk::mesh::MetaData & meta_data , const std::string & s , unsigned n1 )
{
  return stk::mesh::put_field( meta_data.declare_field<CartesianField>(stk::topology::ELEMENT_RANK, s), meta_data.universal_part() , n1 );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Example reader / writer

void example_io_1( stk::ParallelMachine comm,
                   const std::string& in_filename,
                   const std::string& out_filename)
{
  std::cout
 << "========================================================================\n"
 << " Use Case 1: Simple mesh I/O                                            \n"
 << "========================================================================\n";

  // Open, read, filter meta data from the input mesh file:
  // The coordinates field will be set to the correct dimension.

  stk::io::StkMeshIoBroker mesh_data(comm);
  mesh_data.add_mesh_database(in_filename, stk::io::READ_MESH);
  mesh_data.create_input_mesh();
  stk::mesh::MetaData &meta_data = mesh_data.meta_data();

  //----------------------------------
  // Print the parts that were read from the file:

  const stk::mesh::PartVector & io_parts = meta_data.get_parts();
  for (stk::mesh::PartVector::const_iterator i = io_parts.begin();
           i != io_parts.end(); ++i) {
    stk::mesh::Part & part = **i ;

    const CellTopologyData * const top = meta_data.get_cell_topology( part ).getCellTopologyData();

    const unsigned num_subsets = part.subsets().size();

    std::cout << "  Part: " << part.name();
    if ( part.primary_entity_rank() < meta_data.entity_rank_count() ) {
      std::cout << " "
                << meta_data.entity_rank_name( part.primary_entity_rank() );
    }
    if ( top ) { std::cout << " " << top->name ; }
    if ( num_subsets ) { std::cout << " with " << num_subsets << " subsets" ; }
    std::cout << std::endl ;
  }

  //----------------------------------
  // Now I can process what was read from the mesh file to
  // declare more parts and fields.

  meta_data.commit();

  // Read the model (topology, coordinates, attributes, etc)
  // from the mesh-file into the mesh bulk data.
  mesh_data.populate_bulk_data();

  //----------------------------------
  // Create a mesh writer that will simply write out what was read.
  // the parts, attributes, and transient arguments can be different
  // that what was read.
  size_t resultFileIndex = mesh_data.create_output_mesh(out_filename, stk::io::WRITE_RESULTS);
  mesh_data.write_output_mesh(resultFileIndex);
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace stk_examples


