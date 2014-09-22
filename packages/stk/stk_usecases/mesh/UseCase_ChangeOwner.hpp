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

#ifndef TEST_CHANGE_OWNER_HPP
#define TEST_CHANGE_OWNER_HPP

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <stk_mesh/base/CoordinateSystems.hpp>

/*----------------------------------------------------------------------------

A. Generate entire mesh on process #0:

   7-----8-----9
   |     |     |
   |  3  |  4  |
   |     |     |
   4-----5-----6
   |     |     |
   |  1  |  2  |
   |     |     |
   1-----2-----3  ---> X

B. Move Elements { 3 , 4 } and Nodes { 4 , 5 , 6 , 7 , 8 , 9 } to process #1

C. Move Elements { 1 , 2 } and Nodes { 1 , 2 , 3 } to process #1

----------------------------------------------------------------------------*/

typedef stk::mesh::Field<double, stk::mesh::Cartesian> VectorField ;

class Grid2D_Fixture {
public:

  Grid2D_Fixture( stk::ParallelMachine );

  unsigned                    m_spatial_dimension;
  stk::mesh::MetaData         m_fem_meta_data;
  stk::mesh::BulkData         m_bulk_data;
  stk::mesh::Part           & m_quad_part;
  VectorField &               m_coord_field;
  const stk::mesh::EntityRank m_elem_rank;
  const stk::mesh::EntityRank m_node_rank;

  bool test_change_owner( unsigned nx , unsigned ny );

  bool test_change_owner() { return test_change_owner(2,2); }
};

bool test_change_owner_with_constraint( stk::ParallelMachine pm );
bool test_change_owner_2( stk::ParallelMachine pm );
bool test_change_owner_3( stk::ParallelMachine pm );

#endif // TEST_CHANGE_OWNER_HPP

