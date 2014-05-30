#ifndef TEST_CHANGE_OWNER_HPP
#define TEST_CHANGE_OWNER_HPP

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

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

typedef stk_classic::mesh::Field<double, stk_classic::mesh::Cartesian> VectorField ;

class Grid2D_Fixture {
public:

  Grid2D_Fixture( stk_classic::ParallelMachine );

  unsigned                    m_spatial_dimension;
  stk_classic::mesh::fem::FEMMetaData m_fem_meta_data;
  stk_classic::mesh::BulkData         m_bulk_data;
  stk_classic::mesh::Part           & m_quad_part;
  VectorField &               m_coord_field;
  const stk_classic::mesh::EntityRank m_elem_rank;
  const stk_classic::mesh::EntityRank m_node_rank;

  bool test_change_owner( unsigned nx , unsigned ny );

  bool test_change_owner() { return test_change_owner(2,2); }
};

bool test_change_owner_with_constraint( stk_classic::ParallelMachine pm );
bool test_change_owner_2( stk_classic::ParallelMachine pm );
bool test_change_owner_3( stk_classic::ParallelMachine pm );

#endif // TEST_CHANGE_OWNER_HPP

