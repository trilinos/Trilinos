/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/fem/EntityTypes.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace stk {
namespace unit_test {

// Helper function to generate a vector of entity names
std::vector<std::string>  get_entity_type_names ( int rank );

// MetaData wrapper class used to instantiate all the parts as
// references to help with testing of stk::mesh.
class  UnitTestMetaData
{
  protected:
    stk::mesh::MetaData   m_meta_data;
    stk::mesh::Part      &m_test_part;   // A simple part
    stk::mesh::Part      &m_cell_part;   // A part to put cells in

    stk::mesh::Part      &m_part_A_0;
    stk::mesh::Part      &m_part_A_1;
    stk::mesh::Part      &m_part_A_2;
    stk::mesh::Part      &m_part_A_3;

    stk::mesh::Part      &m_part_A_superset;

    stk::mesh::Part      &m_part_B_0;
    stk::mesh::Part      &m_part_B_1;
    stk::mesh::Part      &m_part_B_2;
    stk::mesh::Part      &m_part_B_3;

    stk::mesh::Part      &m_part_B_superset;

    UnitTestMetaData ( const std::vector<std::string> &types );
    UnitTestMetaData ();
   ~UnitTestMetaData () {}

  public:
     enum { MAX_RANK = 3 };
     const stk::mesh::MetaData & meta_data () const { return m_meta_data; }

     stk::mesh::Part &     get_test_part () { return m_test_part; }
     stk::mesh::Part &     get_cell_part () { return m_cell_part; }

     stk::mesh::Part &     get_part_a_0 () { return m_part_A_0; }
     stk::mesh::Part &     get_part_a_1 () { return m_part_A_1; }
     stk::mesh::Part &     get_part_a_2 () { return m_part_A_2; }
     stk::mesh::Part &     get_part_a_3 () { return m_part_A_3; }

     stk::mesh::Part &     get_part_a_superset () { return m_part_A_superset; }

     stk::mesh::Part &     get_part_b_0 () { return m_part_B_0; }
     stk::mesh::Part &     get_part_b_1 () { return m_part_B_1; }
     stk::mesh::Part &     get_part_b_2 () { return m_part_B_2; }
     stk::mesh::Part &     get_part_b_3 () { return m_part_B_3; }

     stk::mesh::Part &     get_part_b_superset () { return m_part_B_superset; }
}
;


class  UnitTestMesh : public UnitTestMetaData
{
  protected:
    stk::mesh::BulkData  m_bulk_data;

    unsigned    m_comm_rank;
    unsigned    m_comm_size;

    stk::mesh::BulkData::BulkDataSyncState     m_previous_state;


  /** Generate a box mesh which is globally ( ngx X ngy X ngz )
   *  elements where:
   *    ngx = root_box[0][1] - root_box[0][0] ;
   *    ngy = root_box[1][1] - root_box[1][0] ;
   *    ngz = root_box[2][1] - root_box[2][0] ;
   *
   *  The box is partitioned via recursive coordinate bisection
   *  and the resulting local box are given in 'local_box'.
   */
     void priv_generate_boxes( stk::mesh::BulkData  & mesh ,
                              const bool  generate_aura ,
                              const int   root_box[][2] ,
                                    int   local_box[][2] );

  /** Generate simple edge-loop mesh with 'nPerProc' edges
   *  on each processor.  Fill the 'node_ids' and 'edge_ids'
   *  with all node and edge ids.  However, each process
   *  will only have a connected arc of the loop.
   *
   *    edge_ids[ nPerProc * p_rank .. nPerProc * ( p_rank + 1 ) - 1 ]
  void priv_generate_loop( stk::mesh::BulkData & mesh ,
                        const PartVector      & edge_parts , 
                        const bool              generate_aura , 
                        const unsigned          nPerProc ,
                        std::vector<EntityId> & node_ids ,
                        std::vector<EntityId> & edge_ids );
   */

     void enter_modification () 
    {
      m_previous_state = m_bulk_data.synchronized_state();
      if ( m_previous_state == stk::mesh::BulkData::SYNCHRONIZED )
        m_bulk_data.modification_begin();
    }

     void exit_modification ()
    {
      if ( m_previous_state == stk::mesh::BulkData::SYNCHRONIZED )
        m_bulk_data.modification_end();
    }

  public:
    UnitTestMesh ( stk::ParallelMachine pm = MPI_COMM_WORLD , unsigned block_size = 1000 );
   ~UnitTestMesh () {}

     const stk::mesh::BulkData & bulk_data () const { return m_bulk_data; } 
     stk::mesh::BulkData       & nonconst_bulk_data () { return m_bulk_data; }

     unsigned  comm_size() const { return m_comm_size; }
     unsigned  comm_rank() const { return m_comm_rank; }

     void  generate_boxes ( bool aura = false );

     stk::mesh::Entity  &get_new_entity ( stk::mesh::EntityType rank , stk::mesh::EntityId parallel_dependent_id )
     {
       return m_bulk_data.declare_entity ( rank , parallel_dependent_id*m_comm_size + m_comm_rank + 1 , std::vector<stk::mesh::Part *> () );
     }
}
;


} // namespace unit_test
} // namespace stk

