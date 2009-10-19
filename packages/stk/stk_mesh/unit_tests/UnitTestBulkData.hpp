
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace stk {
namespace mesh {

class UnitTestBulkData {
public:
  static void testBulkData( ParallelMachine );
  static void testCreateMore_error( ParallelMachine );
  static void testChangeOwner_nodes( ParallelMachine );
  static void testChangeOwner_loop( ParallelMachine );
  static void testChangeOwner_box( ParallelMachine );
  static void testChangeParts( ParallelMachine );
  static void testChangeParts_loop( ParallelMachine );
  static void testDestroy_nodes( ParallelMachine );
  static void testDestroy_loop( ParallelMachine );

  /** Generate simple edge-loop mesh with 'nPerProc' edges
   *  on each processor.  Fill the 'node_ids' and 'edge_ids'
   *  with all node and edge ids.  However, each process
   *  will only have a connected arc of the loop.
   *
   *    edge_ids[ nPerProc * p_rank .. nPerProc * ( p_rank + 1 ) - 1 ]
   */
  static void generate_loop( BulkData & mesh ,
                             const PartVector      & edge_parts , 
                             const bool              generate_aura , 
                             const unsigned          nPerProc ,
                             std::vector<EntityId> & node_ids ,
                             std::vector<EntityId> & edge_ids );

  /** Generate a box mesh which is globally ( ngx X ngy X ngz )
   *  elements where:
   *    ngx = root_box[0][1] - root_box[0][0] ;
   *    ngy = root_box[1][1] - root_box[1][0] ;
   *    ngz = root_box[2][1] - root_box[2][0] ;
   *
   *  The box is partitioned via recursive coordinate bisection
   *  and the resulting local box are given in 'local_box'.
   */
  static void generate_boxes( BulkData  & mesh ,
                              const bool  generate_aura ,
                              const int   root_box[][2] ,
                                    int   local_box[][2] );



  /** Allow unit testing with and without ghosting aura */
  static void modification_end( BulkData & mesh ,bool aura );
};

}
}

