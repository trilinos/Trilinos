#ifndef Stk_Mesh_Use_Cases_UseCase_1_hpp
#define Stk_Mesh_Use_Cases_UseCase_1_hpp

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

/** Assume the following mesh of 4 hex8 elements.
 * 
 *  Global node and element numbering
 * <PRE>
 *      3       7      11      15      19         
 *      +-------+-------+-------+-------+        
 *     /       /       /       /       /|       
 *   4/      8/     12/     16/     20/ |       
 *   +-------+-------+-------+-------+  |       
 *   |       |       |       |       |  +18     
 *   |  e1   |  e2   |  e3   |  e4   | /        
 *   |       |       |       |       |/         
 *   +-------+-------+-------+-------+          
 *   1       5      9       13      17          
 * </PRE>
 * 
 *  Local node numbering
 * <PRE>
 *      8       7
 *      +-------+
 *     /       /| 
 *   5/      6/ | 
 *   +-------+  |
 *   |       |  +3 
 *   |  e1   | /
 *   |       |/ 
 *   +-------+
 *   1       2  
 * </PRE>
 */

namespace stk {
namespace mesh {
namespace use_cases {





// Two part MetaData with four entity types:  
// Node, Edge, Face, Element 
// and two parts (partLeft and partRight)
class UseCase_1_MetaData
{
  protected:
    UseCase_1_MetaData(const std::vector<std::string> & entity_type_names);

    stk::mesh::MetaData m_metaData;

  public:

    const stk::mesh::MetaData & metaData() const
    { return m_metaData; }


    stk::mesh::Part & partLeft;

    stk::mesh::Part & partRight;

};


// Mesh from two part Meta Data above plus 
// field_data_chunk_size = 1000
class UseCase_1_Mesh : public UseCase_1_MetaData 
{
  public:
    enum {
      field_data_chunk_size = 1000
    };

    UseCase_1_Mesh( stk::ParallelMachine comm );

    ~UseCase_1_Mesh();

    const stk::mesh::BulkData & bulkData() const
    { return m_bulkData; }

    stk::mesh::BulkData       & modifiableBulkData()
    { return m_bulkData; }

  private:

    stk::mesh::BulkData m_bulkData;

};

std::vector<std::string> get_entity_type_names_1();

void populate( UseCase_1_Mesh & mesh , unsigned nleft , unsigned nright );
bool verifyMesh( const UseCase_1_Mesh & mesh, unsigned nleft, unsigned nright );

// Helper functions for verifyMesh
bool verifyCellTopology( const UseCase_1_Mesh & mesh );
bool verifyEntityCounts( const UseCase_1_Mesh & mesh, unsigned nleft, unsigned nright );
bool verifyRelations( const UseCase_1_Mesh & mesh, unsigned nleft, unsigned nright );

void get_elem_node_ids_1( stk::mesh::EntityId elem_id , stk::mesh::EntityId node_ids[] );

} //namespace use_cases 
} //namespace mesh 
} //namespace stk 

#endif // Stk_Mesh_Use_Cases_UseCase_1_hpp
