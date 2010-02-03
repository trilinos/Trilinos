#ifndef Stk_Mesh_Use_Cases_UseCase_2_hpp
#define Stk_Mesh_Use_Cases_UseCase_2_hpp

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/fem/FieldTraits.hpp>

/** stk_mesh Use Case 2
 * Note, the basic layout of the mesh is the same as Use Case 1, with
 * the addition of fields.
 *
 * Assume the following mesh of 4 hex8 elements.
 * 
 *  Global node and element numbering
 * <PRE>
 *      3       7      11      15      19         
 *      +-------+-------+-------+-------+        
 *     /       /       /       /       /|       
 *   4/      8/     12/     16/     20/ |       
 *   +-------+-------+-------+-------+  |       
 *   |       |       |       |       |  +18        Z  Y
 *   |  e1   |  e2   |  e3   |  e4   | /           | /
 *   |       |       |       |       |/            |/
 *   +-------+-------+-------+-------+             *--X
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



typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType ;
typedef stk::mesh::Field<double>                      ScalarFieldType ;

// Two part MetaData with four entity types:  
// Node, Edge, Face, Element 
// and two parts (partLeft and partRight)
// and three fields (coordinates, temperature, and volume)
class UseCase_2_MetaData
{
  public:

    const stk::mesh::MetaData & metaData() const
    { return m_metaData; }

    stk::mesh::Part & partLeft() const
    { return m_partLeft; }

    stk::mesh::Part & partRight() const
    { return m_partRight; }

    VectorFieldType & coordinates_field()
    { return m_coordinates_field; }

    const VectorFieldType & const_coordinates_field() const
    { return m_coordinates_field; }

    ScalarFieldType & temperature_field()
    { return m_temperature_field; }

    const ScalarFieldType & const_temperature_field() const
    { return m_temperature_field; }

    ScalarFieldType & volume_field()
    { return m_volume_field; }

    const ScalarFieldType & const_volume_field() const
    { return m_volume_field; }

  protected:
    UseCase_2_MetaData(const std::vector<std::string> & entity_type_names);

    stk::mesh::MetaData m_metaData;

    stk::mesh::Part & m_partLeft;

    stk::mesh::Part & m_partRight;

    VectorFieldType & m_coordinates_field;
    ScalarFieldType & m_temperature_field;
    ScalarFieldType & m_volume_field;

};


// Mesh from two part Meta Data above plus 
// field_data_chunk_size = 1000
class UseCase_2_Mesh : public UseCase_2_MetaData 
{
  public:
    enum {
      field_data_chunk_size = 1000
    };

    UseCase_2_Mesh( stk::ParallelMachine comm );

    ~UseCase_2_Mesh();

    const stk::mesh::BulkData & bulkData() const
    { return m_bulkData; }

    stk::mesh::BulkData       & modifiableBulkData()
    { return m_bulkData; }

  private:

    stk::mesh::BulkData m_bulkData;

};

std::vector<std::string> get_entity_type_names_2();

void populate( UseCase_2_Mesh & mesh , unsigned nleft , unsigned nright );
bool verifyMesh( const UseCase_2_Mesh & mesh, unsigned nleft, unsigned nright );

// Helper functions for verifyMesh
bool verifyCellTopology( const UseCase_2_Mesh & mesh );
bool verifyEntityCounts( const UseCase_2_Mesh & mesh, unsigned nleft, unsigned nright );
bool verifyRelations( const UseCase_2_Mesh & mesh, unsigned nleft, unsigned nright );
bool verifyFields( const UseCase_2_Mesh & mesh );

void get_elem_node_ids_2( stk::mesh::EntityId elem_id , stk::mesh::EntityId node_ids[] );

} //namespace use_cases 
} //namespace mesh 
} //namespace stk 

#endif // Stk_Mesh_Use_Cases_UseCase_2_hpp
