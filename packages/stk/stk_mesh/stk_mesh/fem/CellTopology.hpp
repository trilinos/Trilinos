#ifndef stk_mesh_fem_CellTopology_hpp
#define stk_mesh_fem_CellTopology_hpp

#ifdef HAVE_SHARDS_DEBUG
#define STK_MESH_FEM_CHECK_REQUIRE( S )  S
#else
#define STK_MESH_FEM_CHECK_REQUIRE( S ) /* empty */
#endif

#include <Shards_CellTopologyTraits.hpp>
#include <Shards_CellTopology.hpp>
#include <Shards_CellTopologyData.h>
#include <Shards_BasicTopologies.hpp>

namespace stk {
namespace mesh {
namespace fem {

/** \addtogroup  shards_package_cell_topology
 *  \{
 */

/*------------------------------------------------------------------------*/

class CellTopology ;

/** \brief Overloaded << operator for CellTopologyData objects. */
//std::ostream & operator << ( std::ostream & , const CellTopologyData & );


/** \brief Overloaded << operator for CellTopology objects.  */
std::ostream & operator << ( std::ostream & , const CellTopology & );


/*------------------------------------------------------------------------*/
/** \class CellTopology
 *  \brief Provide input checked access (in debug mode) to
 *         \ref shards::CellTopologyData "cell topology data"
 *         and a procedure to create custom cell topologies. 
 *
 *  Input checking is compiled in when the HAVE_SHARDS_DEBUG macro is defined.
 *
 *  \author Created by P. Bochev, H. C. Edwards and D. Ridzal.
 *
 *  \nosubgrouping
 *
 *  Two kinds of CellTopology objects are used:
 *  (1) wrappers for the predefined basic cell topologies and
 *  (2) temporary custom cell topologies.
 *
 *  A temporary custom cell topology is created by a calling code,
 *  passed to a computational kernel such as those supplied by Intrepid,
 *  and then discarded.  The use case is discretization function
 *  evaluations on a arbitrary polyhedon mesh does not have
 *  standard or consistent cell topologies.
 */
class CellTopology {
public:
  /** \brief  Wrapper for safe access to a raw cell topology data.
   *  \param  cell             [in]  - pointer to raw cell topology
   *          Example use: the statements
   *          \code
   *          CellTopology  myTri_3 ( getCellTopologyData< Triangle<3> >() )
   *          CellTopology  myTet_10( getCellTopologyData< Tetrahedron<10>() );
   *          \endcode
   *          wraps Triangle<3> and Tetrahedron<10> (defined in stk_mesh_fem_BasicTopologies)
   */
  CellTopology( const CellTopologyData * cell )
    : m_cell( cell )
  {}
  
        
  /** \brief Assignment operator <var>*this = right</var>.  */
  CellTopology& operator = (const CellTopology& right) {
    m_cell = right.m_cell;
    return *this;
  }

  /** \brief Copy constructor */
  CellTopology( const CellTopology& right )
    : m_cell(right.m_cell)
  {}

  /** \brief  Default constructor initializes to NULL */
  CellTopology()
    : m_cell(0)
  {}

  /** \brief Destructor */
  ~CellTopology()
  {}
  

  /** \} */

  operator bool() {
    return m_cell != 0;
  }

  /*------------------------------------------------------------------*/
  /** \name  Safe query methods
   *  \{
   */
  
          
  /** \brief  Dimension of this cell topology */
  unsigned getDimension() const
  {
    STK_MESH_FEM_CHECK_REQUIRE( requireCell() );
    return m_cell->dimension ;
  }
  
        
  /** \brief  Unique key for this cell topology;
   *          under certain subcell uniformity conditions.
   *
   *  The key is only guaranteed to be unique for predefined cell topologies.
   */
  unsigned getKey() const
  {
    STK_MESH_FEM_CHECK_REQUIRE( requireCell() );
    return m_cell->key ;
  }

        
  /** \brief  Unique key for this cell's base topology;
   *          under certain subcell uniformity conditions.
   *
   *  The key is only guaranteed to be unique for predefined cell topologies.
   */
  unsigned getBaseKey() const
  {
    STK_MESH_FEM_CHECK_REQUIRE( requireCell() );
    return m_cell->base->key ;
  }
        
        
        
  /** \brief  Unique name for this cell topology;
   *          
   *  The name is only guaranteed to be unique for predefined cell topologies.
   *  A calling code may construct custom cell topologies with arbitrary names.
   */
  const char* getName() const
  {
    STK_MESH_FEM_CHECK_REQUIRE( requireCell() );
    return m_cell->name ;
  }

        
  /** \brief  Unique name for this cell's base topology.
   *          
   */
  const char* getBaseName() const
  {
    STK_MESH_FEM_CHECK_REQUIRE( requireCell() );
    return m_cell->base->name ;
  }
        
        
  /** \brief  Node count of this cell topology */
  unsigned getNodeCount() const
  {
    STK_MESH_FEM_CHECK_REQUIRE( requireCell() );
    return m_cell->node_count ;
  }
  
        
  /** \brief  Vertex count of this cell topology */
  unsigned getVertexCount() const
  {
    STK_MESH_FEM_CHECK_REQUIRE( requireCell() );
    return m_cell->vertex_count ;
  }
  
        
  /** \brief  Edge boundary subcell count of this cell topology */
  unsigned getEdgeCount() const
  {
    STK_MESH_FEM_CHECK_REQUIRE( requireCell() );
    return m_cell->edge_count ;
  }
  
  /** \brief  Face boundary subcell count of this cell topology */
  unsigned getFaceCount() const
  {
    STK_MESH_FEM_CHECK_REQUIRE( requireCell() );
    return m_cell->dimension == 3 ? m_cell->side_count : 0 ;
  }
  
        
  /** \brief  Side boundary subcell count of this cell topology */
  unsigned getSideCount() const
  {
    STK_MESH_FEM_CHECK_REQUIRE( requireCell() );
    return m_cell->side_count ;
  }
  
        
  /** \brief  This cell's raw topology data */
  const CellTopologyData * getTopologyData() const {
    return m_cell ;
  }

        
  /** \brief  This cell's base cell topology's raw topology data */
  const CellTopology getBaseTopology() const
  {
    STK_MESH_FEM_CHECK_REQUIRE( requireCell() );
    return CellTopology(m_cell->base);
  }

        
  /** \brief  Raw cell topology data for a subcell of the
   *          given dimension and ordinal.
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  const CellTopology getTopology( const unsigned subcell_dim , const unsigned subcell_ord ) const
  {
    STK_MESH_FEM_CHECK_REQUIRE( requireCell() );
    STK_MESH_FEM_CHECK_REQUIRE( requireDimension(subcell_dim) );
    STK_MESH_FEM_CHECK_REQUIRE( requireSubcell(subcell_dim,subcell_ord) );
    return CellTopology(m_cell->subcell[subcell_dim][subcell_ord].topology) ;
  }

        
  /** \brief  Raw cell topology data for the base topology of a subcell of the
   *          given dimension and ordinal.
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  const CellTopology getBaseTopology( const unsigned subcell_dim ,
                                            const unsigned subcell_ord ) const
  {
    return getTopology(subcell_dim,subcell_ord).getBaseTopology();
  }

        
  /** \brief  Key of a subcell of the given dimension and ordinal. 
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  unsigned getKey( const unsigned subcell_dim ,
                   const unsigned subcell_ord ) const
  {
    return getTopology(subcell_dim,subcell_ord).getKey() ;
  }


        
  /** \brief  Name of a subcell of the given dimension and ordinal. 
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  const char * getName(const unsigned subcell_dim,   
                       const unsigned subcell_ord) const
  {
    return getTopology(subcell_dim,subcell_ord).getName();
  }
        
        
  /** \brief  Node count of a subcell of the given dimension and ordinal. 
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  unsigned getNodeCount( const unsigned subcell_dim ,
                         const unsigned subcell_ord ) const
  {
    return getTopology(subcell_dim,subcell_ord).getNodeCount() ;
  }

        
  /** \brief  Vertex count of a subcell of the given dimension and ordinal. 
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  unsigned getVertexCount( const unsigned subcell_dim ,
                           const unsigned subcell_ord ) const
  {
    return getTopology(subcell_dim,subcell_ord).getVertexCount() ;
  }

        
  /** \brief  Edge count of a subcell of the given dimension and ordinal. 
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  unsigned getEdgeCount( const unsigned subcell_dim ,
                         const unsigned subcell_ord ) const
  {
    return getTopology(subcell_dim,subcell_ord).getEdgeCount() ;
  }
  
        
  /** \brief  Side count of a subcell of the given dimension and ordinal. 
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */        
  unsigned getSideCount( const unsigned subcell_dim ,
                         const unsigned subcell_ord ) const
  {
    return getTopology(subcell_dim,subcell_ord).getSideCount() ;
  }

        
  /** \brief  Subcell count of subcells of the given dimension. 
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   */        
  unsigned getSubcellCount( const unsigned subcell_dim ) const
  {
    STK_MESH_FEM_CHECK_REQUIRE( requireCell() );
    STK_MESH_FEM_CHECK_REQUIRE( requireDimension(subcell_dim) );
    return m_cell->subcell_count[subcell_dim] ;
  }
  
        
  /** \brief  Query if all subcells of the given dimension
   *          have the same cell topology.
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   */
  bool getSubcellHomogeneity( const unsigned subcell_dim ) const
  {
    STK_MESH_FEM_CHECK_REQUIRE( requireCell() );
    STK_MESH_FEM_CHECK_REQUIRE( requireDimension(subcell_dim) );
    return 0 != m_cell->subcell_homogeneity[subcell_dim] ;
  }
  
        
  /** \brief  Mapping from a subcell's node ordinal to a
   *          node ordinal of this parent cell topology.
   *  \param  subcell_dim      [in]  - spatial dimension of the subcell
   *  \param  subcell_ord      [in]  - subcell ordinal
   *  \param  subcell_node_ord [in]  - node ordinal relative to subcell
   */
  unsigned getNodeMap( const unsigned  subcell_dim ,
                       const unsigned  subcell_ord ,
                       const unsigned  subcell_node_ord ) const
  {
    STK_MESH_FEM_CHECK_REQUIRE( requireCell() );
    STK_MESH_FEM_CHECK_REQUIRE( requireDimension(subcell_dim) );
    STK_MESH_FEM_CHECK_REQUIRE( requireSubcell(subcell_dim,subcell_ord) );
    STK_MESH_FEM_CHECK_REQUIRE( requireNodeMap(subcell_dim,subcell_ord,subcell_node_ord));
    return m_cell->subcell[subcell_dim][subcell_ord].node[subcell_node_ord];
  }


  /** \brief  Number of node permutations defined for this cell */
  unsigned getNodePermutationCount() const
  {
    STK_MESH_FEM_CHECK_REQUIRE(requireCell());
    return m_cell->permutation_count ;
  }

  /** \brief  Permutation of a cell's node ordinals.
   *  \param  permutation_ordinal [in]
   *  \param  node_ordinal        [in]
   */
  unsigned getNodePermutation( const unsigned permutation_ord ,
                               const unsigned node_ord ) const
  {
    STK_MESH_FEM_CHECK_REQUIRE(requireCell());
    STK_MESH_FEM_CHECK_REQUIRE(requireNodePermutation(permutation_ord,node_ord));
    return m_cell->permutation[permutation_ord].node[node_ord];
  }
  
  /** \brief  Permutation of a cell's node ordinals.
   *  \param  permutation_ordinal [in]
   *  \param  node_ordinal        [in]
   */
  unsigned getNodePermutationPolarity( const unsigned permutation_ord ) const
  {
    STK_MESH_FEM_CHECK_REQUIRE(requireCell());
    STK_MESH_FEM_CHECK_REQUIRE(requireNodePermutation(permutation_ord,0));
    return m_cell->permutation[permutation_ord].polarity;
  }
  
  /** \brief  Inverse permutation of a cell's node ordinals.
   *  \param  permutation_ordinal [in]
   *  \param  node_ordinal        [in]
   */
  unsigned getNodePermutationInverse( const unsigned permutation_ord ,
                                      const unsigned node_ord ) const
  {
    STK_MESH_FEM_CHECK_REQUIRE(requireCell());
    STK_MESH_FEM_CHECK_REQUIRE(requireNodePermutation(permutation_ord,node_ord));
    return m_cell->permutation_inverse[permutation_ord].node[node_ord];
  }
  
  /** \} */
  
  /*------------------------------------------------------------------*/
  /** \name Constructors and destructor
   *  \{
   */

private:
  
  /** \brief Throws runtime_error if CellTopology object is null or hase null 
      base topology   
  */
  void requireCell() const ;
  
  
  /** \brief  Throws invalid_argument if subcell dimension exceedes the maximal 
   *          admissible space dimension 3.
   *  \param  subcell_dim    [in]  - spatial dimension of a subcell
   */
  void requireDimension( const unsigned subcellDim ) const ;
  
  
  /** \brief  Throws invalid_argument if subcell_ord exceeds the actual number
   *          of subcells with specified dimension.
   *  \param  subcell_dim    [in]  - spatial dimension of a subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  void requireSubcell( const unsigned subcellDim ,
                       const unsigned subcellOrd ) const ;
  
  
  /** \brief  Throws invalid_argument if node_ord exceeds the actual number 
   *          of nodes in the subcell with specified dimension and ordinal.
   *  \param  subcell_dim    [in]  - spatial dimension of a subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   *  \param  node_ord       [in]  - node ordinal
   */
  void requireNodeMap( const unsigned subcellDim ,
                       const unsigned subcellOrd ,
                       const unsigned nodeOrd ) const ;

  void requireNodePermutation( const unsigned permutationOrd ,
                               const unsigned nodeOrd ) const ;

private:
  const CellTopologyData    * m_cell ;

}; // class CellTopology


template< typename id_type >
int findPermutation( const CellTopology top ,
                     const id_type * const expected_node ,
                     const id_type * const actual_node )
{
  return shards::findPermutation( *top.getTopologyData() , expected_node , actual_node );
}

inline
bool operator==(const CellTopology &left, const CellTopology &right) 
{
  return left.getTopologyData() == right.getTopologyData();
}

inline
bool operator<(const CellTopology &left, const CellTopology &right) 
{
  return left.getTopologyData() < right.getTopologyData();
}

inline
bool operator!=(const CellTopology &left, const CellTopology &right) {
  return !(left == right);
}

/** \} */

} // namespace fem
} // namespace mesh
} // namespace stk

#undef STK_MESH_FEM_CHECK_REQUIRE

#endif // stk_mesh_fem_CellTopology_hpp

