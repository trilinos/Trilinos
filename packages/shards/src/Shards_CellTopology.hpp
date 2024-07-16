// @HEADER
// *****************************************************************************
//                Shards : Shared Discretization Tools
//
// Copyright 2008-2011 NTESS and the Shards contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef Shards_CellTopology_hpp
#define Shards_CellTopology_hpp

#ifdef HAVE_SHARDS_DEBUG
#define SHARDS_REQUIRE( S )  S
#else
#define SHARDS_REQUIRE( S ) /* empty */
#endif

#include <string>
#include <vector>
#include <Shards_CellTopologyData.h>
#include <Shards_BasicTopologies.hpp>

namespace shards {

/** \addtogroup  shards_package_cell_topology
 *  \{
 */

/*------------------------------------------------------------------------*/

class CellTopology ;

/** \brief Overloaded << operator for CellTopologyData objects. */
//std::ostream & operator << ( std::ostream & , const CellTopologyData & );


/** \brief Overloaded << operator for CellTopology objects.  */
std::ostream & operator << ( std::ostream & , const CellTopology & );


/** \enum  Shards::ECellType
    \brief Enumeration of cell types in Shards 
  */
enum ECellType {
  ALL_CELLS = 0,
  STANDARD_CELL,
  NONSTANDARD_CELL
};

inline std::string ECellTypeToString(ECellType cellType) {
  std::string retString;
  switch(cellType){
    case ALL_CELLS:         retString = "All";           break;
    case STANDARD_CELL:     retString = "Standard";      break;
    case NONSTANDARD_CELL:  retString = "Nonstandard";   break;
    default:                retString = "Invalid Cell";
  }
  return retString;
}


/** \enum  Shards::ETopologyType
    \brief Enumeration of topology types in Shards 
  */
enum ETopologyType {
  ALL_TOPOLOGIES,
  BASE_TOPOLOGY,
  EXTENDED_TOPOLOGY
};

inline std::string ETopologyTypeToString(ETopologyType topologyType) {
  std::string retString;
  switch(topologyType){
    case ALL_TOPOLOGIES:      retString = "All";            break;
    case BASE_TOPOLOGY:       retString = "Base";           break;
    case EXTENDED_TOPOLOGY:   retString = "Extended";       break;
    default:                  retString = "Invalid Topology";
  }
  return retString;
}


/** \brief  Returns an std::vector with all cell topologies that meet the specified
            selection flags.
     
    \param  topologies      [out]   - vector with all topologies that 
    \param  cellDim         [in]    - cell dimension; 0, 1, 2, 3, or 4 (default = all dimensions)
    \param  cellType        [in]    - cell type: default = ALL_CELLS  
    \param  topologyType    [in]    - topology type: default = ALL_TOPOLOGIES
  */
void getTopologies(std::vector<shards::CellTopology>& topologies,
                   const unsigned       cellDim   = 4,
                   const ECellType      cellType      = ALL_CELLS,
                   const ETopologyType  topologyType  = ALL_TOPOLOGIES);



/** \brief  Checks if the cell topology is predefined in shards

    \param  cell            [in]  - cell topology
    \return 1   if the cell topology is defined in shards, 
            0   if it is a custom, user-defined cell-topology
*/
int isPredefinedCell(const CellTopology &  cell);



/*------------------------------------------------------------------------*/
/** \class shards::CellTopology
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
  
  const CellTopologyData    * m_cell ;

public:

  /*------------------------------------------------------------------*/
  /** \name  Safe query methods
   *  \{
   */
  
          
  /** \brief  Dimension of this cell topology */
  unsigned getDimension() const
    {
      SHARDS_REQUIRE( requireCell() );
      return m_cell->dimension ;
    }
  
        
  /** \brief  Unique key for this cell topology;
   *          under certain subcell uniformity conditions.
   *
   *  The key is only guaranteed to be unique for predefined cell topologies.
   */
  unsigned getKey() const
    {
      SHARDS_REQUIRE( requireCell() );
      return m_cell->key ;
    }

        
  /** \brief  Unique key for this cell's base topology;
   *          under certain subcell uniformity conditions.
   *
   *  The key is only guaranteed to be unique for predefined cell topologies.
   */
  unsigned getBaseKey() const
    {
      SHARDS_REQUIRE( requireCell() );
      return m_cell->base->key ;
    }
        
        
        
  /** \brief  Unique name for this cell topology;
   *          
   *  The name is only guaranteed to be unique for predefined cell topologies.
   *  A calling code may construct custom cell topologies with arbitrary names.
   */
  const char* getName() const
    {
        SHARDS_REQUIRE( requireCell() );
        return m_cell->name ;
    }

        
  /** \brief  Unique name for this cell's base topology.
   *          
   */
  const char* getBaseName() const
    {
      SHARDS_REQUIRE( requireCell() );
      return m_cell->base->name ;
    }
        
        
  /** \brief  Node count of this cell topology */
  unsigned getNodeCount() const
    {
      SHARDS_REQUIRE( requireCell() );
      return m_cell->node_count ;
    }
  
        
  /** \brief  Vertex count of this cell topology */
  unsigned getVertexCount() const
    {
      SHARDS_REQUIRE( requireCell() );
      return m_cell->vertex_count ;
    }
  
        
  /** \brief  Edge boundary subcell count of this cell topology */
  unsigned getEdgeCount() const
    {
      SHARDS_REQUIRE( requireCell() );
      return m_cell->edge_count ;
    }
  
  /** \brief  Face boundary subcell count of this cell topology */
  unsigned getFaceCount() const
    {
      SHARDS_REQUIRE( requireCell() );
      return m_cell->dimension == 3 ? m_cell->side_count : 0 ;
    }
  
        
  /** \brief  Side boundary subcell count of this cell topology */
  unsigned getSideCount() const
    {
      SHARDS_REQUIRE( requireCell() );
      return m_cell->side_count ;
    }
  
        
  /** \brief  This cell's raw topology data */
  bool isValid() const
    { return m_cell != 0 ; }

        
  /** \brief  This cell's raw topology data */
  const CellTopologyData * getCellTopologyData() const
    { return m_cell ; }

        
  /** \brief  This cell's base cell topology's raw topology data */
  const CellTopologyData * getBaseCellTopologyData() const
    {
      SHARDS_REQUIRE( requireCell() );
      return m_cell->base ;
    }

        
  /** \brief  Raw cell topology data for a subcell of the
   *          given dimension and ordinal.
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  const CellTopologyData * getCellTopologyData( const unsigned subcell_dim ,
                                                const unsigned subcell_ord ) const
    {
      SHARDS_REQUIRE( requireCell() );
      SHARDS_REQUIRE( requireDimension(subcell_dim) );
      SHARDS_REQUIRE( requireSubcell(subcell_dim,subcell_ord) );
      return m_cell->subcell[subcell_dim][subcell_ord].topology ;
    }

        
  /** \brief  Raw cell topology data for the base topology of a subcell of the
   *          given dimension and ordinal.
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  const CellTopologyData * getBaseCellTopologyData( const unsigned subcell_dim ,
                                                    const unsigned subcell_ord ) const
    {
      return getCellTopologyData(subcell_dim,subcell_ord)->base ;
    }

        
  /** \brief  Key of a subcell of the given dimension and ordinal. 
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  unsigned getKey( const unsigned subcell_dim ,
                   const unsigned subcell_ord ) const
    {
      return getCellTopologyData(subcell_dim,subcell_ord)->key ;
    }


        
  /** \brief  Name of a subcell of the given dimension and ordinal. 
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  const char * getName(const unsigned subcell_dim,   
                       const unsigned subcell_ord) const
    {
      return getCellTopologyData(subcell_dim,subcell_ord) -> name;
    }
        
        
  /** \brief  Node count of a subcell of the given dimension and ordinal. 
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  unsigned getNodeCount( const unsigned subcell_dim ,
                         const unsigned subcell_ord ) const
    {
      return getCellTopologyData(subcell_dim,subcell_ord)->node_count ;
    }

        
  /** \brief  Vertex count of a subcell of the given dimension and ordinal. 
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  unsigned getVertexCount( const unsigned subcell_dim ,
                           const unsigned subcell_ord ) const
    {
      return getCellTopologyData(subcell_dim,subcell_ord)->vertex_count ;
    }

        
  /** \brief  Edge count of a subcell of the given dimension and ordinal. 
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  unsigned getEdgeCount( const unsigned subcell_dim ,
                         const unsigned subcell_ord ) const
    {
      return getCellTopologyData(subcell_dim,subcell_ord)->edge_count ;
    }
  
        
  /** \brief  Side count of a subcell of the given dimension and ordinal. 
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */        
  unsigned getSideCount( const unsigned subcell_dim ,
                         const unsigned subcell_ord ) const
    {
      return getCellTopologyData(subcell_dim,subcell_ord)->side_count ;
    }

        
  /** \brief  Subcell count of subcells of the given dimension. 
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   */        
  unsigned getSubcellCount( const unsigned subcell_dim ) const
    {
      SHARDS_REQUIRE( requireCell() );
      SHARDS_REQUIRE( requireDimension(subcell_dim) );
      return m_cell->subcell_count[subcell_dim] ;
    }
  
        
  /** \brief  Query if all subcells of the given dimension
   *          have the same cell topology.
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   */
  bool getSubcellHomogeneity( const unsigned subcell_dim ) const
    {
      SHARDS_REQUIRE( requireCell() );
      SHARDS_REQUIRE( requireDimension(subcell_dim) );
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
      SHARDS_REQUIRE( requireCell() );
      SHARDS_REQUIRE( requireDimension(subcell_dim) );
      SHARDS_REQUIRE( requireSubcell(subcell_dim,subcell_ord) );
      SHARDS_REQUIRE( requireNodeMap(subcell_dim,subcell_ord,subcell_node_ord));
      return m_cell->subcell[subcell_dim][subcell_ord].node[subcell_node_ord];
    }


  /** \brief  Number of node permutations defined for this cell */
  unsigned getNodePermutationCount() const
    {
      SHARDS_REQUIRE(requireCell());
      return m_cell->permutation_count ;
    }

  /** \brief  Permutation of a cell's node ordinals.
   *  \param  permutation_ordinal [in]
   *  \param  node_ordinal        [in]
   */
  unsigned getNodePermutation( const unsigned permutation_ord ,
                               const unsigned node_ord ) const
    {
      SHARDS_REQUIRE(requireCell());
      SHARDS_REQUIRE(requireNodePermutation(permutation_ord,node_ord));
      return m_cell->permutation[permutation_ord].node[node_ord];
    }
  
  /** \brief  Permutation of a cell's node ordinals.
   *  \param  permutation_ordinal [in]
   *  \param  node_ordinal        [in]
   */
  unsigned getNodePermutationPolarity( const unsigned permutation_ord ) const
    {
      SHARDS_REQUIRE(requireCell());
      SHARDS_REQUIRE(requireNodePermutation(permutation_ord,0));
      return m_cell->permutation[permutation_ord].polarity;
    }
  
  /** \brief  Inverse permutation of a cell's node ordinals.
   *  \param  permutation_ordinal [in]
   *  \param  node_ordinal        [in]
   */
  unsigned getNodePermutationInverse( const unsigned permutation_ord ,
                                      const unsigned node_ord ) const
    {
      SHARDS_REQUIRE(requireCell());
      SHARDS_REQUIRE(requireNodePermutation(permutation_ord,node_ord));
      return m_cell->permutation_inverse[permutation_ord].node[node_ord];
    }
  
  /** \} */
  
  /*------------------------------------------------------------------*/
  /** \name Constructors and destructor
   *  \{
   */

  /** \brief  Wrapper for safe access to a raw cell topology data.
   *  \param  cell             [in]  - pointer to raw cell topology
   *          Example use: the statements
   *          \code
   *          CellTopology  myTri_3 ( getCellTopologyData< Triangle<3> >() )
   *          CellTopology  myTet_10( getCellTopologyData< Tetrahedron<10>() );
   *          \endcode
   *          wraps Triangle<3> and Tetrahedron<10> (defined in Shards_BasicTopologies)
   */
  CellTopology( const CellTopologyData * cell )
    : m_cell( cell )
  {}
  
        
  /** \brief  Constructs custom 1-cell (line) with base topology Line<>. 
   *          Example use: the statement
   *          \code
   *          CellTopology  customLine("customLine", 4);    
   *          \endcode
   *           defines custom line with 4 nodes.
   */
  CellTopology( const std::string & name,
                const unsigned      nodeCount);
  
        
  /** \brief  Construct custom 2-cell (polygon) from a list of edges.
   *          The default base topology is the specified custom cell topology.
   *  \param  name             [in]  - descriptive name of the custom 2-cell
   *  \param  vertex_count     [in]  - number of vertices in the custom 2-cell
   *  \param  node_count       [in]  - number of nodes in the custom 2-cell
   *  \param  edges            [in]  - raw CellTopologyData for each edge (can be different!)
   *  \param  edge_node_map    [in]  - flat array with node maps for each edge 
   *  \param  base             [in]  - CellTopologyData of the base topology
   */
  CellTopology( const std::string                             & name,
                const unsigned                                  vertex_count,
                const unsigned                                  node_count,
                const std::vector< const CellTopologyData * > & edges ,
                const std::vector< unsigned >                 & edge_node_map ,
                const CellTopologyData                        * base = NULL );

  
  /** \brief  Construct custom 3-cell (polyhedron) from a list of edges and sides.
   *          The default base topology is the specified custom cell topology.
   *  \param  name             [in]  - descriptive name of the custom 3-cell
   *  \param  vertex_count     [in]  - number of vertices in the custom 3-cell
   *  \param  node_count       [in]  - number of nodes in the custom 3-cell
   *  \param  edges            [in]  - raw CellTopologyData for each edge (can be different!)
   *  \param  edge_node_map    [in]  - flat array with node maps for each edge 
   *  \param  faces            [in]  - raw CellTopologyData for each face (can be different!)
   *  \param  face_node_map    [in]  - flat array with node maps for each face 
   *  \param  base             [in]  - CellTopologyData of the base topology
   */
  CellTopology( const std::string                             & name,
                const unsigned                                  vertex_count,
                const unsigned                                  node_count,
                const std::vector< const CellTopologyData * > & edges ,
                const std::vector< unsigned >                 & edge_node_map ,
                const std::vector< const CellTopologyData * > & faces ,
                const std::vector< unsigned >                 & face_node_map ,
                const CellTopologyData                        * base = NULL );
  
  
  /** \brief Assignment operator <var>*this = right</var>.  */
  CellTopology& operator = (const CellTopology& right);

  /** \brief Copy constructor */
  CellTopology( const CellTopology& right );

  /** \brief  Default constructor initializes to NULL */
  CellTopology();

  /** \brief Destructor */
  ~CellTopology();
  
  /** \} */

}; // class CellTopology

/*------------------------------------------------------------------------*/
/* \brief  Find the permutation from the expected nodes to the actual nodes,
 *
 *  Find permutation 'p' such that:
 *    actual_node[j] == expected_node[ top.permutation[p].node[j] ]
 *  for all vertices.
 */
template< typename id_type >
int findPermutation( const CellTopologyData & top ,
                     const id_type * const expected_node ,
                     const id_type * const actual_node )
{
  const int nv = top.vertex_count ;
  const int np = top.permutation_count ;
  int p = 0 ;
  for ( ; p < np ; ++p ) {
    const unsigned * const perm_node = top.permutation[p].node ;
    int j = 0 ;
    for ( ; j < nv && actual_node[j] == expected_node[ perm_node[j] ] ; ++j );
    if ( nv == j ) break ;
  }
  if ( np == p ) p = -1 ;
  return p ;
}

template< typename id_type >
int findPermutation( const CellTopology & top ,
                     const id_type * const expected_node ,
                     const id_type * const actual_node )
{
  return findPermutation( * top.getCellTopologyData() , expected_node , actual_node );
}

/*------------------------------------------------------------------------*/
/** \brief  Generates detailed message if one or more input parameters are
 *          out of their admissible bounds.
 *  \param  dimension          [in]  - maximum value = 7
 *  \param  face_count         [in]  - maximum value = 63
 *  \param  edge_count         [in]  - maximum value = 63
 *  \param  vertex_count       [in]  - maximum value = 63
 *  \param  node_count         [in]  - maximum value = 1023
 *
 */
void badCellTopologyKey( const unsigned dimension ,
                         const unsigned face_count ,
                         const unsigned edge_count ,
                         const unsigned vertex_count ,
                         const unsigned node_count );
  

/** \brief  Generate integer key from topological dimensions
 *  \param  dimension          [in]  - maximum value = 7    (3 bits)
 *  \param  face_count         [in]  - maximum value = 63   (6 bits)
 *  \param  edge_count         [in]  - maximum value = 63   (6 bits)
 *  \param  vertex_count       [in]  - maximum value = 63   (6 bits)
 *  \param  node_count         [in]  - maximum value = 1023 (10 bits)
 *          The key uses all but the first bit in a 32 bit unsigned. 
 */
inline
unsigned cellTopologyKey( const unsigned dimension ,
                          const unsigned face_count ,
                          const unsigned edge_count ,
                          const unsigned vertex_count ,
                          const unsigned node_count )
{
  const bool bad = ( dimension    >> 3 ) ||
                   ( face_count   >> 6 ) ||
                   ( edge_count   >> 6 ) ||
                   ( vertex_count >> 6 ) ||
                   ( node_count   >> 10 );

  if ( bad ) {
    badCellTopologyKey( dimension , 
                        face_count ,
                        edge_count , 
                        vertex_count , 
                        node_count );
  }

  const unsigned key = ( dimension    << 28  ) |
                       ( face_count   << 22  ) |
                       ( edge_count   << 16  ) |
                       ( vertex_count << 10  ) |
                       ( node_count          ) ;

  return key ;
}

inline
bool operator==(const CellTopology &left, const CellTopology &right)
{
  return left.getCellTopologyData() == right.getCellTopologyData();
}

inline
bool operator<(const CellTopology &left, const CellTopology &right)
{
  return left.getCellTopologyData() < right.getCellTopologyData();
// FIXME: Should is be this?  
//  return left.getKey() < right.getKey(); 
}

inline
bool operator!=(const CellTopology &left, const CellTopology &right) {
  return !(left == right);
}


/** \} */

} // namespace shards

#undef SHARDS_REQUIRE

#endif // Shards_CellTopology_hpp

