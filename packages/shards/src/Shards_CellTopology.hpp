/*------------------------------------------------------------------------*/
/*               shards : Shared Discretization Tools                     */
/*                Copyright (2008) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/* Questions? Contact Pavel Bochev      (pbboche@sandia.gov)              */
/*                    H. Carter Edwards (hcedwar@sandia.gov)              */
/*                    Denis Ridzal      (dridzal@sandia.gov).             */
/*------------------------------------------------------------------------*/

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
class CellTopologyPrivate ;

/** \brief Overloaded << operator for CellTopologyData objects. */
//std::ostream & operator << ( std::ostream & , const CellTopologyData & );

/** \brief Overloaded << operator for CellTopology objects.  */
std::ostream & operator << ( std::ostream & , const CellTopology & );

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
 */
class CellTopology {
private:
  
  /** \brief Deletes m_owned data member
  */
  void deleteOwned();
  
  
  /** \brief Throws runtime_error if CellTopology object is null or hase null 
             base topology   
   */
  void requireCell() const ;
  
  
  /** \brief  Throws invalid_argument if subcell dimension exceedes the maximal 
   *          admissible space dimension 3.
   *  \param  subcell_dim    [in]  - spatial dimension of a subcell
   */
  void requireDimension( const unsigned subcell_dim ) const ;
  
  
  /** \brief  Throws invalid_argument if subcell_ord exceeds the actual number
   *          of subcells with specified dimension.
   *  \param  subcell_dim    [in]  - spatial dimension of a subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  void requireSubcell( const unsigned subcell_dim ,
                       const unsigned subcell_ord ) const ;
  
  
  /** \brief  Throws invalid_argument if node_ord exceeds the actual number 
   *          of nodes in the subcell with specified dimension and ordinal.
   *  \param  subcell_dim    [in]  - spatial dimension of a subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   *  \param  node_ord       [in]  - node ordinal
   */
  void requireNodeMap( const unsigned subcell_dim ,
                       const unsigned subcell_ord ,
                       const unsigned node_ord ) const ;

  
  const CellTopologyData    * m_cell ;
        CellTopologyPrivate * m_owned ;

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
   */
  unsigned getKey() const
    {
      SHARDS_REQUIRE( requireCell() );
      return m_cell->key ;
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
  
        
  /** \brief  Side boundary subcell count of this cell topology */
  unsigned getSideCount() const
    {
      SHARDS_REQUIRE( requireCell() );
      return m_cell->side_count ;
    }
  
        
  /** \brief  This cell's raw topology data */
  const CellTopologyData * getTopology() const
    {
      SHARDS_REQUIRE( requireCell() );
      return m_cell ;
    }

        
  /** \brief  This cell's base cell topology's raw topology data */
  const CellTopologyData * getBaseTopology() const
    {
      SHARDS_REQUIRE( requireCell() );
      return m_cell->base ;
    }

        
  /** \brief  Raw cell topology data for a subcell of the
   *          given dimension and ordinal.
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  const CellTopologyData * getTopology( const unsigned subcell_dim ,
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
  const CellTopologyData * getBaseTopology( const unsigned subcell_dim ,
                                            const unsigned subcell_ord ) const
    {
      return getTopology(subcell_dim,subcell_ord)->base ;
    }

        
  /** \brief  Key of a subcell of the given dimension and ordinal. 
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  unsigned getKey( const unsigned subcell_dim ,
                   const unsigned subcell_ord ) const
    {
      return getTopology(subcell_dim,subcell_ord)->key ;
    }

        
  /** \brief  Node count of a subcell of the given dimension and ordinal. 
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  unsigned getNodeCount( const unsigned subcell_dim ,
                         const unsigned subcell_ord ) const
    {
      return getTopology(subcell_dim,subcell_ord)->node_count ;
    }

        
  /** \brief  Vertex count of a subcell of the given dimension and ordinal. 
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  unsigned getVertexCount( const unsigned subcell_dim ,
                           const unsigned subcell_ord ) const
    {
      return getTopology(subcell_dim,subcell_ord)->vertex_count ;
    }

        
  /** \brief  Edge count of a subcell of the given dimension and ordinal. 
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */
  unsigned getEdgeCount( const unsigned subcell_dim ,
                         const unsigned subcell_ord ) const
    {
      return getTopology(subcell_dim,subcell_ord)->edge_count ;
    }
  
        
  /** \brief  Side count of a subcell of the given dimension and ordinal. 
   *  \param  subcell_dim    [in]  - spatial dimension of the subcell
   *  \param  subcell_ord    [in]  - subcell ordinal
   */        
  unsigned getSideCount( const unsigned subcell_dim ,
                         const unsigned subcell_ord ) const
    {
      return getTopology(subcell_dim,subcell_ord)->side_count ;
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
      return m_cell->subcell_homogeneity[subcell_dim] ;
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
  
  /** \} */
  
  /*------------------------------------------------------------------*/
  /** \name Constructors and destructor
   *  \{
   */

  /** \brief  Default constructor for invalid topology. */
  CellTopology() : m_cell(NULL), m_owned(NULL) {}
  
  
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
    : m_cell( cell ), m_owned( NULL )
    { SHARDS_REQUIRE( requireCell() ); }
  
        
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
  
  
  /** \brief Destructor */
  ~CellTopology() { if ( m_owned ) { deleteOwned(); } }
  
  /** \} */

}; // class CellTopology


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



/** \} */

} // namespace shards

#undef SHARDS_REQUIRE

#endif // Shards_CellTopology_hpp

