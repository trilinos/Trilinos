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

namespace shards {

/** \addtogroup  shards_package_cell_topology
 *  \{
 */

/*------------------------------------------------------------------------*/

class CellTopology ;
class CellTopologyPrivate ;

/** \brief Overloaded << operator for CellTopologyData objects. */
std::ostream & operator << ( std::ostream & , const CellTopologyData & );

/** \brief Overloaded << operator for CellTopology objects.  */
std::ostream & operator << ( std::ostream & , const CellTopology & );

/*------------------------------------------------------------------------*/
/** \class shards::CellTopology
 *  \brief Provide safe access (in debug mode) to cell topology data
 *         and procedure to create custom cell topologies. 
 *  \author Created by P. Bochev, H. C. Edwards and D. Ridzal.
 *  \nosubgrouping
 */
class CellTopology {
private:
  void deleteOwned();
  void requireCell() const ;
  void requireDimension( const unsigned subcell_dim ) const ;
  void requireSubcell( const unsigned subcell_dim ,
                       const unsigned subcell_ord ) const ;
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
   */
  const CellTopologyData * getBaseTopology( const unsigned subcell_dim ,
                                            const unsigned subcell_ord ) const
    {
      return getTopology(subcell_dim,subcell_ord)->base ;
    }

        
   /** \brief  Key of a subcell of the given dimension and ordinal. */
   unsigned getKey( const unsigned subcell_dim ,
                    const unsigned subcell_ord ) const
     {
       return getTopology(subcell_dim,subcell_ord) -> key ;
     }
        
        

  /** \brief  Key of a subcell of the given dimension and ordinal. */
  unsigned getKey( const unsigned subcell_dim ,
                   const unsigned subcell_ord ) const
    {
      return getTopology(subcell_dim,subcell_ord)->key ;
    }

  /** \brief  Node count of a subcell of the given dimension and ordinal. */
  unsigned getNodeCount( const unsigned subcell_dim ,
                         const unsigned subcell_ord ) const
    {
      return getTopology(subcell_dim,subcell_ord)->node_count ;
    }

  /** \brief  Vertex count of a subcell of the given dimension and ordinal. */
  unsigned getVertexCount( const unsigned subcell_dim ,
                           const unsigned subcell_ord ) const
    {
      return getTopology(subcell_dim,subcell_ord)->vertex_count ;
    }

  /** \brief  Edge count of a subcell of the given dimension and ordinal. */
  unsigned getEdgeCount( const unsigned subcell_dim ,
                         const unsigned subcell_ord ) const
    {
      return getTopology(subcell_dim,subcell_ord)->edge_count ;
    }
  
  /** \brief  Side count of a subcell of the given dimension and ordinal. */
  unsigned getSideCount( const unsigned subcell_dim ,
                         const unsigned subcell_ord ) const
    {
      return getTopology(subcell_dim,subcell_ord)->side_count ;
    }

  /** \brief  Subcell count of subcells of the given dimension. */
  unsigned getSubcellCount( const unsigned subcell_dim ) const
    {
      SHARDS_REQUIRE( requireCell() );
      SHARDS_REQUIRE( requireDimension(subcell_dim) );
      return m_cell->subcell_count[subcell_dim] ;
    }
  
  /** \brief  Query if all subcells of the given dimension
   *          have the same cell topology.
   */
  bool getSubcellHomogeneity( const unsigned subcell_dim ) const
    {
      SHARDS_REQUIRE( requireCell() );
      SHARDS_REQUIRE( requireDimension(subcell_dim) );
      return m_cell->subcell_homogeneity[subcell_dim] ;
    }
  
  /** \brief  Node count of subcell of given dimension and ordinal.    */
   unsigned getNodeCount( const unsigned  subcell_dim ,
                          const unsigned  subcell_ord ) const
          {
            SHARDS_REQUIRE( requireCell() );
            SHARDS_REQUIRE( requireSubcell(subcell_dim,subcell_ord) );
            return m_cell -> subcell[subcell_dim][subcell_ord].topology -> node_count;
          }
        
  /** \brief  Mapping from a subcell's node ordinal to a
   *          node ordinal of this parent cell topology.
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
  
  
  /** \brief  Wrapper for safe access to a raw cell topology data. */
  CellTopology( const CellTopologyData * cell )
    : m_cell( cell ), m_owned( NULL )
    { SHARDS_REQUIRE( requireCell() ); }
  
  /** \brief Constructs custom 1-cell (line) with base topology Line<>. */
  CellTopology( const std::string & name,
                const unsigned      nodeCount);
  
  /** \brief  Construct custom 2-cell (polygon) from a list of edges.
   *          The default base topology is the specified custom cell topology.
   */
  CellTopology( const std::string                             & name,
                const unsigned                                  vertex_count,
                const unsigned                                  node_count,
                const std::vector< const CellTopologyData * > & edges ,
                const std::vector< unsigned >                 & edge_node_map ,
                const CellTopologyData                        * base = NULL );

  /** \brief  Construct custom 3-cell (polyhedron) from a list of
   *          edges and sides.
   *          The default base topology is the specified custom cell topology.
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

void badCellTopologyKey( const unsigned dimension ,
                         const unsigned face_count ,
                         const unsigned edge_count ,
                         const unsigned vertex_count ,
                         const unsigned node_count );
  
/** \brief  Generate integer key from topological dimensions
 *  \param  dimension    maximum value = 7
 *  \param  face_count   maximum value = 63
 *  \param  edge_count   maximum value = 63
 *  \param  vertex_count maximum value = 63
 *  \param  node_count   maximum value = 1023
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
    badCellTopologyKey( dimension , face_count ,
                                    edge_count , vertex_count , node_count );
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

