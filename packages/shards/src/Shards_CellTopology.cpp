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

//#define HAVE_SHARDS_DEBUG

#include <stdexcept>
#include <sstream>
#include <Shards_CellTopology.hpp>
//#include <Shards_BasicTopologies.hpp>
#include <iostream>

namespace shards {

typedef CellTopologyData::Subcell Subcell ;

class CellTopologyPrivate {
public:
  CellTopologyData        m_cell ;
  std::string             m_name ;
  std::vector< Subcell >  m_subcell ;
  std::vector< unsigned > m_node_map ;

  /** \brief  1D topology */
  CellTopologyPrivate( const std::string & name , const unsigned nodeCount );

  /** \brief  2D topology */
  CellTopologyPrivate(
    const std::string                             & name,
    const unsigned                                  vertexCount,
    const unsigned                                  nodeCount,
    const std::vector< const CellTopologyData * > & edges ,
    const std::vector< unsigned >                 & edge_node_map ,
    const CellTopologyData                        * base );

  /** \brief  3D topology */
  CellTopologyPrivate(
    const std::string                             & name,
    const unsigned                                  vertexCount,
    const unsigned                                  nodeCount,
    const std::vector< const CellTopologyData * > & edges ,
    const std::vector< unsigned >                 & edge_node_map ,
    const std::vector< const CellTopologyData * > & faces ,
    const std::vector< unsigned >                 & face_node_map ,
    const CellTopologyData                        * base );

private:
  CellTopologyPrivate();
  CellTopologyPrivate( const CellTopologyPrivate & );
  CellTopologyPrivate & operator = ( const CellTopologyPrivate & );
};

//----------------------------------------------------------------------

void CellTopology::deleteOwned()
{ delete m_owned ; }

CellTopology::CellTopology(
  const std::string & name,
  const unsigned      node_count )
  : m_cell(NULL), m_owned(NULL)
{
  m_owned = new CellTopologyPrivate( name , node_count );
  m_cell  = & m_owned->m_cell ;
}



CellTopology::CellTopology(
  const std::string                             & name,
  const unsigned                                  vertex_count ,
  const unsigned                                  node_count,
  const std::vector< const CellTopologyData * > & edges ,
  const std::vector< unsigned >                 & edge_node_map ,
  const CellTopologyData                        * base )
  : m_cell(NULL), m_owned(NULL)
{
  m_owned = new CellTopologyPrivate( name , 
                                     vertex_count , node_count ,
                                     edges , edge_node_map , 
                                     base );
  m_cell  = & m_owned->m_cell ;
}



CellTopology::CellTopology( const std::string                             & name,
                            const unsigned                                  vertex_count,
                            const unsigned                                  node_count,
                            const std::vector< const CellTopologyData * > & edges ,
                            const std::vector< unsigned >                 & edge_node_map ,
                            const std::vector< const CellTopologyData * > & faces ,
                            const std::vector< unsigned >                 & face_node_map ,
                            const CellTopologyData                        * base)
: m_cell(NULL), m_owned(NULL)
{
  m_owned = new CellTopologyPrivate( name , 
                                     vertex_count , node_count ,
                                     edges , edge_node_map , 
                                     faces , face_node_map , 
                                     base );
  m_cell  = & m_owned -> m_cell ;
  
}


//----------------------------------------------------------------------

void CellTopology::requireCell() const
{
  if ( m_cell == NULL || m_cell->base == NULL ) {
    std::string msg ;
    msg.append( "shards::CellTopology::requireCell() : FAILED " );
    if ( m_cell == NULL ) {
      msg.append( "is NULL" );
    }
    else {
      msg.append("'");
      msg.append( m_cell->name );
      msg.append("' has NULL base");
    }
    msg.append( " ) FAILED" );
    throw std::runtime_error( msg );
  }
}

void CellTopology::requireDimension( const unsigned subcellDim ) const
{
  if ( 3 < subcellDim ) {
    std::ostringstream msg ;
    msg << "shards::CellTopology::requireDimension( ERROR: dim = "
        << subcellDim << " > 3 )" ;
    throw std::invalid_argument( msg.str() );
  }
}

void CellTopology::requireSubcell( const unsigned subcellDim ,
                                   const unsigned subcellOrd ) const
{
  if ( m_cell->subcell_count[ subcellDim ] <= subcellOrd ) {
    std::ostringstream msg ;
    msg << "shards::CellTopology::requireSubcell( dim = "
        << subcellDim << " , ERROR: ord = " << subcellOrd
        << " > '" << m_cell->name
        << "'.subcell_count[" << subcellDim
        << "] = " << m_cell->subcell_count[ subcellDim ]
        << " )" ;
    throw std::invalid_argument( msg.str() );
  }
}

void CellTopology::requireNodeMap( const unsigned subcellDim ,
                                   const unsigned subcellOrd ,
                                   const unsigned nodeOrd ) const
{
  const unsigned n =
    m_cell->subcell[subcellDim][subcellOrd].topology->node_count ;

  if ( n <= nodeOrd ) {
    std::ostringstream msg ;
    msg << "shards::CellTopology::requireNodeMap( " 
        << " , " << subcellDim
        << " , " << subcellOrd
        << " , ERROR: " << nodeOrd
        << " >= '"
        << m_cell->name 
        << "'.subcell[" << subcellDim
        << "][" << subcellOrd
        << "].topology->node_count = "
        << n << " )" ;
    throw std::invalid_argument( msg.str() );
  }
}

// 1D --------------------------------------------------------------------

CellTopologyPrivate::CellTopologyPrivate(
  const std::string & name,
  const unsigned      node_count )
  : m_cell(),
    m_name(name),
    m_subcell(1),
    m_node_map()
{
 
  // This 1-subcell is the cell itself, will be assigned to m_cell.subcell[1]
  m_subcell[0].topology = & m_cell ;
  m_subcell[0].node     = index_identity_array();

  m_cell.base = getCellTopologyData< Line<2> >();
  m_cell.name = m_name.c_str();
  m_cell.key  = cellTopologyKey( 1 , 0 , 0 , 2 , node_count );
  m_cell.dimension              = 1 ;
  m_cell.vertex_count           = 2 ;
  m_cell.node_count             = node_count;                                  // PBB 12-03-08
  m_cell.edge_count             = 0 ;
  m_cell.side_count             = 0 ;
  m_cell.subcell_homogeneity[0] = 1 ;
  m_cell.subcell_homogeneity[1] = 0 ;
  m_cell.subcell_homogeneity[2] = 0 ;
  m_cell.subcell_homogeneity[3] = 0 ;
  m_cell.subcell_count[0]       = node_count ;
  m_cell.subcell_count[1]       = 1 ;
  m_cell.subcell_count[2]       = 0 ;
  m_cell.subcell_count[3]       = 0 ;
  m_cell.subcell[0]             = subcell_nodes_array();
  m_cell.subcell[1]             = & m_subcell[0] ;
  m_cell.subcell[2]             = NULL ;
  m_cell.subcell[3]             = NULL ;
  m_cell.side                   = NULL ;
  m_cell.edge                   = NULL ;
}


//2D --------------------------------------------------------------------


CellTopologyPrivate::CellTopologyPrivate(
  const std::string                             & name,
  const unsigned                                  vertex_count ,
  const unsigned                                  node_count,
  const std::vector< const CellTopologyData * > & edges ,
  const std::vector< unsigned >                 & edge_node_map ,
  const CellTopologyData                        * base )
  : m_cell(),
    m_name(name),
    m_subcell(),
    m_node_map()
{

  // Compute size of the edge map & check edge homogeneity (suffices to compare nodes per edge)
  const unsigned edge_count = edges.size();
  unsigned       edge_map_size = 0 ;
  unsigned       node_count_edge0 = edges[0] -> node_count;
  bool           edge_homogeneity = true;
      
  for ( unsigned i = 0 ; i < edge_count ; ++i ) {
    edge_map_size += edges[i]->node_count ;
    if(node_count_edge0 != edges[i] -> node_count ) edge_homogeneity = false;
  }

  const bool error_base = base && (
                          base->base         != base ||
                          base->vertex_count != base->node_count ||
                          base->vertex_count != vertex_count ||
                          base->edge_count   != edges.size() );

  const bool error_base_self = ! base && ( vertex_count != node_count );

  const bool error_edge = edge_map_size != edge_node_map.size();

  if ( error_base || error_base_self || error_edge ) {
    // Throw an error
  }

  m_subcell.resize( 1 + edges.size() ),
    
  // This subcell is the cell itself, will be assigned to m_cell.subcell[2]
  m_subcell[ edge_count ].topology = & m_cell ;
  m_subcell[ edge_count ].node     = index_identity_array();

  m_node_map.resize( edge_map_size );

  for ( unsigned i = 0 ; i < edge_map_size ; ++i ) {
    m_node_map[i] = edge_node_map[i];
  }

  edge_map_size = 0 ;
  for ( unsigned i = 0 ; i < edge_count ; ++i ) {
    m_subcell[i].topology = edges[i] ;
    m_subcell[i].node     = & m_node_map[ edge_map_size ] ;
    edge_map_size += edges[i]->node_count ;
  }

  m_cell.base = (base == NULL ? & m_cell : base) ;                            // PBB 12-03-08
  m_cell.name = m_name.c_str();
  m_cell.key  = cellTopologyKey( 2, 0, edge_count, vertex_count, node_count );
  m_cell.dimension              = 2 ;
  m_cell.vertex_count           = vertex_count ;
  m_cell.node_count             = node_count ;                                // PBB 12-03-08
  m_cell.edge_count             = edge_count ;
  m_cell.side_count             = 0 ;
  m_cell.subcell_homogeneity[0] = 1 ;
  m_cell.subcell_homogeneity[1] = edge_homogeneity ;
  m_cell.subcell_homogeneity[2] = 0 ;
  m_cell.subcell_homogeneity[3] = 0 ;
  m_cell.subcell_count[0]       = node_count ;
  m_cell.subcell_count[1]       = edge_count ;
  m_cell.subcell_count[2]       = 1 ;                                           // PBB 12-03-08
  m_cell.subcell_count[3]       = 0 ;
  m_cell.subcell[0]             = subcell_nodes_array();
  m_cell.subcell[1]             = & m_subcell[0] ;
  m_cell.subcell[2]             = & m_subcell[ edge_count ] ;
  m_cell.subcell[3]             = NULL ;
  m_cell.side                   = m_cell.subcell[1] ;
  m_cell.edge                   = m_cell.subcell[1] ;
}


//3D--------------------------------------------------------------------

CellTopologyPrivate::CellTopologyPrivate(const std::string                             & name,
                                         const unsigned                                  vertex_count,
                                         const unsigned                                  node_count,
                                         const std::vector< const CellTopologyData * > & edges ,
                                         const std::vector< unsigned >                 & edge_node_map ,
                                         const std::vector< const CellTopologyData * > & faces ,
                                         const std::vector< unsigned >                 & face_node_map ,
                                         const CellTopologyData                        * base ) 
: m_cell(),
m_name(name),
m_subcell(),
m_node_map()

{
  const unsigned edge_count = edges.size();
  unsigned edge_map_size = 0 ;
  unsigned node_count_edge0 = edges[0] -> node_count;
  bool edge_homogeneity = true;

  // Compute size of the edge map & check edge homogeneity (suffices to compare nodes per edge)
  for ( unsigned i = 0 ; i < edge_count ; ++i ) {
    edge_map_size += edges[i] -> node_count ;
    if(node_count_edge0 != edges[i] -> node_count ) edge_homogeneity = false;
  }
  
  // Compute size of the face map & check face homogeneity (to do)
  const unsigned face_count = faces.size();
  unsigned face_map_size = 0 ;
  for ( unsigned i = 0 ; i < face_count ; ++i ) {
    face_map_size += faces[i] -> node_count ;
  }
  
  // Set error flags for base, edges and faces
  const bool error_base = base && (base->base         != base             ||
                                   base->vertex_count != base->node_count ||
                                   base->vertex_count != vertex_count     ||
                                   base->edge_count   != edges.size()     ||
                                   base->side_count   != faces.size() );
  
  const bool error_base_self = ! base && ( vertex_count != node_count );
  
  const bool error_edge = edge_map_size != edge_node_map.size();
  
  const bool error_face = face_map_size != face_node_map.size();
  
  if ( error_base || error_base_self || error_edge || error_face) {
    // Throw an error
  }
  
  // Flat array for the 1,2,3-subcells of the custom cell.
  m_subcell.resize( 1 + edges.size() + faces.size() ),
    
  // The last subcell is the cell itself & will be assigned to m_cell.subcell[3]
  m_subcell[ edge_count + face_count].topology = & m_cell ;
  m_subcell[ edge_count + face_count].node     = index_identity_array();
  
  // Flat array with edge nodes followed by face nodes (edges and faces can be non-homogeneous)
  m_node_map.resize( edge_map_size + face_map_size);
  
  // Copy edge nodes followed by face nodes
  for ( unsigned i = 0 ; i < edge_map_size ; ++i ) {
    m_node_map[i] = edge_node_map[i];
  }
  for ( unsigned i = 0 ; i < face_map_size ; ++i ) {
    m_node_map[edge_map_size + i] = face_node_map[i];
  }
  
  
  // Copy edge topologies & list of nodes for each edge to m_subcell:
  edge_map_size = 0 ;
  for ( unsigned i = 0 ; i < edge_count ; ++i ) {
    m_subcell[i].topology = edges[i] ;
    m_subcell[i].node     = & m_node_map[ edge_map_size ] ;
    edge_map_size        += edges[i] -> node_count ;
  }
  
  // Copy face topologies & list of nodes for each face to m_subcell:
  face_map_size = 0;
  for ( unsigned i = 0 ; i < face_count ; ++i ) {
    m_subcell[edge_count + i].topology = faces[i] ;
    m_subcell[edge_count + i].node     = & m_node_map[ edge_map_size + face_map_size ] ;
    face_map_size                     += faces[i] -> node_count ;
  }
  
  // Fill CellTopologyData with custom cell data: default base is the custom cell itself
  m_cell.base = (base == NULL ? & m_cell : base) ;                           
  m_cell.name = m_name.c_str();
  m_cell.key  = cellTopologyKey( 3, 
                                 face_count, 
                                 edge_count, 
                                 vertex_count, 
                                 node_count );
  m_cell.dimension              = 3 ;
  m_cell.vertex_count           = vertex_count ;
  m_cell.node_count             = node_count ;                                
  m_cell.edge_count             = edge_count ;
  m_cell.side_count             = face_count ;
  m_cell.subcell_homogeneity[0] = 1 ;
  m_cell.subcell_homogeneity[1] = edge_homogeneity ;
  m_cell.subcell_homogeneity[2] = 0 ;
  m_cell.subcell_homogeneity[3] = 0 ;
  m_cell.subcell_count[0]       = node_count ;
  m_cell.subcell_count[1]       = edge_count ;
  m_cell.subcell_count[2]       = face_count ;                                           
  m_cell.subcell_count[3]       = 1 ;
  m_cell.subcell[0]             = subcell_nodes_array();
  m_cell.subcell[1]             = & m_subcell[0] ;
  m_cell.subcell[2]             = & m_subcell[ edge_count ] ;
  m_cell.subcell[3]             = & m_subcell[ edge_count + face_count] ;
  m_cell.side                   = m_cell.subcell[2] ;
  m_cell.edge                   = m_cell.subcell[1] ;  
}



//----------------------------------------------------------------------


void badCellTopologyKey( const unsigned dimension ,
                         const unsigned face_count ,
                         const unsigned edge_count ,
                         const unsigned vertex_count ,
                         const unsigned node_count )
{
  const unsigned end_dimension    = 1u << 3 ;
  const unsigned end_face_count   = 1u << 6 ;
  const unsigned end_edge_count   = 1u << 6 ;
  const unsigned end_vertex_count = 1u << 6 ;
  const unsigned end_node_count   = 1u << 10 ;

  std::ostringstream msg ;
  msg << "shards::badCellTopologyKey( " ;
  msg << " dimension = " << dimension ;
  if ( dimension >= end_dimension )
    { msg << " > " << end_dimension << " ERROR"; }
  msg << " , face_count = " << face_count ;
  if ( face_count >= end_face_count )
    { msg << " > " << end_face_count << " ERROR"; }
  msg << " , edge_count = " << edge_count ;
  if ( edge_count >= end_edge_count )
    { msg << " > " << end_edge_count << " ERROR"; }
  msg << " , vertex_count = " << vertex_count ;
  if ( vertex_count >= end_vertex_count )
    { msg << " > " << end_vertex_count << " ERROR"; }
  msg << " , node_count = " << node_count ;
  if ( node_count >= end_node_count )
    { msg << " > " << end_node_count << " ERROR"; }
  msg << " )" ;

  throw std::invalid_argument( msg.str() );
}

std::ostream & operator << ( std::ostream & os, const CellTopology & cell) {
  os << *cell.getTopology();
  return os;
}

} // namespace shards


