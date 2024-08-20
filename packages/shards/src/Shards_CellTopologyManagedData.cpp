// @HEADER
// *****************************************************************************
//                Shards : Shared Discretization Tools
//
// Copyright 2008-2011 NTESS and the Shards contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Shards_CellTopologyManagedData.hpp>

namespace shards {

CellTopologyManagedData *
createCellTopology(
  const std::string & name )
{
  return new CellTopologyManagedData( name );
}



CellTopologyManagedData *
createCellTopology(
  const std::string & name,
  const unsigned      node_count )
{
  return new CellTopologyManagedData( name , node_count );
}



CellTopologyManagedData *
createCellTopology(
  const std::string                             & name,
  const unsigned                                  vertex_count ,
  const unsigned                                  node_count,
  const std::vector< const CellTopologyData * > & edges ,
  const std::vector< unsigned >                 & edge_node_map ,
  const CellTopologyData                        * base )
{
  return new CellTopologyManagedData( name , 
                                      vertex_count , node_count ,
                                      edges , edge_node_map , 
                                      base );
}       

CellTopologyManagedData *
createCellTopology( const std::string                             & name,
                    const unsigned                                  vertex_count,
                    const unsigned                                  node_count,
                    const std::vector< const CellTopologyData * > & edges ,
                    const std::vector< unsigned >                 & edge_node_map ,
                    const std::vector< const CellTopologyData * > & faces ,
                    const std::vector< unsigned >                 & face_node_map ,
                    const CellTopologyData                        * base)
{
  return new CellTopologyManagedData( name , 
                                      vertex_count , node_count ,
                                      edges , edge_node_map , 
                                      faces , face_node_map , 
                                      base );  
}


//----------------------------------------------------------------------

CellTopologyManagedData::CellTopologyManagedData(
  const std::string & name_)
  : m_name(name_),
    m_subcell(1),
    m_node_map()
{
 
  // This 1-subcell is the cell itself, will be assigned to subcell[1]
  m_subcell[0].topology = this;
  m_subcell[0].node     = index_identity_array();

  base = 0;
  name = m_name.c_str();
  key  = 0 ;
  dimension              = 0;
  vertex_count           = 0 ;
  node_count             = 0 ;                                  // PBB 12-03-08
  edge_count             = 0 ;
  side_count             = 0 ;
  subcell_homogeneity[0] = 0 ;
  subcell_homogeneity[1] = 0 ;
  subcell_homogeneity[2] = 0 ;
  subcell_homogeneity[3] = 0 ;
  subcell_count[0]       = 0 ;
  subcell_count[1]       = 0 ;
  subcell_count[2]       = 0 ;
  subcell_count[3]       = 0 ;
  subcell[0]             = 0 ;
  subcell[1]             = 0 ;
  subcell[2]             = NULL ;
  subcell[3]             = NULL ;
  side                   = NULL ;
  edge                   = NULL ;
}


CellTopologyManagedData::CellTopologyManagedData(
  const std::string & name_,
  const unsigned      node_count_ )
  : m_name(name_),
    m_subcell(1),
    m_node_map()
{
 
  // This 1-subcell is the cell itself, will be assigned to subcell[1]
  m_subcell[0].topology = this;
  m_subcell[0].node     = index_identity_array();

  base = getCellTopologyData< Line<2> >();
  name = m_name.c_str();
  key  = cellTopologyKey( 1 , 0 , 0 , 2 , node_count_ );
  dimension              = 1 ;
  vertex_count           = 2 ;
  node_count             = node_count_;                                  // PBB 12-03-08
  edge_count             = 0 ;
  side_count             = 0 ;
  subcell_homogeneity[0] = 1 ;
  subcell_homogeneity[1] = 0 ;
  subcell_homogeneity[2] = 0 ;
  subcell_homogeneity[3] = 0 ;
  subcell_count[0]       = node_count_ ;
  subcell_count[1]       = 1 ;
  subcell_count[2]       = 0 ;
  subcell_count[3]       = 0 ;
  subcell[0]             = subcell_nodes_array();
  subcell[1]             = & m_subcell[0] ;
  subcell[2]             = NULL ;
  subcell[3]             = NULL ;
  side                   = NULL ;
  edge                   = NULL ;
}


//2D --------------------------------------------------------------------


CellTopologyManagedData::CellTopologyManagedData(
  const std::string                             & name_,
  const unsigned                                  vertex_count_ ,
  const unsigned                                  node_count_,
  const std::vector< const CellTopologyData * > & edges ,
  const std::vector< unsigned >                 & edge_node_map ,
  const CellTopologyData                        * base )
  : m_name(name_),
    m_subcell(),
    m_node_map()
{

  // Compute size of the edge map & check edge homogeneity (suffices to compare nodes per edge)
  const unsigned edge_count_ = edges.size();
  unsigned       edge_map_size = 0 ;
  unsigned       node_count_edge0 = edges[0]->node_count;
  bool           edge_homogeneity = true;
      
  for ( unsigned i = 0 ; i < edge_count_ ; ++i ) {
    edge_map_size += edges[i]->node_count ;
    if(node_count_edge0 != edges[i]->node_count ) edge_homogeneity = false;
  }

  const bool error_base = base && (
    base->base         != base ||
    base->vertex_count != base->node_count ||
    base->vertex_count != vertex_count_ ||
    base->edge_count   != edges.size() );

  const bool error_base_self = ! base && ( vertex_count_ != node_count_ );

  const bool error_edge = edge_map_size != edge_node_map.size();

  if ( error_base || error_base_self || error_edge ) {
    // Throw an error
  }

  m_subcell.resize( 1 + edges.size() );
  m_subcell[ edge_count_ ].topology = this ;            // This subcell is the cell itself, will be assigned to subcell[2]
  m_subcell[ edge_count_ ].node     = index_identity_array();

  m_node_map.resize( edge_map_size );

  for ( unsigned i = 0 ; i < edge_map_size ; ++i ) {
    m_node_map[i] = edge_node_map[i];
  }

  edge_map_size = 0 ;
  for ( unsigned i = 0 ; i < edge_count_ ; ++i ) {
    m_subcell[i].topology = edges[i] ;
    m_subcell[i].node     = & m_node_map[ edge_map_size ] ;
    edge_map_size += edges[i]->node_count ;
  }

  base = (base == NULL ? this : base) ;                            // PBB 12-03-08
  name = m_name.c_str();
  key  = cellTopologyKey( 2, 0, edge_count_, vertex_count_, node_count_ );
  dimension              = 2 ;
  vertex_count           = vertex_count_ ;
  node_count             = node_count_ ;                                // PBB 12-03-08
  edge_count             = edge_count_ ;
  side_count             = 0 ;
  subcell_homogeneity[0] = 1 ;
  subcell_homogeneity[1] = edge_homogeneity ;
  subcell_homogeneity[2] = 0 ;
  subcell_homogeneity[3] = 0 ;
  subcell_count[0]       = node_count_ ;
  subcell_count[1]       = edge_count_ ;
  subcell_count[2]       = 1 ;                                           // PBB 12-03-08
  subcell_count[3]       = 0 ;
  subcell[0]             = subcell_nodes_array();
  subcell[1]             = & m_subcell[0] ;
  subcell[2]             = & m_subcell[ edge_count_ ] ;
  subcell[3]             = NULL ;
  side                   = subcell[1] ;
  edge                   = subcell[1] ;
}


//3D--------------------------------------------------------------------

CellTopologyManagedData::CellTopologyManagedData(
  const std::string                             & name_,
  const unsigned                                  vertex_count_,
  const unsigned                                  node_count_,
  const std::vector< const CellTopologyData * > & edges ,
  const std::vector< unsigned >                 & edge_node_map ,
  const std::vector< const CellTopologyData * > & faces ,
  const std::vector< unsigned >                 & face_node_map ,
  const CellTopologyData                        * base ) 
  : m_name(name_),
    m_subcell(),
    m_node_map()

{
  const unsigned edge_count_ = edges.size();
  unsigned edge_map_size = 0 ;
  unsigned node_count_edge0 = edges[0]->node_count;
  bool edge_homogeneity = true;

  // Compute size of the edge map & check edge homogeneity (suffices to compare nodes per edge)
  for ( unsigned i = 0 ; i < edge_count_ ; ++i ) {
    edge_map_size += edges[i]->node_count ;
    if(node_count_edge0 != edges[i]->node_count ) edge_homogeneity = false;
  }
  
  // Compute size of the face map & check face homogeneity (to do)
  const unsigned face_count_ = faces.size();
  unsigned face_map_size = 0 ;
  for ( unsigned i = 0 ; i < face_count_ ; ++i ) {
    face_map_size += faces[i]->node_count ;
  }
  
  // Set error flags for base, edges and faces
  const bool error_base = base && (base->base         != base             ||
                                   base->vertex_count != base->node_count ||
                                   base->vertex_count != vertex_count_    ||
                                   base->edge_count   != edges.size()     ||
                                   base->side_count   != faces.size() );
  
  const bool error_base_self = ! base && ( vertex_count_ != node_count_ );
  
  const bool error_edge = edge_map_size != edge_node_map.size();
  
  const bool error_face = face_map_size != face_node_map.size();
  
  if ( error_base || error_base_self || error_edge || error_face) {
    // Throw an error
  }
  
  // Flat array for the 1,2,3-subcells of the custom cell.
  m_subcell.resize( 1 + edges.size() + faces.size() );
  
  m_subcell[ edge_count_ + face_count_].topology = this ;               // The last subcell is the cell itself & will be assigned to subcell[3]
  m_subcell[ edge_count_ + face_count_].node     = index_identity_array();
  
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
  for ( unsigned i = 0 ; i < edge_count_ ; ++i ) {
    m_subcell[i].topology = edges[i] ;
    m_subcell[i].node     = & m_node_map[ edge_map_size ] ;
    edge_map_size        += edges[i]->node_count ;
  }
  
  // Copy face topologies & list of nodes for each face to m_subcell:
  face_map_size = 0;
  for ( unsigned i = 0 ; i < face_count_ ; ++i ) {
    m_subcell[edge_count_ + i].topology = faces[i] ;
    m_subcell[edge_count_ + i].node     = & m_node_map[ edge_map_size + face_map_size ] ;
    face_map_size                     += faces[i]->node_count ;
  }
  
  // Fill CellTopologyData with custom cell data: default base is the custom cell itself
  base = (base == NULL ? this : base) ;                           
  name = m_name.c_str();
  key  = cellTopologyKey( 3, 
                          face_count_, 
                          edge_count_, 
                          vertex_count_, 
                          node_count_ );
  dimension              = 3 ;
  vertex_count           = vertex_count_ ;
  node_count             = node_count_ ;                                
  edge_count             = edge_count_ ;
  side_count             = face_count_ ;
  subcell_homogeneity[0] = 1 ;
  subcell_homogeneity[1] = edge_homogeneity ;
  subcell_homogeneity[2] = 0 ;
  subcell_homogeneity[3] = 0 ;
  subcell_count[0]       = node_count_ ;
  subcell_count[1]       = edge_count_ ;
  subcell_count[2]       = face_count_ ;                                           
  subcell_count[3]       = 1 ;
  subcell[0]             = subcell_nodes_array();
  subcell[1]             = & m_subcell[0] ;
  subcell[2]             = & m_subcell[ edge_count_ ] ;
  subcell[3]             = & m_subcell[ edge_count_ + face_count_] ;
  side                   = subcell[2] ;
  edge                   = subcell[1] ;  
}

} // namespace shards
