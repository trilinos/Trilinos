// @HEADER
// *****************************************************************************
//                Shards : Shared Discretization Tools
//
// Copyright 2008-2011 NTESS and the Shards contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//#define HAVE_SHARDS_DEBUG

#include <stdexcept>
#include <sstream>
#include <Shards_CellTopology.hpp>
#include <Shards_BasicTopologies.hpp>
#include <iostream>

namespace shards {

typedef CellTopologyData_Subcell Subcell ;

//----------------------------------------------------------------------

CellTopology::~CellTopology()
{}

CellTopology::CellTopology()
  : m_cell( NULL )
{}

CellTopology::CellTopology( const CellTopology & right )
  : m_cell( NULL )
{
  operator=( right );
}

CellTopology & CellTopology::operator = ( const CellTopology & right )
{
  m_cell = right.m_cell ;

  return *this ;
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
        << subcellDim << " , "
        << subcellOrd
        << " , ERROR: " << nodeOrd << " >= '"
        << m_cell->name 
        << "'.subcell[" << subcellDim
        << "][" << subcellOrd
        << "].topology->node_count = "
        << n << " )" ;
    throw std::invalid_argument( msg.str() );
  }
}

void CellTopology::requireNodePermutation( const unsigned permutationOrd ,
                                           const unsigned nodeOrd ) const
{
  const bool bad_p = m_cell->permutation_count <= permutationOrd ;
  const bool bad_n = m_cell->node_count        <= nodeOrd ;
  if ( bad_p || bad_n ) {
    std::ostringstream msg ;
    msg << "shards::CellTopology::requireNodePermutation( " ;
    if ( bad_p ) {
      msg << " ERROR: " << permutationOrd << " >= "
          << m_cell->permutation_count ;
    }
    else {
      msg << permutationOrd ;
    }
    msg << " , " ;
    if ( bad_n ) {
      msg << " ERROR: " << nodeOrd << " >= " << m_cell->node_count ;
    }
    else {
      msg << nodeOrd ;
    }
    msg << " )" ;
    throw std::invalid_argument( msg.str() );
  }
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
  os << *cell.getCellTopologyData();
  return os;
}


void getTopologies(std::vector<shards::CellTopology>& topologies,
                   const unsigned cellDim,
                   const ECellType cellType,
                   const ETopologyType topologyType)
{
 
  if ( 4 < cellDim ) {
    std::ostringstream msg ;
    msg << "shards::CellTopology::getTopologies( ERROR: dim = "
      << cellDim << " > 4 )" ;
    throw std::invalid_argument( msg.str() );
  }

  // clear the vector  
  topologies.clear();
  
  // 0-dimensional cells
  if( (cellDim == 0) || (cellDim == 4) ) {
    if( cellType == STANDARD_CELL || cellType == ALL_CELLS){
      if(topologyType == BASE_TOPOLOGY || topologyType == ALL_TOPOLOGIES) {
        topologies.push_back( CellTopology( shards::getCellTopologyData<Node>() ) ); 
      }
      if(topologyType == EXTENDED_TOPOLOGY || topologyType == ALL_TOPOLOGIES) {
        // No such cells exist
      }
    }
    if( cellType == NONSTANDARD_CELL || cellType == ALL_CELLS){
      if(topologyType == BASE_TOPOLOGY || topologyType == ALL_TOPOLOGIES) {
        // No such cells exist
      }
      if(topologyType == EXTENDED_TOPOLOGY || topologyType == ALL_TOPOLOGIES) {
        // No such cells exist
      }
    }
  } // dim 0
  
  
  // 1-dimensional cells
  if((cellDim == 1) || (cellDim == 4)) {
    
    if( cellType == STANDARD_CELL || cellType == ALL_CELLS){
      if(topologyType == BASE_TOPOLOGY || topologyType == ALL_TOPOLOGIES) {
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Line<2> >() ) ); 
      }
      if(topologyType == EXTENDED_TOPOLOGY || topologyType == ALL_TOPOLOGIES) {
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Line<3> >() ) );
      }
    }    
    if( cellType == NONSTANDARD_CELL || cellType == ALL_CELLS){
      if(topologyType == BASE_TOPOLOGY || topologyType == ALL_TOPOLOGIES) {
        topologies.push_back( CellTopology( shards::getCellTopologyData<Particle>() ) );
      }
      if(topologyType == EXTENDED_TOPOLOGY || topologyType == ALL_TOPOLOGIES) {
        // No such cells exist
      }
    }
  } // dim 1
  
  
  // 2-dimensional cells
  if((cellDim == 2) || (cellDim == 4)) {
    if( cellType == STANDARD_CELL || cellType == ALL_CELLS){
      if(topologyType == BASE_TOPOLOGY || topologyType == ALL_TOPOLOGIES) {
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Triangle<3> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Quadrilateral<4> >() ) );
      }
      if (topologyType == EXTENDED_TOPOLOGY || topologyType == ALL_TOPOLOGIES) {
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Triangle<4> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Triangle<6> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Quadrilateral<8> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Quadrilateral<9> >() ) );
      }      
    }
    if( cellType == NONSTANDARD_CELL || cellType == ALL_CELLS){
      if(topologyType == BASE_TOPOLOGY || topologyType == ALL_TOPOLOGIES) {
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::ShellLine<2> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Beam<2> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Pentagon<5> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Hexagon<6> >() ) );
      }
      if(topologyType == EXTENDED_TOPOLOGY || topologyType == ALL_TOPOLOGIES) {
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::ShellLine<3> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Beam<3> >() ) );  
      }      
    }
  } // dim 2
  
  
  if((cellDim == 3) || (cellDim == 4)) {
    if( cellType == STANDARD_CELL || cellType == ALL_CELLS){
      if(topologyType == BASE_TOPOLOGY || topologyType == ALL_TOPOLOGIES) {
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Tetrahedron<4> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Hexahedron<8> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Pyramid<5> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Wedge<6> >() ) );
      }
      if(topologyType == EXTENDED_TOPOLOGY || topologyType == ALL_TOPOLOGIES) {
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Tetrahedron<8> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Tetrahedron<10> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Tetrahedron<11> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Hexahedron<20> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Hexahedron<27> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Pyramid<13> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Pyramid<14> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Wedge<15> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::Wedge<18> >() ) );
      }      
    }
    if( cellType == NONSTANDARD_CELL || cellType == ALL_CELLS){
      // Predefined Polyhedrons should  go here
      if(topologyType == BASE_TOPOLOGY || topologyType == ALL_TOPOLOGIES) {
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::ShellTriangle<3> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::ShellQuadrilateral<4> >() ) );
      }
      if(topologyType == EXTENDED_TOPOLOGY || topologyType == ALL_TOPOLOGIES) {
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::ShellTriangle<6> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::ShellQuadrilateral<8> >() ) );
        topologies.push_back( CellTopology( shards::getCellTopologyData<shards::ShellQuadrilateral<9> >() ) );
      }      
    }
  } // dim 3    
} // getTopologies


int isPredefinedCell(const CellTopology& cell) {
  
  switch(cell.getKey() ) {
    case Node::key:
    case Particle::key:
    case Line<2>::key:
    case Line<3>::key:
    case ShellLine<2>::key:
    case ShellLine<3>::key:
    case Beam<2>::key:
    case Beam<3>::key:
      
    case Triangle<3>::key:
    case Triangle<4>::key:
    case Triangle<6>::key:
    case ShellTriangle<3>::key:
    case ShellTriangle<6>::key:
      
    case Quadrilateral<4>::key:
    case Quadrilateral<8>::key:
    case Quadrilateral<9>::key:
    case ShellQuadrilateral<4>::key:
    case ShellQuadrilateral<8>::key:
    case ShellQuadrilateral<9>::key:
      
    case Tetrahedron<4>::key:
    case Tetrahedron<8>::key:
    case Tetrahedron<10>::key:
    case Tetrahedron<11>::key:
      
    case Hexahedron<8>::key:
    case Hexahedron<20>::key:
    case Hexahedron<27>::key:
      
    case Pyramid<5>::key:
    case Pyramid<13>::key:
    case Pyramid<14>::key:
      
    case Wedge<6>::key:
    case Wedge<15>::key:
    case Wedge<18>::key:
      
    case Pentagon<5>::key:
    case Hexagon<6>::key:
      return 1;
      
    default:
      return 0;
  }
  
}


} // namespace shards


