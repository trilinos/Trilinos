/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
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
/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   June 2008
 */

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <Dimension.hpp>

namespace phdmesh {

DimTag::~DimTag() {}

const char * DimTag::name() const
{ static const char n[] = "DimTag" ; return n ; }

std::string DimTag::to_string( unsigned size , unsigned index ) const
{
  std::ostringstream tmp ;

  if ( size <= index ) {
    tmp << "DimTag::to_string( " << size << " , " << index << " ) ERROR" ;
    throw std::runtime_error( tmp.str() );
  }

  tmp << index ;

  return tmp.str();
}

unsigned DimTag::to_index( unsigned size , const std::string & label ) const
{
  int index = size ? atoi( label.c_str() ) : 0 ;

  if ( index < 0 || ((int) size ) <= index ) {
    std::ostringstream tmp ;
    tmp << "DimTag::to_index( " << size << " , " << label << " ) ERROR" ;
    throw std::runtime_error( tmp.str() );
  }

  return (unsigned) index ;
}

//----------------------------------------------------------------------

unsigned DimBase::fortran_map( const unsigned * indices ) const
{
  unsigned offset = 0 ;

  switch( m_rank ) {
  case 8 : offset += indices[7] * m_stride[6] ;
  case 7 : offset += indices[6] * m_stride[5] ;
  case 6 : offset += indices[5] * m_stride[4] ;
  case 5 : offset += indices[4] * m_stride[3] ;
  case 4 : offset += indices[3] * m_stride[2] ;
  case 3 : offset += indices[2] * m_stride[1] ;
  case 2 : offset += indices[1] * m_stride[0] ;
  case 1 : offset += indices[0] ;
  }

  return offset ;
}

unsigned DimBase::natural_map( const unsigned * indices ) const
{
  unsigned offset = 0 ;

  switch( m_rank ) {
  case 8 : offset += indices[ m_rank - 8 ] * m_stride[6] ;
  case 7 : offset += indices[ m_rank - 7 ] * m_stride[5] ;
  case 6 : offset += indices[ m_rank - 6 ] * m_stride[4] ;
  case 5 : offset += indices[ m_rank - 5 ] * m_stride[3] ;
  case 4 : offset += indices[ m_rank - 4 ] * m_stride[2] ;
  case 3 : offset += indices[ m_rank - 3 ] * m_stride[1] ;
  case 2 : offset += indices[ m_rank - 2 ] * m_stride[0] ;
  case 1 : offset += indices[ m_rank - 1 ] ;
  }

  return offset ;
}

//----------------------------------------------------------------------

void DimBase::fortran_inv( unsigned offset , unsigned * indices ) const
{
  switch( m_rank ) {
  case 8 : indices[7] = offset / m_stride[6] ; offset %= m_stride[6] ;
  case 7 : indices[6] = offset / m_stride[5] ; offset %= m_stride[5] ;
  case 6 : indices[5] = offset / m_stride[4] ; offset %= m_stride[4] ;
  case 5 : indices[4] = offset / m_stride[3] ; offset %= m_stride[3] ;
  case 4 : indices[3] = offset / m_stride[2] ; offset %= m_stride[2] ;
  case 3 : indices[2] = offset / m_stride[1] ; offset %= m_stride[1] ;
  case 2 : indices[1] = offset / m_stride[0] ; offset %= m_stride[0] ;
  case 1 : indices[0] = offset ;
  }
}

void DimBase::natural_inv( unsigned offset , unsigned * indices ) const
{
  switch( m_rank ) {
  case 8 : indices[ m_rank - 8 ] = offset / m_stride[6]; offset %= m_stride[6];
  case 7 : indices[ m_rank - 7 ] = offset / m_stride[5]; offset %= m_stride[5];
  case 6 : indices[ m_rank - 6 ] = offset / m_stride[4]; offset %= m_stride[4];
  case 5 : indices[ m_rank - 5 ] = offset / m_stride[3]; offset %= m_stride[3];
  case 4 : indices[ m_rank - 4 ] = offset / m_stride[2]; offset %= m_stride[2];
  case 3 : indices[ m_rank - 3 ] = offset / m_stride[1]; offset %= m_stride[1];
  case 2 : indices[ m_rank - 2 ] = offset / m_stride[0]; offset %= m_stride[0];
  case 1 : indices[ m_rank - 1 ] = offset ;
  }
}

//----------------------------------------------------------------------

bool DimBase::fortran_valid( const unsigned * indices ) const
{
  bool result = true ;

  switch( m_rank ) {
  case 8 : result = result && indices[7] * m_stride[6] < m_stride[7] ;
  case 7 : result = result && indices[6] * m_stride[5] < m_stride[6] ;
  case 6 : result = result && indices[5] * m_stride[4] < m_stride[5] ;
  case 5 : result = result && indices[4] * m_stride[3] < m_stride[4] ;
  case 4 : result = result && indices[3] * m_stride[2] < m_stride[3] ;
  case 3 : result = result && indices[2] * m_stride[1] < m_stride[2] ;
  case 2 : result = result && indices[1] * m_stride[0] < m_stride[1] ;
  case 1 : result = result && indices[0]               < m_stride[0] ;
  }

  return result ;
}

bool DimBase::natural_valid( const unsigned * indices ) const
{
  bool result = true ;

  switch( m_rank ) {
  case 8 : result = result && indices[m_rank-8] * m_stride[6] < m_stride[7] ;
  case 7 : result = result && indices[m_rank-7] * m_stride[5] < m_stride[6] ;
  case 6 : result = result && indices[m_rank-6] * m_stride[4] < m_stride[5] ;
  case 5 : result = result && indices[m_rank-5] * m_stride[3] < m_stride[4] ;
  case 4 : result = result && indices[m_rank-4] * m_stride[2] < m_stride[3] ;
  case 3 : result = result && indices[m_rank-3] * m_stride[1] < m_stride[2] ;
  case 2 : result = result && indices[m_rank-2] * m_stride[0] < m_stride[1] ;
  case 1 : result = result && indices[m_rank-1]               < m_stride[0] ;
  }

  return result ;
}

//----------------------------------------------------------------------

void DimBase::fortran_size( unsigned * sizes ) const
{
  switch( m_rank ) {
  case 8 : sizes[7] = m_stride[7] / m_stride[6] ;
  case 7 : sizes[6] = m_stride[6] / m_stride[5] ;
  case 6 : sizes[5] = m_stride[5] / m_stride[4] ;
  case 5 : sizes[4] = m_stride[4] / m_stride[3] ;
  case 4 : sizes[3] = m_stride[3] / m_stride[2] ;
  case 3 : sizes[2] = m_stride[2] / m_stride[1] ;
  case 2 : sizes[1] = m_stride[1] / m_stride[0] ;
  case 1 : sizes[0] = m_stride[0] ;
  }
}

void DimBase::natural_size( unsigned * sizes ) const
{
  switch( m_rank ) {
  case 8 : sizes[m_rank-8] = m_stride[7] / m_stride[6] ;
  case 7 : sizes[m_rank-7] = m_stride[6] / m_stride[5] ;
  case 6 : sizes[m_rank-6] = m_stride[5] / m_stride[4] ;
  case 5 : sizes[m_rank-5] = m_stride[4] / m_stride[3] ;
  case 4 : sizes[m_rank-4] = m_stride[3] / m_stride[2] ;
  case 3 : sizes[m_rank-3] = m_stride[2] / m_stride[1] ;
  case 2 : sizes[m_rank-2] = m_stride[1] / m_stride[0] ;
  case 1 : sizes[m_rank-1] = m_stride[0] ;
  }
}

//----------------------------------------------------------------------

void print( std::ostream & s ,
            const DimBase & dim ,
            const DimTag ** tags ,
            const bool is_dim_natural )
{
  unsigned sizes[ DimMaxRank ];

  if ( is_dim_natural ) { s << "DimNatural<" ; dim.natural_size( sizes ); }
  else                  { s << "DimFortran<" ; dim.fortran_size( sizes ); }

  for ( unsigned i = 0 ; i < dim.rank() ; ++i ) {
    if ( i ) { s << "," ; }
    s << tags[i]->name();
  }

  s << ">(" ;

  for ( unsigned i = 0 ; i < dim.rank() ; ++i ) {
    if ( i ) { s << "," ; }
    s << sizes[i] ;
  }

  s << ")" ;
}


//----------------------------------------------------------------------

}


