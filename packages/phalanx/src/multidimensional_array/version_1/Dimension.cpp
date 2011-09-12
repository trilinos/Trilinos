/*
// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

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


