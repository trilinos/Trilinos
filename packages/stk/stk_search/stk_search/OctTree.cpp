// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <sstream>
#include <ostream>
#include <stdexcept>
#include <stk_search/OctTree.hpp>

namespace std {

ostream & operator << ( ostream & os , const stk::OctTreeKey & otk )
{
  unsigned j = 0 ;

  while ( j < stk::OctTreeKey::MaxDepth ) {
    os << otk.index(++j);
  }
  
  return os ;
}

}

namespace stk {

enum { OK = StaticAssert< OctTreeKey::NWord == 2 >::OK };

namespace {
  void throw_index( const unsigned d , const unsigned i )
  {
    std::ostringstream msg ;
    msg << "OctTree.cpp: index[" << d << "] = " << i << " is out of range [1..8]" ;
    throw std::range_error( msg.str() );
  }
  
  void throw_depth( const unsigned min_depth, const unsigned depth, const unsigned max_depth )
  {
    std::ostringstream msg ;
    msg << "OctTree.cpp: depth = " << depth << " is out of range [" << min_depth << ".." << max_depth << "]" ;
    throw std::range_error( msg.str() );
  }
}
  
unsigned OctTreeKey::depth() const
{
  const unsigned which = m_value[1] ? 1 : 0 ;
  const unsigned val   = m_value[which] ;

  int d = IndexPerWord ;
  while ( d-- && ( val & ( MaskIndex << ( BitsPerIndex * d ) ) ) );
  return ( which + 1 ) * IndexPerWord - ( d + 1 );
}

unsigned OctTreeKey::index( const unsigned Depth ) const
{
  if ( Depth < 1 || MaxDepth < Depth ) { throw_depth( 1, Depth, MaxDepth ); }

  const unsigned which = ( Depth - 1 ) / IndexPerWord ;
  const unsigned shift = BitsPerWord -
                         BitsPerIndex * ( Depth % IndexPerWord ) ;

  return ( m_value[ which ] >> shift ) & MaskIndex ;
}

OctTreeKey & OctTreeKey::clear_index( const unsigned Depth )
{
  if ( Depth < 1 || MaxDepth < Depth ) { throw_depth( 1, Depth, MaxDepth ); }

  const value_type m = MaskIndex ;
  const unsigned which = ( Depth - 1 ) / IndexPerWord ;
  const unsigned shift = BitsPerWord -
                         BitsPerIndex * ( Depth % IndexPerWord ) ;
  const unsigned mask = ~( m << shift );

  m_value[ which ] &= mask ;

  return *this ;
}

OctTreeKey &
OctTreeKey::set_index( const unsigned Depth , const unsigned Index )
{
  if ( Depth < 1 || MaxDepth < Depth ) { throw_depth( 1, Depth, MaxDepth ); }
  if ( Index < 1 || 8        < Index ) { throw_index( Depth , Index ); }

  const value_type m = MaskIndex ;
  const unsigned which = ( Depth - 1 ) / IndexPerWord ;
  const unsigned shift = BitsPerWord -
                         BitsPerIndex * ( Depth % IndexPerWord ) ;

  ( m_value[which] &= ~( m << shift ) ) |= Index << shift ;

  return *this ;
}

OctTreeKey &
OctTreeKey::set_value( const unsigned * val )
{
  Copy<NWord>( m_value , 0u );

  for ( unsigned d = 1 ; d <= MaxDepth ; ++d ) {
    const unsigned which = ( d - 1 ) / IndexPerWord ;
    const unsigned shift = BitsPerWord - BitsPerIndex * ( d % IndexPerWord ) ;
    const unsigned index = ( val[ which ] >> shift ) & MaskIndex ;

    if ( 8 < index ) { throw_index( d , index ); }

    m_value[ which ] |= index << shift ;
  }
  return *this ;
}

bool OctTreeKey::intersect( const OctTreeKey & k ) const
{
  const unsigned which = m_value[0] != k.m_value[0] ? 0 : (
                         m_value[1] != k.m_value[1] ? 1 : 2 );

  bool result = which == 2 ;

  if ( ! result ) {
    const value_type lval =   m_value[which];
    const value_type rval = k.m_value[which];

    result = true ;
    for ( unsigned d = IndexPerWord ; result && d ; ) {
      --d ;
      const value_type mask = MaskIndex << ( BitsPerIndex * d );
      const value_type l = lval & mask ;
      const value_type r = rval & mask ;
      if ( l && r ) { result = l == r ; }
      else          { d = 0 ; }
    }
  }
  return result ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

const unsigned TreeMaxDepth = 10;
const unsigned tree_size[TreeMaxDepth+1] = {
  OctTreeSize< 0>::value ,
  OctTreeSize< 1>::value ,
  OctTreeSize< 2>::value ,
  OctTreeSize< 3>::value ,
  OctTreeSize< 4>::value ,
  OctTreeSize< 5>::value ,
  OctTreeSize< 6>::value ,
  OctTreeSize< 7>::value ,
  OctTreeSize< 8>::value ,
  OctTreeSize< 9>::value ,
  OctTreeSize<10>::value
};

}

unsigned oct_tree_size( const unsigned Depth )
{
  if ( TreeMaxDepth < Depth )
    { throw_depth( 0, Depth, TreeMaxDepth ); }

  return tree_size[ Depth ];
}

unsigned oct_tree_offset( const unsigned Depth , const OctTreeKey & k )
{
  if ( TreeMaxDepth < Depth )
    { throw_depth( 0, Depth, TreeMaxDepth ); }

  unsigned index = 0 ;

  unsigned d = Depth ;

  while ( d ) {
    --d ;
    const unsigned i = k.index( Depth - d );
    if ( i ) { index += 1 + ( i - 1 ) * tree_size[ d ] ; }
    else     { d = 0 ; } // Stop at a zero index
  }

  return index ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

OctTreeKey hsfc3d( const unsigned Depth , const unsigned * const coord )
{
  // Gray & Inverse Gray coding for 3D octree

/*
  const unsigned gray_fwd[8] = { 0 , 1 , 3 , 2 , 6 , 7 , 5 , 4 };
*/
  const unsigned gray_inv[8] = { 0 , 1 , 3 , 2 , 7 , 6 , 4 , 5 };

  // Three dimensional HSFC rotation data

  const unsigned hsfc3d_rotation_perm[8][3] =
    { /* 0 */ { 2 , 1 , 0 } ,
      /* 1 */ { 0 , 2 , 1 } ,
      /* 2 */ { 0 , 1 , 2 } ,
      /* 3 */ { 2 , 0 , 1 } ,
      /* 4 */ { 2 , 0 , 1 } ,
      /* 5 */ { 0 , 1 , 2 } ,
      /* 6 */ { 0 , 2 , 1 } ,
      /* 7 */ { 2 , 1 , 0 } };

  const unsigned hsfc3d_rotation_flip[8][3] =
      { /* 0 */ { 0 , 0 , 0 } ,
        /* 1 */ { 0 , 0 , 0 } ,
        /* 2 */ { 0 , 0 , 0 } ,
        /* 3 */ { 1 , 1 , 0 } ,
        /* 4 */ { 0 , 1 , 1 } ,
        /* 5 */ { 0 , 0 , 0 } ,
        /* 6 */ { 0 , 1 , 1 } ,
        /* 7 */ { 1 , 0 , 1 } };

  OctTreeKey key ;

  // Initial rotation

  unsigned axis[3] ;

  axis[0] = 0 << 1 ;
  axis[1] = 1 << 1 ;
  axis[2] = 2 << 1 ;

  for ( unsigned d = 1 ; d <= Depth ; ++d ) {

    const unsigned s = OctTreeKey::BitsPerWord - d ;

    // The ordinal for this depth

    const unsigned ord = gray_inv[
      (((( coord[ axis[0] >> 1 ] >> s ) ^ axis[0] ) & 01 ) << 0 ) |
      (((( coord[ axis[1] >> 1 ] >> s ) ^ axis[1] ) & 01 ) << 1 ) |
      (((( coord[ axis[2] >> 1 ] >> s ) ^ axis[2] ) & 01 ) << 2 ) ];

    key.set_index( d , ord + 1 );

    // Determine the recursive rotation for the next ordinal
    {
      const unsigned * const p = hsfc3d_rotation_perm[ ord ] ;
      const unsigned * const f = hsfc3d_rotation_flip[ ord ] ;

      unsigned tmp[3] ;

      tmp[0] = axis[0] ;
      tmp[1] = axis[1] ;
      tmp[2] = axis[2] ;

      axis[0] = tmp[ p[0] ] ^ f[0] ;
      axis[1] = tmp[ p[1] ] ^ f[1] ;
      axis[2] = tmp[ p[2] ] ^ f[2] ;
    }
  }

  return key ;
}


}


