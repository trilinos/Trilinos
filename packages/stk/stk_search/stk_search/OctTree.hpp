/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_search_OctTree_hpp
#define stk_search_OctTree_hpp

#include <limits>
#include <iosfwd>
#include <stk_util/util/SimpleArrayOps.hpp>
#include <stk_util/util/StaticAssert.hpp>

namespace stk {

class OctTreeKey ;

}

namespace std {

ostream & operator << ( ostream & , const stk::OctTreeKey & otk);

}

namespace stk {

/** \ingroup util_module
 *  \class  OctTreeKey
 *  \brief  Type for identifying a particular node in a hierarchical oct-tree
 */
class OctTreeKey {
public:
  /** \brief  Integer type for the key components and indices */
  typedef unsigned value_type ;
  enum { MaxDepth     = 16 };
  enum { MaskIndex    = 0x0f };
  enum { BitsPerIndex = 4 };
  enum { BitsPerWord  = std::numeric_limits<value_type>::digits };
  enum { IndexPerWord = BitsPerWord / BitsPerIndex };
  enum { NWord        = MaxDepth / IndexPerWord };
  enum { OKBits = StaticAssert< 0 == BitsPerWord % BitsPerIndex >::OK };
  enum { OKWord = StaticAssert< 0 == MaxDepth    % IndexPerWord >::OK };
private:
  value_type m_value[ NWord ];
public:

  /** \brief  Default construct to the root node */
  OctTreeKey() : m_value()
    { Copy<NWord>( m_value , 0u ); }

  /** \brief  Copy construct */
  OctTreeKey( const OctTreeKey & k ) : m_value()
    { Copy<NWord>( m_value , k.m_value ); }

  /** \brief  Assignment */
  OctTreeKey & operator = ( const OctTreeKey & k )
    { Copy<NWord>( m_value , k.m_value ); return *this ; }

  /** \brief Query depth of this key */
  unsigned depth() const ;

  /** \brief  Index of the key at the depth [1..8]
   *          A zero value indicates it is not defined at that depth.
   */
  unsigned index( const unsigned Depth ) const ;

  /** \brief  Clear index at depth */
  OctTreeKey & clear_index( const unsigned Depth );

  /** \brief  Set index at depth */
  OctTreeKey & set_index( const unsigned Depth , const unsigned Index);

  /** \brief  Query raw value */
  const value_type * value() const { return m_value ; }

  /** \brief  Set raw value */
  OctTreeKey & set_value( const value_type * val);

  /** \brief  Intersects if either key contains the other. */
  bool intersect( const OctTreeKey & k ) const ;
};

/** \brief  Equality test */
inline bool operator == ( const OctTreeKey & l, const OctTreeKey & r ) 
  { return Compare<OctTreeKey::NWord>::equal( l.value() , r.value() ); }

/** \brief  Inequality test */
inline bool operator != ( const OctTreeKey & l, const OctTreeKey & r ) 
  { return Compare<OctTreeKey::NWord>::not_equal( l.value() , r.value() ); }

/** \brief  Comparison using depth-first full ordering */ 
inline bool operator < ( const OctTreeKey & l, const OctTreeKey & r ) 
  { return Compare<OctTreeKey::NWord>::less( l.value() , r.value()) ; }

//----------------------------------------------------------------------

/** \ingroup util_module
 *  \brief  Generate a 3D Hilbert space filling curve oct-tree key from
 *          an integer XYZ coordinate.
 *  \param Depth  Depth for the generated key.
 *  \param coord  XYZ coordinates in the range [0u:~0u].
 */
OctTreeKey hsfc3d( const unsigned Depth , const unsigned * const coord );

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template<unsigned Depth> struct OctTreeSize ;

template<>
struct OctTreeSize<0>
{ enum { value = 1 }; };

template<unsigned Depth>
struct OctTreeSize
{
  enum { MaxDepth = 10 , N = Depth }; // Size representable by an unsigned int

  enum { OK = StaticAssert< N <= MaxDepth >::OK };

  enum { value = 1 + 8 * OctTreeSize<Depth-1>::value };
};

/** \brief  Number of nodes in an oct-tree of a given depth */
unsigned oct_tree_size( const unsigned Depth );

/** \brief  Offset of a oct-tree node in a dense tree of a given depth. */
unsigned oct_tree_offset( const unsigned Depth , const OctTreeKey & k);
}

#endif

