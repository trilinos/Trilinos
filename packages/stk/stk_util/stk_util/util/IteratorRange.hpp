// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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

#ifndef stk_util_util_IteratorRange_hpp
#define stk_util_util_IteratorRange_hpp

#include <utility>
#include <iterator>

namespace stk {


template< class IterType ,
          class IterCategory =
             typename std::iterator_traits< IterType >::iterator_category >
class IteratorRange ;

//----------------------------------------------------------------------
// Specialized for random access iterators, others TBD.

///
/// @addtogroup util_module
/// @{

/** @class IteratorRange
 *  Iterate a span of a container defined by begin and end iterators.
 *  Provides forward iterator and const container-like functionality.
 */
template< class IterType >
class IteratorRange< IterType , std::random_access_iterator_tag > 
  : public std::pair< IterType , IterType >
{
private:
  typedef std::pair< IterType , IterType > Pair ;
  typedef IteratorRange< IterType , std::random_access_iterator_tag > Self ;
  typedef std::iterator_traits< IterType > Traits ;
public:

  //--------------------------------

  typedef          IterType                iterator ;
  typedef typename Traits::value_type      value_type ;
  typedef typename Traits::pointer         pointer ;
  typedef typename Traits::reference       reference ;
  typedef typename Traits::difference_type difference_type ;
  typedef          size_t                  size_type ;

  //--------------------------------

  /** @brief Destructor
  */
  ~IteratorRange() {}

  /** Constructor */
  IteratorRange() : Pair() { Pair::second = Pair::first ; }

  IteratorRange( const Self & rhs ) : Pair( rhs ) {}

  IteratorRange( const Pair & rhs ) : Pair( rhs ) {}

  Self & operator = ( const Self & rhs )
    { Pair::first = rhs.first ; Pair::second = rhs.second ; return *this ; }

  Self & operator = ( const Pair & rhs )
    { Pair::first = rhs.first ; Pair::second = rhs.second ; return *this ; }

  //--------------------------------

  bool operator == ( const Self & rhs ) const
    { return Pair::first == rhs.first && Pair::second == rhs.second ; }

  bool operator != ( const Self & rhs ) const
    { return Pair::first != rhs.first || Pair::second != rhs.second ; }

  bool operator == ( const Pair & rhs ) const
    { return Pair::first == rhs.first && Pair::second == rhs.second ; }

  bool operator != ( const Pair & rhs ) const
    { return Pair::first != rhs.first || Pair::second != rhs.second ; }

  //--------------------------------

  Self & operator ++ () { ++ Pair::first ; return *this ; }

  Self operator ++ (int) { Self tmp(*this); ++ Pair::first ; return tmp ; }

  reference operator * ()  const { return * Pair::first ; }
  pointer   operator -> () const { return & * Pair::first ; }

  //--------------------------------
  // Container-like functionality for random access iterators.

  reference front() const { return * Pair::first ; }
  reference back()  const { return  Pair::second[-1] ; }

  iterator begin() const { return  Pair::first ; }
  iterator end()   const { return  Pair::second ; }

  template<class Iterator>
  IteratorRange( Iterator i , Iterator e ) : Pair(i,e) {}

  template<class Container>
  explicit
  IteratorRange( const Container & c ) : Pair( c.begin() , c.end() ) {}

  template<class Container>
  explicit
  IteratorRange( Container & c ) : Pair( c.begin() , c.end() ) {}

  bool empty () const { return ! ( Pair::first < Pair::second ) ; }

  reference operator [] ( size_t n ) const { return Pair::first[n] ; }

  size_t size() const
    {
      const difference_type d = std::distance( Pair::first , Pair::second );
      return d < 0 ? 0 : static_cast<size_t>(d) ;
    }
};

///
/// @}
///

} // namespace stk

#endif

