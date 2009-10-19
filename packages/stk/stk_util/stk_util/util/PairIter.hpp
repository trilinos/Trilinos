#ifndef stk_util_util_PairIter_hpp
#define stk_util_util_PairIter_hpp

#include <utility>
#include <iterator>

namespace stk {


template< class IterType ,
          class IterCategory =
             typename std::iterator_traits< IterType >::iterator_category >
class PairIter ;

//----------------------------------------------------------------------
// Specialized for random access iterators, others TBD.

///
/// @addtogroup util_module
/// @{

/** @class PairIter
 *  Iterate a span of a container defined by begin and end iterators.
 *  Provides forward iterator and const container-like functionality.
 */
template< class IterType >
class PairIter< IterType , std::random_access_iterator_tag > 
  : public std::pair< IterType , IterType >
{
private:
  typedef std::pair< IterType , IterType > Pair ;
  typedef PairIter< IterType , std::random_access_iterator_tag > Self ;
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
  ~PairIter() {}

  /** Constructor */
  PairIter() : Pair() { Pair::second = Pair::first ; }

  PairIter( const Self & rhs ) : Pair( rhs ) {}

  PairIter( const Pair & rhs ) : Pair( rhs ) {}

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
  PairIter( Iterator i , Iterator e ) : Pair(i,e) {}

  template<class Container>
  explicit
  PairIter( const Container & c ) : Pair( c.begin() , c.end() ) {}

  template<class Container>
  explicit
  PairIter( Container & c ) : Pair( c.begin() , c.end() ) {}

  bool empty () const { return ! ( Pair::first < Pair::second ) ; }

  reference operator [] ( size_t n ) const { return Pair::first[n] ; }

  size_t size() const
    {
      const difference_type d = std::distance( Pair::first , Pair::second );
      return d < 0 ? 0 : (size_t) d ;
    }
};

///
/// @}
///

} // namespace stk

#endif

