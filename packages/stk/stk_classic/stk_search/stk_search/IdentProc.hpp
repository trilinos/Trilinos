/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_search_IdentProc_hpp
#define stk_search_IdentProc_hpp

#include <ostream>
#include <iomanip>

#include <stk_search/SearchTypes.hpp>
#include <stk_util/util/SimpleArrayOps.hpp>

namespace stk_classic {
namespace search {
namespace ident {

/**
 * @brief Class <b>template <class K, class P> IdentProc<K,P></b>.  Template parallel identification used for search, sorts first by identifier and then processor.
 *
 * K and P need to provide a default constructor, copy constructor, assignment operator, less than operator (<), and stream insertion operator (<<).
 */
template <class K = uint64_t, class P = unsigned>
struct IdentProc {

  typedef K Key;
  typedef P Proc;

  /**
   * Creates a new <b>IdentProc</b> instance.
   *
   */
  IdentProc()
    : ident(),
      proc()
  {}

  /**
   * Destroys a <b>IdentProc</b> instance.
   *
   */
  ~IdentProc()
  {}

  /**
   * Creates a new <b>IdentProc</b> instance.
   *
   * @param rhs			an <b>IdentProc</b> const ...
   *
   */
  IdentProc( const IdentProc & rhs )
    : ident(rhs.ident),
      proc(rhs.proc)
  {}

  /**
   * Creates a new <b>IdentProc</b> instance.
   *
   * @param i			an <b>unsigned int</b> ...
   *
   * @param p			an <b>unsigned int</b> ...
   *
   */
  IdentProc( Key i , Proc p )
    : ident(i),
      proc(p)
  {}

  /**
   * @brief Member function <b>=</b> ...
   *
   * @param rhs			an <b>IdentProc</b> const ...
   *
   * @return			an <b>IdentProc</b> ...
   */
  IdentProc & operator = ( const IdentProc & rhs ) {
    ident = rhs.ident ;
    proc = rhs.proc ;
    return *this ;
  }
  /**
  * @brief Member function <b>==</b> ...
  *
  * @param rhs			an <b>IdentProc</b> const ...
  *
  * @return			a <b>bool</b> ...
  */
  inline bool operator==( const IdentProc<K, P> & rhs ) const {
    return !( (ident < rhs.ident || rhs.ident < ident) ||
              (proc  < rhs.proc  || rhs.proc  < proc ) );
  }

  /**
  * @brief Member function <b><</b> ...
  *
  * @param rhs			an <b>IdentProc</b> const ...
  *
  * @return			a <b>bool</b> ...
  */
  inline bool operator<( const IdentProc<K, P> & rhs ) const {
    return ident < rhs.ident || (!(rhs.ident < ident) && proc < rhs.proc);
  }

  /**
  * @brief Member function <b>!=</b> ...
  *
  * @param rhs			an <b>IdentProc</b> const ...
  *
  * @return			a <b>bool</b> ...
  */
  inline bool operator!=( const IdentProc<K, P> &rhs ) const {
    return  ( (ident < rhs.ident || rhs.ident < ident) ||
              (proc  < rhs.proc  || rhs.proc  < proc ) );
  }

  /**
  * @brief Member function <b>></b> ...
  *
  * @param rhs			an <b>IdentProc</b> const ...
  *
  * @return			a <b>bool</b> ...
  */
  inline bool operator>( const IdentProc<K, P> &rhs ) const{
    return rhs < *this;
  }

  /**
  * @brief Member function <b><=</b> ...
  *
  * @param rhs			an <b>IdentProc</b> const ...
  *
  * @return			a <b>bool</b> ...
  */
  inline bool operator<=( const IdentProc<K, P> &rhs ) const {
    return !(rhs < *this);
  }

  /**
  * @brief Member function <b>>=</b> ...
  *
  * @param rhs			an <b>IdentProc</b> const ...
  *
  * @return			a <b>bool</b> ...
  */
  inline bool operator>=(const IdentProc<K, P> &rhs ) const {
    return !(*this < rhs);
  }


  Key           ident;          ///< Identifier
  Proc          proc;           ///< Processor
};


template <class K, class P>
std::ostream& operator<<(std::ostream &dout, const IdentProc<K, P> &ident_proc){
  dout << "id " << ident_proc.ident << ", proc " << ident_proc.proc;
  return dout;
}


} // namespace ident
} // namespace search
} // namespace stk_classic

#endif
