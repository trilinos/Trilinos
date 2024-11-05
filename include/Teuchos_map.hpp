// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_MAP_H
#define TEUCHOS_MAP_H

/*! \file Teuchos_map.hpp
    \brief Provides std::map class for deficient platforms.
*/

#include "Teuchos_ConfigDefs.hpp"
#ifndef TFLOP
#  include <map>
#endif // NOT TFLOP

/*! \class Teuchos::map
    \brief This class creates a basic std::map object for platforms where the std::map is
        deficient, otherwise the std::map is injected into the Teuchos namespace.

    \note
    <ol>
        <li> The std::map is an integral part of Teuchos::ParameterList and Teuchos::CommandLineProcessor.
        <li> Currently the basic std::map class is only used on ASCI Red (Janus).
    </ol>
*/

namespace Teuchos {

#ifdef TFLOP

template<class Key, class T>
class std::map {
public:
  typedef Key key_type;
  typedef T mapped_type;
  typedef std::pair<Key,T>  value_type;
  typedef std::list<value_type>  list_t;
  typedef typename list_t::iterator  iterator;
  typedef typename list_t::const_iterator  const_iterator;

  //! @name Constructor/Destructor.
  //@{

    //! Default Constructor
    std::map() {}

    //! Copy Constructor
    std::map( const std::map<Key,T>& map_in ) : list_( map_in.list_ ) {}

    //! Destructor
    virtual ~std::map() {}
  //@}

  //! @name Accessor methods.
  //@{

    //! Return an iterator that points to the first pair in the std::map.
    iterator begin() { return list_.begin(); }

    //! Return a const iterator that points to the first pair in the std::map.
    const_iterator begin() const { return list_.begin(); }

    //! Return an iterator that points to the last pair in the std::map.
    iterator end() { return list_.end(); }

    //! Return a const iterator that points to the last pair in the std::map.
    const_iterator end() const { return list_.end(); }

    //! Return a reference to the mapped value that belongs to the key \c k.
    /*! \param k - The key for which data should be retrieved.
        If this key doesn't exist then the key is inserted into the std::map and a
        reference to the mapped value is returned.
    */
    mapped_type& operator[]( const key_type& k )
    {
      iterator itr = find(k);
      if(itr != end()) return (*itr).second;
      list_.push_back( value_type( k, T() ) );
      return list_.back().second;
    }
  //@}

  //! @name Search methods.
  //@{

    //! Locate element in the std::map with key_type \c k.
    /*! \param k - The key for which an iterator should be returned.
        \return An iterator that points to the element with key_type \c k, else
        return end().
    */
    iterator find(const key_type& k)
    {
      for( iterator itr = begin(); itr != end(); ++itr ) {
        if( (*itr).first == k ) {
          return itr;
        }
      }
      return end();
    }

    //! Locate element in the std::map with key_type \c k.
    /*! \param k - The key for which a constant iterator should be returned.
        \return A constant iterator that points to the element with key_type \c k, else
        return end().
    */
    const_iterator find(const key_type& k) const
    {
      for( const_iterator itr = begin(); itr != end(); ++itr ) {
        if( (*itr).first == k ) {
          return itr;
        }
      }
      return end();
    }

    bool empty() const { return list_.empty(); }

  //@}
private:
  list_t list_;
};

#else

using std::map;

#endif

} // namespace Teuchos

#endif // TEUCHOS_MAP_H
