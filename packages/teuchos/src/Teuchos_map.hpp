#ifndef TEUCHOS_MAP_H
#define TEUCHOS_MAP_H

#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos {

//#ifdef TFLOP

template<class Key, class T>
class map {
public:
  typedef Key key_type;
  typedef T mapped_type;
  typedef std::pair<Key,T>  value_type;
  typedef std::list<value_type>  list_t;
  typedef typename list_t::iterator  iterator;
  typedef typename list_t::const_iterator  const_iterator;
  iterator begin() { return list_.begin(); }
  iterator end() { return list_.end(); }
  const_iterator begin() const { return list_.begin(); }
  const_iterator end() const { return list_.end(); }
  mapped_type& operator[]( const key_type& k )
    {
      iterator itr = find(k);
      if(itr != end()) return itr->second;
      list_.push_back( value_type( k, T() ) );
      return list_.back().second;
    }
  iterator find(const key_type& k)
    {
      for( iterator itr = begin(); itr != end(); ++itr ) {
        if( itr->first == k ) {
          return itr;
        }
      }
      return end();
    }
  const_iterator find(const key_type& k) const
    {
      for( const_iterator itr = begin(); itr != end(); ++itr ) {
        if( itr->first == k ) {
          return itr;
        }
      }
      return end();
    }
private:
  list_t list_;
};

//#else

//using std::map;

//#endif

} // namespace Teuchos

#endif // TEUCHOS_MAP_H
