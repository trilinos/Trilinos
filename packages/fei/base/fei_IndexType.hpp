/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_IndexType_hpp_
#define _fei_IndexType_hpp_

#include <fei_macros.hpp>

#include <map>
#include <algorithm>

#ifdef HAVE_FEI_BOOST
#include <boost/unordered_map.hpp>
#endif


namespace fei {

template< typename _Key, typename _Value>
class IndexType : public
#ifdef HAVE_FEI_BOOST
   boost::unordered_map<_Key, _Value>
#else
  std::map<_Key, _Value>
#endif
{
  typedef std::map<_Key, _Value> mapRep;
public:
#ifdef HAVE_FEI_BOOST
public:
  typedef boost::unordered_map<_Key, _Value> base;
  bool isStdMap() const {return false;}
  const mapRep& asMap(mapRep &map_rep) const {
    resyncToMap(map_rep);
    return map_rep;
  }
  mapRep& asMap(mapRep &map_rep) {
    resyncToMap(map_rep);
    return map_rep;
  }

  _Key getMaxKey() const {
    if (this->size() == 0 ) 
      return _Value(0);
    _Key k = this->begin()->first;
    typename base::const_iterator iter=++(this->begin()),end=this->end();
    for ( ; iter != end; ++iter) 
      k = std::max(k,iter->first);
    return k;
  }

  _Key getMinKey() const {
    if (this->size() == 0 ) 
      return _Value(0);
    _Key k = this->begin()->first;
    typename base::const_iterator iter=++(this->begin()),end=this->end();
    for ( ; iter != end; ++iter) 
      k = std::min(k,iter->first);
    return k;
  }

#else
  typedef std::map<_Key, _Value> base;
  bool isStdMap() const {return true;}
  mapRep &asMap(mapRep &rep) const{ return *static_cast<base *>(&rep);}
  const mapRep &asMap(const mapRep&rep) const { return *static_cast<const base *>( &rep);}

  _Key getMaxKey() const {
    IndexType<int,int>::const_iterator highest = this->end();
    if (this->size() > 0) 
      return (--highest)->first;
    return _Key(0);
  }

  _Key getMinKey() const {
    if (this->size() > 0) 
      return this->begin()->first;
    return _Key(0);
  }

#endif
  
  bool lightWeightCompare(const mapRep &map_rep) const {
    return (this->size() == map_rep.size());
  }

  bool heavyWeightCompare(const mapRep &map_rep) const {

    if (this->size() != map_rep.size())
      return false;
    typename base::const_iterator iter=this->begin(),end=this->end();
    for ( ; iter != end; ++iter) {
      typename mapRep::const_iterator lookup = map_rep.find(iter->first);
      if ( lookup == map_rep.end())
	return false;
      if (lookup->second  != iter->second)
	return false;
    }
    return true;
  }

  void resyncToMap( mapRep &map_rep) const {
    if ( !heavyWeightCompare(map_rep) ) {
      map_rep.clear();
      typename base::const_iterator iter=this->begin(), end=this->end();
      for ( ; iter != end; ++iter) 
	map_rep.insert(*iter);
    }
  }

  void resyncFromMap(const mapRep &map_rep) {
    if ( !heavyWeightCompare(map_rep) ) {
      this->clear();
      typename mapRep::const_iterator iter=map_rep.begin(), end=map_rep.end();
      for ( ; iter != end; ++iter) 
	this->insert(*iter);
    }
  }  
}; 
}
#endif
