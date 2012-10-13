#ifndef _snl_fei_MapContig_hpp_
#define _snl_fei_MapContig_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

namespace snl_fei {

/** fei implementation class, many std::map characteristics but optimized
 for key-values that are contiguous */
template<typename VAL_TYPE>
class MapContig {
 public:
  /** constructor */
  MapContig(int firstKey, int lastKey);
  /** constructor */
  MapContig(const MapContig<VAL_TYPE>& src);
  /** destructor */
  virtual ~MapContig();

  /** alias for key_type */
  typedef int key_type;

  /** alias for mapped_type */
  typedef VAL_TYPE mapped_type;

  /** alias for value_type */
  typedef typename std::pair<int,VAL_TYPE> value_type;

  /** iterator class */
  class iterator {
   public:
    /** constructor */
    iterator() : offset_(-1), mapPtr_(0) {}

    /** constructor */
    iterator(int offset,
	     MapContig<VAL_TYPE>* mapPtr)
      : offset_(offset), mapPtr_(mapPtr)
    {
    }

    /** destructor */
    virtual ~iterator() {}

    /** operator++ */
    iterator& operator++()
    {
      if (!mapPtr_) return(*this);
      int len = mapPtr_->len_;
      int* keysPtr = mapPtr_->keysPtr_;
      int first = mapPtr_->first_;
      if (offset_ < len) {
        ++offset_;
	while(offset_ < len) {
          if (keysPtr[offset_] >= first) break;
          ++offset_;
        }
      }
      return(*this);
    }

    /** operator== */
    bool operator==(const iterator& rhs)
    {
      return( offset_ == rhs.offset_);
    }

    /** operator!= */
    bool operator!=(const iterator& rhs)
    {
      return( offset_ != rhs.offset_ );
    }

    /** operator* */
    value_type operator*()
    {
      if (!mapPtr_) return(value_type(0,0));

      if (offset_ == mapPtr_->len_) return(value_type(0,0));

      return(value_type(mapPtr_->keysPtr_[offset_],mapPtr_->valuesPtr_[offset_]));
    }

    /** operator= */
    iterator& operator=(const iterator& src)
    {
      offset_ = src.offset_;
      mapPtr_ = src.mapPtr_;
      return(*this);
    }

    /** ugh, public data member... */
    int offset_;
   private:
    MapContig<VAL_TYPE>* mapPtr_;
  };//class iterator

  /** iterator pointing to beginning of MapContig contents */
  iterator begin();
  /** iterator pointing one past the end of MapContig contents */
  iterator& end();

  /** insert */
  std::pair<iterator,bool> insert(value_type val);

  /** insert */
  iterator insert(iterator& pos, value_type val);

  /** find */
  iterator find(int key);

  /** lower_bound */
  iterator lower_bound(int key);

  /** size */
  int size() const;

 private:
  friend class iterator;

  std::vector<int> keys_;
  int* keysPtr_;
  iterator m_end_;
  std::vector<VAL_TYPE> values_;
  VAL_TYPE* valuesPtr_;
  int first_;
  int len_;
}; //class MapContig

template<typename VAL_TYPE>
MapContig<VAL_TYPE>::MapContig(int firstKey, int lastKey)
  : keys_(lastKey-firstKey+1),
    m_end_(),
    values_(lastKey-firstKey+1),
    first_(firstKey),
    len_(lastKey-firstKey+1)
{
  keysPtr_ = &keys_[0];
  for(int i=0; i<len_; ++i) {
    keysPtr_[i] = firstKey+i;
  }
  valuesPtr_ = &values_[0];
  len_ = keys_.size();
  m_end_ = iterator(len_, this);
}

template<typename VAL_TYPE>
MapContig<VAL_TYPE>::MapContig(const MapContig<VAL_TYPE>& src)
 : keys_(src.keys_),
  m_end_(),
  values_(src.values_),
  first_(src.first_),
  len_(src.len_)
{
  keysPtr_ = &keys_[0];
  valuesPtr_ = &values_[0];
  m_end_ = iterator(len_, this);
}

template<typename VAL_TYPE>
MapContig<VAL_TYPE>::~MapContig()
{
}

template<typename VAL_TYPE>
inline typename MapContig<VAL_TYPE>::iterator MapContig<VAL_TYPE>::begin()
{
  return( iterator(0, this) );
}

template<typename VAL_TYPE>
inline typename MapContig<VAL_TYPE>::iterator& MapContig<VAL_TYPE>::end()
{
  return( m_end_ );
}

template<typename VAL_TYPE>
inline std::pair<typename MapContig<VAL_TYPE>::iterator,bool>
 MapContig<VAL_TYPE>::insert(typename MapContig<VAL_TYPE>::value_type val)
{
  int localkey = val.first - first_;
  if (localkey < 0 || localkey >= len_) {
    return( std::pair<iterator,bool>(m_end_, false) );
  }

  valuesPtr_[localkey] = val.second;

  return( std::pair<iterator,bool>(iterator(localkey, this),true));
}

template<typename VAL_TYPE>
inline typename MapContig<VAL_TYPE>::iterator
 MapContig<VAL_TYPE>::insert(typename MapContig<VAL_TYPE>::iterator& pos,
			     typename MapContig<VAL_TYPE>::value_type val)
{
  int offset = pos.offset_;
  if (offset < 0 || offset >=len_ || pos == m_end_) {
    offset = val.first - first_;
    if (offset < 0 || offset >= len_) {
      return(m_end_);
    }
  }

  valuesPtr_[offset] = val.second;

  return( iterator(offset, this) );
}

template<typename VAL_TYPE>
inline typename MapContig<VAL_TYPE>::iterator MapContig<VAL_TYPE>::find(int key)
{
  int localkey = key - first_;
  if (localkey < 0 || localkey >= len_) {
    return( m_end_ );
  }

  return(iterator(localkey, this));
}

template<typename VAL_TYPE>
inline typename MapContig<VAL_TYPE>::iterator MapContig<VAL_TYPE>::lower_bound(int key)
{
  int localkey = key - first_;
  if (localkey < 0 || localkey >= len_) {
    return( m_end_ );
  }

  return(iterator(localkey, this));
}

template<typename VAL_TYPE>
int MapContig<VAL_TYPE>::size() const
{
  return(len_);
}

}//namespace snl_fei

#endif

