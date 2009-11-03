
#ifndef _fei_ctg_set_hpp_
#define _fei_ctg_set_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2006 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_iostream.hpp>

#include <cstring>
#include <vector>
#include <fei_ArrayUtils.hpp>

namespace fei {

const int Set_end_val = -99999999;

/** A specialized container that mimics std::set in many ways.
  This set can only be used for integer types (such as short, int, long).

  This container is optimized for inserting/storing indices that have
  significant contiguous sections or ranges.

  Data is stored in an array in pairs, where each pair represents a
  contiguous range. The first item in a pair is the beginning of the range,
  and the second item in a pair is one greater than the end of the range.
  Example:
    Assume the data to be stored is 0, 1, 2, 3, 6, 7, 8, which forms two
    contiguous ranges 0 - 3 and 6 - 8. This container will store this data
    in an array like this: 0, 4, 6, 9.
*/
template<typename T>
class ctg_set {
  public:
  /** constructor */
  ctg_set(int alloc_incr=32)
    : dataptr_(0), len_(0), highwatermark_(0), alloc_incr_(alloc_incr) {}

  /** constructor */
  ctg_set(const ctg_set<T>& src)
    : dataptr_(0), len_(0),
      highwatermark_(src.highwatermark_), alloc_incr_(src.alloc_incr_)
    {
      if (highwatermark_>0) {
        expand_dataptr(highwatermark_);
        len_ = src.len_;
        for(int i=0; i<len_; ++i) dataptr_[i] = src.dataptr_[i];
      }
    }

  /** destructor */
  virtual ~ctg_set(){clear();}

  /** alias for template parameter */
  typedef T key_type;

  /** const_iterator */
  class const_iterator {
  public:
    /** constructor */
    const_iterator() : set_(0),
      val_(Set_end_val), limit_(Set_end_val), i_(0) {}

    /** constructor */
    const_iterator(const ctg_set<T>* _set,
		   const T& val,
		   int i)
      : set_(_set),
      val_(val), limit_(Set_end_val), i_(i)
      {
	if (set_ != 0) {
	  if (set_->len_ > 0) {
	    limit_ = set_->dataptr_[i+1]-1;
	  }
	}
      }

    /** constructor */
    const_iterator(const const_iterator& src)
      : set_(src.set_),
       val_(src.val_), limit_(src.limit_), i_(src.i_) {}

    /** destructor */
    virtual ~const_iterator(){}

    /** operator* */
    const T& operator*() const { return(val_); }

    /** operator++ */
    const_iterator& operator++()
      {
	if (val_ < limit_) {
	  ++val_;
	}
	else {
	  if (i_ < set_->len_-2) {
	    i_ += 2;
	    val_ = set_->dataptr_[i_];
	    limit_ = set_->dataptr_[i_+1]-1;
	  }
	  else {
	    val_ = Set_end_val;
	    limit_ = Set_end_val;
	  }
	}
	return(*this);
      }

    /** operator= */
    const_iterator& operator=(const const_iterator& src)
      {
	set_ = src.set_;
	i_ = src.i_;
	val_ = src.val_;
	limit_ = src.limit_;
	return(*this);
      }

    /** operator== */
    bool operator==(const const_iterator& rhs)
      {
	return(val_ == rhs.val_);
      }

    /** operator!= */
    bool operator!=(const const_iterator& rhs)
      {
	return(val_ != rhs.val_);
      }

  private:
    const ctg_set<T>* set_;
    T val_;
    T limit_;
    int i_;
  };

  /** iterator pointing to the beginning of ctg_set's contents */
  const_iterator begin() const
    {
      T val = Set_end_val;
      if (len_>0) {
	val = dataptr_[0];
      }
      return(const_iterator(this, val, 0));
    }

  /** iterator pointing one past the end of ctg_set's contents */
  static const_iterator end() { return(const_iterator(0, Set_end_val, 0)); }

  /** assignment operator= */
  ctg_set<T>& operator=(const ctg_set<T>& src)
    {
      highwatermark_ = src.highwatermark_;
      expand_dataptr(highwatermark_);
      len_ = src.len_;

      return(*this);
    }

  /** operator!= */
  bool operator!=(const ctg_set<T>& rhs)
  {
    if (len_ != rhs.len_) {
      return(true);
    }

    for(int i=0; i<len_; ++i) {
      if (dataptr_[i] != rhs.dataptr_[i]) {
        return(true);
      }
    }

    return(false);
  }

  /** Clear the contents of this set. (Get rid of the contents.)
   */
  int clear()
  {
    delete [] dataptr_;
    dataptr_ = 0;
    len_ = 0;
    highwatermark_ = 0;
    return(0);
  }

  /** Insert item in this set, if not already present.
      Return value is a pair with an iterator, and a bool indicating whether
      the insertion was made.
  */
  std::pair<const_iterator,bool> insert(const T& item)
  {
    if (len_ > highwatermark_) {
      FEI_COUT << "error"<<FEI_ENDL;
    }
    if (len_ < 1) {
      highwatermark_ = alloc_incr_;
      expand_dataptr(highwatermark_);
      dataptr_[0] = item;
      dataptr_[1] = item+1;
      len_ = 2;
      return(std::pair<const_iterator,bool>(const_iterator(this, item, 0),true));
    }

    int insertPoint = fei::lowerBound(item, dataptr_, len_);

    if (insertPoint < len_) {

      //dataptr_[insertPoint] is the first entry in dataptr_ that is not
      //less than item.
      //The possibilities are:
      //
      //1. dataptr_[insertPoint] equals item, so:
      //
      //     if (insertPoint is an even number) {
      //       return since item is present
      //     }
      //     else {
      //       expand the range below item to include item
      //       (by incrementing dataptr_[insertPoint])
      //       check whether dataptr_[insertPoint] == dataptr_[insertPoint+1]
      //       if so, merge the range at insertPoint-1 with the
      //       range at insertPoint+1
      //       return
      //     }
      //
      //2. dataptr_[insertPoint] is greater than item, so:
      //
      //     if (insertPoint is an even number) {
      //       if (item == dataptr_[insertPoint]-1) {
      //         expand the range at insertPoint to include item, by
      //         simply decrementing dataptr_[insertPoint]
      //       }
      //       else {
      //         insert a new range at insertPoint
      //       }
      //     }
      //     else {
      //       return since item is already present in the range at
      //       dataptr_[insertPoint-1]
      //     }
      //

      if (dataptr_[insertPoint] == item) {
	if (insertPoint%2 == 0) {
	  //insertPoint is an even number, so return since item is present
	  return(std::pair<const_iterator,bool>(const_iterator(this, item, insertPoint),false));
	}

	//Since insertPoint is an odd number, item lies just outside an existing
	//range so we simply need to add item to the range by incrementing
	//dataptr_[insertPoint].
	++dataptr_[insertPoint];

	//now check whether this newly expanded range should be merged
	//with the range above it
	if (insertPoint < len_-1) {
	  if (dataptr_[insertPoint] == dataptr_[insertPoint+1]) {
	    dataptr_[insertPoint] = dataptr_[insertPoint+2];
	    len_ -= 2;
	    int nmove=len_-insertPoint-1;
	    if (nmove > 0) {
	      T* dest = dataptr_+insertPoint+1;
	      T* src =  dest+2;
	      std::memmove(dest, src, nmove*sizeof(T));
	    }
	  }
	}

	return(std::pair<const_iterator,bool>(const_iterator(this, item, insertPoint-1),true));
      }
      else {
	//dataptr_[insertPoint] is greater than item.

	if (insertPoint%2 == 0) {
	  if (item == dataptr_[insertPoint]-1) {
	    --dataptr_[insertPoint];
	    return(std::pair<const_iterator,bool>(const_iterator(this, item, insertPoint),true));
	  }
	  else {
	    //insert a new range at insertPoint
	    int nmove = len_-insertPoint;
	    if (len_+2 > highwatermark_) {
	      highwatermark_ += alloc_incr_;
	      expand_dataptr(highwatermark_);
	    }
	    len_ += 2;
	    if (nmove > 0) {
	      T* dest = dataptr_+insertPoint+2;
	      T* src = dest - 2;
	      std::memmove(dest, src, nmove*sizeof(T));
	    }
	    dataptr_[insertPoint] = item;
	    dataptr_[insertPoint+1] = item+1;

	    return(std::pair<const_iterator,bool>(const_iterator(this, item, insertPoint),true));
	  }
	}
	else {
	  return(std::pair<const_iterator,bool>(const_iterator(this, item, insertPoint-1),false));
	}
      }
    }

    //If we get to here, insertPoint >= len_, meaning we need to append
    //a new range.
    if (len_+2 > highwatermark_) {
      highwatermark_ += alloc_incr_;
      expand_dataptr(highwatermark_);
    }
    len_ += 2;
    dataptr_[insertPoint] = item;
    dataptr_[insertPoint+1] = item+1;

    return(std::pair<const_iterator,bool>(const_iterator(this, item, insertPoint),true));
  }

  /** insert2 -- power-users only */
  void insert2(const T& item)
  {
    if (len_ < 1) {
      highwatermark_ = alloc_incr_;
      expand_dataptr(highwatermark_);
      dataptr_[0] = item;
      dataptr_[1] = item+1;
      len_ = 2;
      return;
    }

    int insertPoint = fei::lowerBound(item, dataptr_, len_);

    if (insertPoint < len_) {

      //dataptr_[insertPoint] is the first entry in dataptr_ that is not
      //less than item.
      //The possibilities are:
      //
      //1. insertPoint is an even number:
      //   dataptr_[insertPoint] is the start of an existing range.
      //   diff = dataptr_[insertPoint] - item;
      //   switch(diff) {
      //    case 0 : //  item in range, so return
      //    case 1 : //  item just below range, so expand range and return
      //    default: //  insert new range for item
      //   }
      //
      //2. insertPoint is an odd number:
      //   dataptr_[insertPoint] is the end of an existing range
      //   diff = dataptr_[insertPoint] - item;
      //   switch(diff) {
      //    case 0 : {
      //      // item just above range, so expand range
      //      // check whether range should be merged with range above
      //    }
      //    default: // item in range, so return
      //   }
      //

      if (insertPoint%2==0) {
        switch(dataptr_[insertPoint]-item) {
          case 0: break; //item in range
          case 1: {//expand range downwards
            --dataptr_[insertPoint];
            break;
          }
          default: {// insert new range for item
            //insert a new range at insertPoint
            int nmove = len_-insertPoint;
            if (len_+2 > highwatermark_) {
              highwatermark_ += alloc_incr_;
              expand_dataptr(highwatermark_);
            }
            len_ += 2;

            T* dest = dataptr_+insertPoint+2;
            T* src = dest - 2;
            std::memmove(dest, src, nmove*sizeof(T));

            dataptr_[insertPoint] = item;
            dataptr_[insertPoint+1] = item+1;
          }
        }
      }
      else {//insertPoint is odd number
        if (dataptr_[insertPoint] == item) {
          // item just above range, so expand range
          ++dataptr_[insertPoint];

          // check whether range should be merged with range above
          if (insertPoint < len_-1 &&
              dataptr_[insertPoint] == dataptr_[insertPoint+1]) {
            dataptr_[insertPoint] = dataptr_[insertPoint+2];
            len_ -= 2;
            int nmove=len_-insertPoint-1;
            if (nmove > 0) {
              T* dest = dataptr_+insertPoint+1;
              T* src =  dest+2;
              std::memmove(dest, src, nmove*sizeof(T));
            }
          }
        }//end if (dataptr_[insertPoint]==item)...
        //     else do nothing, item is already in existing range
      }

      return;
    } // end if (insertPoint < len_)...

    //If we get to here, insertPoint >= len_, meaning we need to append
    //a new range.
    if (len_+2 > highwatermark_) {
      highwatermark_ += alloc_incr_;
      expand_dataptr(highwatermark_);
    }
    dataptr_[insertPoint] = item;
    dataptr_[insertPoint+1] = item+1;
    len_ += 2;
  }

  /** insert2_dense_group -- not really implemented correctly */
  int insert2_dense_group(const T& starting_index, int group_size)
    {
      for(int i=0; i<group_size; ++i) {
	insert2(starting_index+i);
      }

      return(0);
    }

  /** find specified item in ctg_set. If not found, return iterator==end(). */
  const_iterator find(const T& item)
    {
      if (len_ < 1) {
	return(const_iterator(0, Set_end_val, 0));
      }

      int insertPoint = -1;
      int index = fei::binarySearch(item, dataptr_, len_, insertPoint);

      if (index < 0) {
	if (insertPoint%2==0) {
	  return(end());
	}
	else {
	  return(const_iterator(this, item, insertPoint-1));
	}
      }

      if (index%2==0) {
	return( const_iterator(this, item, index) );
      }

      return(const_iterator(0, Set_end_val, 0));
    }

  /** Copy contents of ctg_set into array of length len.
      The number of items copied into array is min(len, size()).
  */
  int copy_to_array(int len, T* items) const
    {
      const_iterator iter = begin(), iter_end = end();
      int offset = 0;
      for(; iter != iter_end; ++iter) {
	if (offset >= len) {
	  break;
	}

	items[offset++] = *iter;
      }

      return(0);
    }

  /** Copy contents of ctg_set into std::vector. */
  int copy_to_vector(std::vector<T>& items) const
    {
      int sz = size();
      items.resize(sz);
      T* itemsPtr = &(items[0]);
      return(copy_to_array(sz, itemsPtr));
    }

  /** return size of ctg_set. */
  int size() const
    {
      int setsize = 0;
      int offset = 0;
      while(offset<(len_-1)) {
	setsize += dataptr_[offset+1]-dataptr_[offset];
	offset += 2;
      }

      return(setsize);
    }

 private:
  void expand_dataptr(int newlen)
  {
    //on entry to this function, dataptr_ has length 'len_'.
    //we assume newlen is greater than len_.
    //after we create newptr with length 'newlen', we copy
    //len_ positions from dataptr_ to newptr.

    T* newptr = new T[newlen];
    for(int i=0; i<len_; ++i) newptr[i] = dataptr_[i];
    delete [] dataptr_;
    dataptr_ = newptr;
  }

  friend class const_iterator;
  T* dataptr_;
  int len_;
  int highwatermark_;
  int alloc_incr_;
};//class ctg_set

}//namespace fei

#endif

