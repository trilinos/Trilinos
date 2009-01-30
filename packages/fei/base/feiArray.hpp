#ifndef _feiArray_hpp_
#define _feiArray_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <cstring> //include cstring so we can use memcpy
#include "fei_fwd.hpp"
#include "fei_iostream.hpp"

template<class T>
struct lessthan {
  bool operator()(const T& x, const T& y) const { return( x < y ); }
};

namespace snl_fei {

enum Type_enum {
  INT = 0,
  DOUBLE = 1,
  OTHER = 2
};

template<typename T>
struct TypeTraits {

  static Type_enum get_type() { return OTHER; }
};

template<>
struct TypeTraits<int> {
  static Type_enum get_type() { return INT; }
};

template<>
struct TypeTraits<double> {
  static Type_enum get_type() { return DOUBLE; }
};

}//namespace snl_fei

//Note that the second template parameter 'class COMPARE' has a default
//value which is defined in the forward declaration for feiArray in fei_fwd.hpp

/** An array class which provides basic self-description,
   memory allocation/deletion, and searching/inserting/appending.
  
  This class provides management of an array, and keeps track of two lengths:
  an allocated length, and the length which has been requested by the caller.
  The allocated length will always be a multiple of an allocate-increment,
  which is a constructor argument that defaults to 16. The allocated length
  will be automatically expanded as necessary in order to accomodate the
  caller-requested length (requested in calls to resize, as well
  as the constructor argument).
  
  NOTE: upon allocation, this class does NOT initialize the new memory. The
  caller is responsible for that.

  NOTE ALSO: You can declare and use an 'feiArray of feiArrays', e.g.,
  feiArray<feiArray<int> > int_table;
  HOWEVER, you can't use the find, binarySearch or insertSorted functions on
  the first dimension of the table.
*/

template<typename T, class COMPARE >
class feiArray {
 public:
  /** Constructor
      @param length Defaults to 0 if not supplied
      @param incr Allocation increment, the 'chunk-size' by which memory is
      increased when new items are appended. Defaults to 16 if not supplied.
  */
  feiArray(int length=0, int incr=16);

  /** Constructor to wrap an feiArray around existing memory. When this 
      constructor is used, the reAllocate function will always return an error,
      and append/insert etc., will only work when 'length' is less than
      'allocLength'.
      @param length Amount of meaningful data already in place.
      @param allocLength Amount of memory allocated in 'data'.
      @param data Pointer to user-allocated array.
  */
  feiArray(int len, int allocLength, T* data);

  /** Copy constructor. */
  feiArray(const feiArray<T,COMPARE>& src);

  /** Destructor */
  virtual ~feiArray();

  /** Assignment operator. */
  feiArray<T,COMPARE>& operator=(const feiArray<T,COMPARE>& src);

  /** Assignment operator for initialization. Sets the value 'item' 
      throughout the array.
  */
  feiArray<T,COMPARE>& operator=(const T& item);

  /** Comparison operator. */
  bool operator==(const feiArray<T,COMPARE>& rhs);

  /** Not-equal operator. */
  bool operator!=(const feiArray<T,COMPARE>& rhs)
    {
      return(!( (*this) == rhs));
    }

  /** Return the length of the internal data array. This length
      is the amount of memory that has been requested by either the 'length'
      argument to the constructor, or by calls to the resize method. This
      length may be less than the actual allocated length of the array.
      This length will generally indicate how much of the allocated array
      is actually in use by the caller.
  */
  int length() const
    {
      return(length_);
    }

  int size() const
    {
      return(length_);
    }

  /** Return the actual allocated length of the internal data array, which
      will always be equal to or greater than the value returned by length().
  */
  int allocatedLength() const
    {
      return(allocatedLength_);
    }

  /** Return the allocation increment, which is the smallest amount by which
      the allocated length of the internal data array may grow.
  */
  int allocIncrement() const
    {
      return(allocIncrement_);
    }

  /** Alter the amount of memory that's available to the caller. 
      @param newLength If less than the currently available length, memory is
      not re-allocated, but the value returned by length() is adjusted
      downwards. If greater than the current length, then new memory is allocated
      as necessary, and any existing data is copied into the new memory. 
  */
  void resize(int newLength);

  /** Alter the amount of allocated memory. Doesn't change the value
      reported by length().
      @param newLength 
  */
  int reAllocate(int newAllocLength);


  /** A method for accessing the data as if it were a regular
      C array. 
      NOTE: This is bounds-checked access, with a corresponding performance
      penalty.
      @param offset 0-based offset, like regular C array access.
  */
  T& operator[](int offset);

  /** A method for accessing the data as if it were a regular
      C array. 
      NOTE: This is bounds-checked access, with a corresponding performance
      penalty.
      @param offset 0-based offset, like regular C array access.
  */
  const T& operator[](int offset) const;


  /** Return a raw pointer to the internal data array. Provides a high-
      performance method for accessing the data.
  */
  T* dataPtr() {return(data_);};

  /** Return a const raw pointer to the internal data array. Provides a high-
      performance method for accessing the data.
  */
  const T* dataPtr() const {return(data_);};

  /** Insert 'item' at position 'offset', regardless of whether this insertion
      causes the data to no longer be sorted. Expands allocated memory if
      necessary.
      @param item Entity of type T to be inserted.
      @param offset Position at which to insert 'item'. Should be in the
      range [0..length()] (i.e., offset may be 'too big' by one, which is 
      equivalent to calling append(item)).
  */
  void insert(const T& item, int offset);


  /** Append 'item' to the data array, expanding allocated memory if
      necessary. 
      @param item Entity of type T to be appended.
  */
  void append(const T& item);

  /** Perform a linear (exhaustive) search.
      @param item Value to be searched for.
      @return position Offset of the position at which it was found. If it was 
      not found, returns -1.
  */
  int find(const T& item) const;

 private:
  void deleteMemory();

  COMPARE cmp_;
  T* data_;
  int length_;
  int allocatedLength_;
  int allocIncrement_;
  bool userMemory_;
};

#include <fei_macros.hpp>

//------------------------------------------------------------------------------
template<class T,class COMPARE>
inline feiArray<T,COMPARE>& feiArray<T,COMPARE>::operator=(const T& item)
{
  for(int i=0; i<length_; i++) data_[i] = item;

  return(*this);
}

//------------------------------------------------------------------------------
template<typename T, class COMPARE>
void feiArray<T,COMPARE>::resize(int newLength)
{
  //
  //resize re-allocates the internal data array, if newLength is greater than
  //allocatedLength_. If newLength is smaller than allocatedLength_, we just
  //reset length_.
  //
  //If we are re-allocating the array to a larger size, any existing data is
  //copied into the new array.
  //
  if (newLength == length_) return;

  if (newLength < 0) {
    newLength = 0;
  }

  if (newLength <= allocatedLength_) {
    length_ = newLength;
    return;
  }

  if (allocIncrement_ == 0) {
    throw std::runtime_error("feiArray::resize, can't resize, allocIncrement_==0");
  }

  if (allocIncrement_ == 1) {
    allocatedLength_ = newLength;
  }
  else {
    int multiple = newLength/allocIncrement_;
    allocatedLength_ = (multiple+1)*allocIncrement_;
  }

  T* newData = new T[allocatedLength_];

  if (length_>0) {
    snl_fei::Type_enum t_type = snl_fei::TypeTraits<T>::get_type();
    if (t_type != snl_fei::OTHER) {
      std::memcpy(newData, data_, length_*sizeof(T));
    }
    else {
      T* ptr_src = data_;
      T* src_end = data_+length_;
      T* ptr_dest = newData;
      while(ptr_src < src_end) {
	*ptr_dest++ = *ptr_src++;
      }
    }
  }

  length_ = newLength;

  delete [] data_;
  data_ = newData;
}

//------------------------------------------------------------------------------
template<typename T, class COMPARE>
void feiArray<T,COMPARE>::insert(const T& item, int offset)
{
  if (offset < 0 || offset > length_) {
    throw std::runtime_error("feiArray::insert called with invalid offset");
  }

  if (offset == length_) {
    append(item);
    return;
  }

  ++length_;
  T* src = data_;
  bool alloc_new_memory = false;
  T* dest = data_;
  T* newData = NULL;

  if (length_ > allocatedLength_) {
    if (allocIncrement_ < 1) {
      throw std::runtime_error("feiArray::insert, can't resize, allocIncrement_==0");
    }

    while(allocatedLength_ < length_) {
      allocatedLength_ += allocIncrement_;
    }

    newData = new T[allocatedLength_];
    alloc_new_memory = true;
    dest = newData;
  }

  int len1 = offset;
  int len2 = length_-1 - offset;

  int sizeof_T = sizeof(T);
  snl_fei::Type_enum t_type = snl_fei::TypeTraits<T>::get_type();

  if (alloc_new_memory && len1 > 0) {
    if (t_type != snl_fei::OTHER) {
      std::memcpy(dest, src, len1*sizeof_T);
    }
    else {
      for(int i=0; i<len1; ++i) {
	dest[i] = src[i];
      }
    }
  }

  if (len2 > 0) {
    if (t_type != snl_fei::OTHER) {
      if (alloc_new_memory) {
	std::memcpy(dest+offset+1, src+offset, len2*sizeof_T);
      }
      else {
	memmove(dest+offset+1, src+offset, len2*sizeof_T);
      }
    }
    else {
      dest += offset+1;
      src += offset;
      for(int i=len2-1; i>=0; --i) {
	dest[i] = src[i];
      }
    }
  }

  if (alloc_new_memory) {
    if (length_>1) {
      delete [] data_;
    }
    data_ = newData;
  }

  data_[offset] = item;
}

//------------------------------------------------------------------------------
template<typename T, class COMPARE>
void feiArray<T,COMPARE>::append(const T& item)
{
  resize(length_+1);

  data_[length_-1] = item;
}

//------------------------------------------------------------------------------
template<typename T, class COMPARE>
feiArray<T,COMPARE>::feiArray(int len, int incr)
 : cmp_(),
   data_(NULL),
   length_(0),
   allocatedLength_(0),
   allocIncrement_(incr),
   userMemory_(false)
{
   resize(len);
}

//------------------------------------------------------------------------------
template<typename T, class COMPARE>
feiArray<T,COMPARE>::feiArray(int len, int allocLength, T* data)
  : cmp_(),
    data_(data),
    length_(len),
    allocatedLength_(allocLength),
    allocIncrement_(0),
    userMemory_(true)
{
}

//------------------------------------------------------------------------------
template<typename T, class COMPARE>
feiArray<T,COMPARE>::feiArray(const feiArray<T,COMPARE>& src)
 : cmp_(),
   data_(NULL),
   length_(0),
   allocatedLength_(0),
   allocIncrement_(src.allocIncrement_),
   userMemory_(false)
{
  reAllocate(allocatedLength_);
  resize(src.length_);

  for(int i=0; i<src.length_; ++i) {
    data_[i] = src.data_[i];
  }
}

//------------------------------------------------------------------------------
template<typename T, class COMPARE>
feiArray<T,COMPARE>::~feiArray() {
   deleteMemory();
}

//------------------------------------------------------------------------------
template<typename T, class COMPARE>
feiArray<T,COMPARE>& feiArray<T,COMPARE>::operator=(const feiArray<T,COMPARE>& src)
{
  resize(src.length_);

  for(int i=0; i<src.length_; i++) data_[i] = src.data_[i];

  return(*this);
}

//------------------------------------------------------------------------------
template<typename T, class COMPARE>
bool feiArray<T,COMPARE>::operator==(const feiArray<T,COMPARE>& rhs)
{
  if (length_ != rhs.length_) return(false);
  for(int i=0; i<length_; i++) {
    if (cmp_(data_[i], rhs.data_[i])) return(false);
    if (cmp_(rhs.data_[i], data_[i])) return(false);
  }
  return(true);
}

//------------------------------------------------------------------------------
template<typename T, class COMPARE>
int feiArray<T,COMPARE>::reAllocate(int newAllocLength) {
//
//This function alters the amount of allocated memory. It only changes the
//value reported by length() if newAllocLength is less than length().
//
  if (userMemory_) return(-1);

   if (newAllocLength < 0) newAllocLength = 0;

   if (newAllocLength == 0) {
      deleteMemory();
      return(0);
   }

   if (newAllocLength == allocatedLength_) return(0);

   T* newData = new T[newAllocLength];

   if (length_ > newAllocLength) length_ = newAllocLength;

   for(int i=0; i<length_; i++) newData[i] = data_[i];

   allocatedLength_ = newAllocLength;

   delete [] data_;
   data_ = newData;

   return(0);
}

//------------------------------------------------------------------------------
template<typename T, class COMPARE>
T& feiArray<T,COMPARE>::operator[](int offset)
{
  //
  //Bounds-checked array access.
  //
  if (offset < 0 || offset >= length_) {
    FEI_CERR << "feiArray: ERROR, out-of-bounds array access." << FEI_ENDL;
    FEI_CERR << "feiArray: Requested offset " << offset
         << " is outside array bounds [" << 0 << "..." << length_-1 << "]."
         << FEI_ENDL;
  }

  return(data_[offset]);
}

//------------------------------------------------------------------------------
template<typename T, class COMPARE>
const T& feiArray<T,COMPARE>::operator[](int offset) const
{
  //
  //Bounds-checked array access.
  //
  if (offset < 0 || offset >= length_) {
    FEI_CERR << "feiArray: ERROR, out-of-bounds array access." << FEI_ENDL;
    FEI_CERR << "feiArray: Requested offset " << offset
         << " is outside array bounds [" << 0 << "..." << length_-1 << "]."
         << FEI_ENDL;
  }

  return(data_[offset]);
}

//------------------------------------------------------------------------------
template<typename T, class COMPARE>
int feiArray<T,COMPARE>::find(const T& item) const
{
  for(int i=0; i<length_; i++) {
    if (cmp_(data_[i], item)) continue;
    if (cmp_(item, data_[i])) continue;
    return(i);
  }

  return(-1);
}

//------------------------------------------------------------------------------
template<typename T, class COMPARE>
void feiArray<T,COMPARE>::deleteMemory()
{
  if (allocatedLength_ > 0 && !userMemory_) {
    delete [] data_;
  }

  data_ = NULL;
  allocatedLength_ = 0;
  length_ = 0;
}

#endif
