#ifndef _Epetra_Array_h_
#define _Epetra_Array_h_

#include <stdlib.h>

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

  NOTE ALSO: You can declare and use an 'Epetra_Array of Epetra_Arrays', e.g.,
  Epetra_Array<Epetra_Array<int> > int_table;
  HOWEVER, you can't use the find, binarySearch or insertSorted functions on
  the first dimension of the table.
*/

template<class T>
class Epetra_Array {
 public:
 /** Constructor
     @param length Defaults to 0 if not supplied
     @param incr Allocation increment, the 'chunk-size' by which memory is
     increased when new items are appended. Defaults to 16 if not supplied.
  */
   Epetra_Array(int length=0, int incr=16);

   /** Constructor to wrap an Epetra_Array around existing memory. When this 
       constructor is used, the reAllocate function will always return an error,
       and append/insert etc., will only work when 'length' is less than
       'allocLength'.
       @param length Amount of meaningful data already in place.
       @param allocLength Amount of memory allocated in 'data'.
       @param data Pointer to user-allocated array.
    */
   Epetra_Array(int length, int allocLength, T* data);

   /** Copy constructor. */
   Epetra_Array(const Epetra_Array<T>& src);

  /** Destructor */
   ~Epetra_Array();

   /** Assignment operator. */
   inline Epetra_Array<T>& operator=(const Epetra_Array<T>& src);

   /** Assignment operator for initialization. Sets the value 'item' 
    throughout the array.
   */
   inline Epetra_Array<T>& operator=(const T& item);

   /** Comparison operator. */
   bool operator==(const Epetra_Array<T>& rhs);

   /** Not-equal operator. */
   bool operator!=(const Epetra_Array<T>& rhs) {return(!( (*this) == rhs)); }

   /** Less-than operator. */
   bool operator<(const Epetra_Array<T>& rhs) {
     cerr << "Epetra_Array::operator< not implemented." << endl;
     abort();
     return(false);
   }

   /** Greater-than operator. */
   bool operator>(const Epetra_Array<T>& rhs) {
     cerr << "Epetra_Array::operator> not implemented." << endl;
     abort();
     return(false);
   }

   /** Return the length of the internal data array. This length
       is the amount of memory that has been requested by either the 'length'
       argument to the constructor, or by calls to the resize method. This
       length may be less than the actual allocated length of the array.
       This length will generally indicate how much of the allocated array
       is actually in use by the caller.
   */
   int length() const {return(length_);};

   /** Return the actual allocated length of the internal data array, which
     will always be equal to or greater than the value returned by length().
   */
   int allocatedLength() const {return(allocatedLength_);};

   /** Return the size in bytes, of the allocated memory in this Epetra_Array
       instance.
    */
   int allocatedSize() const {return(allocatedLength_*sizeof(int));};

   /** Alter the amount of memory that's available to the caller. 
    @param newLength If less than the currently available length, memory is
    not re-allocated, but the value returned by length() is adjusted
    downwards. If greater than the current length, then new memory is allocated
    as necessary, and any existing data is copied into the new memory. 
    @return err Return value is an error code which is 0 on successful 
    completion, -2 if allocation failed.
   */
   inline int resize(int newLength);

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

   /** Insert 'item' at position 'offset', regardless of whether this insertion
     causes the data to no longer be sorted. Expands allocated memory if
     necessary.
     @param item Entity of type T to be inserted.
     @param offset Position at which to insert 'item'. Should be in the
     range [0..length()] (i.e., offset may be 'too big' by one, which is 
     equivalent to calling append(item)).SThe
     @return err Return value is an error code which is 0 if successful.
     Unsuccessful execution may be caused by allocation failures (returns -2)
     or if offset is out of bounds (returns -1).
   */
   inline int insert(const T& item, int offset);


   /** Append 'item' to the data array, expanding allocated memory if
     necessary. 
     @param item Entity of type T to be appended.
     @return err Return value is an error code, similar to that of the
     insert function (returns -2 if allocation failed).
   */
   int append(const T& item);

   /** Insert 'item' at a position such that the data remains sorted (and only
     if 'item' is not already present).
     @param item Entity of type T to be inserted.
     @return offset If successfully inserted, returns the offset of the position
     at which item was inserted. (Returns -1 if 'item' was already present, -2
     if an allocation failed.)
   */
   inline int insertSorted(const T& item);

   /** Perform a linear (exhaustive) search.
    @param item Value to be searched for.
    @return position Offset of the position at which it was found. If it was 
    not found, returns -1.
   */
   int find(const T& item) const;

   /** Perform a binary search. This
     function assumes that the data is sorted. Its goal is to be faster than
     the above find function, and so does not verify that the data is
     sorted. If the user calls this function when the data is not sorted,
     then it may not perform correctly.
    @param item Value to be searched for.
    @return position Offset of the position at which item was found. If it was
    not found, returns -1.
   */
   inline int binarySearch(const T& item) const;

   /** Perform a binary search but limit the search to a given range.
    @param item Value to be searched for.
    @param start Starting offset of search 'window'.
    @param end Ending offset of search 'window'. end should be less than 
    length().
    @return offset position at which item was found. If not found, returns -1.
    Also returns -1 if start>end, or if start<0 or if end >= length(). 
    (Since 0-based indexing is used, 'end' can't be greater than length()-1.)
   */
   int binarySearch(const T& item, int start, int end) const;

   /** Same as binarySearch(item), but if 'item' is not found, gives the position
  at which it could be inserted with the list still being sorted.
  @param item
  @param insertPoint Not referenced if item is found.
  @return offset If item is found, this is its position. If item is not found,
  returns -1.
   */
   inline int binarySearch(const T& item, int& insertPoint) const;

 private:
   void deleteMemory();

   T* data_;
   int length_;
   int allocatedLength_;
   int allocIncrement_;
   bool userMemory_;
};

#ifdef EPETRA_ESI_INCLUDE_IMPLEMENTATION
#include "Epetra_Array.cpp"
#endif

//------------------------------------------------------------------------------
template<class T>
Epetra_Array<T>::~Epetra_Array() {
   deleteMemory();
}

//------------------------------------------------------------------------------
template<class T>
Epetra_Array<T>& Epetra_Array<T>::operator=(const Epetra_Array<T>& src)
{
  int err = resize(src.length_);
  if (err) cerr << "Epetra_Array::operator=: ERROR, resize failed." << endl;

  for(int i=0; i<src.length_; i++) data_[i] = src.data_[i];

  return(*this);
}

//------------------------------------------------------------------------------
template<class T>
Epetra_Array<T>& Epetra_Array<T>::operator=(const T& item)
{
  for(int i=0; i<length_; i++) data_[i] = item;

  return(*this);
}

//------------------------------------------------------------------------------
template<class T>
bool Epetra_Array<T>::operator==(const Epetra_Array<T>& rhs)
{
  if (length_ != rhs.length_) return(false);
  for(int i=0; i<length_; i++) if (data_[i] != rhs.data_[i]) return(false);
  return(true);
}

//------------------------------------------------------------------------------
template<class T>
int Epetra_Array<T>::resize(int newLength) {
//
//resize re-allocates the internal data array, if newLength is greater than
//allocatedLength_. If newLength is smaller than allocatedLength_, we just
//reset length_.
//
//If we are re-allocating the array to a larger size, any existing data is
//copied into the new array.
//
   if (newLength <= 0) {
      newLength = 0;
   }

   if (newLength <= allocatedLength_) {
      length_ = newLength;
      return(0);
   }

   if (newLength == 0 && allocatedLength_ == 0) return(0);

   int newAllocSize = allocatedLength_;

   if (newAllocSize < newLength && allocIncrement_ == 0) return(-1);

   while(newAllocSize < newLength) newAllocSize += allocIncrement_;

   T* newData = NULL;
   newData = new T[newAllocSize];

   if (!newData) {
      return(-2);
   }

   for(int i=0; i<length_; i++) {
      newData[i] = data_[i];
   }

   length_ = newLength;
   allocatedLength_ = newAllocSize;

   delete [] data_;
   data_ = newData;

   return(0);
}

//------------------------------------------------------------------------------
template<class T>
int Epetra_Array<T>::reAllocate(int newAllocLength) {
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

   if (newData == NULL) return(-2);

   if (length_ > newAllocLength) length_ = newAllocLength;

   for(int i=0; i<length_; i++) newData[i] = data_[i];

   allocatedLength_ = newAllocLength;

   delete [] data_;
   data_ = newData;

   return(0);
}

//------------------------------------------------------------------------------
template<class T>
T& Epetra_Array<T>::operator[](int offset) {
//
//Bounds-checked array access.
//
   if (offset < length_ && offset >= 0) return(data_[offset]);
   else {
      cerr << "Epetra_Array: ERROR, out-of-bounds array access." << endl;
      cerr << "Epetra_Array: Requested offset " << offset 
           << " is outside array bounds [" << 0 << "..." << length_-1 << "]."
           << endl;
      abort();
      return(data_[0]);
   }
}

//------------------------------------------------------------------------------
template<class T>
const T& Epetra_Array<T>::operator[](int offset) const {
//
//Bounds-checked array access.
//
   if (offset < length_ && offset >= 0) return(data_[offset]);
   else {
      cerr << "Epetra_Array: ERROR, out-of-bounds array access." << endl;
      cerr << "Epetra_Array: Requested offset " << offset
           << " is outside array bounds [" << 0 << "..." << length_-1 << "]."
           << endl;
      abort();
      return(data_[0]);
   }
}

//------------------------------------------------------------------------------
template<class T>
int Epetra_Array<T>::insert(const T& item, int offset) {
   if (offset < 0 || offset > length_) return(-1);

   if (offset == length_) return(append(item));

   int oldLength = length_;
   int err = resize(oldLength+1);
   if (err) return(err);

   int i;
   for(i=oldLength; i>offset; i--) {
      data_[i] = data_[i-1];
   }

   data_[offset] = item;

   return(0);
}

//------------------------------------------------------------------------------
template<class T>
int Epetra_Array<T>::append(const T& item) {
   int err = resize(length_+1);
   if (err) return(err);

   data_[length_-1] = item;

   return(0);
}

//------------------------------------------------------------------------------
template<class T>
int Epetra_Array<T>::insertSorted(const T& item) {
   int insertPoint = -1;
   int position = binarySearch(item, insertPoint);

   if (position >= 0) return(-1);
   else {
      int err = insert(item, insertPoint);
      if (err) return(err);
      else return(insertPoint);
   }
}

//------------------------------------------------------------------------------
template<class T>
int Epetra_Array<T>::find(const T& item) const {

   for(int i=0; i<length_; i++) {
      if (data_[i] == item) return(i);
   }

   return(-1);
}

//------------------------------------------------------------------------------
template<class T>
int Epetra_Array<T>::binarySearch(const T& item) const {

   int insertPoint = -1;

   return( binarySearch(item, insertPoint) );
}

//------------------------------------------------------------------------------
template<class T>
int Epetra_Array<T>::binarySearch(const T& item, int start, int end) const {
//
//Credit for this binary-search code goes to Ray Tuminaro, who wrote the
//Aztec function AZ_find_index.
//   
   if (length_ == 0) {
      return(-1);
   }

   if (start > end || start < 0 || end >= length_) return(-1);

   if (length_ == 1) {
      if (data_[0] == item) return(0);
      else {
         return(-1);
      }
   }

   unsigned ustart = start, uend = end;

   while(uend - ustart > 1) {
      int mid = (ustart + uend) >> 1;
      if (data_[mid] < item) ustart = mid;
      else uend = mid;
   }

   if (data_[ustart] == item) return((int)ustart);
   if (data_[uend] == item) return((int)uend);

   return(-1);
}

//------------------------------------------------------------------------------
template<class T>
int Epetra_Array<T>::binarySearch(const T& item, int& insertPoint) const {

   if (length_ == 0) {
      insertPoint = 0;
      return(-1);
   }

   if (length_ == 1) {
      if (data_[0] == item) return(0);
      else {
         if (data_[0] < item) insertPoint = 1;
         else insertPoint = 0;
         return(-1);
      }
   }

   unsigned start = 0;
   unsigned end = length_ - 1;

   while(end - start > 1) {
     unsigned mid = (start + end) >> 1;
     if (data_[mid] < item) start = mid;
     else end = mid;
   }

   if (data_[start] == item) return((int)start);
   if (data_[end] == item) return((int)end);

   if (data_[end] < item) insertPoint = (int)end+1;
   else if (data_[start] < item) insertPoint = (int)end;
   else insertPoint = (int)start;

   return(-1);
}

//------------------------------------------------------------------------------
template<class T>
void Epetra_Array<T>::deleteMemory() {

   if (allocatedLength_ > 0 && !userMemory_) {
      delete [] data_;
   }

   data_ = NULL;
   allocatedLength_ = 0;
   length_ = 0;
}

#endif
