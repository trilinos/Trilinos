#ifndef TEUCHOS_ARRAY_H
#define TEUCHOS_ARRAY_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Error.hpp"
#include "Teuchos_Utils.hpp"

namespace Teuchos
{

  /**
   * Array is a templated container class similar to the STL vector, but with
   * a more lightweight implementation and API as well as optional bounds-checking. 
   */
  template<class T>
  class Array
  {
  public:
    /** Empty ctor */
    Array();

    /** Allocate an array with n elements */
    Array(int n);

    /** Allocate n elements, and fill with value */
    Array(int n, const T& t);

    /** dtor deletes the contents */
    ~Array();

    /** copy ctor makes a deep copy of all elements */
    Array(const Array<T>& other);

    /** assignment operator makes a deep copy of all elements */
    const Array<T>& operator=(const Array<T>& other);

    /**
     * return number of elements.
     */
    int length() const ;

    /** change size */
    void resize(int newN);
    /** preallocate space */
    void reserve(int n);
    /** get preallocated space */
    int  capacity() const ;

    /** Stick a new entry at the end of the array. Resize to allow
     * space for the new entry
     */
    Array<T>& append(const T& entry);

    /**
     * Read/Write access to a the i-th element.
     */
    T& operator[](int i);
    /**
     * Read-only access to a the i-th element.
     */
    const T& operator[](int i) const;

    /**
     * Remove the i-th element. Subsequent elements are bumped back
     * to fill in the hole.
     */
    void remove(int i);

    /** Write as a string */
    std::string toString() const ;

    /** 
     * indicate whether Array has been compiled with boundschecking on 
     */
    static bool hasBoundsChecking();

  private:
    /** C array of contents */
    T* data_;
    /** number of active elements */
    int len_;
    /** number of elements that have been allocated */
    int capacity_;

    /** check for a bounds violation if HAVE_ARRAY_BOUNDSCHECK has been
     * defined as 1. */
    void indexCheckCrash(int i) const;
  };

  /** \relates Array */
  template<class T> std::ostream& operator<<(std::ostream& os, const Array<T>& array);

  /** \relates Array */
  template<class T> int hashCode(const Array<T>& array);

  /** \relates Array */
  template<class T> std::string toString(const Array<T>& array);

  template<class T> inline Array<T>::Array()
    : data_(0),
      len_(0),
      capacity_(0)
  {
  }



  template<class T> inline Array<T>::Array(int n)
    : data_(0),
      len_(n),
      capacity_(n)
  {
    if (len_ < 0)
      Error::raise("Negative length passed to Array<T>::Array(int n)");
    if (len_ > 0)
      {
        data_ = new T [len_];
        if (!data_)
          Error::raise("Array constructor out of memory");
      }
  }


  template<class T> inline Array<T>::Array(int n, const T& t)
    : data_(0),
      len_(n),
      capacity_(n)
  {
    if (len_ < 0)
      Error::raise("Negative length passed to Array<T>::Array(int n)");
    if (len_ > 0)
      {
        data_ = new T [len_];
        if (!data_)
          Error::raise("Array constructor out of memory");
      }
    for (int i = 0; i < len_; ++ i)
      data_[i] = t;
  }


  template<class T> inline  Array<T>::Array(const Array<T>& arr)
    : data_(0),
      len_(arr.len_),
      capacity_(arr.len_)
  {
    if (len_ > 0)
      {
        data_ = new T [capacity_];
        if (!data_)
          Error::raise("Array constructor out of memory");
        for (int i = 0; i < len_; ++i)
          data_[i] = arr.data_[i];
      }
  }

  template<class T> inline Array<T>::~Array()
  {
    delete [] data_;
  }

  template<class T> inline
  const Array<T>& Array<T>::operator=(const Array<T>& arr)
  {
    if (this != &arr)
      { // don't bother to assign if they're already identical
        if (capacity_ < arr.len_)
          { //If the reserved space is too small to hold arr
            delete [] data_;
            data_ = 0;
            capacity_ = arr.len_;
            if (capacity_ > 0) {
              data_ = new T[capacity_];
              if (!data_)
                Error::raise("Array constructor out of memory");
            }
          }
        len_ = arr.len_;
        for (int i = 0; i < len_; ++i)
          data_[i] = arr[i];
      }
    return *this;
  }


  /* bump up array size by factors of two until the desired size is reached */
  inline int bump(int start, int finish)
  {
    if (start == 0)
      start = 1;
    while (start < finish)
      start *= 2;
    return start;
  }


  template<class T> inline
  void Array<T>::resize(int newN) {
    if (len_ != newN) { // do not do anything if new size is not different
      if (newN < 0)
        Error::raise("Negative length passed to Array<T>::resize(int newN)");
      if(newN > capacity_)
        reserve(bump(capacity_, newN));
      len_ = newN;
    }
  }




  template<class T> inline
  void Array<T>::reserve(int N){
    if(capacity_ != N){
      if(N < 0){
        Error::raise("Negative length passed to Array<T>::reserve(int N)");
      }
      if(N < len_){ len_ = N;}
      capacity_ = N;
      T* oldData = data_;
      data_ = 0;
      data_ = new T [capacity_];
      if (!data_)
        Error::raise("Array<T>::reserve(int N) out of memory");
      for (int i = 0; i < len_; i++)
        data_[i] = oldData[i];
      delete [] oldData;
    }
  }

  template<class T> inline
  int Array<T>::capacity() const{
    return capacity_;
  }



  template<class T> inline
  Array<T>& Array<T>::append(const T& rhs)
  {
    resize(len_+1);
    data_[len_-1] = rhs;
    return *this;
  }



  template<class T> inline
  T& Array<T>::operator[](int i) {
#ifdef TEUCHOS_HAVE_ARRAY_BOUNDSCHECK
    indexCheckCrash(i);
#endif
    return data_[i];
  }

  template<class T> inline
  const T& Array<T>::operator[](int i) const {
#ifdef TEUCHOS_HAVE_ARRAY_BOUNDSCHECK
    indexCheckCrash(i);
#endif
    return data_[i];
  }

  template<class T> inline
  void Array<T>::remove(int i)
  {
#ifdef TEUCHOS_HAVE_ARRAY_BOUNDSCHECK
    indexCheckCrash(i);
#endif
    for (int j=i+1; j<length(); j++)
      {
        data_[j-1] = data_[j];
      }
    data_[len_-1] = T();
    len_--;
  }

  template<class T> inline
  int Array<T>::length() const {
    return len_;
  }

  template<class T> inline
  bool Array<T>::hasBoundsChecking()
  {
#ifdef TEUCHOS_HAVE_ARRAY_BOUNDSCHECK  
    return true;
#else
    return false;
#endif
  }

  template<class T> inline
  void Array<T>::indexCheckCrash(int i) const
  {
    if (i<0 || i>=len_)
      Error::boundsError("Array accessor", i, 0, len_);
  }

  template<class T>
  Array<T> sliceArray(const Array<Array<T> >& array, int index)
  {
    Array<T> rtn(array.length());

    for (int i=0; i<array.length(); i++)
      {
        rtn[i] = array[i][index];
      }
    return rtn;
  }


  // print in form (), (1), or (1,2)
  template<class T> inline ostream& operator<<(ostream& os, const Array<T>& array)
  {
    return os << toString(array);
  }

  template<class T> inline int hashCode(const Array<T>& array)
  {
    int rtn = hashCode(array.length());
    for (int i=0; i<array.length(); i++)
      {
        rtn += hashCode(array[i]);
      }
    return rtn;
  }

  template<class T> inline std::string Array<T>::toString() const
  {
    std::string rtn = "{";

    for (int i=0; i<length(); i++)
      {
        rtn += Teuchos::toString(data_[i]);
        if (i<length()-1) rtn += ", ";
      }
    rtn += "}";

    return rtn;
  }

  template<class T> inline std::string toString(const Array<T>& array)
  {
    return array.toString();
  }


  /** \relates Array create an array with one entry */
  template<class T> inline
  Array<T> tuple(const T& a)
  {
    Array<T> rtn(1, a);
    return rtn;
  }

  /** \relates Array create an array with two entries */
  template<class T> inline
  Array<T> tuple(const T& a, const T& b)
  {
    Array<T> rtn(2);
    rtn[0] = a;
    rtn[1] = b;
    return rtn;
  }

  /** \relates Array create an array with three entries */
  template<class T> inline
  Array<T> tuple(const T& a, const T& b, const T& c)
  {
    Array<T> rtn(3);
    rtn[0] = a;
    rtn[1] = b;
    rtn[2] = c;
    return rtn;
  }

  /** \relates Array create an array with four entries */
  template<class T> inline
  Array<T> tuple(const T& a, const T& b, const T& c, const T& d)
  {
    Array<T> rtn(4);
    rtn[0] = a;
    rtn[1] = b;
    rtn[2] = c;
    rtn[3] = d;
    return rtn;
  }

  /** \relates Array create an array with five entries */
  template<class T> inline
  Array<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e)
  {
    Array<T> rtn(5);
    rtn[0] = a;
    rtn[1] = b;
    rtn[2] = c;
    rtn[3] = d;
    rtn[4] = e;
    return rtn;
  }


  /** \relates Array create an array with six entries */
  template<class T> inline
  Array<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
                 const T& f)
  {
    Array<T> rtn(6);
    rtn[0] = a;
    rtn[1] = b;
    rtn[2] = c;
    rtn[3] = d;
    rtn[4] = e;
    rtn[5] = f;
    return rtn;
  }

  /** \relates Array create an array with seven entries */
  template<class T> inline
  Array<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
                 const T& f, const T& g)
  {
    Array<T> rtn(7);
    rtn[0] = a;
    rtn[1] = b;
    rtn[2] = c;
    rtn[3] = d;
    rtn[4] = e;
    rtn[5] = f;
    rtn[6] = g;
    return rtn;
  }

  /** \relates Array create an array with eight entries */
  template<class T> inline
  Array<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
                 const T& f, const T& g, const T& h)
  {
    Array<T> rtn(8);
    rtn[0] = a;
    rtn[1] = b;
    rtn[2] = c;
    rtn[3] = d;
    rtn[4] = e;
    rtn[5] = f;
    rtn[6] = g;
    rtn[7] = h;
    return rtn;
  }

  /** \relates Array create an array with nine entries */
  template<class T> inline
  Array<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
                 const T& f, const T& g, const T& h, const T& i)
  {
    Array<T> rtn(9);
    rtn[0] = a;
    rtn[1] = b;
    rtn[2] = c;
    rtn[3] = d;
    rtn[4] = e;
    rtn[5] = f;
    rtn[6] = g;
    rtn[7] = h;
    rtn[8] = i;
    return rtn;
  }


  /** \relates Array create an array with ten entries */
  template<class T> inline
  Array<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
                 const T& f, const T& g, const T& h, const T& i, const T& j)
  {
    Array<T> rtn(10);
    rtn[0] = a;
    rtn[1] = b;
    rtn[2] = c;
    rtn[3] = d;
    rtn[4] = e;
    rtn[5] = f;
    rtn[6] = g;
    rtn[7] = h;
    rtn[8] = i;
    rtn[9] = j;
    return rtn;
  }
}

#endif

