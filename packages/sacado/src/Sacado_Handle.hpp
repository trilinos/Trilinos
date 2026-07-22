// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_HANDLE_HPP
#define SACADO_HANDLE_HPP

namespace Sacado {

  /*!
   * \brief A generic handle class
   */
  template <typename T>
  class Handle {
  public:

    //! Create new handle from pointer \c p
    Handle(T* p) : rep(p), count(new int(1)) {}

    //! Copy handle
    Handle(const Handle& h) : rep(h.rep), count(h.count) { (*count)++; }

    //! Destroy handle
    ~Handle() { decrementRef(); }

    //! Return pointer to underlying data
    T* get() { return rep; }

    //! Return pointer to underlying data
    const T* get() const { return rep; }

    //! Assign to handle \c h as its own copy
    void Assign(const Handle& h) {
      decrementRef();
      rep = new T(*(h.rep));
      count = new int(1);
    }

    //! Make handle have its own copy of \c rep
    void makeOwnCopy() {
      T *tmp;
      if (*count > 1) {
	tmp = rep;
	(*count)--;
	rep = new T(*tmp);
	count = new int(1);
      }
    }
    
    //! Assignment operator
    Handle& operator = (const Handle& h) {
      if (this != &h) {
	decrementRef();
	rep = h.rep; 
	count = h.count;
	(*count)++;
      }
      return *this;
    }

    //! Dereference
    T* operator -> () const { return rep; }

    //! Dereference
    const T& operator * () const { return *rep; }

    //! Dereference
    T& operator * () { return *rep; }
    
  private:

    //! Pointer to data
    T *rep;

    //! Reference count
    int *count;
    
    //! Decrement reference 
    void decrementRef() {
      (*count)--;
      if (*count == 0) {
	delete rep;
	delete count;
      }
    }

  }; // class Handle

  //! Compare two handles
  template <typename T>
  bool operator==(const Handle<T>& h1, const Handle<T>& h2) {
    return h1.get() == h2.get();
  }

} // namespace Sacado

#endif // SACADO_HANDLE_HPP
