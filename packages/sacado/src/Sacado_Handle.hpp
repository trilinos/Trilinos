// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
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
