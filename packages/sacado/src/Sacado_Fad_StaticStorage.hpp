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

#ifndef SACADO_FAD_STATICSTORAGE_HPP
#define SACADO_FAD_STATICSTORAGE_HPP

namespace Sacado {

  namespace Fad {

    /*!
     * \brief Static array allocation class that works for any type
     */
    template <typename T, bool isScalar = IsScalarType<T>::value>
    struct ss_array {

      //! Copy array from \c src to \c dest of length \c sz
      static inline void copy(const T* src, T*  dest, int sz) {
	for (int i=0; i<sz; ++i)
	  *(dest++) = *(src++);
      }

      //! Zero out array \c dest of length \c sz
      static inline void zero(T* dest, int sz) {
	for (int i=0; i<sz; ++i)
	  *(dest++) = T(0.);
      }
    };

    /*!
     * \brief Static array allocation class that is specialized for scalar
     * i.e., fundamental or built-in types (float, double, etc...).
     */
//     template <typename T>
//     struct ss_array<T,true> {

//       //! Copy array from \c src to \c dest of length \c sz
//       static inline void copy(const T* src, T* dest, int sz) {
// 	memcpy(dest,src,sz*sizeof(T));
//       }

//       //! Zero out array \c dest of length \c sz
//       static inline void zero(T* dest, int sz) {
// 	memset(dest,0,sz*sizeof(T));
//       }
//     };

    //! Derivative array storage class using static memory allocation
    /*!
     * This class uses a statically allocated array whose dimension is fixed
     * by the template parameter \c Num.
     */
    template <typename T, int Num> 
    class StaticStorage {

    public:

      //! Default constructor
      StaticStorage(const T & x) : val_(x), sz_(0) {}

      //! Constructor with size \c sz
      /*!
       * Initializes derivative array 0 of length \c sz
       */
      StaticStorage(const int sz, const T & x) : val_(x), sz_(sz) { 
#ifdef SACADO_DEBUG
	if (sz > Num)
	  throw "StaticStorage::StaticStorage() Error:  Supplied derivative dimension exceeds maximum length.";
#endif
	ss_array<T>::zero(dx_, sz_); 
      }

      //! Copy constructor
      StaticStorage(const StaticStorage& x) : 
	val_(x.val_), sz_(x.sz_) { ss_array<T>::copy(x.dx_, dx_, sz_); }

      //! Destructor
      ~StaticStorage() {}

      //! Assignment
      StaticStorage& operator=(const StaticStorage& x) {
	val_ = x.val_;
	sz_ = x.sz_;
	ss_array<T>::copy(x.dx_, dx_, sz_);
	return *this;
      }

      //! Returns number of derivative components
      int size() const { return sz_;}

      //! Resize the derivative array to sz
      void resize(int sz) { 
#ifdef SACADO_DEBUG
	if (sz > Num)
	  throw "StaticStorage::resize() Error:  Supplied derivative dimension exceeds maximum length.";
#endif
	sz_ = sz; 
      }

      //! Zero out derivative array
      void zero() { ss_array<T>::zero(dx_, sz_); }

    public:

      //! Value
      T val_;

      //! Derivative array
      T dx_[Num];

      //! Size of derivative array
      int sz_;

    }; // class StaticStorage

  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_STATICSTORAGE_HPP
