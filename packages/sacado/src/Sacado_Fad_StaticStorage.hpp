// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

#include "Sacado_ConfigDefs.h"
#include "Sacado_StaticArrayTraits.hpp"

namespace Sacado {

  namespace Fad {

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
	val_(x.val_), sz_(x.sz_) { 
	//ss_array<T>::copy(x.dx_, dx_, sz_); 
	for (int i=0; i<sz_; i++)
	  dx_[i] = x.dx_[i];
      }

      //! Destructor
      ~StaticStorage() {}

      //! Assignment
      StaticStorage& operator=(const StaticStorage& x) {
	val_ = x.val_;
	sz_ = x.sz_;
	//ss_array<T>::copy(x.dx_, dx_, sz_);
	for (int i=0; i<sz_; i++)
	  dx_[i] = x.dx_[i];
	return *this;
      }

      //! Returns number of derivative components
      int size() const { return sz_;}

      //! Returns array length
      int length() const { return Num; }

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
