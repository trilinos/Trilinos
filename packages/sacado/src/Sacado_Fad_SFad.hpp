// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2004) Sandia Corporation
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

#ifndef SACADO_FAD_SFAD_HPP
#define SACADO_FAD_SFAD_HPP

#include "Sacado_Fad_GeneralFad.hpp"

namespace Sacado {

  namespace Fad {

    // Forward declaration
    template <typename T, int Num> 
    class StaticStorage;

    //! Forward-mode AD class using static memory allocation
    /*!
     * This is the user-level class for forward mode AD with static
     * memory allocation, and is appropriate for whenever the number
     * of derivative components is known at compile time.  The largest size
     * of the derivative array is fixed by the template parameter \c Num
     * while the actual size used is set by the \c sz argument to the 
     * constructor or the \c n argument to diff().  The user
     * interface is split between the Sacado::Fad::GeneralFad and 
     * and Sacado::Fad::Implementation classes.
     */
    template <typename T, int Num>
    class SFad : public GeneralFad<T,StaticStorage<T,Num> > {

    public:

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      SFad() : GeneralFad<T,StaticStorage<T,Num> >() {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      SFad(const T & x) : GeneralFad<T,StaticStorage<T,Num> >(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      SFad(const int sz, const T & x) : 
	GeneralFad<T,StaticStorage<T,Num> >(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      SFad(const int sz, const int i, const T & x) : 
	GeneralFad<T,StaticStorage<T,Num> >(sz,i,x) {}

      //! Copy constructor
      SFad(const SFad& x) : GeneralFad<T,StaticStorage<T,Num> >(x) {}

      //! Copy constructor from any Expression object
      template <typename S> SFad(const Expr<S>& x) : 
	GeneralFad<T,StaticStorage<T,Num> >(x) {}

      //@}

      //! Destructor
      ~SFad() {}

    }; // class SFad<T>

    //! Derivative array storage class using static memory allocation
    /*!
     * This class uses a statically allocated array whose dimension is fixed
     * by the template parameter \c Num.
     */
    template <typename T, int Num> 
    class StaticStorage {

    public:

      //! Default constructor
      StaticStorage() : sz_(0) {}

      //! Constructor with size \c sz
      /*!
       * Initializes derivative array 0 of length \c sz
       */
      StaticStorage(const int sz) : sz_(sz) { memset(dx_,0,sz_*sizeof(T)); }

      //! Copy constructor
      StaticStorage(const StaticStorage& x) : 
	sz_(x.sz_) { for (int i=0; i<sz_; i++) dx_[i] = x.dx_[i]; }

      //! Destructor
      ~StaticStorage() {}

      //! Assignment
      StaticStorage& operator=(const StaticStorage& x) {
	sz_ = x.sz_;
	for (int i=0; i<sz_; i++) 
	  dx_[i] = x.dx_[i]; 
	return *this;
      }

      //! Returns number of derivative components
      int size() const { return sz_;}

      //! Resize the derivative array to sz
      void resize(int sz) { sz_ = sz; }

      //! Zero out derivative array
      void zero() { memset(dx_,0,sz_*sizeof(T)); }

    protected:

      //! Size of derivative array
      int sz_;

      //! Derivative array
      T dx_[Num];

    }; // class DynamicStorage

  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_SFAD_HPP
