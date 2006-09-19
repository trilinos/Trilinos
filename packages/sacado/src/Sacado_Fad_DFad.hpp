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

#ifndef SACADO_FAD_DFAD_HPP
#define SACADO_FAD_DFAD_HPP

#include <valarray>
#include "Sacado_Fad_GeneralFad.hpp"

namespace Sacado {

  namespace Fad {

    // Forward declaration
    template <typename T> 
    class DynamicStorage;

    //! Forward-mode AD class using dynamic memory allocation
    /*!
     * This is the user-level class for forward mode AD with dynamic
     * memory allocation, and is appropriate for whenever the number
     * of derivative components is not known at compile time.  The user
     * interface is split between the Sacado::Fad::GeneralFad and 
     * and Sacado::Fad::Implementation classes.
     */
    template <typename T>
    class DFad : public GeneralFad<T,DynamicStorage<T> > {

    public:

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      DFad() :  GeneralFad<T,DynamicStorage<T> >() {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      DFad(const T & x) : GeneralFad<T,DynamicStorage<T> >(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      DFad(const int sz, const T & x) : 
	GeneralFad<T,DynamicStorage<T> >(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      DFad(const int sz, const int i, const T & x) : 
	GeneralFad<T,DynamicStorage<T> >(sz,i,x) {}

      //! Copy constructor
      DFad(const DFad& x) : GeneralFad<T,DynamicStorage<T> >(x) {}

      //! Copy constructor from any Expression object
      template <typename S> DFad(const Expr<S>& x) : 
	GeneralFad<T,DynamicStorage<T> >(x) {}

      //@}

      //! Destructor
      ~DFad() {}

    }; // class DFad<T>

    //! Derivative array storage class using dynamic memory allocation
    /*!
     * This class dynamically allocates the derivative array using
     * std::valarray.
     *
     * We might want to allow the user to provide a custom allocator.
     */
    template <typename T> 
    class DynamicStorage {

    public:

      //! Default constructor
      DynamicStorage() : dx_() {}

      //! Constructor with size \c sz
      /*!
       * Initializes derivative array 0 of length \c sz
       */
      DynamicStorage(const int sz) : dx_(T(0),sz) {}

      //! Copy constructor
      DynamicStorage(const DynamicStorage& x) : 
	dx_(x.dx_) {}

      //! Destructor
      ~DynamicStorage() {}

      //! Assignment
      DynamicStorage& operator=(const DynamicStorage& x) { 
	if (dx_.size() != x.dx_.size())
	  dx_.resize(x.dx_.size());
	dx_ = x.dx_; 
	return *this; 
      } 

      //! Returns number of derivative components
      int size() const { return dx_.size();}

      //! Resize the derivative array to sz
      void resize(int sz) { dx_.resize(sz); }

      //! Zero out derivative array
      void zero() { dx_ = T(0.); }

    protected:

      //! Derivative array
      std::valarray<T> dx_;

    }; // class DynamicStorage

  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_DFAD_HPP
