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

#ifndef SACADO_FAD_DVFAD_HPP
#define SACADO_FAD_DVFAD_HPP

#include "Sacado_Fad_GeneralFadExpr.hpp"
#include "Sacado_Fad_VectorDynamicStorage.hpp"
#include "Sacado_Fad_DVFadTraits.hpp"
#include "Sacado_dummy_arg.hpp"

namespace Sacado {

  namespace Fad {

    /*! 
     * \brief Forward-mode AD class using dynamic memory allocation and
     * expression templates.
     */
    /*!
     * This is the user-level class for forward mode AD with dynamic
     * memory allocation, and is appropriate for whenever the number
     * of derivative components is not known at compile time.  The user
     * interface is provided by Sacado::Fad::GeneralVFad.
     */
    template <typename ValueT>
    class DVFad : 
      public Expr< GeneralFad<ValueT,VectorDynamicStorage<ValueT> > > {

    public:

      //! Typename of scalar's (which may be different from ValueT)
      typedef typename ScalarType<ValueT>::type ScalarT;

      //! Turn DVFad into a meta-function class usable with mpl::apply
      template <typename T> 
      struct apply {
	typedef DVFad<T> type;
      };

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      DVFad() : 
	Expr< GeneralFad< ValueT,VectorDynamicStorage<ValueT> > >() {}

      //! Constructor with supplied value \c x of type ValueT
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      DVFad(const ValueT& x) : 
	Expr< GeneralFad< ValueT,VectorDynamicStorage<ValueT> > >(x) {}

      //! Constructor with supplied value \c x of type ScalarT
      /*!
       * Initializes value to \c ValueT(x) and derivative array is empty.
       * Creates a dummy overload when ValueT and ScalarT are the same type.
       */
      DVFad(const typename dummy<ValueT,ScalarT>::type& x) : 
	Expr< GeneralFad< ValueT,VectorDynamicStorage<ValueT> > >(ValueT(x)) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      DVFad(const int sz, const ValueT& x) : 
	Expr< GeneralFad< ValueT,VectorDynamicStorage<ValueT> > >(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      DVFad(const int sz, const int i, const ValueT & x) : 
	Expr< GeneralFad< ValueT,VectorDynamicStorage<ValueT> > >(sz,i,x) {}

      //! Constructor with supplied memory
      /*!
       * Initializes value to point to \c x and derivative array to point 
       * to\c dx.  Derivative array is zero'd out if \c zero_out is true.
       */
      DVFad(const int sz, ValueT* x, ValueT* dx, bool zero_out = false) : 
	Expr< GeneralFad< ValueT,VectorDynamicStorage<ValueT> > >(sz,x,dx,zero_out) {}

      //! Constructor with supplied memory and index \c i
      /*!
       * Initializes value to point to \c x and derivative array to point 
       * to\c dx.  Initializes derivative array row \c i of the identity matrix,
       * i.e., sets derivative component \c i to 1 and all other's to zero.
       */
      DVFad(const int sz, const int i, ValueT* x, ValueT* dx) : 
	Expr< GeneralFad< ValueT,VectorDynamicStorage<ValueT> > >(sz,i,x,dx) {}

      //! Copy constructor
      DVFad(const DVFad& x) : 
	Expr< GeneralFad< ValueT,VectorDynamicStorage<ValueT> > >(x) {}

      //! Copy constructor from any Expression object
      template <typename S> DVFad(const Expr<S>& x) : 
	Expr< GeneralFad< ValueT,VectorDynamicStorage<ValueT> > >(x) {}

      //@}

      //! Destructor
      ~DVFad() {}

      //! Assignment operator with constant right-hand-side
      DVFad& operator=(const ValueT& v) {
	GeneralFad< ValueT,VectorDynamicStorage<ValueT> >::operator=(v);
	return *this;
      }

      //! Assignment operator with constant right-hand-side
      /*!
       * Creates a dummy overload when ValueT and ScalarT are the same type.
       */
      DVFad& operator=(const typename dummy<ValueT,ScalarT>::type& v) {
	GeneralFad< ValueT,VectorDynamicStorage<ValueT> >::operator=(ValueT(v));
	return *this;
      }

      //! Assignment operator with DVFad right-hand-side
      DVFad& operator=(const DVFad& x) {
	GeneralFad< ValueT,VectorDynamicStorage<ValueT> >::operator=(static_cast<const GeneralFad< ValueT,VectorDynamicStorage<ValueT> >&>(x));
	return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S> DVFad& operator=(const Expr<S>& x) 
      {
	GeneralFad< ValueT,VectorDynamicStorage<ValueT> >::operator=(x);
	return *this;
      }
	
    }; // class DVFad<ValueT>

  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_DFAD_HPP
