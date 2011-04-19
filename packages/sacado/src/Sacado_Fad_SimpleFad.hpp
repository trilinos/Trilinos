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

#ifndef SACADO_FAD_SIMPLEFAD_HPP
#define SACADO_FAD_SIMPLEFAD_HPP

#include "Sacado_Fad_SimpleFadTraits.hpp"
#include "Sacado_Fad_DynamicStorage.hpp"
#include "Sacado_dummy_arg.hpp"

namespace Sacado {

  namespace Fad {

    /*! 
     * \brief Forward-mode AD class using dynamic memory allocation but no
     * expression templates.
     */
    /*!
     * This is the user-level class for forward mode AD with dynamic
     * memory allocation, and is appropriate for whenever the number
     * of derivative components is not known at compile time.  The user
     * interface is provided by Sacado::Fad::GeneralFad.
     */
    template <typename ValueT>
    class SimpleFad : public GeneralFad<ValueT,DynamicStorage<ValueT> >  {

    public:

      //! Typename of scalar's (which may be different from ValueT)
      typedef typename ScalarType<ValueT>::type ScalarT;

      //! Turn SimpleFad into a meta-function class usable with mpl::apply
      template <typename T> 
      struct apply {
	typedef SimpleFad<T> type;
      };

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      SimpleFad() : 
	GeneralFad< ValueT,DynamicStorage<ValueT> >() {}

      //! Constructor with supplied value \c x of type ValueT
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      SimpleFad(const ValueT& x) : 
	GeneralFad< ValueT,DynamicStorage<ValueT> >(x) {}

      //! Constructor with supplied value \c x of type ScalarT
      /*!
       * Initializes value to \c ValueT(x) and derivative array is empty.
       * Creates a dummy overload when ValueT and ScalarT are the same type.
       */
      SimpleFad(const typename dummy<ValueT,ScalarT>::type& x) : 
	GeneralFad< ValueT,DynamicStorage<ValueT> >(ValueT(x)) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      SimpleFad(const int sz, const ValueT& x) : 
	GeneralFad< ValueT,DynamicStorage<ValueT> >(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      SimpleFad(const int sz, const int i, const ValueT & x) : 
	GeneralFad< ValueT,DynamicStorage<ValueT> >(sz,i,x) {}

      //! Copy constructor
      SimpleFad(const SimpleFad& x) : 
	GeneralFad< ValueT,DynamicStorage<ValueT> >(x) {}

      //! Tangent copy constructor
      SimpleFad(const SimpleFad& x, const ValueT& v, const ValueT& partial) :
	GeneralFad< ValueT,DynamicStorage<ValueT> >(x.size(), v) {
	for (int i=0; i<this->size(); i++)
	  this->fastAccessDx(i) = x.fastAccessDx(i)*partial;
      }

      //@}

      //! Destructor
      ~SimpleFad() {}

      //! Assignment operator with constant right-hand-side
      SimpleFad& operator=(const ValueT& v) {
	GeneralFad< ValueT,DynamicStorage<ValueT> >::operator=(v);
	return *this;
      }

      //! Assignment operator with constant right-hand-side
      /*!
       * Creates a dummy overload when ValueT and ScalarT are the same type.
       */
      SimpleFad& operator=(const typename dummy<ValueT,ScalarT>::type& v) {
	GeneralFad< ValueT,DynamicStorage<ValueT> >::operator=(ValueT(v));
	return *this;
      }

      //! Assignment operator with SimpleFad right-hand-side
      SimpleFad& operator=(const SimpleFad& x) {
	GeneralFad< ValueT,DynamicStorage<ValueT> >::operator=(static_cast<const GeneralFad< ValueT,DynamicStorage<ValueT> >&>(x));
	return *this;
      }

      //! Addition-assignment operator with Expr right-hand-side
      SimpleFad& operator += (const SimpleFad& x);

      //! Subtraction-assignment operator with Expr right-hand-side
      SimpleFad& operator -= (const SimpleFad& x);
  
      //! Multiplication-assignment operator with Expr right-hand-side
      SimpleFad& operator *= (const SimpleFad& x);

      //! Division-assignment operator with Expr right-hand-side
      SimpleFad& operator /= (const SimpleFad& x);
	
    }; // class SimpleFad<ValueT>

  } // namespace Fad

} // namespace Sacado

// Include method definitions
#include "Sacado_Fad_SimpleFadImp.hpp"

// Include elementary operation overloads
#include "Sacado_Fad_SimpleFadOps.hpp"

#endif // SACADO_FAD_SIMPLEFAD_HPP
