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

#ifndef SACADO_FAD_DMFAD_HPP
#define SACADO_FAD_DMFAD_HPP

#include "Sacado_Fad_GeneralFadExpr.hpp"
#include "Sacado_Fad_MemPoolStorage.hpp"
#include "Sacado_Fad_DMFadTraits.hpp"

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
     * interface is provided by Sacado::Fad::GeneralFad.
     *
     * The class is templated on two types, \c ValueT and \c ScalarT.  Type
     * \c ValueT is the type for values the derivative class holds, while
     * type \c ScalarT is the type of basic scalars in the code being
     * differentiated (usually \c doubles).  When computing first derivatives, 
     * these two types are generally the same,  However when computing
     * higher derivatives, \c ValueT may be DFad<double> while \c ScalarT will
     * still be \c double.  Usually \c ScalarT does not need to be explicitly
     * specified since it can be deduced from \c ValueT through the template
     * metafunction ScalarValueType.
     */
    template <typename ValueT, 
	      typename ScalarT = typename ScalarValueType<ValueT>::type >
    class DMFad : public Expr< GeneralFad<ValueT,MemPoolStorage<ValueT> > > {

    public:

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      DMFad() : 
	Expr< GeneralFad< ValueT,MemPoolStorage<ValueT> > >() {}

      //! Constructor with supplied value \c x of type ValueT
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      DMFad(const ValueT& x) : 
	Expr< GeneralFad< ValueT,MemPoolStorage<ValueT> > >(x) {}

      //! Constructor with supplied value \c x of type ScalarT
      /*!
       * Initializes value to \c ValueT(x) and derivative array is empty
       */
      DMFad(const ScalarT& x) : 
	Expr< GeneralFad< ValueT,MemPoolStorage<ValueT> > >(ValueT(x)) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      DMFad(const int sz, const ValueT& x) : 
	Expr< GeneralFad< ValueT,MemPoolStorage<ValueT> > >(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      DMFad(const int sz, const int i, const ValueT & x) : 
	Expr< GeneralFad< ValueT,MemPoolStorage<ValueT> > >(sz,i,x) {}

      //! Copy constructor
      DMFad(const DMFad& x) : 
	Expr< GeneralFad< ValueT,MemPoolStorage<ValueT> > >(x) {}

      //! Copy constructor from any Expression object
      template <typename S> DMFad(const Expr<S>& x) : 
	Expr< GeneralFad< ValueT,MemPoolStorage<ValueT> > >(x) {}

      //@}

      //! Destructor
      ~DMFad() {}

      //! Assignment operator with any expression right-hand-side
      template <typename S> DMFad& operator=(const Expr<S>& x) 
      {
	GeneralFad< ValueT,MemPoolStorage<ValueT> >::operator=(x);
	return *this;
      }

      //! Set the default memory pool for new objects
      static void setDefaultPool(MemPool* pool) {
	MemPoolStorage<ValueT>::defaultPool_ = pool;
      }
	
    }; // class DMFad<ScalarT,ValueT>

    /*! 
     * \brief Forward-mode AD class using dynamic memory allocation and
     * expression templates.
     */
    /*!
     * This is the specialization of DMFad<ValueT,ScalarT> for when
     * \c ValueT and \c ScalarT are the same type.  It removes an extra
     * constructor that would be duplicated in this case.
     */
    template <typename ValueT>
    class DMFad<ValueT,ValueT> : 
      public Expr< GeneralFad<ValueT,MemPoolStorage<ValueT> > > {

    public:

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      DMFad() : 
	Expr< GeneralFad< ValueT,MemPoolStorage<ValueT> > >() {}

      //! Constructor with supplied value \c x of type ValueT
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      DMFad(const ValueT& x) : 
	Expr< GeneralFad< ValueT,MemPoolStorage<ValueT> > >(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      DMFad(const int sz, const ValueT& x) : 
	Expr< GeneralFad< ValueT,MemPoolStorage<ValueT> > >(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      DMFad(const int sz, const int i, const ValueT& x) : 
	Expr< GeneralFad< ValueT,MemPoolStorage<ValueT> > >(sz,i,x) {}

      //! Copy constructor
      DMFad(const DMFad& x) : 
	Expr< GeneralFad< ValueT,MemPoolStorage<ValueT> > >(x) {}

      //! Copy constructor from any Expression object
      template <typename S> DMFad(const Expr<S>& x) : 
	Expr< GeneralFad< ValueT,MemPoolStorage<ValueT> > >(x) {}

      //@}

      //! Destructor
      ~DMFad() {}

      //! Assignment operator with any expression right-hand-side
      template <typename S> DMFad& operator=(const Expr<S>& x) 
      {
	GeneralFad< ValueT,MemPoolStorage<ValueT> >::operator=(x);
	return *this;
      }

      //! Set the default memory pool for new objects
      static void setDefaultPool(MemPool* pool) {
	MemPoolStorage<ValueT>::defaultPool_ = pool;
      }

    }; // class DMFad<ValueT>

  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_DMFAD_HPP
