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

#ifndef SACADO_FAD_GENERALFAD_HPP
#define SACADO_FAD_GENERALFAD_HPP

#include "Sacado_ConfigDefs.h"
#include "Sacado_Fad_Expression.hpp"

namespace Sacado {

  //! Namespace for forward-mode AD classes
  namespace Fad {

    //! Forward-mode AD class templated on the storage for the derivative array
    /*!
     * This class provides a general forward mode AD implementation for any
     * type of derivative array storage.  It does not incorporate expression
     * templates.
     */
    template <typename T, typename Storage> 
    class GeneralFad {

    public:

      //! Typename of values
      typedef T value_type;

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor
      GeneralFad() : val_( T(0)), s_() {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      GeneralFad(const T & x) : val_(x), s_() {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      GeneralFad(const int sz, const T & x) : val_(x), s_(sz) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      GeneralFad(const int sz, const int i, const T & x) : 
	val_(x), s_(sz) { 
	s_.dx_[i]=1.; 
      }

      //! Copy constructor
      GeneralFad(const GeneralFad& x) : 
	val_(x.val_), s_(x.s_) {}

      //! Copy constructor from any Expression object
      template <typename S> GeneralFad(const Expr<S>& x);

      //! Destructor
      ~GeneralFad() {}

      //! Set %GeneralFad object as the \c ith independent variable
      /*!
       * Sets the derivative array of length \c n to the \c ith row of the
       * identity matrix and has the same affect as the 
       * Implementation(const int sz, const int i, const T & x) 
       * constructor.
       */
      void diff(const int ith, const int n);

      //! Resize derivative array to length \c sz
      /*!
       * This method does not (re)initialize the derivative components, so any
       * previous values may be lost.  Also any pointers to derivative
       * components may be invalid.
       */
      void resize(int sz) { s_.resize(sz); }

      //! Zero out the derivative array
      void zero() { s_.zero(); }

      //@}

      /*!
       * @name Value accessor methods
       */
      //@{

      //! Returns value
      const T& val() const { return val_;}

      //! Returns value
      T& val() { return val_;}

      //@}

      /*!
       * @name Derivative accessor methods
       */
      //@{

      //! Returns number of derivative components
      int size() const { return s_.size();}

      //! Returns true if derivative array is not empty
      bool hasFastAccess() const { return s_.size()!=0;}

      //! Returns true if derivative array is empty
      bool isPassive() const { return s_.size()!=0;}
      
      //! Set whether variable is constant
      void setIsConstant(bool is_const) { 
	if (is_const && s_.size()!=0)
	  s_.resize(0);
      }

      //! Returns derivative array
      const T* dx() const { return &(s_.dx_[0]);}

      //! Returns derivative component \c i with bounds checking
      T dx(int i) const { 
	return s_.size() ? s_.dx_[i] : T(0); }
    
      //! Returns derivative component \c i without bounds checking
      T& fastAccessDx(int i) { return s_.dx_[i];}

      //! Returns derivative component \c i without bounds checking
      T fastAccessDx(int i) const { return s_.dx_[i];}
    
      //@}

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      GeneralFad& operator=(const T& val);

      //! Assignment with Expr right-hand-side
      GeneralFad& 
      operator=(const GeneralFad& x);

      //! Assignment operator with any expression right-hand-side
      template <typename S> 
      GeneralFad& operator=(const Expr<S>& x); 

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Addition-assignment operator with constant right-hand-side
      GeneralFad& operator += (const T& x);

      //! Subtraction-assignment operator with constant right-hand-side
      GeneralFad& operator -= (const T& x);

      //! Multiplication-assignment operator with constant right-hand-side
      GeneralFad& operator *= (const T& x);

      //! Division-assignment operator with constant right-hand-side
      GeneralFad& operator /= (const T& x);

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S> 
      GeneralFad& operator += (const Expr<S>& x);

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S> 
      GeneralFad& operator -= (const Expr<S>& x);
  
      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S> 
      GeneralFad& operator *= (const Expr<S>& x);

      //! Division-assignment operator with Expr right-hand-side
      template <typename S> 
      GeneralFad& operator /= (const Expr<S>& x);

      //@}

    protected:

      //! Value
      T val_;

      //! Derivatives
      Storage s_;

    }; // class GeneralFad

  } // namespace Fad

} // namespace Sacado

#include "Sacado_Fad_GeneralFadImp.hpp"

#endif // SACADO_FAD_GENERALFAD_HPP
