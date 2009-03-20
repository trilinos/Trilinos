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
//
// The forward-mode AD classes in Sacado are a derivative work of the
// expression template classes in the Fad package by Nicolas Di Cesare.  
// The following banner is included in the original Fad source code:
//
// ************ DO NOT REMOVE THIS BANNER ****************
//
//  Nicolas Di Cesare <Nicolas.Dicesare@ann.jussieu.fr>
//  http://www.ann.jussieu.fr/~dicesare
//
//            CEMRACS 98 : C++ courses, 
//         templates : new C++ techniques 
//            for scientific computing 
// 
//********************************************************
//
//  A short implementation ( not all operators and 
//  functions are overloaded ) of 1st order Automatic
//  Differentiation in forward mode (FAD) using
//  EXPRESSION TEMPLATES.
//
//********************************************************
// @HEADER

#ifndef SACADO_FAD_GENERALVFAD_HPP
#define SACADO_FAD_GENERALVFAD_HPP

#include "Sacado_Fad_Expression.hpp"

namespace Sacado {

  //! Namespace for forward-mode AD classes
  namespace Fad {

    /*! 
     * \brief Forward-mode AD class templated on the storage for the 
     * derivative array that supports contiguous allocation of values and
     * derivative components (i.e., vector mode).
     */
    /*!
     * This class provides a general forward mode AD implementation for any
     * type of derivative array storage.  It does not incorporate expression
     * templates.  It is different from GeneralFad in that it adds a new
     * constructor providing pointers to the value and derivative components.
     * The value is then accessed in the storage through dereferencing a 
     * pointer (which adds extra cost), and the derivative array can have
     * a non-unit stride (which also adds additional cost).  
     */
    template <typename T, typename Storage> 
    class GeneralVFad {

    public:

      //! Typename of values
      typedef T value_type;

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor
      GeneralVFad() : s_(T(0.)) {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      GeneralVFad(const T & x) : s_(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      GeneralVFad(const int sz, const T & x) : s_(sz, x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      GeneralVFad(const int sz, const int i, const T & x) : 
	s_(sz, x) { 
	s_.dx_[i*s_.stride_]=1.; 
      }

      //! Constructor with supplied memory
      /*!
       * Initializes value to point to \c x and derivative array to point 
       * to\c dx.  Stride of the derivative array is specified through
       * \c stride.  Derivative array is zero'd out if \c zero_out is true.
       */
      GeneralVFad(const int sz, T* x, T* dx_, const int stride = 1,
		  bool zero_out = false) 
	: s_(sz, x, dx_, stride, zero_out) {}

      //! Constructor with supplied memory and index \c i
      /*!
       * Initializes value to point to \c x and derivative array to point 
       * to\c dx.  Stride of the derivative array is specified through
       * \c stride.  Initializes derivative array row \c i of the identity 
       * matrix, i.e., sets derivative component \c i to 1 and all other's to 
       * zero.
       */
      GeneralVFad(const int sz, const int i, T* x, T* dx_, const int stride = 1)
	: 
	s_(sz, x, dx_, stride, true) { 
	s_.dx_[i*s_.stride_]=1.; 
      }

      //! Copy constructor
      GeneralVFad(const GeneralVFad& x) : 
	s_(x.s_) {}

      //! Copy constructor from any Expression object
      template <typename S> GeneralVFad(const Expr<S>& x);

      //! Destructor
      ~GeneralVFad() {}

      //! Set %GeneralVFad object as the \c ith independent variable
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
       * components may be invalid.  However if the memory for the derivative
       * array is supplied through the constructor, the array cannot be
       * resized beyond its original length.
       */
      void resize(int sz) { s_.resize(sz); }

      //! Zero out the derivative array
      void zero() { s_.zero(); }

      //! Set value/derivative memory
      /*!
       * Any previous values in value/derivative array storage will be lost
       */
      void setMemory(const int sz, T* x, T* dx_, const int stride = 1) { 
	s_.setMemory(sz,x,dx_,stride); }

      //@}

      /*!
       * @name Value accessor methods
       */
      //@{

      //! Returns value
      const T& val() const { return *s_.val_;}

      //! Returns value
      T& val() { return *s_.val_;}

      //@}

      /*!
       * @name Derivative accessor methods
       */
      //@{

      //! Returns number of derivative components
      int size() const { return s_.size(); }

      /*! 
       * \brief Returns number of derivative components that can be stored 
       * without reallocation
       */
      int availableSize() const { return s_.length(); }

      //! Returns true if derivative array is not empty
      bool hasFastAccess() const { return s_.size()!=0; }

      //! Returns true if derivative array is empty
      bool isPassive() const { return s_.size()!=0; }
      
      //! Set whether variable is constant
      void setIsConstant(bool is_const) { 
	if (is_const && s_.size()!=0)
	  s_.resize(0);
      }

      //! Returns derivative array
      const T* dx() const { return &(s_.dx_[0]);}

      //! Returns derivative component \c i with bounds checking
      T dx(int i) const { 
	return s_.size() ? s_.dx_[i*s_.stride_] : T(0.); }
    
      //! Returns derivative component \c i without bounds checking
      T& fastAccessDx(int i) { return s_.dx_[i*s_.stride_];}

      //! Returns derivative component \c i without bounds checking
      const T& fastAccessDx(int i) const { return s_.dx_[i*s_.stride_];}
    
      //@}

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      GeneralVFad& operator=(const T& val);

      //! Assignment with Expr right-hand-side
      GeneralVFad& 
      operator=(const GeneralVFad& x);

      //! Assignment operator with any expression right-hand-side
      template <typename S> 
      GeneralVFad& operator=(const Expr<S>& x); 

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Addition-assignment operator with constant right-hand-side
      GeneralVFad& operator += (const T& x);

      //! Subtraction-assignment operator with constant right-hand-side
      GeneralVFad& operator -= (const T& x);

      //! Multiplication-assignment operator with constant right-hand-side
      GeneralVFad& operator *= (const T& x);

      //! Division-assignment operator with constant right-hand-side
      GeneralVFad& operator /= (const T& x);

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S> 
      GeneralVFad& operator += (const Expr<S>& x);

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S> 
      GeneralVFad& operator -= (const Expr<S>& x);
  
      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S> 
      GeneralVFad& operator *= (const Expr<S>& x);

      //! Division-assignment operator with Expr right-hand-side
      template <typename S> 
      GeneralVFad& operator /= (const Expr<S>& x);

      //@}

    protected:

      //! Value & Derivatives
      Storage s_;

    }; // class GeneralVFad

  } // namespace Fad

} // namespace Sacado

#include "Sacado_Fad_GeneralVFadImp.hpp"

#endif // SACADO_FAD_GENERALVFAD_HPP
