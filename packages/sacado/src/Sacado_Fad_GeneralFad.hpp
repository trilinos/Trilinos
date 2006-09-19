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

#ifndef SACADO_FAD_IMPLEMENTATION_HPP
#define SACADO_FAD_IMPLEMENTATION_HPP

#include "Sacado_ConfigDefs.h"
#include "Sacado_Fad_Expression.hpp"

// forward decalarations
namespace Sacado {
  namespace Fad {
    template <class ExprT> class UnaryPlusOp;
    template <class ExprT> class UnaryMinusOp;
  }
}

namespace Sacado {

  //! Namespace for forward-mode AD classes
  namespace Fad {

    //! Forward-mode AD class templated on the storage for the deriavtive array
    /*!
     * This class provides the implementation of the Fad object required
     * for expression templating.  Class GeneralFad provides the complete
     * user inteface.
     */
    template <typename T, typename Storage> 
    class Implementation : protected Storage {

    public:

      //! Typename of values
      typedef T value_type;

      //! Default constructor
      Implementation() : Storage(), val_( T(0)) {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      Implementation(const T & x) : Storage(), val_(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      Implementation(const int sz, const T & x) : Storage(sz), val_(x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      Implementation(const int sz, const int i, const T & x) : 
	Storage(sz), val_(x) { 
	this->dx_[i]=1.; 
      }

      //! Copy constructor
      Implementation(const Implementation& x) : 
	Storage(x), val_(x.val_) {}

      //! Destructor
      ~Implementation() {}

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
      int size() const { return Storage::size();}

      //! Returns true if derivative array is not empty
      bool hasFastAccess() const { return Storage::size()!=0;}

      //! Returns true if derivative array is empty
      bool isPassive() const { return Storage::size()!=0;}
      
      //! Set whether variable is constant
      void setIsConstant(bool is_const) { 
	if (is_const && Storage::size()!=0)
	  Storage::resize(0);
      }

      //! Returns derivative array
      const T* dx() const { return &(this->dx_[0]);}

      //! Returns derivative component \c i with bounds checking
      T dx(int i) const { 
	return Storage::size() ? this->dx_[i] : T(0); }
    
      //! Returns derivative component \c i without bounds checking
      T& fastAccessDx(int i) { return this->dx_[i];}

      //! Returns derivative component \c i without bounds checking
      T fastAccessDx(int i) const { return this->dx_[i];}
    
      //@}

    protected:

      //! Value
      T val_;

    }; // class Implementation

    //! GeneralFad expression template specialization
    /*!
     * This template class represents a simple GeneralFad expression.
     */
    template <typename T, typename Storage> 
    class Expr< Implementation<T,Storage> > : 
      public Implementation<T,Storage> {

    public:

      //! Default constructor
      Expr() : Implementation<T,Storage>() {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      Expr(const T & x) : Implementation<T,Storage>(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      Expr(const int sz, const T & x) : Implementation<T,Storage>(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      Expr(const int sz, const int i, const T & x) : 
	Implementation<T,Storage>(sz,i,x) {}

      //! Copy constructor
      Expr(const Expr& x) : Implementation<T,Storage>(x) {}

    }; // class Expr< Implementation<T,Storage> >

    //! Forward-mode AD class templated on the storage for the derivative array
    /*!
     * This class provides the user interface of the GeneralFad object.  Class
     * Implementation provides the implementation.
     *
     * We should probably do more checking on derivative array dimensions.  
     */
    template <typename T, typename Storage>
    class GeneralFad : public Expr< Implementation<T,Storage> > {

    public:

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      GeneralFad() : Expr< Implementation<T,Storage> >() {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      GeneralFad(const T & x) : Expr< Implementation<T,Storage> >(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      GeneralFad(const int sz, const T & x) : 
	Expr< Implementation<T,Storage> >(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      GeneralFad(const int sz, const int i, const T & x) : 
	Expr< Implementation<T,Storage> >(sz,i,x) {}

      //! Copy constructor
      GeneralFad(const GeneralFad& x) : Expr< Implementation<T,Storage> >(x) {}

      //! Copy constructor from any Expression object
      template <typename S> GeneralFad(const Expr<S>& x);

      //! Set %GeneralFad object as the \c ith independent variable
      /*!
       * Sets the derivative array of length \c n to the \c ith row of the
       * identity matrix and has the same affect as the 
       * Implementation(const int sz, const int i, const T & x) 
       * constructor.
       */
      void diff(const int ith, const int n);

      //@}

      //! Destructor
      ~GeneralFad() {}

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      GeneralFad<T,Storage>& operator=(const T& val);

      //! Assignment with GeneralFad right-hand-side
      GeneralFad<T,Storage>& operator=(const GeneralFad<T,Storage>& x);

      //! Assignment operator with any expression right-hand-side
      template <typename S> GeneralFad<T,Storage>& operator=(const Expr<S>& x); 

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Unary-plus operator
      inline Expr< UnaryExpr< GeneralFad<T,Storage>, UnaryPlusOp > >
      operator + () const {
	typedef UnaryExpr< GeneralFad<T,Storage>, UnaryPlusOp > expr_t;
	return Expr<expr_t>(expr_t(*this));
      }

      //! Unary-minus operator
      inline Expr< UnaryExpr< GeneralFad<T,Storage>, UnaryMinusOp > >
      operator - () const {
	typedef UnaryExpr< GeneralFad<T,Storage>, UnaryMinusOp > expr_t;
	return Expr<expr_t>(expr_t(*this));
      }

      //! Addition-assignment operator with constant right-hand-side
      GeneralFad<T,Storage>& operator += (const T& x);

      //! Subtraction-assignment operator with constant right-hand-side
      GeneralFad<T,Storage>& operator -= (const T& x);

      //! Multiplication-assignment operator with constant right-hand-side
      GeneralFad<T,Storage>& operator *= (const T& x);

      //! Division-assignment operator with constant right-hand-side
      GeneralFad<T,Storage>& operator /= (const T& x);

      //! Addition-assignment operator with GeneralFad right-hand-side
      template <typename S> GeneralFad<T,Storage>& operator += (const S& x);

      //! Subtraction-assignment operator with GeneralFad right-hand-side
      template <typename S> GeneralFad<T,Storage>& operator -= (const S& x);
  
      //! Multiplication-assignment operator with GeneralFad right-hand-side
      template <typename S> GeneralFad<T,Storage>& operator *= (const S& x);

      //! Division-assignment operator with GeneralFad right-hand-side
      template <typename S> GeneralFad<T,Storage>& operator /= (const S& x);

      //@}

    }; // class GeneralFad<T,Storage>

  } // namespace Fad

} // namespace Sacado

#include "Sacado_Fad_GeneralFadTraits.hpp"
#include "Sacado_Fad_GeneralFadImp.hpp"
#include "Sacado_Fad_Ops.hpp"

#endif // SACADO_FAD_DFAD_HPP
