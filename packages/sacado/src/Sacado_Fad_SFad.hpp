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

#ifndef SACADO_FAD_SFAD_HPP
#define SACADO_FAD_SFAD_HPP

#include "Sacado_Fad_SFadTraits.hpp"
#include "Sacado_Fad_Expression.hpp"
#include "Sacado_StaticArrayTraits.hpp"

namespace Sacado {

  //! Namespace for forward-mode AD classes
  namespace Fad {

    //! A tag for specializing Expr for SFad expressions
    template <typename T, int Num> 
    struct SFadExprTag {};

    /*!
     * \brief Expression template forward-mode AD class with static memory
     * allocation.
     */
    /*!
     * This classes specializes Expr to SFad expressions.
     */
    template <typename T, int Num> 
    class Expr< SFadExprTag<T,Num> > {

    public:

      //! Typename of values
      typedef T value_type;

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor
      Expr() : val_( T(0.)) { ss_array<T>::zero(dx_, Num); }

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      Expr(const T & x) : val_(x)  { ss_array<T>::zero(dx_, Num); }

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      Expr(const int sz, const T & x);

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      Expr(const int sz, const int i, const T & x);

      //! Copy constructor
      Expr(const Expr& x);

      //! Copy constructor from any Expression object
      template <typename S> Expr(const Expr<S>& x);

      //! Destructor
      ~Expr() {}

      //! Set %Fad object as the \c ith independent variable
      /*!
       * Sets the derivative array of length \c n to the \c ith row of the
       * identity matrix and has the same affect as the 
       * Implementation(const int sz, const int i, const T & x) 
       * constructor.
       */
      void diff(const int ith, const int n);

      //! Resize derivative array to length \c sz
      /*!
       * Since the derivative array length is not dynamic, this method
       * throws an error if compiled with SACADO_DEBUG defined.
       */
      void resize(int sz);

      //! Zero out the derivative array
      void zero() { ss_array<T>::zero(dx_, Num); }

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
      int size() const { return Num;}

      /*! 
       * \brief Returns number of derivative components that can be stored 
       * without reallocation
       */
      int availableSize() const { return Num; }

      //! Returns true if derivative array is not empty
      bool hasFastAccess() const { return true; }

      //! Returns true if derivative array is empty
      bool isPassive() const { return false; }
      
      //! Set whether variable is constant
      void setIsConstant(bool is_const) {}

      //! Returns derivative array
      const T* dx() const { return &(dx_[0]);}

      //! Returns derivative component \c i with bounds checking
      const T& dx(int i) const { return dx_[i]; }
    
      //! Returns derivative component \c i without bounds checking
      T& fastAccessDx(int i) { return dx_[i];}

      //! Returns derivative component \c i without bounds checking
      const T& fastAccessDx(int i) const { return dx_[i];}
    
      //@}

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      Expr< SFadExprTag<T,Num> >& operator=(const T& val);

      //! Assignment with Expr right-hand-side
      Expr< SFadExprTag<T,Num> >& 
      operator=(const Expr< SFadExprTag<T,Num> >& x);

      //! Assignment operator with any expression right-hand-side
      template <typename S> 
      Expr< SFadExprTag<T,Num> >& operator=(const Expr<S>& x); 

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Addition-assignment operator with constant right-hand-side
      Expr< SFadExprTag<T,Num> >& operator += (const T& x);

      //! Subtraction-assignment operator with constant right-hand-side
      Expr< SFadExprTag<T,Num> >& operator -= (const T& x);

      //! Multiplication-assignment operator with constant right-hand-side
      Expr< SFadExprTag<T,Num> >& operator *= (const T& x);

      //! Division-assignment operator with constant right-hand-side
      Expr< SFadExprTag<T,Num> >& operator /= (const T& x);

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S> 
      Expr< SFadExprTag<T,Num> >& operator += (const Expr<S>& x);

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S> 
      Expr< SFadExprTag<T,Num> >& operator -= (const Expr<S>& x);
  
      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S> 
      Expr< SFadExprTag<T,Num> >& operator *= (const Expr<S>& x);

      //! Division-assignment operator with Expr right-hand-side
      template <typename S> 
      Expr< SFadExprTag<T,Num> >& operator /= (const Expr<S>& x);

      //@}

    protected:

      //! Value
      T val_;

      //! Derivatives
      T dx_[Num];

    }; // class Expr<SFadExprTag>

    //! Forward-mode AD class using static memory allocation
    /*!
     * This is the user-level class for forward mode AD with static
     * memory allocation, and is appropriate for whenever the number
     * of derivative components is known at compile time.  The size
     * of the derivative array is fixed by the template parameter \c Num.  
     *
     * The class is templated on two types, \c ValueT and \c ScalarT.  Type
     * \c ValueT is the type for values the derivative class holds, while
     * type \c ScalarT is the type of basic scalars in the code being
     * differentiated (usually \c doubles).  When computing first derivatives, 
     * these two types are generally the same,  However when computing
     * higher derivatives, \c ValueT may be SFad<double> while \c ScalarT will
     * still be \c double.  Usually \c ScalarT does not need to be explicitly
     * specified since it can be deduced from \c ValueT through the template
     * metafunction ScalarValueType.
     */
    template <typename ValueT, int Num,
	      typename ScalarT = typename ScalarValueType<ValueT>::type >
    class SFad : 
      public Expr< SFadExprTag<ValueT,Num > > {

    public:

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      SFad() : 
	Expr< SFadExprTag< ValueT,Num > >() {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      SFad(const ValueT & x) : 
	Expr< SFadExprTag< ValueT,Num > >(x) {}

      //! Constructor with supplied value \c x of type ScalarT
      /*!
       * Initializes value to \c ValueT(x) and derivative array is empty
       */
      SFad(const ScalarT& x) : 
	Expr< SFadExprTag< ValueT,Num > >(ValueT(x)) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      SFad(const int sz, const ValueT & x) : 
	Expr< SFadExprTag< ValueT,Num > >(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      SFad(const int sz, const int i, const ValueT & x) : 
	Expr< SFadExprTag< ValueT,Num > >(sz,i,x) {}

      //! Copy constructor
      SFad(const SFad& x) : 
	Expr< SFadExprTag< ValueT,Num > >(x) {}

      //! Copy constructor from any Expression object
      template <typename S> SFad(const Expr<S>& x) : 
	Expr< SFadExprTag< ValueT,Num > >(x) {}

      //@}

      //! Destructor
      ~SFad() {}

      //! Assignment operator with constant right-hand-side
      SFad& operator=(const ValueT& v) {
	Expr< SFadExprTag< ValueT,Num > >::operator=(v);
	return *this;
      }

      //! Assignment operator with constant right-hand-side
      SFad& operator=(const ScalarT& v) {
	Expr< SFadExprTag< ValueT,Num > >::operator=(ValueT(v));
	return *this;
      }

      //! Assignment operator with DFad right-hand-side
      SFad& operator=(const SFad& x) {
	Expr< SFadExprTag< ValueT,Num > >::operator=(static_cast<const Expr< SFadExprTag< ValueT,Num > >&>(x));
	return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S> SFad& operator=(const Expr<S>& x) 
      {
	Expr< SFadExprTag< ValueT,Num > >::operator=(x);
	return *this;
      }

    }; // class SFad<ValueT,Num,ScalarT>

    //! Forward-mode AD class using static memory allocation
    /*!
     * This is the specialization of SFad<ValueT,ScalarT> for when
     * \c ValueT and \c ScalarT are the same type.  It removes an extra
     * constructor that would be duplicated in this case.
     */
    template <typename ValueT, int Num>
    class SFad<ValueT,Num,ValueT> : 
      public Expr< SFadExprTag<ValueT,Num >  >{

    public:

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      SFad() : Expr< SFadExprTag< ValueT,Num > >() {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      SFad(const ValueT & x) : 
	Expr< SFadExprTag< ValueT,Num > >(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      SFad(const int sz, const ValueT & x) : 
	Expr< SFadExprTag< ValueT,Num > >(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      SFad(const int sz, const int i, const ValueT & x) : 
	Expr< SFadExprTag< ValueT,Num > >(sz,i,x) {}

      //! Copy constructor
      SFad(const SFad& x) : 
	Expr< SFadExprTag< ValueT,Num > >(x) {}

      //! Copy constructor from any Expression object
      template <typename S> SFad(const Expr<S>& x) : 
	Expr< SFadExprTag< ValueT,Num > >(x) {}

      //@}

      //! Destructor
      ~SFad() {}

      //! Assignment operator with constant right-hand-side
      SFad& operator=(const ValueT& v) {
	Expr< SFadExprTag< ValueT,Num > >::operator=(v);
	return *this;
      }

      //! Assignment operator with DFad right-hand-side
      SFad& operator=(const SFad& x) {
	Expr< SFadExprTag< ValueT,Num > >::operator=(static_cast<const Expr< SFadExprTag< ValueT,Num > >&>(x));
	return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S> SFad& operator=(const Expr<S>& x) 
      {
	Expr< SFadExprTag< ValueT,Num > >::operator=(x);
	return *this;
      }

    }; // class SFad<ValueT,Num>

  } // namespace Fad

} // namespace Sacado

#include "Sacado_Fad_SFadImp.hpp"
#include "Sacado_Fad_Ops.hpp"

#endif // SACADO_FAD_SFAD_HPP
