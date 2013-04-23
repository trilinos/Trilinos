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

#ifndef SACADO_ELRFAD_GENERALFAD_HPP
#define SACADO_ELRFAD_GENERALFAD_HPP

#include "Sacado_ELRFad_Expression.hpp"
#include "Sacado_dummy_arg.hpp"
#include<ostream>

namespace Sacado {

  //! Namespace for expression-level reverse forward-mode AD classes
  namespace ELRFad {

    //! Forward-mode AD class templated on the storage for the derivative array
    /*!
     * This class provides a general forward mode AD implementation for any
     * type of derivative array storage.  It does not incorporate expression
     * templates.
     */
    template <typename T, typename Storage> 
    class GeneralFad : public Storage {

    public:

      //! Typename of values
      typedef T value_type;

      //! Typename of scalar's (which may be different from T)
      typedef typename ScalarType<T>::type scalar_type;

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor
      GeneralFad() : Storage(T(0.)), update_val_(true) {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      GeneralFad(const T & x) : Storage(x), update_val_(true) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      GeneralFad(const int sz, const T & x) : 
	Storage(sz, x), update_val_(true) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      GeneralFad(const int sz, const int i, const T & x) : 
	Storage(sz, x), update_val_(true) { 
	this->fastAccessDx(i)=1.; 
      }

      //! Copy constructor
      GeneralFad(const GeneralFad& x) : 
	Storage(x), update_val_(x.update_val_) {}

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

      //! Set whether this Fad object should update values
      void setUpdateValue(bool update_val) { update_val_ = update_val; }

      //! Return whether this Fad object has an updated value
      bool updateValue() const { return update_val_; }

      //! Returns whether two Fad objects have the same values
      template <typename S>
      bool isEqualTo(const Expr<S>& x) const {
	typedef IsEqual<value_type> IE;
	if (x.size() != this->size()) return false;
	bool eq = IE::eval(x.val(), this->val());
	for (int i=0; i<this->size(); i++)
	  eq = eq && IE::eval(x.dx(i), this->dx(i));
	return eq;
      }

      //@}

      /*!
       * @name Derivative accessor methods
       */
      //@{

      /*! 
       * \brief Returns number of derivative components that can be stored 
       * without reallocation
       */
      int availableSize() const { return this->length(); }

      //! Returns true if derivative array is not empty
      bool hasFastAccess() const { return this->size()!=0;}

      //! Returns true if derivative array is empty
      bool isPassive() const { return this->size()==0;}
      
      //! Set whether variable is constant
      void setIsConstant(bool is_const) { 
	if (is_const && this->size()!=0)
	  this->resize(0);
      }
    
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

      //! Addition-assignment operator with constant right-hand-side
      /*!
       * Creates a dummy overload when value_type and scalar_type are the 
       * same type.
       */
      GeneralFad& operator += (const typename dummy<value_type,scalar_type>::type& x);

      //! Subtraction-assignment operator with constant right-hand-side
      /*!
       * Creates a dummy overload when value_type and scalar_type are the 
       * same type.
       */
      GeneralFad& operator -= (const typename dummy<value_type,scalar_type>::type& x);

      //! Multiplication-assignment operator with constant right-hand-side
      /*!
       * Creates a dummy overload when value_type and scalar_type are the 
       * same type.
       */
      GeneralFad& operator *= (const typename dummy<value_type,scalar_type>::type& x);

      //! Division-assignment operator with constant right-hand-side
      /*!
       * Creates a dummy overload when value_type and scalar_type are the 
       * same type.
       */
      GeneralFad& operator /= (const typename dummy<value_type,scalar_type>::type& x);

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

      //! Update value
      bool update_val_;

      // Functor for mpl::for_each to compute the local accumulation
      // of a tangent derivative
      template <typename ExprT>
      struct FastLocalAccumOp {
	typedef typename ExprT::value_type value_type;
	static const int N = ExprT::num_args;
	const ExprT& x;
	mutable value_type t;
	value_type partials[N];
	const typename ExprT::base_expr_type* args[N];
	int i;
	inline FastLocalAccumOp(const ExprT& x_) : x(x_) { 
	  x.computePartials(value_type(1.), partials); 
	  for (int j=0; j<N; j++)
	    args[j] = &(x.getArg(j));
	}
	template <typename ArgT>
	inline void operator () (ArgT arg) const {
	  const int Arg = ArgT::value;
	  t += partials[Arg] * args[Arg]->fastAccessDx(i);
	}
      };

      template <typename ExprT>
      struct SlowLocalAccumOp : FastLocalAccumOp<ExprT> {
	inline SlowLocalAccumOp(const ExprT& x_) : 
	  FastLocalAccumOp<ExprT>(x_) {}
	template <typename ArgT>
	inline void operator () (ArgT arg) const {
	  const int Arg = ArgT::value;
	  if (this->x.template isActive<Arg>())
	    this->t += this->partials[Arg] * this->args[Arg]->fastAccessDx(this->i);
	}
      };

    }; // class GeneralFad

    
    template <typename T, typename Storage>
    std::ostream& operator << (std::ostream& os, 
                               const GeneralFad<T,Storage>& x) {
      os << x.val() << " [";
      
      for (int i=0; i< x.size(); i++) {
        os << " " << x.dx(i);
      }

      os << " ]";
      return os;
    }

  } // namespace ELRFad

} // namespace Sacado

#include "Sacado_ELRFad_GeneralFadImp.hpp"

#endif // SACADO_ELRFAD_GENERALFAD_HPP
