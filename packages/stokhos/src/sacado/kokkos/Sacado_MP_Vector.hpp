// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_MP_VECTOR_HPP
#define SACADO_MP_VECTOR_HPP

#include "Stokhos_ConfigDefs.h"

#ifdef HAVE_STOKHOS_SACADO

#include <ostream>	// for std::ostream

#include "Stokhos_mpl_for_each.hpp"

namespace Sacado {

  //! Namespace for multipoint classes
  namespace MP {

    //! Wrapper for a generic expression template
    /*!
     * This class is used to limit the overload set for building up 
     * expressions.  Each expression object should derive from this
     * using CRTP:
     *
     * \code
     * class T : public Expr<T> { ... };
     * \endcode
     *
     * In this case the default implementation here should be correct for
     * any expression class.  If not, an expression class is free to change
     * the implementation through partial specialization.
     */
    template <typename T, typename node> class Expr {
    public:

      //! Node type
      typedef node node_type;

      //! Typename of derived object, returned by derived()
      /*!
       * This assumes a CRTP pattern where T is infact derived from
       * Expr<T>
       */
      typedef T derived_type;

      //! Return derived object
      /*!
       * This assumes a CRTP pattern where T is infact derived from
       * Expr<T>.  This will only compile if this infact the case.
       */
      const derived_type& derived() const;

    };

    //! Vectorized evaluation class
    template <typename Storage, typename Node> 
    class Vector : public Expr< Vector<Storage,Node>, Node > {
    public:

      //! Typename of storage class
      typedef Storage storage_type;

      //! Node type
      typedef Node node_type;

      typedef typename storage_type::value_type value_type;
      typedef typename storage_type::ordinal_type ordinal_type;
      typedef typename storage_type::pointer pointer;
      typedef typename storage_type::const_pointer const_pointer;
      typedef typename storage_type::reference reference;
      typedef typename storage_type::const_reference const_reference;

      //! Typename of scalar's (which may be different from T)
      typedef typename ScalarType<value_type>::type scalar_type;     

      //! Turn Vector into a meta-function class usable with mpl::apply
      template <typename S> struct apply {};

      //! Number of arguments
      static const int num_args = 1;

      //! Default constructor
      /*!
       * Sets size to 1 and first coefficient to 0 (represents a constant).
       */
      Vector();

      //! Constructor with supplied value \c x
      /*!
       * Sets size to 1 and first coefficient to x (represents a constant).
       */
      Vector(const value_type& x);

      //! Constructor with specified size \c sz
      /*!
       * Creates array of size \c sz and initializes coeffiencts to 0.
       */
      Vector(ordinal_type sz, const value_type& x);

      //! Copy constructor
      Vector(const Vector& x);

      //! Copy constructor from any Expression object
      template <typename S> Vector(const Expr<S,node_type>& x);

      //! Destructor
      ~Vector() {}

      //! Initialize coefficients to value
      void init(const value_type& v);

      //! Initialize coefficients to an array of values
      void init(const value_type* v);

      //! Initialize coefficients from an Vector with different storage
      template <typename S, typename N>
      void init(const Vector<S,N>& v);

      //! Load coefficients to an array of values
      void load(value_type* v);

      //! Load coefficients into an Vector with different storage
      template <typename S, typename N>
      void load(Vector<S,N>& v);

      //! Reset size
      /*!
       * Coefficients are preserved.  
       */
      void reset(ordinal_type sz_new);

      //! Prepare vector for writing 
      /*!
       * This method prepares the vector for writing through coeff() and 
       * fastAccessCoeff() member functions.  It ensures the handle for the
       * coefficients is not shared among any other vector.  
       * If the handle is not shared it does nothing, so there
       * is no cost in calling this method in this case.  If the handle is 
       * shared and this method is not called, any changes to the coefficients
       * by coeff() or fastAccessCoeff() may change other vector objects.
       */
      void copyForWrite();

      //! Returns whether two ETV objects have the same values
      template <typename S>
      bool isEqualTo(const Expr<S,node_type>& xx) const;

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      Vector& operator=(const value_type& val);

      //! Assignment with Vector right-hand-side
      Vector& operator=(const Vector& x);

      //! Assignment with any expression right-hand-side
      template <typename S> 
      Vector& operator=(const Expr<S,node_type>& x);

      //@}

      /*!
       * Accessor methods
       */

      /*!
       * @name Value accessor methods
       */
      //@{

      //! Returns storage object
      const storage_type& storage() const;

      //! Returns storage object
      storage_type& storage();

      //! Returns value
      const_reference val() const;

      //! Returns value
      reference val();

      //@}

      /*!
       * @name Coefficient accessor methods
       */
      //@{

      //! Returns size of polynomial
      ordinal_type size() const;

      //! Returns true if polynomial has size >= sz
      bool hasFastAccess(ordinal_type sz) const;

      //! Returns Hermite coefficient array
      const_pointer coeff() const;

      //! Returns Hermite coefficient array
      pointer coeff();

      //! Returns degree \c i term with bounds checking
      value_type coeff(ordinal_type i) const;
    
      //! Returns degree \c i term without bounds checking
      reference fastAccessCoeff(ordinal_type i);

      //! Returns degree \c i term without bounds checking
      value_type fastAccessCoeff(ordinal_type i) const;
    
      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Addition-assignment operator with constant right-hand-side
      Vector& operator += (const value_type& x);

      //! Subtraction-assignment operator with constant right-hand-side
      Vector& operator -= (const value_type& x);

      //! Multiplication-assignment operator with constant right-hand-side
      Vector& operator *= (const value_type& x);

      //! Division-assignment operator with constant right-hand-side
      Vector& operator /= (const value_type& x);

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S> 
      Vector& operator += (const Expr<S,node_type>& x);

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S> 
      Vector& operator -= (const Expr<S,node_type>& x);
  
      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S> 
      Vector& operator *= (const Expr<S,node_type>& x);

      //! Division-assignment operator with Expr right-hand-side
      template <typename S> 
      Vector& operator /= (const Expr<S,node_type>& x);

      //@}

      std::string name() const;

    }; // class Vector

    //! Type for storing nodes in expression graph
    /*!
     * Since expression nodes are returned by value in the overloaded
     * operators, we can't store them by reference in general.
     */
    template <typename T> struct const_expr_ref {
      typedef const T type;
    };

    //! Type for storing nodes in expression graph
    /*!
     * Specialization for leaf-nodes, which can be stored by reference
     * since they are an argument to the expression.
     */
    template <typename S, typename N> struct const_expr_ref< Vector<S,N> > {
      typedef const Vector<S,N>& type;
    };

    template <typename T, typename node> class UnaryPlusOp {};
    template <typename T, typename node> class UnaryMinusOp {};
    template <typename T, typename node> class ExpOp {};
    template <typename T, typename node> class LogOp {};
    template <typename T, typename node> class Log10Op {};
    template <typename T, typename node> class SqrtOp {};
    template <typename T, typename node> class CosOp {};
    template <typename T, typename node> class SinOp {};
    template <typename T, typename node> class TanOp {};
    template <typename T, typename node> class ACosOp {};
    template <typename T, typename node> class ASinOp {};
    template <typename T, typename node> class ATanOp {};
    template <typename T, typename node> class CoshOp {};
    template <typename T, typename node> class SinhOp {};
    template <typename T, typename node> class TanhOp {};
    template <typename T, typename node> class ACoshOp {};
    template <typename T, typename node> class ASinhOp {};
    template <typename T, typename node> class ATanhOp {};
    template <typename T, typename node> class AbsOp {};
    template <typename T, typename node> class FAbsOp {};

    template <typename T1, typename T2, typename node> class AdditionOp {};
    template <typename T1, typename T2, typename node> class SubtractionOp {};
    template <typename T1, typename T2, typename node> class MultiplicationOp {};
    template <typename T1, typename T2, typename node> class DivisionOp {};
    template <typename T1, typename T2, typename node> class Atan2Op {};
    template <typename T1, typename T2, typename node> class PowerOp {};
    template <typename T1, typename T2, typename node> class MaxOp {};
    template <typename T1, typename T2, typename node> class MinOp {};

  } // namespace MP

} // namespace Sacado

#include "Sacado_MP_VectorTraits.hpp"

// Host specialization
#include "KokkosArray_Host.hpp"
#include "KokkosArray_Host_macros.hpp"
#include "Sacado_MP_Vector_impl.hpp"
#include "KokkosArray_Clear_macros.hpp"

// Cuda specialization
#include "KokkosArray_Cuda.hpp"
#include "KokkosArray_Cuda_macros.hpp"
#include "Sacado_MP_Vector_impl.hpp"
#include "KokkosArray_Clear_macros.hpp"

#endif // HAVE_STOKHOS_SACADO

#endif // SACADO_MP_VECTOR_HPP
