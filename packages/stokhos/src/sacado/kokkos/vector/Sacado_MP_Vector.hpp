// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef SACADO_MP_VECTOR_HPP
#define SACADO_MP_VECTOR_HPP

#include "Stokhos_ConfigDefs.h"

#ifdef HAVE_STOKHOS_SACADO

#include <ostream>      // for std::ostream

#include "Kokkos_Macros.hpp"

#include "Sacado_MP_ExpressionTraits.hpp"
#include "Sacado_MP_VectorTraits.hpp"
#include "Sacado_Traits.hpp"
#include "Sacado_mpl_apply.hpp"
#include "Sacado_mpl_range_c.hpp"
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
    template <typename T>
    class Expr {
    public:

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
      KOKKOS_INLINE_FUNCTION
      const derived_type& derived() const {
        return static_cast<const derived_type&>(*this);
      }

    };

    //! Vectorized evaluation class
    template <typename Storage>
    class Vector : public Expr< Vector<Storage> > {
    public:

      //! Typename of storage class
      typedef Storage storage_type;

      typedef typename storage_type::value_type value_type;
      typedef typename storage_type::ordinal_type ordinal_type;
      typedef typename storage_type::pointer pointer;
      typedef typename storage_type::const_pointer const_pointer;
      typedef typename storage_type::reference reference;
      typedef typename storage_type::const_reference const_reference;

      //! Typename of scalar's (which may be different from T)
      typedef typename ScalarType<value_type>::type scalar_type;

      //! Turn Vector into a meta-function class usable with mpl::apply
      template < class NewStorageType >
      struct apply {
        typedef Vector< NewStorageType > type;
      };

      //! Number of arguments
      static const int num_args = 1;

      // A temporary hack to allow taking the address of a temporary
      // Vector with ViewStorage.  A better approach would be to return
      // a VectorViewStoragePtr with overloaded * to return a new
      // Vector<ViewStorage>
      KOKKOS_INLINE_FUNCTION
      Vector* operator&() { return this; }
      KOKKOS_INLINE_FUNCTION
      const Vector* operator&() const { return this; }

      //! Default constructor
      /*!
       * Sets size to 1 and first coefficient to 0 (represents a constant).
       */
      KOKKOS_INLINE_FUNCTION
      Vector() : s(1) {}

      //! Constructor with supplied value \c x
      /*!
       * Sets size to 1 and first coefficient to x (represents a constant).
       */
      KOKKOS_INLINE_FUNCTION
      Vector(const value_type& x) : s(1) { s.init(x); }

      //! Constructor with specified size \c sz
      /*!
       * Creates array of size \c sz and initializes coeffiencts to 0.
       */
      KOKKOS_INLINE_FUNCTION
      Vector(ordinal_type sz, const value_type& x) : s(sz,x) {}

      //! Constructor with supplied storage
      KOKKOS_INLINE_FUNCTION
      Vector(const storage_type& ss) : s(ss) {}

      //! Copy constructor
      KOKKOS_INLINE_FUNCTION
      Vector(const Vector& x) : s(x.s) {}

      //! Copy constructor from any Expression object
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Vector(const Expr<S>& xx) :
        s(xx.derived().size()) {
        typedef typename Expr<S>::derived_type expr_type;
        const expr_type& x = xx.derived();

        if (x.hasFastAccess(s.size())) {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] = x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] = x.coeff(i);
        }
      }

      //! Destructor
      KOKKOS_INLINE_FUNCTION
      ~Vector() {}

      //! Initialize coefficients to value
      KOKKOS_INLINE_FUNCTION
      void init(const value_type& v) { s.init(v); }

      //! Initialize coefficients to an array of values
      KOKKOS_INLINE_FUNCTION
      void init(const value_type* v) { s.init(v); }

      //! Initialize coefficients from an Vector with different storage
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      void init(const Vector<S>& v) {
        s.init(v.s.coeff(), v.s.size());
      }

      //! Load coefficients to an array of values
      KOKKOS_INLINE_FUNCTION
      void load(value_type* v) { s.load(v); }

      //! Load coefficients into an Vector with different storage
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      void load(Vector<S>& v) { s.load(v.s.coeff()); }

      //! Reset size
      /*!
       * Coefficients are preserved.
       */
      KOKKOS_INLINE_FUNCTION
      void reset(ordinal_type sz_new) {
        ordinal_type sz = this->size();
        s.resize(sz_new);
        if (sz == 1 && sz_new > sz)
          for (ordinal_type i=1; i<sz_new; i++)
            s[i] = s[0];
      }

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
      KOKKOS_INLINE_FUNCTION
      void copyForWrite() {  }

      //! Returns whether two ETV objects have the same values
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      bool isEqualTo(const Expr<S>& xx) const {
        const typename Expr<S>::derived_type& x = xx.derived();
        typedef IsEqual<value_type> IE;
        if (x.size() != this->size()) return false;
        bool eq = true;
        for (ordinal_type i=0; i<this->size(); i++)
          eq = eq && IE::eval(x.coeff(i), this->coeff(i));
        return eq;
      }

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator=(const value_type& x) {
        s.init(x);
        return *this;
      }

      //! Assignment with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator=(const Vector& x) {
        s = x.s;
        return *this;
      }

      //! Assignment with any expression right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Vector& operator=(const Expr<S>& xx) {
        typedef typename Expr<S>::derived_type expr_type;
        const expr_type& x = xx.derived();

        this->reset(x.size());
        if (x.hasFastAccess(s.size())) {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] = x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] = x.coeff(i);
        }
        return *this;
      }

      //! Assignment operator only valid for view storage
      template< typename S >
      KOKKOS_INLINE_FUNCTION
      typename Kokkos::Impl::enable_if<( ! Kokkos::Impl::is_same<S,void>::value &&
                                         Stokhos::is_ViewStorage<Storage>::value
                                       ), Vector >
        ::type const & operator = ( const Expr<S> & xx ) const
      {
        const typename Expr<S>::derived_type & x = xx.derived();

        for ( ordinal_type i = 0 ; i < s.size() ; ++i ) { s[i] = x.coeff(i); }

        return *this ;
      }

      //@}

      /*!
       * Accessor methods
       */

      //! Returns storage object
      KOKKOS_INLINE_FUNCTION
      const storage_type& storage() const { return s; }

      //! Returns storage object
      KOKKOS_INLINE_FUNCTION
      storage_type& storage() { return s; }

      /*!
       * @name Value accessor methods
       */
      //@{

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      const_reference val() const { return s[0]; }

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      reference val() { return s[0]; }

      //@}

      /*!
       * @name Coefficient accessor methods
       */
      //@{

      //! Returns size of polynomial
      KOKKOS_INLINE_FUNCTION
      ordinal_type size() const { return s.size();}

      //! Returns true if polynomial has size >= sz
      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess(ordinal_type sz) const { return s.size()>=sz;}

      //! Returns Hermite coefficient array
      KOKKOS_INLINE_FUNCTION
      const_pointer coeff() const { return s.coeff();}

      //! Returns Hermite coefficient array
      KOKKOS_INLINE_FUNCTION
      pointer coeff() { return s.coeff();}

      //! Returns degree \c i term with bounds checking
      KOKKOS_INLINE_FUNCTION
      value_type coeff(ordinal_type i) const {
        return i<s.size() ? s[i] : s[0]; }

      //! Returns degree \c i term with bounds checking
      KOKKOS_INLINE_FUNCTION
      value_type & coeff(ordinal_type i) {
        return i<s.size() ? s[i] : s[0]; }

      //! Returns degree \c i term without bounds checking
      KOKKOS_INLINE_FUNCTION
      reference fastAccessCoeff(ordinal_type i) { return s[i];}

      //! Returns degree \c i term without bounds checking
      KOKKOS_INLINE_FUNCTION
      value_type fastAccessCoeff(ordinal_type i) const { return s[i];}

      template <int i>
      KOKKOS_INLINE_FUNCTION
      value_type getCoeff() const { return s.template getCoeff<i>(); }

      template <int i>
      KOKKOS_INLINE_FUNCTION
      reference getCoeff() { return s.template getCoeff<i>(); }

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Addition-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator += (const value_type& x) {
        for (ordinal_type i=0; i<s.size(); i++)
          s[i] += x;
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator -= (const value_type& x) {
        for (ordinal_type i=0; i<s.size(); i++)
          s[i] -= x;
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator *= (const value_type& x) {
        for (ordinal_type i=0; i<s.size(); i++)
          s[i] *= x;
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator /= (const value_type& x) {
        for (ordinal_type i=0; i<s.size(); i++)
          s[i] /= x;
        return *this;
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Vector& operator += (const Expr<S>& x) {
        *this = *this + x;
        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Vector& operator -= (const Expr<S>& x) {
        *this = *this - x;
        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Vector& operator *= (const Expr<S>& x) {
        *this = *this * x;
        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Vector& operator /= (const Expr<S>& x) {
        *this = *this / x;
        return *this;
      }

      //@}

      KOKKOS_INLINE_FUNCTION
      std::string name() const { return "x"; }

    protected:

      Storage s;

      template <typename expr_type>
      struct StaticOp {
        storage_type& s;
        const expr_type& x;

        KOKKOS_INLINE_FUNCTION
        StaticOp(storage_type& s_, const expr_type& x_) : s(s_), x(x_) {}

        template <typename ArgT>
        KOKKOS_INLINE_FUNCTION
        void operator() (ArgT arg) const {
          const int Arg = ArgT::value;
          s.template getCoeff<Arg>() = x.template getCoeff<Arg>();
        }

      };

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
    template <typename S> struct const_expr_ref< Vector<S> > {
      typedef const Vector<S>& type;
    };

    template <typename Storage>
    std::ostream&
    operator << (std::ostream& os,
                 const Vector<Storage>& a)
    {
      typedef typename Vector<Storage>::ordinal_type ordinal_type;

      os << "[ ";

      for (ordinal_type i=0; i<a.size(); i++) {
        os << a.coeff(i) << " ";
      }

      os << "]\n";
      return os;
    }

  } // namespace MP

  //! Trait class to determine if a scalar type is a Vector
  template <typename T> struct is_mp_vector {
    static const bool value = false;
  };
  template <typename S> struct is_mp_vector< MP::Vector<S> > {
    static const bool value = true;
  };

} // namespace Sacado

#include "Sacado_MP_Vector_ops.hpp"

#endif // HAVE_STOKHOS_SACADO

#endif // SACADO_MP_VECTOR_HPP
