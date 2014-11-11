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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
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

#ifndef SACADO_CACHEFAD_SFAD_HPP
#define SACADO_CACHEFAD_SFAD_HPP

#include "Sacado_CacheFad_SFadTraits.hpp"
#include "Sacado_CacheFad_Expression.hpp"
#include "Sacado_StaticArrayTraits.hpp"
#include "Sacado_dummy_arg.hpp"

namespace Sacado {

  //! Namespace for forward-mode AD classes
  namespace CacheFad {

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
      typedef typename RemoveConst<T>::type value_type;

      //! Typename of scalar's (which may be different from T)
      typedef typename ScalarType<value_type>::type scalar_type;

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor
      KOKKOS_INLINE_FUNCTION
      Expr() : val_( T(0.)), update_val_(true) { ss_array<T>::zero(dx_, Num); }

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      KOKKOS_INLINE_FUNCTION
      Expr(const T & x) : val_(x), update_val_(true)  {
        ss_array<T>::zero(dx_, Num); }

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      KOKKOS_INLINE_FUNCTION
      Expr(const int sz, const T & x);

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      KOKKOS_INLINE_FUNCTION
      Expr(const int sz, const int i, const T & x);

      //! Copy constructor
      KOKKOS_INLINE_FUNCTION
      Expr(const Expr& x);

      //! Copy constructor from any Expression object
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Expr(const Expr<S>& x);

      //! Destructor
      KOKKOS_INLINE_FUNCTION
      ~Expr() {}

      //! Set %Fad object as the \c ith independent variable
      /*!
       * Sets the derivative array of length \c n to the \c ith row of the
       * identity matrix and has the same affect as the
       * Implementation(const int sz, const int i, const T & x)
       * constructor.
       */
      KOKKOS_INLINE_FUNCTION
      void diff(const int ith, const int n);

      //! Resize derivative array to length \c sz
      /*!
       * Since the derivative array length is not dynamic, this method
       * throws an error if compiled with SACADO_DEBUG defined.
       */
      KOKKOS_INLINE_FUNCTION
      void resize(int sz);

      //! Expand derivative array to size sz
      /*!
       * Since the derivative array length is not dynamic, this method
       * throws an error if compiled with SACADO_DEBUG defined.
       */
      KOKKOS_INLINE_FUNCTION
      void expand(int sz) { resize(sz); }

      //! Zero out the derivative array
      KOKKOS_INLINE_FUNCTION
      void zero() { ss_array<T>::zero(dx_, Num); }

      //! Set whether this Fad object should update values
      KOKKOS_INLINE_FUNCTION
      void setUpdateValue(bool update_val) { update_val_ = update_val; }

      //! Return whether this Fad object has an updated value
      KOKKOS_INLINE_FUNCTION
      bool updateValue() const { return update_val_; }

      //! Cache values
      KOKKOS_INLINE_FUNCTION
      void cache() const {}

      //! Returns whether two Fad objects have the same values
      template <typename S>
      KOKKOS_INLINE_FUNCTION
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
       * @name Value accessor methods
       */
      //@{

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      const T& val() const { return val_;}

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      T& val() { return val_;}

      //@}

      /*!
       * @name Derivative accessor methods
       */
      //@{

      //! Returns number of derivative components
      KOKKOS_INLINE_FUNCTION
      int size() const { return Num;}

      /*!
       * \brief Returns number of derivative components that can be stored
       * without reallocation
       */
      KOKKOS_INLINE_FUNCTION
      int availableSize() const { return Num; }

      //! Returns true if derivative array is not empty
      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess() const { return true; }

      //! Returns true if derivative array is empty
      KOKKOS_INLINE_FUNCTION
      bool isPassive() const { return false; }

      //! Set whether variable is constant
      KOKKOS_INLINE_FUNCTION
      void setIsConstant(bool is_const) {}

      //! Returns derivative array
      KOKKOS_INLINE_FUNCTION
      const T* dx() const { return &(dx_[0]);}

      //! Returns derivative component \c i with bounds checking
      KOKKOS_INLINE_FUNCTION
      const T& dx(int i) const { return dx_[i]; }

      //! Returns derivative component \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      T& fastAccessDx(int i) { return dx_[i];}

      //! Returns derivative component \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      const T& fastAccessDx(int i) const { return dx_[i];}

      //@}

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      Expr< SFadExprTag<T,Num> >& operator=(const T& val);

      //! Assignment with Expr right-hand-side
      KOKKOS_INLINE_FUNCTION
      Expr< SFadExprTag<T,Num> >&
      operator=(const Expr< SFadExprTag<T,Num> >& x);

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Expr< SFadExprTag<T,Num> >& operator=(const Expr<S>& x);

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Addition-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      Expr< SFadExprTag<T,Num> >& operator += (const T& x);

      //! Subtraction-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      Expr< SFadExprTag<T,Num> >& operator -= (const T& x);

      //! Multiplication-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      Expr< SFadExprTag<T,Num> >& operator *= (const T& x);

      //! Division-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      Expr< SFadExprTag<T,Num> >& operator /= (const T& x);

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Expr< SFadExprTag<T,Num> >& operator += (const Expr<S>& x);

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Expr< SFadExprTag<T,Num> >& operator -= (const Expr<S>& x);

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Expr< SFadExprTag<T,Num> >& operator *= (const Expr<S>& x);

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Expr< SFadExprTag<T,Num> >& operator /= (const Expr<S>& x);

      //@}

    protected:

      //! Value
      T val_;

      //! Derivatives
      T dx_[Num];

      //! Update value
      bool update_val_;

    }; // class Expr<SFadExprTag>

    /*!
     * \brief Forward-mode AD class using static memory allocation and
     * caching expression templates.
     */
    /*!
     * This is the user-level class for forward mode AD with static
     * memory allocation, and is appropriate for whenever the number
     * of derivative components is known at compile time.  The size
     * of the derivative array is fixed by the template parameter \c Num.
     * It is similar to Sacado::Fad::SFad, except it uses the
     * caching expression templates that cache the results of val()
     * calculations for later dx() calculations.
     */
    template <typename ValueT, int Num>
    class SFad :
      public Expr< SFadExprTag<ValueT,Num > > {

    public:

      //! Typename of scalar's (which may be different from ValueT)
      typedef typename ScalarType<ValueT>::type ScalarT;

      //! Turn SFad into a meta-function class usable with mpl::apply
      template <typename T>
      struct apply {
        typedef SFad<T,Num> type;
      };

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      KOKKOS_INLINE_FUNCTION
      SFad() :
        Expr< SFadExprTag< ValueT,Num > >() {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      KOKKOS_INLINE_FUNCTION
      SFad(const ValueT & x) :
        Expr< SFadExprTag< ValueT,Num > >(x) {}

      //! Constructor with supplied value \c x of type ScalarT
      /*!
       * Initializes value to \c ValueT(x) and derivative array is empty.
       * Creates a dummy overload when ValueT and ScalarT are the same type.
       */
      KOKKOS_INLINE_FUNCTION
      SFad(const typename dummy<ValueT,ScalarT>::type& x) :
        Expr< SFadExprTag< ValueT,Num > >(ValueT(x)) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      KOKKOS_INLINE_FUNCTION
      SFad(const int sz, const ValueT & x) :
        Expr< SFadExprTag< ValueT,Num > >(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      KOKKOS_INLINE_FUNCTION
      SFad(const int sz, const int i, const ValueT & x) :
        Expr< SFadExprTag< ValueT,Num > >(sz,i,x) {}

      //! Copy constructor
      KOKKOS_INLINE_FUNCTION
      SFad(const SFad& x) :
        Expr< SFadExprTag< ValueT,Num > >(static_cast<const Expr< SFadExprTag< ValueT,Num > >&>(x)) {}

      //! Copy constructor from any Expression object
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SFad(const Expr<S>& x) :
        Expr< SFadExprTag< ValueT,Num > >(x) {}

      //@}

      //! Destructor
      KOKKOS_INLINE_FUNCTION
      ~SFad() {}

      //! Assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      SFad& operator=(const ValueT& v) {
        Expr< SFadExprTag< ValueT,Num > >::operator=(v);
        return *this;
      }

      //! Assignment operator with constant right-hand-side
      /*!
       * Creates a dummy overload when ValueT and ScalarT are the same type.
       */
      KOKKOS_INLINE_FUNCTION
      SFad& operator=(const typename dummy<ValueT,ScalarT>::type& v) {
        Expr< SFadExprTag< ValueT,Num > >::operator=(ValueT(v));
        return *this;
      }

      //! Assignment operator with DFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      SFad& operator=(const SFad& x) {
        Expr< SFadExprTag< ValueT,Num > >::operator=(static_cast<const Expr< SFadExprTag< ValueT,Num > >&>(x));
        return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SFad& operator=(const Expr<S>& x)
      {
        Expr< SFadExprTag< ValueT,Num > >::operator=(x);
        return *this;
      }

    }; // class SFad<ValueT,Num>

  } // namespace CacheFad

} // namespace Sacado

#include "Sacado_CacheFad_SFadImp.hpp"
#include "Sacado_CacheFad_Ops.hpp"

//
// Classes needed for Kokkos::View< SFad<...> ... > specializations
//
// Users can disable these view specializations either at configure time or
// by defining SACADO_DISABLE_FAD_VIEW_SPEC in their code.
//

#if defined(HAVE_SACADO_KOKKOSCORE) && defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#include "impl/Kokkos_AnalyzeShape.hpp"
#include "Kokkos_AnalyzeSacadoShape.hpp"
#include "Sacado_CacheFad_ViewFad.hpp"

namespace Kokkos {
namespace Impl {

// Forward declarations
struct ViewSpecializeSacadoFad;
template <typename T,unsigned,unsigned> struct ViewFadType;

//! The View Fad type associated with this type
template< class ValueType, int N, unsigned length, unsigned stride >
struct ViewFadType< Sacado::CacheFad::SFad< ValueType, N >, length, stride > {
  typedef Sacado::CacheFad::ViewFad<ValueType,length,stride> type;
};

//! The View Fad type associated with this type
template< class ValueType, int N, unsigned length, unsigned stride >
struct ViewFadType< const Sacado::CacheFad::SFad< ValueType, N >, length, stride > {
  typedef Sacado::CacheFad::ViewFad<const ValueType,length,stride> type;
};

/** \brief  Analyze the array shape of a Sacado::CacheFad::SFad<T,N>.
 *
 *  This specialization is required so that the array shape of
 *  Kokkos::View< Sacado::CacheFad::SFad<T,N>, ... >
 *  can be determined at compile-time.
 *
 *  We add one to the SFad dimension (N) to store the value component.
 */
template< class ValueType, int N >
struct AnalyzeShape< Sacado::CacheFad::SFad< ValueType, N > >
  : Shape< sizeof(Sacado::CacheFad::SFad< ValueType, N >) , 0 > // Treat as a scalar
{
public:

  typedef ViewSpecializeSacadoFad specialize ;

  typedef Shape< sizeof(Sacado::CacheFad::SFad< ValueType, N >) , 0 > shape ;

  typedef       Sacado::CacheFad::SFad< ValueType, N >        array_intrinsic_type ;
  typedef const Sacado::CacheFad::SFad< ValueType, N >  const_array_intrinsic_type ;
  typedef array_intrinsic_type non_const_array_intrinsic_type ;

  typedef       Sacado::CacheFad::SFad< ValueType, N >  type ;
  typedef const Sacado::CacheFad::SFad< ValueType, N >  const_type ;
  typedef       Sacado::CacheFad::SFad< ValueType, N >  non_const_type ;

  typedef       Sacado::CacheFad::SFad< ValueType, N >  value_type ;
  typedef const Sacado::CacheFad::SFad< ValueType, N >  const_value_type ;
  typedef       Sacado::CacheFad::SFad< ValueType, N >  non_const_value_type ;
};

/** \brief  Analyze the array shape of a Sacado::CacheFad::SFad<T,N>.
 *
 *  This specialization is required so that the array shape of
 *  Kokkos::View< Sacado::CacheFad::SFad<T,N>, ... >
 *  can be determined at compile-time.
 *
 *  We add one to the SFad dimension (N) to store the value component.
 */
template< class ValueType, class Layout, int N >
struct AnalyzeSacadoShape< Sacado::CacheFad::SFad< ValueType, N >, Layout >
  : ShapeInsert< typename AnalyzeSacadoShape< ValueType, Layout >::shape , N+1 >::type
{
private:

  typedef AnalyzeSacadoShape< ValueType, Layout > nested ;

public:

  typedef ViewSpecializeSacadoFad specialize ;

  typedef typename ShapeInsert< typename nested::shape , N+1 >::type shape ;

  typedef typename nested::array_intrinsic_type         array_intrinsic_type [N+1];
  typedef typename nested::const_array_intrinsic_type   const_array_intrinsic_type [N+1] ;
  typedef array_intrinsic_type non_const_array_intrinsic_type ;

  typedef       Sacado::CacheFad::SFad< ValueType, N >  type ;
  typedef const Sacado::CacheFad::SFad< ValueType, N >  const_type ;
  typedef       Sacado::CacheFad::SFad< ValueType, N >  non_const_type ;

  typedef       Sacado::CacheFad::SFad< ValueType, N >  value_type ;
  typedef const Sacado::CacheFad::SFad< ValueType, N >  const_value_type ;
  typedef       Sacado::CacheFad::SFad< ValueType, N >  non_const_value_type ;

  typedef typename nested::type           flat_array_type ;
  typedef typename nested::const_type     const_flat_array_type ;
  typedef typename nested::non_const_type non_const_flat_array_type ;
};

} // namespace Impl
} // namespace Kokkos

#endif

#endif // SACADO_CACHEFAD_SFAD_HPP
