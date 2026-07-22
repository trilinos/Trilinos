// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SACADO_MP_VECTOR_HPP
#define SACADO_MP_VECTOR_HPP

#include "Stokhos_ConfigDefs.h"

#ifdef HAVE_STOKHOS_SACADO

#include <ostream>      // for std::ostream
#include <initializer_list>

#include "Kokkos_Macros.hpp"

#include "Sacado_MP_ExpressionTraits.hpp"
#include "Sacado_MP_VectorTraits.hpp"
#include "Sacado_Traits.hpp"
#include "Sacado_mpl_apply.hpp"
#include "Sacado_mpl_range_c.hpp"
#include "Stokhos_mpl_for_each.hpp"
#include "Stokhos_MemoryTraits.hpp"
#include "Stokhos_Is_Constant.hpp"
#include "Stokhos_ViewStorage.hpp"

#include "Kokkos_View_Utils.hpp"

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

      //! Return derived object
      /*!
       * This assumes a CRTP pattern where T is infact derived from
       * Expr<T>.  This will only compile if this infact the case.
       */
      KOKKOS_INLINE_FUNCTION
      const volatile derived_type& derived() const volatile {
        return static_cast<const volatile derived_type&>(*this);
      }

      // Allow explicit casting to integral types, since we don't have an
      // integral ensemble type.
      template <typename U, typename Enabled = typename std::enable_if<std::is_integral<U>::value>::type>
      KOKKOS_INLINE_FUNCTION
      explicit operator U() const { return static_cast<U>(derived().val()); }

    };

    //! Vectorized evaluation class
    template <typename Storage>
    class Vector : public Expr< Vector<Storage> > {
    public:

      //! Typename of storage class
      typedef Storage storage_type;

      typedef typename storage_type::value_type value_type;
      typedef typename storage_type::ordinal_type ordinal_type;
      typedef typename storage_type::execution_space execution_space;
      typedef typename storage_type::pointer pointer;
      typedef typename storage_type::volatile_pointer volatile_pointer;
      typedef typename storage_type::const_pointer const_pointer;
      typedef typename storage_type::const_volatile_pointer const_volatile_pointer;
      typedef typename storage_type::reference reference;
      typedef typename storage_type::volatile_reference volatile_reference;
      typedef typename storage_type::const_reference const_reference;
      typedef typename storage_type::const_volatile_reference const_volatile_reference;

      typedef typename execution_space::memory_space memory_space;
      typedef typename Stokhos::MemoryTraits<memory_space> MemTraits;

      //! Typename of scalar's (which may be different from value_type)
      typedef typename ScalarType<value_type>::type scalar_type;

      typedef Vector base_expr_type;

      //! Turn Vector into a meta-function class usable with mpl::apply
      template < class NewStorageType >
      struct apply {
        typedef Vector< NewStorageType > type;
      };

      //! Number of arguments
      static const int num_args = 1;

#if STOKHOS_ALIGN_MEMORY
      KOKKOS_INLINE_FUNCTION
      static void* operator new(std::size_t sz) {
        return MemTraits::alloc(sz);
      }
      KOKKOS_INLINE_FUNCTION
      static void* operator new[](std::size_t sz) {
        return MemTraits::alloc(sz);
      }
      KOKKOS_INLINE_FUNCTION
      static void* operator new(std::size_t sz, void* ptr) {
        return ptr;
      }
      KOKKOS_INLINE_FUNCTION
      static void* operator new[](std::size_t sz, void* ptr) {
        return ptr;
      }
      KOKKOS_INLINE_FUNCTION
      static void operator delete(void* ptr) {
        MemTraits::free(ptr);
      }
      KOKKOS_INLINE_FUNCTION
      static void operator delete[](void* ptr) {
        MemTraits::free(ptr);
      }
      KOKKOS_INLINE_FUNCTION
      static void operator delete(void* ptr, void*) {
        MemTraits::free(ptr);
      }
      KOKKOS_INLINE_FUNCTION
      static void operator delete[](void* ptr, void*) {
        MemTraits::free(ptr);
      }
#endif

      //! Default constructor
      /*!
       * May not intialize the coefficient array.
       */
      KOKKOS_DEFAULTED_FUNCTION
      Vector() = default;

      //! Constructor with supplied value \c x
      /*!
       * Sets size to 1 and first coefficient to x (represents a constant).
       */
      KOKKOS_INLINE_FUNCTION
      Vector(const value_type& x) : s(1,x) {}

      //! View constructor
      /*!
       * Creates vector with pre-allocated data.  Set \c owned = true
       * if this Vector should take over management of the data.
       */
      KOKKOS_INLINE_FUNCTION
      Vector(ordinal_type sz, pointer v, bool owned) : s(sz,v,owned) {}

      //! Constructor for creating a view out of pre-allocated memory
      /*!
       * This does not do any initialization of the coefficients.
       */
      KOKKOS_INLINE_FUNCTION
      Vector(ordinal_type sz, const value_type& x) : s(sz,x) {}

      //! Constructor with supplied storage
      KOKKOS_INLINE_FUNCTION
      Vector(const storage_type& ss) : s(ss) {}

      //! Copy constructor
      KOKKOS_DEFAULTED_FUNCTION
      Vector(const Vector& x) = default;

      //! Copy constructor
      KOKKOS_INLINE_FUNCTION
      Vector(const volatile Vector& x) : s(x.s) {}

      //! Copy constructor from any Expression object
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Vector(const Expr<S>& xx) :
        s(xx.derived().size()) {
        typedef typename Expr<S>::derived_type expr_type;
        const expr_type& x = xx.derived();

#ifdef STOKHOS_DEBUG
        if (s.size() != x.size())
          Kokkos::Impl::raise_error("Vector():  Mismatched sizes");
#endif

        if (x.hasFastAccess(s.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] = x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] = x.coeff(i);
        }
      }

      //! Intialize from initializer_list
      /*!
       * When static fixed storage is used:
       *    * passing an empty initializer list initializes all samples in the ensemble to zero;
       *    * passing an initializer list of length 1 initializes all samples in the ensemble to the passed value;
       *    * all other mismatches in size between the initializer list and the ensemble produce an abort.
       */
      KOKKOS_INLINE_FUNCTION
      Vector(std::initializer_list<value_type> l) : s(l.size(), l.begin()) {
        if constexpr (Storage::is_static) {
          const auto         lsz = static_cast<ordinal_type>(l.size());
          const ordinal_type  sz = this->size();
          if (lsz < sz) {
            if (lsz > 1) {
                Kokkos::abort("Size mismatch in list initialization of MP Vector with static fixed storage.");
            } 
            else {
              const value_type v = lsz > 0 ? *l.begin() : value_type(0);
              s.init(v);
            }
          }
        }
      }

      //! Destructor
      KOKKOS_DEFAULTED_FUNCTION
      ~Vector() = default;

      //! Initialize coefficients to value
      KOKKOS_INLINE_FUNCTION
      void init(const value_type& v) { s.init(v); }

      //! Initialize coefficients to value
      KOKKOS_INLINE_FUNCTION
      void init(const value_type& v) volatile { s.init(v); }

      //! Initialize coefficients to an array of values
      KOKKOS_INLINE_FUNCTION
      void init(const value_type* v) { s.init(v); }

      //! Initialize coefficients to an array of values
      KOKKOS_INLINE_FUNCTION
      void init(const value_type* v) volatile { s.init(v); }

      //! Initialize coefficients from an Vector with different storage
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      void init(const Vector<S>& v) {
        s.init(v.s.coeff(), v.s.size());
      }

      //! Initialize coefficients from an Vector with different storage
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      void init(const Vector<S>& v) volatile {
        s.init(v.s.coeff(), v.s.size());
      }

      //! Load coefficients to an array of values
      KOKKOS_INLINE_FUNCTION
      void load(value_type* v) { s.load(v); }

      //! Load coefficients to an array of values
      KOKKOS_INLINE_FUNCTION
      void load(value_type* v) volatile { s.load(v); }

      //! Load coefficients into an Vector with different storage
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      void load(Vector<S>& v) { s.load(v.s.coeff()); }

      //! Load coefficients into an Vector with different storage
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      void load(Vector<S>& v) volatile { s.load(v.s.coeff()); }

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

      //! Reset size
      /*!
       * Coefficients are preserved.
       */
      KOKKOS_INLINE_FUNCTION
      void reset(ordinal_type sz_new) volatile {
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
      void copyForWrite() volatile {  }

      //! Returns whether two MP objects have the same values
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

      //! Returns whether two MP objects have the same values
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      bool isEqualTo(const Expr<S>& xx) const volatile {
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

      //! Assignment from initializer_list
      /*!
       * No KOKKOS_INLINE_FUNCTION as it is not callable from the device
       */
      Vector& operator=(std::initializer_list<value_type> l) {
        const ordinal_type lsz = l.size();
        if (lsz != s.size())
          s.resize(lsz);
        s.init(l.begin(), lsz);
        return *this;
      }

      //! Assignment from initializer_list
      /*!
       * No KOKKOS_INLINE_FUNCTION as it is not callable from the device
       */
      /*volatile*/ Vector&
      operator=(std::initializer_list<value_type> l) volatile {
        const ordinal_type lsz = l.size();
        if (lsz != s.size())
          s.resize(lsz);
        s.init(l.begin(), lsz);
        return const_cast<Vector&>(*this);
      }

      //! Assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator=(const value_type& x) {
        s.init(x);
        return *this;
      }

      //! Assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator=(const value_type& x) volatile {
        s.init(x);
        return const_cast<Vector&>(*this);
      }

      //! Assignment with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator=(const Vector& x) {
        if (this != &x) {
          s = x.s;

          // For DyamicStorage as a view (is_owned=false), we need to set
          // the trailing entries when assigning a constant vector (because
          // the copy constructor in this case doesn't reset the size of this)
          //
          // Note:  supporting this technically makes the Vector non-POD, even
          // with a static storage type where this branch will get optimized
          // out.  We would have to remove DynamicStorage-as-a view as an option
          // or partial specialize on StaticFixedStorage to fix this.  However
          // the volatile operator=() and copy constructor overloads make
          // Vector non-POD anyway.
          if (s.size() > x.s.size())
            for (ordinal_type i=x.s.size(); i<s.size(); i++)
              s[i] = s[0];
        }

        return *this;
      }

      //! Assignment with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator=(const volatile Vector& x) {
        if (this != &x) {
          s = x.s;

          // For DyamicStorage as a view (is_owned=false), we need to set
          // the trailing entries when assigning a constant vector (because
          // the copy constructor in this case doesn't reset the size of this)
          //
          // Note:  supporting this technically makes the Vector non-POD, even
          // with a static storage type where this branch will get optimized
          // out.  We would have to remove DynamicStorage-as-a view as an option
          // or partial specialize on StaticFixedStorage to fix this.  However
          // the volatile operator=() and copy constructor overloads make
          // Vector non-POD anyway.
          if (s.size() > x.s.size())
            for (ordinal_type i=x.s.size(); i<s.size(); i++)
              s[i] = s[0];
        }

        return *this;
      }

      //! Assignment with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator=(const Vector& x) volatile {
        if (this != &x) {
          s = x.s;

          // For DyamicStorage as a view (is_owned=false), we need to set
          // the trailing entries when assigning a constant vector (because
          // the copy constructor in this case doesn't reset the size of this)
          //
          // Note:  supporting this technically makes the Vector non-POD, even
          // with a static storage type where this branch will get optimized
          // out.  We would have to remove DynamicStorage-as-a view as an option
          // or partial specialize on StaticFixedStorage to fix this.  However
          // the volatile operator=() and copy constructor overloads make
          // Vector non-POD anyway.
          if (s.size() > x.s.size())
            for (ordinal_type i=x.s.size(); i<s.size(); i++)
              s[i] = s[0];
        }

        return const_cast<Vector&>(*this);
      }

      //! Assignment with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator=(const volatile Vector& x) volatile {
        if (this != &x) {
          s = x.s;

          // For DyamicStorage as a view (is_owned=false), we need to set
          // the trailing entries when assigning a constant vector (because
          // the copy constructor in this case doesn't reset the size of this)
          //
          // Note:  supporting this technically makes the Vector non-POD, even
          // with a static storage type where this branch will get optimized
          // out.  We would have to remove DynamicStorage-as-a view as an option
          // or partial specialize on StaticFixedStorage to fix this.  However
          // the volatile operator=() and copy constructor overloads make
          // Vector non-POD anyway.
          if (s.size() > x.s.size())
            for (ordinal_type i=x.s.size(); i<s.size(); i++)
              s[i] = s[0];
        }

        return const_cast<Vector&>(*this);
      }

      //! Assignment with any expression right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Vector& operator=(const Expr<S>& xx) {
        typedef typename Expr<S>::derived_type expr_type;
        const expr_type& x = xx.derived();

        this->reset(x.size());

#ifdef STOKHOS_DEBUG
        if (s.size() != x.size())
          Kokkos::Impl::raise_error("Vector::operator=():  Mismatched sizes");
#endif

        if (x.hasFastAccess(s.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] = x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] = x.coeff(i);
        }
        return *this;
      }

      //! Assignment with any expression right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator=(const Expr<S>& xx) volatile {
        typedef typename Expr<S>::derived_type expr_type;
        const expr_type& x = xx.derived();

        this->reset(x.size());

#ifdef STOKHOS_DEBUG
        if (s.size() != x.size())
          Kokkos::Impl::raise_error("Vector::operator=():  Mismatched sizes");
#endif

        if (x.hasFastAccess(s.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] = x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] = x.coeff(i);
        }
        return const_cast<Vector&>(*this);
      }

      //! Assignment operator only valid for view storage
      template< typename S >
      KOKKOS_INLINE_FUNCTION
      typename std::enable_if<( ! std::is_same<S,void>::value &&
                                         Stokhos::is_ViewStorage<Storage>::value
                                       ), Vector >
        ::type const & operator = ( const Expr<S> & xx ) const
      {
        const typename Expr<S>::derived_type & x = xx.derived();

#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for ( ordinal_type i = 0 ; i < s.size() ; ++i ) { s[i] = x.coeff(i); }

        return *this ;
      }

      //! Assignment operator only valid for view storage
      template< typename S >
      KOKKOS_INLINE_FUNCTION
      volatile
      typename std::enable_if<( ! std::is_same<S,void>::value &&
                                         Stokhos::is_ViewStorage<Storage>::value
                                       ), Vector >
        ::type const & operator = ( const Expr<S> & xx ) const volatile
      {
        const typename Expr<S>::derived_type & x = xx.derived();

#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for ( ordinal_type i = 0 ; i < s.size() ; ++i ) { s[i] = x.coeff(i); }

        return *this ;
      }

      //@}

      /*!
       * Accessor methods
       */

      //! Returns storage object
      KOKKOS_INLINE_FUNCTION
      const volatile storage_type& storage() const volatile { return s; }

      //! Returns storage object
      KOKKOS_INLINE_FUNCTION
      const storage_type& storage() const { return s; }

      //! Returns storage object
      KOKKOS_INLINE_FUNCTION
      volatile storage_type& storage() volatile { return s; }

      //! Returns storage object
      KOKKOS_INLINE_FUNCTION
      storage_type& storage() { return s; }

      /*!
       * @name Value accessor methods
       */
      //@{

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      const_volatile_reference val() const volatile { return s[0]; }

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      const_reference val() const { return s[0]; }

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      volatile_reference val() volatile { return s[0]; }

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      reference val() { return s[0]; }

      //@}

      /*!
       * @name Coefficient accessor methods
       */
      //@{

      //! Returns size of vector
      KOKKOS_INLINE_FUNCTION
      ordinal_type size() const { return s.size();}

      //! Returns size of vector
      KOKKOS_INLINE_FUNCTION
      ordinal_type size() const volatile { return s.size();}

      //! Returns true if vector has size >= sz
      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess(ordinal_type sz) const { return s.size()>=sz;}

      //! Returns true if vector has size >= sz
      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess(ordinal_type sz) const volatile { return s.size()>=sz;}

      //! Returns coefficient array
      KOKKOS_INLINE_FUNCTION
      const_pointer coeff() const { return s.coeff();}

      //! Returns coefficient array
      KOKKOS_INLINE_FUNCTION
      const_volatile_pointer coeff() const volatile { return s.coeff();}

      //! Returns coefficient array
      KOKKOS_INLINE_FUNCTION
      volatile_pointer coeff() volatile { return s.coeff();}

      //! Returns coefficient array
      KOKKOS_INLINE_FUNCTION
      pointer coeff() { return s.coeff();}

      //! Returns term \c i with bounds checking
      KOKKOS_INLINE_FUNCTION
      value_type coeff(ordinal_type i) const volatile {
        return i<s.size() ? s[i] : s[0]; }

      //! Returns term \c i with bounds checking
      KOKKOS_INLINE_FUNCTION
      value_type coeff(ordinal_type i) const {
        return i<s.size() ? s[i] : s[0]; }

      //! Returns term \c i with bounds checking
      KOKKOS_INLINE_FUNCTION
      value_type coeff(ordinal_type i) volatile {
        return i<s.size() ? s[i] : s[0]; }

      //! Returns term \c i with bounds checking
      KOKKOS_INLINE_FUNCTION
      value_type coeff(ordinal_type i) {
        return i<s.size() ? s[i] : s[0]; }

      //! Returns term \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      const_volatile_reference fastAccessCoeff(ordinal_type i) const volatile {
        return s[i];}

      //! Returns term \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      const_reference fastAccessCoeff(ordinal_type i) const {
        return s[i];}

      //! Returns term \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      volatile_reference fastAccessCoeff(ordinal_type i) volatile {
        return s[i];}

      //! Returns term \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      reference fastAccessCoeff(ordinal_type i) {
        return s[i];}

      //! Returns term \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      const_volatile_reference operator[](ordinal_type i) const volatile {
        return s[i];}

      //! Returns term \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      const_reference operator[](ordinal_type i) const {
        return s[i];}

      //! Returns term \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      volatile_reference operator[](ordinal_type i) volatile {
        return s[i];}

      //! Returns term \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      reference operator[](ordinal_type i) {
        return s[i];}

      template <int i>
      KOKKOS_INLINE_FUNCTION
      value_type getCoeff() const volatile {
        return s.template getCoeff<i>(); }

      template <int i>
      KOKKOS_INLINE_FUNCTION
      value_type getCoeff() const {
        return s.template getCoeff<i>(); }

      template <int i>
      KOKKOS_INLINE_FUNCTION
      volatile_reference getCoeff() volatile {
        return s.template getCoeff<i>(); }

      template <int i>
      KOKKOS_INLINE_FUNCTION
      reference getCoeff() {
        return s.template getCoeff<i>(); }

      //! Return iterator to first element of coefficient array
      KOKKOS_INLINE_FUNCTION
      pointer begin() { return s.coeff(); }

      //! Return iterator to first element of coefficient array
      KOKKOS_INLINE_FUNCTION
      const_pointer begin() const { return s.coeff(); }

      //! Return iterator to first element of coefficient array
      KOKKOS_INLINE_FUNCTION
      volatile_pointer begin() volatile { return s.coeff(); }

      //! Return iterator to first element of coefficient array
      KOKKOS_INLINE_FUNCTION
      const_volatile_pointer begin() const volatile { return s.coeff(); }

      //! Return iterator to first element of coefficient array
      KOKKOS_INLINE_FUNCTION
      const_pointer cbegin() const { return s.coeff(); }

      //! Return iterator to first element of coefficient array
      KOKKOS_INLINE_FUNCTION
      const_volatile_pointer cbegin() const volatile { return s.coeff(); }

      //! Return iterator following last element of coefficient array
      KOKKOS_INLINE_FUNCTION
      pointer end() { return s.coeff() + s.size(); }

      //! Return iterator following last element of coefficient array
      KOKKOS_INLINE_FUNCTION
      const_pointer end() const { return s.coeff() + s.size(); }

      //! Return iterator following last element of coefficient array
      KOKKOS_INLINE_FUNCTION
      volatile_pointer end() volatile { return s.coeff() + s.size(); }

      //! Return iterator following last element of coefficient array
      KOKKOS_INLINE_FUNCTION
      const_volatile_pointer end() const volatile { return s.coeff() + s.size(); }

      //! Return iterator following last element of coefficient array
      KOKKOS_INLINE_FUNCTION
      const_pointer cend() const { return s.coeff()+ s.size(); }

      //! Return iterator following last element of coefficient array
      KOKKOS_INLINE_FUNCTION
      const_volatile_pointer cend() const volatile { return s.coeff()+ s.size(); }

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Addition-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator += (const value_type& x) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (ordinal_type i=0; i<s.size(); i++)
          s[i] += x;
        return *this;
      }

      //! Addition-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator += (const volatile value_type& x) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (ordinal_type i=0; i<s.size(); i++)
          s[i] += x;
        return *this;
      }

      //! Addition-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator += (const value_type& x) volatile {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (ordinal_type i=0; i<s.size(); i++)
          s[i] += x;
        return const_cast<Vector&>(*this);
      }

      //! Addition-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator += (const volatile value_type& x) volatile {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (ordinal_type i=0; i<s.size(); i++)
          s[i] += x;
        return const_cast<Vector&>(*this);
      }

      //! Subtraction-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator -= (const value_type& x) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (ordinal_type i=0; i<s.size(); i++)
          s[i] -= x;
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator -= (const volatile value_type& x) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (ordinal_type i=0; i<s.size(); i++)
          s[i] -= x;
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator -= (const value_type& x) volatile {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (ordinal_type i=0; i<s.size(); i++)
          s[i] -= x;
        return const_cast<Vector&>(*this);
      }

      //! Subtraction-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator -= (const volatile value_type& x) volatile {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (ordinal_type i=0; i<s.size(); i++)
          s[i] -= x;
        return const_cast<Vector&>(*this);
      }

      //! Multiplication-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator *= (const value_type& x) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (ordinal_type i=0; i<s.size(); i++)
          s[i] *= x;
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator *= (const volatile value_type& x) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (ordinal_type i=0; i<s.size(); i++)
          s[i] *= x;
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator *= (const value_type& x) volatile {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (ordinal_type i=0; i<s.size(); i++)
          s[i] *= x;
        return const_cast<Vector&>(*this);
      }

      //! Multiplication-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator *= (const volatile value_type& x) volatile {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (ordinal_type i=0; i<s.size(); i++)
          s[i] *= x;
        return const_cast<Vector&>(*this);
      }

      //! Division-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator /= (const value_type& x) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (ordinal_type i=0; i<s.size(); i++)
          s[i] /= x;
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator /= (const volatile value_type& x) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (ordinal_type i=0; i<s.size(); i++)
          s[i] /= x;
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator /= (const value_type& x) volatile {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (ordinal_type i=0; i<s.size(); i++)
          s[i] /= x;
        return const_cast<Vector&>(*this);
      }

      //! Division-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator /= (const volatile value_type& x) volatile {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (ordinal_type i=0; i<s.size(); i++)
          s[i] /= x;
        return const_cast<Vector&>(*this);
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Vector& operator += (const Expr<S>& xx) {
        //*this = *this + x;
        typedef typename Expr<S>::derived_type expr_type;
        const expr_type& x = xx.derived();

        if (x.size() > s.size())
          this->reset(x.size());

#ifdef STOKHOS_DEBUG
        if (s.size() < x.size())
          Kokkos::Impl::raise_error("Vector::operator+=():  Mismatched sizes");
#endif

        if (x.hasFastAccess(s.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] += x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] += x.coeff(i);
        }
        return *this;
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Vector& operator += (const volatile Expr<S>& xx) {
        //*this = *this + x;
        typedef typename Expr<S>::derived_type expr_type;
        const volatile expr_type& x = xx.derived();

        if (x.size() > s.size())
          this->reset(x.size());

#ifdef STOKHOS_DEBUG
        if (s.size() < x.size())
          Kokkos::Impl::raise_error("Vector::operator+=():  Mismatched sizes");
#endif

        if (x.hasFastAccess(s.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] += x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] += x.coeff(i);
        }
        return *this;
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator += (const Expr<S>& xx) volatile {
        //*this = *this + x;
        typedef typename Expr<S>::derived_type expr_type;
        const expr_type& x = xx.derived();

        if (x.size() > s.size())
          this->reset(x.size());

#ifdef STOKHOS_DEBUG
        if (s.size() < x.size())
          Kokkos::Impl::raise_error("Vector::operator+=():  Mismatched sizes");
#endif

        if (x.hasFastAccess(s.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] += x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] += x.coeff(i);
        }
        return const_cast<Vector&>(*this);
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator += (const volatile Expr<S>& xx) volatile {
        //*this = *this + x;
        typedef typename Expr<S>::derived_type expr_type;
        const volatile expr_type& x = xx.derived();

        if (x.size() > s.size())
          this->reset(x.size());

#ifdef STOKHOS_DEBUG
        if (s.size() < x.size())
          Kokkos::Impl::raise_error("Vector::operator+=():  Mismatched sizes");
#endif

        if (x.hasFastAccess(s.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] += x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] += x.coeff(i);
        }
        return const_cast<Vector&>(*this);
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Vector& operator -= (const Expr<S>& xx) {
        //*this = *this - x;
        typedef typename Expr<S>::derived_type expr_type;
        const expr_type& x = xx.derived();

        if (x.size() > s.size())
          this->reset(x.size());

#ifdef STOKHOS_DEBUG
        if (s.size() < x.size())
          Kokkos::Impl::raise_error("Vector::operator-=():  Mismatched sizes");
#endif

        if (x.hasFastAccess(s.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] -= x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] -= x.coeff(i);
        }
        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Vector& operator -= (const volatile Expr<S>& xx) {
        //*this = *this - x;
        typedef typename Expr<S>::derived_type expr_type;
        const volatile expr_type& x = xx.derived();

        if (x.size() > s.size())
          this->reset(x.size());

#ifdef STOKHOS_DEBUG
        if (s.size() < x.size())
          Kokkos::Impl::raise_error("Vector::operator-=():  Mismatched sizes");
#endif

        if (x.hasFastAccess(s.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] -= x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] -= x.coeff(i);
        }
        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator -= (const Expr<S>& xx) volatile {
        //*this = *this - x;
        typedef typename Expr<S>::derived_type expr_type;
        const expr_type& x = xx.derived();

        if (x.size() > s.size())
          this->reset(x.size());

#ifdef STOKHOS_DEBUG
        if (s.size() < x.size())
          Kokkos::Impl::raise_error("Vector::operator-=():  Mismatched sizes");
#endif

        if (x.hasFastAccess(s.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] -= x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] -= x.coeff(i);
        }
        return const_cast<Vector&>(*this);
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator -= (const volatile Expr<S>& xx) volatile {
        //*this = *this - x;
        typedef typename Expr<S>::derived_type expr_type;
        const volatile expr_type& x = xx.derived();

        if (x.size() > s.size())
          this->reset(x.size());

#ifdef STOKHOS_DEBUG
        if (s.size() < x.size())
          Kokkos::Impl::raise_error("Vector::operator-=():  Mismatched sizes");
#endif

        if (x.hasFastAccess(s.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] -= x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] -= x.coeff(i);
        }
        return const_cast<Vector&>(*this);
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Vector& operator *= (const Expr<S>& xx) {
        //*this = *this * x;
        typedef typename Expr<S>::derived_type expr_type;
        const expr_type& x = xx.derived();

        if (x.size() > s.size())
          this->reset(x.size());

#ifdef STOKHOS_DEBUG
        if (s.size() < x.size())
          Kokkos::Impl::raise_error("Vector::operator*=():  Mismatched sizes");
#endif

        if (x.hasFastAccess(s.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] *= x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] *= x.coeff(i);
        }
        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Vector& operator *= (const volatile Expr<S>& xx) {
        //*this = *this * x;
        typedef typename Expr<S>::derived_type expr_type;
        const volatile expr_type& x = xx.derived();

        if (x.size() > s.size())
          this->reset(x.size());

#ifdef STOKHOS_DEBUG
        if (s.size() < x.size())
          Kokkos::Impl::raise_error("Vector::operator*=():  Mismatched sizes");
#endif

        if (x.hasFastAccess(s.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] *= x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] *= x.coeff(i);
        }
        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator *= (const Expr<S>& xx) volatile {
        //*this = *this * x;
        typedef typename Expr<S>::derived_type expr_type;
        const expr_type& x = xx.derived();

        if (x.size() > s.size())
          this->reset(x.size());

#ifdef STOKHOS_DEBUG
        if (s.size() < x.size())
          Kokkos::Impl::raise_error("Vector::operator*=():  Mismatched sizes");
#endif

        if (x.hasFastAccess(s.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] *= x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] *= x.coeff(i);
        }
        return const_cast<Vector&>(*this);
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator *= (const volatile Expr<S>& xx) volatile {
        //*this = *this * x;
        typedef typename Expr<S>::derived_type expr_type;
        const volatile expr_type& x = xx.derived();

        if (x.size() > s.size())
          this->reset(x.size());

#ifdef STOKHOS_DEBUG
        if (s.size() < x.size())
          Kokkos::Impl::raise_error("Vector::operator*=():  Mismatched sizes");
#endif

        if (x.hasFastAccess(s.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] *= x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] *= x.coeff(i);
        }
        return const_cast<Vector&>(*this);
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Vector& operator /= (const Expr<S>& xx) {
        //*this = *this / x;
        typedef typename Expr<S>::derived_type expr_type;
        const expr_type& x = xx.derived();

        if (x.size() > s.size())
          this->reset(x.size());

#ifdef STOKHOS_DEBUG
        if (s.size() < x.size())
          Kokkos::Impl::raise_error("Vector::operator/=():  Mismatched sizes");
#endif

        if (x.hasFastAccess(s.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] /= x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] /= x.coeff(i);
        }
        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Vector& operator /= (const volatile Expr<S>& xx) {
        //*this = *this / x;
        typedef typename Expr<S>::derived_type expr_type;
        const volatile expr_type& x = xx.derived();

        if (x.size() > s.size())
          this->reset(x.size());

#ifdef STOKHOS_DEBUG
        if (s.size() < x.size())
          Kokkos::Impl::raise_error("Vector::operator/=():  Mismatched sizes");
#endif

        if (x.hasFastAccess(s.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] /= x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] /= x.coeff(i);
        }
        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator /= (const Expr<S>& xx) volatile {
        //*this = *this / x;
        typedef typename Expr<S>::derived_type expr_type;
        const expr_type& x = xx.derived();

        if (x.size() > s.size())
          this->reset(x.size());

#ifdef STOKHOS_DEBUG
        if (s.size() < x.size())
          Kokkos::Impl::raise_error("Vector::operator/=():  Mismatched sizes");
#endif

        if (x.hasFastAccess(s.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] /= x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] /= x.coeff(i);
        }
        return const_cast<Vector&>(*this);
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator /= (const volatile Expr<S>& xx) volatile {
        //*this = *this / x;
        typedef typename Expr<S>::derived_type expr_type;
        const volatile expr_type& x = xx.derived();

        if (x.size() > s.size())
          this->reset(x.size());

#ifdef STOKHOS_DEBUG
        if (s.size() < x.size())
          Kokkos::Impl::raise_error("Vector::operator/=():  Mismatched sizes");
#endif

        if (x.hasFastAccess(s.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] /= x.fastAccessCoeff(i);
        }
        else {
          for (ordinal_type i=0; i<s.size(); i++)
            s[i] /= x.coeff(i);
        }
        return const_cast<Vector&>(*this);
      }

      //! Prefix ++
      KOKKOS_INLINE_FUNCTION
      Vector& operator++() {
        for (ordinal_type i=0; i<s.size(); i++)
          ++(s[i]);
        return *this;
      }

      //! Prefix ++
      KOKKOS_INLINE_FUNCTION
      volatile Vector& operator++() volatile {
        for (ordinal_type i=0; i<s.size(); i++)
          ++(s[i]);
        return *this;
      }

      //! Postfix ++
      KOKKOS_INLINE_FUNCTION
      Vector operator++(int) {
        Vector tmp(*this);
        ++(*this);
        return tmp;
      }

      //! Postfix ++
      KOKKOS_INLINE_FUNCTION
      Vector operator++(int) volatile {
        Vector tmp(*this);
        ++(*this);
        return tmp;
      }

      //! Prefix --
      KOKKOS_INLINE_FUNCTION
      Vector& operator--() {
        for (ordinal_type i=0; i<s.size(); i++)
          --(s[i]);
        return *this;
      }

      //! Prefix --
      KOKKOS_INLINE_FUNCTION
      volatile Vector& operator--() volatile {
        for (ordinal_type i=0; i<s.size(); i++)
          --(s[i]);
        return *this;
      }

      //! Postfix --
      KOKKOS_INLINE_FUNCTION
      Vector operator--(int) {
        Vector tmp(*this);
        --(*this);
        return tmp;
      }

      //! Postfix --
      KOKKOS_INLINE_FUNCTION
      Vector operator--(int) volatile {
        Vector tmp(*this);
        --(*this);
        return tmp;
      }

      //@}

      KOKKOS_INLINE_FUNCTION
      std::string name() const volatile { return "x"; }

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
  }
}

#if STOKHOS_USE_MP_VECTOR_SFS_SPEC
#include "Sacado_MP_Vector_SFS.hpp"
#endif

namespace Sacado{
  namespace MP {

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

    //! Type for storing nodes in expression graph
    /*!
     * Specialization for leaf-nodes, which can be stored by reference
     * since they are an argument to the expression.
     */
    template <typename S> struct const_expr_ref< volatile Vector<S> > {
      typedef const volatile Vector<S>& type;
    };

    //! Traits class for removing volatile from type
    template <typename T> struct remove_volatile {
      typedef T type;
    };
    template <typename T> struct remove_volatile<volatile T> {
      typedef T type;
    };

    //! Traits class for adding volatile to type
    template <typename T> struct add_volatile {
      typedef volatile T type;
    };
    template <typename T> struct add_volatile<volatile T> {
      typedef volatile T type;
    };

    template <typename Storage>
    std::ostream&
    operator << (std::ostream& os, const Vector<Storage>& a)
    {
      typedef typename Vector<Storage>::ordinal_type ordinal_type;

      os << "[ ";

      for (ordinal_type i=0; i<a.size(); i++) {
        os << a.coeff(i) << " ";
      }

      os << "]";
      return os;
    }

    template <typename Storage>
    std::ostream&
    operator << (std::ostream& os, const volatile Vector<Storage>& a)
    {
      typedef typename Vector<Storage>::ordinal_type ordinal_type;

      os << "[ ";

      for (ordinal_type i=0; i<a.size(); i++) {
        os << a.coeff(i) << " ";
      }

      os << "]";
      return os;
    }

    template <typename Storage>
    std::istream&
    operator >> (std::istream& is, Vector<Storage>& a)
    {
      typedef typename Vector<Storage>::ordinal_type ordinal_type;
      typedef typename Vector<Storage>::value_type value_type;

      //
      // Need to check all of this for errors, end-of-line, etc...
      //

      char b = 0;
      if (Storage::is_static) {
        is >> b; // "["
        for (ordinal_type i=0; i<a.size(); i++) {
          is >> a.fastAccessCoeff(i);
        }
        is >> b; // "]";
      }
      else {
        std::vector<value_type> c;
        value_type v;
        is >> b; // "["
        while (is >> b && b != ']') {
          is.putback(b);
          is >> v;
          c.push_back(v);
        }
        ordinal_type n = c.size();
        a.reset(n);
        for (ordinal_type i=0; i<n; ++i)
          a.fastAccessCoeff(i) = c[i];
      }

      return is;
    }

    //------------------------------------------------------------------------
    //------------------------------------------------------------------------

    /**\brief  Define a partition of a View of Sacado::MP::Vector type */
    template <unsigned Size = 0>
    struct VectorPartition {
      static const unsigned PartitionSize = Size;
      unsigned begin ;
      unsigned end ;

      template< typename iType0 , typename iType1 >
      KOKKOS_INLINE_FUNCTION
      VectorPartition( const iType0 & i0 , const iType1 & i1 ) :
        begin(i0), end(i1) {
      }
    };

    template <typename T>
    struct is_vector_partition {
      static const bool value = false;
    };

    template <unsigned Size>
    struct is_vector_partition< VectorPartition<Size> > {
      static const bool value = true;
    };

  } // namespace MP

  template <typename T>
  struct IsExpr< MP::Expr<T> > {
    static const bool value = true;
  };

  template <typename T>
  struct BaseExprType< MP::Expr<T> > {
    typedef typename MP::Expr<T>::derived_type derived_type;
    typedef typename derived_type::base_expr_type type;
  };

  template <typename S>
  struct IsExpr< MP::Vector<S> > {
    static const bool value = true;
  };

  template <typename S>
  struct BaseExprType< MP::Vector<S> > {
    typedef MP::Vector<S> type;
  };

  //! Trait class to determine if a scalar type is a Vector
  template <typename T> struct is_mp_vector {
    static const bool value = false;
  };
  template <typename S> struct is_mp_vector< MP::Vector<S> > {
    static const bool value = true;
  };
  template <typename T> struct is_mp_vector< const T > {
    static const bool value = is_mp_vector<T>::value;
  };
  template <typename T> struct is_mp_vector< T* > {
    static const bool value = is_mp_vector<T>::value;
  };
  template <typename T> struct is_mp_vector< T[] > {
    static const bool value = is_mp_vector<T>::value;
  };
  template <typename T, unsigned N> struct is_mp_vector< T[N] > {
    static const bool value = is_mp_vector<T>::value;
  };

  // Utility function to see if a MP::Vector is really a constant
  template <typename Storage>
  bool is_constant(const Sacado::MP::Vector<Storage>& x)
  {
    typedef typename Storage::ordinal_type ordinal_type;
    typedef typename Storage::value_type value_type;

    // All size-1 vectors are constants
    const ordinal_type sz = x.size();
    if (sz == 1) return true;

    // Maybe use a tolerance????
    const value_type val = x.fastAccessCoeff(0);
    for (ordinal_type i=1; i<sz; ++i)
      if (x.fastAccessCoeff(i) != val) return false;

    return true;
  }

} // namespace Sacado

#include "Sacado_MP_Vector_ops.hpp"
#include "Stokhos_MP_Vector_MaskTraits.hpp"

#if STOKHOS_ALIGN_MEMORY

#include <memory>

namespace std {

template <typename Storage>
class allocator< Sacado::MP::Vector< Storage > >
  : public Stokhos::aligned_allocator< Sacado::MP::Vector< Storage > > {
public:
  typedef Sacado::MP::Vector<Storage>    T;
  typedef Stokhos::aligned_allocator<T>  Base;
  typedef typename Base::value_type      value_type;
  typedef typename Base::pointer         pointer;
  typedef typename Base::const_pointer   const_pointer;
  typedef typename Base::reference       reference;
  typedef typename Base::const_reference const_reference;
  typedef typename Base::size_type       size_type;
  typedef typename Base::difference_type difference_type;

  template <class U> struct rebind { typedef allocator<U> other; };
  allocator() {}
  template <class U> allocator(const allocator<U>&) {}
};

template <typename Storage>
class allocator< const Sacado::MP::Vector< Storage > >
  : public Stokhos::aligned_allocator< const Sacado::MP::Vector< Storage > > {
public:
  typedef Sacado::MP::Vector<Storage>    T;
  typedef Stokhos::aligned_allocator<const T> Base;
  typedef typename Base::value_type      value_type;
  typedef typename Base::pointer         pointer;
  typedef typename Base::const_pointer   const_pointer;
  typedef typename Base::reference       reference;
  typedef typename Base::const_reference const_reference;
  typedef typename Base::size_type       size_type;
  typedef typename Base::difference_type difference_type;

  template <class U> struct rebind { typedef allocator<U> other; };
  allocator() {}
  template <class U> allocator(const allocator<U>&) {}
};

}

#endif

#include "Kokkos_NumericTraits.hpp"

namespace Kokkos {

template <typename Storage>
struct reduction_identity< Sacado::MP::Vector<Storage> > {
  typedef Sacado::MP::Vector<Storage> Vector;
  typedef typename Storage::value_type scalar;
  typedef reduction_identity<scalar> RIS;
  KOKKOS_FORCEINLINE_FUNCTION constexpr static Vector sum()  {
    return Vector(RIS::sum());
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static Vector prod() {
    return Vector(RIS::prod());
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static Vector max()  {
    return Vector(RIS::max());
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static Vector min()  {
    return Vector(RIS::min());
  }
};

namespace Impl {
  template <typename Storage>
  struct promote<Sacado::MP::Vector<Storage>,false> {
    using type = typename Sacado::MP::Vector<Storage>;
  };
}

}

#endif // HAVE_STOKHOS_SACADO

#endif // SACADO_MP_VECTOR_HPP
