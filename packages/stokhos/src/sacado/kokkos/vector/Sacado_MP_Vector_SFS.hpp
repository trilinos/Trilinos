// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SACADO_MP_VECTOR_SFS_HPP
#define SACADO_MP_VECTOR_SFS_HPP

#include "Stokhos_ConfigDefs.h"

#ifdef HAVE_STOKHOS_SACADO

//#include "Sacado_MP_Vector.hpp"
#include "Stokhos_StaticFixedStorage.hpp"

namespace Sacado {

  namespace MP {

    //! Vectorized evaluation class
    template <typename ordinal_t, typename value_t, int Num, typename device_t>
    class Vector< Stokhos::StaticFixedStorage<ordinal_t,value_t,Num,device_t> >{
    public:

      //! Typename of storage class
      typedef Stokhos::StaticFixedStorage<ordinal_t,value_t,Num,device_t> Storage;
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
       * Sets size to 1 and first coefficient to 0 (represents a constant).
       */
      KOKKOS_INLINE_FUNCTION
      Vector() : s(1) {}

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
      KOKKOS_INLINE_FUNCTION
      Vector(const Vector& x) : s(x.s) {}

      //! Copy constructor
      KOKKOS_INLINE_FUNCTION
      Vector(const volatile Vector& x) : s(x.s) {}

      //! Intialize from initializer_list
      /*!
       * No KOKKOS_INLINE_FUNCTION as it is not callable from the device
       */
      Vector(std::initializer_list<value_type> l) : s(l.size(), l.begin()) {}

      //! Destructor
      KOKKOS_INLINE_FUNCTION
      ~Vector() {}

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
      void reset(ordinal_type sz_new) {}

      //! Reset size
      /*!
       * Coefficients are preserved.
       */
      KOKKOS_INLINE_FUNCTION
      void reset(ordinal_type sz_new) volatile {}

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
      KOKKOS_INLINE_FUNCTION
      bool isEqualTo(const Vector& x) const {
        typedef IsEqual<value_type> IE;
        bool eq = true;
        for (ordinal_type i=0; i<this->size(); i++)
          eq = eq && IE::eval(x.coeff(i), this->coeff(i));
        return eq;
      }

      //! Returns whether two MP objects have the same values
      KOKKOS_INLINE_FUNCTION
      bool isEqualTo(const Vector& x) const volatile {
        typedef IsEqual<value_type> IE;
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
        }

        return *this;
      }

      //! Assignment with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator=(const volatile Vector& x) {
        if (this != &x) {
          s = x.s;
        }

        return *this;
      }

      //! Assignment with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator=(const Vector& x) volatile {
        if (this != &x) {
          s = x.s;
        }

        return const_cast<Vector&>(*this);
      }

      //! Assignment with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator=(const volatile Vector& x) volatile {
        if (this != &x) {
          s = x.s;
        }

        return const_cast<Vector&>(*this);
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
      static ordinal_type size() { return storage_type::size();}

      //! Returns true if vector has size >= sz
      KOKKOS_INLINE_FUNCTION
      static bool hasFastAccess(ordinal_type sz) { return true; }

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
      value_type coeff(ordinal_type i) const volatile { return s[i]; }

      //! Returns term \c i with bounds checking
      KOKKOS_INLINE_FUNCTION
      value_type coeff(ordinal_type i) const { return s[i]; }

      //! Returns term \c i with bounds checking
      KOKKOS_INLINE_FUNCTION
      value_type coeff(ordinal_type i) volatile { return s[i]; }

      //! Returns term \c i with bounds checking
      KOKKOS_INLINE_FUNCTION
      value_type coeff(ordinal_type i) { return s[i]; }

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

      //! Addition-assignment operator with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator += (const Vector& x) {
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
        return *this;
      }

      //! Addition-assignment operator with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator += (const volatile Vector& x) {
        *this = *this + x;
        return *this;
      }

      //! Addition-assignment operator with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator += (const Vector& x) volatile {
        *this = *this + x;
        return const_cast<Vector&>(*this);
      }

      //! Addition-assignment operator with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator += (const volatile Vector& x) volatile {
        *this = *this + x;
        return const_cast<Vector&>(*this);
      }

      //! Subtraction-assignment operator with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator -= (const Vector& x) {
        *this = *this - x;
        return *this;
      }

      //! Subtraction-assignment operator with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator -= (const volatile Vector& x) {
        *this = *this - x;
        return *this;
      }

      //! Subtraction-assignment operator with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator -= (const Vector& x) volatile {
        *this = *this - x;
        return const_cast<Vector&>(*this);
      }

      //! Subtraction-assignment operator with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator -= (const volatile Vector& x) volatile {
        *this = *this - x;
        return const_cast<Vector&>(*this);
      }

      //! Multiplication-assignment operator with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator *= (const Vector& x) {
        *this = *this * x;
        return *this;
      }

      //! Multiplication-assignment operator with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator *= (const volatile Vector& x) {
        *this = *this * x;
        return *this;
      }

      //! Multiplication-assignment operator with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator *= (const Vector& x) volatile {
        *this = *this * x;
        return const_cast<Vector&>(*this);
      }

      //! Multiplication-assignment operator with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator *= (const volatile Vector& x) volatile {
        *this = *this * x;
        return const_cast<Vector&>(*this);
      }

      //! Division-assignment operator with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator /= (const Vector& x) {
        *this = *this / x;
        return *this;
      }

      //! Division-assignment operator with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      Vector& operator /= (const volatile Vector& x) {
        *this = *this / x;
        return *this;
      }

      //! Division-assignment operator with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator /= (const Vector& x) volatile {
        *this = *this / x;
        return const_cast<Vector&>(*this);
      }

      //! Division-assignment operator with Vector right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ Vector& operator /= (const volatile Vector& x) volatile {
        *this = *this / x;
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

    }; // class Vector

    //------------------------------------------------------------------------
    //------------------------------------------------------------------------

  } // namespace MP

} // namespace Sacado

#include "Sacado_MP_Vector_SFS_ops.hpp"

#endif // HAVE_STOKHOS_SACADO

#endif // SACADO_MP_VECTOR_SFS_HPP
