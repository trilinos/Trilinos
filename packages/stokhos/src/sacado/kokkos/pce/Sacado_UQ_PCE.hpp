// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SACADO_UQ_PCE_HPP
#define SACADO_UQ_PCE_HPP

#include "Stokhos_ConfigDefs.h"

#ifdef HAVE_STOKHOS_SACADO

#include "Kokkos_Macros.hpp"

#include "Sacado_Traits.hpp"
#include "Sacado_mpl_apply.hpp"
#include "Stokhos_Is_Constant.hpp"

#include "Stokhos_CrsProductTensor.hpp"

#include <cmath>
#include <algorithm>    // for std::min and std::max
#include <ostream>      // for std::ostream
#include <initializer_list>

namespace Sacado {

  // Forward declaration
  template <typename Storage>
  KOKKOS_INLINE_FUNCTION
  bool is_constant(const Sacado::UQ::PCE<Storage>& x);

  //! Namespace for UQ classes classes
  namespace UQ {

    //! Generalized polynomial chaos expansion class
    /*!
     * Uses a handle and a "copy-on-write" strategy for efficient copying, but
     * no expression templating.
     */
    template <typename Storage >
    class PCE {

      template <class> friend class PCE;

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

      //! Typename of scalar's (which may be different from value_type)
      typedef typename ScalarType<value_type>::type scalar_type;

      //! Basis type
      //typedef Stokhos::OrthogPolyBasis<ordinal_type,value_type> basis_type;

      //! Cijk type
      typedef Stokhos::CrsProductTensor<value_type, execution_space, Kokkos::MemoryUnmanaged> my_cijk_type;
      typedef Stokhos::CrsProductTensor<value_type, execution_space> cijk_type;

      //! Turn PCE into a meta-function class usable with mpl::apply
      template <typename S>
      struct apply {
        typedef PCE<S> type;
      };

      //! Default constructor
      /*!
       * May not intialize the coefficient array.
       */
      KOKKOS_DEFAULTED_FUNCTION 
      PCE() = default;

      //! Constructor with supplied value \c x
      /*!
       * Sets size to 1 and first coefficient to x (represents a constant).
       */
      KOKKOS_INLINE_FUNCTION
      PCE(const value_type& x) : cijk_(), s_(1, x) {}

      //! Constructor with Cijk \c cijk (General case)
      /*!
       * Creates array of correct size and initializes coeffiencts to 0.
       */
      template <typename M>
      KOKKOS_INLINE_FUNCTION
      PCE(const Stokhos::CrsProductTensor<value_type, execution_space, M>& cijkVal) :
        cijk_(cijkVal), s_(cijk_.dimension()) {}

      //! Constructor with Cijk \c cijk and specified size \c sz
      /*!
       * Creates array of size \c sz and initializes coeffiencts to 0.
       */
      template <typename M>
      KOKKOS_INLINE_FUNCTION
      PCE(const Stokhos::CrsProductTensor<value_type, execution_space, M>& cijkVal,
          ordinal_type sz) : cijk_(cijkVal), s_(sz) {}

      //! View constructor
      /*!
       * Creates PCE with pre-allocated data.  Set \c owned = true
       * if this PCE should take over management of the data.
       */
      template <typename M>
      KOKKOS_INLINE_FUNCTION
      PCE(const Stokhos::CrsProductTensor<value_type, execution_space, M>& cijkVal,
          ordinal_type sz, pointer v, bool owned) :
        cijk_(cijkVal), s_(sz,v,owned) {}

      //! Copy constructor
      KOKKOS_INLINE_FUNCTION
      PCE(const PCE& x) : cijk_(x.cijk_), s_(1,x.fastAccessCoeff(0)) {
        if (x.size() > 1 && !is_constant(x))
          s_ = x.s_;
      }

      //! Copy constructor
      KOKKOS_INLINE_FUNCTION
      PCE(const volatile PCE& x) :
        cijk_(const_cast<const my_cijk_type&>(x.cijk_)),
        s_(1,value_type(x.fastAccessCoeff(0))) {
        if (x.size() > 1 && !is_constant(x))
          s_ = x.s_;
      }

      //! Intialize from initializer_list
      /*!
       * No KOKKOS_INLINE_FUNCTION as it is not callable from the device
       */
      PCE(std::initializer_list<value_type> l) :
        cijk_(), s_(l.size(), l.begin()) {}

      //! Destructor
      KOKKOS_DEFAULTED_FUNCTION
      ~PCE() = default;

      //! Initialize coefficients to value
      KOKKOS_INLINE_FUNCTION
      void init(const_reference v) { s_.init(v); }

      //! Initialize coefficients to value
      KOKKOS_INLINE_FUNCTION
      void init(const_reference v) volatile { s_.init(v); }

      //! Initialize coefficients to an array of values
      KOKKOS_INLINE_FUNCTION
      void init(const_pointer v) { s_.init(v); }

      //! Initialize coefficients to an array of values
      KOKKOS_INLINE_FUNCTION
      void init(const_pointer v) volatile { s_.init(v); }

      //! Initialize coefficients from an PCE with different storage
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      void init(const PCE<S>& v) { s_.init(v.coeff()); }

      //! Initialize coefficients from an PCE with different storage
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      void init(const PCE<S>& v) volatile { s_.init(v.coeff()); }

      //! Load coefficients to an array of values
      KOKKOS_INLINE_FUNCTION
      void load(pointer v) { s_.load(v); }

      //! Load coefficients to an array of values
      KOKKOS_INLINE_FUNCTION
      void load(pointer v) volatile { s_.load(v); }

      //! Load coefficients into an PCE with different storage
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      void load(PCE<S>& v) { s_.load(v.coeff()); }

      //! Load coefficients into an PCE with different storage
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      void load(PCE<S>& v) volatile { s_.load(v.coeff()); }

      //! Reset expansion
      /*!
       * May change size of array.  Coefficients are preserved.
       */
      template <typename M>
      KOKKOS_INLINE_FUNCTION
      void reset(const Stokhos::CrsProductTensor<value_type, execution_space, M>& cijkVal) {
        cijk_ = cijkVal;
        s_.resize(cijk_.dimension());
      }

      //! Reset expansion
      /*!
       * May change size of array.  Coefficients are preserved.
       */
      template <typename M>
      KOKKOS_INLINE_FUNCTION
      void reset(const Stokhos::CrsProductTensor<value_type, execution_space, M>& cijkVal) volatile {
        cijk_ = cijkVal;
        s_.resize(cijk_.dimension());
      }

      //! Reset expansion and size
      /*!
       * Coefficients are preserved.
       */
      template <typename M>
      KOKKOS_INLINE_FUNCTION
      void reset(const Stokhos::CrsProductTensor<value_type, execution_space, M>& cijkVal, ordinal_type sz) {
        cijk_ = cijkVal;
        s_.resize(sz);
      }

      //! Reset expansion and size
      /*!
       * Coefficients are preserved.
       */
      template <typename M>
      KOKKOS_INLINE_FUNCTION
      void reset(const Stokhos::CrsProductTensor<value_type, execution_space, M>& cijkVal, ordinal_type sz) volatile {
        cijk_ = cijkVal;
        s_.resize(sz);
      }

      //! Prepare polynomial for writing
      /*!
       * This method prepares the polynomial for writing through coeff() and
       * fastAccessCoeff() member functions.  It ensures the handle for the
       * %PCE coefficients is not shared among any other %PCE polynomial
       * objects.  If the handle is not shared it does nothing, so there
       * is no cost in calling this method in this case.  If the handle is
       * shared and this method is not called, any changes to the coefficients
       * by coeff() or fastAccessCoeff() may change other polynomial objects.
       */
      KOKKOS_INLINE_FUNCTION
      void copyForWrite() { }

      /*
      //! Evaluate polynomial approximation at a point
      KOKKOS_INLINE_FUNCTION
      value_type evaluate(const Teuchos::Array<value_type>& point) const;

      //! Evaluate polynomial approximation at a point with given basis values
      KOKKOS_INLINE_FUNCTION
      value_type evaluate(const Teuchos::Array<value_type>& point,
                          const Teuchos::Array<value_type>& bvals) const;
      */

      //! Compute mean of expansion
      KOKKOS_INLINE_FUNCTION
      value_type mean() const { return this->fastAccessCoeff(0); }

      //! Compute standard deviation of expansion
      KOKKOS_INLINE_FUNCTION
      value_type standard_deviation() const;

      //! Compute the two-norm of expansion
      KOKKOS_INLINE_FUNCTION
      value_type two_norm() const {
        return std::sqrt(this->two_norm_squared());
      }

      //! Compute the squared two-norm of expansion
      KOKKOS_INLINE_FUNCTION
      value_type two_norm_squared() const;

      //! Compute the L2 inner product of 2 PCEs
      KOKKOS_INLINE_FUNCTION
      value_type inner_product(const PCE& x) const;

      //! Returns whether two PCE objects have the same values
      KOKKOS_INLINE_FUNCTION
      bool isEqualTo(const PCE& x) const;

      //! Returns whether two PCE objects have the same values
      KOKKOS_INLINE_FUNCTION
      bool isEqualTo(const PCE& x) const volatile;

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      PCE<Storage>& operator=(const value_type val);

      //! Assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ PCE<Storage>& operator=(const value_type val) volatile;

      //! Assignment with PCE right-hand-side
      KOKKOS_INLINE_FUNCTION
      PCE<Storage>& operator=(const PCE<Storage>& x);

      //! Assignment with PCE right-hand-side
      KOKKOS_INLINE_FUNCTION
      PCE<Storage>& operator=(const volatile PCE<Storage>& x);

      //! Assignment with PCE right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ PCE<Storage>& operator=(const PCE<Storage>& x) volatile;

      //! Assignment with PCE right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ PCE<Storage>& operator=(const volatile PCE<Storage>& x) volatile;

      //! Assignment with PCE right-hand-side and different storage
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      PCE<Storage>& operator=(const PCE<S>& x) {
        // Don't copy cijk as it may be on a different device
        const ordinal_type sz_new = x.size();
        s_.resize(sz_new);
        for (ordinal_type i=0; i<sz_new; i++)
          s_[i] = x.s_[i];

        // For DyamicStorage as a view (is_owned=false), we need to set
        // the trailing entries when assigning a constant vector (because
        // the copy constructor in this case doesn't reset the size of this)
        const ordinal_type sz = s_.size();
        if (sz > sz_new)
          for (ordinal_type i=sz_new; i<sz; i++)
            s_[i] = value_type(0);

        return *this;
      }

      //! Assignment from initializer_list
      /*!
       * No KOKKOS_INLINE_FUNCTION as it is not callable from the device
       */
      PCE& operator=(std::initializer_list<value_type> l) {
        const ordinal_type lsz = l.size();
        if (lsz != s_.size())
          s_.resize(lsz);
        s_.init(l.begin(), lsz);
        return *this;
      }

      //! Assignment from initializer_list
      /*!
       * No KOKKOS_INLINE_FUNCTION as it is not callable from the device
       */
      /*volatile*/ PCE&
      operator=(std::initializer_list<value_type> l) volatile {
        const ordinal_type lsz = l.size();
        if (lsz != s_.size())
          s_.resize(lsz);
        s_.init(l.begin(), lsz);
        return const_cast<PCE&>(*this);
      }

      //@}

      /*!
       * Accessor methods
       */
      //@{

      /*
      //! Get basis
      KOKKOS_INLINE_FUNCTION
      Teuchos::RCP<const basis_type> basis() const { return s_.basis(); }
      */

      //! Get cijk
      KOKKOS_INLINE_FUNCTION
      my_cijk_type cijk() const { return cijk_; }

      //! Get cijk
      KOKKOS_INLINE_FUNCTION
      my_cijk_type cijk() const volatile {
        return const_cast<const my_cijk_type&>(cijk_);
      }

      //@}

      /*!
       * @name Value accessor methods
       */
      //@{

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      const_volatile_reference val() const volatile { return s_[0]; }

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      const_reference val() const { return s_[0]; }

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      volatile_reference val() volatile { return s_[0]; }

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      reference val() { return s_[0]; }

      //@}

      /*!
       * @name Coefficient accessor methods
       */
      //@{

      //! Returns size of polynomial
      KOKKOS_INLINE_FUNCTION
      ordinal_type size() const { return s_.size();}

      //! Returns size of polynomial
      KOKKOS_INLINE_FUNCTION
      ordinal_type size() const volatile { return s_.size();}

      //! Returns true if polynomial has size >= sz
      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess(ordinal_type sz) const { return s_.size()>=sz;}

      //! Returns true if polynomial has size >= sz
      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess(ordinal_type sz) const volatile {
        return s_.size()>=sz;
      }

      //! Returns coefficient array
      KOKKOS_INLINE_FUNCTION
      const_pointer coeff() const { return s_.coeff();}

      //! Returns coefficient array
      KOKKOS_INLINE_FUNCTION
      const_volatile_pointer coeff() const volatile { return s_.coeff();}

      //! Returns coefficient array
      KOKKOS_INLINE_FUNCTION
      volatile_pointer coeff() volatile { return s_.coeff();}

      //! Returns coefficient array
      KOKKOS_INLINE_FUNCTION
      pointer coeff() { return s_.coeff();}

      //! Returns degree \c i term with bounds checking
      KOKKOS_INLINE_FUNCTION
      value_type coeff(ordinal_type i) const {
        value_type tmp= i<s_.size() ? s_[i]:value_type(0.); return tmp;}

      //! Returns degree \c i term with bounds checking
      KOKKOS_INLINE_FUNCTION
      value_type coeff(ordinal_type i) const volatile {
        value_type tmp= i<s_.size() ? s_[i]:value_type(0.); return tmp;}

      //! Returns degree \c i term without bounds checking
      KOKKOS_INLINE_FUNCTION
      const_volatile_reference fastAccessCoeff(ordinal_type i) const volatile {
        return s_[i];}

      //! Returns degree \c i term without bounds checking
      KOKKOS_INLINE_FUNCTION
      const_reference fastAccessCoeff(ordinal_type i) const {
        return s_[i];}

      //! Returns degree \c i term without bounds checking
      KOKKOS_INLINE_FUNCTION
      volatile_reference fastAccessCoeff(ordinal_type i) volatile {
        return s_[i];}

      //! Returns degree \c i term without bounds checking
      KOKKOS_INLINE_FUNCTION
      reference fastAccessCoeff(ordinal_type i) {
        return s_[i];}

      /*
      //! Get coefficient term for given dimension and order
      KOKKOS_INLINE_FUNCTION
      reference term(ordinal_type dimension, ordinal_type order) {
        return s_.term(dimension, order); }

      //! Get coefficient term for given dimension and order
      KOKKOS_INLINE_FUNCTION
      const_reference term(ordinal_type dimension, ordinal_type order) const {
        return s_.term(dimension, order); }

      //! Get orders for a given term
      KOKKOS_INLINE_FUNCTION
      Teuchos::Array<ordinal_type> order(ordinal_type term) const {
        return s_.order(term); }
      */

      //! Return iterator to first element of coefficient array
      KOKKOS_INLINE_FUNCTION
      pointer begin() { return s_.coeff(); }

      //! Return iterator to first element of coefficient array
      KOKKOS_INLINE_FUNCTION
      const_pointer begin() const { return s_.coeff(); }

      //! Return iterator to first element of coefficient array
      KOKKOS_INLINE_FUNCTION
      volatile_pointer begin() volatile { return s_.coeff(); }

      //! Return iterator to first element of coefficient array
      KOKKOS_INLINE_FUNCTION
      const_volatile_pointer begin() const volatile { return s_.coeff(); }

      //! Return iterator to first element of coefficient array
      KOKKOS_INLINE_FUNCTION
      const_pointer cbegin() const { return s_.coeff(); }

      //! Return iterator to first element of coefficient array
      KOKKOS_INLINE_FUNCTION
      const_volatile_pointer cbegin() const volatile { return s_.coeff(); }

      //! Return iterator following last element of coefficient array
      KOKKOS_INLINE_FUNCTION
      pointer end() { return s_.coeff() + s_.size(); }

      //! Return iterator following last element of coefficient array
      KOKKOS_INLINE_FUNCTION
      const_pointer end() const { return s_.coeff() + s_.size(); }

      //! Return iterator following last element of coefficient array
      KOKKOS_INLINE_FUNCTION
      volatile_pointer end() volatile { return s_.coeff() + s_.size(); }

      //! Return iterator following last element of coefficient array
      KOKKOS_INLINE_FUNCTION
      const_volatile_pointer end() const volatile { return s_.coeff() + s_.size(); }

      //! Return iterator following last element of coefficient array
      KOKKOS_INLINE_FUNCTION
      const_pointer cend() const { return s_.coeff()+ s_.size(); }

      //! Return iterator following last element of coefficient array
      KOKKOS_INLINE_FUNCTION
      const_volatile_pointer cend() const volatile { return s_.coeff()+ s_.size(); }

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Unary-plus operator
      KOKKOS_INLINE_FUNCTION
      PCE operator + () const { return *this; }

      //! Unary-plus operator
      KOKKOS_INLINE_FUNCTION
      PCE operator + () const volatile { return *this; }

      //! Unary-minus operator
      KOKKOS_INLINE_FUNCTION
      PCE operator - () const;

      //! Unary-minus operator
      KOKKOS_INLINE_FUNCTION
      PCE operator - () const volatile;

      //! Addition-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      PCE& operator += (const value_type x) {
        s_[0] += x;
        return *this;
      }

      //! Addition-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ PCE& operator += (const value_type x) volatile {
        s_[0] += x;
        return const_cast<PCE&>(*this);
      }

      //! Subtraction-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      PCE& operator -= (const value_type x) {
        s_[0] -= x;
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ PCE& operator -= (const value_type x) volatile {
        s_[0] -= x;
        return const_cast<PCE&>(*this);
      }

      //! Multiplication-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      PCE& operator *= (const value_type x);

      //! Multiplication-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ PCE& operator *= (const value_type x) volatile;

      //! Division-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      PCE& operator /= (const value_type x);

      //! Division-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ PCE& operator /= (const value_type x) volatile;

      //! Addition-assignment operator with PCE right-hand-side
      KOKKOS_INLINE_FUNCTION
      PCE& operator += (const PCE& x);

      //! Addition-assignment operator with PCE right-hand-side
      KOKKOS_INLINE_FUNCTION
      PCE& operator += (const volatile PCE& x) {
        return *this += const_cast<const PCE&>(x);
      }

      //! Addition-assignment operator with PCE right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ PCE& operator += (const PCE& x) volatile {
        return const_cast<PCE&>(*this) += x;
      }

      //! Addition-assignment operator with PCE right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ PCE& operator += (const volatile PCE& x) volatile {
        return const_cast<PCE&>(*this) += const_cast<const PCE&>(x);
      }

      //! Subtraction-assignment operator with PCE right-hand-side
      KOKKOS_INLINE_FUNCTION
      PCE& operator -= (const PCE& x);

      //! Subtraction-assignment operator with PCE right-hand-side
      KOKKOS_INLINE_FUNCTION
      PCE& operator -= (const volatile PCE& x) {
        return *this -= const_cast<const PCE&>(x);
      }

      //! Subtraction-assignment operator with PCE right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ PCE& operator -= (const PCE& x) volatile {
        return const_cast<PCE&>(*this) -= x;
      }

      //! Subtraction-assignment operator with PCE right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ PCE& operator -= (const volatile PCE& x) volatile {
        return const_cast<PCE&>(*this) -= const_cast<const PCE&>(x);
      }

      //! Multiplication-assignment operator with PCE right-hand-side
      KOKKOS_INLINE_FUNCTION
      PCE& operator *= (const PCE& x);

      //! Multiplication-assignment operator with PCE right-hand-side
      KOKKOS_INLINE_FUNCTION
      PCE& operator *= (const volatile PCE& x) {
        return *this *= const_cast<const PCE&>(x);
      }

      //! Multiplication-assignment operator with PCE right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ PCE& operator *= (const PCE& x) volatile {
        return const_cast<PCE&>(*this) *= x;
      }

      //! Multiplication-assignment operator with PCE right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ PCE& operator *= (const volatile PCE& x) volatile {
        return const_cast<PCE&>(*this) *= const_cast<const PCE&>(x);
      }

      //! Division-assignment operator with PCE right-hand-side
      KOKKOS_INLINE_FUNCTION
      PCE& operator /= (const PCE& x);

      //! Division-assignment operator with PCE right-hand-side
      KOKKOS_INLINE_FUNCTION
      PCE& operator /= (const volatile PCE& x) {
        return *this /= const_cast<const PCE&>(x);
      }

      //! Division-assignment operator with PCE right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ PCE& operator /= (const PCE& x) volatile {
        return const_cast<PCE&>(*this) /= x;
      }

      //! Division-assignment operator with PCE right-hand-side
      KOKKOS_INLINE_FUNCTION
      /*volatile*/ PCE& operator /= (const volatile PCE& x) volatile {
        return const_cast<PCE&>(*this) /= const_cast<const PCE&>(x);
      }

      //! Prefix ++
      KOKKOS_INLINE_FUNCTION
      PCE& operator++() {
        ++(s_[0]);
        return *this;
      }

      //! Prefix ++
      KOKKOS_INLINE_FUNCTION
      volatile PCE& operator++() volatile {
        ++(s_[0]);
        return *this;
      }

      //! Postfix ++
      KOKKOS_INLINE_FUNCTION
      PCE operator++(int) {
        PCE tmp(*this);
        ++(*this);
        return tmp;
      }

      //! Postfix ++
      KOKKOS_INLINE_FUNCTION
      PCE operator++(int) volatile {
        PCE tmp(*this);
        ++(*this);
        return tmp;
      }

      //! Prefix --
      KOKKOS_INLINE_FUNCTION
      PCE& operator--() {
        --(s_[0]);
        return *this;
      }

      //! Prefix --
      KOKKOS_INLINE_FUNCTION
      volatile PCE& operator--() volatile {
        --(s_[0]);
        return *this;
      }

      //! Postfix --
      KOKKOS_INLINE_FUNCTION
      PCE operator--(int) {
        PCE tmp(*this);
        --(*this);
        return tmp;
      }

      //! Postfix --
      KOKKOS_INLINE_FUNCTION
      PCE operator--(int) volatile {
        PCE tmp(*this);
        --(*this);
        return tmp;
      }



      //@}

    protected:

      typedef typename my_cijk_type::size_type cijk_size_type;

      //! Cijk class
      my_cijk_type cijk_;

      Storage s_;

    private:

      //PCE division using CG with diagonal preconditioning
//      KOKKOS_INLINE_FUNCTION
//      void CG_divide(const PCE& a, const PCE& b, PCE& c);


    }; // class PCE

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    void CG_divide(const PCE<Storage>& a, const PCE<Storage>& b, PCE<Storage>& c);

    // Operations
    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    operator+(const PCE<Storage>& a, const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    operator+(const volatile PCE<Storage>& a, const PCE<Storage>& b) {
      return const_cast<const PCE<Storage>&>(a) + b;
    }

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    operator+(const PCE<Storage>& a, const volatile PCE<Storage>& b) {
      return a + const_cast<const PCE<Storage>&>(b);
    }

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    operator+(const volatile PCE<Storage>& a, const volatile PCE<Storage>& b) {
      return const_cast<const PCE<Storage>&>(a) +
             const_cast<const PCE<Storage>&>(b);
    }

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    operator+(const typename PCE<Storage>::value_type& a,
              const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    operator+(const PCE<Storage>& a,
              const typename PCE<Storage>::value_type& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    operator-(const PCE<Storage>& a, const PCE<Storage>& b);

    template <typename Storage> PCE<Storage>
    KOKKOS_INLINE_FUNCTION
    operator-(const typename PCE<Storage>::value_type& a,
              const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    operator-(const PCE<Storage>& a,
              const typename PCE<Storage>::value_type& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    operator*(const PCE<Storage>& a, const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    operator*(const typename PCE<Storage>::value_type& a,
              const PCE<Storage>& b);

    template <typename Storage> PCE<Storage>
    KOKKOS_INLINE_FUNCTION
    operator*(const PCE<Storage>& a,
              const typename PCE<Storage>::value_type& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    operator/(const PCE<Storage>& a, const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    operator/(const typename PCE<Storage>::value_type& a,
              const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    operator/(const PCE<Storage>& a,
              const typename PCE<Storage>::value_type& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    exp(const PCE<Storage>& a);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    log(const PCE<Storage>& a);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    void
    log(PCE<Storage>& c, const PCE<Storage>& a);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    log10(const PCE<Storage>& a);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    sqrt(const PCE<Storage>& a);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    cbrt(const PCE<Storage>& a);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    pow(const PCE<Storage>& a, const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    pow(const typename PCE<Storage>::value_type& a, const PCE<Storage>& b);

    template <typename Storage> PCE<Storage>
    KOKKOS_INLINE_FUNCTION
    pow(const PCE<Storage>& a,
        const typename PCE<Storage>::value_type& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    cos(const PCE<Storage>& a);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    sin(const PCE<Storage>& a);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    tan(const PCE<Storage>& a);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    cosh(const PCE<Storage>& a);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    sinh(const PCE<Storage>& a);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    tanh(const PCE<Storage>& a);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    acos(const PCE<Storage>& a);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    asin(const PCE<Storage>& a);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    atan(const PCE<Storage>& a);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    atan2(const PCE<Storage>& a, const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    atan2(const typename PCE<Storage>::value_type& a,
          const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    atan2(const PCE<Storage>& a,
          const typename PCE<Storage>::value_type& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    acosh(const PCE<Storage>& a);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    asinh(const PCE<Storage>& a);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    atanh(const PCE<Storage>& a);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    abs(const PCE<Storage>& a);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    fabs(const PCE<Storage>& a);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    ceil(const PCE<Storage>& a);

    // template <typename Storage>
    // KOKKOS_INLINE_FUNCTION
    // PCE<Storage>
    // max(const PCE<Storage>& a, const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    max(const typename PCE<Storage>::value_type& a,
        const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    max(const PCE<Storage>& a,
        const typename PCE<Storage>::value_type& b);

    // template <typename Storage>
    // KOKKOS_INLINE_FUNCTION
    // PCE<Storage>
    // min(const PCE<Storage>& a, const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    min(const typename PCE<Storage>::value_type& a,
        const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    PCE<Storage>
    min(const PCE<Storage>& a,
        const typename PCE<Storage>::value_type& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    bool
    operator==(const PCE<Storage>& a, const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    bool
    operator==(const typename PCE<Storage>::value_type& a,
               const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    bool
    operator==(const PCE<Storage>& a,
               const typename PCE<Storage>::value_type& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    bool
    operator!=(const PCE<Storage>& a, const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    bool
    operator!=(const typename PCE<Storage>::value_type& a,
               const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    bool
    operator!=(const PCE<Storage>& a,
               const typename PCE<Storage>::value_type& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    bool
    operator<=(const PCE<Storage>& a,
               const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    bool
    operator<=(const typename PCE<Storage>::value_type& a,
               const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    bool
    operator<=(const PCE<Storage>& a,
               const typename PCE<Storage>::value_type& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    bool
    operator>=(const PCE<Storage>& a, const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    bool
    operator>=(const typename PCE<Storage>::value_type& a,
               const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    bool
    operator>=(const PCE<Storage>& a,
               const typename PCE<Storage>::value_type& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    bool
    operator<(const PCE<Storage>& a, const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    bool
    operator<(const typename PCE<Storage>::value_type& a,
              const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    bool
    operator<(const PCE<Storage>& a,
              const typename PCE<Storage>::value_type& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    bool
    operator>(const PCE<Storage>& a, const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    bool
    operator>(const typename PCE<Storage>::value_type& a,
              const PCE<Storage>& b);

    template <typename Storage>
    KOKKOS_INLINE_FUNCTION
    bool
    operator>(const PCE<Storage>& a,
              const typename PCE<Storage>::value_type& b);

    template <typename Storage>
    std::ostream&
    operator << (std::ostream& os, const PCE<Storage>& a);

    template <typename Storage>
    std::istream&
    operator >> (std::istream& os, PCE<Storage>& a);

    /**\brief  Define a partition of a View of Sacado::UQ::PCE type */
    struct PCEPartition {
      unsigned begin ;
      unsigned end ;

      template< typename iType0 , typename iType1 >
      KOKKOS_INLINE_FUNCTION
      PCEPartition( const iType0 & i0 , const iType1 & i1 ) :
        begin(i0), end(i1) {}
    };

/*   template <typename Storage>
   KOKKOS_INLINE_FUNCTION
   void
   PCE<Storage>::
   CG_divide(const PCE<Storage>& a, const PCE<Storage>& b, PCE<Storage>& c) {
//     typedef typename PCE<Storage>::ordinal_type ordinal_type;
//     typedef typename PCE<Storage>::value_type value_type;

     const ordinal_type size = c.size();
     //Needed scalars
     value_type alpha, beta, rTz, rTz_old, resid;

     //Needed temporary PCEs
     PCE<Storage> r(a.cijk(),size);
     PCE<Storage> p(a.cijk(),size);
     PCE<Storage> bp(a.cijk(),size);
     PCE<Storage> z(a.cijk(),size);

     //compute residual = a - b*c
     r =  a - b*c;
     z = r/b.coeff(0);
     p = z;
     resid = r.two_norm();
     //Compute <r,z>=rTz (L2 inner product)
     rTz = r.inner_product(z);
     ordinal_type k = 0;
     value_type tol = 1e-6;
     while ( resid > tol && k < 100){
       //Compute b*p
       bp = b*p;

       //Compute alpha = <r,z>/<p,b*p>
       alpha = rTz/p.inner_product(bp);

       //Check alpha is positive!
       //Update the solution c = c + alpha*p
       c = c + alpha*p;
       rTz_old = rTz;

       //Compute the new residual r = r - alpha*b*p
       r = r - alpha*bp;
       resid = r.two_norm();

       //Compute beta = rTz_new/ rTz_old and new p
       z = r/b.coeff(0);
       rTz = r.inner_product(z);
       beta = rTz/rTz_old;
       p = z + beta*p;
       k++;
      }
//  return c;
   }
*/

    template <typename S>
    void memcpy(PCE<S>* dst, const PCE<S>* src, const size_t sz) {
      const size_t n = sz / sizeof(PCE<S>);
      for (size_t i=0; i<n; ++i)
        dst[i] = src[i];
    }

  } // namespace PCE

  //! Trait class to determine if a scalar type is a PCE
  template <typename T> struct is_uq_pce {
    static const bool value = false;
  };
  template <typename S> struct is_uq_pce< UQ::PCE<S> > {
    static const bool value = true;
  };
  template <typename T> struct is_uq_pce< const T > {
    static const bool value = is_uq_pce<T>::value;
  };
  template <typename T> struct is_uq_pce< T* > {
    static const bool value = is_uq_pce<T>::value;
  };
  template <typename T> struct is_uq_pce< T[] > {
    static const bool value = is_uq_pce<T>::value;
  };
  template <typename T, unsigned N> struct is_uq_pce< T[N] > {
    static const bool value = is_uq_pce<T>::value;
  };

  // Utility function to see if a PCE is really a constant
  template <typename Storage>
  KOKKOS_INLINE_FUNCTION
  bool is_constant(const Sacado::UQ::PCE<Storage>& x)
  {
    typedef typename Storage::ordinal_type ordinal_type;
    typedef typename Storage::value_type value_type;

    // All size-1 expansions are constants
    const ordinal_type sz = x.size();
    if (sz == 1) return true;

    // Maybe use a tolerance????
    const value_type zero = 0;
    for (ordinal_type i=1; i<sz; ++i)
      if (x.fastAccessCoeff(i) != zero) return false;

    return true;
  }

} // namespace Sacado

#include "Sacado_UQ_PCE_Traits.hpp"
#include "Sacado_UQ_PCE_Imp.hpp"

#include "Kokkos_NumericTraits.hpp"

namespace Kokkos {

template <typename Storage>
struct reduction_identity< Sacado::UQ::PCE<Storage> > {
  typedef Sacado::UQ::PCE<Storage> pce;
  typedef typename Storage::value_type scalar;
  typedef reduction_identity<scalar> RIS;
  KOKKOS_FORCEINLINE_FUNCTION constexpr static pce sum()  {
    return pce(RIS::sum());
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static pce prod() {
    return pce(RIS::prod());
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static pce max()  {
    return pce(RIS::max());
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static pce min()  {
    return pce(RIS::min());
  }
};

namespace Impl {
  template <typename Storage>
  struct promote<Sacado::UQ::PCE<Storage>,false> {
    using type = typename Sacado::UQ::PCE<Storage>;
  };
}

}

#endif // HAVE_STOKHOS_SACADO

#endif // SACADO_PCE_ORTHOGPOLY_HPP
