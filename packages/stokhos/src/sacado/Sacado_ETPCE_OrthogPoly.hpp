// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SACADO_ETPCE_ORTHOGPOLY_HPP
#define SACADO_ETPCE_ORTHOGPOLY_HPP

#include "Stokhos_ConfigDefs.h"

#ifdef HAVE_STOKHOS_SACADO

#include "Teuchos_RCP.hpp"

#include "Sacado_Traits.hpp"
#include "Sacado_Handle.hpp"

#include "Stokhos_OrthogPolyExpansion.hpp"
#include "Stokhos_QuadOrthogPolyExpansion.hpp"
#include "Stokhos_ConstantOrthogPolyExpansion.hpp"
#include "Stokhos_OrthogPolyApprox.hpp"

#include "Sacado_mpl_apply.hpp"

#include <ostream>      // for std::ostream

#ifdef HAVE_STOKHOS_THRUST
#include "thrust/tuple.h"
#endif

#ifndef KERNEL_PREFIX
#define KERNEL_PREFIX
#endif

namespace Sacado {

  //! Namespace for expression templated polynomial chaos expansion classes
  namespace ETPCE {

    template <int k, typename T> KERNEL_PREFIX T&
    get(T* a) { return a[k]; }
    template <int k, typename T> KERNEL_PREFIX const T&
    get(const T* a) { return a[k]; }

    template <int k, int N, typename T> KERNEL_PREFIX T&
    get(T a[N]) { return a[k]; }
    template <int k, int N, typename T> KERNEL_PREFIX const T&
    get(const T a[N]) { return a[k]; }

    //! Wrapper for a generic expression template
    /*!
     * This template class serves as a wrapper for all expression
     * template classes.
     */
    template <typename ExprT> class Expr {};

    // Forward declaration
    template <typename T, typename S> class OrthogPoly;

    //! Generalized polynomial chaos expansion class implementation
    template <typename T, typename Storage >
    class OrthogPolyImpl {
    public:

      //! Typename of values
      typedef T value_type;

      //! Typename of scalar's (which may be different from T)
      typedef typename ScalarType<T>::type scalar_type;

      //! Typename of ordinals
      typedef int ordinal_type;

      //! Typename of storage class
      typedef Storage storage_type;

      //! Basis type
      typedef Stokhos::OrthogPolyBasis<ordinal_type,T> basis_type;

      //! Quad Expansion type
      typedef Stokhos::OrthogPolyExpansion<ordinal_type,T,Storage> expansion_type;

      //! Quad Expansion type
      typedef Stokhos::QuadOrthogPolyExpansion<ordinal_type,T,Storage> quad_expansion_type;

      //! Stokhos approximation type
      typedef Stokhos::OrthogPolyApprox<ordinal_type,T,Storage> approx_type;

      typedef typename approx_type::pointer pointer;
      typedef typename approx_type::const_pointer const_pointer;
      typedef typename approx_type::reference reference;
      typedef typename approx_type::const_reference const_reference;


      //! Default constructor
      /*!
       * Sets size to 1 and first coefficient to 0 (represents a constant).
       */
      OrthogPolyImpl();

      //! Constructor with supplied value \c x
      /*!
       * Sets size to 1 and first coefficient to x (represents a constant).
       */
      OrthogPolyImpl(const value_type& x);

      //! Constructor with expansion \c expansion (General case)
      /*!
       * Creates array of correct size and initializes coeffiencts to 0.
       */
      OrthogPolyImpl(const Teuchos::RCP<expansion_type>& expansion);

      //! Constructor with expansion \c expansion and specified size \c sz
      /*!
       * Creates array of size \c sz and initializes coeffiencts to 0.
       */
      OrthogPolyImpl(const Teuchos::RCP<expansion_type>& expansion,
                     ordinal_type sz);

      //! Copy constructor
      OrthogPolyImpl(const OrthogPolyImpl& x);

      //! Copy constructor from any Expression object
      template <typename S> OrthogPolyImpl(const Expr<S>& x);

      //! Destructor
      ~OrthogPolyImpl() {}

      //! Initialize coefficients to value
      void init(const T& v) { th_->init(v); }

      //! Initialize coefficients to an array of values
      void init(const T* v) { th_->init(v); }

      //! Initialize coefficients from an OrthogPolyImpl with different storage
      template <typename S>
      void init(const OrthogPolyImpl<T,S>& v) { th_->init(v.getOrthogPolyApprox()); }

      //! Load coefficients to an array of values
      void load(T* v) { th_->load(v); }

      //! Load coefficients into an OrthogPolyImpl with different storage
      template <typename S>
      void load(OrthogPolyImpl<T,S>& v) { th_->load(v.getOrthogPolyApprox()); }

      //! Reset expansion
      /*!
       * May change size of array.  Coefficients are preserved.
       */
      void reset(const Teuchos::RCP<expansion_type>& expansion);

      //! Reset expansion and size
      /*!
       * Coefficients are preserved.
       */
      void reset(const Teuchos::RCP<expansion_type>& expansion,
                 ordinal_type sz);

      //! Prepare polynomial for writing
      /*!
       * This method prepares the polynomial for writing through coeff() and
       * fastAccessCoeff() member functions.  It ensures the handle for the
       * %Hermite coefficients is not shared among any other %Hermite polynomial
       * objects.  If the handle is not shared it does nothing, so there
       * is no cost in calling this method in this case.  If the handle is
       * shared and this method is not called, any changes to the coefficients
       * by coeff() or fastAccessCoeff() may change other polynomial objects.
       */
      void copyForWrite() { th_.makeOwnCopy(); }

      //! Evaluate polynomial approximation at a point
      value_type evaluate(const Teuchos::Array<value_type>& point) const;

      //! Evaluate polynomial approximation at a point with given basis values
      value_type evaluate(const Teuchos::Array<value_type>& point,
                          const Teuchos::Array<value_type>& bvals) const;

      //! Compute mean of expansion
      value_type mean() const {return th_->mean(); }

      //! Compute standard deviation of expansion
      value_type standard_deviation() const { return th_->standard_deviation(); }

      //! Compute the two-norm of expansion
      value_type two_norm() const { return th_->two_norm(); }

      //! Compute the squared two-norm of expansion
      value_type two_norm_squared() const { return th_->two_norm_squared(); }

      //! Compute the L2 inner product of 2 PCEs
      value_type inner_product(const OrthogPolyImpl& b) const {
        return th_->inner_product(b.getOrthogPolyApprox()); }

      //! Print approximation in basis
      std::ostream& print(std::ostream& os) const { return th_->print(os); }

      //! Returns whether two PCE objects have the same values
      template <typename S> bool isEqualTo(const Expr<S>& x) const;

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      OrthogPolyImpl& operator=(const value_type& val);

      //! Assignment with OrthogPolyImpl right-hand-side
      OrthogPolyImpl& operator=(const OrthogPolyImpl& x);

      //! Assignment with any expression right-hand-side
      template <typename S>
      OrthogPolyImpl& operator=(const Expr<S>& x);

      //@}

      /*!
       * Accessor methods
       */
      //@{

      //! Get basis
      Teuchos::RCP<basis_type> basis() const { return th_->basis(); }

      //! Get expansion
      Teuchos::RCP<expansion_type> expansion() const { return expansion_; }

      //! Get quad expansion
      Teuchos::RCP<quad_expansion_type> quad_expansion() const { return quad_expansion_; }

      //@}

      /*!
       * @name Value accessor methods
       */
      //@{

      //! Returns value
      const_reference val() const { return (*th_)[0]; }

      //! Returns value
      reference val() { return (*th_)[0]; }

      //@}

      /*!
       * @name Coefficient accessor methods
       */
      //@{

      //! Returns size of polynomial
      ordinal_type size() const { return th_->size();}

      //! Returns true if polynomial has size >= sz
      bool hasFastAccess(ordinal_type sz) const { return th_->size()>=sz;}

      //! Returns Hermite coefficient array
      const_pointer coeff() const { return th_->coeff();}

      //! Returns Hermite coefficient array
      pointer coeff() { return th_->coeff();}

      //! Returns degree \c i term with bounds checking
      value_type coeff(ordinal_type i) const {
        return i<th_->size() ? (*th_)[i]:value_type(0.); }

      //! Returns degree \c i term without bounds checking
      reference fastAccessCoeff(ordinal_type i) { return (*th_)[i];}

      //! Returns degree \c i term without bounds checking
      value_type fastAccessCoeff(ordinal_type i) const { return (*th_)[i];}

      //! Get coefficient term for given dimension and order
      reference term(ordinal_type dimension, ordinal_type order) {
        return th_->term(dimension, order); }

      //! Get coefficient term for given dimension and order
      const_reference term(ordinal_type dimension, ordinal_type order) const {
        return th_->term(dimension, order); }

      //! Get orders for a given term
      Teuchos::Array<ordinal_type> order(ordinal_type term) const {
        return th_->order(term); }

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Addition-assignment operator with constant right-hand-side
      OrthogPolyImpl& operator += (const value_type& x);

      //! Subtraction-assignment operator with constant right-hand-side
      OrthogPolyImpl& operator -= (const value_type& x);

      //! Multiplication-assignment operator with constant right-hand-side
      OrthogPolyImpl& operator *= (const value_type& x);

      //! Division-assignment operator with constant right-hand-side
      OrthogPolyImpl& operator /= (const value_type& x);

      //@}

      //! Get underlying Stokhos::OrthogPolyApprox
      const approx_type& getOrthogPolyApprox() const { return *th_; }

      //! Get underlying Stokhos::OrthogPolyApprox
      approx_type& getOrthogPolyApprox() { return *th_; }

    protected:

      //! Copy from Expression object
      template <typename S> void expressionCopy(const Expr<S>& x);

    protected:

      //! Expansion class
      Teuchos::RCP<expansion_type> expansion_;

      //! Cast of expansion class to QuadExpansion
      Teuchos::RCP<quad_expansion_type> quad_expansion_;

      //! Constant expansion class for constants
      Teuchos::RCP<expansion_type> const_expansion_;

      //! Handle to underlying OrthogPolyApprox
      Sacado::Handle< Stokhos::OrthogPolyApprox<int,value_type,Storage> > th_;

    }; // class OrthogPolyImpl

    //! OrthogPolyImpl expression template specialization
    /*!
     * This template class represents a simple OrthogPolyImpl expression and
     * mixes-in the OrthogPolyImpl interface and the expression template
     * interface.
     */
    template <typename T, typename Storage>
    class Expr< OrthogPolyImpl<T,Storage> > :
        public OrthogPolyImpl<T,Storage> {

    public:

      //! Typename of values
      typedef typename OrthogPolyImpl<T,Storage>::value_type value_type;
      typedef typename OrthogPolyImpl<T,Storage>::scalar_type scalar_type;
      typedef typename OrthogPolyImpl<T,Storage>::approx_type approx_type;
      typedef typename OrthogPolyImpl<T,Storage>::storage_type storage_type;
      typedef typename OrthogPolyImpl<T,Storage>::const_reference const_reference;

      typedef OrthogPoly<T,Storage> base_expr_type;

      //! Number of arguments
      static const int num_args = 1;

      //! Default constructor
      Expr() :
        OrthogPolyImpl<T,Storage>() {}

      //! Constructor with supplied value \c x
      Expr(const T & x) :
        OrthogPolyImpl<T,Storage>(x) {}

      //! Constructor with expansion \c expansion (General case)
      Expr(const Teuchos::RCP<typename OrthogPolyImpl<T,Storage>::expansion_type>& expansion) :
        OrthogPolyImpl<T,Storage>(expansion) {}

      //! Constructor with expansion \c expansion and specified size \c sz
      Expr(const Teuchos::RCP<typename OrthogPolyImpl<T,Storage>::expansion_type>& expansion,
           typename OrthogPolyImpl<T,Storage>::ordinal_type sz) :
        OrthogPolyImpl<T,Storage>(expansion, sz) {}

      //! Copy constructor
      Expr(const Expr& x) :
        OrthogPolyImpl<T,Storage>(static_cast<const OrthogPolyImpl<T,Storage>&>(x)) {}

      //! Copy constructor
      Expr(const OrthogPolyImpl<T,Storage>& x) :
        OrthogPolyImpl<T,Storage>(x) {}

      //! Copy constructor from any Expression object
      template <typename S> Expr(const Expr<S>& x) :
        OrthogPolyImpl<T,Storage>(x) {}

      //! Destructor
      ~Expr() {}

      const approx_type& getArg(int i) const {
        return this->getOrthogPolyApprox(); }

      bool has_fast_access(int sz) const { return this->size() >= sz; }

      bool has_nonconst_expansion() const {
        return this->expansion_ != this->const_expansion_;
      }

      int order() const { return this->size() == 1 ? 0 : 1; }

      value_type fast_higher_order_coeff(int i) const {
        return this->fastAccessCoeff(i);
      }

      value_type higher_order_coeff(int i) const {
        return this->coeff(i);
      }

      template <int offset, typename tuple_type>
      KERNEL_PREFIX
      value_type eval_sample(tuple_type x) const {
        return get<offset>(x);
      }

      std::string name() const { return "x"; }

    }; // class Expr<OrthogPolyImpl>

  } // namespace ETPCE

} // namespace Sacado

#include "Sacado_ETPCE_ExpressionTraits.hpp"
#include "Sacado_ETPCE_OrthogPolyTraits.hpp"
#include "Sacado_ETPCE_OrthogPolyImp.hpp"
#include "Sacado_ETPCE_OrthogPolyOps.hpp"

namespace Sacado {

  namespace ETPCE {

    //! Generalized polynomial chaos expansion class
    template <typename T, typename Storage>
    class OrthogPoly : public Expr< OrthogPolyImpl<T,Storage> > {
    public:

      //! Typename of values
      typedef typename OrthogPolyImpl<T,Storage>::value_type value_type;

      //! Typename of scalars
      typedef typename OrthogPolyImpl<T,Storage>::scalar_type scalar_type;

      //! Typename of ordinals
      typedef typename OrthogPolyImpl<T,Storage>::ordinal_type ordinal_type;

      //! Typename of storage class
      typedef typename OrthogPolyImpl<T,Storage>::storage_type storage_type;

      //! Basis type
      typedef typename OrthogPolyImpl<T,Storage>::basis_type basis_type;

      //! Expansion type
      typedef typename OrthogPolyImpl<T,Storage>::expansion_type expansion_type;

      //! Stokhos approximation type
      typedef typename OrthogPolyImpl<T,Storage>::approx_type approx_type;

      typedef typename OrthogPolyImpl<T,Storage>::pointer pointer;
      typedef typename OrthogPolyImpl<T,Storage>::const_pointer const_pointer;
      typedef typename OrthogPolyImpl<T,Storage>::reference reference;
      typedef typename OrthogPolyImpl<T,Storage>::const_reference const_reference;

      //! Turn OrthogPoly into a meta-function class usable with mpl::apply
      template <typename S>
      struct apply {
        typedef typename Sacado::mpl::apply<Storage,ordinal_type,S>::type storage_type;
        typedef OrthogPoly<S,storage_type> type;
      };

      //! Default constructor
      /*!
       * Sets size to 1 and first coefficient to 0 (represents a constant).
       */
      OrthogPoly() :
        Expr< OrthogPolyImpl<T,Storage> >() {}

      //! Constructor with supplied value \c x
      /*!
       * Sets size to 1 and first coefficient to x (represents a constant).
       */
      OrthogPoly(const value_type& x) :
        Expr< OrthogPolyImpl<T,Storage> >(x) {}

      //! Constructor with expansion \c expansion (General case)
      /*!
       * Creates array of correct size and initializes coeffiencts to 0.
       */
      OrthogPoly(const Teuchos::RCP<expansion_type>& expansion) :
        Expr< OrthogPolyImpl<T,Storage> >(expansion) {}

      //! Constructor with expansion \c expansion and specified size \c sz
      /*!
       * Creates array of size \c sz and initializes coeffiencts to 0.
       */
      OrthogPoly(const Teuchos::RCP<expansion_type>& expansion,
                     ordinal_type sz) :
        Expr< OrthogPolyImpl<T,Storage> >(expansion, sz) {}

      //! Copy constructor
      OrthogPoly(const OrthogPoly& x) :
        Expr< OrthogPolyImpl<T,Storage> >(x) {}

      //! Copy constructor from any Expression object
      template <typename S> OrthogPoly(const Expr<S>& x) :
        Expr< OrthogPolyImpl<T,Storage> >(x) {}

      //! Destructor
      ~OrthogPoly() {}

      //! Assignment operator with constant right-hand-side
      OrthogPoly& operator=(const value_type& val) {
        OrthogPolyImpl<T,Storage>::operator=(val);
        return *this;
      }

      //! Assignment with OrthogPoly right-hand-side
      OrthogPoly& operator=(const OrthogPoly& x) {
        OrthogPolyImpl<T,Storage>::operator=(static_cast<const OrthogPolyImpl<T,Storage>&>(x));
        return *this;
      }

      //! Assignment with Expr< OrthogPolyImpl > right-hand-side
      OrthogPoly& operator=(const Expr< OrthogPolyImpl<T,Storage> >& x) {
        OrthogPolyImpl<T,Storage>::operator=(static_cast<const OrthogPolyImpl<T,Storage>&>(x));
        return *this;
      }

      //! Assignment with any expression right-hand-side
      template <typename S>
      OrthogPoly& operator=(const Expr<S>& x) {
        OrthogPolyImpl<T,Storage>::operator=(x);
        return *this;
      }

      //@{

      //! Addition-assignment operator with constant right-hand-side
      OrthogPoly& operator += (const value_type& x) {
        OrthogPolyImpl<T,Storage>::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      OrthogPoly& operator -= (const value_type& x) {
        OrthogPolyImpl<T,Storage>::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      OrthogPoly& operator *= (const value_type& x) {
        OrthogPolyImpl<T,Storage>::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      OrthogPoly& operator /= (const value_type& x) {
        OrthogPolyImpl<T,Storage>::operator/=(x);
        return *this;
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      OrthogPoly& operator += (const Expr<S>& x) {
        *this = *this + x;
        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      OrthogPoly& operator -= (const Expr<S>& x) {
        *this = *this - x;
        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      OrthogPoly& operator *= (const Expr<S>& x) {
        *this = *this * x;
        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      OrthogPoly& operator /= (const Expr<S>& x) {
        *this = *this / x;
        return *this;
      }

      //@}

    }; // class OrthogPoly

  } // namespace ETPCE

  template <typename T>
  struct IsExpr< ETPCE::Expr<T> > {
    static const bool value = true;
  };

  template <typename T>
  struct BaseExprType< ETPCE::Expr<T> > {
    typedef typename ETPCE::Expr<T>::base_expr_type type;
  };

  template <typename T, typename S>
  struct IsExpr< ETPCE::OrthogPoly<T,S> > {
    static const bool value = true;
  };

  template <typename T, typename S>
  struct BaseExprType< ETPCE::OrthogPoly<T,S> > {
    typedef ETPCE::OrthogPoly<T,S> type;
  };

} // namespace Sacado

#endif // HAVE_STOKHOS_SACADO

#endif // SACADO_ETPCE_ORTHOGPOLY_HPP
