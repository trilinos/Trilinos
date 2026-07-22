// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_CACHEFAD_GENERALFADEXPR_HPP
#define SACADO_CACHEFAD_GENERALFADEXPR_HPP

#include "Sacado_CacheFad_GeneralFad.hpp"

namespace Sacado {

  namespace CacheFad {

    //! GeneralFad expression template specialization
    /*!
     * This template class represents a simple GeneralFad expression and
     * mixes-in the GeneralFad interface and the expression template
     * interface using the caching expression templates.
     */
    template <typename T, typename Storage>
    class Expr< GeneralFad<T,Storage> > : public GeneralFad<T,Storage> {

    public:

      //! Typename of values
      typedef typename GeneralFad<T,Storage>::value_type value_type;

      //! Typename of scalar's (which may be different from value_type)
      typedef typename GeneralFad<T,Storage>::scalar_type scalar_type;

      //! Typename of base-expressions
      typedef typename BaseExpr< GeneralFad<T,Storage> >::type base_expr_type;

      //! Default constructor
      SACADO_INLINE_FUNCTION
      Expr() :
        GeneralFad<T,Storage>() {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      template <typename S>
      SACADO_INLINE_FUNCTION
      Expr(const S & x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        GeneralFad<T,Storage>(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      SACADO_INLINE_FUNCTION
      Expr(const int sz, const T & x, const DerivInit zero_out = InitDerivArray) :
        GeneralFad<T,Storage>(sz,x,zero_out) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      SACADO_INLINE_FUNCTION
      Expr(const int sz, const int i, const T & x) :
        GeneralFad<T,Storage>(sz,i,x) {}

      //! Constructor with supplied storage \c s
      SACADO_INLINE_FUNCTION
      Expr(const Storage& s) :
        GeneralFad<T,Storage>(s) {}

      //! Copy constructor
      SACADO_INLINE_FUNCTION
      Expr(const Expr& x) :
        GeneralFad<T,Storage>(static_cast<const GeneralFad<T,Storage>&>(x)) {}

      //! Copy constructor from any Expression object
      template <typename S>
      SACADO_INLINE_FUNCTION
      Expr(const Expr<S>& x, SACADO_ENABLE_EXPR_CTOR_DECL) :
        GeneralFad<T,Storage>(x) {}

      //! Destructor
      SACADO_INLINE_FUNCTION
      ~Expr() {}

    }; // class Expr<GeneralFad>

  } // namespace CacheFad

} // namespace Sacado

#include "Sacado_CacheFad_Ops.hpp"

#endif // SACADO_CACHEFAD_GENERALFADEXPR_HPP
