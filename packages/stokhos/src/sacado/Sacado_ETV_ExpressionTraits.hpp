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
// @HEADER

#ifndef SACADO_ETV_EXPRESSIONTRAITS_HPP
#define SACADO_ETV_EXPRESSIONTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace ETV {
    template <typename T> class Expr;
  }
}

namespace Sacado {

  // We don't specialize Promote because otherwise we can get ambiguous
  // partial specializations with the Vector classes.

  //! Specialization of %ScalarType to Expr types
  template <typename T>
  struct ScalarType< ETV::Expr<T> > {
    typedef typename ScalarType< typename ETV::Expr<T>::value_type >::type type;
  };

  //! Specialization of %ValueType to Expr types
  template <typename T>
  struct ValueType< ETV::Expr<T> > {
    typedef typename ETV::Expr<T>::value_type type;
  };

  //! Specialization of %IsADType to Expr types
  template <typename T>
  struct IsADType< ETV::Expr<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to Expr types
  template <typename T>
  struct IsScalarType< ETV::Expr<T> > {
    static const bool value = false;
  };

  //! Specialization of %Value to Expr types
  template <typename T>
  struct Value< ETV::Expr<T> > {
    typedef typename ValueType< ETV::Expr<T> >::type value_type;
    static const value_type& eval(const ETV::Expr<T>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to Expr types
  template <typename T>
  struct ScalarValue< ETV::Expr<T> > {
    typedef typename ValueType< ETV::Expr<T> >::type value_type;
    typedef typename ScalarType< ETV::Expr<T> >::type scalar_type;
    static const scalar_type& eval(const ETV::Expr<T>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

} // namespace Sacado

#endif // SACADO_ETV_EXPRESSIONTRAITS_HPP
