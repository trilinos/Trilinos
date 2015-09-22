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

#ifndef SACADO_ETPCE_EXPRESSIONTRAITS_HPP
#define SACADO_ETPCE_EXPRESSIONTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace ETPCE {
    template <typename T> class Expr;
  }
}

namespace Sacado {

  // We don't specialize Promote because otherwise we can get ambiguous
  // partial specializations with the OrthogPoly classes.

  //! Specialization of %ScalarType to Expr types
  template <typename T>
  struct ScalarType< ETPCE::Expr<T> > {
    typedef typename ScalarType< typename ETPCE::Expr<T>::value_type >::type type;
  };

  //! Specialization of %ValueType to Expr types
  template <typename T>
  struct ValueType< ETPCE::Expr<T> > {
    typedef typename ETPCE::Expr<T>::value_type type;
  };

  //! Specialization of %IsADType to Expr types
  template <typename T>
  struct IsADType< ETPCE::Expr<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to Expr types
  template <typename T>
  struct IsScalarType< ETPCE::Expr<T> > {
    static const bool value = false;
  };

  //! Specialization of %Value to Expr types
  template <typename T>
  struct Value< ETPCE::Expr<T> > {
    typedef typename ValueType< ETPCE::Expr<T> >::type value_type;
    static const value_type& eval(const ETPCE::Expr<T>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to Expr types
  template <typename T>
  struct ScalarValue< ETPCE::Expr<T> > {
    typedef typename ValueType< ETPCE::Expr<T> >::type value_type;
    typedef typename ScalarType< ETPCE::Expr<T> >::type scalar_type;
    static const scalar_type& eval(const ETPCE::Expr<T>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

} // namespace Sacado

#endif // SACADO_ETPCE_EXPRESSIONTRAITS_HPP
