// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_EXP_EXTENDER_HPP
#define SACADO_FAD_EXP_EXTENDER_HPP

#include "Sacado_Fad_Expression.hpp"

namespace Sacado {

  namespace Fad {
  namespace Exp {

    //! Extension class for extending interface of its argument
    template <typename T, typename Enabled = void>
    class Extender : public T {
    public:
      using T::T;

      // Default expression template specialization
      typedef ExprSpecDefault expr_spec_type;
    };

  } // namespace Exp
  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_EXP_EXTENDER_HPP
