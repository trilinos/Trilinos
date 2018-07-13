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
// @HEADER

#ifndef SACADO_FAD_DFAD_HPP
#define SACADO_FAD_DFAD_HPP

#include "Sacado_ConfigDefs.h"

#ifdef SACADO_NEW_FAD_DESIGN_IS_DEFAULT

#include "Sacado_Fad_Exp_DFad.hpp"

namespace Sacado {
  namespace Fad {
    template <typename T>
    using DFad = Exp::GeneralFad< Exp::DynamicStorage<T> >;
  }
}

#else

#include "Sacado_Fad_GeneralFadExpr.hpp"
#include "Sacado_Fad_DFadTraits.hpp"
#include "Sacado_Fad_DynamicStorage.hpp"

#define FAD_NS Fad
#include "Sacado_Fad_DFad_tmpl.hpp"
#undef FAD_NS

#endif // SACADO_NEW_FAD_DESIGN_IS_DEFAULT

#include "Sacado_Fad_ViewFad.hpp"

#endif // SACADO_FAD_DFAD_HPP
