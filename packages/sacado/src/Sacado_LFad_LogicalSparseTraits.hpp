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

#ifndef SACADO_LFAD_LOGICALSPARSETRAITS_HPP
#define SACADO_LFAD_LOGICALSPARSETRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace LFad {
    template <typename T1, typename T2> class LogicalSparse;
  }
}

namespace Sacado {

  namespace LFad {
    template <typename T> class Expr;
  }

  //! Specialization of %Promote to LogicalSparse types
  SACADO_AD_PROMOTE_SPEC2( LFad, LogicalSparse )

  //! Specialization of %ScalarType to LogicalSparse types
  template <typename ValT, typename LogT>
  struct ScalarType< LFad::LogicalSparse<ValT,LogT> > {
    typedef typename ScalarType< typename LFad::LogicalSparse<ValT,LogT>::value_type >::type type;
  };

  //! Specialization of %ValueType to LogicalSparse types
  template <typename ValT, typename LogT>
  struct ValueType< LFad::LogicalSparse<ValT,LogT> > {
    typedef typename LFad::LogicalSparse<ValT,LogT>::value_type type;
  };

  //! Specialization of %IsADType to LogicalSparse types
  template <typename ValT, typename LogT>
  struct IsADType< LFad::LogicalSparse<ValT,LogT> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to LogicalSparse types
  template <typename ValT, typename LogT>
  struct IsScalarType< LFad::LogicalSparse<ValT,LogT> > {
    static const bool value = false;
  };

  //! Specialization of %Value to LogicalSparse types
  template <typename ValT, typename LogT>
  struct Value< LFad::LogicalSparse<ValT,LogT> > {
    typedef typename ValueType< LFad::LogicalSparse<ValT,LogT> >::type value_type;
    static const value_type& eval(const LFad::LogicalSparse<ValT,LogT>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to DFad types
  template <typename ValT, typename LogT>
  struct ScalarValue< LFad::LogicalSparse<ValT,LogT> > {
    typedef typename ValueType< LFad::LogicalSparse<ValT,LogT> >::type value_type;
    typedef typename ScalarType< LFad::LogicalSparse<ValT,LogT> >::type scalar_type;
    static const scalar_type& eval(const LFad::LogicalSparse<ValT,LogT>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to DFad types
  template <typename ValT, typename LogT>
  struct StringName< LFad::LogicalSparse<ValT,LogT> > {
    static std::string eval() {
      return std::string("Sacado::LFad::LoginalSparse< ") +
        StringName<ValT>::eval() + ", " +
        StringName<LogT>::eval() + " >"; }
  };

  //! Specialization of %IsEqual to DFad types
  template <typename ValT, typename LogT>
  struct IsEqual< LFad::LogicalSparse<ValT,LogT> > {
    static bool eval(const LFad::LogicalSparse<ValT,LogT>& x,
                     const LFad::LogicalSparse<ValT,LogT>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to DFad types
  template <typename ValT, typename LogT>
  struct IsStaticallySized< LFad::LogicalSparse<ValT,LogT> > {
    static const bool value = false;
  };

} // namespace Sacado

// Define Teuchos traits classes
#ifdef HAVE_SACADO_TEUCHOS
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Sacado_Fad_ScalarTraitsImp.hpp"
#include "Teuchos_SerializationTraits.hpp"

namespace Teuchos {

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename ValT, typename LogT>
  struct PromotionTraits< Sacado::LFad::LogicalSparse<ValT,LogT>,
                          Sacado::LFad::LogicalSparse<ValT,LogT> > {
    typedef typename Sacado::Promote< Sacado::LFad::LogicalSparse<ValT,LogT>,
                                      Sacado::LFad::LogicalSparse<ValT,LogT> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename ValT, typename LogT, typename R>
  struct PromotionTraits< Sacado::LFad::LogicalSparse<ValT,LogT>, R > {
    typedef typename Sacado::Promote< Sacado::LFad::LogicalSparse<ValT,LogT>, R >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename L, typename ValT, typename LogT>
  struct PromotionTraits< L, Sacado::LFad::LogicalSparse<ValT,LogT> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::LFad::LogicalSparse<ValT,LogT> >::type
    promote;
  };

  //
  // These specialization implementations don't work for LFad because
  // the value type (e.g., double) is different from the dx() type (e.g., bool)
  //

  // //! Specializtion of %Teuchos::ScalarTraits
  // template <typename ValT, typename LogT>
  // struct ScalarTraits< Sacado::LFad::LogicalSparse<ValT,LogT> > :
  //   public Sacado::Fad::ScalarTraitsImp< Sacado::LFad::LogicalSparse<ValT,LogT> >
  // {};

  // //! Specialization of %Teuchos::SerializationTraits
  // template <typename Ordinal, typename ValT, typename LogT>
  // struct SerializationTraits<Ordinal, Sacado::LFad::LogicalSparse<ValT,LogT> > :
  //   public Sacado::Fad::SerializationTraitsImp< Ordinal,
  //                                            Sacado::LFad::LogicalSparse<ValT,LogT> >
  // {};
}
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_LFAD_LOGICALSPARSETRAITS_HPP
