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

#ifndef SACADO_ETPCE_ORTHOGPOLYTRAITS_HPP
#define SACADO_ETPCE_ORTHOGPOLYTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace ETPCE {
    template <typename T, typename S> class OrthogPoly;
  }
}

namespace Sacado {

  //! Specialization of %Promote to OrthogPoly types
  template <typename T, typename S>
  struct Promote< ETPCE::OrthogPoly<T,S>,
                  ETPCE::OrthogPoly<T,S> > {
    typedef ETPCE::OrthogPoly<T,S> type;
  };
  template <typename T, typename S>
  struct Promote< ETPCE::OrthogPoly<T,S>,
                  typename ETPCE::OrthogPoly<T,S>::value_type > {
    typedef ETPCE::OrthogPoly<T,S> type;
  };
  template <typename T, typename S>
  struct Promote< typename ETPCE::OrthogPoly<T,S>::value_type,
                  ETPCE::OrthogPoly<T,S> > {
    typedef ETPCE::OrthogPoly<T,S> type;
  };
  template <typename T, typename S>
  struct Promote< ETPCE::OrthogPoly<T,S>,
                  typename dummy< T,
                                  typename ETPCE::OrthogPoly<T,S>::scalar_type
                                  >::type > {
    typedef ETPCE::OrthogPoly<T,S> type;
  };
  template <typename T, typename S>
  struct Promote< typename dummy< T,
                                  typename ETPCE::OrthogPoly<T,S>::scalar_type
                                  >::type,
                  ETPCE::OrthogPoly<T,S> > {
    typedef ETPCE::OrthogPoly<T,S> type;
  };

  //! Specialization of %ScalarType to OrthogPoly types
  template <typename T, typename S>
  struct ScalarType< ETPCE::OrthogPoly<T,S> > {
    typedef typename ScalarType<typename ETPCE::OrthogPoly<T,S>::value_type>::type type;
  };

  //! Specialization of %ValueType to OrthogPoly types
  template <typename T, typename S>
  struct ValueType< ETPCE::OrthogPoly<T,S> > {
    typedef typename ETPCE::OrthogPoly<T,S>::value_type type;
  };

  //! Specialization of %IsADType to OrthogPoly types
  template <typename T, typename S>
  struct IsADType< ETPCE::OrthogPoly<T,S> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to OrthogPoly types
  template <typename T, typename S>
  struct IsScalarType< ETPCE::OrthogPoly<T,S> > {
    static const bool value = false;
  };

  //! Specialization of %Value to OrthogPoly types
  template <typename T, typename S>
  struct Value< ETPCE::OrthogPoly<T,S> > {
    typedef typename ValueType< ETPCE::OrthogPoly<T,S> >::type value_type;
    static const value_type& eval(const ETPCE::OrthogPoly<T,S>& x) { 
      return x.val(); }
  };

  //! Specialization of %ScalarValue to OrthogPoly types
  template <typename T, typename S>
  struct ScalarValue< ETPCE::OrthogPoly<T,S> > {
    typedef typename ValueType< ETPCE::OrthogPoly<T,S> >::type value_type;
    typedef typename ScalarType< ETPCE::OrthogPoly<T,S> >::type scalar_type;
    static const scalar_type& eval(const ETPCE::OrthogPoly<T,S>& x) { 
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to OrthogPoly types
  template <typename T, typename S>
  struct StringName< ETPCE::OrthogPoly<T,S> > {
    static std::string eval() { 
      return std::string("Sacado::ETPCE::OrthogPoly< ") + 
	StringName<T>::eval() + " >"; }
  };

  //! Specialization of %IsEqual to OrthogPoly types
  template <typename T, typename S>
  struct IsEqual< ETPCE::OrthogPoly<T,S> > {
    static bool eval(const ETPCE::OrthogPoly<T,S>& x, 
		     const ETPCE::OrthogPoly<T,S>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to OrthogPoly types
  template <typename T, typename S>
  struct IsStaticallySized< ETPCE::OrthogPoly<T,S> > {
    static const bool value = false;
  };

} // namespace Sacado

// Define Teuchos traits classes
#ifdef HAVE_SACADO_TEUCHOS
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Sacado_PCE_ScalarTraitsImp.hpp"
#include "Teuchos_SerializationTraits.hpp"

namespace Teuchos {

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename T, typename S>
  struct PromotionTraits< Sacado::ETPCE::OrthogPoly<T,S>, 
			  Sacado::ETPCE::OrthogPoly<T,S> > {
    typedef typename Sacado::Promote< Sacado::ETPCE::OrthogPoly<T,S>,
				      Sacado::ETPCE::OrthogPoly<T,S> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename T, typename S, typename R>
  struct PromotionTraits< Sacado::ETPCE::OrthogPoly<T,S>, R > {
    typedef typename Sacado::Promote< Sacado::ETPCE::OrthogPoly<T,S>, R >::type 
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename L, typename T, typename S>
  struct PromotionTraits< L, Sacado::ETPCE::OrthogPoly<T,S> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::ETPCE::OrthogPoly<T,S> >::type 
    promote;
  };

  //! Specializtion of %Teuchos::ScalarTraits
  template <typename T, typename S>
  struct ScalarTraits< Sacado::ETPCE::OrthogPoly<T,S> > :
    public Sacado::PCE::ScalarTraitsImp< Sacado::ETPCE::OrthogPoly<T,S> > {};

  //! Specialization of %Teuchos::SerializationTraits
  template <typename Ordinal, typename T, typename S>
  struct SerializationTraits<Ordinal, Sacado::ETPCE::OrthogPoly<T,S> > :
    public Sacado::PCE::SerializationTraitsImp< Ordinal, 
						Sacado::ETPCE::OrthogPoly<T,S> >
  {};

  //! Specialization of %Teuchos::ValueTypeSerializer
  template <typename Ordinal, typename T, typename S>
  struct ValueTypeSerializer<Ordinal, Sacado::ETPCE::OrthogPoly<T,S> > :
    public Sacado::PCE::SerializerImp< Ordinal, 
				       Sacado::ETPCE::OrthogPoly<T,S>,
				       ValueTypeSerializer<Ordinal,T> > 
  {
    typedef Sacado::ETPCE::OrthogPoly<T,S> PCEType;
    typedef ValueTypeSerializer<Ordinal,T> ValueSerializer;
    typedef Sacado::PCE::SerializerImp< Ordinal,PCEType,ValueSerializer> Base;
    typedef typename Base::expansion_type expansion_type;
    ValueTypeSerializer(const Teuchos::RCP<expansion_type>& expansion,
			const Teuchos::RCP<const ValueSerializer>& vs) :
      Base(expansion,vs) {}
  };
}
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_ETPCE_ORTHOGPOLYTRAITS_HPP
