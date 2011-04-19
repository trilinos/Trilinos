// $Id$ 
// $Source$ 
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

  //! Specialization of %Promote to LogicalSparse types
  template <typename ValT, typename LogT>
  struct Promote< LFad::LogicalSparse<ValT,LogT>, 
		  LFad::LogicalSparse<ValT,LogT> > {
    typedef LFad::LogicalSparse<ValT,LogT> type;
  };

  //! Specialization of %Promote to LogicalSparse types
  template <typename ValT, typename LogT, typename R>
  struct Promote< LFad::LogicalSparse<ValT,LogT>, R > {
    typedef typename ValueType< LFad::LogicalSparse<ValT,LogT> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef LFad::LogicalSparse<value_type,LogT> type;
  };

  //! Specialization of %Promote to LogicalSparse types
  template <typename L, typename ValT, typename LogT>
  struct Promote< L, LFad::LogicalSparse<ValT, LogT> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< LFad::LogicalSparse<ValT,LogT> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef LFad::LogicalSparse<value_type,LogT> type;
  };

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

} // namespace Sacado

#endif // SACADO_LFAD_LOGICALSPARSETRAITS_HPP
