// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#ifndef SACADO_CACHEFAD_DFADTRAITS_HPP
#define SACADO_CACHEFAD_DFADTRAITS_HPP

#include "Sacado_ConfigDefs.h"
#include "Sacado_ADTraits.hpp"
#include "Sacado_Promote.hpp"

// Forward declarations
namespace Sacado {
  namespace CacheFad {
    template <typename T> class DFad;
  }
}

namespace Sacado {

  //! Specialization of ADTraits to DFad types
  template <typename T>
  class ADTraits< CacheFad::DFad<T>, CacheFad::DFad<T> > {
  public:

    typedef CacheFad::DFad<T> promote;
  };

  //! Specialization of ADTraits to DFad types
  template <typename L, typename R>
  class ADTraits< CacheFad::DFad<L>, R > {
  public:

    typedef typename CacheFad::DFad<L>::value_type value_type_l;
    typedef typename Promote<R>::value_type value_type_r;
    typedef typename ADTraits<value_type_l,value_type_r>::promote value_type;

    typedef CacheFad::DFad<value_type> promote;
  };

  //! Specialization of ADTraits to DFad types
  template <typename L, typename R>
  class ADTraits< L, CacheFad::DFad<R> > {
  public:

    typedef typename Promote<L>::value_type value_type_l;
    typedef typename CacheFad::DFad<R>::value_type value_type_r;
    typedef typename ADTraits<value_type_l,value_type_r>::promote value_type;

    typedef CacheFad::DFad<value_type> promote;
  };

} // namespace Sacado

#endif // SACADO_CACHEFAD_DFADTRAITS_HPP
