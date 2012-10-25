// $Id$
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_MPL_FOR_EACH_HPP
#define STOKHOS_MPL_FOR_EACH_HPP

#include "Sacado_mpl_begin.hpp"
#include "Sacado_mpl_end.hpp"
#include "Sacado_mpl_next.hpp"
#include "Sacado_mpl_deref.hpp"

#include "KokkosArray_Macros.hpp"

namespace Stokhos {

  namespace mpl {

    template <class Seq, 
	      class node_t, 
	      class Iter1 = typename Sacado::mpl::begin<Seq>::type, 
	      class Iter2 = typename Sacado::mpl::end<Seq>::type>
    struct for_each {
      typedef node_t node_type;
      template <typename Op>
      KOKKOSARRAY_INLINE_DEVICE_FUNCTION
      for_each(const Op& op) {
	op(typename Sacado::mpl::deref<Iter1>::type());
	for_each<Seq, node_type, typename Sacado::mpl::next<Iter1>::type, Iter2> f(op);
      }
    };

    template <class Seq, class node_t, class Iter1>
    struct for_each<Seq, node_t, Iter1, Iter1> {
      template <typename Op>
      KOKKOSARRAY_INLINE_DEVICE_FUNCTION
      for_each(const Op& op) {}
    };

  }

}

#endif // STOKHOS_MPL_FOR_EACH_HPP
