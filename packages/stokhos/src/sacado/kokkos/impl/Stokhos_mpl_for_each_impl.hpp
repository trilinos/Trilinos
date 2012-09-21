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

#if ! defined(KOKKOSARRAY_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOSARRAY_MACRO_DEVICE)                  || \
    ! defined(KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION)

#error "Including <Stokhos_mpl_for_each_impl.hpp> without macros defined"

#else

namespace Stokhos {

  namespace mpl {

    template <class Seq, class Iter1, class Iter2>
    struct for_each<Seq, KOKKOSARRAY_MACRO_DEVICE, Iter1, Iter2> {
      typedef KOKKOSARRAY_MACRO_DEVICE node_type;
      template <typename Op>
      KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
      for_each(const Op& op) {
	op(typename Sacado::mpl::deref<Iter1>::type());
	for_each<Seq, node_type, typename Sacado::mpl::next<Iter1>::type, Iter2> f(op);
      }
    };

    template <class Seq, class Iter1>
    struct for_each<Seq, KOKKOSARRAY_MACRO_DEVICE, Iter1, Iter1> {
      template <typename Op>
      KOKKOSARRAY_MACRO_DEVICE_AND_HOST_FUNCTION
      for_each(const Op& op) {}
    };

  }

}

#endif
