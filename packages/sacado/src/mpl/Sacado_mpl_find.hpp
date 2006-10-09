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

#ifndef SACADO_MPL_FIND_HPP
#define SACADO_MPL_FIND_HPP

#include "Sacado_mpl_none.hpp"
#include "Sacado_mpl_begin.hpp"
#include "Sacado_mpl_end.hpp"
#include "Sacado_mpl_deref.hpp"
#include "Sacado_mpl_next.hpp"
#include "Sacado_mpl_is_same.hpp"
#include "Sacado_mpl_if.hpp"

namespace Sacado {

  namespace mpl {

    template <class Seq, class T>
    class TypeSequenceDoesNotContainType {};

    template <class Seq, 
	      class T,
	      class Iter1 = typename mpl::begin<Seq>::type, 
	      class Iter2 = typename mpl::end<Seq>::type>
    struct find {
      static const int value = 
        mpl::mpl_if< mpl::is_same<typename mpl::deref<Iter1>::type, T>::value, 
		     Iter1,
		     find<Seq, T, typename mpl::next<Iter1>::type, 
			  Iter2> >::value;
    };

    template <class Seq, class T, class Iter1>
    struct find<Seq, T, Iter1, Iter1> {
      static const int value = TypeSequenceDoesNotContainType<Seq,T>::value;
    };

  }

}

#endif // SACADO_MPL_FIND_HPP
