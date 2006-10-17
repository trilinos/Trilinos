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

#ifndef SACADO_MPL_VECTOR_HPP
#define SACADO_MPL_VECTOR_HPP

#include "Sacado_mpl_none.hpp"
#include "Sacado_mpl_size.hpp"
#include "Sacado_mpl_begin.hpp"
#include "Sacado_mpl_end.hpp"
#include "Sacado_mpl_next.hpp"
#include "Sacado_mpl_push_back.hpp"
#include "Sacado_mpl_at.hpp"
#include "Sacado_mpl_deref.hpp"

// Specializations needed for various mpl operations
#include "Sacado_mpl_vector_size_spec.hpp"
#include "Sacado_mpl_vector_at_spec.hpp"

namespace Sacado {

  namespace mpl {

    // vector tag for mpl operations
    struct vector_tag {};

    // vector
    template <class T0 = mpl::none, 
	      class T1 = mpl::none, 
	      class T2 = mpl::none, 
	      class T3 = mpl::none, 
	      class T4 = mpl::none, 
	      class T5 = mpl::none, 
	      class T6 = mpl::none, 
	      class T7 = mpl::none, 
	      class T8 = mpl::none, 
	      class T9 = mpl::none>
    struct vector : vector_size<T0,T1,T2,T3,T4,T5,T6,T7,T8,T9> {
      typedef vector_tag tag;
      typedef vector type;
      typedef T0 t0;
      typedef T1 t1;
      typedef T2 t2;
      typedef T3 t3;
      typedef T4 t4;
      typedef T5 t5;
      typedef T6 t6;
      typedef T7 t7;
      typedef T8 t8;
      typedef T9 t9;
    };
    
    // iterator
    template <class Vector, int Pos>
    struct vector_iterator {
      static const int value = Pos;
    };

    // size
    template <> 
    struct size_impl<vector_tag> {
      template <class Vector>
      struct apply {
	static const int value = Vector::sz;
      };
    };

    // begin
    template <> 
    struct begin_impl<vector_tag> {
      template <class Vector>
      struct apply {
	typedef vector_iterator<Vector,0> type;
      };
    };

    // end
    template <> 
    struct end_impl<vector_tag> {
      template <class Vector>
      struct apply {
	typedef vector_iterator<Vector,Vector::sz> type;
      };
    };

    // next
    template <class Vector, int Pos>
    struct next< vector_iterator<Vector,Pos> > {
      typedef vector_iterator<Vector,Pos+1> type;
    };

    

    // at
    template <int Pos> 
    struct at_impl<vector_tag, Pos> {
      template <class Vector>
      struct apply : vector_at<Vector,Pos> {};
    };

    // deref
    template <class Vector, int Pos>
    struct deref< vector_iterator<Vector,Pos> > : mpl::at<Vector,Pos> {};

  }

}

#include "Sacado_mpl_vector_push_back_spec.hpp"

namespace Sacado {
  namespace mpl {

    // push_back
    template <> 
    struct push_back_impl<vector_tag> {
      template <class Vector, class T>
      struct apply : vector_push_back<Vector, T, mpl::size<Vector>::value> {};
    };
    
  }
}

#endif
