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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SCADO_MPL_VECTOR_PUSH_BACK_SPEC_HPP
#define SCADO_MPL_VECTOR_PUSH_BACK_SPEC_HPP

namespace Sacado {

  namespace mpl {

    template <class Vector, class T, int N>
    struct vector_push_back {};

    template <class Vector, class T> struct vector_push_back<Vector,T,0> :
      mpl::vector<T> {};
    template <class Vector, class T> struct vector_push_back<Vector,T,1> : 
      mpl::vector<typename Vector::t0,T> {};
    template <class Vector, class T> struct vector_push_back<Vector,T,2> : 
      mpl::vector<typename Vector::t0,
		  typename Vector::t1,T> {};
    template <class Vector, class T> struct vector_push_back<Vector,T,3> : 
      mpl::vector<typename Vector::t0,
		  typename Vector::t1,
		  typename Vector::t2,T> {};
    template <class Vector, class T> struct vector_push_back<Vector,T,4> : 
      mpl::vector<typename Vector::t0,
		  typename Vector::t1,
		  typename Vector::t2,
		  typename Vector::t3,T> {};
    template <class Vector, class T> struct vector_push_back<Vector,T,5> : 
      mpl::vector<typename Vector::t0,
		  typename Vector::t1,
		  typename Vector::t2,
		  typename Vector::t3,
		  typename Vector::t4,T> {};
    template <class Vector, class T> struct vector_push_back<Vector,T,6> : 
      mpl::vector<typename Vector::t0,
		  typename Vector::t1,
		  typename Vector::t2,
		  typename Vector::t3,
		  typename Vector::t4,
		  typename Vector::t5,T> {};
    template <class Vector, class T> struct vector_push_back<Vector,T,7> : 
      mpl::vector<typename Vector::t0,
		  typename Vector::t1,
		  typename Vector::t2,
		  typename Vector::t3,
		  typename Vector::t4,
		  typename Vector::t5,
		  typename Vector::t6,T> {};
    template <class Vector, class T> struct vector_push_back<Vector,T,8> : 
      mpl::vector<typename Vector::t0,
		  typename Vector::t1,
		  typename Vector::t2,
		  typename Vector::t3,
		  typename Vector::t4,
		  typename Vector::t5,
		  typename Vector::t6,
		  typename Vector::t7,T> {};
    template <class Vector, class T> struct vector_push_back<Vector,T,9> : 
      mpl::vector<typename Vector::t0,
		  typename Vector::t1,
		  typename Vector::t2,
		  typename Vector::t3,
		  typename Vector::t4,
		  typename Vector::t5,
		  typename Vector::t6,
		  typename Vector::t7,
		  typename Vector::t8,T> {};

  }

}

#endif // SCADO_MPL_VECTOR_PUSH_BACK_SPEC_HPP
