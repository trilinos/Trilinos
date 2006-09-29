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

#ifndef SCADO_MPL_VECTOR_HPP
#define SCADO_MPL_VECTOR_HPP

namespace Sacado {

  namespace mpl {

    struct none {}; // tag type to denote no element
    struct vector_tag {};

    template <class T0, class T1, class T2>
    struct vector_size { 
      static const int sz = 3; 
    };
    template <class T0, class T1>
    struct vector_size<T0,T1,none> { 
      static const int sz = 2; 
    };
    template <class T0>
    struct vector_size<T0,none,none> { 
      static const int sz = 1; 
    };
    template <>
    struct vector_size<none,none,none> { 
      static const int sz = 0; 
    };

    template <class T0 = none, 
	      class T1 = none, 
	      class T2 = none>
    struct vector : vector_size<T0,T1,T2> {
      typedef vector_tag tag;
      typedef vector type;
      typedef T0 t0;
      typedef T1 t1;
      typedef T2 t2;
    };

    template <class Vector, int Pos> struct vector_at {};
    template <class Vector>
    struct vector_at<Vector,0> {
      typedef typename Vector::t0 type;
    };
    template <class Vector>
    struct vector_at<Vector,1> {
      typedef typename Vector::t1 type;
    };
    template <class Vector>
    struct vector_at<Vector,2> {
      typedef typename Vector::t2 type;
    };

    template <class Vector, int Pos>
    struct vector_iterator {
      static const int value = Pos;
    };

    template <class T> struct next {};
    template <class Vector, int Pos>
    struct next< vector_iterator<Vector,Pos> > {
      typedef vector_iterator<Vector,Pos+1> type;
    };

    template <class T> struct begin_impl {};

    template <> 
    struct begin_impl<vector_tag> {
      template <class Vector>
      struct apply {
	typedef vector_iterator<Vector,0> type;
      };
    };

    template <class T>
    struct begin : 
      begin_impl<typename T::tag>:: template apply<T> {};

    template <class T> struct end_impl {};

    template <> 
    struct end_impl<vector_tag> {
      template <class Vector>
      struct apply {
	typedef vector_iterator<Vector,Vector::sz> type;
      };
    };

    template <class T>
    struct end : 
      end_impl<typename T::tag>:: template apply<T> {};

    template <class T> struct size_impl {};

    template <> 
    struct size_impl<vector_tag> {
      template <class Vector>
      struct apply {
	static const int value = Vector::sz;
      };
    };

    template <class T>
    struct size : 
      size_impl<typename T::tag>:: template apply<T> {};

    template <class T, int Pos> struct at_impl {};

    template <int Pos> 
    struct at_impl<vector_tag, Pos> {
      template <class Vector>
      struct apply : vector_at<Vector,Pos> {};
    };

    template <class T, int Pos>
    struct at : 
      at_impl<typename T::tag,Pos>:: template apply<T> {};

    template <class T> struct deref {};
    template <class Vector, int Pos>
    struct deref< vector_iterator<Vector,Pos> > : at<Vector,Pos> {};

    template <class T1, class T2>
    struct is_same {
      static const bool value = false;
    };

    template <class T>
    struct is_same<T,T> {
      static const bool value = true;
    };

    template <bool cond, class T1, class T2> struct mpl_if {};
    template <class T1, class T2> struct mpl_if<true,T1,T2> : T1 {};
    template <class T1, class T2> struct mpl_if<false,T1,T2> : T2 {};

    template <class T1, class T2, class T>
    struct find_it {};

    template <class Vector, int Pos1, int Pos2, class T>
    struct find_it<vector_iterator<Vector,Pos1>,
		   vector_iterator<Vector,Pos2>,
		   T> {
      static const int value = 
      mpl_if< is_same< typename deref< vector_iterator<Vector,Pos1> >::type,T>::value,
	      vector_iterator<Vector,Pos1>,
	      find_it<typename next< vector_iterator<Vector,Pos1> >::type,
		      vector_iterator<Vector,Pos2>,
		      T> >::value;
    };

    template <class Vector, int Pos, class T>
    struct find_it<vector_iterator<Vector,Pos>,
		   vector_iterator<Vector,Pos>,
		   T> {
      static const int value = vector_iterator<Vector,Pos>::value;
    };

    template <class Seq, class T>
    struct find : 
      find_it<typename begin<Seq>::type, typename end<Seq>::type, T> {};

    template <class T1, class T2>
    struct for_each_it {
      template <typename Op>
      for_each_it(const Op& op) {}
    };

    template <class Vector, int Pos1, int Pos2>
    struct for_each_it< vector_iterator<Vector,Pos1>,
		        vector_iterator<Vector,Pos2> > {
      template <typename Op>
      for_each_it(const Op& op) {
	op(typename deref< vector_iterator<Vector,Pos1> >::type());
	for_each_it< typename next< vector_iterator<Vector,Pos1> >::type,
	  vector_iterator<Vector,Pos2> > f(op);
      }
    };

    template <class Vector, int Pos>
    struct for_each_it< vector_iterator<Vector,Pos>,
		        vector_iterator<Vector,Pos> > {
      template <typename Op>
      for_each_it(const Op& op) {}
    };

    template <class Seq>
    struct for_each {
      template <typename Op>
      for_each(const Op& op) {
	for_each_it<typename begin<Seq>::type, typename end<Seq>::type> f(op);
      }
    };

  }

}

#endif
