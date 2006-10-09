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

#ifndef SCADO_MPL_VECTOR_AT_SPEC_HPP
#define SCADO_MPL_VECTOR_AT_SPEC_HPP

namespace Sacado {

  namespace mpl {

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
    template <class Vector>
    struct vector_at<Vector,3> {
      typedef typename Vector::t3 type;
    };
    template <class Vector>
    struct vector_at<Vector,4> {
      typedef typename Vector::t4 type;
    };
    template <class Vector>
    struct vector_at<Vector,5> {
      typedef typename Vector::t5 type;
    };
    template <class Vector>
    struct vector_at<Vector,6> {
      typedef typename Vector::t6 type;
    };
    template <class Vector>
    struct vector_at<Vector,7> {
      typedef typename Vector::t7 type;
    };
    template <class Vector>
    struct vector_at<Vector,8> {
      typedef typename Vector::t8 type;
    };
    template <class Vector>
    struct vector_at<Vector,9> {
      typedef typename Vector::t9 type;
    };

  }

}

#endif // SCADO_MPL_VECTOR_AT_SPEC_HPP
