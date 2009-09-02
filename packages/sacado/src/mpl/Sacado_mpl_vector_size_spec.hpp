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

#ifndef SCADO_MPL_VECTOR_SIZE_SPEC_HPP
#define SCADO_MPL_VECTOR_SIZE_SPEC_HPP

namespace Sacado {

  namespace mpl {

    template <class T0, class T1, class T2, class T3, class T4, class T5,
	      class T6, class T7, class T8, class T9>
    struct vector_size { 
      static const int sz = 10; 
    };
    template <class T0, class T1, class T2, class T3, class T4, class T5,
	      class T6, class T7, class T8>
    struct vector_size<T0,T1,T2,T3,T4,T5,T6,T7,T8,mpl::none> { 
      static const int sz = 9; 
    };
    template <class T0, class T1, class T2, class T3, class T4, class T5,
	      class T6, class T7>
    struct vector_size<T0,T1,T2,T3,T4,T5,T6,T7,mpl::none,mpl::none> { 
      static const int sz = 8; 
    };
    template <class T0, class T1, class T2, class T3, class T4, class T5,
	      class T6>
    struct vector_size<T0,T1,T2,T3,T4,T5,T6,mpl::none,mpl::none,mpl::none> { 
      static const int sz = 7; 
    };
    template <class T0, class T1, class T2, class T3, class T4, class T5>
    struct vector_size<T0,T1,T2,T3,T4,T5,
		       mpl::none,mpl::none,mpl::none,mpl::none> { 
      static const int sz = 6; 
    };
    template <class T0, class T1, class T2, class T3, class T4>
    struct vector_size<T0,T1,T2,T3,T4,mpl::none,
		       mpl::none,mpl::none,mpl::none,mpl::none> { 
      static const int sz = 5; 
    };
    template <class T0, class T1, class T2, class T3>
    struct vector_size<T0,T1,T2,T3,mpl::none,mpl::none,
		       mpl::none,mpl::none,mpl::none,mpl::none> { 
      static const int sz = 4; 
    };
    template <class T0, class T1, class T2>
    struct vector_size<T0,T1,T2,mpl::none,mpl::none,mpl::none,
		       mpl::none,mpl::none,mpl::none,mpl::none> { 
      static const int sz = 3; 
    };
    template <class T0, class T1>
    struct vector_size<T0,T1,mpl::none,mpl::none,mpl::none,mpl::none,
		       mpl::none,mpl::none,mpl::none,mpl::none> { 
      static const int sz = 2; 
    };
    template <class T0>
    struct vector_size<T0,mpl::none,mpl::none,mpl::none,mpl::none,mpl::none,
		       mpl::none,mpl::none,mpl::none,mpl::none> { 
      static const int sz = 1; 
    };
    template <>
    struct vector_size<mpl::none,mpl::none,mpl::none,mpl::none,mpl::none,
		       mpl::none,mpl::none,mpl::none,mpl::none,mpl::none> { 
      static const int sz = 0; 
    };

  }

}

#endif // SCADO_MPL_VECTOR_SIZE_SPEC_HPP
