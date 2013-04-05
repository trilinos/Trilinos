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

#ifndef SACADO_MPL_DISABLE_IF_HPP
#define SACADO_MPL_DISABLE_IF_HPP

namespace Sacado {

  namespace mpl {

    template <bool, typename T = void>
    struct disable_if_c {};

    template <typename T>
    struct disable_if_c<false, T> {
      typedef T type;
    };

    template <class Cond, typename T = void>
    struct disable_if
      : disable_if_c<Cond::value, T> {};

    template <bool, typename T = void>
    struct lazy_disable_if_c {};

    template <typename T>
    struct lazy_disable_if_c<false, T> {
      typedef typename T::type type;
    };

    template <class Cond, typename T = void>
    struct lazy_disable_if
      : lazy_disable_if_c<Cond::value, T> {};

  }

}

#endif // SACADO_MPL_DISABLE_IF_HPP
