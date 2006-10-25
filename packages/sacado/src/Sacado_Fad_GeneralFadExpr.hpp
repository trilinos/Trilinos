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

#ifndef SACADO_FAD_GENERALFADEXPR_HPP
#define SACADO_FAD_GENERALFADEXPR_HPP

#include "Sacado_Fad_GeneralFad.hpp"

namespace Sacado {

  namespace Fad {

    //! GeneralFad expression template specialization
    /*!
     * This template class represents a simple GeneralFad expression and
     * mixes-in the GeneralFad interface and the expression template
     * interface.
     */
    template <typename T, typename Storage> 
    class Expr< GeneralFad<T,Storage> > : public GeneralFad<T,Storage> {

    public:

      //! Default constructor
      Expr() : 
	GeneralFad<T,Storage>() {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      Expr(const T & x) : 
	GeneralFad<T,Storage>(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      Expr(const int sz, const T & x) : 
	GeneralFad<T,Storage>(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      Expr(const int sz, const int i, const T & x) : 
	GeneralFad<T,Storage>(sz,i,x) {}

      //! Copy constructor
      Expr(const Expr& x) : 
	GeneralFad<T,Storage>(x) {}

      //! Copy constructor from any Expression object
      template <typename S> Expr(const Expr<S>& x) :
	GeneralFad<T,Storage>(x) {}

      //! Destructor
      ~Expr() {}

    }; // class Expr<GeneralFad>

  } // namespace Fad

} // namespace Sacado

#include "Sacado_Fad_Ops.hpp"

#endif // SACADO_FAD_GENERALFADEXPR_HPP
