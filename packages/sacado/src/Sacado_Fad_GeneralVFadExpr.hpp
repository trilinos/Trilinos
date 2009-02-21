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

#ifndef SACADO_FAD_GENERALVFADEXPR_HPP
#define SACADO_FAD_GENERALVFADEXPR_HPP

#include "Sacado_Fad_GeneralVFad.hpp"

namespace Sacado {

  namespace Fad {

    //! GeneralVFad expression template specialization
    /*!
     * This template class represents a simple GeneralVFad expression and
     * mixes-in the GeneralVFad interface and the expression template
     * interface.
     */
    template <typename T, typename Storage> 
    class Expr< GeneralVFad<T,Storage> > : public GeneralVFad<T,Storage> {

    public:

      //! Default constructor
      Expr() : 
	GeneralVFad<T,Storage>() {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      Expr(const T & x) : 
	GeneralVFad<T,Storage>(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      Expr(const int sz, const T & x) : 
	GeneralVFad<T,Storage>(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      Expr(const int sz, const int i, const T & x) : 
	GeneralVFad<T,Storage>(sz,i,x) {}

      //! Constructor with supplied memory
      /*!
       * Initializes value to point to \c x and derivative array to point 
       * to\c dx.  Derivative array is zero'd out if \c zero_out is true.
       */
      Expr(const int sz, T* x, T* dx, bool zero_out = false) : 
	GeneralVFad<T,Storage>(sz, x, dx, zero_out) {}

      //! Constructor with supplied memory and index \c i
      /*!
       * Initializes value to point to \c x and derivative array to point 
       * to\c dx.  Initializes derivative array row \c i of the identity matrix,
       * i.e., sets derivative component \c i to 1 and all other's to zero.
       */
      Expr(const int sz, const int i, T* x, T* dx) : 
	GeneralVFad<T,Storage>(sz, i, x, dx) {}

      //! Copy constructor
      Expr(const Expr& x) : 
	GeneralVFad<T,Storage>(x) {}

      //! Copy constructor from any Expression object
      template <typename S> Expr(const Expr<S>& x) :
	GeneralVFad<T,Storage>(x) {}

      //! Destructor
      ~Expr() {}

    }; // class Expr<GeneralVFad>

  } // namespace Fad

} // namespace Sacado

#include "Sacado_Fad_Ops.hpp"

#endif // SACADO_FAD_GENERALVFADEXPR_HPP
