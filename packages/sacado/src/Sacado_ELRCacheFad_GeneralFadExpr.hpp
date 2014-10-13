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

#ifndef SACADO_ELRCACHEFAD_GENERALFADEXPR_HPP
#define SACADO_ELRCACHEFAD_GENERALFADEXPR_HPP

#include "Sacado_ELRCacheFad_GeneralFad.hpp"

namespace Sacado {

  namespace ELRCacheFad {

    //! GeneralFad expression template specialization
    /*!
     * This template class represents a simple GeneralFad expression and
     * mixes-in the GeneralFad interface and the expression template
     * interface using the caching expression-level-revese-mode
     * expression templates.
     */
    template <typename T, typename Storage>
    class Expr< GeneralFad<T,Storage> > : public GeneralFad<T,Storage> {

    public:

      typedef typename GeneralFad<T,Storage>::value_type value_type;
      typedef typename GeneralFad<T,Storage>::scalar_type scalar_type;

      //! Typename of base-expressions
      typedef GeneralFad<T,Storage> base_expr_type;

      //! Number of arguments
      static const int num_args = 1;

      //! Is expression linear
      static const bool is_linear = true;

      //! Default constructor
      KOKKOS_INLINE_FUNCTION
      Expr() :
        GeneralFad<T,Storage>() {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      KOKKOS_INLINE_FUNCTION
      Expr(const T & x) :
        GeneralFad<T,Storage>(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      KOKKOS_INLINE_FUNCTION
      Expr(const int sz, const T & x) :
        GeneralFad<T,Storage>(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      KOKKOS_INLINE_FUNCTION
      Expr(const int sz, const int i, const T & x) :
        GeneralFad<T,Storage>(sz,i,x) {}

      //! Constructor with supplied storage \c s
      KOKKOS_INLINE_FUNCTION
      Expr(const Storage& s) :
        GeneralFad<T,Storage>(s) {}

      //! Copy constructor
      KOKKOS_INLINE_FUNCTION
      Expr(const Expr& x) :
        GeneralFad<T,Storage>(static_cast<const GeneralFad<T,Storage>&>(x)) {}

      //! Copy constructor from any Expression object
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Expr(const Expr<S>& x) :
        GeneralFad<T,Storage>(x) {}

      //! Destructor
      KOKKOS_INLINE_FUNCTION
      ~Expr() {}

      //! Return partials w.r.t. arguments
      KOKKOS_INLINE_FUNCTION
      void computePartials(const T& bar, value_type partials[]) const {
        partials[0] = bar;
      }

      //! Return tangent component \c i of arguments
      KOKKOS_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        if (i<this->size())
          dots[0] = this->fastAccessDx(i);
        else
          dots[0] = 0.0;
      }

      //! Return whether argument is active
      template <int Arg>
      KOKKOS_INLINE_FUNCTION
      bool isActive() const { return this->size() > 0; }

      //! Return whether expression is linear
      KOKKOS_INLINE_FUNCTION
      bool isLinear() const { return true; }

      //! Return tangent component \c i of argument \c Arg
      template <int Arg>
      KOKKOS_INLINE_FUNCTION
      T getTangent(int i) const { return this->fastAccessDx(i); }

      KOKKOS_INLINE_FUNCTION
      const base_expr_type& getArg(int j) const { return *this; }

    }; // class Expr<GeneralFad>

    //! Specialization of %ExprPromote to GeneralFad types
    template <typename T, typename S>
    struct ExprPromote< GeneralFad<T,S>,
                        typename GeneralFad<T,S>::value_type > {
      typedef GeneralFad<T,S> type;
    };

    //! Specialization of %ExprPromote to GeneralFad types
    template <typename T, typename S>
    struct ExprPromote< typename GeneralFad<T,S>::value_type,
                        GeneralFad<T,S> > {
      typedef GeneralFad<T,S> type;
    };

    //! Specialization of %ExprPromote to GeneralFad types
    template <typename T, typename S>
    struct ExprPromote< GeneralFad<T,S>,
                        typename dummy<typename GeneralFad<T,S>::value_type,
                                       typename GeneralFad<T,S>::scalar_type
                                       >::type > {
      typedef GeneralFad<T,S> type;
    };

    //! Specialization of %ExprPromote to GeneralFad types
    template <typename T, typename S>
    struct ExprPromote< typename dummy<typename GeneralFad<T,S>::value_type,
                                       typename GeneralFad<T,S>::scalar_type
                                       >::type,
                        GeneralFad<T,S> > {
      typedef GeneralFad<T,S> type;
    };

  } // namespace ELRCacheFad

} // namespace Sacado

#include "Sacado_ELRCacheFad_Ops.hpp"

#endif // SACADO_ELRCACHEFAD_GENERALFADEXPR_HPP
