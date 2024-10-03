// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_KL_ONE_D_EXPONENTIAL_EIGENPAIR_HPP
#define STOKHOS_KL_ONE_D_EXPONENTIAL_EIGENPAIR_HPP

#include <iostream>
#include <cmath>

#include "Kokkos_Core.hpp"

namespace Stokhos {

  //! Namespace for analytic %KL expansions
  namespace KL {

    //! Container for one-dimensional eigenfunction and eigenvalue
    template <typename eigen_function_type>
    struct OneDEigenPair {
      typedef typename eigen_function_type::value_type value_type;
      eigen_function_type eig_func;
      value_type eig_val;
    }; // struct OneDEigenPair

    //! One-dimensional eigenfunction for exponential covariance function
    /*!
     * Represents an eigenfunction of the form \f$A \sin(\omega (x-(b+a)/2))\f$
     * or  \f$A \cos(\omega (x-(b+a)/2))\f$ over the domain \f$[a,b]\f$ where
     * \f[
     *   A = \frac{1}{\sqrt{\frac{b-a}{2} \pm \frac{\sin(\omega(b-a)}{2\omega}}}
     * \f]
     * for \f$\cos\f$, \f$\sin\f$ respectively.
     */
    template <typename Value>
    class ExponentialOneDEigenFunction {
    public:

      typedef Value value_type;

      //! Enum identifying the type of eigenfunction
      enum TYPE {
        SIN, ///< A*sin(omega*(x-b))
        COS  ///< A*cos(omega*(x-b))
      };

      //! Default Constructor
      KOKKOS_INLINE_FUNCTION
      ExponentialOneDEigenFunction() :
        type(SIN), a(0), b(0), A(0), omega(0), dim_name(0) {}

      //! Constructor
      KOKKOS_INLINE_FUNCTION
      ExponentialOneDEigenFunction(TYPE type_, const value_type& a_,
                                   const value_type& b_,
                                   const value_type& omega_,
                                   const int dim_name_) :
        type(type_), a((b_-a_)/2.0), b((b_+a_)/2.0), omega(omega_),
        dim_name(dim_name_) {
        if (type == SIN)
          A = 1.0/std::sqrt(a - std::sin(2.*omega*a)/(2.*omega));
        else
          A = 1.0/std::sqrt(a + std::sin(2.*omega*a)/(2.*omega));
      }

      //! Destructor
      KOKKOS_INLINE_FUNCTION
      ~ExponentialOneDEigenFunction() {}

      //! Evaluate eigenfunction
      template <typename point_type>
      KOKKOS_INLINE_FUNCTION
      point_type
      evaluate(const point_type& x) const {
        if (type == SIN)
          return A*sin(omega*(x-b));
        return A*cos(omega*(x-b));
      }

      //! Print eigenfunction
      void print(std::ostream& os) const {
        os << A << " * ";
        if (type == SIN)
          os << "sin(";
        else
          os << "cos(";
        os << omega << " * (x_" << dim_name << " - " << b << "))";
      }

      //! Return type
      KOKKOS_INLINE_FUNCTION
      TYPE getType() const { return type; }

      //! Return frequency
      KOKKOS_INLINE_FUNCTION
      value_type getFrequency() const { return omega; }

      //! Return multiplier
      KOKKOS_INLINE_FUNCTION
      value_type getMultiplier() const { return A; }

      //! Get shift
      KOKKOS_INLINE_FUNCTION
      value_type getShift() const { return b; }

    protected:

      //! Type of eigenfunction (sin or cos)
      TYPE type;

      //! Domain length
      value_type a;

      //! Domain center
      value_type b;

      //! Multiplier for eigenfunction
      value_type A;

      //! Frequency of eigenfunction
      value_type omega;

      //! Dimesion name (e.g., x_1) for printing eigenfunction
      int dim_name;
    };

  } // namespace KL

} // namespace Stokhos

#endif // STOKHOS_KL_ONE_D_EXPONENTIALEIGENPAIR_HPP
