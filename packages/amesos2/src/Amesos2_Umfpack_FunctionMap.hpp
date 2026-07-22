// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef AMESOS2_UMFPACK_FUNCTIONMAP_HPP
#define AMESOS2_UMFPACK_FUNCTIONMAP_HPP

#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include "Amesos2_FunctionMap.hpp"
#include "Amesos2_Umfpack_TypeMap.hpp"

extern "C"
{
  #include "umfpack.h"
}

namespace Amesos2 {

  /* ==================== Specializations ====================
   *
   * \cond Umfpack_function_specializations
   */

  template <>
  struct FunctionMap<Umfpack,double>
  {
    typedef TypeMap<Umfpack,double> type_map;

    /**
     * \brief Binds to the appropriate Umfpack solver driver based on data type
     */

    static int umfpack_solve(
    int sys,
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ],
    double X [ ],
    const double B [ ],
    void *Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO])
    {
      return ::umfpack_di_solve(sys, Ap, Ai, Ax, X, B, Numeric, Control, Info);
    }

    static int umfpack_numeric(
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ],
    void *Symbolic,
    void **Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO])
    {
      return ::umfpack_di_numeric(Ap, Ai, Ax, Symbolic, Numeric, Control, Info);
    }

    static int umfpack_symbolic(
    int n_row,
    int n_col,
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO])
    {
      return ::umfpack_di_symbolic(n_row, n_col, Ap, Ai, Ax, Symbolic, Control, Info);
    }

    static void umfpack_defaults(
    double Control [UMFPACK_CONTROL])
    {
      ::umfpack_di_defaults(Control);
    }

    static void umfpack_free_numeric(void **Numeric)
    {
      return ::umfpack_di_free_numeric(Numeric);
    }

    static void umfpack_free_symbolic(void **Symbolic)
    {
      return ::umfpack_di_free_symbolic(Symbolic);
    }
  };


#ifdef HAVE_TEUCHOS_COMPLEX

  template <>
  struct FunctionMap<Umfpack,Kokkos::complex<double>>
  {
    typedef TypeMap<Umfpack,Kokkos::complex<double>> type_map;

    /**
     * \brief Binds to the appropriate Umfpack solver driver based on data type
     */

    static double * stdComplexToUmfpackDoubleConversion(
    const Kokkos::complex<double> v [ ])
    {
      return (double*)(&v[0]);
    }

    static int umfpack_solve(
    int sys,
    const int Ap [ ],
    const int Ai [ ],
    const Kokkos::complex<double> Ax [ ],
    Kokkos::complex<double> X [ ],
    const Kokkos::complex<double> B [ ],
    void *Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO])
    {
      return ::umfpack_zi_solve(sys, Ap, Ai,
        stdComplexToUmfpackDoubleConversion(Ax), NULL,
        stdComplexToUmfpackDoubleConversion(X), NULL,
        stdComplexToUmfpackDoubleConversion(B), NULL,
        Numeric, Control, Info);
    }

    static int umfpack_numeric(
    const int Ap [ ],
    const int Ai [ ],
    const Kokkos::complex<double> Ax [ ],
    void *Symbolic,
    void **Numeric,
    const double Control[UMFPACK_CONTROL],
    double Info[UMFPACK_INFO])
    {
      return ::umfpack_zi_numeric(Ap, Ai, stdComplexToUmfpackDoubleConversion(Ax), NULL, Symbolic, Numeric, Control, Info);
    }

    static int umfpack_symbolic(
    int n_row,
    int n_col,
    const int Ap [ ],
    const int Ai [ ],
    const Kokkos::complex<double> Ax [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO])
    {
      return ::umfpack_zi_symbolic(n_row, n_col, Ap, Ai, stdComplexToUmfpackDoubleConversion(Ax), NULL, Symbolic, Control, Info);
    }

    static void umfpack_defaults(double Control [UMFPACK_CONTROL])
    {
      ::umfpack_zi_defaults(Control);
    }

    static void umfpack_free_numeric(void **Numeric)
    {
      ::umfpack_zi_free_numeric(Numeric);
    }

    static void umfpack_free_symbolic(void **Symbolic)
    {
      ::umfpack_zi_free_symbolic(Symbolic);
    }
  };

#endif	// HAVE_TEUCHOS_COMPLEX

  /* \endcond Umfpack_function_specializations */


} // end namespace Amesos2

#endif  // AMESOS2_UMFPACK_FUNCTIONMAP_HPP
