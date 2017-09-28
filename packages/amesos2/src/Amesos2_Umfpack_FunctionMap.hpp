// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Sivasankaran Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#ifndef AMESOS2_UMFPACK_FUNCTIONMAP_HPP
#define AMESOS2_UMFPACK_FUNCTIONMAP_HPP

#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include "Amesos2_FunctionMap.hpp"
#include "Amesos2_Umfpack_TypeMap.hpp"

namespace UMFPACK {
  #include "umfpack.h"
} // end namespace UMFPACK

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
      return UMFPACK::umfpack_di_solve(sys, Ap, Ai, Ax, X, B, Numeric, Control, Info);
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
      return UMFPACK::umfpack_di_numeric(Ap, Ai, Ax, Symbolic, Numeric, Control, Info);
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
      return UMFPACK::umfpack_di_symbolic(n_row, n_col, Ap, Ai, Ax, Symbolic, Control, Info);
    }

    static void umfpack_defaults(
    double Control [UMFPACK_CONTROL])
    {
      UMFPACK::umfpack_di_defaults(Control);
    }

    static void umfpack_free_numeric(void **Numeric)
    {
      return UMFPACK::umfpack_di_free_numeric(Numeric);
    }

    static void umfpack_free_symbolic(void **Symbolic)
    {
      return UMFPACK::umfpack_di_free_symbolic(Symbolic);
    }
  };


#ifdef HAVE_TEUCHOS_COMPLEX

  template <>
  struct FunctionMap<Umfpack,std::complex<double>>
  {
    typedef TypeMap<Umfpack,std::complex<double>> type_map;

    /**
     * \brief Binds to the appropriate Umfpack solver driver based on data type
     */

    static double * stdComplexToUmfpackDoubleConversion(
    const std::complex<double> v [ ])
    {
      return (double*)(&v[0]);
    }

    static int umfpack_solve(
    int sys,
    const int Ap [ ],
    const int Ai [ ],
    const std::complex<double> Ax [ ],
    std::complex<double> X [ ],
    const std::complex<double> B [ ],
    void *Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO])
    {
      return UMFPACK::umfpack_zi_solve(sys, Ap, Ai,
        stdComplexToUmfpackDoubleConversion(Ax), NULL,
        stdComplexToUmfpackDoubleConversion(X), NULL,
        stdComplexToUmfpackDoubleConversion(B), NULL,
        Numeric, Control, Info);
    }

    static int umfpack_numeric(
    const int Ap [ ],
    const int Ai [ ],
    const std::complex<double> Ax [ ],
    void *Symbolic,
    void **Numeric,
    const double Control[UMFPACK_CONTROL],
    double Info[UMFPACK_INFO])
    {
      return UMFPACK::umfpack_zi_numeric(Ap, Ai, stdComplexToUmfpackDoubleConversion(Ax), NULL, Symbolic, Numeric, Control, Info);
    }

    static int umfpack_symbolic(
    int n_row,
    int n_col,
    const int Ap [ ],
    const int Ai [ ],
    const std::complex<double> Ax [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO])
    {
      return UMFPACK::umfpack_zi_symbolic(n_row, n_col, Ap, Ai, stdComplexToUmfpackDoubleConversion(Ax), NULL, Symbolic, Control, Info);
    }

    static void umfpack_defaults(double Control [UMFPACK_CONTROL])
    {
      UMFPACK::umfpack_zi_defaults(Control);
    }

    static void umfpack_free_numeric(void **Numeric)
    {
      UMFPACK::umfpack_zi_free_numeric(Numeric);
    }

    static void umfpack_free_symbolic(void **Symbolic)
    {
      UMFPACK::umfpack_zi_free_symbolic(Symbolic);
    }
  };

#endif	// HAVE_TEUCHOS_COMPLEX

  /* \endcond Umfpack_function_specializations */


} // end namespace Amesos2

#endif  // AMESOS2_UMFPACK_FUNCTIONMAP_HPP
