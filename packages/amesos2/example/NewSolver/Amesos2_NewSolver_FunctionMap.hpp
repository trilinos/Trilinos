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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/**
   \file   Amesos2_NewSolver_FunctionMap.hpp
   \author John Doe <jd@sandia.gov>
   \date   Fri Jul  9 14:48:16 CDT 2010
   
   \brief  Template for providing a mechanism to map function calls to the
           correct Solver function based on the scalar type of Matrices and
           MultiVectors
*/

#ifndef AMESOS2_NEWSOLVER_FUNCTIONMAP_HPP
#define AMESOS2_NEWSOLVER_FUNCTIONMAP_HPP

#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include "Amesos2_FunctionMap.hpp"
#include "Amesos2_MatrixAdapter.hpp"
#include "Amesos2_NewSolver_TypeMap.hpp"


namespace Amesos2 {

  /*
   * External declarations of the NewSolver functions in namespace New_Solver
   */
  namespace New_Solver {

    extern "C" {

      /*
       * Include solver headers for different data types in different namespaces.
       */
      namespace S {
	// include here single-precision real definition headers
      }

      namespace D {
	// include here double-precision real definition headers
      }

      namespace C {
	// include here single-precision complex definition headers
      }

      namespace Z {
	// include here double-precision complex definition headers
      }

    } // end extern "C"

  } // end namespace New_Solver


  /** 
   * Helper class which passes on function calls to the appropriate
   * NewSolver function based on the type of its scalar template
   * argument.
   *
   * Some solver libraries have solver and matrix builder functions
   * defined based on data type.  One function for complex, one for
   * double precision complex, another for \c float , and yet another
   * for \c double.  To work elegantly with the Amesos::NewSolver
   * interface we want to be able to perform a single function call
   * which is appropriate for the scalar type of the Matrix and
   * MultiVectors that we are working with.  The \c FunctionMap class
   * provides that capability.
   *
   * The class template is specialized for each data type that
   * NewSolver supports.
   *
   * Provide a specialization for the template class which binds each
   * function to the appropriate NewSolver package function call for
   * the given data type.
   * 
   */
  template <>
  struct FunctionMap<NewSolver,float>
  {
    static void solve( /* args */ )
    {
      New_Solver::S::single_solve( /* args */ );
    }

    static void symbolic_fact( /* args */ )
    {
      New_Solver::S::single_symbolic_fact( /* args */ );
    }

    static void numeric_fact( /* args */ )
    {
      New_Solver::S::single_numeric_fact( /* args */ );
    }

    static void create_CompCol_Matrix( /* args */ )
    {
      New_Solver::S::single_create_CompCol_Matrix( /* args */ );
    }

    static void create_CompRow_Matrix( /* args */ )
    {
      New_Solver::S::single_create_CompRow_Matrix( /* args */ );
    }

    static void create_Dense_Matrix( /* args */ )
    {
      New_Solver::S::single_create_Dense_Matrix( /* args */ );
    }
  };


  template <>
  struct FunctionMap<NewSolver,double>
  {
    static void solve( /* args */ )
    {
      New_Solver::D::double_solve( /* args */ );
    }

    static void symbolic_fact( /* args */ )
    {
      New_Solver::D::double_symbolic_fact( /* args */ );
    }

    static void numeric_fact( /* args */ )
    {
      New_Solver::D::double_numeric_fact( /* args */ );
    }

    static void create_CompCol_Matrix( /* args */ )
    {
      New_Solver::D::double_create_CompCol_Matrix( /* args */ );
    }

    static void create_CompRow_Matrix( /* args */ )
    {
      New_Solver::D::double_create_CompRow_Matrix( /* args */ );
    }

    static void create_Dense_Matrix( /* args */ )
    {
      New_Solver::D::double_create_Dense_Matrix( /* args */ );
    }
  };


  /* The specializations for Teuchos::as<> for New_Solver::complex and
   * New_Solver::doublecomplex are provided in Amesos2_NewSolver_Type.hpp
   */
  template <>
  struct FunctionMap<NewSolver,std::complex<float> >
  {
    static void solve( /* args */ )
    {
      New_Solver::C::complex_solve( /* args */ );
    }

    static void symbolic_fact( /* args */ )
    {
      New_Solver::C::complex_symbolic_fact( /* args */ );
    }

    static void numeric_fact( /* args */ )
    {
      New_Solver::C::complex_numeric_fact( /* args */ );
    }

    static void create_CompCol_Matrix( /* args */ )
    {
      New_Solver::C::complex_create_CompCol_Matrix( /* args */ );
    }

    static void create_CompRow_Matrix( /* args */ )
    {
      New_Solver::C::complex_create_CompRow_Matrix( /* args */ );
    }

    static void create_Dense_Matrix( /* args */ )
    {
      New_Solver::C::complex_create_Dense_Matrix( /* args */ );
    }
  };


  template <>
  struct FunctionMap<NewSolver,std::complex<double> >
  {
    static void solve( /* args */ )
    {
      New_Solver::Z::doublecomplex_solve( /* args */ );
    }

    static void symbolic_fact( /* args */ )
    {
      New_Solver::Z::doublecomplex_symbolic_fact( /* args */ );
    }

    static void numeric_fact( /* args */ )
    {
      New_Solver::Z::doublecomplex_numeric_fact( /* args */ );
    }

    static void create_CompCol_Matrix( /* args */ )
    {
      New_Solver::Z::doublecomplex_create_CompCol_Matrix( /* args */ );
    }

    static void create_CompRow_Matrix( /* args */ )
    {
      New_Solver::Z::doublecomplex_create_CompRow_Matrix( /* args */ );
    }

    static void create_Dense_Matrix( /* args */ )
    {
      New_Solver::Z::doublecomplex_create_Dense_Matrix( /* args */ );
    }
  };

  /*
   * etc. for as many types as need supporting
   */

} // end namespace Amesos2

#endif  // AMESOS2_NEWSOLVER_FUNCTIONMAP_HPP
