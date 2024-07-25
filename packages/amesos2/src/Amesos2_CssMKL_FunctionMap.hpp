// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/**
   \file   Amesos2_NewSolver_FunctionMap.hpp
   \author Eric Bavier <etbavie@sandia.gov>
   \date   Wed Jul 27 12:53:10 MDT 2011

   \brief  Template for providing a mechanism to map function calls to the
           correct Solver function based on the scalar type of Matrices and
           MultiVectors
*/

#ifndef AMESOS2_CSSMKL_FUNCTIONMAP_HPP
#define AMESOS2_CSSMKL_FUNCTIONMAP_HPP

#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include "Amesos2_FunctionMap.hpp"
#include "Amesos2_CssMKL_TypeMap.hpp"


namespace Amesos2 {

  namespace PMKL {
    #ifdef __MKL_PARDISO_H
      #undef __MKL_PARDISO_H
    #endif
    #include "mkl_pardiso.h"
    #ifdef __MKL_CLUSTER_SPARSE_SOLVER_H
      #undef __MKL_CLUSTER_SPARSE_SOLVER_H
    #endif
    #include "mkl_cluster_sparse_solver.h"
  }

  /** \internal
   * 
   * For Pardiso we bind to the library functions based on the local
   * ordinal type.  If the local ordinal type is bigger than int, then
   * we use pardiso_64 instead.  The void* arrays are interpreted by
   * the function based on the value of mtype and iparm(28) as
   * single/double and complex/real.
   */
  template <>
  struct FunctionMap<CssMKL,PMKL::_INTEGER_t>
  {
    static void cluster_sparse_solver( void* pt,
			 PMKL::_INTEGER_t* maxfct, PMKL::_INTEGER_t* mnum,
			 PMKL::_INTEGER_t* mtype , PMKL::_INTEGER_t* phase,
			 PMKL::_INTEGER_t* n     , void* a, PMKL::_INTEGER_t* ia,
			 PMKL::_INTEGER_t* ja    , PMKL::_INTEGER_t* perm,
			 PMKL::_INTEGER_t* nrhs  , PMKL::_INTEGER_t* iparm,
			 PMKL::_INTEGER_t* msglvl, void* b, void* x,
			 const MPI_Fint  * comm  , PMKL::_INTEGER_t* error)
    {
      PMKL::cluster_sparse_solver(pt, maxfct, mnum, mtype, phase, n, a, ia, ja,
		                  perm, nrhs, iparm, msglvl, b, x, comm, error);
    }
  };


  template <>
  struct FunctionMap<CssMKL,long long int>
  {
    static void cluster_sparse_solver( void* pt,
			 long long int*  maxfct, long long int* mnum,
			 long long int*  mtype , long long int* phase,
			 long long int*  n     , void* a, long long int* ia,
			 long long int*  ja    , long long int* perm,
			 long long int*  nrhs  , long long int* iparm,
			 long long int*  msglvl, void* b, void* x,
			 const MPI_Fint* comm  , long long int* error)
    {
      PMKL::cluster_sparse_solver_64(pt, maxfct, mnum, mtype, phase, n, a, ia, ja,
		                     perm, nrhs, iparm, msglvl, b, x, comm, error);
    }
  };
} // end namespace Amesos2

#endif  // AMESOS2_NEWSOLVER_FUNCTIONMAP_HPP
