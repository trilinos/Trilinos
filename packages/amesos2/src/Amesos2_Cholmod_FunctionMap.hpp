// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/**
   \file   Amesos2_Cholmod_FunctionMap.hpp
   \author Kevin Deweese <kdewees@sandia.gov>
   \date   Tue Aug 6 12:53:10 MDT 2013

   \brief  Template for providing a mechanism to map function calls to the
           correct Solver function based on the scalar type of Matrices and
           MultiVectors
*/

#ifndef AMESOS2_CHOLMOD_FUNCTIONMAP_HPP
#define AMESOS2_CHOLMOD_FUNCTIONMAP_HPP

#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include "Amesos2_FunctionMap.hpp"
#include "Amesos2_Cholmod_TypeMap.hpp"

#include "cholmod.h"

namespace Amesos2 {
  
  template <>
  struct FunctionMap<Cholmod,double>
  {
    
    static void cholmod_init_sparse(size_t nrow, size_t ncol, size_t nzmax,
				    int sorted, void *p, void *x, void *i,
				    cholmod_sparse *sparse, int cholmod_itype)
    {
      sparse->nrow = nrow;
      sparse->ncol = ncol;
      sparse->nzmax = nzmax;
      sparse->stype = 1;
      sparse->itype = cholmod_itype;
      sparse->sorted = 0;
      sparse->packed = 1;
      sparse->xtype = CHOLMOD_REAL;
      sparse->dtype = CHOLMOD_DOUBLE;
      sparse->x = x;
      sparse->p = p;
      sparse->i = i;
    }

    static void cholmod_init_dense(int nrow, int ncol, int d, void *x,
				   cholmod_dense *dense)
    {
      dense->nrow = nrow;
      dense->ncol = ncol;
      dense->d = d;
      dense->xtype = CHOLMOD_REAL;
      dense->dtype = CHOLMOD_DOUBLE;
      dense->x = x;
      dense->nzmax = 0;
      dense->z = NULL;
    }
  };

  template <>
  struct FunctionMap<Cholmod,float> // Cholmod does not support float yet
  {
    static void cholmod_init_sparse(size_t nrow, size_t ncol, size_t nzmax,
				    int sorted, void *p, void *x, void*i,
				    cholmod_sparse* sparse, int cholmod_itype)
    {
      sparse->nrow = nrow;
      sparse->ncol = ncol;
      sparse->nzmax = nzmax;
      sparse->stype = 1;
      sparse->itype = cholmod_itype;
      sparse->sorted = 0;
      sparse->packed = 1;
      sparse->xtype = CHOLMOD_REAL;
      sparse->dtype = CHOLMOD_SINGLE;
      sparse->x = x;
      sparse->p = p;
      sparse->i = i;
    }


    static void cholmod_init_dense(int nrow, int ncol, int d, void *x,
				   cholmod_dense *dense)
    {
      dense->nrow = nrow;
      dense->ncol = ncol;
      dense->d = d;
      dense->xtype = CHOLMOD_REAL;
      dense->dtype = CHOLMOD_SINGLE;
      dense->x = x;
      dense->nzmax = 0;
      dense->z = NULL;
    }
  };

#ifdef HAVE_TEUCHOS_COMPLEX
  template <>
  struct FunctionMap<Cholmod,Kokkos::complex<double>>
  {

    static void cholmod_init_sparse(size_t nrow, size_t ncol, size_t nzmax,
				    int sorted, void *p, void *x, void *i,
				    cholmod_sparse* sparse, int cholmod_itype)
    {
      sparse->nrow = nrow;
      sparse->ncol = ncol;
      sparse->nzmax = nzmax;
      sparse->stype = 1;
      sparse->itype = cholmod_itype;
      sparse->sorted = 0;
      sparse->packed = 1;
      sparse->xtype = CHOLMOD_COMPLEX;
      sparse->dtype = CHOLMOD_DOUBLE;
      sparse->x = x;
      sparse->p = p;
      sparse->i = i;

    }
  
    static void cholmod_init_dense(int nrow, int ncol, int d, void *x,
				   cholmod_dense *dense)
    {
      dense->nrow = nrow;
      dense->ncol = ncol;
      dense->d = d;
      dense->xtype = CHOLMOD_COMPLEX;
      dense->dtype = CHOLMOD_DOUBLE;
      dense->x = x;
      dense->nzmax = 0;
      dense->z = NULL;
    }
  };

  template <>
  struct FunctionMap<Cholmod,Kokkos::complex<float>>
  {

    static void cholmod_init_sparse(size_t nrow, size_t ncol, size_t nzmax,
            int sorted, void *p, void *x, void *i,
            cholmod_sparse* sparse, int cholmod_itype)
    {
      sparse->nrow = nrow;
      sparse->ncol = ncol;
      sparse->nzmax = nzmax;
      sparse->stype = 1;
      sparse->itype = cholmod_itype;
      sparse->sorted = 0;
      sparse->packed = 1;
      sparse->xtype = CHOLMOD_COMPLEX;
      sparse->dtype = CHOLMOD_SINGLE;
      sparse->x = x;
      sparse->p = p;
      sparse->i = i;

    }

    static void cholmod_init_dense(int nrow, int ncol, int d, void *x,
           cholmod_dense *dense)
    {
      dense->nrow = nrow;
      dense->ncol = ncol;
      dense->d = d;
      dense->xtype = CHOLMOD_COMPLEX;
      dense->dtype = CHOLMOD_SINGLE;
      dense->x = x;
      dense->nzmax = 0;
      dense->z = NULL;
    }
  };
#endif

} // end namespace Amesos2

#endif  // AMESOS2_CHOLMOD_FUNCTIONMAP_HPP
