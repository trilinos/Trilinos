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
				    cholmod_sparse *sparse)
    {
      sparse->nrow = nrow;
      sparse->ncol = ncol;
      sparse->nzmax = nzmax;
      sparse->stype = 1;
      sparse->itype = CHOLMOD_LONG;
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
				    cholmod_sparse* sparse)
    {
      sparse->nrow = nrow;
      sparse->ncol = ncol;
      sparse->nzmax = nzmax;
      sparse->stype = 1;
      sparse->itype = CHOLMOD_LONG;
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
				    cholmod_sparse* sparse)
    {
      sparse->nrow = nrow;
      sparse->ncol = ncol;
      sparse->nzmax = nzmax;
      sparse->stype = 1;
      sparse->itype = CHOLMOD_LONG;
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
            cholmod_sparse* sparse)
    {
      sparse->nrow = nrow;
      sparse->ncol = ncol;
      sparse->nzmax = nzmax;
      sparse->stype = 1;
      sparse->itype = CHOLMOD_LONG;
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
