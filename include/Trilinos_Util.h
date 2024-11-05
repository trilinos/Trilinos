// @HEADER
// ***********************************************************************
//
//                 TriUtils: Trilinos Utilities Package
//                 Copyright (2011) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// @HEADER

#ifndef _TRILINOS_UTIL_H_
#define _TRILINOS_UTIL_H_

#if defined(Triutils_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Triutils package is deprecated"
#endif
#endif

#define Trilinos_Util_max(x,y) (( x > y ) ? x : y)     /* max function  */
#define Trilinos_Util_min(x,y) (( x < y ) ? x : y)     /* min function */

#ifndef TRILINOS_NO_CONFIG_H
/*
 * The macros PACKAGE, PACKAGE_NAME, etc, get defined for each package and need to
 * be undef'd here to avoid warnings when this file is included from another package.
 * Fix from KL
 */
#ifdef PACKAGE
#undef PACKAGE
#endif

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif

#ifdef VERSION
#undef VERSION
#endif

#include "Triutils_config.h"

#ifdef HAVE_DEBUG
#ifndef DEBUG
#define DEBUG
#endif
#endif

#ifndef JANUS_STLPORT
#ifdef HAVE_CSTDLIB
#include <cstdlib>
using std::calloc;
using std::free;
using std::exit;
using std::rand;
using std::abort;
using std::malloc;
using std::free;
#else /* HAVE_STDLIB_H */
#include <stdlib.h>
#endif
#else /* JANUS_STLPORT */
#include <stdlib.h>
#endif /* JANUS_STLPORT */

#ifdef HAVE_CSTDIO
#include <cstdio>
#ifndef AVOID_BUG1472
using std::fopen;
using std::fclose;
using std::FILE;
using std::fprintf;
using std::fscanf;
using std::printf;
using std::perror;
using std::feof;
#endif
#else
#include <stdio.h>
#endif

#ifdef HAVE_IOSTREAM
#include <iostream>
#else
#include <iostream.h>
#endif

#ifdef HAVE_STRING
#include <string>
#else
#include <string.h>
#endif

#ifndef JANUS_STLPORT
#ifdef HAVE_CMATH
#include <cmath>
using std::fabs;
using std::sqrt;
#else
#include <math.h>
#endif
#else /* JANUS_STLPORT */
#include <math.h>
#endif /* JANUS_STLPORT */

#ifdef HAVE_MATH_H
#include <math.h>
#endif

#ifdef HAVE_MAP
#include <map>
#else
#include <map.h>
#endif

#else /*ndef TRILINOS_NO_CONFIG_H*/

#include <iostream.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <math.h>

#endif /*ndef TRILINOS_NO_CONFIG_H*/

void Trilinos_Util_read_hb(const char *data_file, int MyPID,
    int *N_global, int *n_nonzeros,
    double **val, int **bindx,
    double **x, double **b, double **xexact);

void Trilinos_Util_read_hb(const char *data_file, int MyPID,
    int *N_global, int *n_nonzeros,
    double **val, int **bindx);

void Trilinos_Util_read_coo(const char *data_file, int MyPID,
    int *N_global, int *n_nonzeros,
    double **val, int **bindx,
    double **x, double **b, double **xexact);

double Trilinos_Util_smsrres (int m, int n,
    double *val, int *indx,
    double *xlocal, double *x, double *b);

double Trilinos_Util_scscres (int isym, int m, int n,
    double *val, int *indx, int *pntr,
    double *x, double *b);

void  Trilinos_Util_scscmv (int isym, int m, int n,
    double *val, int *indx, int *pntr,
    double *x, double *b);

double Trilinos_Util_svbrres (int m, int n, int m_blk,
    double *val, int *indx, int *bindx, int *rpntr,
    int *cpntr, int *bpntrb, int *bpntre,
    double *x, double *b);

void Trilinos_Util_msr2vbr(
    double val[], int indx[], int rnptr[],
    int cnptr[], int bnptr[],
    int bindx[], int msr_bindx[], double msr_val[],
    int total_blk_rows, int total_blk_cols, int blk_space,
    int nz_space, int blk_type);

int Trilinos_Util_find_block_col(int cnptr[], int column, int
    max_blocks, int blk_size);

int Trilinos_Util_find_block_in_row(int bindx[], int bnptr[], int
    blk_row, int blk_col,
    int indx[], int no_elements, double val[],
    int blk_space, int nz_space);


void Trilinos_Util_add_new_ele(int cnptr[], int col, int blk_row, int
    bindx[], int bnptr[],
    int indx[], double val[], int row, double new_ele,
    int maxcols, int blk_space, int nz_space, int blk_type);

int Trilinos_Util_find_closest_not_larger(int key, int list[], int length);

void Trilinos_Util_convert_values_to_ptrs(int array[], int length, int start);

int Trilinos_Util_csrcsc(int n, int n2, int job, int ipos, double * a,
    int *ja, int *ia, double *ao, int *jao, int *iao);

int Trilinos_Util_csrmsr( int n, double *a, int *ja, int *ia, double *ao,
    int *jao, double *wk, int *iwk);

int Trilinos_Util_ssrcsr( int job, int value2, int nrow, double *a,
    int *ja, int *ia, int nzmax,
    double *ao, int *jao, int *iao, int *indu,
    int *iwk);

int Trilinos_Util_coocsr( int nrow, int nnz, double *a, int *ir, int *jc,
    double *ao, int *jao, int *iao);


struct SPBLASMAT_STRUCT {
  int n;
  double *val;
  int *indx;
  int *bindx;
  int *rpntr;
  int *cpntr;
  int *bpntrb;
  int *bpntre;
  int buffersize;
  int bufferstride;
  double *buffer;
  int *ncolvec;
  double nops_per_rhs;
  int minblocksize;
  int maxblocksize;
};
typedef struct SPBLASMAT_STRUCT SPBLASMAT;

#define MAXNRHS 1

void  Trilinos_Util_duscr_vbr(
    int n, double *val, int *indx, int *bindx,
    int *rpntr, int *cpntr, int *bpntrb, int *bpntre,
    SPBLASMAT *A
    );

void Trilinos_Util_dusmm(
    int m, int nrhs, int k, double alpha, SPBLASMAT *A,
    double *x, int xstride, double beta, double *b, int bstride
    );
void  Trilinos_Util_dusds_vbr( SPBLASMAT *A);

void Trilinos_Util_write_vec(const char *filename, int n_equations, double *x);

void Trilinos_Util_read_vec(const char *filename, int n_equations, double *x);

// Start of Epetra dependencies
#ifdef HAVE_TRIUTILS_EPETRA

class Epetra_Comm;
class Epetra_Vector;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_CrsMatrix;
class Epetra_VbrMatrix;
class Epetra_MultiVector;

#include "Epetra_ConfigDefs.h"
#include "Trilinos_Util_ReadMatrixMarket2Epetra.h"

void Trilinos_Util_ReadHb2Epetra_internal(const char *data_file,
    const Epetra_Comm  &comm,
    Epetra_Map *& map,
    Epetra_CrsMatrix *& A,
    Epetra_Vector *& x,
    Epetra_Vector *& b,
    Epetra_Vector *&xexact,
    bool FakeLongLong);

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES

void Trilinos_Util_ReadHb2Epetra(const char *data_file,
    const Epetra_Comm  &comm,
    Epetra_Map *& map,
    Epetra_CrsMatrix *& A,
    Epetra_Vector *& x,
    Epetra_Vector *& b,
    Epetra_Vector *&xexact);

#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

void Trilinos_Util_ReadHb2Epetra64(const char *data_file,
    const Epetra_Comm  &comm,
    Epetra_Map *& map,
    Epetra_CrsMatrix *& A,
    Epetra_Vector *& x,
    Epetra_Vector *& b,
    Epetra_Vector *&xexact);

#endif

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES

void Trilinos_Util_ReadHpc2Epetra(const char *data_file,
    const Epetra_Comm  &comm,
    Epetra_Map *& map,
    Epetra_CrsMatrix *& A,
    Epetra_Vector *& x,
    Epetra_Vector *& b,
    Epetra_Vector *&xexact);

#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

void Trilinos_Util_ReadHpc2Epetra64(const char *data_file,
    const Epetra_Comm  &comm,
    Epetra_Map *& map,
    Epetra_CrsMatrix *& A,
    Epetra_Vector *& x,
    Epetra_Vector *& b,
    Epetra_Vector *&xexact);

#endif

// CJ TODO FIXME: Trilinos_Util_ReadHb2EpetraVbr available only if 32 bit GIDs available.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES

void Trilinos_Util_ReadHb2EpetraVbr(const char *data_file, const char * partitioning,
    const Epetra_Comm  &comm,
    Epetra_BlockMap *& map,
    Epetra_VbrMatrix *& A,
    Epetra_Vector *& x,
    Epetra_Vector *& b,
    Epetra_Vector *&xexact);

#endif

// CJ TODO FIXME: Trilinos_Util_distrib_msr_matrix available only if 32 bit GIDs available.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES

void Trilinos_Util_distrib_msr_matrix(
    const Epetra_Comm & Comm,
    int *N_global, int *n_nonzeros,
    int *N_update, int **update,
    double **val, int **bindx,
    double **x, double **b, double **xexact
    );

void Trilinos_Util_distrib_msr_matrix(
    const Epetra_Comm & Comm,
    int *N_global, int *n_nonzeros,
    int *N_update, int **update,
    double **val, int **bindx
    );

#endif

// CJ TODO FIXME: Trilinos_Util_distrib_vbr_matrix available only if 32 bit GIDs available.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES

void Trilinos_Util_distrib_vbr_matrix(const Epetra_Comm & Comm,
    int *N_global, int *N_blk_global,
    int *n_nonzeros,  int *n_blk_nonzeros,
    int *N_update, int **update,
    double **val, int **indx, int **rpntr, int **cpntr,
    int **bpntr, int **bindx,
    double **x, double **b, double **xexact);

#endif

void Trilinos_Util_create_vbr(const Epetra_Comm & Comm, const char *part_file,
    int *N_global, int *N_blk_global,
    int *n_nonzeros, int *n_blk_nonzeros,
    int *N_update, int **update,
    int *bindx_msr, double *val_msr,
    double **val, int **indx, int **rpntr, int **cpntr,
    int **bpntr, int **bindx);

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES

void Trilinos_Util_GenerateCrsProblem(
    int nx, int ny, int npoints, int * xoff, int * yoff,
    const Epetra_Comm  &comm,
    Epetra_Map *& map,
    Epetra_CrsMatrix *& A,
    Epetra_Vector *& x,
    Epetra_Vector *& b,
    Epetra_Vector *&xexact, int indexBase = 0);

void Trilinos_Util_GenerateCrsProblem(
    int nx, int ny, int npoints, int * xoff, int * yoff, int nrhs,
    const Epetra_Comm  &comm,
    Epetra_Map *& map,
    Epetra_CrsMatrix *& A,
    Epetra_MultiVector *& x,
    Epetra_MultiVector *& b,
    Epetra_MultiVector *&xexact, int indexBase = 0);

#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

void Trilinos_Util_GenerateCrsProblem64(
    int nx, int ny, int npoints, int * xoff, int * yoff,
    const Epetra_Comm  &comm,
    Epetra_Map *& map,
    Epetra_CrsMatrix *& A,
    Epetra_Vector *& x,
    Epetra_Vector *& b,
    Epetra_Vector *&xexact, long long indexBase = 0);

void Trilinos_Util_GenerateCrsProblem64(
    int nx, int ny, int npoints, int * xoff, int * yoff, int nrhs,
    const Epetra_Comm  &comm,
    Epetra_Map *& map,
    Epetra_CrsMatrix *& A,
    Epetra_MultiVector *& x,
    Epetra_MultiVector *& b,
    Epetra_MultiVector *&xexact, long long indexBase = 0);

#endif

// CJ TODO FIXME: Trilinos_Util_GenerateVbrProblem available only if 32 bit GIDs available.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES

void Trilinos_Util_GenerateVbrProblem(
    int nx, int ny, int npoints, int * xoff, int * yoff,
    int nsizes, int * sizes,
    const Epetra_Comm  &comm,
    Epetra_BlockMap *& map,
    Epetra_VbrMatrix *& A,
    Epetra_Vector *& x,
    Epetra_Vector *& b,
    Epetra_Vector *&xexact);

void Trilinos_Util_GenerateVbrProblem(
    int nx, int ny, int npoints, int * xoff, int * yoff,
    int nsizes, int * sizes, int nrhs,
    const Epetra_Comm  &comm,
    Epetra_BlockMap *& map,
    Epetra_VbrMatrix *& A,
    Epetra_MultiVector *& x,
    Epetra_MultiVector *& b,
    Epetra_MultiVector *&xexact);

#endif

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES

int Trilinos_Util_ReadTriples2Epetra( const char *data_file,
    bool symmetric,
    const Epetra_Comm  &comm,
    Epetra_Map *& map,
    Epetra_CrsMatrix *& A,
    Epetra_Vector *& x,
    Epetra_Vector *& b,
    Epetra_Vector *&xexact,
    bool NonUniformMap=false,
    bool TimDavisHeader=false,
    bool ZeroBased=false ) ;

#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

int Trilinos_Util_ReadTriples2Epetra64(
    const char *data_file,
    bool symmetric,
    const Epetra_Comm  &comm,
    Epetra_Map *& map,
    Epetra_CrsMatrix *& A,
    Epetra_Vector *& x,
    Epetra_Vector *& b,
    Epetra_Vector *&xexact,
    bool NonUniformMap=false,
    bool TimDavisHeader=false,
    bool ZeroBased=false );

#endif

#endif //ifdef HAVE_TRIUTILS_EPETRA

#endif /* _TRILINOS_UTIL_H_ */
