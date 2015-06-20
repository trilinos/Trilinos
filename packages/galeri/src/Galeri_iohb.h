// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
// @HEADER

#include<stdio.h>
#include<stdlib.h>

#ifdef HAVE_MALLOC_H
#include<malloc.h>
#endif

int readHB_info(const char* filename, int* M, int* N, int* nz, char** Type,
                                                      int* Nrhs);

int readHB_header(FILE* in_file, char* Title, char* Key, char* Type,
                    int* Nrow, int* Ncol, int* Nnzero, int* Nrhs,
                    char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt,
                    int* Ptrcrd, int* Indcrd, int* Valcrd, int* Rhscrd,
                    char *Rhstype);

int readHB_mat_double(const char* filename, int colptr[], int rowind[],
                                                                 double val[]);

int readHB_newmat_double(const char* filename, int* M, int* N, int* nonzeros,
                         int** colptr, int** rowind, double** val);

int readHB_aux_double(const char* filename, const char AuxType, double b[]);

int readHB_newaux_double(const char* filename, const char AuxType, double** b);

int writeHB_mat_double(const char* filename, int M, int N,
                        int nz, const int colptr[], const int rowind[],
                        const double val[], int Nrhs, const double rhs[],
                        const double guess[], const double exact[],
                        const char* Title, const char* Key, const char* Type,
                        char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt,
                        const char* Rhstype);

int readHB_mat_char(const char* filename, int colptr[], int rowind[],
                                           char val[], char* Valfmt);

int readHB_newmat_char(const char* filename, int* M, int* N, int* nonzeros, int** colptr,
                          int** rowind, char** val, char** Valfmt);

int readHB_aux_char(const char* filename, const char AuxType, char b[]);

int readHB_newaux_char(const char* filename, const char AuxType, char** b, char** Rhsfmt);

int writeHB_mat_char(const char* filename, int M, int N,
                        int nz, const int colptr[], const int rowind[],
                        const char val[], int Nrhs, const char rhs[],
                        const char guess[], const char exact[],
                        const char* Title, const char* Key, const char* Type,
                        char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt,
                        const char* Rhstype);

int ParseIfmt(char* fmt, int* perline, int* width);

int ParseRfmt(char* fmt, int* perline, int* width, int* prec, int* flag);

void IOHBTerminate(const char* message);
