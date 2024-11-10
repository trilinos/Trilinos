// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ConfigDefs.hpp"

#define GAUSSQ_F77 F77_FUNC(gaussq,GAUSSQ)

#ifdef __cplusplus
extern "C" {
#endif

extern void GAUSSQ_F77(int *kind, int *n, double *alpha, double *beta,
		       int *kpts, double *endpts, double *b, double *x, double *w);

#ifdef __cplusplus
}
#endif
