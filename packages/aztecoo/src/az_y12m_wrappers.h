/*@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef _AZ_Y12M_WRAPPERS_H_
#define _AZ_Y12M_WRAPPERS_H_

#include "az_f77func.h"

#   define Y12MBF_F77                 F77_FUNC(y12mbf,Y12MBF)
#   define Y12MCF_F77                 F77_FUNC(y12mcf,Y12MCF)
#   define Y12MDF_F77                 F77_FUNC(y12mdf,Y12MDF)

#ifdef __cplusplus
extern "C" {
#include <stdio.h>
#endif

void PREFIX Y12MBF_F77(int *n, int *z, double val[], int snr[], int *nn,
                      int rnr[], int *nn1, int ha[], int *iha, double aflag[],
                      int iflag[], int *ifail);

void PREFIX Y12MCF_F77(int *n, int *z, double val[], int snr[], int *nn,
                      int rnr[], int *nn1, double pivot[], double b[], int ha[],
                      int *iha, double aflag[], int iflag[], int *ifail);

void PREFIX Y12MDF_F77(int *n, double val[], int *nn, double b[], double pivot[],
                      int snr[], int ha[], int *iha, int iflag[], int *ifail);

#ifdef __cplusplus
}
#endif

#endif /* _AZ_Y12M_WRAPPERS_H_ */
