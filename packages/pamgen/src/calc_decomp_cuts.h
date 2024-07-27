// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef calc_decomp_cutsH
#define calc_decomp_cutsH
namespace PAMGEN_NEVADA {

long long dom_decomp_2d(const long long Nx, const long long Ny,
		  const long long Np, long long *pNGx, long long *pNGy);

long long dom_decomp_3d(const long long Nx, const long long Ny, const long long Nz,
		  const long long Np, long long *pNGx, long long *pNGy, long long *pNGz);
}
#endif
