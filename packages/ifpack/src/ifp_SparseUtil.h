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

#ifndef _IFP_SPARSEUTIL_H_
#define _IFP_SPARSEUTIL_H_

void shell_sort(
  const int n,
  int x[]);

void allocate_ilu(
  const int levfill,
  const int n,
  int *nzl, int *nzu,
  const int ia[], const int ja[],
  int *ial[], int *jal[],
  int *iau[], int *jau[],
  int growthl, int growthu);

int symbolic_ilu(
  const int levinc,
  const int n,
  int *nzl,
  int *nzu,
  const int ia[], const int ja[],
  int ial[], int jal[],
  int iau[], int jau[]);

#endif // _IFP_SPARSEUTIL_H_
