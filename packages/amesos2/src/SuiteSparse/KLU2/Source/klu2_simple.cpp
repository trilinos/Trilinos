// @HEADER
// ***********************************************************************
//
//                   KLU2: A Direct Linear Solver package
//                    Copyright 2011 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, with Sandia Corporation, the 
// U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Mike A. Heroux (maherou@sandia.gov)
//
// KLU2 is derived work from KLU, licensed under LGPL, and copyrighted by
// University of Florida. The Authors of KLU are Timothy A. Davis and
// Eka Palamadai. See Doc/KLU_README.txt for the licensing and copyright
// information for KLU.
//
// ***********************************************************************
// @HEADER

/* klu_simple: a simple KLU demo; solution is x = (1,2,3,4,5) */
/* TODO : Need to go to an example directory */

#include <stdio.h>
#include <iostream>
#include <complex>
#include "klu2_defaults.hpp"
#include "klu2_analyze.hpp"
#include "klu2_factor.hpp"
#include "klu2_solve.hpp"
#include "klu2_free_symbolic.hpp"
#include "klu2_free_numeric.hpp"

int    n = 5 ;
int    Ap [ ] = {0, 2, 5, 9, 10, 12} ;
int    Ai [ ] = { 0,  1,  0,   2,  4,  1,  2,  3,   4,  2,  1,  4} ;
double Ax [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;
double b [ ] = {8., 45., -3., 3., 19.} ;

std::complex<double> Ax_cd [] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;
std::complex<double> b_cd [ ] = {8., 45., -3., 3., 19.} ;

int main (void)
{
    klu_symbolic<double, int> *Symbolic ;
    klu_numeric<double, int> *Numeric ;
    klu_common<double, int> Common ;
    int i ;
    klu_defaults<double, int> (&Common) ;
    Symbolic = klu_analyze<double, int> (n, Ap, Ai, &Common) ;
    Numeric = klu_factor<double, int> (Ap, Ai, Ax, Symbolic, &Common) ;
    klu_solve<double, int> (Symbolic, Numeric, 5, 1, b, &Common) ;
    klu_free_symbolic<double, int> (&Symbolic, &Common) ;
    klu_free_numeric<double, int> (&Numeric, &Common) ;
    for (i = 0 ; i < n ; i++) printf ("x [%d] = %g\n", i, b [i]) ;

    std::cout << "Complex double case " << std::endl ;
    typedef std::complex<double> ComplexD ;
    klu_symbolic<ComplexD, int> *Symbolic_C ;
    klu_numeric<ComplexD, int> *Numeric_C ;
    klu_common<ComplexD, int> Common_C ;
    klu_defaults<ComplexD, int> (&Common_C) ;
    Symbolic_C = klu_analyze<ComplexD, int> (n, Ap, Ai, &Common_C) ;
    Numeric_C = klu_factor<ComplexD, int> (Ap, Ai, Ax_cd, Symbolic_C, &Common_C) ;
    klu_solve<ComplexD, int> (Symbolic_C, Numeric_C, 5, 1, b_cd, &Common_C) ;
    klu_free_symbolic<ComplexD, int> (&Symbolic_C, &Common_C) ;
    klu_free_numeric<ComplexD, int> (&Numeric_C, &Common_C) ;
    for (i = 0 ; i < n ; i++) 
        std::cout << "x [" << i << "] = "<< b_cd[i] << std::endl ;

    return (0) ;
}

