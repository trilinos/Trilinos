AMD Version 2.2, Copyright (c) 2007 by Timothy A.
Davis, Patrick R. Amestoy, and Iain S. Duff.  All Rights Reserved.
AMD is available under alternate licences; contact T. Davis for details.

AMD:  a set of routines for permuting sparse matrices prior to
    factorization.  Includes a version in C, a version in Fortran, and a MATLAB
    mexFunction.

Requires UFconfig, in the ../UFconfig directory relative to this directory.

Quick start (Unix, or Windows with Cygwin):

    To compile, test, and install AMD, you may wish to first configure the
    installation by editting the ../UFconfig/UFconfig.mk file.  Next, cd to this
    directory (AMD) and type "make" (or "make lib" if you do not have MATLAB).
    To compile and run a demo program for the Fortran version, type
    "make fortran".  When done, type "make clean" to remove unused *.o files
    (keeps the compiled libraries and demo programs).  See the User Guide
    (Doc/AMD_UserGuide.pdf), or ../UFconfig/UFconfig.mk for more details.

Quick start (for MATLAB users);

    To compile, test, and install the AMD mexFunction, cd to the
    AMD/MATLAB directory and type amd_make at the MATLAB prompt.

If you have MATLAB 7.2 or earlier and use "make mex", you must first edit
UFconfig/UFconfig.h to remove the "-largeArrayDims" option from the MEX command
(or just use amd_make.m inside MATLAB).

-------------------------------------------------------------------------------

AMD License:

    Your use or distribution of AMD or any modified version of
    AMD implies that you agree to this License.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301
    USA

    Permission is hereby granted to use or copy this program under the
    terms of the GNU LGPL, provided that the Copyright, this License,
    and the Availability of the original version is retained on all copies.
    User documentation of any code that uses this code or any modified
    version of this code must cite the Copyright, this License, the
    Availability note, and "Used by permission." Permission to modify
    the code and to distribute modified code is granted, provided the
    Copyright, this License, and the Availability note are retained,
    and a notice that the code was modified is included.

Availability:

    http://www.cise.ufl.edu/research/sparse/amd

-------------------------------------------------------------------------------

This is the AMD README file.  It is a terse overview of AMD.
Refer to the User Guide (Doc/AMD_UserGuide.pdf) for how to install
and use AMD.

Description:

    AMD is a set of routines for pre-ordering sparse matrices prior to Cholesky
    or LU factorization, using the approximate minimum degree ordering
    algorithm.  Written in ANSI/ISO C with a MATLAB interface, and in
    Fortran 77.

Authors:

    Timothy A. Davis (davis at cise.ufl.edu), University of Florida.
    Patrick R. Amestory, ENSEEIHT, Toulouse, France.
    Iain S. Duff, Rutherford Appleton Laboratory, UK.

Acknowledgements:

    This work was supported by the National Science Foundation, under
    grants DMS-9504974, DMS-9803599, and CCR-0203270.

    Portions of this work were done while on sabbatical at Stanford University
    and Lawrence Berkeley National Laboratory (with funding from the SciDAC
    program).  I would like to thank Gene Golub, Esmond Ng, and Horst Simon
    for making this sabbatical possible.

-------------------------------------------------------------------------------
Files and directories in the AMD distribution:
-------------------------------------------------------------------------------

    ---------------------------------------------------------------------------
    Subdirectories of the AMD directory:
    ---------------------------------------------------------------------------

    Source	primary source code
    Include	include file for use in your code that calls AMD

    ---------------------------------------------------------------------------
    Files in the AMD directory:
    ---------------------------------------------------------------------------

    README.txt	this file

    ---------------------------------------------------------------------------
    Source directory:
    ---------------------------------------------------------------------------

    amd_order.c			user-callable, primary AMD ordering routine
    amd_control.c		user-callable, prints the control parameters
    amd_defaults.c		user-callable, sets default control parameters
    amd_info.c			user-callable, prints the statistics from AMD

    amd_1.c			non-user-callable, construct A+A'
    amd_2.c			user-callable, primary ordering kernel
				(a C version of amd.f and amdbar.f, with
				post-ordering added)
    amd_aat.c			non-user-callable, computes nnz (A+A')
    amd_dump.c			non-user-callable, debugging routines
    amd_postorder.c		non-user-callable, postorder
    amd_post_tree.c		non-user-callable, postorder just one tree
    amd_valid.c			non-user-callable, verifies a matrix
    amd_preprocess.c		non-user-callable, computes A', removes duplic

    amd.f			user-callable Fortran 77 version
    amdbar.f			user-callable Fortran 77 version

    ---------------------------------------------------------------------------
    Include directory:
    ---------------------------------------------------------------------------

    amd.h			include file for C programs that use AMD
    amd_internal.h		non-user-callable, include file for AMD

