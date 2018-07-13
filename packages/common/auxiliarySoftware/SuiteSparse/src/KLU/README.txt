KLU Version 1.0, May 31, 2007, by Timothy A. Davis and Ekanathan Palamadai.
Copyright (C) 2004-2007, University of Florida
KLU is also available under other licenses; contact authors for details.
http://www.cise.ufl.edu/research/sparse

Requires the AMD, COLAMD, and BTF libraries, in ../AMD, ../COLAMD, and ../BTF,
respectively.  Requires the ../UFconfig/UFconfig.mk configuration file.
Optionally uses CHOLMOD (KLU/User example ordering).  The Tcov tests and
the Demo both require CHOLMOD.

To compile the libklu.a library, type "make".  The compiled library is located
in KLU/Lib/libklu.a.  Compile code that uses KLU with -IKLU/Include.

Type "make clean" to remove all but the compiled library, and "make distclean"
to remove all files not in the original distribution.

--------------------------------------------------------------------------------

KLU is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This Module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this Module; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

--------------------------------------------------------------------------------

A full text of the license is in Doc/lesser.txt.

Files in this distribution:

    Include		include files
    README.txt		this file
    Source		source code

./Include:
    klu.h		user include file
    klu_internal.h	internal include file, not needed by the user
    klu_version.h	internal include file, not needed by the user

./Source:
    klu_analyze.c	klu_analyze and supporting functions
    klu_analyze_given.c	klu_analyze_given and supporting functions
    klu.c		kernel factor/solve functions, not user-callable
    klu_defaults.c	klu_defaults function
    klu_diagnostics.c	klu_rcond, klu_condest, klu_rgrowth, kluflops
    klu_dump.c		debugging functions
    klu_extract.c	klu_extract
    klu_factor.c	klu_factor and supporting functions
    klu_free_numeric.c	klu_free_numeric function
    klu_free_symbolic.c	klu_free_symbolic function
    klu_kernel.c	kernel factor functions, not user-callable
    klu_memory.c	klu_malloc, klu_free, klu_realloc, and supporing func.
    klu_refactor.c	klu_refactor function
    klu_scale.c		klu_scale function
    klu_solve.c		klu_solve function
    klu_sort.c		klu_sort and supporting functions
    klu_tsolve.c	klu_tsovle function
