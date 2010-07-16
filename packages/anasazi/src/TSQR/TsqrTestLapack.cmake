# @HEADER
# ***********************************************************************
#
#                 Anasazi: Block Eigensolvers Package
#                 Copyright (2010) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
#
# ***********************************************************************
# @HEADER

message (STATUS "Checking LAPACK QR factorization routine names")

# Use the CheckFortranFunctionExists module to search the LAPACK library 
# for {D,S,Z,C}LARFP.  Also, search for DLARFG, just as a sanity check.
# The LARFP functions were added to LAPACK in version 3.2; they produce 
# a Householder reflector with the BETA output always nonnegative.  The 
# LARFG functions may produce a negative BETA output (in the real case).
#
# Note that CheckFortranFunctionExists requires support for compiling
# Fortran (77) files.  That's OK because the ParallelTSQR project
# requires a Fortran compiler to build anyway.
include (CheckFortranFunctionExists)

# Remember the current CMAKE_REQUIRED_LIBRARIES so we can restore it
# after we are done testing.
set (CMAKE_REQUIRED_LIBRARIES_SAVE ${CMAKE_REQUIRED_LIBRARIES})
set (CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${TPL_LAPACK_LIBRARIES} ${TPL_BLAS_LIBRARIES})

#message(STATUS "\nCMAKE_REQUIRED_LIBRARIES == ${CMAKE_REQUIRED_LIBRARIES}\n")

# The check for DLARFG is just a sanity check to make sure that CMake is 
# finding the LAPACK library.  All modern versions of LAPACK should have 
# DLARFG.
check_fortran_function_exists ("DLARFG" HAVE_LAPACK_DLARFG)
check_fortran_function_exists ("SLARFG" HAVE_LAPACK_SLARFG)
check_fortran_function_exists ("ZLARFG" HAVE_LAPACK_ZLARFG)
check_fortran_function_exists ("CLARFG" HAVE_LAPACK_CLARFG)

# LAPACK versions 3.2 and 3.2.1 have _LARFP, which is a variant of
# DLARFG that always produces a nonnegative BETA output.  In mid-June
# 2010, these routines will be renamed to _LARFGP.  As of 21 June
# 2010, this has already happened in LAPACK's SVN repository.  The
# changes will go out in the next bugfix release of LAPACK.
# 
# This rename means that we have to check for three possibilities.  If
# the _LARFGP routines exist, we use them preferentially.  Otherwise
# we use the _LARFP routines if they exist.  If neither exist, we've
# chosen for benchmarking purposes (even though it may result in
# negative BETA and therefore an R factor in the QR factorization with
# negative diagonal entries) to substitute _LARFG.
#
# FIXME (mfh 22 June 2010) If none of these routines exist, we should
# provide a fallback implementation that does extra work to ensure
# that the diagonal entries of the R factor in the QR factorization
# are nonnegative.  This has to be taken care of in TSQR itself,
# though: we can't just flip the signs of the diagonal entries and the
# signs of their corresponding TAU scaling factors from GEQR2.  The
# problem could be solved with an extra write of the output matrix in
# explicit_Q() or apply(), by keeping the appropriate +/-1 scaling
# factors.  I haven't done that yet.  
#
# Probably the right thing to do in the TSQR interface would be to
# make producing an R factor with nonnegative diagonal entries an
# option, which is enabled by default if _LARFGP or _LARFP exist, and
# disabled by default if neither of these routines exist.
check_fortran_function_exists ("DLARFP" HAVE_LAPACK_DLARFP)
check_fortran_function_exists ("SLARFP" HAVE_LAPACK_SLARFP)
check_fortran_function_exists ("ZLARFP" HAVE_LAPACK_ZLARFP)
check_fortran_function_exists ("CLARFP" HAVE_LAPACK_CLARFP)
check_fortran_function_exists ("DLARFGP" HAVE_LAPACK_DLARFGP)
check_fortran_function_exists ("SLARFGP" HAVE_LAPACK_SLARFGP)
check_fortran_function_exists ("ZLARFGP" HAVE_LAPACK_ZLARFGP)
check_fortran_function_exists ("CLARFGP" HAVE_LAPACK_CLARFGP)

# The same LAPACK change mentioned above from _LARFP to _LARFGP
# changed the name of the QR factorization routines that promise an R
# factor with nonnegative diagonal entries, from _GEQRF (which was the
# default QR factorization used when solving least-squares problems)
# to _GEQRFP (which is not the default).
#
# For the implemenetation of our wrapper function LAPACK::GEQRF, we
# use _GEQRFP preferentially if it exists.  Otherwise, we fall back to
# _GEQRF.  The latter only promises an R factor with nonnegative
# diagonal entries in LAPACK versions 3.2 and 3.2.1.  In LAPACK
# releases before and after those, it may produce an R factor with
# negative diagonal entries.
check_fortran_function_exists ("DGEQRFP" HAVE_LAPACK_DGEQRFP)
check_fortran_function_exists ("SGEQRFP" HAVE_LAPACK_SGEQRFP)
check_fortran_function_exists ("ZGEQRFP" HAVE_LAPACK_ZGEQRFP)
check_fortran_function_exists ("CGEQRFP" HAVE_LAPACK_CGEQRFP)

# Similar situation as mentioned above with _GEQRF.
check_fortran_function_exists ("DGEQR2P" HAVE_LAPACK_DGEQR2P)
check_fortran_function_exists ("SGEQR2P" HAVE_LAPACK_SGEQR2P)
check_fortran_function_exists ("ZGEQR2P" HAVE_LAPACK_ZGEQR2P)
check_fortran_function_exists ("CGEQR2P" HAVE_LAPACK_CGEQR2P)

# Write out a header file with the appropriate #defines.
configure_file ("TSQR/Tsqr_Config.hpp.in" "TSQR/Tsqr_Config.hpp")

# Restore the original value of CMAKE_REQUIRED_LIBRARIES.
set (CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES_SAVE})

