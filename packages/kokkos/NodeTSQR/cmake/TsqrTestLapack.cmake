# @HEADER
# ***********************************************************************
#  
#           Kokkos: Node API and Parallel Node Kernels
#               Copyright (2009) Sandia Corporation
#  
#  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
#  license for use of this work by or on behalf of the U.S. Government.
#  
#  This library is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation; either version 2.1 of the
#  License, or (at your option) any later version.
#   
#  This library is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#   
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
#  USA
#  Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
#  
# ***********************************************************************
# @HEADER

IF (NOT ${PACKAGE_NAME}_FINISHED_FIRST_CONFIGURE)
  message (STATUS "Checking LAPACK QR factorization routine names")
ENDIF()

# Use TsqrCheckLapackRoutine to search the LAPACK library for
# {D,S,Z,C}LARFP.  Also, search for DLARFG, just as a sanity check.
# The LARFP routines were added to LAPACK in version 3.2; they produce
# a Householder reflector with the BETA output always nonnegative.  In
# LAPACK version 3.2.2, the LARFP routines were renamed to LARFGP.
# The (original) LARFG routines may produce a negative BETA output (in
# the real case).

# The check for DLARFG is just a sanity check to make sure that CMake is 
# finding the LAPACK library.  All modern versions of LAPACK should have 
# DLARFG.
TSQR_CHECK_LAPACK_ROUTINE ("DLARFG" HAVE_LAPACK_DLARFG)
TSQR_CHECK_LAPACK_ROUTINE ("SLARFG" HAVE_LAPACK_SLARFG)
TSQR_CHECK_LAPACK_ROUTINE ("ZLARFG" HAVE_LAPACK_ZLARFG)
TSQR_CHECK_LAPACK_ROUTINE ("CLARFG" HAVE_LAPACK_CLARFG)

# LAPACK versions 3.2 and 3.2.1 have _LARFP, which is a variant of
# DLARFG that always produces a nonnegative BETA output.  LAPACK
# version 3.2.2 (released ~ July 2010) renamed these routines to
# _LARFGP.
# 
# This rename means that we have to check for three possibilities.  If
# the _LARFGP routines exist, we use them preferentially.  Otherwise
# we use the _LARFP routines if they exist.  If neither exist, we've
# chosen for benchmarking purposes (even though it may result in
# negative BETA and therefore an R factor in the QR factorization with
# negative diagonal entries) to substitute _LARFG.  We will print out
# a big warning message if that happens, because it may break the
# assumptions of some of the orthogonalization routines.
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
TSQR_CHECK_LAPACK_ROUTINE ("DLARFP" HAVE_LAPACK_DLARFP)
TSQR_CHECK_LAPACK_ROUTINE ("SLARFP" HAVE_LAPACK_SLARFP)
TSQR_CHECK_LAPACK_ROUTINE ("ZLARFP" HAVE_LAPACK_ZLARFP)
TSQR_CHECK_LAPACK_ROUTINE ("CLARFP" HAVE_LAPACK_CLARFP)
TSQR_CHECK_LAPACK_ROUTINE ("DLARFGP" HAVE_LAPACK_DLARFGP)
TSQR_CHECK_LAPACK_ROUTINE ("SLARFGP" HAVE_LAPACK_SLARFGP)
TSQR_CHECK_LAPACK_ROUTINE ("ZLARFGP" HAVE_LAPACK_ZLARFGP)
TSQR_CHECK_LAPACK_ROUTINE ("CLARFGP" HAVE_LAPACK_CLARFGP)

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
TSQR_CHECK_LAPACK_ROUTINE ("DGEQRFP" HAVE_LAPACK_DGEQRFP)
TSQR_CHECK_LAPACK_ROUTINE ("SGEQRFP" HAVE_LAPACK_SGEQRFP)
TSQR_CHECK_LAPACK_ROUTINE ("ZGEQRFP" HAVE_LAPACK_ZGEQRFP)
TSQR_CHECK_LAPACK_ROUTINE ("CGEQRFP" HAVE_LAPACK_CGEQRFP)

# Similar situation as mentioned above with _GEQRF.
TSQR_CHECK_LAPACK_ROUTINE ("DGEQR2P" HAVE_LAPACK_DGEQR2P)
TSQR_CHECK_LAPACK_ROUTINE ("SGEQR2P" HAVE_LAPACK_SGEQR2P)
TSQR_CHECK_LAPACK_ROUTINE ("ZGEQR2P" HAVE_LAPACK_ZGEQR2P)
TSQR_CHECK_LAPACK_ROUTINE ("CGEQR2P" HAVE_LAPACK_CGEQR2P)
