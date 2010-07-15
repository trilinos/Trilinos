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

message (STATUS "Generating LAPACK wrappers for TSQR")

# Generate a Fortran wrapper module for the LAPACK routines that
# compute Householder vectors with nonnegative real BETA output.  We
# wrap the routines, because the LAPACK routines' name depends on the
# LAPACK release version.  We wrap in Fortran because it's harder to
# call C from Fortran than vice versa.  We generate the Fortran code
# directly because it's easier to use Fortran modules than Fortran
# include files.
#
# Note that if your LAPACK version is < 3.2, then LAPACK has no
# Householder reflector routines that promise nonnegative real BETA
# output.  DLARFG and SLARFG may produce negative BETA output.  In
# LAPACK versions 3.2 and 3.2.1, the desired routines are called
# DLARFP resp. SLARFP.  In LAPACK version 3.2.2 (latest as of writing
# this, on 30 Jun 2010), the desired routines are called DLARFGP
# resp. SLARFGP.
set (_FORTRAN_MODULE_FILENAME "Tsqr_HouseholderReflector.f90")
file (WRITE "${_FORTRAN_MODULE_FILENAME}" "module TsqrHouseholderReflector
  implicit none

  contains
")
foreach (_LAPACK_PREFIX "D" "S" "Z" "C")
  set (_HAVE_LARFGP "HAVE_LAPACK_${_LAPACK_PREFIX}LARFGP")
  set (_HAVE_LARFP "HAVE_LAPACK_${_LAPACK_PREFIX}LARFP")
  set (_HAVE_LARFG "HAVE_LAPACK_${_LAPACK_PREFIX}LARFG")

  if ("${_HAVE_LARFGP}")
    set (_LAPACK_ROUTINE "${_LAPACK_PREFIX}LARFGP")
  elseif ("${_HAVE_LARFP}")
    set (_LAPACK_ROUTINE "${_LAPACK_PREFIX}LARFP")
  elseif ("${_HAVE_LARFG}")
    set (_LAPACK_ROUTINE "${_LAPACK_PREFIX}LARFG")
  else ()
    message (FATAL_ERROR "Failed to link with LAPACK routine ${_LAPACK_PREFIX}LARFG, which should exist")
  endif ()

  if ("${_LAPACK_PREFIX}" STREQUAL "D")
    set (_DATATYPE "real (8)")
  elseif ("${_LAPACK_PREFIX}" STREQUAL "S")
    set (_DATATYPE "real (4)")
  elseif ("${_LAPACK_PREFIX}" STREQUAL "Z")
    set (_DATATYPE "complex (8)")
  elseif ("${_LAPACK_PREFIX}" STREQUAL "C")
    set (_DATATYPE "complex (4)")
  endif ()
  
  message (STATUS "Detected LAPACK routine ${_LAPACK_ROUTINE}")

  file (APPEND "${_FORTRAN_MODULE_FILENAME}" 
"  subroutine ${_LAPACK_PREFIX}LARFP_wrapper( n, alpha, x, incx, tau )
     integer, intent(in)     :: n, incx
     ${_DATATYPE}, intent(inout) :: alpha
     ${_DATATYPE}, intent(out)   :: tau
     ${_DATATYPE}, intent(inout) :: x(*)

     interface
       subroutine ${_LAPACK_ROUTINE}( n, alpha, x, incx, tau )
         integer, intent(in)     :: n, incx
         ${_DATATYPE}, intent(inout) :: alpha
         ${_DATATYPE}, intent(out)   :: tau
         ${_DATATYPE}, intent(inout) :: x(*)
       end subroutine ${_LAPACK_ROUTINE}
     end interface

     call ${_LAPACK_ROUTINE}( n, alpha, x, incx, tau )
  end subroutine ${_LAPACK_PREFIX}LARFP_wrapper

")
endforeach ()
file (APPEND "${_FORTRAN_MODULE_FILENAME}" "end module TsqrHouseholderReflector")

