!-----------------------------------------------------------------------------!
! \file   kba_equations/Parameters.f90
! \author Thomas M. Evans
! \date   Thur Jan 24 13:59:01 2008
! \brief  Parameters for FORTRAN kernels in KBA.
! \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC
!-----------------------------------------------------------------------------!
! $Id: Parameters.f90,v 1.2 2009/07/23 18:45:16 9te Exp $
!-----------------------------------------------------------------------------!

module Parameters

  implicit none

  ! Single precision control.
  integer, parameter :: SHORT  = selected_int_kind(3)
  integer, parameter :: SINGLE = selected_real_kind(5)
  integer, parameter :: DOUBLE = selected_real_kind(14)

  ! Indices
  integer, parameter :: I = 1
  integer, parameter :: J = 2
  integer, parameter :: K = 3

  ! Fractions.
  real(DOUBLE), parameter :: zero        = 0.0_DOUBLE
  real(DOUBLE), parameter :: one         = 1.0_DOUBLE
  real(DOUBLE), parameter :: two         = 2.0_DOUBLE
  real(DOUBLE), parameter :: three       = 3.0_DOUBLE
  real(DOUBLE), parameter :: four        = 4.0_DOUBLE
  real(DOUBLE), parameter :: one_half    = one / 2.0_DOUBLE
  real(DOUBLE), parameter :: one_third   = one / 3.0_DOUBLE
  real(DOUBLE), parameter :: one_fourth  = one / 4.0_DOUBLE
  real(DOUBLE), parameter :: one_fifth   = one / 5.0_DOUBLE
  real(DOUBLE), parameter :: one_sixth   = one / 6.0_DOUBLE
  real(DOUBLE), parameter :: one_seventh = one / 7.0_DOUBLE
  real(DOUBLE), parameter :: one_eighth  = one / 8.0_DOUBLE

  ! Constants.
  real(DOUBLE), parameter :: pi          = 3.1415926535897931159979635_DOUBLE
  real(DOUBLE), parameter :: two_pi      = two * pi
  real(DOUBLE), parameter :: four_pi     = four * pi
  real(DOUBLE), parameter :: inv_two_pi  = one / two_pi
  real(DOUBLE), parameter :: inv_four_pi = one / four_pi

  ! Huge constant.
  real(DOUBLE), parameter :: huge = 1.0e14

  ! Size of mesh in each direction.
  integer, save :: im, jm, km

end module Parameters

!-----------------------------------------------------------------------------!
! Set the dimensionality of the mesh in each direction.

subroutine set_dimensions(im_in, jm_in, km_in)
  
  use Parameters, only : im, jm, km
  implicit none

  integer, intent(in) :: im_in, jm_in, km_in

  im = im_in
  jm = jm_in
  km = km_in

end subroutine set_dimensions

!-----------------------------------------------------------------------------!
!                           end of Parameters.f90
!-----------------------------------------------------------------------------!
