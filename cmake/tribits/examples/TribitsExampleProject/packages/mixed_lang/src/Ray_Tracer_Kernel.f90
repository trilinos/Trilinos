!-----------------------------------------------------------------------------!
! \file   kba/Block_Tracer.f90
! \author Thomas M. Evans
! \date   Thur Jan 24 14:00 2008
! \brief  FORTRAN kernel for tracing a ray through a block of mesh.
! \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC
!-----------------------------------------------------------------------------!
! $Id: Ray_Tracer_Kernel.f90,v 1.1 2008/02/07 21:27:48 9te Exp $
!-----------------------------------------------------------------------------!

module Ray_Tracer

  use Parameters, only : DOUBLE
  implicit none
  public

  ! Distance-to-surface arrays used by the ray tracer.
  real(DOUBLE), dimension(:), allocatable, save :: dx, dy, dz

  ! Ray direction.
  real(DOUBLE), dimension(3), save :: omega

contains
  
  !---------------------------------------------------------------------------!
  ! Calculate distance to boundary arrays
  !---------------------------------------------------------------------------!
  
  subroutine calc_distances(start, edges, r, p, angle, d)

    implicit none
    
    ! >>> INPUT/OUTPUT DATA
    
    ! Starting cell position in edge array.
    integer, intent(in) :: start

    ! Edges in the chosen direction.
    real(DOUBLE), dimension(:), intent(in) :: edges

    ! Starting point position in the chosen direction.
    real(DOUBLE), intent(in) :: r

    ! Target (ending) point position in the chosen direction.
    real(DOUBLE), intent(in) :: p

    ! Inverse of the angle-of-flight.
    real(DOUBLE), intent(in) :: angle

    ! Edges in the chosen direction.
    real(DOUBLE), dimension(:), intent(out) :: d

    ! >>> LOCAL DATA
    
    ! Indices.
    integer :: n, index

    ! Number of edges.
    integer :: length

    ! Inverse angle.
    real(DOUBLE) :: inv_angle
    
    ! >>> BODY
    
    ! intialize d
    d = 1.0e14

    ! calculate the size of the edges array
    length = size(edges)

    ! initialize starting index
    n     = start
    index = 1

    ! step through edges in the positive direction
    if (angle > 0.0) then

       ! calculate inverse of angle
       inv_angle = 1.0_DOUBLE / angle

       ! continue until we hit the point or we get to the edge of the mesh
       do while (edges(n) < p .and. n < length)
          
          ! assign distance to boundary
          d(index) = (edges(n + 1) - r) * inv_angle

          ! update the indices
          n     = n + 1
          index = index + 1

       end do
    
    else if (angle < 0.0) then

       ! calculate inverse of angle
       inv_angle = 1.0_DOUBLE / angle

       ! continue until we hit the point or we get to the edge of the mesh
       do while (edges(n + 1) > p .and. n > 0)
          
          ! assign distance to boundary
          d(index) = (edges(n) - r) * inv_angle

          ! update the indices
          n     = n - 1
          index = index + 1

       end do
   
    end if
       
  end subroutine calc_distances
  
  !---------------------------------------------------------------------------!
  ! CONSTRUCT
  !---------------------------------------------------------------------------!
  
  subroutine construct_ray_tracer

    use Parameters, only : im, jm, km
    implicit none

    if (.not. allocated(dx)) allocate(dx(im + 1))
    if (.not. allocated(dy)) allocate(dy(jm + 1))
    if (.not. allocated(dz)) allocate(dz(km + 1))

  end subroutine construct_ray_tracer

  !---------------------------------------------------------------------------!
  ! DESTRUCT
  !---------------------------------------------------------------------------!

  subroutine destruct_ray_tracer

    implicit none

    if (allocated(dx)) deallocate(dx)
    if (allocated(dy)) deallocate(dy)
    if (allocated(dz)) deallocate(dz)

  end subroutine destruct_ray_tracer

end module Ray_Tracer

!-----------------------------------------------------------------------------!
! Setup routines for ray-tracing
!-----------------------------------------------------------------------------!

subroutine setup_ray_tracing

  use Ray_Tracer, only : construct_ray_tracer

  call construct_ray_tracer

end subroutine setup_ray_tracing

!-----------------------------------------------------------------------------!

subroutine cleanup_ray_tracing
  
  use Ray_Tracer, only : destruct_ray_tracer
  
  call destruct_ray_tracer

end subroutine cleanup_ray_tracing

!-----------------------------------------------------------------------------!
! Ray-tracing through a block mesh
!-----------------------------------------------------------------------------!

! Trace through a mesh block.

subroutine block_tracer(x, y, z, sigma, r, p, ijk, tau, terminator)
  
  use Parameters, only : SHORT, DOUBLE, im, jm, km, I, J, K, one
  use Ray_Tracer, only : dx, dy, dz, omega, calc_distances
  implicit none

  ! >>> INPUT DATA

  ! Edge arrays.
  real(DOUBLE), dimension(im + 1), intent(in) :: x
  real(DOUBLE), dimension(jm + 1), intent(in) :: y
  real(DOUBLE), dimension(km + 1), intent(in) :: z

  ! Total cross sections in each cell for all groups in (1/cm).
  real(DOUBLE), dimension(im * jm * km), intent(in) :: sigma

  ! Ray target (ending) point.
  real(DOUBLE), dimension(3), intent(in) :: p

  ! >>> IN/OUT DATA

  ! (i, j, k) position of ray in mesh.
  integer, dimension(3), intent(inout) :: ijk

  ! Ray position.
  real(DOUBLE), dimension(3), intent(inout) :: r

  ! Tau = sum_{i} sigma_i * l_i where l is the length of each step.
  real(DOUBLE), intent(inout) :: tau

  ! >>> OUT DATA

  ! Termination condition.
  !   0 - hit target point
  !   1 - exited low-X face of block
  !   2 - exited high-X face of block
  !   3 - exited low-Y face of block
  !   4 - exited high-Y face of block
  !   5 - exited low-Z face of block
  !   6 - exited high-Z face of block
  !   7 - negative distance-to-boundary
  integer, intent(out) :: terminator

  ! >>> LOCAL DATA

  ! Distance to the source point: norm = 1 / ||(r - rp)||.
  real(DOUBLE) :: dtarget, dtravel, dstep, d
  
  ! Boolean flag for ray-tracing.
  logical :: done

  ! Indices into distance arrays.
  integer :: ii, jj, kk

  ! Canonical cell index.
  integer :: cell, step_i, step_j, step_k
  
  ! >>> BODY

  ! add 1 to ijk to adjust to FORTRAN-style indexing
  ijk  = ijk + 1

  ! calculate the distance to the target and the direction of the ray
  dtarget  = sqrt(dot_product(r - p, r - p))
  omega    = (p - r) / dtarget
  
  ! calculate the distances to each surface
  call calc_distances(ijk(I), x, r(I), p(I), omega(I), dx)
  call calc_distances(ijk(J), y, r(J), p(J), omega(J), dy)
  call calc_distances(ijk(K), z, r(K), p(K), omega(K), dz)

  ! initialize indices into distance arrays
  ii = 1
  jj = 1
  kk = 1

  ! steps
  step_i = int(sign(one, omega(I)))
  step_j = int(sign(one, omega(J)))
  step_k = int(sign(one, omega(K)))

  ! initialize loop
  done    = .false.
  dtravel = 0.0_DOUBLE

  ! ray-trace to point or until we leave the block
  do while (done .eqv. .false.)
     
     ! calculate the canonical cell index
     cell = ijk(I) + (ijk(J) - 1) * im + (ijk(K) - 1) * im * jm
     
     ! get the distance
     d = min(dx(ii), min(dy(jj), dz(kk)))

     ! calculate the step
     dstep = d - dtravel

     ! push to next boundary or to target point
     if (dstep < dtarget) then

        ! subtract the distance from the target and add to the travel
        dtarget = dtarget - dstep
        dtravel = d

        ! add up tau
        tau = tau + sigma(cell) * dstep

        ! step to next face
        if (d == dx(ii)) then

           ii     = ii + 1
           ijk(I) = ijk(I) + step_i
           
           ! terminating conditions
           if (ijk(I) == 0) then
              terminator = 1
              done       = .true.
           else if (ijk(I) == im + 1) then
              terminator = 2
              done       = .true.
           end if

        else if (d == dy(jj)) then

           jj     = jj + 1
           ijk(J) = ijk(J) + step_j

           ! terminating conditions
           if (ijk(J) == 0) then
              terminator = 3
              done       = .true.
           else if (ijk(J) == jm + 1) then
              terminator = 4
              done       = .true.
           end if

        else if (d == dz(kk)) then

           kk     = kk + 1
           ijk(K) = ijk(K) + step_k

           ! terminating conditions
           if (ijk(K) == 0) then
              terminator = 5
              done       = .true.
           else if (ijk(K) == km + 1) then
              terminator = 6
              done       = .true.
           end if

        end if

     else

        ! add to the travel
        dtravel = dtravel + dtarget

        ! add up tau
        tau = tau + sigma(cell) * dtarget

        ! we are finished
        done       = .true.
        terminator = 0
        
     end if

  end do

  ! calculate final location of particle
  r = r + dtravel * omega

  ! subtract 1 to ijk to adjust to C-style indexing
  ijk = ijk - 1

end subroutine block_tracer

!-----------------------------------------------------------------------------!
!                          end of Block_Tracer.f90
!-----------------------------------------------------------------------------!
