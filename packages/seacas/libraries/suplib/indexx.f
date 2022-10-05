C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

c$$$      program test
c$$$      parameter (n=20)
c$$$      dimension a(n)
c$$$      integer indx(n)
c$$$
c$$$      integer*4 seed
c$$$      seed=0
c$$$
c$$$      call srand(time())
c$$$      do i=1, n
c$$$        a(i) = rand(seed)
c$$$      end do
c$$$
c$$$      do i=1,n
c$$$        write (*,*) a(i)
c$$$      end do
c$$$      write (*,*) '-------'
c$$$
c$$$      call indexx(a, indx, n, .true.)
c$$$
c$$$      do i=1,n
c$$$        write (*,*) i, indx(i), a(indx(i))
c$$$      end do
c$$$
c$$$      end

C------------------------------------------------------------------------
C     SUBROUTINE INDEXX: Indexes an array ARRAY, that is
C           it outputs an array INDX such that ARRAY(INDX(J)) is in
C           ascending order for J=1,2,...,N.  The input quantities N and
C           ARRAY are not changed.

C     ARRAY (*)        -  Array to be sorted
C     INDX  (modified) -  Sorted order of ARRAY
C     N                -  Number of elements in ARRAY
C     INIT             -  .FALSE. if INDX already setup
C                         .TRUE.  if INDX must be initialized
C------------------------------------------------------------------------

      subroutine indexx(a, indx, n, init)

      dimension a(*)
      integer indx(0:*)
      integer n
      logical init

      integer start, bottom
      integer temp

      if (init) then
        do i=1, n
            indx(i-1) = i
          end do
      end if

      do start = (n-2)/2, 0, -1
        call my_siftdown(a, indx, start, n)
      end do

      do bottom = n-1, 1, -1
        temp = indx(0)
        indx(0) = indx(bottom)
        indx(bottom) = temp
        call my_siftdown(a, indx, 0, bottom)
      end do
      end

      subroutine my_siftdown(a, indx, start, bottom)

      real a(*)
      integer indx(0:*)

      integer start, bottom
      integer child, root
      integer temp

      root = start
      do while(root*2 + 1 .lt. bottom)
        child = root * 2 + 1

        if ((child + 1 .lt. bottom) .and.
     *    (a(indx(child)) .lt. a(indx(child+1)))) then
          child = child + 1
        end if

        if (a(indx(root)) .lt. a(indx(child))) then
          temp = indx(child)
          indx(child) = indx(root)
          indx(root) = temp
          root = child
        else
          return
        end if
      end do
      return
      end

