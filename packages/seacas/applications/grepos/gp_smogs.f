C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SMOGS (X, Y, Z, numelb, nlink, idelb, link,
     $     isbnd, nelblk, NUMNP, NIT, EPS, R0,
     $     xscr, yscr, zscr, nscr, ndim)
C=======================================================================
C***********************************************************************

C  SUBROUTINE SMOGS  =  MESH SMOOTHING BY LAPLACE-S USING GAUSS-SEIDEL

C***********************************************************************

C  VARIABLES USED:
C     NIT   =  THE MAX NUMBER OF ITERATIONS TO DO.
C     EPS   =  MINIMUM DISTANCE NODES MUST MOVE TO CONTINUE ITERATIONS
C     RO    =  AN UNDER- OR OVER-RELAXATION FACTOR (NORMALLY 1.0)

C***********************************************************************

      real x(*), y(*), z(*)
      integer numelb(*), nlink(*), idelb(*), link(*)
      logical isbnd(*)
      real xscr(*), yscr(*), zscr(*)
      integer nscr(*)
      LOGICAL BIG

      IF (R0 .LT. 0.01) R0 = 1.
      EPS2 = (EPS*R0)**2

C  ITERATION LOOP

      DO 120 IT = 1, NIT
         call inirea (numnp, 0.0, xscr)
         call inirea (numnp, 0.0, yscr)
         if (ndim .eq. 3) call inirea (numnp, 0.0, zscr)
         call iniint (numnp, 0, nscr)
         ielnk = 0
         BIG = .FALSE.

C  NODE LOOP

         do 110 iblk = 1, nelblk
            islnk = ielnk+1
            ielnk = ielnk + nlink(iblk) * numelb(iblk)

C  SKIP CONTINUATION AND BOUNDARY LINES

            if (ndim .eq. 2) then
               call smogs2 (xscr, yscr, isbnd, x, y, numelb(iblk),
     $              link(islnk), nlink(iblk), numnp, nscr)
            else
               call smogs3 (xscr, yscr, zscr, isbnd, x, y, z,
     $              numelb(iblk), link(islnk), nlink(iblk),
     $              numnp, nscr)
            end if
 110     continue

         delmax = 0.0
         do 130 i=1, numnp
            if (.not. isbnd(i)) then
               SUMX = xscr(i)/dble(nscr(i))
               SUMY = yscr(i)/dble(nscr(i))
               XDEL = r0*(SUMX - X(i))
               YDEL = r0*(SUMY - Y(i))
               X(I) = X(I) + XDEL
               Y(I) = Y(I) + YDEL
               if (ndim .eq. 3) then
                  sumz = zscr(i)/dble(nscr(i))
                  zdel = r0*(sumz - z(i))
                  z(i) = z(i) + zdel
               else
                  zdel = 0.0
               end if
               del = xdel*xdel + ydel*ydel + zdel*zdel
               if (del .gt. delmax) then
                  delmax = del
                  idel = i
               end if

C  CHECK FOR CONVERGENCE

            end if

 130     continue
         IF (delmax .GT. EPS2) BIG = .TRUE.

C  IF NO SIGNIFICANT MOVEMENTS OCCURRED,  RETURN

         write (*,*) 'Iteration = ', it, sqrt(delmax)/r0, idel
         IF (.NOT.BIG) RETURN
  120 CONTINUE
      call prterr ('INFO', 'Iterations exceeded in smogs')
      RETURN
      END
