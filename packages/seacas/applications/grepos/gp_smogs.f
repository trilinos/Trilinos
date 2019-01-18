C Copyright(C) 2011-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C           
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C                         
C * Neither the name of NTESS nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
C                                                 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C=======================================================================
      SUBROUTINE SMOGS (X, Y, Z, numelb, nlink, idelb, link,
     $     isbnd, nelblk, NUMNP, NIT, EPS, R0,
     $     xscr, yscr, zscr, nscr, ndim)
C=======================================================================
C***********************************************************************
C
C  SUBROUTINE SMOGS  =  MESH SMOOTHING BY LAPLACE-S USING GAUSS-SEIDEL
C
C***********************************************************************
C
C  VARIABLES USED:
C     NIT   =  THE MAX NUMBER OF ITERATIONS TO DO.
C     EPS   =  MINIMUM DISTANCE NODES MUST MOVE TO CONTINUE ITERATIONS
C     RO    =  AN UNDER- OR OVER-RELAXATION FACTOR (NORMALLY 1.0)
C
C***********************************************************************
C
      real x(*), y(*), z(*)
      integer numelb(*), nlink(*), idelb(*), link(*)
      logical isbnd(*)
      real xscr(*), yscr(*), zscr(*)
      integer nscr(*)
      LOGICAL BIG
C
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
               SUMX = xscr(i)/float(nscr(i))
               SUMY = yscr(i)/float(nscr(i))
               XDEL = r0*(SUMX - X(i))
               YDEL = r0*(SUMY - Y(i))
               X(I) = X(I) + XDEL
               Y(I) = Y(I) + YDEL
               if (ndim .eq. 3) then
                  sumz = zscr(i)/float(nscr(i))
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

C
C  IF NO SIGNIFICANT MOVEMENTS OCCURRED,  RETURN
C
         write (*,*) 'Iteration = ', it, sqrt(delmax)/r0, idel
         IF (.NOT.BIG) RETURN
  120 CONTINUE
      call prterr ('INFO', 'Iterations exceeded in smogs')
      RETURN
      END
