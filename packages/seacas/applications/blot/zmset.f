C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ZMSET (XZMMIN, XZMMAX, YZMMIN, YZMMAX, XN, YN,
     &   LINSET, IPSET, NPART)
C=======================================================================

C   --*** ZMSET *** (MESH) Delete partial lines outside zoom window
C   --   Written by Amy Gilkey - revised 10/02/86
C   --
C   --ZMSET deletes partial lines with both nodes outside the zoom window,
C   --and the line not crossing the zoom window.
C   --
C   --Parameters:
C   --   XZMMIN, XZMMAX, YZMMIN, YZMMAX - IN - the box enclosing the
C   --      zoom window
C   --   XN, YN - IN - the nodal coordinates
C   --   LINSET - IN/OUT - the sorted line set
C   --   IPSET - IN/OUT - the indices of the partial line set
C   --   NPART - IN/OUT - the number of lines in the partial line set
C   --
C   --Common Variables:
C   --   Uses LLNSET of /D3NUMS/

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'd3nums.blk'

      REAL XN(*), YN(*)
      INTEGER LINSET(LLNSET,*)
      INTEGER IPSET(*)

      YLN = YZMMAX - YZMMIN

      nhid = 0
      IP = 1
  100 CONTINUE
      IF (IP .LE. NPART) THEN

C      --Delete partial line if both nodes are outside zoom window

         N1 = LINSET(1,IPSET(IP))
         X1 = XN(N1)
         Y1 = YN(N1)
         IF ((X1 .LE. XZMMIN) .OR. (X1 .GE. XZMMAX)
     &      .OR. (Y1 .LE. YZMMIN) .OR. (Y1 .GE. YZMMAX)) THEN
            N2 = LINSET(2,IPSET(IP))
            X2 = XN(N2)
            Y2 = YN(N2)
            IF ((X2 .LE. XZMMIN) .OR. (X2 .GE. XZMMAX)
     &         .OR. (Y2 .LE. YZMMIN) .OR. (Y2 .GE. YZMMAX)) THEN
               IF (((X1 .LE. XZMMIN) .AND. (X2 .LE. XZMMIN)) .OR.
     &            ((X1 .GE. XZMMAX) .AND. (X2 .GE. XZMMAX)) .OR.
     &            ((Y1 .LE. YZMMIN) .AND. (Y2 .LE. YZMMIN)) .OR.
     &            ((Y1 .GE. YZMMAX) .AND. (Y2 .GE. YZMMAX))) THEN
                  IPSET(IP) = IPSET(NPART)
                  NPART = NPART - 1
                  IP = IP - 1
                  nhid = nhid + 1

               ELSE
C               --Calculate the intersection of the window and the partial line
C               --Solve the simultaneous equations:
C               --   X = X2 + (X1 - X2) * TVH = XLIM
C               --   Y = Y2 + (Y1 - Y2) * TVH = YZMMIN + (YZMMAX - YZMMIN) * TLN

                  IF (X1 .LE. XZMMIN) THEN
                     XLIM = XZMMIN
                  ELSE
                     XLIM = XZMMAX
                  END IF

                  XVH = X1 - X2
                  YVH = Y1 - Y2
                  XLH = XLIM - X2
                  YLH = YZMMIN - Y2
                  IF ((XVH .NE. 0.0) .AND. (YLN .NE. 0.0)) THEN
                     TLN = - (-YVH * XLH + XVH * YLH) / (XVH * YLN)
                     IF ((TLN .LE. 0) .OR. (TLN .GE. 1)) THEN
                        IPSET(IP) = IPSET(NPART)
                        NPART = NPART - 1
                        IP = IP - 1
                        nhid = nhid + 1
                     END IF
                  END IF
               END IF
            END IF
         END IF

         IP = IP + 1
         GOTO 100
      END IF
      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &   write (*, '(1x,a,i5)') 'partial lines outside zoom =', nhid

      RETURN
      END
