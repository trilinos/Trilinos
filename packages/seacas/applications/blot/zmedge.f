C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ZMEDGE (XZMMIN, XZMMAX, YZMMIN, YZMMAX, XN, YN,
     &   IEDSET, NEDGES)
C=======================================================================

C   --*** ZMEDGE *** (MESH) Delete edges outside zoom window
C   --   Written by Amy Gilkey - revised 10/02/86
C   --
C   --ZMEDGE deletes edges with both nodes outside the zoom window,
C   --and the line not crossing the zoom window.
C   --
C   --Parameters:
C   --   XZMMIN, XZMMAX, YZMMIN, YZMMAX - IN - the box enclosing the
C   --      zoom window
C   --   XN, YN - IN - the nodal coordinates
C   --   IEDSET - IN/OUT - the edge line set;
C   --      (0) = face defining edge; 0 to delete edge
C   --   NEDGES - IN - the number of lines in the edge set

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      REAL XN(*), YN(*)
      INTEGER IEDSET(0:2,*)

      nhid = 0
      DO 100 IEDG = 1, NEDGES
         IF (IEDSET(0,IEDG) .EQ. 0) GOTO 100

C      --Delete edge if both nodes are outside zoom window

         N1 = IEDSET(1,IEDG)
         X1 = XN(N1)
         Y1 = YN(N1)
         IF ((X1 .LE. XZMMIN) .OR. (X1 .GE. XZMMAX)
     &      .OR. (Y1 .LE. YZMMIN) .OR. (Y1 .GE. YZMMAX)) THEN
            N2 = IEDSET(2,IEDG)
            X2 = XN(N2)
            Y2 = YN(N2)
            IF (((X1 .LE. XZMMIN) .AND. (X2 .LE. XZMMIN)) .OR.
     &         ((X1 .GE. XZMMAX) .AND. (X2 .GE. XZMMAX)) .OR.
     &         ((Y1 .LE. YZMMIN) .AND. (Y2 .LE. YZMMIN)) .OR.
     &         ((Y1 .GE. YZMMAX) .AND. (Y2 .GE. YZMMAX))) THEN
               IEDSET(0,IEDG) = 0
               nhid = nhid + 1
            END IF
         END IF

  100 CONTINUE
      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &   write (*, '(1x,a,i5)') 'edges outside zoom window =', nhid

      RETURN
      END
