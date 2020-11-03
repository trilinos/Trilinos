C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE EXPXYZ (XN, YN, ZN, ICONOD, XEXPL, YEXPL, ZEXPL,
     &   MATMAP, NELBLK, IDELB, NUMELB, NUMLNK, LINK, NUMNP, NDIM,
     $     MODE)
C=======================================================================

C   --*** EXPXYZ *** (GREPOS) Modify coordinates for each block separately
C   --   Written by Amy Gilkey - revised 05/09/88
C   --   Modified by Greg Sjaardema - 02/06/89
C   --
C   --EXPXYZ modifies the coordinate array for the database.
C   --       each block is treated separately if not connected
C   --
C   --Parameters:
C   --   XN, YN, ZN - OUT - the coordinates
C   --   IDELB - IN - the element block IDs for each block
C   --   NUMELB - IN - the number of elements in each block
C   --   NUMLNK - IN - the number of nodes per element in each block
C   --   LINK -   IN - the connectivity for each block
C   --   NUMNP - IN - Number of nodes
C   --   NDIM  - IN - Number of spatial dimensions
C   --   MODE  - 1 = Explode
C   --           2 = Scale
C   --           3 = Randomize

      CHARACTER*80 STRING
      REAL XN(*), YN(*), ZN(*)
      INTEGER ICONOD(NUMNP, NELBLK)
      REAL XEXPL(*), YEXPL(*), ZEXPL(*)
      INTEGER MATMAP(NELBLK,NELBLK)
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER LINK(*)
      INTEGER MODE

C   --Determine material block numbers of each node

      call iniint(numnp*nelblk, 0, iconod)

      IELNK = 0
      DO 10 IBLK = 1, NELBLK
         ISLNK = IELNK + 1
         IELNK = IELNK + NUMLNK(IBLK) * NUMELB(IBLK)

         CALL EXFCON (IBLK, NUMELB(IBLK), NUMLNK(IBLK), LINK(ISLNK),
     &      ICONOD(1,IBLK))
   10 CONTINUE

C ... Setup MATMAP which identifies connected material blocks.

      DO 30 IMAT = 1, NELBLK
         DO 20 IMAT2 = 1, NELBLK
            MATMAP(IMAT,IMAT2) = 0
   20    CONTINUE
         MATMAP(IMAT,IMAT) = 1
   30 CONTINUE

      DO 60 INOD = 1, NUMNP
         DO 50 IBLK = 1, NELBLK
            IF (ICONOD(INOD, IBLK) .NE. 0) THEN
               DO 40 IRST = IBLK+1, NELBLK
                  IF (ICONOD(INOD, IRST) .NE. 0) THEN
                     MATMAP(IBLK, IRST) = 1
                  END IF
   40          CONTINUE
            END IF
   50    CONTINUE
   60 CONTINUE

C ... If column J of row I of MATMAP is nonzero, then set
C     _EXPL(J) equal to _EXPL(I) and print message that material blocks
C     are connected.

      DO 80 IROW = 1, NELBLK
         DO 70 ICOL = IROW+1, NELBLK
            IF (MATMAP(IROW,ICOL) .NE. 0) THEN
               XEXPL(ICOL) = XEXPL(IROW)
               YEXPL(ICOL) = YEXPL(IROW)
               IF (NDIM .EQ. 3) THEN
                  ZEXPL(ICOL) = ZEXPL(IROW)
               END IF
               WRITE (STRING, 90) IDELB(IROW), IDELB(ICOL)
               CALL SQZSTR (STRING, LSTR)
               CALL PRTERR ('CMDSPEC', STRING(:LSTR))
            END IF
   70    CONTINUE
   80 CONTINUE
   90 FORMAT ('EXPXYZ -- Material ',I10,
     &        ' is connected to material ',I10)

C ... Collapse ICONOD array so first column is material ---
C      if multiple materials attached to one node, use first

      DO 120 INOD = 1, NUMNP
         IF (ICONOD(INOD, 1) .EQ. 0) THEN
            DO 100 IBLK = 2, NELBLK
               IF (ICONOD(INOD, IBLK) .NE. 0) THEN
                  ICONOD(INOD, 1) = IBLK
                  GO TO 110
               END IF
  100       CONTINUE
            CALL PRTERR ('ERROR', 'Node not connected to any elements')
  110       CONTINUE
         END IF
  120 CONTINUE

C ... Offset or scale each node
      IF (MODE .EQ. 1) THEN
         DO 130 INOD = 1, NUMNP
            IBLK = ICONOD(INOD,1)
            XN(INOD) = XN(INOD) + XEXPL(IBLK)
            YN(INOD) = YN(INOD) + YEXPL(IBLK)
            IF (NDIM .EQ. 3) THEN
               ZN(INOD) = ZN(INOD) + ZEXPL(IBLK)
            END IF
  130    CONTINUE
C ... Scale Each Node
      ELSE IF (MODE .EQ. 2) THEN
         DO 140 INOD = 1, NUMNP
            IBLK = ICONOD(INOD,1)
            XN(INOD) = XN(INOD) * XEXPL(IBLK)
            YN(INOD) = YN(INOD) * YEXPL(IBLK)
            IF (NDIM .EQ. 3) THEN
               ZN(INOD) = ZN(INOD) * ZEXPL(IBLK)
            END IF
  140    CONTINUE
C ... Randomize Each Node
      ELSE IF (MODE .EQ. 3) THEN
         IDUM = 1
         DO 150 INOD = 1, NUMNP
            IBLK = ICONOD(INOD,1)
            XN(INOD) = (2.0*RAN1(IDUM)-1.0) * XEXPL(IBLK) + XN(INOD)
            YN(INOD) = (2.0*RAN1(IDUM)-1.0) * YEXPL(IBLK) + YN(INOD)
            IF (NDIM .EQ. 3) THEN
               ZN(INOD) = (2.0*RAN1(IDUM)-1.0) * ZEXPL(IBLK) + ZN(INOD)
            END IF
  150    CONTINUE
      END IF
      RETURN
      END
