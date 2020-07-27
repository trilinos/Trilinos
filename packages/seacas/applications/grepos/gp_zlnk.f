C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ZLNK (NELB, IDELB, NUMELB, NUMLNK, LINK, NEWLNK,
     $     NUMATR, ATRIB, NUMEL, MAP, MAPNOD, BLKTYP)
C=======================================================================

C   --***        *** (GREPOS)
C   --   Written by Greg Sjaardema - revised 03/02/90
C   --
C   --Parameters:
C   --   NELB - IN - the number of element blocks
C   --   IDELB - IN - the element block IDs for each block
C   --   NUMELB - IN - the number of elements in each block
C   --   NUMLNK - IN - the number of nodes per element in each block
C   --   LINK - IN - the connectivity for each block
C   --   NEWLNK - OUT - the new connectivity for each block
C   --   MAP - OUT - list of active nodes
C   --               map(i) > 0 for active node
C   --               map(i) = length of map()+1 for deleted node

      include 'gp_params.blk'
      INTEGER IDELB(NELB)
      INTEGER NUMELB(NELB)
      INTEGER NUMLNK(NELB)
      INTEGER LINK(*), NEWLNK(*)
      INTEGER NUMATR(NELB)
      REAL    ATRIB(*)
      INTEGER MAP(NUMEL), MAPNOD(*)
      CHARACTER*(MXSTLN) BLKTYP(*)

      CALL INIINT (NUMEL, 0, MAP)

      IOFF  = 0
      IOFFA = 0
      IELNK = 0
      IEATR = 0
      IEELE = 0
      DO 40 IELB = 1, NELB
         ISLNK = IELNK + 1
         IELNK = IELNK + NUMLNK(IELB) * NUMELB(IELB)
         ISATR = IEATR + 1
         IEATR = IEATR + NUMATR(IELB) * NUMELB(IELB)
         ISELE = IEELE + 1
         IEELE = IEELE + NUMELB(IELB)

         IF (IDELB(IELB) .EQ. 0) THEN
C -- deleted element block
            IOFF = IOFF + NUMLNK(IELB) * NUMELB(IELB)
            IOFFA = IOFFA + NUMATR(IELB) * NUMELB(IELB)
         ELSE
C??            IF (IOFF .NE. 0) THEN
               DO 10 I = ISLNK, IELNK
                  NEWLNK(I-IOFF) = MAPNOD(LINK(I))
   10          CONTINUE
C??            END IF
            IF (IOFFA .NE. 0) THEN
               DO 20 I = ISATR, IEATR
                  ATRIB(I-IOFFA) = ATRIB(I)
   20          CONTINUE
            END IF
            DO 30 I = ISELE, IEELE
               MAP(I) = 1
   30       CONTINUE
         END IF
   40 CONTINUE

C .. Set up map of old number to new number.  If old id = I, then
C     new id = MAP(I)
      NEWID = 0
      DO 50 IEL = 1, NUMEL
         IF (MAP(IEL) .GT. 0) THEN
            NEWID = NEWID + 1
            MAP(IEL) = NEWID
         ELSE
            MAP(IEL) = NUMEL+1
         END IF
   50 CONTINUE
      NUMEL = NEWID

C .. Update element block counter arrays
      IAELB = 0
      DO 60 IELB = 1, NELB
         IF (IDELB(IELB) .NE. 0) THEN
            IAELB = IAELB + 1
            IDELB(IAELB) = IDELB(IELB)
            NUMELB(IAELB) = NUMELB(IELB)
            NUMLNK(IAELB) = NUMLNK(IELB)
            NUMATR(IAELB) = NUMATR(IELB)
            BLKTYP(IAELB) = BLKTYP(IELB)
         END IF
   60 CONTINUE
      NELB = IAELB

      RETURN
      END
