C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE DBPELB (OPTION, NELBLK, IDELB, NUMELB, NUMLNK, NUMATR,
     &   BLKTYP, NVAREL, NAMEEV, ISEVOK, LISEV)
C=======================================================================

C   --*** DBPELB *** (EXOLIB) Print database element block summary
C   --   Written by Amy Gilkey - revised 01/19/88
C   --
C   --DBPELB displays the database element block summary (block ID,
C   --number of elements, element block name).
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --      ' ' to print block summary
C   --      'N' to print element block name
C   --      'V' to print element variable truth table
C   --   NELBLK - IN - the number of element blocks
C   --   IDELB - IN - the element block ID for each block
C   --   NUMELB - IN - the number of elements for each block
C   --   NUMLNK - IN - the number of nodes per element for each block
C   --   NUMATR - IN - the number of attributes for each block
C   --   BLKTYP - IN - the names of the element block types
C   --   NVAREL - IN - the number of element variables
C   --   NAMEEV - IN - the names of the element variables
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   LISEV - SCRATCH - size = NVAREL (if 'V' in OPTION)

      include 'exodusII.inc'
      include 'ag_namlen.blk'

      CHARACTER*(*) OPTION
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      CHARACTER*(MXSTLN) BLKTYP(*)
      CHARACTER*(namlen) NAMEEV(*)
      LOGICAL ISEVOK(*)
      INTEGER LISEV(*)

      LOGICAL ISABRT
      LOGICAL DONAM, DOVTBL
      CHARACTER*20 STRA

      DONAM  = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0))
      DOVTBL = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'V') .GT. 0))

      WRITE (*, *)

      WRITE (STRA, 10000, IOSTAT=IDUM) NELBLK
10000  FORMAT ('(#', I5, ')')
      CALL PCKSTR (1, STRA)
      LSTRA = LENSTR (STRA)

      DO 110 IELB = 1, NELBLK
         IF (ISABRT ()) RETURN

         WRITE (STRA, 10000, IOSTAT=IDUM) IELB
         CALL PCKSTR (1, STRA)
         WRITE (*, 10010, IOSTAT=IDUM) IDELB(IELB), STRA(:LSTRA),
     &      NUMELB(IELB), NUMLNK(IELB), NUMATR(IELB)

         IF (DONAM) THEN
            WRITE (*, 10020) BLKTYP(IELB)
         END IF

         IF (DOVTBL) THEN
            NSEL = 0
            DO 100 I = 1, NVAREL
               IF (ISEVOK( (IELB-1)*NELBLK+I )) THEN
                  NSEL = NSEL + 1
                  LISEV(NSEL) = I
               END IF
  100       CONTINUE

            WRITE (*, 10030, IOSTAT=IDUM) (NAMEEV(LISEV(I)), I=1,NSEL)
         END IF
  110 CONTINUE

      RETURN

10010  FORMAT (1X, 'Block', I8, 1X, A, ':',
     &   I10, ' elements',
     &   I10, '-node', I8, ' attributes')
10020  FORMAT (4X, 'Element block type = "', A, '"')
10030  FORMAT (4X, 'Defined variables:', :, 3X, 4 (2X, A), :, /
     &   (4X, 1X, 6 (2X, A)))
      END
