C Copyright(C) 1999-2020, 2025 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBPELB (OPTION, NELBLK, IDELB, NUMELB, NUMLNK, NUMATR,
     &   BLKTYP, EBNAME, ATNAME, NVAREL, NAMEEV, ISEVOK, LISEV, STAT)
C=======================================================================

C   --*** DBPELB *** (EXOLIB) Print database element block summary
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

      include 'gp_params.blk'
      include 'gp_namlen.blk'
      CHARACTER*(*) OPTION
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      INTEGER STAT(*)
      CHARACTER*(MXSTLN) BLKTYP(*)
      CHARACTER*(maxnam) EBNAME(*)
      CHARACTER*(maxnam) NAMEEV(*)
      CHARACTER*(maxnam) ATNAME(*)
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

      isatn = 1
      DO IELB = 1, NELBLK
         IF (ISABRT ()) RETURN

         WRITE (STRA, 10000, IOSTAT=IDUM) IELB
         CALL PCKSTR (1, STRA)
            if (stat(ielb) .lt. 0) then
               WRITE (*, 10010, IOSTAT=IDUM) IDELB(IELB), STRA(:LSTRA),
     &              NUMELB(IELB), NUMLNK(IELB), NUMATR(IELB), 'DELETED'
            else
               WRITE (*, 10010, IOSTAT=IDUM) IDELB(IELB), STRA(:LSTRA),
     &              NUMELB(IELB), NUMLNK(IELB), NUMATR(IELB), ''
            end if
         IF (DONAM) THEN
               WRITE (*, 10020) EBNAME(IELB)(:LENSTR(EBNAME(IELB))),
     $              BLKTYP(IELB)(:LENSTR(BLKTYP(IELB)))
         END IF

         if (numatr(ielb) .gt. 0) then
           mxnam = 0
           DO I = 1, numatr(ielb)
             LNAM = lenstr(atname(isatn+i-1))
             mxnam = max(lnam, mxnam)
           END DO
           if (mxnam .gt. 1) then
             WRITE (*, 10090, IOSTAT=IDUM)
     &         (ATNAME(I)(:mxnam),
     *         I=ISATN,ISATN+NUMATR(IELB)-1)
           end if
           isatn = isatn + numatr(ielb)
         end if

         IF (DOVTBL) THEN
            NSEL = 0
            DO I = 1, NVAREL
               IF (ISEVOK( (IELB-1)*NELBLK+I )) THEN
                  NSEL = NSEL + 1
                  LISEV(NSEL) = I
               END IF
            END DO

            WRITE (*, 10030, IOSTAT=IDUM) (NAMEEV(LISEV(I)), I=1,NSEL)
         END IF
      END DO

      RETURN

10010  FORMAT (1X, 'Block', I11, 1X, A, ':',
     &   I10, ' elements',
     &   I10, '-node', I8, ' attributes',2X,A)
10020  FORMAT (4X, 'Element block name = "',A,
     $      '", Element type = "', A, '"')
10030  FORMAT (4X, 'Defined variables:', :, 3X, 4 (2X, A), :, /
     &   (4X, 1X, 6 (2X, A)))
10090 FORMAT (4X, 'Attributes: ', 10(2X, A))
      END
