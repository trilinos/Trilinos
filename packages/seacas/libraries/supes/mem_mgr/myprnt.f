C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE MYPRNT (NAME1, UNIT, MYCV, OFFSET, TOFFST,
     *   DICT, DPOINT, LDICT, NNAMES, CHRNUM,
     *   CHRCOL, NCOLP, WRDSIZ, LASTER)

      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'

C***********************************************************************

C     NAME1    Name of array to be printed
               CHARACTER*8 NAME1
C     UNIT     Output unit number.
C     MYCV     Internal character array
               CHARACTER*1 MYCV(*)
C     OFFSET   Offset between internal reference and users reference
C              string.
C     TOFFST   Offset between internal reference and internal character
C              array.
C     DICT     Dictionary name table
C     DPOINT   Dictionary pointer table
C     LDICT    Dimension of dictionary
C     NNAMES   Number of names in the dictionary
               CHARACTER*8 DICT(LDICT,CHRCOL)
               DIMENSION DPOINT(LDICT,CHRCOL,3), NNAMES(2)
C     CHRCOL   Number of column for character names.
C     NCOLP    Number of print columns
C     WRDSIZ   Number of characters to group together in printing.
C     LASTER   Error return

C***********************************************************************

C     Check worklength

      IF (WRDSIZ .LT. 1 .OR. WRDSIZ+2+11 .GT. NCOLP) THEN
         LASTER = BADLEN
         RETURN
      END IF
C     FIND NAME1 IN DICTIONARY

      CALL MYFIND (NAME1, DICT, DPOINT, LDICT, NNAMES,
     *   CHRCOL, LASTER, ROW)
      IF (LASTER .NE. SUCESS) RETURN
      DELTA = (DPOINT(ROW,CHRCOL,1) - 1) * CHRNUM + 1 + OFFSET
      WRITE(UNIT,10000) DICT(ROW,CHRCOL),
     *   DELTA,
     *   DPOINT(ROW,CHRCOL,3)
      IF (DPOINT(ROW,CHRCOL,2) .LT. 0) THEN
         LASTER = DEFRON
         WRITE (UNIT, *) 'THIS VECTOR WAS RESERVED IN THE DEFERRED '//
     *      'MODE AND IS NOT YET RESOLVED.'
         RETURN
      END IF

      DELTA = (DPOINT(ROW,CHRCOL,1) - 1) * CHRNUM + TOFFST
      NCOL = (NCOLP - 11) / (WRDSIZ + 2)
      NROW = (DPOINT(ROW,CHRCOL,3) + WRDSIZ * NCOL - 1)
     *   / (WRDSIZ * NCOL)
      NGRP = (NROW + 9) / 10
      DO 110 IGRP = 1, NGRP
         WRITE(UNIT,10010)
         NPRT = (IGRP - 1) * 10 * NCOL * WRDSIZ
         NREM = DPOINT(ROW,CHRCOL,3) - NPRT
         NROW = (NREM + WRDSIZ * NCOL - 1) / (NCOL * WRDSIZ)
         NROW = MIN(10, NROW)
         J = NPRT + 1
         DO 100 IROW = 1, NROW
            WRITE (UNIT, 10020) J,
     *         ((MYCV(K), K=J+DELTA+WRDSIZ*(IWRD-1),
     *         MIN(DELTA+DPOINT(ROW,CHRCOL,3),
     *         J+DELTA+WRDSIZ*IWRD-1)),
     *         ' ', ' ', IWRD = 1, NCOL)
            J = J + NCOL * WRDSIZ
  100    CONTINUE
  110 CONTINUE

      LASTER = SUCESS
      RETURN
10000 FORMAT('0'/'0ARRAY NAME = ',A,3X,'LOCATION = ',I16,3X,
     *  'LENGTH = ',I8)
10010 FORMAT(' ')
10020 FORMAT(1X,I6,':',132A1)
      END
