C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE MXPRNT (NAME1, UNIT, NAME2, MYV, RMYV, OFFSET,
     *   DICT, DPOINT, LDICT, NNAMES, CHRCOL, NCOLP, LASTER)

      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'

C***********************************************************************

C     NAME1    Name of array to be printed
               CHARACTER*8 NAME1
C     UNIT     Output unit number.
C     NAME2    Type of array to be printed
               CHARACTER*(*) NAME2
C     MYV      Internal integer array
               INTEGER MYV(*)
C     RMYV     Internal real array
               REAL RMYV(*)
C     DICT     Dictionary name table
C     DPOINT   Dictionary pointer table
C     LDICT    Dimension of dictionary
C     NNAMES   Number of names in the dictionary
               CHARACTER*8 DICT(LDICT,CHRCOL)
               DIMENSION DPOINT(LDICT,CHRCOL,2), NNAMES(2)
C     CHRCOL   Number of column for character names.
C     NCOLP    Number of print columns
C     LASTER   Error return

C***********************************************************************

C     FIND NAME1 IN DICTIONARY

      CALL MXFIND (NAME1, DICT, DPOINT, LDICT, NNAMES,
     *   CHRCOL, LASTER, ROW)
      IF (LASTER .NE. SUCESS) RETURN
      DELTA = DPOINT(ROW,1,1) - 1
      WRITE(UNIT,10000) DICT(ROW,1), DPOINT(ROW,1,1)+OFFSET,
     *   DPOINT(ROW,1,2)
      IF (DPOINT(ROW,1,2) .LT. 0) THEN
         LASTER = DEFRON
         WRITE (UNIT, *) 'THIS VECTOR WAS RESERVED IN THE DEFERRED '//
     *      'MODE AND IS NOT YET RESOLVED.'
         RETURN
      END IF
      IF (NAME2(1:1) .EQ. 'R') THEN

C        VECTOR IS REAL

         NCOL=(NCOLP-11)/13
         NROW=DPOINT(ROW,1,2)/NCOL+1
         NGRP=NROW/10+1
         DO 110 IGRP=1,NGRP
            WRITE(UNIT,10010)
            NPRT=(IGRP-1)*10*NCOL
            NREM=DPOINT(ROW,1,2)-NPRT
            NROW=(NREM+NCOL-1)/NCOL
            NROW=MIN0(10,NROW)
            DO 100 IROW=1,NROW
               J=NPRT+1+(IROW-1)*NCOL
               KU=MIN0(DPOINT(ROW,1,2),J+NCOL-1)
               WRITE(UNIT,10020)J,(RMYV(K),K=J+DELTA,KU+DELTA)
  100       CONTINUE
  110    CONTINUE

      ELSE IF (NAME2(1:1) .EQ. 'I') THEN

C        VECTOR IS INTEGER

         NCOL=(NCOLP-11)/10
         NROW=DPOINT(ROW,1,2)/NCOL+1
         NGRP=NROW/10+1
         DO 130 IGRP=1,NGRP
            WRITE(UNIT,10010)
            NPRT=(IGRP-1)*10*NCOL
            NREM=DPOINT(ROW,1,2)-NPRT
            NROW=(NREM+NCOL-1)/NCOL
            NROW=MIN0(10,NROW)
            DO 120 IROW=1,NROW
               J=NPRT+1+(IROW-1)*NCOL
               KU=MIN0(DPOINT(ROW,1,2),J+NCOL-1)
               WRITE(UNIT,10030)J,(MYV(K),K=J+DELTA,KU+DELTA)
  120       CONTINUE
  130    CONTINUE

      ELSE

C        TYPE IS UNKNOWN

         LASTER = BDTYPE
         RETURN

      END IF
      LASTER = SUCESS
      RETURN
10000 FORMAT(//' ARRAY NAME = ',A,3X,'LOCATION = ',I16,3X,
     *  'LENGTH = ',I8)
10010 FORMAT(' ')
10020 FORMAT(1X,I6,':',9(2X,1PE11.4))
10030 FORMAT(1X,I6,':',12(2X,I8))
      END
