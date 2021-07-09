C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE MYNSRT (NAME1, NEWLOC, NUMLEN, CLEN,
     *   DICT, DPOINT, LDICT,
     *   NNAMES, CHRCOL, LASTER)

      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'

C     This routine updates the dictionary with a new name (if it is new)
C     and updates the location and length tables.  The length of the
C     dictionary is checked before the new name is added.  If LASTER is
C     not returned with a value of SUCESS, the tables are unchanged.

C***********************************************************************

C     NAME1    Name to be inserted
               CHARACTER*8 NAME1
C     NEWLOC   Location of storage
C     NUMLEN   Numeric length of storage
C     CLEN     Character length of storage
C     DICT     Dictionary name table
C     DPOINT   Dictionary pointer table
C     LDICT    Dimensioned size of dictionary
C     NNAMES   Number of entries in the dictionary
               CHARACTER*8 DICT(LDICT,CHRCOL)
               DIMENSION DPOINT(LDICT,CHRCOL,3), NNAMES(2)
C     CHRCOL   Number of column for character names.
C     LASTER   Error return

C***********************************************************************

C     IS THERE ROOM IN THE DICTIONARY?

      IF (NNAMES(CHRCOL) .GE. LDICT) THEN
         LASTER = DFULL
         RETURN
      END IF

C     FIND NAME1 IN DICTIONARY

      CALL MYFIND (NAME1, DICT, DPOINT, LDICT, NNAMES,
     *   CHRCOL, LASTER, ROW)
      IF (LASTER .EQ. WRTYPE) THEN
         RETURN
      ELSE IF (LASTER .EQ. SUCESS) THEN
         LASTER = BDNAME
         RETURN
      ELSE IF (LASTER .EQ. NONAME) THEN
         LASTER = SUCESS
      END IF

C     UPDATE DICTIONARY.

      CALL SHFTC (DICT(1,CHRCOL), CHRCOL*LDICT, ROW, NNAMES(CHRCOL), -1)
      CALL SHFTI (DPOINT(1,CHRCOL,1), CHRCOL*LDICT, 3, ROW,
     *   NNAMES(CHRCOL), -1)
      NNAMES(CHRCOL) = NNAMES(CHRCOL) + 1
      DICT(ROW,CHRCOL) = NAME1
      DPOINT(ROW,CHRCOL,1) = NEWLOC
      DPOINT(ROW,CHRCOL,2) = NUMLEN
      DPOINT(ROW,CHRCOL,3) = CLEN
      RETURN
      END
