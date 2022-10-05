C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE MXNSRT (NAME1, NEWLOC, NEWLEN,
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
C     NEWLEN   Length of storage
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

      IF (NNAMES(1) .GE. LDICT) THEN
         LASTER = DFULL
         RETURN
      END IF

C     FIND NAME1 IN DICTIONARY

      CALL MXFIND (NAME1, DICT, DPOINT, LDICT, NNAMES,
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

      CALL SHFTC (DICT, CHRCOL*LDICT, ROW, NNAMES(1), -1)
      CALL SHFTI (DPOINT, CHRCOL*LDICT, 3, ROW, NNAMES(1), -1)
      NNAMES(1) = NNAMES(1) + 1
      DICT(ROW,1) = NAME1
      DPOINT(ROW,1,1) = NEWLOC
      DPOINT(ROW,1,2) = NEWLEN
      DPOINT(ROW,1,3) = -1
      RETURN
      END
