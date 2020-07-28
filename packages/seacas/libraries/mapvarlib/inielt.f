C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE INIELT(SOLEB,IVAR,TIMES,ISTP,IDBLK,CENTER,DUME)

C *********************************************************************

C INIELT initializes element variable values based on TIME, ELEMENT
C BLOCK, VARIABLE NAME, COORDINATE, etc. By default, element variable
C values are set  to zero. It is intended that the user rewrite this
C subroutine to provide values that are appropriate to the problem
C being solved. This is the preferred method to handle element variable
C assignment for recipient mesh nodes that lie outside the boundary
C of the donor mesh.

C Called by INTRPE, SINTPE, TRANAB, STRAN

C *********************************************************************

C SOLEB   REAL  Array of element variable values
C               (1:numebb,1:nvarel)
C TIMES   REAL  Array of times (1:ntimes)
C IDBLK   INT   The element block I. D.
C CENTER  REAL  Array of element centroid coordinates
C               (1;numebb,1:3)

C *********************************************************************

      include 'exodusII.inc'

      include 'aexds1.blk'
      include 'ebbyeb.blk'
      include 'inival.blk'

      DIMENSION SOLEB(NUMEBB,NVAREL), TIMES(*), CENTER(NUMEBB,*)
      DIMENSION DUME(*)

C *********************************************************************

C Code to help you find some potentially useful stuff
C The actual time (real number)
C     TIME = TIMES(ISTP)

C The pointer into VARNAM to get the variable name being processed
C     INAM = IVAR + NVARGP

C The name of the variable (character) being processed
C     NAME = NAMVAR(INAM)

C The coordinates of the point (element centroid)

C XP = CENTER(IELT,1)
C YP = CENTER(IELT,2)
C ZP = CENTER(IELT,3)

C By default, set value to 0.
C User to replace this with whatever code he wishes.

      DO 10 IELT = 1, NUMEBB
        SOLEB(IELT,IVAR) = VALINI
   10 CONTINUE

      RETURN
      END
