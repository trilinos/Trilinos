C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE BLDPTE(X,Y,Z,ICON,CENTER)

C***********************************************************************

C BLDPTE LOADS THE ELEMENT CENTROIDS INTO THE XYZSRF ARRAY

C Calls subroutine CNTR

C Called by MAPVAR

C***********************************************************************

      include 'bmesh.blk'
      include 'ebbyeb.blk'
      include 'tapes.blk'

      DIMENSION X(*),Y(*),Z(*),CENTER(NUMEBB,*),ICON(NELNDB,*)
      DIMENSION XX(27), YY(27), ZZ(27)

C      NNODES = NNELM(ITYPE)
      NNODES = NELNDB
      IF (ITYPE .EQ. 6) NNODES = 4
      DO 10 IEL = 1, NUMEBB
        DO 20 INOD = 1, NNODES
          XX(INOD) = X(ICON(INOD,IEL))
          YY(INOD) = Y(ICON(INOD,IEL))
          IF (NDIMB .EQ. 3)THEN
            ZZ(INOD) = Z(ICON(INOD,IEL))
          ELSE
            ZZ(INOD) = 0.
          END IF
   20   CONTINUE
        CALL CNTR(ITYPE,XX,YY,ZZ,CENTER(IEL,1),CENTER(IEL,2),
     &            CENTER(IEL,3))
   10 CONTINUE
      RETURN
      END
