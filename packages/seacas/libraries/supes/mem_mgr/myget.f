C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE MYGET (MYLOC, MNGET, VOID, LVOID, NVOIDS,
     *   CHRCOL, MAXSIZ, LASTER, VROW)

      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'

C     This subroutine returns the location (row number) of a void with
C     sufficient space for the memory request.  If necessary, memory is
C     requested from the system.  The memory is contiguous.
C     This routine is to be used only if CHRCOL = 2.

C***********************************************************************

C     MYLOC    Address of internal reference array
C     MNGET    Memory request in character storage units
C     VOID     Void table
C     LVOID    Dimension of void table
C     NVOIDS   Number of voids
               DIMENSION VOID(LVOID,CHRCOL,2), NVOIDS(2)
C     CHRCOL   Column for character tables (must be 2)
C     MAXSIZ   Dimension of static character array.
C     LASTER   Error return
C     VROW     Row number of void which satisfies the memory request

C***********************************************************************

C     IS THE MEMORY REQUEST SENSIBLE?

      IF (MNGET .LT. 0) THEN
         LASTER = BADLEN
         RETURN
      ELSE IF (MNGET .EQ. 0) THEN
         LASTER = SUCESS
         RETURN
      END IF

      CALL MXLOOK (MNGET, VOID(1,CHRCOL,1), CHRCOL*LVOID,
     *   NVOIDS(CHRCOL), VROW, LASTER)
      IF (LASTER .EQ. SUCESS) RETURN

C     CALL EXTENSION LIBRARY ROUTINE TO GET SPACE FROM SYSTEM.

      CALL MYMEMY (MNGET, LOC, MEMRET, MAXSIZ)
      LOC = LOC - MYLOC + 1

      IF (MEMRET .LT. 0) THEN

C        ILLEGAL MEMORY BLOCK SIZE.

         LASTER = ILBLK
         RETURN

      END IF

C     UPDATE VOID TABLE.

      CALL VTABLE (LOC, MEMRET, VOID(1,CHRCOL,1), LVOID,
     *   NVOIDS(CHRCOL), CHRCOL, LASTER)
      IF (LASTER .NE. SUCESS) RETURN

      CALL MXLOOK (MNGET, VOID(1,CHRCOL,1), CHRCOL*LVOID,
     *   NVOIDS(CHRCOL), VROW, LASTER)

      RETURN
      END
