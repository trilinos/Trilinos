C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE MXGET (MYLOC, MNGET, VOID, LVOID, NVOIDS,
     *   CHRCOL, LASTER, VROW)

      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'

C     This subroutine returns the location (row number) of a void with
C     sufficient space for the memory request.  If necessary, memory is
C     acquired from the system.  The memory is contiguous.

C***********************************************************************

C     MYLOC    Address of internal reference array
C     MNGET    Memory request in numerical storage units
C     VOID     Void table
C     LVOID    Dimension of void table
C     NVOIDS   Number of voids
               DIMENSION VOID(LVOID,CHRCOL,2), NVOIDS(2)
C     CHRCOL   Column for character tables.
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

      CALL MXLOOK (MNGET, VOID, CHRCOL*LVOID, NVOIDS(1), VROW, LASTER)
      IF (LASTER .EQ. SUCESS) RETURN

C     CALL EXTENSION LIBRARY ROUTINE TO GET SPACE FROM SYSTEM.

      CALL EXMEMY (MNGET, LOC, MEMRET)
      LOC = LOC - MYLOC + 1

c  On return from exmemy, memret is set equal to -1 on an invalid
c  memory request (at least that's the plan under the new C code
c  extension library).  Therefore, I've made the change that should
c  test the appropriate condition.

      IF (MEMRET .LT. 0) THEN

C        ILLEGAL MEMORY REQUEST.

         LASTER = NOGET
         RETURN

      END IF

C     UPDATE VOID TABLE.

      CALL VTABLE (LOC, MEMRET, VOID, LVOID, NVOIDS(1), CHRCOL, LASTER)
      IF (LASTER .NE. SUCESS) RETURN

      CALL MXLOOK (MNGET, VOID, CHRCOL*LVOID, NVOIDS(1), VROW, LASTER)

      RETURN
      END
