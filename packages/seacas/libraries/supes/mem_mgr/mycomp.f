C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE MYCOMP (MYV, VOID, LVOID,
     *   NVOIDS, DPOINT, LDICT, NNAMES, CHRCOL, LASTER)

      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'

C     THIS ROUTINE PERFORMS THE NUMERIC DATA COMPRESSION OPERATION.

C************************************************************************

C     MYV      Reference array
C     VOID     Void table
C     LVOID    Dimension of VOID
C     NVOIDS   Number of voids
C     DPOINT   Dictionary pointer table
C     LDICT    Dimension of DPOINT
C     NNAMES   Number of names
C     CHRCOL   Column for character tables
C     LASTER   Error return code

      DIMENSION DPOINT(LDICT,CHRCOL,2), VOID(LVOID,CHRCOL,2)
      DIMENSION NNAMES(2), NVOIDS(2)
      CHARACTER*1 MYV(1)

C************************************************************************

      LASTER = SUCESS

C     The basic strategy is to look for an array in the dictionary
C     which is immediately preceded by a void.  If found, a shift
C     is performed, and the void table is updated.

      IVOID = 0
  100 CONTINUE
      IVOID = IVOID + 1
  110 IF (IVOID .GT. NVOIDS(2)) GO TO 130
         VADDR = VOID(IVOID,2,1) + VOID(IVOID,2,2)
         DO 120 IDICT = 1, NNAMES(2)
            DADDR = DPOINT(IDICT,2,1)
            IF (VADDR .EQ. DADDR .AND. DPOINT(IDICT,2,2) .GT. 0) THEN

C              Perform data shift and update void table.

               CALL SHFTC (MYV, LDICT,
     *            DADDR, DADDR+DPOINT(IDICT,2,2)-1, VOID(IVOID,2,2))
               DPOINT(IDICT,2,1) = VOID(IVOID,2,1)
               VOID(IVOID,2,1) = DPOINT(IDICT,2,1) + DPOINT(IDICT,2,2)
               CALL VTABLE (0, 0, VOID(1,2,1), LVOID, NVOIDS(2),
     *            CHRCOL, LASTER)
               IF (LASTER .NE. SUCESS) RETURN
               GO TO 110

            END IF
  120    CONTINUE
      GO TO 100
  130 CONTINUE
      RETURN
      END
