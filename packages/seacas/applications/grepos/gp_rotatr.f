C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ROTATR (NELBLK, NDIM, IDELB, BLKTYP, NUMATR,
     *  NUMELB, ATRIB)
C=======================================================================

      include 'exodusII.inc'
      include 'gp_attrot.blk'

      CHARACTER*(MXSTLN) BLKTYP(*)
      INTEGER NUMATR(*)
      INTEGER NUMELB(*)
      INTEGER IDELB(*)
      REAL ATRIB(*)

      IEATR = 0
      IAT   = 1

      DO IELB = 1, NELBLK
        ISATR = IEATR + 1
        IEATR = IEATR + NUMATR(IELB) * NUMELB(IELB)
        if (numatr(ielb) .ge. ATTIND+NDIM-1) THEN
          if (rotall .or.
     *      (.NOT. ROTTYP .AND. IDELB(IELB) .EQ. ATTBLK)  .OR.
     *      (ROTTYP .AND. BLKTYP(IELB) .EQ. ROTBLK)) then
            CALL ROTAT1 (NDIM, NUMELB(IELB), NUMATR(IELB), ATTIND,
     *        ATRIB(ISATR), ROTATT)
          end if
        end if
        IAT = IAT + NUMATR(IELB)
      END DO

      RETURN
      END

      SUBROUTINE ROTAT1(NDIM, NUMEL, NUMATR, IDX, ATRIB, ROTMAT)
      REAL ATRIB(*)
      REAL ROTMAT(3,3)

      IBEG = 0
      IF (NDIM .EQ. 3) THEN
        DO IEL = 1, NUMEL
          X = ATRIB(IDX+0 + IBEG)
          Y = ATRIB(IDX+1 + IBEG)
          Z = ATRIB(IDX+2 + IBEG)
          XN = X*ROTMAT(1,1) + Y*ROTMAT(2,1) + Z*ROTMAT(3,1)
          YN = X*ROTMAT(1,2) + Y*ROTMAT(2,2) + Z*ROTMAT(3,2)
          ZN = X*ROTMAT(1,3) + Y*ROTMAT(2,3) + Z*ROTMAT(3,3)
          ATRIB(IDX+0 + IBEG) = XN
          ATRIB(IDX+1 + IBEG) = YN
          ATRIB(IDX+2 + IBEG) = ZN
          IBEG = IBEG + NUMATR
        END DO
      ELSE IF (NDIM .EQ. 2) THEN
        DO IEL = 1, NUMEL
          X = ATRIB(IDX+0 + IBEG)
          Y = ATRIB(IDX+1 + IBEG)
          XN = X*ROTMAT(1,1) + Y*ROTMAT(2,1)
          YN = X*ROTMAT(1,2) + Y*ROTMAT(2,2)
          ATRIB(IDX+0 + IBEG) = XN
          ATRIB(IDX+1 + IBEG) = YN
          IBEG = IBEG + NUMATR
        END DO
      END IF

      RETURN
      END

