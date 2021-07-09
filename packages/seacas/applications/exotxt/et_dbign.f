C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBIGN (NDB, NUMESS, IDESS, NNESS, IXNESS,
     &                  LTNESS, LTNNN, IOERR)
C=======================================================================
C   --*** DBIGN *** Get node from the side sets
C   --   Written 9/10/95 for ExodusIIv2
C   --
C   --Parameters:
C   --   NDB    - IN  - the database number
C   --   NUMESS - IN  - number of side sets
C   --   IDESS  - OUT - array of side set IDS
C   --   NNESS  - OUT - array of the number of nodes for each side set
C   --   IXNESS - OUT - array of indices into LTNESS - 1st node each set
C   --   LTNESS - OUT - array of nodes for all side sets
C   --   LTNNN  - OUT _ array of number of nodes for each side in a side sets
C   --   IOERR  - OUT - error flag

      INTEGER NDB
      INTEGER NUMESS
      INTEGER IDESS(*)
      INTEGER NNESS(*)
      INTEGER IXNESS(*)
      INTEGER LTNESS(*)
      INTEGER LTNNN(*)
      INTEGER IOERR
      IOERR = 0

C     Offset into element list of current side set
      ISOFF  = 0
C     Node count for current side set
      NODCNT = 0
      DO 100 I = 1, NUMESS
C        Set index of the first node for each side set in LTNESS
         IXNESS(I) = NODCNT + 1
         CALL EXGSP(NDB, IDESS(I), NSIDE, NDIST, IOERR)
C        NSIDE - number of sides in side set IDESS(I)
C        NDIST - number of distribution factors in side set IDESS(I)
         IF (IOERR .EQ. 1) RETURN

         CALL exgssn(NDB, IDESS(I), LTNNN(ISOFF+1),
     &               LTNESS(NODCNT+1), IOERR)
C        LTNNN(ISOFF+1) - number of nodes for each side in side set IDESS(I)
C        LTNESS(NODCNT+1) - nodes for current set
         IF (IOERR .EQ. 1) RETURN
C        Calculate node count sum for the current side set
         NCSUM = 0
         DO 90 J = 0, NSIDE-1
            NCSUM = NCSUM + LTNNN(ISOFF+1+J)
  90     CONTINUE
         NNESS(I) = NCSUM
         NODCNT = NODCNT + NCSUM
         ISOFF = ISOFF + NSIDE
 100  CONTINUE

      RETURN
      END
