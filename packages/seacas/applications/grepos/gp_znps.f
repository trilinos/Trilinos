C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ZNPS (IDNPS, NNNPS, IXNNPS, LTNNPS, FACNPS, NUMNPS,
     &   LNPSNL, MAXVAL)
C=======================================================================
C   --   NUMNPS - IN - the number of node sets
C   --   LNPSNL - IN - the length of the node sets node list
C   --   IDNPS - OUT - the node set ID for each set
C   --   NNNPS - OUT - the number of nodes for each set
C   --   IXNNPS - OUT - the index of the first node for each set
C   --   LTNNPS - OUT - the nodes for all sets
C   --   FACNPS - OUT - the distribution factors for all sets

      INTEGER NUMNPS, LNPSNL
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER IXNNPS(*)
      INTEGER LTNNPS(*)
      REAL FACNPS(*)

      IOFF = 0
      ISET = 0
      DO 20 I=1, NUMNPS
         IBEG = IXNNPS(I)
         IEND = IBEG + NNNPS(I) - 1
         IF (IDNPS(I) .NE. 0) THEN
C ... This set may have been partially or totally deleted
C ... Count number of deleted nodes (node number = MAXVAL)
            NUMDEL = NUMEQI (MAXVAL, NNNPS(I), LTNNPS(IBEG))
            IF (NUMDEL .EQ. NNNPS(I)) THEN
C ... This set has been totally deleted
               IDNPS(I) = 0
            ELSE
C ... This set has been partially deleted, NNNPS(I)-NUMDEL nodes left
               ISET = ISET + 1
               IDNPS(ISET)  = IDNPS(I)
               NNNPS(ISET)  = NNNPS(I) - NUMDEL
               IXNNPS(ISET) = IOFF + 1
               DO 10 INOD = IBEG, IEND
                  IF (LTNNPS(INOD) .LT. MAXVAL) THEN
                     IOFF = IOFF + 1
                     LTNNPS(IOFF) = LTNNPS(INOD)
                     FACNPS(IOFF) = FACNPS(INOD)
                  END IF
   10          CONTINUE
            END IF
         END IF
         IF (IDNPS(I) .EQ. 0) THEN
C ... This set has been totally deleted
         END IF
   20 CONTINUE

      NUMNPS = ISET
      LNPSNL = IOFF

      RETURN
      END
