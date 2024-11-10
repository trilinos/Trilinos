C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MUNNPS (NUMNPS, ISTAT, LNPSNL,
     &   IDNPS, NNNPS, IXNNPS, LTNNPS, FACNPS,
     &   LTNX, FACX, IXNPS, NNX, ISCR, NODSCR,
     $   NAMNS, NAMSC, NUMNP)
C=======================================================================

C   --*** MUNNPS *** (GJOIN) Compress and rearrange nodal point sets
C   --   Written by Amy Gilkey - revised 02/25/88
C   --
C   --MUNNPS processes the nodal point sets according to the set status.
C   --Sets may be combined or deleted.
C   --
C   --Parameters:
C   --   NUMNPS - IN/OUT - the number of nodal point sets
C   --   ISTAT - IN - the status of each set:
C   --      0 = same
C   --      - = delete
C   --      n = combine with set n
C   --   LNPSNL - IN/OUT - the length of the nodal point sets node list
C   --   IDNPS - IN/OUT - the nodal point set ID for each set
C   --   NNNPS - IN/OUT - the number of nodes for each set
C   --   IXNNPS - IN/OUT - the index of the first node for each set
C   --   LTNNPS - IN/OUT - the nodes for all sets
C   --   FACNPS - IN/OUT - the distribution factors for all sets
C   --   LTNX - SCRATCH - sized to hold the set nodes
C   --   FACX - SCRATCH - sized to hold the set factors
C   --   IXNPS - SCRATCH - size = NUMNPS
C   --   NNX - SCRATCH - size = NUMNPS
C   --   ISCR - SCRATCH - size = NUMNPS
C   --   NUMNP - IN -- number of nodes in model
C   --   NODSCR - SCRATCH - size = NUMNOD

      include 'gj_namlen.blk'

      INTEGER ISTAT(*)
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER IXNNPS(*)
      INTEGER LTNNPS(*), LTNX(*)
      REAL FACNPS(*), FACX(*)
      INTEGER IXNPS(*)
      INTEGER NNX(*)
      INTEGER ISCR(*)
      INTEGER NODSCR(*)
      character*(namlen) namns(*)
      character*(namlen) namsc(*)
      
      IF (NUMNPS .LE. 0) RETURN

      JNPS = 0
      JNN0 = 0
      DO 120 INPS = 1, NUMNPS

         IF (ISTAT(INPS) .EQ. 0) THEN
            NINSET = 1
            ISCR(NINSET) = INPS
         ELSE IF (ISTAT(INPS) .EQ. INPS) THEN
            CALL GETALL (INPS, NUMNPS, ISTAT, NINSET, ISCR)
         ELSE
            NINSET = 0
         END IF

         IF (NINSET .GT. 0) THEN
            JNPS = JNPS + 1
            IXNPS(JNPS) = INPS
            NNX(JNPS) = 0
            JOLD = JNN0+1
         END IF

         DO 110 ISET = 1, NINSET
            N = ISCR(ISET)
            JLEN = JNN0+1 - JOLD
            INN0 = IXNNPS(N) - 1

            call iniint(numnp, 0, nodscr)
            do 90 i=jold, jold+jlen-1
              nodscr(ltnx(i)) = 1
 90         continue

            NNEW = 0
            DO 100 I = 1, NNNPS(N)
 2             IF (nodscr(ltnnps(inn0+i)) .eq. 0) then
                  NNEW = NNEW + 1
                  LTNX(JNN0+NNEW) = LTNNPS(INN0+I)
                  FACX(JNN0+NNEW) = FACNPS(INN0+I)
               END IF
  100       CONTINUE
            NNX(JNPS) = NNX(JNPS) + NNEW
            JNN0 = JNN0 + NNEW
  110    CONTINUE
  120 CONTINUE

      CALL ORDIX (JNPS, IXNPS, NUMNPS, IDNPS, ISCR, IDNPS)
      CALL ORDNAM (JNPS, IXNPS, NUMNPS, NAMNS, NAMSC, NAMNS)
      CALL MOVINT (JNPS, NNX, NNNPS)
      NUMNPS = JNPS
      JNN = 1
      DO 130 INPS = 1, NUMNPS
         IXNNPS(INPS) = JNN
         JNN = JNN + NNNPS(INPS)
  130 CONTINUE
      LNPSNL = JNN - 1

      CALL MOVINT (LNPSNL, LTNX, LTNNPS)
      CALL MOVREA (LNPSNL, FACX, FACNPS)

      RETURN
      END
