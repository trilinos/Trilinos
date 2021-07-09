C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE CKNPS (NUMNPS, LNPSNL, NUMNP,
     &   IDNPS, NNNPS, IXNNPS, LTNNPS, FACNPS, ICHECK)
C=======================================================================

C   --*** CKNPS *** (EXPLORE) Check database nodal point sets
C   --
C   --CKNPS checks the nodal point set information.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NUMNPS - IN - the number of nodal point sets
C   --   LNPSNL - IN - the number of nodes for all sets
C   --   NUMNP - IN - the number of nodes
C   --   IDNPS - IN - the nodal point set ID for each set
C   --   NNNPS - IN - the number of nodes for each set
C   --   IXNNPS - IN - the index of the first node for each set
C   --   LTNNPS - IN - the nodes for all sets
C   --   FACNPS - IN - the distribution factors for all sets
C   --   ICHECK - SCRATCH - size = MAX (NUMNPS, LNPSNL)

      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER IXNNPS(*)
      INTEGER LTNNPS(*)
      REAL FACNPS(*)
      INTEGER ICHECK(*)

      CHARACTER*5 STRA

C   --Check for unique identifier

      DO 100 INPS = 1, NUMNPS
         IF (LOCINT (IDNPS(INPS), INPS-1, IDNPS) .GT. 0) THEN
            CALL INTSTR (1, 0, IDNPS(INPS), STRA, LSTRA)
            CALL PRTERR ('CMDSPEC', 'Nodal point set ID '
     &         // STRA(:LSTRA) // ' is not unique')
         END IF
  100 CONTINUE

C   --Check number of nodes in nodal point sets

      NNPS = 0
      DO 110 INPS = 1, NUMNPS
         NNPS = MAX (NNPS, IXNNPS(INPS) + NNNPS(INPS) - 1)
  110 CONTINUE

      IF (NNPS .NE. LNPSNL) THEN
         CALL PRTERR ('WARNING', 'Maximum node index' //
     &      ' in all nodal point sets does not match total')
      END IF

C   --Check all nodes in node point sets are within node range

      CALL CHKRNG (LTNNPS, LNPSNL, NUMNP, NZERO, NERR)
      IF (NERR .GT. 0) THEN
        CALL PRTERR ('FATAL',
     &    'Nodal point set node ids are out of range')
      END IF
      IF (NZERO .GT. 0) THEN
        CALL PRTERR ('FATAL',
     &    'Nodal point set node ids are zero')
      END IF

      RETURN
      END
