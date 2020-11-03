C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PRNPS (OPTION, NOUT, NUMNPS, LISNPS, LNPSNL,
     &   IDNPS, NNNPS, IXNNPS, LTNNPS, FACNPS, NSNAME, MAPND)
C=======================================================================

C   --*** PRNPS *** (BLOT) Display database node set
C   --   Written by Amy Gilkey - revised 01/18/88
C   --
C   --PRNPS displays the node sets.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --      ' ' - set summary (number of nodes, etc)
C   --      'N' - nodes in set
C   --      'F' - distribution factors for set
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NUMNPS - IN - the number of node sets
C   --   LISNPS - IN - the indices of the selected node sets
C   --   LNPSNL - IN - the number of nodes for all sets
C   --   IDNPS - IN - the node set ID for each set
C   --   NNNPS - IN - the number of nodes for each set
C   --   IXNNPS - IN - the index of the first node for each set
C   --   LTNNPS - IN - the nodes for all sets
C   --   FACNPS - IN - the distribution factors for all sets

      CHARACTER*(*) OPTION
      INTEGER LISNPS(0:*)
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER IXNNPS(*)
      INTEGER LTNNPS(*)
      REAL FACNPS(*)
      CHARACTER*(*) NSNAME(*)
      INTEGER MAPND(*)

      LOGICAL ISABRT
      LOGICAL DONOD, DOFAC
      CHARACTER*20 STRA, STRB

      DONOD = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0))
      DOFAC = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'F') .GT. 0))

      IF (NOUT .GT. 0) THEN
         IF (DONOD .AND. DOFAC) THEN
            WRITE (NOUT, 10020) 'Node List and Distribution Factors'
         ELSE IF (DONOD) THEN
            WRITE (NOUT, 10020) 'Node List'
         ELSE IF (DOFAC) THEN
            WRITE (NOUT, 10020) 'Distribution Factors'
         ELSE
            WRITE (NOUT, 10020)
         END IF
      END IF

      IF (NOUT .GT. 0) THEN
         WRITE (NOUT, *)
      ELSE
         WRITE (*, *)
      END IF

      WRITE (STRA, 10000, IOSTAT=IDUM) NUMNPS
10000  FORMAT ('(#', I4, ')')
      CALL PCKSTR (1, STRA)
      LSTRA = LENSTR (STRA)
      WRITE (STRB, 10010, IOSTAT=IDUM) LNPSNL
10010  FORMAT ('(index=', I9, ')')
      CALL PCKSTR (1, STRB)
      LSTRB = LENSTR (STRB)

      DO 100 IX = 1, LISNPS(0)
         IF (ISABRT ()) RETURN
         INPS = LISNPS(IX)
         WRITE (STRA, 10000, IOSTAT=IDUM) INPS
         CALL PCKSTR (1, STRA)
         WRITE (STRB, 10010, IOSTAT=IDUM) IXNNPS(INPS)
         CALL PCKSTR (1, STRB)
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10030, IOSTAT=IDUM)
     &         IDNPS(INPS), STRA(:LSTRA),
     &         NNNPS(INPS), STRB(:LSTRB),
     $         NSNAME(INPS)(:LENSTR(NSNAME(INPS)))
         ELSE
            WRITE (*, 10030, IOSTAT=IDUM)
     &         IDNPS(INPS), STRA(:LSTRA),
     &         NNNPS(INPS), STRB(:LSTRB),
     $         NSNAME(INPS)(:LENSTR(NSNAME(INPS)))
         END IF

         IF (DONOD .AND. (NNNPS(INPS) .GT. 0)) THEN
            IS = IXNNPS(INPS)
            IE = IS + NNNPS(INPS) - 1
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10040, IOSTAT=IDUM)
     &          (MAPND(LTNNPS(I)), I=IS,IE)
            ELSE
               WRITE (*, 10040, IOSTAT=IDUM)
     &          (MAPND(LTNNPS(I)), I=IS,IE)
            END IF
         END IF

         IF (DOFAC .AND. (NNNPS(INPS) .GT. 0)) THEN
            IS = IXNNPS(INPS)
            IE = IS + NNNPS(INPS) - 1
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10050, IOSTAT=IDUM)
     &            (FACNPS(I), I=IS,IE)
            ELSE
               WRITE (*, 10050, IOSTAT=IDUM)
     &            (FACNPS(I), I=IS,IE)
            END IF
         END IF
  100 CONTINUE

      RETURN

10020  FORMAT (/, 1X, 'Node Sets (Global Node Ids)', :, ' - ', A)
10030  FORMAT (1X, 'Set', I9, 1X, A, ':',
     &   I9, ' nodes', 1X, A,' name = "',A,'"')
10040  FORMAT ((1X, 8I11))
10050  FORMAT ((1X, 6 (2X, E11.4)))
      END
