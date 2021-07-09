C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PRESS (OPTION, NOUT, NUMESS, LISESS, LESSEL, LESSNL,
     &   IDESS, NEESS, NNESS, IXEESS, IXNESS, LTEESS, LTNESS, FACESS,
     $   SSNAME, MAPEL, MAPND)
C=======================================================================

C   --*** PRESS *** (BLOT) Display database side set
C   --   Written by Amy Gilkey - revised 01/18/88
C   --
C   --PRESS displays the side sets.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --      ' ' - set summary (number of nodes, etc)
C   --      'E' - elements in set (may not be combined with 'N' or 'F')
C   --      'N' - nodes in set (may not be combined with 'E')
C   --      'F' - distribution factors for set (may not be combined with 'E')
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NUMESS - IN - the number of side sets
C   --   LISESS - IN - the indices of the selected side sets
C   --   LESSEL - IN - the number of elements for all sets
C   --   LESSNL - IN - the number of nodes for all sets
C   --   IDESS - IN - the side set ID for each set
C   --   NEESS - IN - the number of elements for each set
C   --   NNESS - IN - the number of nodes for each set
C   --   IXEESS - IN - the index of the first element for each set
C   --   IXNESS - IN - the index of the first node for each set
C   --   LTEESS - IN - the elements for all sets
C   --   LTNESS - IN - the nodes for all sets
C   --   FACESS - IN - the distribution factors for all sets

      CHARACTER*(*) OPTION
      INTEGER LISESS(0:*)
      INTEGER IDESS(*)
      INTEGER NEESS(*)
      INTEGER NNESS(*)
      INTEGER IXEESS(*)
      INTEGER IXNESS(*)
      INTEGER LTEESS(*)
      INTEGER LTNESS(*)
      REAL FACESS(*)
      CHARACTER*(*) SSNAME(*)
      INTEGER MAPEL(*)
      INTEGER MAPND(*)

      LOGICAL ISABRT
      LOGICAL DOELE, DONOD, DOFAC
      CHARACTER*20 STRA, STRB, STRC

      DOELE = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'E') .GT. 0))
      DONOD = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0))
      DOFAC = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'F') .GT. 0))

      IF (NOUT .GT. 0) THEN
         IF (DOELE) THEN
            WRITE (NOUT, 10020) 'ELEMENT LIST'
         ELSE IF (DONOD .AND. DOFAC) THEN
            WRITE (NOUT, 10020) 'NODE LIST AND DISTRIBUTION FACTORS'
         ELSE IF (DONOD) THEN
            WRITE (NOUT, 10020) 'NODE LIST'
         ELSE IF (DOFAC) THEN
            WRITE (NOUT, 10020) 'DISTRIBUTION FACTORS'
         ELSE
            WRITE (NOUT, 10020)
         END IF
      END IF

      IF (NOUT .GT. 0) THEN
         WRITE (NOUT, *)
      ELSE
         WRITE (*, *)
      END IF

      WRITE (STRA, 10000, IOSTAT=IDUM) NUMESS
10000  FORMAT ('(#', I4, ')')
      CALL PCKSTR (1, STRA)
      LSTRA = LENSTR (STRA)
      WRITE (STRB, 10010, IOSTAT=IDUM) LESSEL
10010  FORMAT ('(index=', I9, ')')
      CALL PCKSTR (1, STRB)
      LSTRB = LENSTR (STRB)
      WRITE (STRC, 10010, IOSTAT=IDUM) LESSNL
      CALL PCKSTR (1, STRC)
      LSTRC = LENSTR (STRC)

      DO 100 IX = 1, LISESS(0)
         IF (ISABRT ()) RETURN
         IESS = LISESS(IX)

         WRITE (STRA, 10000, IOSTAT=IDUM) IESS
         CALL PCKSTR (1, STRA)
         WRITE (STRB, 10010, IOSTAT=IDUM) IXEESS(IESS)
         CALL PCKSTR (1, STRB)
         WRITE (STRC, 10010, IOSTAT=IDUM) IXNESS(IESS)
         CALL PCKSTR (1, STRC)
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10030, IOSTAT=IDUM)
     &         IDESS(IESS), STRA(:LSTRA),
     &         NEESS(IESS), STRB(:LSTRB), NNESS(IESS), STRC(:LSTRC),
     $         SSNAME(IESS)(:LENSTR(SSNAME(IESS)))
         ELSE
            WRITE (*, 10030, IOSTAT=IDUM)
     &         IDESS(IESS), STRA(:LSTRA),
     &         NEESS(IESS), STRB(:LSTRB), NNESS(IESS), STRC(:LSTRC),
     $         SSNAME(IESS)(:LENSTR(SSNAME(IESS)))
         END IF

         IF (DOELE .AND. (NEESS(IESS) .GT. 0)) THEN
            IS = IXEESS(IESS)
            IE = IS + NEESS(IESS) - 1
            IF (NOUT .GT. 0) THEN
              WRITE (NOUT, 10040, IOSTAT=IDUM)
     &          (MAPEL(LTEESS(I)), I=IS,IE)
            ELSE
              WRITE (*, 10040, IOSTAT=IDUM)
     &          (MAPEL(LTEESS(I)), I=IS,IE)
            END IF
          END IF

          IF (DONOD .AND. (NNESS(IESS) .GT. 0)) THEN
            IS = IXNESS(IESS)
            IE = IS + NNESS(IESS) - 1
            IF (NOUT .GT. 0) THEN
              WRITE (NOUT, 10040, IOSTAT=IDUM)
     &          (MAPND(LTNESS(I)), I=IS,IE)
            ELSE
              WRITE (*, 10040, IOSTAT=IDUM)
     &          (MAPND(LTNESS(I)), I=IS,IE)
            END IF
         END IF

         IF (DOFAC .AND. (NNESS(IESS) .GT. 0)) THEN
            IS = IXNESS(IESS)
            IE = IS + NNESS(IESS) - 1
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10050, IOSTAT=IDUM)
     &            (FACESS(I), I=IS,IE)
            ELSE
               WRITE (*, 10050, IOSTAT=IDUM)
     &            (FACESS(I), I=IS,IE)
            END IF
         END IF
  100 CONTINUE

      RETURN

10020  FORMAT (/, 1X, 'SIDE SETS', :, ' - ', A)
10030  FORMAT (1X, 'Set', I9, 1X, A, ':',
     &   I9, ' elements', 1X, A,
     &   I9, ' nodes', 1X, A,' name = "',A,'"')
10040  FORMAT ((1X, 8I11))
10050  FORMAT ((1X, 6 (2X, E11.4)))
      END
