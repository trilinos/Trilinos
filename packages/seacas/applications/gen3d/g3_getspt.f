C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE GETSPT(A)
      DIMENSION A(*)
      INCLUDE 'g3_splxyz.blk'

      PARAMETER (BINGO = 1.0E38)
      PARAMETER (MAXFLD = 10)
      PARAMETER (NDEFS = 64)
      CHARACTER*8 WORD, VERB
      INTEGER     INTYP(MAXFLD)
      CHARACTER*8 CFIELD(MAXFLD)
      INTEGER     IFIELD(MAXFLD)
      REAL        RFIELD(MAXFLD)

      LOGICAL MATSTR, HELP, ISHELP

      CHARACTER*8 CMDTBL(6)
      SAVE CMDTBL
C     --CMDTBL - the valid commands table

C     --Command table follows.  Remember to change the dimensioned size when
C     --changing the table.
      DATA CMDTBL /
     1     'SLOPE   ', 'END     ', 'EXIT    ', 'LIST    ', 'HELP    ',
     3     '        ' /

      CALL SHOCMD ('COMMANDS', CMDTBL)

C     ... Initialize default values

      IPT = 0
      SLBOT(1) = BINGO
      SLTOP(1) = BINGO
      SLBOT(2) = BINGO
      SLTOP(2) = BINGO
      ISPL = 1

C     ... Allocate arrays, guess on amount and increase if more points entered.

      CALL MDRSRV ('ZSPL',  KZSPL,  NDEFS)
      CALL MDRSRV ('XSPL',  KXSPL,  NDEFS)
      CALL MDRSRV ('YSPL',  KYSPL,  NDEFS)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) THEN
         CALL MEMERR
         STOP 'Memory Allocation in GETSPL'
      END IF

      NTOT = NDEFS

   10 CONTINUE

C     --Read command line

      WRITE (*, *)
      CALL FREFLD (0, 0, '   Spline Option > ', MAXFLD,
     &     IOSTAT, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)
      IF (IOSTAT .LT. 0) GOTO 20
      IF (NUMFLD .EQ. 0) GOTO 10
      INTYP(MIN(MAXFLD,NUMFLD)+1) = -999

      IFLD = 1
      CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
      CALL ABRSTR (VERB, WORD, CMDTBL)
      IF (VERB .EQ. ' ') VERB = WORD

      IF (VERB .EQ. 'EXIT' .OR. VERB .EQ. 'END') GO TO 20

      IF (VERB .EQ. 'HELP') THEN
         ISHELP = HELP ('GEN3D', 'COMMANDS', CFIELD(IFLD))
         IF (.NOT. ISHELP) CALL SHOCMD ('COMMANDS', CMDTBL)
         VERB = ' '
      END IF

      IF (INTYP(1) .EQ. 0) THEN
         IF ( VERB .EQ. 'SLOPE') THEN
            CALL FFCHAR (IFLD, INTYP, CFIELD, 'TOP', VERB)
            IF (MATSTR(VERB, 'TOP', 1) .OR.
     &           MATSTR(VERB, 'FRONT', 1)) THEN
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &              'X slope at bottom end', 0.0, SLTOP(1), *10)
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &              'Y slope at bottom end', 0.0, SLTOP(2), *10)
            END IF

         ELSE IF (VERB .EQ. 'LIST') THEN
            CALL SHOCMD ('COMMANDS', CMDTBL)
         END IF
      ELSE IF (NUMFLD .GE. 2) THEN
         IPT = IPT + 1
         IF (IPT .GT. NTOT) THEN
            NTOT = NTOT + NDEFS
            CALL MDLONG ('ZSPL', KZSPL, NTOT)
            CALL MDLONG ('XSPL', KXSPL, NTOT)
            CALL MDLONG ('YSPL', KYSPL, NTOT)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) THEN
               CALL MEMERR
               STOP 'Memory Allocation in GETSPL'
            END IF
         END IF
         A(KZSPL + IPT - 1) = RFIELD(1)
         A(KXSPL + IPT - 1) = RFIELD(2)
         A(KYSPL + IPT - 1) = RFIELD(3)
      END IF
      GO TO 10

C     ... Done reading data, compress arrays and allocate remaining arrays

   20 CONTINUE
      CALL MDLONG ('ZSPL',  KZSPL, IPT)
      CALL MDLONG ('XSPL',  KXSPL, IPT)
      CALL MDLONG ('YSPL',  KYSPL, IPT)

      CALL MDRSRV ('XSPL2', KXSPL2, IPT)
      CALL MDRSRV ('YSPL2', KYSPL2, IPT)
      CALL MDRSRV ('SCR',   KSCR,   IPT)

      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) THEN
         CALL MEMERR
         STOP 'Memory Allocation in GETSPL'
      END IF

      NSPLT = IPT

      RETURN
      END
