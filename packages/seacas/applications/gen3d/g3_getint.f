C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GETINT (TYPE, IFLD, INTYP, IFIELD, RFIELD,
     *  NBLK, NRTRAN, D3TRAN, ZGRAD, NEREPL, NNREPL, DIM3,
     *  MAXBLK, *)
C=======================================================================

C   --*** GETINT *** (GEN3D) Input number of intervals, translation/
C   --                       rotation distances and gradients
C   --
C   --   Written by Greg Sjaardema - 5/11/89
C   --
C   --Parameters:
C   --   TYPE    -  IN -  Type of transformation: 'rotation' 'translation'
C   --   IFLD    - I/O -  Index of the field to be read
C   --   INTYP   -  IN -  the free-field reader type array
C   --   IFIELD  -  IN -  the free-field reader integer array
C   --   RFIELD  -  IN -  the free-field reader real array
C   --   NBLK    - OUT -  number of transformation blocks
C   --   NRTRAN  - OUT -  number of transformation increments
C   --   D3TRAN  - OUT -  distance of transformation increments
C   --   ZGRAD   - OUT -  gradient of transformation increments
C   --   NEREPL  - OUT -  total number of 3D elements / 2D element
C   --   NNREPL  - OUT -  total number of 3D nodes / 2D nodes
C   --   DIM3    - OUT -  total transformation distance
C   --   MAXBLK  -  IN -  maximum number of increments to read
C   --   *  --- Return statement for quit

      CHARACTER*(*) TYPE
      CHARACTER*80  PRMPTA, PRMPTB
      INTEGER   IFLD, INTYP(*), IFIELD(*)
      REAL      RFIELD(*), D3TRAN(*), ZGRAD(*), DIM3
      INTEGER   NBLK, NRTRAN(*), NEREPL, NNREPL
      LOGICAL   FFEXST

      PRMPTA = 'Expected number of ' // TYPE // 's'
      LA = LENSTR(PRMPTA)
      PRMPTB = 'total ' // TYPE
      LB = LENSTR(PRMPTB)

      IF (TYPE(:1) .EQ. 'R' .OR. TYPE(:1) .EQ. 'r') THEN
        TRDEF = 360.0
      ELSE IF (TYPE(:1) .EQ. 'T' .OR. TYPE(:1) .EQ. 't') THEN
        TRDEF =   1.0
      ELSE
        TRDEF =   1.0
      END IF

      NBLK = 0
 10   CONTINUE
      IF ((NBLK .EQ. 0) .OR. FFEXST (IFLD, INTYP)) THEN

        CALL FFINTG (IFLD, INTYP, IFIELD,
     &    PRMPTA(10:LA), 1, NTRAN, *20)
        IF (NTRAN .LT. 1) THEN
          CALL PRTERR ('CMDERR', PRMPTA(:LA))
          GOTO 20
        END IF

        CALL FFREAL (IFLD, INTYP, RFIELD,
     &    PRMPTB(:LB), TRDEF, TRNAMT, *20)
        TRNAMT = ABS (TRNAMT)

C         IDEGR = NINT (DEGR/NROT)
C         IF (IDEGR .GE. 180) THEN
C            CALL PRTERR ('CMDERR',
C     &         'Single rotation cannot cover 180 degrees')
C            RETURN 1
C         END IF

        CALL FFREAL (IFLD, INTYP, RFIELD,
     &    'gradient', 1.0, GRADNT, *20)
        GRADNT = ABS (GRADNT)

        NBLK = NBLK + 1
        if (nblk .gt. maxblk) then
          CALL PRTERR('CMDERR',
     *      'Too many intervals specified. Reduce and rerun')
          stop 'Interval Error'
        else
          NRTRAN(NBLK) = NTRAN
          D3TRAN(NBLK) = TRNAMT
          IF (GRADNT .LE. 0.0) GRADNT = 1.0
          ZGRAD(NBLK) = GRADNT
        end if
        if (nblk .lt. maxblk) GOTO 10
      END IF

 20   CONTINUE
      IF (NBLK .LE. 0) RETURN 1

      NEREPL = 0
      DIM3 = 0
      DO 30 IBLK = 1, NBLK
        NEREPL = NEREPL + NRTRAN(IBLK)
        DIM3 = DIM3 + D3TRAN(IBLK)
 30   CONTINUE
      NNREPL = NEREPL + 1

      RETURN
      END
