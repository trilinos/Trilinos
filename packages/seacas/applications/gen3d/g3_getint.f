C Copyright(C) 2011-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C
C * Neither the name of NTESS nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
C
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
C
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
