C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C============================================================================
      SUBROUTINE LINTHK (INLINE, IFLD, INTYP, IFIELD, RFIELD, CFIELD,
     &   RESET)
C============================================================================

C   --*** LINTHK ***  (BLOT) Process LINETHICKNESS command
C   --   Written by John Glick - 1/13/88
C   --
C   --Parameters:
C   --   INLINE - IN/OUT - the parsed input lines for the log file
C   --   IFLD, INTYP, IFIELD, CFIELD, - IN/OUT - the free-field reader
C   --          index and character field.
C   --   RESET - IN - =.TRUE. if call is only to reset the linethickness
C   --                        parameter.
C   --                 .FALSE. if call is to set the parameters.
      CHARACTER*(*) INLINE(*)
      INTEGER IFLD, INTYP(*), IFIELD(*)
      REAL RFIELD(*)
      CHARACTER*(*) CFIELD(*)
      LOGICAL RESET

      LOGICAL FFEXST

      INTEGER NUMSP
      CHARACTER*8 LINIDA, LINIDF
      REAL TKSPCR
      CHARACTER*8 TKSPCC, TKSPCF

      COMMON /LINTHC/ MSHBND, BLKBND, ELEBND, THKNSS
      REAL MSHBND, BLKBND, ELEBND
C      --      Line thickness specification for lines appearing
C      --      on mesh plots.  Specification is a real value in the
C      --      range 0. - 1000., with 0. being the thinest line and
C      --      1000. being the thickest.
C      -- MSHBND - Thickness of lines forming the mesh boundary.
C      -- BLKBND - Thickness of lines forming element block boundaries.
C      -- ELEBND - Thickness of lines forming element boundaries.
      REAL THKNSS(3)
C     --       Line thickness specifications for THICK, MEDIUM,
C     --       and THIN keywords.

      REAL THKRES(3)
      CHARACTER*8 LINLST(5)
      CHARACTER*8 THKNAM(4)
      DATA LINLST /'MESH    ', 'BLOCK   ', 'ELEMENT ',
     &   'ALL     ', '        '/
      DATA THKNAM /'THIN    ', 'MEDIUM  ', 'THICK   ',
     &   '        '/
      DATA THKRES /280., 200., 160./

C *****************************************************************

      IF (RESET) THEN
         MSHBND = THKRES(1)
         BLKBND = THKRES(2)
         ELEBND = THKRES(3)
         THKNSS(1) = 280.
         THKNSS(2) = 200.
         THKNSS(3) = 160.
      ELSE
         NUMSP = 0

C           Check for existence of a field.
  100    CONTINUE
         IF (FFEXST (IFLD, INTYP)) THEN

C           Check that next field has characters in it.
            IF (INTYP(IFLD) .GE. 0) THEN
C              Check for 'MESH', 'BLOCK', 'ELEMENT', or, 'ALL'
C              field.
               LINIDA = CFIELD(IFLD)
               IFLD = IFLD + 1
               CALL ABRSTR (LINIDF, LINIDA, LINLST)
               IF (LINIDF .EQ. ' ') THEN
                  WRITE (*, 10000) LINIDA
10000              FORMAT (1X, A, ' not a valid line type identifier.'/,
     &               1x, 'Expected ''MESH'', ''BLOCK'', ''ELEMENT'',',
     &               ' or ''ALL'' '/1x,
     &               'Rest of command not processed.')
                  GO TO 120
               ELSE

C                 Check for existence of another field
                  IF (FFEXST (IFLD, INTYP)) THEN
                     IF (INTYP(IFLD) .GE. 1) THEN
C                       Real value specified for line thickness
                        TKSPCR = RFIELD(IFLD)
                        IFLD = IFLD + 1
                        IF ((TKSPCR .LT. 0.0)  .OR.
     &                     (TKSPCR .GT. 1000.)) THEN
                           WRITE (*, 10010) TKSPCR
10010                       FORMAT (1X, E12.5, ' not a valid line ',
     &                        'thickness specification.  It must',
     &                        ' be ''THIN'','/1x, '''MEDIUM'',',
     &                        ' ''THICK'', or a value between',
     &                        ' 0.0 and 1000.')
                           GO TO 110
                        ELSE
                           IF ((LINIDF .EQ. 'MESH')  .OR.
     &                        (LINIDF .EQ. 'ALL')) MSHBND = TKSPCR
                           IF ((LINIDF .EQ. 'BLOCK')  .OR.
     &                        (LINIDF .EQ. 'ALL')) BLKBND = TKSPCR
                           IF ((LINIDF .EQ. 'ELEMENT')  .OR.
     &                        (LINIDF .EQ. 'ALL')) ELEBND = TKSPCR
                           NUMSP = NUMSP + 1
                           CALL FFADDC (LINIDF, INLINE(1))
                           CALL FFADDR (TKSPCR, INLINE(1))
                        ENDIF
                     ELSE
C                       Character string specified for line thickness
                        TKSPCC = CFIELD(IFLD)
                        IFLD = IFLD + 1
                        CALL ABRSTR (TKSPCF, TKSPCC, THKNAM)
                        IF (TKSPCF .EQ. ' ') THEN
                           WRITE (*, 10020) TKSPCC
10020                       FORMAT (1X, A, ' not a valid line ',
     &                        'thickness specification.  It must be',
     &                        ' ''THIN'','/1x, '''MEDIUM'',',
     &                        ' ''THICK'', or a value between',
     &                        ' 0.0 and 1000.')
                           GO TO 110
                        ELSE
                           IF (TKSPCF .EQ. 'THIN') THEN
                              TKSPCR = THKNSS(3)
                           ELSE IF (TKSPCF .EQ. 'MEDIUM') THEN
                              TKSPCR = THKNSS(2)
                           ELSE IF (TKSPCF .EQ. 'THICK') THEN
                              TKSPCR = THKNSS(1)
                           ENDIF
                           IF ((LINIDF .EQ. 'MESH')  .OR.
     &                        (LINIDF .EQ. 'ALL')) MSHBND = TKSPCR
                           IF ((LINIDF .EQ. 'BLOCK')  .OR.
     &                        (LINIDF .EQ. 'ALL')) BLKBND = TKSPCR
                           IF ((LINIDF .EQ. 'ELEMENT')  .OR.
     &                        (LINIDF .EQ. 'ALL')) ELEBND = TKSPCR
                           NUMSP = NUMSP + 1
                           CALL FFADDC (LINIDF, INLINE(1))
                           CALL FFADDC (TKSPCF, INLINE(1))
                        ENDIF
                     ENDIF

                  ELSE
                     WRITE (*, 10030) LINIDF
10030                 FORMAT (1X, ' No line thickness specification',
     &                  ' to go with ', A, ' keyword')
                     GO TO 120
                  ENDIF
               ENDIF
            ELSE
               WRITE (*, 10040)
10040           FORMAT(1X,' Valid line type identifier not specified.'/,
     &            1x, 'Expected ''MESH'', ''BLOCK'', ''ELEMENT'',',
     &            ' or ''ALL'' '/1x,
     &            'Rest of command not processed.')
               GO TO 120
            ENDIF

  110       CONTINUE
            GO TO 100
         ENDIF

  120    CONTINUE
         IF (NUMSP .EQ. 0) INLINE(1) = ' '
      ENDIF
      RETURN
      END
