C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CMDVAR (MODDET, MODTYP, NAMES,
     &   INLINE, IFLD, INTYP, CFIELD,
     &   IDTVAR, *)
C=======================================================================

C   --*** CMDVAR *** (DETOUR) Read the command variables
C   --   Written by Amy Gilkey - revised 03/03/88
C   --
C   --CMDVAR reads the variables associated with the command and checks
C   --that they are appropriate for the command.
C   --
C   --Parameters:
C   --   MODDET - IN - the display mode for the command (as in /DETOPT/)
C   --   MODTYP - IN/OUT - the mode type for the command (set for vector mode)
C   --      (as in /DETOPT/)
C   --   NAMES - IN - the variable names
C   --   INLINE - IN/OUT - the parsed input line for the log file
C   --   IFLD - IN/OUT - the free-field current field number
C   --   INTYP - IN - the free-field type
C   --   CFIELD - IN - the free-field character field
C   --   IDTVAR - IN/OUT - the current variables
C   --   * - return statement if an error is encountered
C   --
C   --Common Variables:
C   --   Uses NDIM, NVARNP, NVAREL, NSTEPW of /DBNUMS/

      include 'params.blk'
      include 'dbnums.blk'

      INTEGER GETVNM

      CHARACTER*(*) INLINE
      CHARACTER*(*) MODDET, MODTYP
      CHARACTER*(*) NAMES(*)
      INTEGER INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER IDTVAR(4)

      LOGICAL FFEXST

      CHARACTER*(MXNAME) WORD
      LOGICAL GIVEN1
      INTEGER NVNEW(4)

      CHARACTER TYP1, TYP2

      LOGICAL SHRTHD

      NEEDED = 0
      IF (MODDET .EQ. 'CONTOUR') THEN
         NEEDED = 1
      ELSE IF (MODDET .EQ. 'ELEMCONT') THEN
         NEEDED = 1
      ELSE IF (MODDET .EQ. 'VECTOR') THEN
         IF ((MODTYP .NE. 'SIGMAX') .AND. (MODTYP .NE. 'SIGMIN')) THEN
            NEEDED = NDIM
         ELSE
            NEEDED = 3
         END IF
      ELSE IF (MODDET .EQ. 'SYMBOL') THEN
         NEEDED = 1
      ELSE IF (MODDET .EQ. 'GAUSS') THEN
         NEEDED = 4
      END IF

      IF (NEEDED .GT. 0) THEN
         IF (NSTEPW .LE. 0) THEN
            CALL PRTERR ('CMDERR',
     &         'Display mode is invalid with no whole time steps')
            GOTO 110
         END IF

         NMIN = 9999
         NMAX = 0
         GIVEN1 = (FFEXST (IFLD, INTYP))
         SHRTHD = .FALSE.
         DO 100 IV = 1, NEEDED
            IF (.NOT. SHRTHD) THEN
               IF (GIVEN1) THEN
                  IF (.NOT. FFEXST (IFLD, INTYP)) THEN
                     CALL PRTERR ('CMDERR',
     &                  'All variables must be given')
                     GOTO 110
                  END IF
               END IF
               CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
               IF (WORD .NE. ' ') CALL FFADDC (WORD, INLINE)
            ENDIF

            IF ((MODDET .EQ. 'VECTOR') .AND. (WORD .EQ. '0')) THEN
                  IVAR = 0
            ELSE IF (WORD .NE. ' ') THEN
               IVAR = LOCSTR (WORD, NVARHI+NVARGL+NVARNP+NVAREL, NAMES)
               IF ((MODDET .EQ. 'VECTOR') .AND.(IVAR .LE. 0)) THEN
                  IF (IV .EQ. 1) THEN
                    ivar = getvnm(word(:lenstr(word)), 'X',
     &                      NVARHI+NVARGL+NVARNP+NVAREL, NAMES)
                    IF (IVAR .GT. 0) SHRTHD = .TRUE.
                  ELSE IF (IV .EQ. 2) THEN
                    ivar = getvnm(word(:lenstr(word)), 'Y',
     &                      NVARHI+NVARGL+NVARNP+NVAREL, NAMES)
                  ELSE IF (IV .EQ. 3) THEN
                    ivar = getvnm(word(:lenstr(word)), 'Z',
     &                      NVARHI+NVARGL+NVARNP+NVAREL, NAMES)
                  END IF
               END IF
               IF (IVAR .LE. 0) THEN
                  CALL PRTERR ('CMDERR',
     &               'Variable "' // WORD(:LENSTR(WORD))
     &               // '" does not exist')
                  GOTO 110
               END IF
            ELSE
               IVAR = IDTVAR(IV)
               IF (IVAR .LE. 0) THEN
                  CALL PRTERR ('CMDERR',
     &               'Variable must be given')
                  GOTO 110
               END IF
            END IF
            NVNEW(IV) = IVAR

            IF (IVAR .NE. 0) THEN
               NMIN = MIN (NMIN, IVAR)
               NMAX = MAX (NMAX, IVAR)
            END IF
  100    CONTINUE
         SHRTHD = .FALSE.
         IF (NMIN .GT. NMAX) NMIN = NMAX

         CALL DBVTYP_BL (NMIN, TYP1, IDUM)
         CALL DBVTYP_BL (NMAX, TYP2, IDUM)

         IF (MODDET .EQ. 'CONTOUR') THEN
            IF ((TYP1 .NE. 'N') .AND. (TYP1 .NE. 'E')) THEN
               CALL PRTERR ('CMDERR',
     &            'Contour variable must be nodal or element variable')
               GOTO 110
            END IF
         ELSE IF (MODDET .EQ. 'ELEMCONT') THEN
            IF (TYP1 .NE. 'E') THEN
               CALL PRTERR ('CMDERR',
     &            'Contour variable must be an element variable')
               GOTO 110
            END IF
         ELSE IF (MODDET .EQ. 'VECTOR') THEN
            IF (MODTYP .EQ. ' ') THEN
               CALL DBVTYP_BL (NMIN, TYP1, IDUM)
               IF ((TYP1 .EQ. 'N') .OR. (TYP1 .EQ. ' ')) THEN
                  MODTYP = 'NODE'
               ELSE
                  MODTYP = 'ELEMENT'
               END IF
            END IF
            IF ((MODTYP .EQ. 'NODE') .OR. (MODTYP .EQ. 'ELEMENT')) THEN
               IF (((TYP1 .NE. 'N') .AND. (TYP1 .NE. 'E'))
     &            .OR. (TYP1 .NE. TYP2)) THEN
                  CALL PRTERR ('CMDERR',
     &               'Vector variables must'
     &               // ' all be nodal or all be element')
                  GOTO 110
               END IF
            ELSE IF ((MODTYP .EQ. 'SIGMAX')
     &         .OR. (MODTYP .EQ. 'SIGMIN')) THEN
               IF (TYP1 .NE. 'E') THEN
                  CALL PRTERR ('CMDERR',
     &               'Stress variables must be element variables')
                  GOTO 110
               END IF
            END IF
         ELSE IF (MODDET .EQ. 'SYMBOL') THEN
            IF (TYP1 .NE. 'E') THEN
               CALL PRTERR ('CMDERR',
     &            'Symbol variable must be an element variable')
               GOTO 110
            END IF
         ELSE IF (MODDET .EQ. 'GAUSS') THEN
            IF ((TYP1 .NE. 'E') .AND. (TYP2 .NE. 'E')) THEN
               CALL PRTERR ('CMDERR',
     &            'Gauss variables must be element variables')
               GOTO 110
            END IF
         END IF

         CALL CPYINT (NEEDED, NVNEW, IDTVAR)
      END IF

      RETURN

  110 CONTINUE
      RETURN 1
      END

      integer function getvnm(word, suffix, icnt, names)
      include 'params.blk'

      CHARACTER*(MXNAME) WORD
      character*1 suffix
      CHARACTER*(*) NAMES(*)
      CHARACTER*(MXNAME) TEMP

      temp = word(:lenstr(word))//suffix
      getvnm = locstr(temp, icnt, names)

      if (getvnm .le. 0) then
        temp = word(:lenstr(word))//'_'// suffix
        getvnm = locstr(temp, icnt, names)
      endif

      if (getvnm .le. 0) then
        temp = word(:lenstr(word))//'.'// suffix
        getvnm = locstr(temp, icnt, names)
      endif

      return
      end
