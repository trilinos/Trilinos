C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CHKVAR (MODDET, MODTYP, IVIEW, IDTVAR, NVOLD,
     &   NNDVAR, NEDVAR, OK)
C=======================================================================

C   --*** CHKVAR *** (DETOUR) Check the variables by their mode
C   --   Written by Amy Gilkey - revised 03/03/88
C   --
C   --CHKVAR counts the number of nodal and element database variables
C   --needed for the display modes.  It then checks that no variables have
C   --been changed for other modes.  The modes are then checked for
C   --consistency.  Mode types which are dependent on the variable type
C   --may be changed.
C   --
C   --Parameters:
C   --   MODDET - IN - the modes for all views (as in /DETOPT/)
C   --   MODTYP - IN/OUT - the mode types for all views (as in /DETOPT/)
C   --   IVIEW - IN - the view that has just been changed, <0 for final check
C   --   IDTVAR - IN - the current variables
C   --   NVOLD - IN - the old variables
C   --   NNDVAR - OUT - the number of nodal variables needed
C   --   NEDVAR - OUT - the number of element variables needed
C   --   OK - OUT - true iff the modes and variables are consistent
C   --
C   --Common Variables:
C   --   Uses NDIM, NVARNP, NVAREL of /DBNUMS/

      include 'dbnums.blk'

      CHARACTER*(*) MODDET(4), MODTYP(4)
      INTEGER IDTVAR(4), NVOLD(4)
      LOGICAL OK

      INTEGER NDEFVW, IXVW
      CHARACTER TYP1, TYP2

      NNDVAR = 0
      NEDVAR = 0
      NOVAR = 0
      DO 110 IVW = 1, NDEFVW (.FALSE.)
         I = IXVW (.FALSE., IVW)
         IF (MODDET(I) .EQ. 'CONTOUR') THEN
            N = 1
            NNDVAR = MAX (NNDVAR, N)
            CALL DBVTYP_BL (IDTVAR(1), TYP1, IDUM)
            IF (TYP1 .EQ. 'E') NEDVAR = MAX (NEDVAR, N)
         ELSE IF (MODDET(I) .EQ. 'ELEMCONT') THEN
            N = 1
            NEDVAR = MAX (NEDVAR, N)
         ELSE IF (MODDET(I) .EQ. 'VECTOR') THEN
            NMIN = 9999
            NMAX = 0
            DO 100 IV = 1, NDIM
               IF (IDTVAR(IV) .NE. 0) THEN
                  NMIN = MIN (NMIN, IDTVAR(IV))
                  NMAX = MAX (NMAX, IDTVAR(IV))
               END IF
  100       CONTINUE
            IF (NMIN .GT. NMAX) NMIN = NMAX
            IF ((MODTYP(I) .NE. 'SIGMAX')
     &         .OR. (MODTYP(I) .NE. 'SIGMIN')) THEN
               CALL DBVTYP_BL (NMIN, TYP1, IDUM)
               IF ((TYP1 .EQ. 'N') .OR. (TYP1 .EQ. ' ')) THEN
                  MODTYP(I) = 'NODE'
               ELSE
                  MODTYP(I) = 'ELEMENT'
               END IF
            END IF
            IF (MODTYP(I) .EQ. 'NODE') THEN
               N = NDIM
               NNDVAR = MAX (NNDVAR, N)
            ELSE IF (MODTYP(I) .EQ. 'ELEMENT') THEN
               N = NDIM
               NEDVAR = MAX (NEDVAR, N)
            ELSE IF ((MODTYP(I) .EQ. 'SIGMAX')
     &         .OR. (MODTYP(I) .EQ. 'SIGMIN')) THEN
               N = 3
               NEDVAR = MAX (NEDVAR, N)
            END IF
         ELSE IF (MODDET(I) .EQ. 'SYMBOL') THEN
            N = 1
            NEDVAR = MAX (NEDVAR, N)
         ELSE IF (MODTYP(I) .EQ. 'GAUSS') THEN
            N = 4
            NEDVAR = MAX (NEDVAR, N)
         ELSE
            N = 0
         END IF
         IF (I .NE. IVIEW) NOVAR = MAX (NOVAR, N)
  110 CONTINUE

      IF (IVIEW .NE. 0) THEN
         DO 120 I = 1, NOVAR
            IF (IDTVAR(I) .NE. NVOLD(I)) THEN
               CALL PRTERR ('CMDWARN',
     &            'Variable changes for all views')
               GOTO 130
            END IF
  120    CONTINUE
  130    CONTINUE
      END IF

      OK = .TRUE.
      DO 150 IVW = 1, NDEFVW (.FALSE.)
         I = IXVW (.FALSE., IVW)
         IF (MODDET(I) .EQ. 'CONTOUR') THEN
            IF (IDTVAR(1) .EQ. 0) THEN
               CALL PRTERR ('CMDERR',
     &            'Contour variable is undefined')
               OK = .FALSE.
               GOTO 150
            END IF
         ELSE IF (MODDET(I) .EQ. 'ELEMCONT') THEN
            IF (IDTVAR(1) .EQ. 0) THEN
               CALL PRTERR ('CMDERR',
     &            'Contour variable is undefined')
               OK = .FALSE.
               GOTO 150
            END IF
            CALL DBVTYP_BL (IDTVAR(1), TYP1, IDUM)
            IF (TYP1 .NE. 'E') THEN
               CALL PRTERR ('CMDERR',
     &            'Contour variable must be an element variable')
               OK = .FALSE.
               GOTO 150
            END IF
         ELSE IF (MODDET(I) .EQ. 'VECTOR') THEN
            IF ((MODTYP(I) .EQ. 'NODE')
     &         .OR. (MODTYP(I) .EQ. 'ELEMENT')) THEN
               CALL DBVTYP_BL (NMIN, TYP1, IDUM)
               CALL DBVTYP_BL (NMAX, TYP2, IDUM)
               IF (TYP1 .NE. TYP2) THEN
                  CALL PRTERR ('CMDERR',
     &               'Vector variables must'
     &               // ' all be nodal or all be element')
                  OK = .FALSE.
                  GOTO 150
               END IF
            ELSE IF ((MODTYP(I) .EQ. 'SIGMAX')
     &         .OR. (MODTYP(I) .EQ. 'SIGMIN')) THEN
               CALL DBVTYP_BL (NMIN, TYP1, IDUM)
               IF (TYP1 .NE. 'E') THEN
                  CALL PRTERR ('CMDERR',
     &               'Stress variables must be element variables')
                  OK = .FALSE.
                  GOTO 150
               END IF
            END IF
         ELSE IF (MODDET(I) .EQ. 'SYMBOL') THEN
            IF (MODTYP(I) .NE. 'SPHERE'  .AND.
     &         MODTYP(I) .NE. 'FSPHERE') THEN
               IF (IDTVAR(1) .EQ. 0) THEN
                  CALL PRTERR ('CMDERR',
     &               'Symbol variable is undefined')
                  OK = .FALSE.
                  GOTO 150
               END IF
               CALL DBVTYP_BL (IDTVAR(1), TYP1, IDUM)
               IF (TYP1 .NE. 'E') THEN
                  CALL PRTERR ('CMDERR',
     &               'Symbol variable must be an element variable')
                  OK = .FALSE.
                  GOTO 150
               END IF
            ENDIF
         ELSE IF (MODTYP(I) .EQ. 'GAUSS') THEN
            DO 140 J = 1, 4
               IF (IDTVAR(J) .EQ. 0) THEN
                  CALL PRTERR ('CMDERR',
     &               'Gauss variable is undefined')
                  OK = .FALSE.
                  GOTO 150
               END IF
               CALL DBVTYP_BL (IDTVAR(J), TYP1, IDUM)
               IF (TYP1 .NE. 'E') THEN
                  CALL PRTERR ('CMDERR',
     &               'Gauss variables must be element variables')
                  OK = .FALSE.
                  GOTO 150
               END IF
  140       CONTINUE
         END IF
  150 CONTINUE

      RETURN
      END
