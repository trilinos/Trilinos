C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SCAELE (A, IVAR,
     &   LENE, ISEVOK, NALIVR, ALIVE, VAR, WHOTIM, XE, YE, ZE,
     &   VALMN, NUMMN, XYZMN, ISTMN, VALMX, NUMMX, XYZMX, ISTMX)
C=======================================================================

C   --*** SCAELE *** (BLOT) Scale element variable by element block
C   --   Written by Amy Gilkey - revised 04/01/88
C   --
C   --SCAELE reads the values for the element variable from the database
C   --and finds the minimum and maximum value for each element block.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   IVAR - IN - the variable index (for GETVAR)
C   --   LENE - IN - the cumulative element counts by element block
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   NALIVR - IN - the element birth/death variable
C   --   ALIVE - SCRATCH - the element selection scratch array
C   --   VAR - SCRATCH - the variable array
C   --   WHOTIM - IN - true iff time step is a whole (versus history) time step
C   --   XE, YE, ZE - IN - the element coordinates
C   --   VALMN, VALMX - IN/OUT - the minimum and maximum value for each
C   --      element block
C   --   NUMMN, NUMMX - IN/OUT - the element number of the minimum and
C   --      maximum value for each element block
C   --   XYZMN, XYZMX - IN/OUT - the coordinates of NUMMN, NUMMX
C   --   ISTMN, ISTMX - IN/OUT - the step number of the minimum and maximum
C   --      value for each element block
C   --
C   --Common Variables:
C   --   Uses NUMEL, NELBLK, NVAREL, NSTEPS of /DBNUMS/

      include 'dbnums.blk'
      include 'mshopt.blk'

      DIMENSION A(*)
      INTEGER LENE(0:*)
      LOGICAL ISEVOK(NELBLK,NVAREL)
      LOGICAL ALIVE(*)
      REAL VAR(*)
      LOGICAL WHOTIM(*)
      REAL XE(*), YE(*), ZE(*)
      REAL VALMN(0:NELBLK), VALMX(0:NELBLK)
      INTEGER NUMMN(0:NELBLK), NUMMX(0:NELBLK)
      REAL XYZMN(3,0:NELBLK), XYZMX(3,0:NELBLK)
      INTEGER ISTMN(0:NELBLK), ISTMX(0:NELBLK)

      LOGICAL INIT
      CHARACTER CDUM

C   --Transfer birth/death variable to random file (for efficiency)

      IF (NALIVR .GT. 0) THEN
         MXSTEP = 0
         DO 100 ISTEP = 1, NSTEPS
            IF (WHOTIM(ISTEP)) MXSTEP = ISTEP
  100    CONTINUE
         CALL GETVAR (A, NALIVR, -1, -MXSTEP, NUMEL, VAR)
      END IF

      CALL DBVTYP_BL (IVAR, CDUM, IXVAR)

C   --Initialize the birth/death settings

      IF (NALIVR .LE. 0) THEN
         CALL INILOG (NUMEL, .TRUE., ALIVE)
      END IF

      DO 130 ISTEP = 1, NSTEPS
         IF (.NOT. WHOTIM(ISTEP)) GOTO 130

C      --Read the new birth/death settings (if any) and the variable values
C      --for this time step

         IF (NALIVR .GT. 0) THEN
C         --Order to optimize read from a sequential file
            IF (NALIVR .LT. IVAR) THEN
               CALL GETALV (A, NALIVR, ALIVAL, ISTEP, LENE, ISEVOK,
     &            ALIVE, ALIVE)
               CALL GETVAR (A, IVAR, -1, ISTEP, NUMEL, VAR)
            ELSE
               CALL GETVAR (A, IVAR, -1, ISTEP, NUMEL, VAR)
               CALL GETALV (A, NALIVR, ALIVAL, ISTEP, LENE, ISEVOK,
     &            ALIVE, ALIVE)
            END IF
         ELSE
            CALL GETVAR (A, IVAR, -1, ISTEP, NUMEL, VAR)
         END IF

C      --Find minimum and maximum variable values for element,
C      --by element block

         DO 120 IELB = 1, NELBLK
            IF (ISEVOK(IELB,IXVAR)) THEN

               INIT = (ISTEP .EQ. 1) .OR. (ISTMN(IELB) .LE. 0)
               IF (.NOT. INIT) THEN
C               --Reinitialize for this element block
                  VALMIN = VALMN(IELB)
                  NUMMIN = -999
                  VALMAX = VALMX(IELB)
                  NUMMAX = -999
               END IF

               DO 110 IEL = LENE(IELB-1)+1, LENE(IELB)
                  IF (ALIVE(IEL)) THEN
                     IF (INIT) THEN
                        VALMIN = VAR(IEL)
                        NUMMIN = IEL
                        VALMAX = VAR(IEL)
                        NUMMAX = IEL
                        INIT = .FALSE.
                     ELSE IF (VALMIN .GT. VAR(IEL)) THEN
                        VALMIN = VAR(IEL)
                        NUMMIN = IEL
                     ELSE IF (VALMAX .LT. VAR(IEL)) THEN
                        VALMAX = VAR(IEL)
                        NUMMAX = IEL
                     END IF
                  END IF
  110          CONTINUE

               IF (INIT) THEN
C               --Initialize for the case where no element has the
C               --element block
                  VALMN(IELB) = 0.0
                  NUMMN(IELB) = 0
                  ISTMN(IELB) = 0
                  VALMX(IELB) = 0.0
                  NUMMX(IELB) = 0
                  ISTMX(IELB) = 0
               ELSE
                  IF (NUMMIN .GT. 0) THEN
                     VALMN(IELB) = VALMIN
                     NUMMN(IELB) = NUMMIN
                     ISTMN(IELB) = ISTEP
                  END IF
                  IF (NUMMAX .GT. 0) THEN
                     VALMX(IELB) = VALMAX
                     NUMMX(IELB) = NUMMAX
                     ISTMX(IELB) = ISTEP
                  END IF
               END IF
            END IF
  120    CONTINUE
  130 CONTINUE

C   --Fill in coordinates for all element blocks

      DO 160 IELB = 1, NELBLK
         DO 140 I = 1, 3
            XYZMN(I,IELB) = 0.0
  140    CONTINUE
         DO 150 I = 1, 3
            XYZMX(I,IELB) = 0.0
  150    CONTINUE
         IF (ISTMN(IELB) .GT. 0) THEN
            IF (NDIM .GE. 1) XYZMN(1,IELB) = XE(NUMMN(IELB))
            IF (NDIM .GE. 2) XYZMN(2,IELB) = YE(NUMMN(IELB))
            IF (NDIM .GE. 3) XYZMN(3,IELB) = ZE(NUMMN(IELB))
            IF (NDIM .GE. 1) XYZMX(1,IELB) = XE(NUMMX(IELB))
            IF (NDIM .GE. 2) XYZMX(2,IELB) = YE(NUMMX(IELB))
            IF (NDIM .GE. 3) XYZMX(3,IELB) = ZE(NUMMX(IELB))
         END IF
  160 CONTINUE

C   --Calculate min/max for element by element block

      INIT = .TRUE.
      DO 170 IELB = 1, NELBLK
         IF (ISTMN(IELB) .GT. 0) THEN
            IF (INIT) THEN
               IMIN = IELB
               IMAX = IELB
               INIT = .FALSE.
            ELSE
               IF (VALMN(IMIN) .GT. VALMN(IELB)) THEN
                  IMIN = IELB
               END IF
               IF (VALMX(IMAX) .LT. VALMX(IELB)) THEN
                  IMAX = IELB
               END IF
            END IF
         END IF
  170 CONTINUE

C   --Store minimum and maximum variable values for all element blocks

      IF (INIT) THEN
         VALMN(0) = 0.0
         NUMMN(0) = 0
         DO 180 I = 1, 3
            XYZMN(I,0) = 0.0
  180    CONTINUE
         ISTMN(0) = 0
         VALMX(0) = 0.0
         NUMMX(0) = 0
         DO 190 I = 1, 3
            XYZMX(I,0) = 0.0
  190    CONTINUE
         ISTMX(0) = 0
      ELSE
         VALMN(0) = VALMN(IMIN)
         NUMMN(0) = NUMMN(IMIN)
         DO 200 I = 1, 3
            XYZMN(I,0) = XYZMN(I,IMIN)
  200    CONTINUE
         ISTMN(0) = ISTMN(IMIN)
         VALMX(0) = VALMX(IMAX)
         NUMMX(0) = NUMMX(IMAX)
         DO 210 I = 1, 3
            XYZMX(I,0) = XYZMX(I,IMAX)
  210    CONTINUE
         ISTMX(0) = ISTMX(IMAX)
      END IF

      RETURN
      END
