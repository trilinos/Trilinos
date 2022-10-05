C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE VECTOR (ISYTYP, VAR, NUM, LENF, NLNKF, HIDENE,
     &   XNE, YNE, ZNE, IN2ELB, ISVOK, VECMAX, BLKCOL, IDELB, *)
C=======================================================================

C   --*** VECTOR *** (DETOUR) Plot nodal or element vector
C   --   Written by Amy Gilkey - revised 03/10/88
C   --   D. P. Flanagan, 12/08/83
C   --
C   --VECTOR draws a vector representing 2 or 3 variables for each node
C   --or element.  It processes each node or element by element block.
C   --Only nodes or elements of selected element blocks are drawn.
C   --
C   --Parameters:
C   --   ISYTYP - IN - the vector type (as in MODTYP of /DETOPT/)
C   --   VAR - IN - the vector variable values
C   --   NUM - IN - the number of nodes (if nodal) or faces (if element)
C   --   LENF - IN - the cumulative face counts by element block
C   --      (only if element)
C   --   NLNKF - IN - the number of nodes per face (only if element)
C   --   HIDENE(i) - IN - true iff node i or face i is hidden (3D only)
C   --   XNE, YNE, ZNE - IN - the nodal coordinates (if nodal)
C   --      or the face center coordinates (if element)
C   --   IN2ELB - IN - the element block for each node; (only if nodal)
C   --      <0 if not in any selected element block
C   --      =0 if in more than one selected element block
C   --   ISVOK - IN - ISVOK(i) is true iff the vector variables are defined
C   --      for element block i (only if element)
C   --   VECMAX - IN - the maximum vector variable value, scaled
C   --   BLKCOL - IN - the user selected colors of the element blocks.
C   --                    BLKCOL(0) = 1 if the user defined material
C   --                                colors should be used in mesh plots.
C   --                              = -1 if program selected colors should
C   --                                be used.
C   --                    BLKCOL(i) = the user selected color of element
C   --                               block i:
C   --                                  -2 - no color selected by user.
C   --                                  -1 - black
C   --                                   0 - white
C   --                                   1 - red
C   --                                   2 - green
C   --                                   3 - yellow
C   --                                   4 - blue
C   --                                   5 - cyan
C   --                                   6 - magenta
C   --   * - return statement if cancel function active
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses IS3DIM, NUMNPF of /D3NUMS/
C   --   Uses VECSCL of /ETCOPT/
C   --   Uses ROTMAT of /ROTOPT/
C   --   Uses DTW, VWSCL of /DEVDAT/

      include 'dbnums.blk'
      include 'd3nums.blk'
      include 'etcopt.blk'
      include 'rotopt.blk'
      include 'devdat.blk'

      CHARACTER*(*) ISYTYP
      REAL VAR(NUM,*)
      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      LOGICAL HIDENE(NUM)
      REAL XNE(NUM), YNE(NUM), ZNE(NUM)
      INTEGER IN2ELB(NUMNPF)
      LOGICAL ISVOK(NELBLK)
      REAL VECMAX
      INTEGER BLKCOL(0:NELBLK)
      INTEGER IDELB(*)

      LOGICAL GRABRT
      REAL ZERO3(3)
      SAVE ZERO3

      DATA ZERO3 / 0.0, 0.0, 0.0 /

      IF (VECSCL .EQ. 0.0) THEN
         VSCL = 1.0
      ELSE
         VSCL = VECSCL * 0.05 * VECMAX
      END IF
      ASCL = VWSCL

      IF (ISYTYP .EQ. 'NODE') THEN

         DO 110 IELB = 0, NELBLK

C         --Set the vector color
c            IF (IELB .GT. 0) THEN
c               ITEMP = IDELB(IELB)
c            ELSE
               ITEMP = IELB
c            ENDIF
            CALL UGRCOL (ITEMP, BLKCOL)

            DO 100 INP = 1, NUMNPF

               IF (IS3DIM) THEN
                  IF (HIDENE(INP)) GOTO 100
               END IF

               IF (IN2ELB(INP) .EQ. IELB) THEN
                  IF (GRABRT ()) RETURN 1

C               --Call vector routine

                  IF (.NOT. IS3DIM) THEN
                     CALL VEC (IS3DIM, XNE(INP), YNE(INP), 0.0,
     &                  VAR(INP,1), VAR(INP,2), 0.0, VSCL, ASCL)
                  ELSE
                     XVAR = VAR(INP,1)
                     YVAR = VAR(INP,2)
                     ZVAR = VAR(INP,3)
                     CALL BL_ROTATE (1, 1, ROTMAT, ZERO3,
     &                  XVAR, YVAR, ZVAR, XVAR, YVAR, ZVAR)
                     CALL VEC (IS3DIM, XNE(INP), YNE(INP), ZNE(INP),
     &                  XVAR, YVAR, ZVAR, VSCL, ASCL)
                  END IF
               END IF

  100       CONTINUE

            CALL PLTFLU
  110    CONTINUE

      ELSE IF (ISYTYP .EQ. 'ELEMENT') THEN

         DO 130 IELB = 1, NELBLK
            IF (ISVOK(IELB)) THEN

C            --Set the vector color
c               CALL UGRCOL (IDELB(IELB), BLKCOL)
               CALL UGRCOL (IELB, BLKCOL)

               DO 120 IFAC = LENF(IELB-1)+1, LENF(IELB)
                  IF (IS3DIM) THEN
                     IF (HIDENE(IFAC)) GOTO 120
                  END IF

C               --Call vector routine

                  IF (GRABRT ()) RETURN 1

                  IF (.NOT. IS3DIM) THEN
                     CALL VEC (IS3DIM, XNE(IFAC), YNE(IFAC), 0.0,
     &                  VAR(IFAC,1), VAR(IFAC,2), 0.0, VSCL, ASCL)
                  ELSE
                     XVAR = VAR(IFAC,1)
                     YVAR = VAR(IFAC,2)
                     ZVAR = VAR(IFAC,3)
                     CALL BL_ROTATE (1, 1, ROTMAT, ZERO3,
     &                  XVAR, YVAR, ZVAR, XVAR, YVAR, ZVAR)
                     CALL VEC (IS3DIM, XNE(IFAC), YNE(IFAC), ZNE(IFAC),
     &                  XVAR, YVAR, ZVAR, VSCL, ASCL)
                  END IF

  120          CONTINUE

               CALL PLTFLU
            END IF
  130    CONTINUE

      ELSE IF ((ISYTYP .EQ. 'SIGMAX') .OR. (ISYTYP .EQ. 'SIGMIN')) THEN

         DEG180 = 2.0 * ATAN (1.0)

         DO 150 IELB = 1, NELBLK
            IF (ISVOK(IELB)) THEN

C            --Set the vector color
c               CALL UGRCOL (IDELB(IELB), BLKCOL)
               CALL UGRCOL (IELB, BLKCOL)

               DO 140 IFAC = LENF(IELB-1)+1, LENF(IELB)
                  IF (IS3DIM) THEN
                     IF (HIDENE(IFAC)) GOTO 140
                  END IF

                  SIGXX = VAR(IFAC,1)
                  SIGYY = VAR(IFAC,2)
                  TAUXY = VAR(IFAC,3)
                  Q1 = 0.5 * (SIGXX + SIGYY)
                  Q2 = 0.5 * (SIGXX - SIGYY)
                  Q3 = SQRT (Q2 * Q2 + TAUXY * TAUXY)
                  SIGMAX = Q1 + Q3
                  SIGMIN = Q1 - Q3
                  IF ((TAUXY .EQ. 0.0) .AND. (SIGMAX-SIGYY .EQ. 0.0))
     &               GOTO 140
                  THETA = ATAN2 (TAUXY, SIGMAX - SIGYY)
                  IF (ISYTYP .EQ. 'SIGMAX') THEN
                     XSIG = SIGMAX * COS (THETA)
                     YSIG = SIGMAX * SIN (THETA)
                  ELSE
                     XSIG = SIGMIN * COS (THETA + DEG180)
                     YSIG = SIGMIN * SIN (THETA + DEG180)
                  END IF

C               --Call vector routine

                  IF (GRABRT ()) RETURN 1

                  IF (.NOT. IS3DIM) THEN
                     CALL VEC (IS3DIM, XNE(IFAC), YNE(IFAC), 0.0,
     &                  XSIG, YSIG, 0.0, VSCL, ASCL)
                  ELSE
                     ZSIG = 0.0
                     CALL BL_ROTATE (1, 1, ROTMAT, ZERO3,
     &                  XSIG, YSIG, ZSIG, XSIG, YSIG, ZSIG)
                     CALL VEC (IS3DIM, XNE(IFAC), YNE(IFAC), ZNE(IFAC),
     &                  XSIG, YSIG, ZNE(IFAC), VSCL, ASCL)
                  END IF

  140          CONTINUE

               CALL PLTFLU
            END IF
  150    CONTINUE
      END IF

      RETURN
      END

C=======================================================================
      SUBROUTINE VECTORN (ISYTYP, VAR, NUM, HIDENE,
     &   XNE, YNE, ZNE, IN2ELB, VECMAX, BLKCOL, IDELB, *)
C=======================================================================

C   --*** VECTOR *** (DETOUR) Plot nodal or element vector
C   --   Written by Amy Gilkey - revised 03/10/88
C   --   D. P. Flanagan, 12/08/83
C   --
C   --VECTOR draws a vector representing 2 or 3 variables for each node
C   --or element.  It processes each node or element by element block.
C   --Only nodes or elements of selected element blocks are drawn.
C   --
C   --Parameters:
C   --   ISYTYP - IN - the vector type (as in MODTYP of /DETOPT/)
C   --   VAR - IN - the vector variable values
C   --   NUM - IN - the number of nodes (if nodal) or faces (if element)
C   --   LENF - IN - the cumulative face counts by element block
C   --      (only if element)
C   --   NLNKF - IN - the number of nodes per face (only if element)
C   --   HIDENE(i) - IN - true iff node i or face i is hidden (3D only)
C   --   XNE, YNE, ZNE - IN - the nodal coordinates (if nodal)
C   --      or the face center coordinates (if element)
C   --   IN2ELB - IN - the element block for each node; (only if nodal)
C   --      <0 if not in any selected element block
C   --      =0 if in more than one selected element block
C   --   VECMAX - IN - the maximum vector variable value, scaled
C   --   BLKCOL - IN - the user selected colors of the element blocks.
C   --                    BLKCOL(0) = 1 if the user defined material
C   --                                colors should be used in mesh plots.
C   --                              = -1 if program selected colors should
C   --                                be used.
C   --                    BLKCOL(i) = the user selected color of element
C   --                               block i:
C   --                                  -2 - no color selected by user.
C   --                                  -1 - black
C   --                                   0 - white
C   --                                   1 - red
C   --                                   2 - green
C   --                                   3 - yellow
C   --                                   4 - blue
C   --                                   5 - cyan
C   --                                   6 - magenta
C   --   * - return statement if cancel function active
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses IS3DIM, NUMNPF of /D3NUMS/
C   --   Uses VECSCL of /ETCOPT/
C   --   Uses ROTMAT of /ROTOPT/
C   --   Uses DTW, VWSCL of /DEVDAT/

      include 'dbnums.blk'
      include 'd3nums.blk'
      include 'etcopt.blk'
      include 'rotopt.blk'
      include 'devdat.blk'

      CHARACTER*(*) ISYTYP
      REAL VAR(NUM,*)
      LOGICAL HIDENE(NUM)
      REAL XNE(NUM), YNE(NUM), ZNE(NUM)
      INTEGER IN2ELB(NUMNPF)
      LOGICAL ISVOK(NELBLK)
      REAL VECMAX
      INTEGER BLKCOL(0:NELBLK)

      LOGICAL GRABRT
      REAL ZERO3(3)
      SAVE ZERO3

      DATA ZERO3 / 0.0, 0.0, 0.0 /

      IF (VECSCL .EQ. 0.0) THEN
         VSCL = 1.0
      ELSE
         VSCL = VECSCL * 0.05 * VECMAX
      END IF
      ASCL = VWSCL

      IF (ISYTYP .EQ. 'NODE') THEN

         DO 110 IELB = 0, NELBLK

C         --Set the vector color
            ITEMP = IELB
            CALL UGRCOL (ITEMP, BLKCOL)

            DO 100 INP = 1, NUMNPF

               IF (IS3DIM) THEN
                  IF (HIDENE(INP)) GOTO 100
               END IF

               IF (IN2ELB(INP) .EQ. IELB) THEN
                  IF (GRABRT ()) RETURN 1

C               --Call vector routine

                  IF (.NOT. IS3DIM) THEN
                     CALL VEC (IS3DIM, XNE(INP), YNE(INP), 0.0,
     &                  VAR(INP,1), VAR(INP,2), 0.0, VSCL, ASCL)
                  ELSE
                     XVAR = VAR(INP,1)
                     YVAR = VAR(INP,2)
                     ZVAR = VAR(INP,3)
                     CALL BL_ROTATE (1, 1, ROTMAT, ZERO3,
     &                  XVAR, YVAR, ZVAR, XVAR, YVAR, ZVAR)
                     CALL VEC (IS3DIM, XNE(INP), YNE(INP), ZNE(INP),
     &                  XVAR, YVAR, ZVAR, VSCL, ASCL)
                  END IF
               END IF

  100       CONTINUE

            CALL PLTFLU
  110    CONTINUE

      END IF

      RETURN
      END

