C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DTREAD (A, ISTEP, IDTVAR, NNDVAR, NEDVAR,
     &   LENF, IF2EL, IELBST, ISEVOK, VARNP, VARFAC, VAR, LVARF)
C=======================================================================

C   --*** DTREAD *** (DETOUR) Read variables for time step
C   --   Written by Amy Gilkey - revised 05/11/88
C   --
C   --DTREAD reads the nodal and element variables needed for the time step.
C   --The element variables are converted to face variables.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   ISTEP - IN - the selected time step number
C   --   IDTVAR - IN - the variable numbers
C   --   NNDVAR - IN - the number of nodal variables needed
C   --   NEDVAR - IN - the number of element variables needed
C   --   LENF - IN - the cumulative face counts by element block
C   --   IF2EL - IN - the element number of each face
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   VARNP - OUT - the nodal variable values
C   --      (sized for all needed nodal variables)
C   --   VARFAC - OUT - the face variable values
C   --      (sized for all needed element variables)
C   --   VAR - SCRATCH - array to hold an element variable
C   --
C   --Common Variables:
C   --   Uses NUMNP, NUMEL, NELBLK of /DBNUMS/

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      DIMENSION A(*)
      INTEGER IDTVAR(*)
      INTEGER LENF(0:NELBLK)
      INTEGER IF2EL(*)
      INTEGER IELBST(NELBLK)
      LOGICAL ISEVOK(NELBLK,*)
      REAL VARNP(NUMNP,*)
      REAL VARFAC(LVARF,*)
      REAL VAR(NUMEL)

      CHARACTER TYP

      DO 160 IVAR = MAX (NNDVAR, NEDVAR), 1, -1

C      --Get the variable type

         CALL DBVTYP_BL (IDTVAR(IVAR), TYP, ID)

C      --Read in variable

         IF (TYP .EQ. ' ') THEN

C         --Zero out nodal variable and element variable if needed

            IF (NNDVAR .GE. IVAR) THEN
               DO 110 INP = 1, NUMNP
                  VARNP(INP,IVAR) = 0.0
  110          CONTINUE
            END IF

            IF (NEDVAR .GE. IVAR) THEN
               DO 130 IELB = 1, NELBLK
                  IF (IELBST(IELB) .GT. 0) THEN
                     DO 120 IFAC = LENF(IELB-1)+1, LENF(IELB)
                        VARFAC(IFAC,IVAR) = 0.0
  120                CONTINUE
                  END IF
  130          CONTINUE
            END IF

         ELSE IF (TYP .EQ. 'N') THEN

C         --Get nodal variable

            CALL GTMVAR (A, IDTVAR(IVAR), -999, ISTEP, NUMNPF,
     &         VARNP(1,IVAR))

         ELSE IF (TYP .EQ. 'E') THEN

C            --Get element variable, change to face variable

               CALL GETVAR (A, IDTVAR(IVAR), -1, ISTEP, NUMEL, VAR)

               DO 150 IELB = 1, NELBLK
                  IF ((IELBST(IELB) .GT. 0) .AND. ISEVOK(IELB,ID)) THEN
                     DO 140 IFAC = LENF(IELB-1)+1, LENF(IELB)
                        VARFAC(IFAC,IVAR) = VAR(IF2EL(IFAC))
  140                CONTINUE
                  END IF
  150          CONTINUE

         END IF
  160 CONTINUE

      RETURN
      END
