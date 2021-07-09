C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SCAGLO (A, VAR, WHOTIM,
     &   VALMIN, ISTMIN, VALMAX, ISTMAX)
C=======================================================================

C   --*** SCAGLO *** (BLOT) Scale all global variables
C   --   Written by Amy Gilkey - revised 04/01/88
C   --
C   --SCAGLO reads the values for the global variables from the database
C   --and finds the minimum and maximum values.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   IVAR - IN - the variable index (for GETVAR)
C   --   VAR - SCRATCH - the variable array
C   --   WHOTIM - IN - true iff time step is a whole (versus history) time step
C   --   VALMIN, VALMAX - OUT - the minimum and maximum value for each variable
C   --   ISTMIN, ISTMAX - OUT - the step number of the minimum and maximum
C   --      value for each variable
C   --
C   --Common Variables:
C   --   Uses NVARGL, NSTEPS of /DBNUMS/

      include 'dbnums.blk'

      DIMENSION A(*)
      REAL VAR(NVARGL)
      LOGICAL WHOTIM(*)
      REAL VALMIN(NVARGL), VALMAX(NVARGL)
      INTEGER ISTMIN(NVARGL), ISTMAX(NVARGL)

      CALL DBVIX_BL ('G', 1, IVAR)

      DO 110 ISTEP = 1, NSTEPS
         IF (.NOT. WHOTIM(ISTEP)) GOTO 110

C      --Read the variables

         CALL GETVAR (A, IVAR, -999, ISTEP, NVARGL, VAR)

C      --Find minimum and maximum variable values for variable

         DO 100 IXVAR = 1, NVARGL

            IF (ISTEP .EQ. 1) THEN
               VALMIN(IXVAR) = VAR(IXVAR)
               ISTMIN(IXVAR) = ISTEP
               VALMAX(IXVAR) = VAR(IXVAR)
               ISTMAX(IXVAR) = ISTEP
            ELSE IF (VALMIN(IXVAR) .GT. VAR(IXVAR)) THEN
               VALMIN(IXVAR) = VAR(IXVAR)
               ISTMIN(IXVAR) = ISTEP
            ELSE IF (VALMAX(IXVAR) .LT. VAR(IXVAR)) THEN
               VALMAX(IXVAR) = VAR(IXVAR)
               ISTMAX(IXVAR) = ISTEP
            END IF

  100    CONTINUE
  110 CONTINUE

      RETURN
      END
