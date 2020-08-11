C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SCANOD (A, IVAR, VAR, WHOTIM, XN, YN, ZN,
     &   VALMIN, NUMMIN, XYZMIN, ISTMIN, VALMAX, NUMMAX, XYZMAX, ISTMAX)
C=======================================================================

C   --*** SCANOD *** (BLOT) Scale nodal variable
C   --   Written by Amy Gilkey - revised 04/01/88
C   --
C   --SCANOD reads the values for the nodal variable from the database
C   --and finds the minimum and maximum value.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   IVAR - IN - the variable index (for GETVAR)
C   --   VAR - SCRATCH - the variable array
C   --   WHOTIM - IN - true iff time step is a whole (versus history) time step
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   VALMIN, VALMAX - OUT - the minimum and maximum value
C   --   NUMMIN, NUMMAX - OUT - the node number of the minimum and maximum value
C   --   XYZMIN, XYZMAX - OUT - the coordinates of NUMMIN, NUMMAX
C   --   ISTMIN, ISTMAX - OUT - the step number of the minimum and maximum value
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP, NVARNP, NSTEPS of /DBNUMS/

      include 'dbnums.blk'

      DIMENSION A(*)
      REAL VAR(*)
      LOGICAL WHOTIM(*)
      REAL XN(*), YN(*), ZN(*)
      REAL VALMIN, VALMAX
      INTEGER NUMMIN, NUMMAX
      REAL XYZMIN(3), XYZMAX(3)
      INTEGER ISTMIN, ISTMAX

      DO 110 ISTEP = 1, NSTEPS
         IF (.NOT. WHOTIM(ISTEP)) GOTO 110

C      --Read the variables

         CALL GETVAR (A, IVAR, -999, ISTEP, NUMNP, VAR)

C      --Find minimum and maximum variable values for variable

         IF (ISTEP .EQ. 1) THEN
            INP = 1
            VALMIN = VAR(INP)
            NUMMIN = INP
            ISTMIN = ISTEP
            VALMAX = VAR(INP)
            NUMMAX = INP
            ISTMAX = ISTEP
         END IF

         DO 100 INP = 1, NUMNP
            IF (VALMIN .GT. VAR(INP)) THEN
               VALMIN = VAR(INP)
               NUMMIN = INP
               ISTMIN = ISTEP
            ELSE IF (VALMAX .LT. VAR(INP)) THEN
               VALMAX = VAR(INP)
               NUMMAX = INP
               ISTMAX = ISTEP
            END IF
  100    CONTINUE
  110 CONTINUE

      DO 120 I = 1, 3
         XYZMIN(I) = 0.0
  120 CONTINUE
      IF (NDIM .GE. 1) XYZMIN(1) = XN(NUMMIN)
      IF (NDIM .GE. 2) XYZMIN(2) = YN(NUMMIN)
      IF (NDIM .GE. 3) XYZMIN(3) = ZN(NUMMIN)
      DO 130 I = 1, 3
         XYZMAX(I) = 0.0
  130 CONTINUE
      IF (NDIM .GE. 1) XYZMAX(1) = XN(NUMMAX)
      IF (NDIM .GE. 2) XYZMAX(2) = YN(NUMMAX)
      IF (NDIM .GE. 3) XYZMAX(3) = ZN(NUMMAX)

      RETURN
      END
