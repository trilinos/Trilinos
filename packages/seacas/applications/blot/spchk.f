C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SPCHK (PRTOK, OKAY)
C=======================================================================

C   --*** SPCHK *** (SPLOT) Check that SPLOT is a valid program to run
C   --   Written by Amy Gilkey - revised 11/03/87
C   --
C   --SPCHK checks that SPLOT is a valid program.  This is true if
C   --a mesh is defined (nodes, elements, connectivity, etc.) and
C   --there are nodal or element variables in the database.
C   --
C   --If the program is valid, a welcome message is printed.
C   --
C   --Parameters:
C   --   PRTOK - IN - if true, print error messages and welcome messages,
C   --      if false, check but do not print anything
C   --   OKAY - OUT - true iff the program is valid
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP, NUMEL, NVARNP, NVAREL, NSTEPW of /DBNUMS/

      include 'dbnums.blk'

      LOGICAL PRTOK
      LOGICAL OKAY

C   --Check header information from database

      IF ((NUMNP .LE. 0) .OR. (NUMEL .LE. 0) .OR. (NDIM .LE. 0)) THEN
         IF (PRTOK) CALL PRTERR ('CMDERR',
     &      'No mesh is defined')
         GOTO 100
      END IF
      IF ((NDIM .NE. 2) .AND. (NDIM .NE. 3)) THEN
         IF (PRTOK) CALL PRTERR ('CMDERR',
     &      'SPLOT can only process 2D and 3D meshes')
         GOTO 100
      END IF
      IF ((NVARNP + NVAREL) .LE. 0) THEN
         IF (PRTOK) CALL PRTERR ('CMDERR',
     &      'No nodal or element variables are defined')
         GOTO 100
      END IF
      IF (NSTEPW .LE. 0) THEN
         IF (PRTOK) CALL PRTERR ('CMDERR',
     &      'No time steps with nodal or element variables are defined')
         GOTO 100
      END IF

C   --Write greeting message
      IF (PRTOK) WRITE (*, 10000)

      OKAY = .TRUE.
      RETURN

  100 CONTINUE
      OKAY = .FALSE.
      RETURN
10000  FORMAT (' SPLOT - a distance-versus-variable curve plot program')
      END
