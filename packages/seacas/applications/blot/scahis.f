C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: scahis.f,v $
C Revision 1.3  2009/03/25 12:36:47  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  2009/01/22 21:34:22  gdsjaar
C There were several inline dbnums common blocks. Replaced with the
C include so they all have the same size with the added variable types.
C
C Added minor support for nodeset and sideset variables.
C
C It can print the count and the names, but that is all
C at this time.
C
C Revision 1.1  1994/04/07 20:10:42  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:56:54  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SCAHIS (A, VAR, WHOTIM,
     &   VALMIN, ISTMIN, VALMAX, ISTMAX)
C=======================================================================

C   --*** SCAHIS *** (BLOT) Scale all history variables
C   --   Written by Amy Gilkey - revised 04/01/88
C   --
C   --SCAHIS reads the values for the history variables from the database
C   --and finds the minimum and maximum values.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   VAR - SCRATCH - the variable array
C   --   WHOTIM - IN - true iff time step is a whole (versus history) time step
C   --   VALMIN, VALMAX - OUT - the minimum and maximum value for each variable
C   --   ISTMIN, ISTMAX - OUT - the step number of the minimum and maximum
C   --      value for each variable
C   --
C   --Common Variables:
C   --   Uses NVARHI, NSTEPS of /DBNUMS/

      include 'dbnums.blk'

      DIMENSION A(*)
      REAL VAR(NVARHI)
      LOGICAL WHOTIM(*)
      REAL VALMIN(NVARHI), VALMAX(NVARHI)
      INTEGER ISTMIN(NVARHI), ISTMAX(NVARHI)

      CALL DBVIX_BL ('H', 1, IVAR)

      DO 110 ISTEP = 1, NSTEPS

C      --Read the variables

         CALL GETVAR (A, IVAR, -999, ISTEP, NVARHI, VAR)

C      --Find minimum and maximum variable values for each variable

         DO 100 IXVAR = 1, NVARHI

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
