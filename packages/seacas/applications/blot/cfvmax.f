C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: cfvmax.f,v $
C Revision 1.3  2009/03/25 12:36:42  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  1997/03/20 21:23:52  caforsy
C Update Imakefile for Imake 6.1
C
C Revision 1.1  1994/04/07 19:55:26  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:47:58  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE CFVMAX (NUMFAC, VARFAC, NMIN, NMAX, FMIN, FMAX)
C=======================================================================

C   --*** CFVMAX *** (DETOUR) Calculate face variable min/max and count number
C   --   Written by Amy Gilkey - revised 03/09/88
C   --
C   --CFVMAX calculates the minimum and maximum value of a face variable
C   --and counts the number of occurrences of the minimum and maximum.
C   --
C   --Parameters:
C   --   NUMFAC - IN - the number of surface faces
C   --   VARFAC - IN - the face variable values
C   --   NMIN, NMAX - OUT - the number of variables values matching the
C   --      minimum and the maximum
C   --   FMIN, FMAX - OUT - the face variable minimums and maximums

      REAL VARFAC(*)

      LOGICAL INIT

      NMIN = 0
      NMAX = 0
      FMIN = 0.0
      FMAX = 0.0
      INIT = .TRUE.
      DO 100 IFAC = 1, NUMFAC
         IF (INIT) THEN
            FMIN = VARFAC(IFAC)
            NMIN = 1
            FMAX = VARFAC(IFAC)
            NMAX = 1
            INIT = .FALSE.
         ELSE
            IF (FMIN .GE. VARFAC(IFAC)) THEN
               IF (FMIN .EQ. VARFAC(IFAC)) THEN
                  NMIN = NMIN + 1
               ELSE
                  FMIN = VARFAC(IFAC)
                  NMIN = 1
               END IF
            END IF
            IF (FMAX .LE. VARFAC(IFAC)) THEN
               IF (FMAX .EQ. VARFAC(IFAC)) THEN
                  NMAX = NMAX + 1
               ELSE
                  FMAX = VARFAC(IFAC)
                  NMAX = 1
               END IF
            END IF
         END IF
  100 CONTINUE

      RETURN
      END

