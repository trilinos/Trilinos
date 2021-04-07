C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CNVMAX (NUMNPF, VARNP, IN2ELB, NMIN, NMAX, FMIN, FMAX)
C=======================================================================

C   --*** CNVMAX *** (DETOUR) Calculate nodal variable min/max and count number
C   --   Written by Amy Gilkey - revised 03/09/88
C   --
C   --CNVMAX calculates the minimum and maximum value of a nodal variable
C   --and counts the number of occurrences of the minimum and maximum.
C   --
C   --Parameters:
C   --   NUMNPF - IN - the number of nodes
C   --   VARNP - IN - the nodal variable values
C   --   IN2ELB - IN - the element block for each node;
C   --      <0 if not in any selected element block
C   --      =0 if in more than one selected element block
C   --   NMIN, NMAX - OUT - the number of variables values matching the
C   --      minimum and the maximum
C   --   FMIN, FMAX - OUT - the nodal variable minimums and maximums

      REAL VARNP(*)
      INTEGER IN2ELB(*)

      LOGICAL INIT

      NMIN = 0
      NMAX = 0
      FMIN = 0.0
      FMAX = 0.0
      INIT = .TRUE.
      DO 100 INP = 1, NUMNPF
         IF (IN2ELB(INP) .GE. 0) THEN
            IF (INIT) THEN
               FMIN = VARNP(INP)
               NMIN = 1
               FMAX = VARNP(INP)
               NMAX = 1
               INIT = .FALSE.
            ELSE
               IF (FMIN .GE. VARNP(INP)) THEN
                  IF (FMIN .EQ. VARNP(INP)) THEN
                     NMIN = NMIN + 1
                  ELSE
                     FMIN = VARNP(INP)
                     NMIN = 1
                  END IF
               END IF
               IF (FMAX .LE. VARNP(INP)) THEN
                  IF (FMAX .EQ. VARNP(INP)) THEN
                     NMAX = NMAX + 1
                  ELSE
                     FMAX = VARNP(INP)
                     NMAX = 1
                  END IF
               END IF
            END IF
         END IF
  100 CONTINUE

      RETURN
      END
