C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE NXMMAX (ISTEP, MMNUM, IVAR, NVAR, NUM, VAR,
     &   OLDMIN, OLDMAX, XMIN, XMAX,
     &   MINSTE, MAXSTE, MINNE, MAXNE, NMIN, NMAX)
C=======================================================================

C   --*** NXMMAX *** (EXPLORE) Find the min/max for this time step
C   --
C   --NXMMAX finds an incremental minimum and maximum of a variable
C   --for this time step.  It does not reset the min/max if not on the
C   --first step.  It may find the min/max that are greater/less than
C   --the previous values.
C   --
C   --Parameters:
C   --   ISTEP - IN - the current step number, previous XMIN and XMAX
C   --      are overwritten if <=1
C   --   MMNUM - IN - number of sequential min/max requests for this
C   --      variable, if >1 then find next min/max
C   --   IVAR - IN - min/max variable number
C   --   NUM - IN - the number of nodes or elements or 1 for global
C   --   VAR - IN - the variables for current time step
C   --   OLDMIN, OLDMAX - IN/OUT - the last minimum and maximum values,
C   --      only if MMNUM > 1
C   --   XMIN, XMAX - IN/OUT - the minimum and maximum values
C   --   MINSTE, MAXSTE - IN/OUT - the step number where the minimum and
C   --      maximum were found
C   --   MINNE, MAXNE - IN/OUT - the node or element number where the
C   --      minimum and maximum were found
C   --   NMIN, NMAX - IN/OUT - the number of values equal to the minimum
C   --      and maximum were found

      REAL VAR(NUM, NVAR)

      IF (ISTEP .LE. 1) THEN
         XMIN = 1.0E36
         XMAX = - 1.0E36
         IF (MMNUM .LE. 1) THEN
            OLDMIN = - 1.0E36
            OLDMAX = 1.0E36
         ELSE
            OLDMIN = XMAX
            OLDMAX = XMIN
         END IF
      END IF

      DO 100 I = 1, NUM
         IF ((XMIN .GT. VAR(I,IVAR))
     &      .AND. (VAR(I,IVAR) .GT. OLDMIN)) THEN
            XMIN = VAR(I,IVAR)
            MINSTE = ISTEP
            MINNE = I
            NMIN = 1
         ELSE IF (XMIN .EQ. VAR(I,IVAR)) THEN
            NMIN = NMIN + 1
         END IF

         IF ((XMAX .LT. VAR(I,IVAR))
     &      .AND. (VAR(I,IVAR) .LT. OLDMAX)) THEN
            XMAX = VAR(I,IVAR)
            MAXSTE = ISTEP
            MAXNE = I
            NMAX = 1
         ELSE IF (XMAX .EQ. VAR(I,IVAR)) THEN
            NMAX = NMAX + 1
         END IF
  100 CONTINUE

      RETURN
      END
