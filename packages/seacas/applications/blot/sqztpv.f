C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SQZTPV (NPTIMS, IPTIMS, WHOTIM, NPTS, PLTVAL)
C=======================================================================

C   --*** SQZTPV *** (TPLOT) Compress curves for whole time steps
C   --   Written by Amy Gilkey - revised 11/11/87
C   --
C   --SQZTPV compresses the values for TPLOT curves so that only the values
C   --for whole time steps are in the curve.  Curves of history variables
C   --only are left intact.
C   --
C   --Parameters:
C   --   NPTIMS - IN - the number of selected time steps
C   --   IPTIMS - IN - the selected time steps
C   --   WHOTIM - IN - true iff whole (versus history) time step
C   --   NPTS - OUT - the number of points on each curve
C   --   PLTVAL - IN/OUT - the plot data;
C   --      PLTVAL(x,NTPVAR+1) holds the times if TIMPLT
C   --      PLTVAL(x,NTPVAR+2) holds the compressed times if TIMPLT and needed
C   --
C   --Common Variables:
C   --   Uses NTPCRV, NTPVAR, TIMPLT of /TPVARS/

      include 'params.blk'
      include 'tpvars.blk'

      INTEGER IPTIMS(*)
      LOGICAL WHOTIM(*)
      INTEGER NPTS(NTPVAR+2)
      REAL PLTVAL(NPTIMS,NTPVAR+2)

      CHARACTER TYPX, TYPY

      NSQZ = 0
      NPTIMW = NWHSEL (NPTIMS, IPTIMS, WHOTIM)

      IF (TIMPLT) TYPX = 'H'
      N = 0
      DO 100 NP = 1, NTPCRV
         IF (.NOT. TIMPLT) THEN
            N = N + 1
            CALL DBVTYP_BL (ITVID(N), TYPX, IDUM)
         END IF
         N = N + 1
         CALL DBVTYP_BL (ITVID(N), TYPY, IDUM)
         IF ((NPTIMS .NE. NPTIMW) .AND.
     &      ((TYPX .NE. 'H') .OR. (TYPY .NE. 'H'))) THEN
            NSQZ = NSQZ + 1
            IF (.NOT. TIMPLT) THEN
               CALL SQZLGV (NPTIMS, IPTIMS, WHOTIM,
     &            PLTVAL(1,N-1), NPTS(N-1), PLTVAL(1,N-1))
            END IF
            CALL SQZLGV (NPTIMS, IPTIMS, WHOTIM,
     &         PLTVAL(1,N), NPTS(N), PLTVAL(1,N))
         ELSE
            IF (.NOT. TIMPLT) THEN
               NPTS(N-1) = NPTIMS
            END IF
            NPTS(N) = NPTIMS
         END IF
  100 CONTINUE

      IF ((NSQZ .GT. 0) .AND. TIMPLT) THEN
         CALL SQZLGV (NPTIMS, IPTIMS, WHOTIM,
     &      PLTVAL(1,NTPVAR+1), IDUM, PLTVAL(1,NTPVAR+2))
      END IF

      RETURN
      END
