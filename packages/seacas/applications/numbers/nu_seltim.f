C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SELTIM (TIMES, ITMSEL)
      DIMENSION TIMES(*)
      LOGICAL ITMSEL(*)
      CHARACTER*16 ENGNOT, STRA, STRB
      CHARACTER*80 STRTMP
      EXTERNAL ENGNOT

      include 'nu_ptim.blk'

C ... TOLER is the tolerance for matching a timestep.  If the difference
C     is less than TOLER, then a match occurs.  TOLER1 = TOLER + 1

      PARAMETER (TOLER = 1.0E-3)
      TOLERP1 = 1.0 + TOLER
      TOLERM1 = 1.0 - TOLER

      CALL INILOG (NSTEP, .FALSE., ITMSEL)

      NLAST  = 0
      LSTSEL = NSTEP
      TIMGET = STMIN
      NUMSEL = 0

      IF (STDEL .EQ. 0.0) THEN
   10    CONTINUE
         NLAST = NLAST + 1
         IF (NLAST .GT. NSTEP) GO TO 30
         IF (TIMES(NLAST) .LE. TOLERP1 * STMAX .AND.
     *       TIMES(NLAST) .GE. TOLERM1 * STMIN ) THEN
            ITMSEL(NLAST) = .TRUE.
            NUMSEL = NUMSEL + 1
            LSTSEL = NLAST
         ELSE IF (TIMES(NLAST) .GT. TOLERP1 * STMAX) THEN
            LSTSEL = NLAST - 1
            GO TO 30
         END IF
         GO TO 10
      ELSE
   20    CONTINUE
         NLAST = NLAST + 1
         IF (NLAST .GT. NSTEP) GO TO 30
         TDELT = ABS(TIMES(NLAST) - TIMGET)
         IF (TDELT .LT. (TOLER * TMAX) .AND.
     *      TIMES(NLAST) .LE. TOLERP1 * STMAX) THEN
            ITMSEL(NLAST) = .TRUE.
            NUMSEL = NUMSEL + 1
            TIMGET = MIN ( TIMES(NLAST) + STDEL, TMAX)
            LSTSEL = NLAST
         ELSE IF (TIMES(NLAST) .GT. TOLERP1 * STMAX) THEN
            LSTSEL = NLAST - 1
            GO TO 30
         END IF
         GO TO 20
      END IF
   30 CONTINUE

      IF (NUMSEL .EQ. 0) THEN
         CALL PRTERR ('WARNING', 'No time steps selected.')
      ELSE
        STRA = ENGNOT(STMIN,2)
        STRB = ENGNOT(STMAX,2)
        WRITE (STRTMP, 40) NUMSEL, STRA, STRB
        CALL SQZSTR(STRTMP, LSTR)
        WRITE (*, 50) STRTMP(:LSTR)
      END IF
      RETURN
   40 FORMAT (I5,' Steps Selected from ',A16,' to ',A16)
   50 FORMAT (/,5X,A)
      END
