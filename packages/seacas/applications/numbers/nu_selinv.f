C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SELINV (TIMES, ITMSEL, NINTV)
      DIMENSION TIMES(*)
      LOGICAL ITMSEL(*)
      CHARACTER*16 ENGNOT, STRA, STRB
      CHARACTER*80 STRTMP

      include 'nu_ptim.blk'

      NLAST  = 0
      NUMSEL = 0

C      IFIRST = LOCRL (STMIN, NSTEP, ITMSEL, TIMES)
C      ILAST  = LOCRL (STMAX, NSTEP, ITMSEL, TIMES)
      IFIRST = LOCREA (STMIN, NSTEP, TIMES)
      ILAST  = LOCREA (STMAX, NSTEP, TIMES)
      NBETWN = ILAST - IFIRST
      CALL INILOG (NSTEP, .FALSE., ITMSEL)

      IF (NINTV .EQ. 0) THEN
         CALL PRTERR ('WARNING', 'No time steps selected.')
         RETURN
      ELSE IF (NINTV .LT. 0) THEN
C - Include tmin step
         II = -NINTV - 1
         II = MAX(1, II)
         RINC = DBLE(NBETWN) / DBLE(II)
         IB = IFIRST
         IE = ILAST
      ELSE
         II = NINTV
         RINC = DBLE(NBETWN) / DBLE(II)
         IB = IFIRST + INT(RINC + 0.5)
         IE = ILAST
      END IF

      NUMSEL = 0
      RTIM = DBLE(IB)
  100 CONTINUE
         ITMSEL(INT(RTIM+0.5)) = .TRUE.
         NUMSEL = NUMSEL + 1
         RTIM = MIN( RTIM + RINC, DBLE(IE))
         IF (NUMSEL .LT. ABS(NINTV)) GO TO 100
      ITMSEL(ILAST) = .TRUE.
      STRA = ENGNOT(STMIN,2)
      STRB = ENGNOT(STMAX,2)
      WRITE (STRTMP, 40) NUMSEL, STRA, STRB
      CALL SQZSTR(STRTMP, LSTR)
      WRITE (*, 50) STRTMP(:LSTR)
      RETURN
   40 FORMAT (I5,' Steps Selected from ',A16,' to ',A16)
   50 FORMAT (/,5X,A)
      END
