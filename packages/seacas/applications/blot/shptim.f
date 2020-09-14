C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SHPTIM (WHONLY, NPTIMS, IPTIMS, TIMES, WHOTIM)
C=======================================================================

C   --*** SHPTIM *** (TIMSEL) Display time request
C   --   Written by Amy Gilkey - revised 02/01/88
C   --
C   --SHPTIM displays the selected times.
C   --
C   --Parameters:
C   --   WHONLY - IN - true iff only whole times may be selected
C   --   NPTIMS - IN - the number of selected times
C   --   IPTIMS - IN - the selected time steps
C   --   TIMES - IN - the database times
C   --   WHOTIM - IN - true iff TIMES(i) is whole (versus history) time step

C   --Routines Called:
C   --   LENSTR - (STRLIB) Find string length
C   --   NUMSTR - (STRLIB) Convert numbers to engineering notation
C   --   SQZSTR - (STRLIB) Delete extra blanks from string

      LOGICAL WHONLY
      INTEGER NPTIMS
      INTEGER IPTIMS(*)
      REAL TIMES(*)
      LOGICAL WHOTIM(*)

      CHARACTER*5 ISTRA
      CHARACTER*20 RSTR(3)
      REAL RNUM(3)

      WRITE (*, *)

      IF (.NOT. WHONLY) THEN
         NPTIMW = NWHSEL (NPTIMS, IPTIMS, WHOTIM)
      ELSE
         NPTIMW = NPTIMS
      END IF
      NPTIMH = NPTIMS - NPTIMW

      CALL INTSTR (1, 0, NPTIMS, ISTRA, LSTRA)
      WRITE (*, 10010)
     &  ISTRA(:LSTRA), ' selected time steps'

      IF (NPTIMS .GT. 0) THEN
         MXSTEP = 0
         DO 100 I = 1, NPTIMS
            MXSTEP = MAX (MXSTEP, IPTIMS(I))
  100    CONTINUE
         CALL INTSTR (1, 0, MXSTEP, ISTRA, LSTRA)
         RNUM(2) = TIMES(IPTIMS(1))
         IF ((RNUM(2) .EQ. 0.0) .AND. (NPTIMS .GT. 1))
     &      RNUM(2) = TIMES(IPTIMS(2))
         RNUM(3) = TIMES(IPTIMS(NPTIMS))

         DO 110 I = 1, NPTIMS
            CALL INTSTR (1, LSTRA, IPTIMS(I), ISTRA, L)
            RNUM(1) = TIMES(IPTIMS(I))
            CALL NUMSTR (3, 4, RNUM, RSTR, LSTR)
            WRITE (*, 10000, IOSTAT=IDUM) I,
     &        ISTRA(:L), RSTR(1)(:LSTR)
  110    CONTINUE
      END IF

      RETURN
10000 FORMAT (1X, I5, ')', 3X, '(step ', A, ')', 3X, A)
10010 FORMAT (1X, 5A)
      END
