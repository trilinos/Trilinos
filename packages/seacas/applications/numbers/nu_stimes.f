C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE STIMES (OPTION, ALLPRT, ALLTIM, NSTEPS, TIMES, SELTIM)
C=======================================================================

C   --*** STIMES *** Print selected database steps and min/max times
C   --   Modified from DBPTIM Written by Amy Gilkey - revised 12/18/87
C   --
C   --STIMES displays the number of selected time steps and the minimum
C   -- and maximum selected time on the database.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --      'N' to print the number of time steps
C   --      'M' to print the minimum and maximum time step times
C   --      'T' to print the time step times
C   --   ALLPRT - IN - the type of time steps to print:
C   --      true to print both selected and non-selected information;
C   --      false to print selected time steps
C   --   ALLTIM - IN - the type of the input time steps:
C   --      true if selected and non-selected information;
C   --      false if selected time steps only
C   --   NSTEPS - IN - the number of time steps
C   --   TIMES - IN - the database time steps
C   --   SELTIM - IN - true iff TIMES(i) is selected time step
C   --      (only if ALLTIM = '*')

C   --Routines Called:
C   --   INTSTR - (STRLIB) Convert integers to strings
C   --   NUMSTR - (STRLIB) Convert numbers to engineering notation
C   --   SQZSTR - (STRLIB) Delete extra blanks from string

      CHARACTER*(*) OPTION
      LOGICAL ALLPRT, ALLTIM
      INTEGER NSTEPS
      REAL TIMES(*)
      LOGICAL SELTIM(*)

      LOGICAL ISABRT
      CHARACTER*80 STRING
      CHARACTER*5 ISTRA
      CHARACTER*20 RSTR(3)
      REAL RNUM(3)

      WRITE (*, *)

      IF (ALLTIM) THEN
         NSTEPW = NUMEQL (.TRUE., NSTEPS, SELTIM)
      ELSE
         NSTEPW = NSTEPS
      END IF
      NSTEPH = NSTEPS - NSTEPW
      IF (ALLPRT) THEN
         NSTEPX = NSTEPS
      ELSE
         NSTEPX = NSTEPW
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0)) THEN
         IF (ALLPRT) THEN
            IF (ALLTIM) THEN
               WRITE (STRING, 10000) NSTEPS, NSTEPH
            ELSE
               WRITE (STRING, 10000) NSTEPS
            END IF
10000        FORMAT ('Number of time steps = ', I5, :,
     &         ' (including ', I5, ' non-selected)')
         ELSE
            WRITE (STRING, 10010) NSTEPW
10010        FORMAT ('Number of selected time steps = ', I5)
         END IF
         CALL SQZSTR (STRING, LSTR)
         WRITE (*, 10030) STRING(:LSTR)
      END IF

      IF (((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'M') .GT. 0))
     &   .AND. (NSTEPX .GT. 0)) THEN
         IF (ALLTIM .AND. (.NOT. ALLPRT)) THEN
            CALL MINMXL (NSTEPS, SELTIM, TIMES, TIMMIN, TIMMAX)
         ELSE
            CALL MINMAX (NSTEPS, TIMES, TIMMIN, TIMMAX)
         END IF
         IF (NSTEPX .EQ. 1) THEN
            CALL NUMSTR1 (4, TIMMIN, RSTR, LSTR)
            IF (ALLPRT) THEN
               WRITE (*, 10030) '   Time = ', RSTR(1)(:LSTR)
            ELSE
               WRITE (*, 10030) '   Selected Time = ', RSTR(1)(:LSTR)
            END IF
         ELSE
            RNUM(1) = TIMMIN
            RNUM(2) = TIMMAX
            CALL NUMSTR (2, 4, RNUM, RSTR, LSTR)
            IF (ALLPRT) THEN
               WRITE (*, 10030)
     &            '   Minimum time = ', RSTR(1)(:LSTR)
               WRITE (*, 10030)
     &            '   Maximum time = ', RSTR(2)(:LSTR)
            ELSE
               WRITE (*, 10030)
     &            '   Minimum selected time = ', RSTR(1)(:LSTR)
               WRITE (*, 10030)
     &            '   Maximum selected time = ', RSTR(2)(:LSTR)
            END IF
         END IF
      END IF

      IF (((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'T') .GT. 0))
     &   .AND. (NSTEPX .GT. 0)) THEN
         CALL INTSTR (1, 0, NSTEPS, ISTRA, LISTR)
         RNUM(2) = TIMES(1)
         IF ((RNUM(2) .EQ. 0.0) .AND. (NSTEPS .GT. 1))
     &      RNUM(2) = TIMES(2)
         RNUM(3) = TIMES(NSTEPS)
         IF (ALLTIM .AND. (.NOT. ALLPRT)) THEN
            DO 100 I = 1, NSTEPS
               IF (SELTIM(I)) THEN
                  IF (ISABRT ()) RETURN
                  CALL INTSTR (1, LISTR, I, ISTRA, LI)
                  RNUM(1) = TIMES(I)
                  CALL NUMSTR (3, 4, RNUM, RSTR, LSTR)
                  WRITE (*, 10020, IOSTAT=IDUM)
     &               ISTRA(:LI), RSTR(1)(:LSTR)
               END IF
  100       CONTINUE
         ELSE
            DO 110 I = 1, NSTEPS
               CALL INTSTR (1, LISTR, I, ISTRA, LI)
               RNUM(1) = TIMES(I)
               CALL NUMSTR (3, 4, RNUM, RSTR, LSTR)
               IF (.NOT. ALLTIM) THEN
                  WRITE (*, 10020, IOSTAT=IDUM)
     &               ISTRA(:LI), RSTR(1)(:LSTR)
10020              FORMAT (1X, 3X, 'Step ', A, ')', 3X, A, :, 3X, A)
               ELSE IF (SELTIM(I)) THEN
                  WRITE (*, 10020, IOSTAT=IDUM)
     &               ISTRA(:LI), RSTR(1)(:LSTR)
               ELSE
                  WRITE (*, 10020, IOSTAT=IDUM)
     &               ISTRA(:LI), RSTR(1)(:LSTR), '(non-selected)'
               END IF
  110       CONTINUE
         END IF
      END IF

      RETURN
10030  FORMAT (1X, 5A)
      END
