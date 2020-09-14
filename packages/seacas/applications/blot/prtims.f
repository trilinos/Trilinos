C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PRTIMS (OPTION, NOUT, ALLPRT, ALLTIM,
     &   NSTEPS, TIMES, WHOTIM)
C=======================================================================

C   --*** PRTIMS *** (BLOT) Display database time step times
C   --   Written by Amy Gilkey - revised 12/23/87
C   --
C   --PRTIMS displays the time for all time steps.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --      'N' to print the number of time steps
C   --      'M' to print the minimum and maximum time step times
C   --      'T' to print the time step times
C   --   NOUT - IN - the output file, <=0 for standard
C   --   ALLPRT - IN - the type of time steps to print:
C   --      true to print both whole and history information;
C   --      false to print whole time steps
C   --   ALLTIM - IN - the type of the input time steps:
C   --      true if whole and history information;
C   --      false if whole time steps only
C   --   NSTEPS - IN - the number of time steps
C   --   TIMES - IN - the times for each time step
C   --   WHOTIM - IN - true iff TIMES(i) is whole (versus history) time step

      CHARACTER*(*) OPTION
      LOGICAL ALLPRT, ALLTIM
      REAL TIMES(*)
      LOGICAL WHOTIM(*)

      LOGICAL ISABRT
      CHARACTER*80 STRING
      CHARACTER*5 ISTRA
      CHARACTER*20 RSTR(3)
      REAL RNUM(3)

      IF (NOUT .GT. 0) WRITE (NOUT, 10030)

      IF (NOUT .GT. 0) THEN
         WRITE (NOUT, *)
      ELSE
         WRITE (*, *)
      END IF

      IF (ALLTIM) THEN
         NSTEPW = NUMEQL (.TRUE., NSTEPS, WHOTIM)
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
        WRITE (STRING, 10000, IOSTAT=IDUM) NSTEPS
10000     FORMAT ('Number of time steps = ', I5)
         CALL SQZSTR (STRING, LSTR)
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10050) STRING(:LSTR)
         ELSE
            WRITE (*, 10050) STRING(:LSTR)
         END IF
      END IF

      IF (((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'M') .GT. 0))
     &   .AND. (NSTEPX .GT. 0)) THEN
         IF (ALLTIM .AND. (.NOT. ALLPRT)) THEN
            CALL MINMXL (NSTEPS, WHOTIM, TIMES, TIMMIN, TIMMAX)
         ELSE
            CALL MINMAX (NSTEPS, TIMES, TIMMIN, TIMMAX)
         END IF
         IF (NSTEPX .EQ. 1) THEN
            CALL NUMSTR1(8, TIMMIN, RSTR, LSTR)
            IF (ALLPRT) THEN
               WRITE (STRING, 10040) 'Time = ', RSTR(1)(:LSTR)
            ELSE
               WRITE (STRING, 10040) 'Whole Time = ', RSTR(1)(:LSTR)
            END IF
            LSTR = LENSTR (STRING)
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10060) STRING(:LSTR)
            ELSE
               WRITE (*, 10060) STRING(:LSTR)
            END IF
         ELSE
            RNUM(1) = TIMMIN
            RNUM(2) = TIMMAX
            CALL NUMSTR (2, 8, RNUM, RSTR, LSTR)
            IF (ALLPRT) THEN
               WRITE (STRING, 10040)
     &            'Minimum time = ', RSTR(1)(:LSTR)
            ELSE
               WRITE (STRING, 10040)
     &            'Minimum whole time = ', RSTR(1)(:LSTR)
            END IF
            LSTR = LENSTR (STRING)
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10060) STRING(:LSTR)
            ELSE
               WRITE (*, 10060) STRING(:LSTR)
            END IF
            IF (ALLPRT) THEN
              LSTR = LENSTR(RSTR(2))
              WRITE (STRING, 10040)
     &          'Maximum time = ', RSTR(2)(:LSTR)
            ELSE
              LSTR = LENSTR(RSTR(2))
              WRITE (STRING, 10040)
     &          'Maximum whole time = ', RSTR(2)(:LSTR)
            END IF
            LSTR = LENSTR (STRING)
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10060) STRING(:LSTR)
            ELSE
               WRITE (*, 10060) STRING(:LSTR)
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
               IF (WHOTIM(I)) THEN
                  IF (ISABRT ()) RETURN
                  CALL INTSTR (1, LISTR, I, ISTRA, LI)
                  RNUM(1) = TIMES(I)
                  CALL NUMSTR (3, 8, RNUM, RSTR, LSTR)
                  WRITE (STRING, 10020)
     &               ISTRA(:LI), RSTR(1)(:LSTR)
10020              FORMAT ('Step ', A, ')', 3X, A, :, 3X, A)
                  LSTR = LENSTR (STRING)
                  IF (NOUT .GT. 0) THEN
                     WRITE (NOUT, 10060) STRING(:LSTR)
                  ELSE
                     WRITE (*, 10060) STRING(:LSTR)
                  END IF
               END IF
  100       CONTINUE
         ELSE
            DO 110 I = 1, NSTEPS
               CALL INTSTR (1, LISTR, I, ISTRA, LI)
               RNUM(1) = TIMES(I)
               CALL NUMSTR (3, 8, RNUM, RSTR, LSTR)
               IF (.NOT. ALLTIM) THEN
                  WRITE (STRING, 10020)
     &               ISTRA(:LI), RSTR(1)(:LSTR)
               ELSE IF (WHOTIM(I)) THEN
                  WRITE (STRING, 10020)
     &               ISTRA(:LI), RSTR(1)(:LSTR)
               ELSE
                  WRITE (STRING, 10020)
     &               ISTRA(:LI), RSTR(1)(:LSTR), '(history-only)'
               END IF
               LSTR = LENSTR (STRING)
               IF (NOUT .GT. 0) THEN
                  WRITE (NOUT, 10060) STRING(:LSTR)
               ELSE
                  WRITE (*, 10060) STRING(:LSTR)
               END IF
  110       CONTINUE
         END IF
      END IF

      RETURN

10030  FORMAT (/, 1X, 'TIME STEP TIMES')
10040  FORMAT (5A)
10050  FORMAT (1X, 5A)
10060  FORMAT (4X, 5A)
      END
