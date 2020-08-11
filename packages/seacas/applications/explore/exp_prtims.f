C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE PRTIMS (OPTION, NOUT, NSTEPS, TIMES)
C=======================================================================

C   --*** PRTIMS *** (EXPLORE) Display database time step times
C   --
C   --PRTIMS displays the time for all time steps.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --      'N' to print the number of time steps
C   --      'M' to print the minimum and maximum time step times
C   --      'T' to print the time step times
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NSTEPS - IN - the number of time steps
C   --   TIMES - IN - the times for each time step

      CHARACTER*(*) OPTION
      REAL TIMES(*)

      CHARACTER*80 STRING
      CHARACTER*5 ISTRA
      CHARACTER*64 RSTR(3)
      REAL RNUM(3)
      INTEGER LSTR, NPREC
      INTEGER GETPRC

      IF (NOUT .GT. 0) WRITE (NOUT, 10030)

      IF (NOUT .GT. 0) THEN
         WRITE (NOUT, *)
      ELSE
         WRITE (*, *)
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0)) THEN
         WRITE (STRING, 10000, IOSTAT=IDUM) NSTEPS
10000    FORMAT ('Number of time steps = ', I12)
         CALL SQZSTR (STRING, LSTR)
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10050) STRING(:LSTR)
         ELSE
            WRITE (*, 10050) STRING(:LSTR)
         END IF
      END IF

      NPREC=GETPRC()

      IF (((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'M') .GT. 0))
     &   .AND. (NSTEPS .GT. 0)) THEN
         CALL MINMAX (NSTEPS, TIMES, TIMMIN, TIMMAX)
         IF (NSTEPS .EQ. 1) THEN
            CALL NUMSTR1 (NPREC, TIMMIN, RSTR, LSTR)
            WRITE (STRING, 10040) 'Time = ', RSTR(1)(:LSTR)
            LSTR = LENSTR (STRING)
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10060) STRING(:LSTR)
            ELSE
               WRITE (*, 10060) STRING(:LSTR)
            END IF
         ELSE
            RNUM(1) = TIMMIN
            RNUM(2) = TIMMAX
            CALL NUMSTR (2, NPREC, RNUM, RSTR, LSTR)
            WRITE (STRING, 10040)
     &           'Minimum time = ', RSTR(1)(:LSTR)
            LSTR = LENSTR (STRING)
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10060) STRING(:LSTR)
            ELSE
               WRITE (*, 10060) STRING(:LSTR)
            END IF
            WRITE (STRING, 10040)
     &           'Maximum time = ', RSTR(2)(:LSTR)
            LSTR = LENSTR (STRING)
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10060) STRING(:LSTR)
            ELSE
               WRITE (*, 10060) STRING(:LSTR)
            END IF
         END IF
      END IF

      IF (((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'T') .GT. 0))
     &   .AND. (NSTEPS .GT. 0)) THEN
         CALL INTSTR (1, 0, NSTEPS, ISTRA, LISTR)
         RNUM(2) = TIMES(1)
         IF ((RNUM(2) .EQ. 0.0) .AND. (NSTEPS .GT. 1))
     &      RNUM(2) = TIMES(2)
         RNUM(3) = TIMES(NSTEPS)
         DO 110 I = 1, NSTEPS
           CALL INTSTR (1, LISTR, I, ISTRA, LI)
           RNUM(1) = TIMES(I)
           CALL NUMSTR (3, NPREC, RNUM, RSTR, LSTR)
           WRITE (STRING, 10020) ISTRA(:LI), RSTR(1)(:LSTR)
10020      FORMAT ('Step ', A, ')', 3X, A, :, 3X, A)
           LSTR = LENSTR (STRING)
           IF (NOUT .GT. 0) THEN
              WRITE (NOUT, 10060) STRING(:LSTR)
           ELSE
              WRITE (*, 10060) STRING(:LSTR)
           END IF
 110     CONTINUE
      END IF

      RETURN

10030  FORMAT (/, 1X, 'TIME STEP TIMES')
10040  FORMAT (5A)
10050  FORMAT (1X, 5A)
10060  FORMAT (4X, 5A)
      END
