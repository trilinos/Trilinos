C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
C         
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C=======================================================================
      SUBROUTINE DBPTIM (OPTION, ALLPRT, ALLTIM, NSTEPS, TIMES, WHOTIM)
C=======================================================================
C$Id: dbptim.f,v 1.3 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dbptim.f,v $
CRevision 1.3  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.2  1997/03/20 19:40:14  caforsy
CUpdated Imakefile for Imake 6.1.  Changed printing routines to handle
Clarger problems.
C
CRevision 1.1.1.1  1990/08/14 16:13:56  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:13:55  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:20  gdsjaar
c Initial revision
c 

C   --*** DBPTIM *** (EXOLIB) Print database steps and min/max times
C   --   Written by Amy Gilkey - revised 12/18/87
C   --
C   --DBPTIM displays the number of time steps and the minimum and maximum
C   --time on the database.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --      'N' to print the number of time steps
C   --      'M' to print the minimum and maximum time step times
C   --      'T' to print the time step times
C   --   ALLPRT - IN - the type of time steps to print:
C   --      true to print both whole and history information;
C   --      false to print whole time steps
C   --   ALLTIM - IN - the type of the input time steps:
C   --      true if whole and history information;
C   --      false if whole time steps only
C   --   NSTEPS - IN - the number of time steps
C   --   TIMES - IN - the database time steps
C   --   WHOTIM - IN - true iff TIMES(i) is whole (versus history) time step
C   --      (only if ALLTIM = '*')

C   --Routines Called:
C   --   INTSTR - (STRLIB) Convert integers to strings
C   --   NUMSTR - (STRLIB) Convert numbers to engineering notation
C   --   SQZSTR - (STRLIB) Delete extra blanks from string

      CHARACTER*(*) OPTION
      LOGICAL ALLPRT, ALLTIM
      INTEGER NSTEPS
      REAL TIMES(*)
      LOGICAL WHOTIM(*)

      LOGICAL ISABRT
      CHARACTER*80 STRING
      CHARACTER*5 ISTRA
      CHARACTER*20 RSTR(3)
      REAL RNUM(3)

      WRITE (*, *)

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
         IF (ALLPRT) THEN
            IF (ALLTIM) THEN
               WRITE (STRING, 10000) NSTEPS, NSTEPH
            ELSE
               WRITE (STRING, 10000) NSTEPS
            END IF
10000        FORMAT ('Number of time steps = ', I8, :,
     &         ' (including ', I8, ' history-only)')
         ELSE
            WRITE (STRING, 10010) NSTEPW
10010        FORMAT ('Number of whole time steps = ', I8)
         END IF
         CALL SQZSTR (STRING, LSTR)
         WRITE (*, 10030) STRING(:LSTR)
      END IF

      IF (((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'M') .GT. 0))
     &   .AND. (NSTEPX .GT. 0)) THEN
         IF (ALLTIM .AND. (.NOT. ALLPRT)) THEN
            CALL MINMXL (NSTEPS, WHOTIM, TIMES, TIMMIN, TIMMAX)
         ELSE
            CALL MINMAX (NSTEPS, TIMES, TIMMIN, TIMMAX)
         END IF
         IF (NSTEPX .EQ. 1) THEN
            CALL NUMSTR (1, 4, TIMMIN, RSTR, LSTR)
            IF (ALLPRT) THEN
               WRITE (*, 10030) '   Time = ', RSTR(1)(:LSTR)
            ELSE
               WRITE (*, 10030) '   Whole Time = ', RSTR(1)(:LSTR)
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
     &            '   Minimum whole time = ', RSTR(1)(:LSTR)
               WRITE (*, 10030)
     &            '   Maximum whole time = ', RSTR(2)(:LSTR)
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
               ELSE IF (WHOTIM(I)) THEN
                  WRITE (*, 10020, IOSTAT=IDUM)
     &               ISTRA(:LI), RSTR(1)(:LSTR)
               ELSE
                  WRITE (*, 10020, IOSTAT=IDUM)
     &               ISTRA(:LI), RSTR(1)(:LSTR), '(history-only)'
               END IF
  110       CONTINUE
         END IF
      END IF

      RETURN
10030  FORMAT (1X, 5A)
      END
