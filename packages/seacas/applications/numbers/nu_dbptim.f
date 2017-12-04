C    Copyright(C) 1988-2017 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C              
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C                            
C    * Neither the name of NTESS nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C                                                    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C=======================================================================
      SUBROUTINE DBPTIM (OPTION, NSTEPS, TIMES)
C=======================================================================
C$Id: dbptim.f,v 1.2 2004/06/29 18:05:32 gdsjaar Exp $
C$Log: dbptim.f,v $
CRevision 1.2  2004/06/29 18:05:32  gdsjaar
CGeneral cleanup. Remove unused labels and variables.
C
CRevision 1.1  1999/02/16 21:37:59  gdsjaar
CConverted to read exodusII database format.  Somewhat tested, not
Cready for production yet.
C
CRevision 1.2  1997/03/20 20:55:05  caforsy
CUpdated Imakefile for Imake 6.1. Added changed to amod in order to
Cport to tflop machine
C
CRevision 1.1  1995/10/03 21:43:39  mksmith
CAdding files new to algII for algebra2
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
C   --                 'N' to print the number of time steps
C   --                 'M' to print the minimum and maximum time step times
C   --                 'T' to print the time step times
C   --   NSTEPS - IN - the number of time steps
C   --   TIMES  - IN - the database time steps

C   --Routines Called:
C   --   INTSTR - (STRLIB) Convert integers to strings
C   --   NUMSTR - (STRLIB) Convert numbers to engineering notation
C   --   SQZSTR - (STRLIB) Delete extra blanks from string

      CHARACTER*(*) OPTION
      INTEGER       NSTEPS
      REAL          TIMES(*)
      CHARACTER*80  STRING
      CHARACTER*5   ISTRA
      CHARACTER*20  RSTR(3)
      REAL          RNUM(3)
C     Logical flags for options
      LOGICAL ALL

      ALL = OPTION .EQ. '*'
      WRITE (*, *)

      IF (ALL .OR. (INDEX (OPTION, 'N') .GT. 0)) THEN
         WRITE (STRING, 10000) NSTEPS
10000    FORMAT ('Number of time steps = ', I8)
         CALL SQZSTR (STRING, LSTR)
         WRITE (*, 10030) STRING(:LSTR)
      END IF

      IF (ALL .OR. (INDEX (OPTION, 'M') .GT. 0)
     &    .AND. (NSTEPS .GT. 0)) THEN
         CALL MINMAX (NSTEPS, TIMES, TIMMIN, TIMMAX)
         
         RNUM(1) = TIMMIN
         RNUM(2) = TIMMAX
         CALL NUMSTR (2, 4, RNUM, RSTR, LSTR)
         WRITE (*, 10030)
     &      '   Minimum time = ', RSTR(1)(:LSTR)
         WRITE (*, 10030)
     &      '   Maximum time = ', RSTR(2)(:LSTR)
      END IF

      IF (ALL .OR. (INDEX (OPTION, 'T') .GT. 0)
     &    .AND. (NSTEPS .GT. 0)) THEN
         CALL INTSTR (1, 0, NSTEPS, ISTRA, LISTR)
         RNUM(2) = TIMES(1)
         IF ((RNUM(2) .EQ. 0.0) .AND. (NSTEPS .GT. 1))
     &      RNUM(2) = TIMES(2)
         RNUM(3) = TIMES(NSTEPS)
         DO 110 I = 1, NSTEPS
            CALL INTSTR (1, LISTR, I, ISTRA, LI)
            RNUM(1) = TIMES(I)
            CALL NUMSTR (3, 4, RNUM, RSTR, LSTR)
            WRITE (*, 10020) ISTRA(:LI), RSTR(1)(:LSTR)
10020       FORMAT (1X, 3X, 'Step ', A, ')', 3X, A, :, 3X, A)
  110       CONTINUE
      END IF

      RETURN
10030  FORMAT (1X, 5A)
      END
