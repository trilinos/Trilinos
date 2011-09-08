C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
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
C    * Neither the name of Sandia Corporation nor the names of its
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
C    
C=======================================================================
      SUBROUTINE TOSTEP (NSTEP, NUMELB, IDELB, ISEVOK,
     &   TIME, VARGL, VARNP, VAREL)
C=======================================================================

C   --*** TOSTEP *** (GROPE) Move database to given time step number
C   --
C   --TOSTEP reads the given time step variables (NSTEP), setting
C   --the time step to the new current time step.  A message is
C   --displayed if the time step is too large.
C   --
C   --Parameters:
C   --   NSTEP - IN/OUT - the number of the time step to be read;
C   --      the returned number of the time step actually read
C   --   NUMELB - IN - the number of elements per block
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(i,j) is NOT 0
C   --   TIME - OUT - the time step time
C   --   VARGL - IN/OUT - the global variables for the time step
C   --   VARNP - IN/OUT - the nodal variables for the time step
C   --   VAREL - IN/OUT - the element variables for the time step
C   --
C   --Common Variables:
C   --   Sets and uses NCSTEP of /DBASE/
C   --   Uses NELBLK, NQAREC, NINFO, NVARHI, NVARGL, NVARNP, NVAREL, NSTEPS
C   --      of /DBNUMS/

      INCLUDE 'dbase.blk'
      INCLUDE 'dbnums.blk'

      INTEGER NUMELB(*)
      INTEGER IDELB(*)
      INTEGER ISEVOK(NVAREL,*)
      REAL VARGL(*), VARNP(*), VAREL(*)

      CHARACTER*80 ERRMSG
      CHARACTER*5 STRA

      LOGICAL FIRST
      SAVE FIRST

      INTEGER NCSTEP
      REAL CTIME
      save ctime
C      --NCSTEP - the current step number
C      --CTIME - the current time

      DATA FIRST / .TRUE. /

      IF (FIRST) THEN
         NCSTEP = -999
         FIRST = .FALSE.
      END IF

      IF (NSTEP .LE. 0) THEN
         CALL INTSTR (1, 0, NSTEP, STRA, LSTRA)
         CALL PRTERR ('CMDERR',
     &      'Invalid time step number ' // STRA(:LSTRA))
         GOTO 150
      END IF

      if (nstep .ne. ncstep) then
         CALL RDSTEP (NDB, NSTEP,
     &        NUMELB, IDELB, ISEVOK,
     &        CTIME, VARGL, VARNP, VAREL,
     $        max(1,nvarel))
         ncstep = nstep
      END IF
      TIME = CTIME
      
      IF (NSTEP .GT. NSTEPS) THEN
         IF (NCSTEP .EQ. 0) THEN
            WRITE (*, 10010) 'There are no time steps'
         ELSE
            WRITE (*, 10010)
     &         'All time steps have been read from the database'
         END IF
      END IF

  150 CONTINUE
      NSTEP = NCSTEP
      RETURN

  160 CONTINUE
      WRITE (ERRMSG, 10000) 'DATABASE HEADER'
      GOTO 180
  170 CONTINUE
      WRITE (ERRMSG, 10000) 'DATABASE TIME STEPS'
      GOTO 180
  180 CONTINUE
      CALL WDBERR (IERR, ERRMSG)
      NSTEP = NCSTEP
      RETURN

10000  FORMAT (A)
10010  FORMAT (1X, 5A)
      END
