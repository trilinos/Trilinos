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
      SUBROUTINE RDSTEP (NDB, NCSTEP,
     &     NUMELB, IDELB, ISEVOK,
     &     TIME, VARGL, VARNP, VAREL,
     $     neldm)
C=======================================================================

C   --*** RDSTEP *** (GROPE) Read current database variables
C   --
C   --RDSTEP reads the variable data for the next time step from the
C   --database.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   NCSTEP - IN/OUT - the current step number
C   --   NVARGL - IN - the number of global variables
C   --   NVARNP - IN - the number of nodal variables
C   --   NVAREL - IN - the number of element variables
C   --   NUMELB - IN - the number of elements per block
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(i,j) is NOT 0
C   --   TIME - OUT - the time step time
C   --   VARGL - OUT - the global variables for the time step
C   --   VARNP - OUT - the nodal variables for the time step
C   --   VAREL - OUT - the element variables for the time step

      include 'exodusII.inc'
      include 'dbnums.blk'

      INTEGER NUMELB(*)
      INTEGER IDELB(*)
      INTEGER ISEVOK(nvarel,*)
      REAL VARGL(*)
      REAL VARNP(numnp, *)
      REAL VAREL(numel, *)

      CHARACTER*80 ERRMSG

C ... Initialize variables to known value in case error on read
      TIME = -999.9
      CALL INIREA (NVARGL, -999.9, VARGL)
      CALL INIREA (NVARNP * NUMNP, -999.9, VARNP)
      CALL INIREA (NVAREL * NUMEL, -999.9, VAREL)

C ... Read the time for the specified step
      call exgtim(ndb, ncstep, time, ierr)
      if (ierr .ne. 0) go to 150

C ... Read the global variables (if any)
      if (nvargl .gt. 0) then
         call exggv (ndb, ncstep, nvargl, vargl, ierr)
         if (ierr .ne. 0) go to 170
      end if
      
C ... Read the nodal variables (if any)
      DO 120 I = 1, NVARNP
        call exgnv (ndb, ncstep, i, numnp, varnp(1,i), ierr)
        if (ierr .ne. 0) go to 180
 120  CONTINUE
      
C ... Read the element variables (if any)
      DO 140 I = 1, NVAREL
        IEL0 = 1
        DO 130 IELB = 1, NELBLK
          IF (ISEVOK(I,IELB) .NE. 0 .AND. numelb(ielb) .GT. 0) THEN
             call exgev (ndb, ncstep, i, idelb(ielb),
     &            numelb(ielb), varel(iel0,i), ierr)
             if (ierr .ne. 0) go to 190
          END IF
          IEL0 = IEL0 + NUMELB(IELB)
 130    CONTINUE
 140  CONTINUE
      
      RETURN

  150 CONTINUE
      WRITE (ERRMSG, 10000, IOSTAT=IDUM)
     &   'TIME on step', NCSTEP
      GOTO 200
  170 CONTINUE
      WRITE (ERRMSG, 10000, IOSTAT=IDUM)
     &   'GLOBAL VARIABLES on step', NCSTEP
      GOTO 200
  180 CONTINUE
      WRITE (ERRMSG, 10000, IOSTAT=IDUM)
     &   'NODAL VARIABLE', I, ' on step', NCSTEP
      GOTO 200
  190 CONTINUE
      WRITE (ERRMSG, 10000, IOSTAT=IDUM)
     &   'ELEMENT VARIABLE', I, ' in block', IELB, ' on step', NCSTEP
      GOTO 200
  200 CONTINUE
      CALL WDBERR (IERR, ERRMSG)
  220 CONTINUE
      RETURN

10000  FORMAT (5 (A, I10))
      END
