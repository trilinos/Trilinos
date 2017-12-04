C Copyright (C) 2009-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C 
C     * Neither the name of NTESS nor the names of its
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
C 

C=======================================================================
      SUBROUTINE DBOSTE (NDB, ISTEP,
     &   NVARHI, NVARGL, NVARNP, NUMNP, NVAREL, NELBLK, NUMELB, ISEVOK,
     &   TIME, WHOTIM, VARHI, VARGL, VARNP, VAREL)
C=======================================================================

C   --*** DBOSTE *** (EXOLIB) Write database variables for one time step
C   --      Removed MAX from dimension statements, added routine
C   --      DBOST1 to do all work.
C   --
C   --DBOSTE writes the database history, global, nodal, and element variables
C   --for one time step.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   ISTEP - IN - the time step number
C   --   NVARHI - IN - the number of history variables
C   --   NVARGL - IN - the number of global variables
C   --   NVARNP - IN - the number of nodal variables
C   --   NUMNP - IN - the number of nodes
C   --   NVAREL - IN - the number of element variables
C   --   NELBLK - IN - the number of element blocks
C   --   NUMELB - IN - the number of elements per block
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(i,j)
C   --   TIME - IN - the time step time
C   --   WHOTIM - IN - true iff whole (versus history) time step
C   --   VARHI - IN - the history variables for the time step
C   --   VARGL - IN - the global variables for the time step
C   --   VARNP - IN - the nodal variables for the time step
C   --   VAREL - IN - the element variables for the time step
C   --
C   --Database must be positioned in front of time step upon entry;
C   --upon exit positioned after time step.

      INTEGER NDB
      INTEGER ISTEP
      INTEGER NVARHI, NVARGL, NVARNP, NUMNP, NVAREL, NELBLK
      INTEGER NUMELB(*)
c      LOGICAL ISEVOK(nelblk,*)
      integer ISEVOK(nvarel,*)
      REAL TIME
      LOGICAL WHOTIM
      REAL VARHI(*)
      REAL VARGL(*)
      REAL VARNP(*)
      REAL VAREL(*)

C   --Write step time

      IF (WHOTIM) THEN
         HISTFL = 0.0
      ELSE
         HISTFL = -1.0
      END IF

      WRITE (NDB) TIME, HISTFL

C   --Write history variables

      IF (NVARHI .GT. 0) THEN
         WRITE (NDB) (VARHI(I), I=1,NVARHI)
      ELSE
         WRITE (NDB) 0
      END IF

      IF (WHOTIM) THEN

C      --Write global variables

         IF (NVARGL .GT. 0) THEN
            WRITE (NDB) (VARGL(I), I=1,NVARGL)
         ELSE
            WRITE (NDB) 0
         END IF

C      --Write nodal variables

         indo=0
         DO 100 I = 1, NVARNP
            IF (NUMNP .GT. 0) THEN
               WRITE (NDB) (VARNP(indo+inp), INP=1,NUMNP)
               indo=indo+numnp
            ELSE
               WRITE (NDB) 0
            END IF
  100    CONTINUE

C      --Write element variables

         IELO = 0
         DO 120 IELB = 1, NELBLK
            DO 110 I = 1, NVAREL
               IF (ISEVOK(I,IELB) .ne. 0) THEN
                  IF (NUMELB(IELB) .GT. 0) THEN
                     WRITE (NDB) (VAREL(IELO+N), N=1,NUMELB(IELB))
                     ielo=ielo+numelb(ielb)
                  ELSE
                     WRITE (NDB) 0
                  END IF
               END IF
  110       CONTINUE
  120    CONTINUE
      END IF

      RETURN
      END
