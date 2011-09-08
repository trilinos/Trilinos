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
      SUBROUTINE DBOSTE (NDB, ISTEP,
     &   NVARHI, NVARGL, NVARNP, NUMNP, NVAREL, NELBLK, NUMELB, ISEVOK,
     &   TIME, WHOTIM, VARHI, VARGL, VARNP, VAREL)
C=======================================================================
C$Id: dboste.f,v 1.3 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dboste.f,v $
CRevision 1.3  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.2  1990/10/02 12:21:48  gdsjaar
CFixed problem with MAX() in dimension statement.  Previous fix did not
Cwork, fixed by adding DBOST1 subroutine which is passed the MAX()'d
Cdimensions. 
C
c Revision 1.1.1.1  90/08/14  16:13:40  gdsjaar
c Testing
c 
c Revision 1.1  90/08/14  16:13:38  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:17  gdsjaar
c Initial revision
c 

C   --*** DBOSTE *** (EXOLIB) Write database variables for one time step
C   --   Written by Amy Gilkey - revised 10/14/87
C   --   Modified by Greg Sjaardema - 10/2/90
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
C   --      variable i of block j exists iff ISEVOK(j,i)
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
      LOGICAL ISEVOK(*)
      REAL TIME
      LOGICAL WHOTIM
      REAL VARHI(*)
      REAL VARGL(*)
      REAL VARNP(*)
      REAL VAREL(*)

      CALL DBOST1 (NDB, ISTEP,
     $   NVARHI, NVARGL, NVARNP, NUMNP, NVAREL, NELBLK, NUMELB, ISEVOK, 
     &   TIME, WHOTIM, VARHI, VARGL, VARNP, VAREL,
     $   MAX(1, NVARNP), MAX(1, NVAREL))

      RETURN
      END

C=======================================================================
      SUBROUTINE DBOST1 (NDB, ISTEP,
     &   NVARHI, NVARGL, NVARNP, NUMNP, NVAREL, NELBLK, NUMELB, ISEVOK,
     &   TIME, WHOTIM, VARHI, VARGL, VARNP, VAREL,
     &   NVNPDM, NVELDM)
C=======================================================================

      INTEGER NDB
      INTEGER ISTEP
      INTEGER NVARHI, NVARGL, NVARNP, NUMNP, NVAREL, NELBLK
      INTEGER NUMELB(NELBLK)
      LOGICAL ISEVOK(NELBLK,*)
      REAL TIME
      LOGICAL WHOTIM
      REAL VARHI(*)
      REAL VARGL(*)
      REAL VARNP(NVNPDM,*)
      REAL VAREL(NVELDM,*)

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

         DO 100 I = 1, NVARNP
            IF (NUMNP .GT. 0) THEN
               WRITE (NDB) (VARNP(I,INP), INP=1,NUMNP)
            ELSE
               WRITE (NDB) 0
            END IF
  100    CONTINUE

C      --Write element variables

         IEL0 = 0
         DO 120 IELB = 1, NELBLK
            DO 110 I = 1, NVAREL
               IF (ISEVOK(IELB,I)) THEN
                  IF (NUMELB(IELB) .GT. 0) THEN
                     WRITE (NDB) (VAREL(I,IEL0+N), N=1,NUMELB(IELB))
                  ELSE
                     WRITE (NDB) 0
                  END IF
               END IF
  110       CONTINUE
            IEL0 = IEL0 + NUMELB(IELB)
  120    CONTINUE
      END IF

      RETURN
      END
