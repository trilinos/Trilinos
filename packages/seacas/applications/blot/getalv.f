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

C $Log: getalv.f,v $
C Revision 1.4  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.3  2009/01/22 21:34:21  gdsjaar
C There were several inline dbnums common blocks. Replaced with the
C include so they all have the same size with the added variable types.
C
C Added minor support for nodeset and sideset variables.
C
C It can print the count and the names, but that is all
C at this time.
C
C Revision 1.2  2007/11/14 20:14:53  gdsjaar
C Added optional 'alive value' to the death on variable command.  The
C default value is 0.0, but you can now specify a different value to
C indicate aliveness (for example, the presto DEATH_DUMMY_VAR treats 1.0
C as the alive value).
C
C Example: DEATH ON DEATH_DUMMY_VAR 1
C
C Revision 1.1  1994/04/07 20:01:31  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
CRevision 1.2  1990/12/14  08:51:02  gdsjaar
CAdded RCS Id and Log to all files
C
C=======================================================================
      SUBROUTINE GETALV (A, NALVAR, ALIVAL, ISTEP, LENE, ISEVOK,
     *  ALIVE, VAR)
C=======================================================================

C   --*** GETALV *** (MESH) Read birth/death variable
C   --   Written by Amy Gilkey - revised 10/28/87
C   --
C   --GETALV reads the values for the requested birth/death variable and
C   --returns the element state.
C   --
C   --The element is alive iff the input variable is 0.0.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   NALVAR - IN - the variable sequence number
C   --   ALIVAL - IN - the value to indicate element is fully alive
C   --   ISTEP - IN - the time step number
C   --   LENE - IN - the cumulative element counts by element block
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   ALIVE - OUT - true iff the element i is alive
C   --   VAR - SCRATCH - the birth/death variable array; may be ALIVE
C   --
C   --Common Variables:
C   --   Uses NUMEL, NELBLK of /DBNUMS/

      include 'dbnums.blk'

      DIMENSION A(*)
      INTEGER LENE(0:*)
      LOGICAL ISEVOK(NELBLK,NVAREL)
      LOGICAL ALIVE(NUMEL)
      REAL VAR(NUMEL)

      CHARACTER CDUM

      CALL GETVAR (A, NALVAR, -1, ISTEP, NUMEL, VAR)

      CALL DBVTYP (NALVAR, CDUM, IDALV)
      DO 120 IELB = 1, NELBLK
         IF (ISEVOK(IELB,IDALV)) THEN
            DO 100 IEL = LENE(IELB-1)+1, LENE(IELB)
               ALIVE(IEL) = (VAR(IEL) .EQ. ALIVAL)
  100       CONTINUE
         ELSE
            DO 110 IEL = LENE(IELB-1)+1, LENE(IELB)
               ALIVE(IEL) = .TRUE.
  110       CONTINUE
         END IF
  120 CONTINUE

      RETURN
      END
