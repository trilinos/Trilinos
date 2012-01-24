C Copyright(C) 2011 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C           
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C                         
C * Neither the name of Sandia Corporation nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
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
      SUBROUTINE NEWMAP (MAPEL, MAPEL3, IXEL, INCEL, NREL, IELCOL)
C=======================================================================

C   $Id: newmap.f,v 1.1 1990/08/20 12:22:22 gdsjaar Exp $
C   $Log: newmap.f,v $
C   Revision 1.1  1990/08/20 12:22:22  gdsjaar
C   Initial revision
C

C   --*** NEWMAP *** (GEN3D) Calculate 3D element order map
C   --   Written by Amy Gilkey - revised 05/05/86
C   --
C   --NEWMAP calculates the element order map for the 3D database.
C   --
C   --Parameters:
C   --   MAPEL - IN - the 2D element order map
C   --   MAPEL3 - OUT - the 3D element order map
C   --   IXEL - IN - the new index for each element
C   --   INCEL - IN - the increment for each element, needed for blocks
C   --      that become multiple blocks
C   --   NREL - IN - the number of new elements generated for each element
C   --   IELCOL - IN - the row number for each element, 0 if not needed
C   --
C   --Common Variables:
C   --   Uses NUMEL of /DBNUMS/
C   --   Uses NUMNP3 of /DBNUM3/
C   --   Uses NEREPL of /PARAMS/

      INCLUDE 'gs_dbnums.blk'
      INCLUDE 'gs_dbnum3.blk'
      INCLUDE 'gs_params.blk'

      INTEGER MAPEL(NUMEL), MAPEL3(NUMEL3)
      INTEGER IXEL(*), INCEL(*), NREL(*), IELCOL(*)

C   --Element order map - add on elements for each plate/slice

      N = 0
      DO 20 I = 1, NUMEL
         IEL = MAPEL(I)
         JEL = IXEL(IEL)
         DO 10 NR = 1, NREL(IEL)
            N = N + 1
            MAPEL3(N) = JEL
            JEL = JEL + INCEL(IEL)
   10    CONTINUE
   20 CONTINUE

      RETURN
      END
