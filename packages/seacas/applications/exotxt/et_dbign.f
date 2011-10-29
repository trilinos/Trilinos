C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
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
C 

C=======================================================================
      SUBROUTINE DBIGN (NDB, NUMESS, IDESS, NNESS, IXNESS,
     &                  LTNESS, LTNNN, IOERR)
C=======================================================================
C$Id: dbign.f,v 1.4 2007/10/17 18:46:09 gdsjaar Exp $
C$Log: dbign.f,v $
CRevision 1.4  2007/10/17 18:46:09  gdsjaar
CAdded copyright notice to all files.
C
Cexotxt2 is licensed under the BSD license
C
CRevision 1.3  1996/05/21 16:52:17  caforsy
CAdded read/write for property data.  Cleaned up exodusII error checks
C
CRevision 1.2  1995/11/07 15:01:25  gdsjaar
CInitial checkin of ACCESS/translate/exotxt2
C

C   --*** DBIGN *** Get node from the side sets
C   --   Written 9/10/95 for ExodusIIv2
C   --
C   --Parameters:
C   --   NDB    - IN  - the database number
C   --   NUMESS - IN  - number of side sets
C   --   IDESS  - OUT - array of side set IDS
C   --   NNESS  - OUT - array of the number of nodes for each side set
C   --   IXNESS - OUT - array of indices into LTNESS - 1st node each set 
C   --   LTNESS - OUT - array of nodes for all side sets
C   --   LTNNN  - OUT _ array of number of nodes for each side in a side sets
C   --   IOERR  - OUT - error flag

      INTEGER NDB
      INTEGER NUMESS
      INTEGER IDESS(*)
      INTEGER NNESS(*)
      INTEGER IXNESS(*)
      INTEGER LTNESS(*)
      INTEGER LTNNN(*)
      INTEGER IOERR
      IOERR = 0

C     Offset into element list of current side set
      ISOFF  = 0
C     Node count for current side set
      NODCNT = 0
      DO 100 I = 1, NUMESS
C        Set index of the first node for each side set in LTNESS
         IXNESS(I) = NODCNT + 1
         CALL EXGSP(NDB, IDESS(I), NSIDE, NDIST, IOERR)
C        NSIDE - number of sides in side set IDESS(I)
C        NDIST - number of distribution factors in side set IDESS(I)
         IF (IOERR .EQ. 1) RETURN

         CALL exgssn(NDB, IDESS(I), LTNNN(ISOFF+1),
     &               LTNESS(NODCNT+1), IOERR)
C        LTNNN(ISOFF+1) - number of nodes for each side in side set IDESS(I)
C        LTNESS(NODCNT+1) - nodes for current set
         IF (IOERR .EQ. 1) RETURN
C        Calculate node count sum for the current side set
         NCSUM = 0
         DO 90 J = 0, NSIDE-1
            NCSUM = NCSUM + LTNNN(ISOFF+1+J)
  90     CONTINUE
         NNESS(I) = NCSUM
         NODCNT = NODCNT + NCSUM
         ISOFF = ISOFF + NSIDE
 100  CONTINUE   

      RETURN
      END
