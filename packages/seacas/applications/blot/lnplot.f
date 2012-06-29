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
      SUBROUTINE LNPLOT (NPTIMS, NPTIMW, XLN, YLN, ZLN, IXNODE, *)
C=======================================================================

C   --*** LNPLOT *** (PATHLN) Plot the pathlines
C   --   Written by Amy Gilkey - revised 05/27/88
C   --
C   --LNPLOT plots the pathlines.  For 3D meshes, the pathlines are rotated.
C   --
C   --Parameters:
C   --   NPTIMS - IN - the number of points on a history pathline
C   --   NPTIMW - IN - the number of points on a non-history pathline
C   --   XLN, YLN, ZLN - IN - the pathline data
C   --   IXNODE - SCRATCH - size = NPTIMS (3D only)
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NLNCRV, ILVID of /LNVARS/
C   --   Uses ROTMAT, ROTCEN of /ROTOPT/

      PARAMETER (NUMSYM = 6, NUMLIN = 6)

      include 'lnvars.blk'
      include 'd3nums.blk'
      include 'rotopt.blk'

      REAL XLN(NPTIMS,NLNCRV), YLN(NPTIMS,NLNCRV), ZLN(NPTIMS,NLNCRV)
      INTEGER IXNODE(NPTIMS)

      LOGICAL GRABRT
      LOGICAL NUMCRV
      CHARACTER TYP
      CHARACTER*8 LABSID

      lintyp = 1
      isytyp = 0
      numcrv = .false.
      labsid = 'NONE'
      
      IF (IS3DIM) THEN
         DO 100 I = 1, NPTIMS
            IXNODE(I) = I
  100    CONTINUE
      END IF

      DO 110 NP = 1, NLNCRV

         IF (GRABRT()) RETURN 1
         CALL GRCOLR (NP)
         CALL GRSYMB (LINTYP, ISYTYP, NP)

         CALL DBVTYP (ILVID(1,NP), TYP, IDUM)
         IF (TYP .EQ. 'H') THEN
            NPTS = NPTIMS
         ELSE
            NPTS = NPTIMW
         END IF

C      --Rotate the 3D pathline

         IF (IS3DIM) THEN
            CALL ROTATE (NPTS, IXNODE, ROTMAT, ROTCEN,
     &         XLN, YLN, ZLN, XLN, YLN, ZLN)
         END IF

C      --Plot pathline

         IF (GRABRT()) RETURN 1
         CALL PLTCUR (XLN(1,NP), YLN(1,NP), NPTS)

         IF (NUMCRV) THEN
            IF (GRABRT()) RETURN 1
            CALL GRNCRV (LABSID, NP, NPTS,
     &         XLN(1,NP), YLN(1,NP), (LINTYP .EQ. 0))
         END IF

C      --Finish plot

         CALL PLTFLU
  110 CONTINUE

      RETURN
      END
