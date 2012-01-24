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
      SUBROUTINE PROXYZ (XN, YN, XN3, YN3, ZN3, ATRIB)
C=======================================================================

C   $Id: proxyz.f,v 1.2 1991/01/09 12:59:30 gdsjaar Exp $
C   $Log: proxyz.f,v $
C   Revision 1.2  1991/01/09 12:59:30  gdsjaar
C   Initial conversion from GEN3D to GENSHELL, no BC yet
C
c Revision 1.1.1.1  90/08/20  12:22:35  gdsjaar
c Gen3D Mesh Generation Program
c 
c Revision 1.1  90/08/20  12:22:34  gdsjaar
c Initial revision
c 

C   --*** PROXYZ *** (GEN3D) Calculate 3D coordinates for experimental
C   --   Modified by Greg Sjaardema - 02/06/89
C   --
C   --Parameters:
C   --   XN, YN - IN - the 2D coordinates, destroyed
C   --   XN3, YN3, ZN3 - OUT - the 3D coordinates
C   --   IXNP - IN - the new index for each node
C   --   NRNP - IN - the number of new nodes generated for each node
C   --   ZCORD - SCRATCH - size = NNREPL, holds z coordinate for transformations
C   --   SINANG, COSANG - SCRATCH - size = NNREPL, holds sin and cos of
C   --      angles for rotations
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP of /DBNUMS/
C   --   Uses NDIM3, NUMNP3 of /DBNUM3/
C   --   Uses DOTRAN, NNREPL, DIM3, NRTRAN, D3TRAN, XXGRAD,
C   --      CENTER, NUMCOL, NUMROW of /PARAMS/
C   --   Uses XOFFS, YOFFS, ZOFFS of /XYZOFF/
C   --   Uses ROT3D, ROTMAT of /XYZROT/

      INCLUDE 'gs_dbnums.blk'
      INCLUDE 'gs_dbnum3.blk'
      INCLUDE 'gs_params.blk'
      INCLUDE 'gs_xxxxx.blk'

      REAL XN(NUMNP), YN(NUMNP),
     &   XN3(NUMNP3), YN3(NUMNP3), ZN3(NUMNP3)
      REAL ATRIB(NUMEL)

C      --Initialize the parametric interval distance

      ZTOT = 0.0
      DO 10 IBLK = 1, NBLK
         IF (NRTRAN(IBLK) .GT. 0) THEN
            ZTOT = ZTOT + D3TRAN(IBLK)
         ELSE
            CALL PRTERR ('PROGRAM',
     *         'Zero translations found')
            NBLK = IBLK
            GO TO 20
         END IF
   10 CONTINUE

      IF (ZTOT .EQ. 0.0) THEN
         CALL PRTERR ('CMDERR', 'Total translation distance is zero')
         STOP
      END IF

      NXTNR = 1
      ZEND = 0.0
   20 CONTINUE
      DO 30 IBLK = 1, NBLK
         ZBEG = ZEND
         ZEND = ZBEG + D3TRAN(IBLK)
         CALL INIGRD (ZBEG/ZTOT, ZEND/ZTOT, ZGRAD(IBLK),
     *      NRTRAN(IBLK), NRTRAN(IBLK)+1, ZCORD(NXTNR) )
         NXTNR = NXTNR + NRTRAN(IBLK)
   30 CONTINUE

C      --Project bottom surface onto a plane

      IF (.NOT. ISXWRP) THEN
         DO 50 INP = 1, NUMNP
            XB = XN(INP)
            YB = YN(INP)
            ZB = 0.0
            XT = (XN(INP) - XXSCL0) * XXSCAL + XXSCL0
            YT = (YN(INP) - XYSCL0) * XYSCAL + XYSCL0
            ZT =  - ZTOT + (XXA * XT + XXB * YT) / XXC
            XT =  XT + XXOFFS
            YT =  YT + XYOFFS

            JNP = IXNP(INP) - 1
            DO 40 NR = 1, NNREPL
               XN3(JNP+NR) = XB + (XT - XB) * ZCORD(NR)
               YN3(JNP+NR) = YB + (YT - YB) * ZCORD(NR)
               ZN3(JNP+NR) = ZB + (ZT - ZB) * ZCORD(NR)
   40       CONTINUE
   50    CONTINUE

C ... Warp of translated surface

      ELSE
         IF (CONVEX) THEN
            ZCEN = ZTOT - XWARP
            RMULT = 1.0
         ELSE
            ZCEN = ZTOT + XWARP
            RMULT = -1.0
         END IF
         DO 70 INP = 1, NUMNP
            XB = XN(INP)
            YB = YN(INP)
            ZB = 0.0
            XT = (XN(INP) - XXSCL0) * XXSCAL + XXSCL0
            YT = (YN(INP) - XYSCL0) * XYSCAL + XYSCL0

C ... Note: ZT must be negative

            ZT =  - (ZCEN + RMULT * SQRT(XWARP**2 - XB**2 - YB**2))
            XT =  XT + XXOFFS
            YT =  YT + XYOFFS

            JNP = IXNP(INP) - 1
            DO 60 NR = 1, NNREPL
               XN3(JNP+NR) = XB + (XT - XB) * ZCORD(NR)
               YN3(JNP+NR) = YB + (YT - YB) * ZCORD(NR)
               ZN3(JNP+NR) = ZB + (ZT - ZB) * ZCORD(NR)
   60       CONTINUE
   70    CONTINUE
      END IF

      RETURN
      END
