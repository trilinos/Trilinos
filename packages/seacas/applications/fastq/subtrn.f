C    Copyright(C) 2014-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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
C

C $Id: subtrn.f,v 1.1 1990/11/30 11:16:58 gdsjaar Exp $
C $Log: subtrn.f,v $
C Revision 1.1  1990/11/30 11:16:58  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]SUBTRN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SUBTRN (NPER, NEWPER, IP, X, Y, NID, XSUB, YSUB,
     &   NIDSUB, I1, I2, I3, I4, I5, I6, I7, I8, XCEN1, YCEN1, XCEN2,
     &   YCEN2, XMID1, YMID1, XMID2, YMID2, CCW, ERR)
C***********************************************************************
C
C  SUBROUTINE SUBTRN = PUTS A TRANSITION'S SUBREGION'S PERIMETER INTO
C                      THE NPERIM ARRAYS
C
C***********************************************************************
C
      DIMENSION X (NPER), Y (NPER), NID (NPER)
      DIMENSION XSUB (NPER), YSUB (NPER), NIDSUB (NPER)
C
C  PUT THE CORRECT PORTION OF THE PERIMETER IN XSUB,  YSUB,  AND NIDSUB
C  BASED ON THE VALUE OF IP (WHICH OF THE 6 SUBREGIONS ARE NEEDED)
C
      KOUNT = 0
C
C  SUBREGION 1  -  SIDE 1
C
      IF (IP .EQ. 1) THEN
         DO 100 I = I1, I2
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = X (I)
            YSUB (KOUNT) = Y (I)
            NIDSUB (KOUNT) = NID (I)
  100    CONTINUE
C
C  SUBREGION 1  -  SIDE 2
C
         XDIF = XCEN2 - X (I2)
         YDIF = YCEN2 - Y (I2)
         XINT = XDIF / FLOAT (NPER - I8 + 1)
         YINT = YDIF / FLOAT (NPER - I8 + 1)
         DO 110 I = 1, NPER - I8
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
            YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
            NIDSUB (KOUNT) = 300000 + NPER - I8 - I + 2
  110    CONTINUE
         KOUNT = KOUNT + 1
         XSUB (KOUNT) = XCEN2
         YSUB (KOUNT) = YCEN2
         NIDSUB (KOUNT) = 200000
C
C  SUBREGION 1  -  SIDE 3
C
         XDIF = X (I8) - XCEN2
         YDIF = Y (I8) - YCEN2
         XINT = XDIF / FLOAT (I2 - I1)
         YINT = YDIF / FLOAT (I2 - I1)
         DO 120 I = 1, I2 - I1 - 1
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
            YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
            NIDSUB (KOUNT) = 100000 + I + 1
  120    CONTINUE
C
C  SUBREGION 1  -  SIDE 4
C
         DO 130 I = I8, NPER
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = X (I)
            YSUB (KOUNT) = Y (I)
            NIDSUB (KOUNT) = NID (I)
  130    CONTINUE
C
C  SUBREGION 2  -  SIDE 1
C
      ELSEIF (IP .EQ. 2) THEN
         DO 140 I = I7, I8
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = X (I)
            YSUB (KOUNT) = Y (I)
            NIDSUB (KOUNT) = NID (I)
  140    CONTINUE
C
C  SUBREGION 2  -  SIDE 2
C
         XDIF = XCEN2 - X (I8)
         YDIF = YCEN2 - Y (I8)
         XINT = XDIF / FLOAT (I2 - I1)
         YINT = YDIF / FLOAT (I2 - I1)
         DO 150 I = 1, I2 - I1 - 1
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
            YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
            NIDSUB (KOUNT) = 100000 + I2 - I1 - I + 1
  150    CONTINUE
         KOUNT = KOUNT + 1
         XSUB (KOUNT) = XCEN2
         YSUB (KOUNT) = YCEN2
         NIDSUB (KOUNT) = 200000
C
C  SUBREGION 2  -  SIDE 3
C
         XDIF = XMID2 - XCEN2
         YDIF = YMID2 - YCEN2
         XINT = XDIF / FLOAT (I3 - I2)
         YINT = YDIF / FLOAT (I3 - I2)
         DO 160 I = 1, I3 - I2 - 1
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
            YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
            NIDSUB (KOUNT) = 200000 + I + 1
  160    CONTINUE
         KOUNT = KOUNT + 1
         XSUB (KOUNT) = XMID2
         YSUB (KOUNT) = YMID2
         NIDSUB (KOUNT) = 700000 + NPER - I8 + 2
C
C  SUBREGION 2  -  SIDE 4
C
         XDIF = X (I7) - XMID2
         YDIF = Y (I7) - YMID2
         XINT = XDIF / FLOAT (I2 - I1)
         YINT = YDIF / FLOAT (I2 - I1)
         DO 170 I = 1, I2 - I1 - 1
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
            YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
            NIDSUB (KOUNT) = NIDSUB (KOUNT - 1) + 1
  170    CONTINUE
C
C  SUBREGION 3  -  SIDE 1
C
      ELSEIF (IP .EQ. 3) THEN
         KOUNT = KOUNT + 1
         XSUB (KOUNT) = X (I3)
         YSUB (KOUNT) = Y (I3)
         NIDSUB (KOUNT) = NID (I3)
         XDIF = XMID2 - X (I3)
         YDIF = YMID2 - Y (I3)
         XINT = XDIF / FLOAT (NPER - I8 + 1)
         YINT = YDIF / FLOAT (NPER - I8 + 1)
         DO 180 I = 2, NPER - I8 + 1
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
            YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
            NIDSUB (KOUNT) = 700000 + I
  180    CONTINUE
         KOUNT = KOUNT + 1
         XSUB (KOUNT) = XMID2
         YSUB (KOUNT) = YMID2
         NIDSUB (KOUNT) = 700000 + NPER - I8 + 2
C
C  SUBREGION 3  -  SIDE 2
C
         XDIF = XCEN2 - XMID2
         YDIF = YCEN2 - YMID2
         XINT = XDIF / FLOAT (I3 - I2)
         YINT = YDIF / FLOAT (I3 - I2)
         DO 190 I = 1, I3 - I2 - 1
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
            YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
            NIDSUB (KOUNT) = 200000 + I3 - I2 - I + 1
  190    CONTINUE
         KOUNT = KOUNT + 1
         XSUB (KOUNT) = XCEN2
         YSUB (KOUNT) = YCEN2
         NIDSUB (KOUNT) = 200000
C
C  SUBREGION 3  -  SIDE 3
C
         XDIF = X (I2) - XCEN2
         YDIF = Y (I2) - YCEN2
         XINT = XDIF / FLOAT (NPER - I8 + 1)
         YINT = YDIF / FLOAT (NPER - I8 + 1)
         DO 200 I = 1, NPER - I8
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
            YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
            NIDSUB (KOUNT) = 300000 + I + 1
  200    CONTINUE
C
C  SUBREGION 3  -  SIDE 4
C
         DO 210 I = I2, I3 - 1
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = X (I)
            YSUB (KOUNT) = Y (I)
            NIDSUB (KOUNT) = NID (I)
  210    CONTINUE
C
C  SUBREGION 4  -  SIDE 1 AND 2
C
      ELSEIF (IP .EQ. 4) THEN
         DO 220 I = I4, I6
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = X (I)
            YSUB (KOUNT) = Y (I)
            NIDSUB (KOUNT) = NID (I)
  220    CONTINUE
C
C  SUBREGION 4  -  SIDE 3
C
         XDIF = XCEN1 - X (I6)
         YDIF = YCEN1 - Y (I6)
         XINT = XDIF / FLOAT (I5 - I4)
         YINT = YDIF / FLOAT (I5 - I4)
         DO 230 I = 1, I5 - I4 - 1
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
            YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
            NIDSUB (KOUNT) = 400000 + I5 - I4 - I + 1
  230    CONTINUE
         KOUNT = KOUNT + 1
         XSUB (KOUNT) = XCEN1
         YSUB (KOUNT) = YCEN1
         NIDSUB (KOUNT) = 100000
C
C  SUBREGION 4  -  SIDE 4
C
         XDIF = X (I4) - XCEN1
         YDIF = Y (I4) - YCEN1
         XINT = XDIF / FLOAT (I6 - I5)
         YINT = YDIF / FLOAT (I6 - I5)
         DO 240 I = 1, I6 - I5 - 1
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
            YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
            NIDSUB (KOUNT) = 600000 + I + 1
  240    CONTINUE
C
C  SUBREGION 5  -  SIDE 1
C
      ELSEIF (IP .EQ. 5) THEN
         DO 250 I = I6, I7
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = X (I)
            YSUB (KOUNT) = Y (I)
            NIDSUB (KOUNT) = NID (I)
  250    CONTINUE
C
C  SUBREGION 5  -  SIDE 2
C
         XDIF = XMID1 - X (I7)
         YDIF = YMID1 - Y (I7)
         XINT = XDIF / FLOAT (I5 - I4)
         YINT = YDIF / FLOAT (I5 - I4)
         DO 260 I = 1, I5 - I4 - 1
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
            YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
            NIDSUB (KOUNT) = 700000 + I6 - I4 - I + 1
  260    CONTINUE
         KOUNT = KOUNT + 1
         XSUB (KOUNT) = XMID1
         YSUB (KOUNT) = YMID1
         NIDSUB (KOUNT) = 700000 + I6 - I5 + 1
C
C  SUBREGION 5  -  SIDE 3
C
         XDIF = XCEN1 - XMID1
         YDIF = YCEN1 - YMID1
         XINT = XDIF / FLOAT (I4 - I3)
         YINT = YDIF / FLOAT (I4 - I3)
         DO 270 I = 1, I4 - I3 - 1
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
            YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
            NIDSUB (KOUNT) = 500000 + I4 - I3 - I + 1
  270    CONTINUE
         KOUNT = KOUNT + 1
         XSUB (KOUNT) = XCEN1
         YSUB (KOUNT) = YCEN1
         NIDSUB (KOUNT) = 100000
C
C  SUBREGION 5  -  SIDE 4
C
         XDIF = X (I6) - XCEN1
         YDIF = Y (I6) - YCEN1
         XINT = XDIF / FLOAT (I5 - I4)
         YINT = YDIF / FLOAT (I5 - I4)
         DO 280 I = 1, I5 - I4 - 1
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
            YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
            NIDSUB (KOUNT) = 400000 + I + 1
  280    CONTINUE
C
C  SUBREGION 6  -  SIDE 1
C
      ELSEIF (IP .EQ. 6) THEN
         KOUNT = KOUNT + 1
         XSUB (KOUNT) = XMID1
         YSUB (KOUNT) = YMID1
         NIDSUB (KOUNT) = 700000 + I6 - I5 + 1
         XDIF = X (I3) - XMID1
         YDIF = Y (I3) - YMID1
         XINT = XDIF / FLOAT (I6 - I5)
         YINT = YDIF / FLOAT (I6 - I5)
         DO 290 I = 1, I6 - I5 - 1
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
            YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
            NIDSUB (KOUNT) = 700000 + I6 - I5 - I + 1
  290    CONTINUE
C
C  SUBREGION 6  -  SIDE 2
C
         DO 300 I = I3, I4
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = X (I)
            YSUB (KOUNT) = Y (I)
            NIDSUB (KOUNT) = NID (I)
  300    CONTINUE
C
C  SUBREGION 6  -  SIDE 3
C
         XDIF = XCEN1 - X (I4)
         YDIF = YCEN1 - Y (I4)
         XINT = XDIF / FLOAT (I6 - I5)
         YINT = YDIF / FLOAT (I6 - I5)
         DO 310 I = 1, I6 - I5 - 1
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
            YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
            NIDSUB (KOUNT) = 600000 + I6 - I5 - I + 1
  310    CONTINUE
         KOUNT = KOUNT + 1
         XSUB (KOUNT) = XCEN1
         YSUB (KOUNT) = YCEN1
         NIDSUB (KOUNT) = 100000
C
C  SUBREGION 6  -  SIDE 4
C
         XDIF = XMID1 - XCEN1
         YDIF = YMID1 - YCEN1
         XINT = XDIF / FLOAT (I4 - I3)
         YINT = YDIF / FLOAT (I4 - I3)
         DO 320 I = 1, I4 - I3 - 1
            KOUNT = KOUNT + 1
            XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
            YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
            NIDSUB (KOUNT) = 500000 + I + 1
  320    CONTINUE
      ENDIF
C
      NEWPER = KOUNT
      RETURN
C
      END
