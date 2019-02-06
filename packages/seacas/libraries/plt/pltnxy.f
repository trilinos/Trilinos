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
      SUBROUTINE PLTNXY(X,Y,NUM,XLAB,XUNIT,YLAB,YUNIT)
      REAL DEVCAP(23)
      REAL DEFOUT(7)
      COMMON /STATUS/DEVCAP,DEFOUT
      REAL DEVP(5)
      COMMON /DEVICE/DEVP
      REAL COLP(3)
      REAL PALETT(3,16)
      COMMON /COLOR/COLP,PALETT
      REAL TEXTP(40)
      COMMON /TEXT/TEXTP
      REAL VECTP(5)
      REAL XCUR
      REAL YCUR
      COMMON /VECTRC/VECTP,XCUR,YCUR
      INTEGER IDEX(200,2)
      INTEGER NVECT(200,2)
      REAL XSIZE(200,2)
      REAL YSIZE(200,2)
      REAL X0(2300,2)
      REAL Y0(2300,2)
      REAL X1(2300,2)
      REAL Y1(2300,2)
      COMMON /FONT/IDEX,NVECT,XSIZE,YSIZE,X0,Y0,X1,Y1
      REAL GRAPHP(100)
      COMMON /GRAPH/GRAPHP
      COMMON /MAPPAR/MAPP(11)
      REAL MAPP
      COMMON /STORAG/MEMORY(1000)
      CHARACTER*(*) XLAB,XUNIT,YLAB,YUNIT
      DIMENSION X(1),Y(1)
      REAL INTERX,INTERY

      XLENT = GRAPHP(3)
      YLENT = GRAPHP(4)
      CALL VECRGS(IABS(NUM),X,XMAX,XMIN)
      CALL VECRGS(IABS(NUM),Y,YMAX,YMIN)
      IF (GRAPHP(22).EQ.1.) THEN
         CALL PLTINO(XMIN,XMAX,FNLOWX,FNUPPX,INTERX,IEXPX,NMINX)
         XSTART = FNLOWX
         XEND = FNUPPX
         TNEXPX = 10.**IEXPX
         CALL PLTINO(YMIN,YMAX,FNLOWY,FNUPPY,INTERY,IEXPY,NMINY)
         YSTART = FNLOWY
         YEND = FNUPPY
         TNEXPY = 10.**IEXPY
         GRAPHP(24) = FNLOWX*TNEXPX
         GRAPHP(25) = FNUPPX*TNEXPX
         GRAPHP(26) = (FNUPPX-FNLOWX)/INTERX
         GRAPHP(27) = NMINX
         GRAPHP(28) = FNLOWY*TNEXPY
         GRAPHP(29) = FNUPPY*TNEXPY
         GRAPHP(30) = (FNUPPY-FNLOWY)/INTERY
         GRAPHP(31) = NMINY
         GRAPHP(78) = FNLOWX*TNEXPX
         GRAPHP(79) = GRAPHP(78)
         GRAPHP(80) = FNUPPX*TNEXPX
         GRAPHP(81) = INTERX*TNEXPX
         GRAPHP(82) = NMINX
         GRAPHP(83) = FNLOWY*TNEXPY
         GRAPHP(84) = GRAPHP(83)
         GRAPHP(85) = FNUPPY*TNEXPY
         GRAPHP(86) = INTERY*TNEXPY
         GRAPHP(87) = NMINY

      ELSE IF (GRAPHP(22).EQ.2.) THEN
         DELX = XMAX - XMIN
         DELY = YMAX - YMIN
         CALL PLTINO(XMIN,XMAX,FNLOWX,FNUPPX,INTERX,IEXPX,NMINX)
         CALL PLTINO(YMIN,YMAX,FNLOWY,FNUPPY,INTERY,IEXPY,NMINY)
         IF (DELX.GT.DELY) THEN
            IF (IEXPX.NE.IEXPY) THEN
               CALL PLTNCF(YMIN,'u',FNLOWY,IEXPX)
               IEXPY = IEXPX
               INTERY = INTERX
               IF (INTERY.EQ.2. .AND. AMOD(ABS(FNLOWY),2.).EQ.1.) THEN
                  FNLOWY = FNLOWY - 1.
               END IF

               NMINY = NMINX
            END IF

            FNUPPY = FNLOWY + FNUPPX - FNLOWX
            IF (FNUPPY*10.**IEXPY.LT.YMAX) THEN
               FNUPPY = FNUPPY + 1.
               FNUPPX = FNUPPX + 1.
            END IF

         END IF

         IF (DELX.LT.DELY) THEN
            IF (IEXPX.NE.IEXPY) THEN
               CALL PLTNCF(XMIN,'u',FNLOWX,IEXPY)
               IEXPX = IEXPY
               INTERX = INTERY
               IF (INTERX.EQ.2. .AND. AMOD(ABS(FNLOWX),2.).EQ.1.) THEN
                  FNLOWX = FNLOWX - 1.
               END IF

               NMINX = NMINY
            END IF

            FNUPPX = FNLOWX + FNUPPY - FNLOWY
            IF (FNUPPX*10.**IEXPX.LT.XMAX) THEN
               FNUPPX = FNUPPX + 1.
               FNUPPY = FNUPPY + 1.
            END IF

         END IF

         IF (GRAPHP(3).NE.GRAPHP(4)) THEN
            XLENT = AMIN1(GRAPHP(3),GRAPHP(4))
            YLENT = XLENT
         END IF

         TNEXPX = 10.**IEXPX
         TNEXPY = 10.**IEXPY
         XSTART = FNLOWX
         XEND = FNUPPX
         YSTART = FNLOWY
         YEND = FNUPPY
         GRAPHP(24) = FNLOWX*TNEXPX
         GRAPHP(25) = FNUPPX*TNEXPX
         GRAPHP(26) = (FNUPPX-FNLOWX)/INTERX
         GRAPHP(27) = NMINX
         GRAPHP(28) = FNLOWY*TNEXPY
         GRAPHP(29) = FNUPPY*TNEXPY
         GRAPHP(30) = (FNUPPY-FNLOWY)/INTERY
         GRAPHP(31) = NMINY
         GRAPHP(78) = FNLOWX*TNEXPX
         GRAPHP(79) = GRAPHP(78)
         GRAPHP(80) = FNUPPX*TNEXPX
         GRAPHP(81) = INTERX*TNEXPX
         GRAPHP(82) = NMINX
         GRAPHP(83) = FNLOWY*TNEXPY
         GRAPHP(84) = GRAPHP(83)
         GRAPHP(85) = FNUPPY*TNEXPY
         GRAPHP(86) = INTERY*TNEXPY
         GRAPHP(87) = NMINY

      ELSE IF (GRAPHP(22).EQ.3.) THEN
         TINT = (GRAPHP(25)-GRAPHP(24))/GRAPHP(26)
         IEXPX = NINT(LOG10(ABS(TINT)))
         TNEXPX = 10.**IEXPX
         FNLOWX = GRAPHP(24)/TNEXPX
         FNUPPX = GRAPHP(25)/TNEXPX
         INTERX = (FNUPPX-FNLOWX)/NINT(GRAPHP(26))
         NMINX = INT(GRAPHP(27))
         XSTART = FNLOWX
         XEND = FNUPPX
         TINT = (GRAPHP(29)-GRAPHP(28))/GRAPHP(30)
         IEXPY = NINT(LOG10(ABS(TINT)))
         TNEXPY = 10.**IEXPY
         FNLOWY = GRAPHP(28)/TNEXPY
         FNUPPY = GRAPHP(29)/TNEXPY
         INTERY = (FNUPPY-FNLOWY)/NINT(GRAPHP(30))
         NMINY = INT(GRAPHP(31))
         YSTART = FNLOWY
         YEND = FNUPPY
         GRAPHP(78) = XSTART*TNEXPX
         GRAPHP(79) = GRAPHP(78)
         GRAPHP(80) = XEND*TNEXPX
         GRAPHP(81) = INTERX*TNEXPX
         GRAPHP(82) = NMINX
         GRAPHP(83) = YSTART*TNEXPY
         GRAPHP(84) = GRAPHP(83)
         GRAPHP(85) = YEND*TNEXPY
         GRAPHP(86) = INTERY*TNEXPY
         GRAPHP(87) = NMINY

      ELSE IF (GRAPHP(22).EQ.4.) THEN
         IEXPX = NINT(LOG10(ABS(GRAPHP(81))))
         TNEXPX = 10.**IEXPX
         XSTART = GRAPHP(78)/TNEXPX
         XEND = GRAPHP(80)/TNEXPX
         FNLOWX = GRAPHP(79)/TNEXPX
         INTERX = GRAPHP(81)/TNEXPX
         NMINX = INT(GRAPHP(82))
         IEXPY = NINT(LOG10(ABS(GRAPHP(86))))
         TNEXPY = 10.**IEXPY
         YSTART = GRAPHP(83)/TNEXPY
         YEND = GRAPHP(85)/TNEXPY
         FNLOWY = GRAPHP(84)/TNEXPY
         INTERY = GRAPHP(86)/TNEXPY
         NMINY = INT(GRAPHP(87))
         GRAPHP(24) = XSTART*TNEXPX
         GRAPHP(25) = XEND*TNEXPX
         GRAPHP(26) = (XSTART-XEND)/INTERX
         GRAPHP(27) = NMINX
         GRAPHP(28) = YSTART*TNEXPY
         GRAPHP(29) = YEND*TNEXPY
         GRAPHP(30) = (YSTART-YEND)/INTERY
         GRAPHP(31) = NMINY
      END IF

      IF (GRAPHP(91).NE.-999999.) THEN
         FAC = 10.** (IEXPX-GRAPHP(91))
         IEXPX = INT(GRAPHP(91))
         TNEXPX = 10.**IEXPX
         XSTART = XSTART*FAC
         XEND = XEND*FAC
         FNLOWX = FNLOWX*FAC
         INTERX = INTERX*FAC
      END IF

      IF (GRAPHP(90).NE.-999999.) THEN
         FAC = 10.** (IEXPY-GRAPHP(90))
         IEXPY = INT(GRAPHP(90))
         TNEXPY = 10.**IEXPY
         YSTART = YSTART*FAC
         YEND = YEND*FAC
         FNLOWY = FNLOWY*FAC
         INTERY = INTERY*FAC
      END IF

      IF (GRAPHP(40).EQ.1.) THEN
         YSTART = YSTART*TNEXPY
         YEND = YEND*TNEXPY
         FNLOWY = FNLOWY*TNEXPY
         FNUPPY = FNUPPY*TNEXPY
         INTERY = INTERY*TNEXPY
         IEXPY = 0
         TNEXPY = 1.
      END IF

      IF (GRAPHP(40).EQ.2.) THEN
         XSTART = XSTART*TNEXPX
         XEND = XEND*TNEXPX
         FNLOWX = FNLOWX*TNEXPX
         FNUPPX = FNUPPX*TNEXPX
         INTERX = INTERX*TNEXPX
         IEXPX = 0
         TNEXPX = 1.
      END IF

      IF (GRAPHP(40).EQ.4.) THEN
         XSTART = XSTART*TNEXPX
         XEND = XEND*TNEXPX
         FNLOWX = FNLOWX*TNEXPX
         FNUPPX = FNUPPX*TNEXPX
         INTERX = INTERX*TNEXPX
         IEXPX = 0
         TNEXPX = 1.
         YSTART = YSTART*TNEXPY
         YEND = YEND*TNEXPY
         FNLOWY = FNLOWY*TNEXPY
         FNUPPY = FNUPPY*TNEXPY
         INTERY = INTERY*TNEXPY
         IEXPY = 0
         TNEXPY = 1.
      END IF

      CALL PLTGM2(XSTART*TNEXPX,XEND*TNEXPX,YSTART*TNEXPY,
     *            YEND*TNEXPY,GRAPHP(1),GRAPHP(1)+XLENT,GRAPHP(2),
     *            GRAPHP(2)+YLENT,GRAPHP(7))
      CALL PLTUWN(GRAPHP(7))
      CALL PLTAXS(GRAPHP(1),GRAPHP(2),XLENT,YLENT,'x',XSTART,XEND,
     *            FNLOWX,INT(GRAPHP(41)),INTERX,NMINX,XLAB,XUNIT,IEXPX)

      CALL PLTAXS(GRAPHP(1),GRAPHP(2),XLENT,YLENT,'y',YSTART,YEND,
     *            FNLOWY,INT(GRAPHP(42)),INTERY,NMINY,YLAB,YUNIT,IEXPY)

      CALL PLTCUR(X,Y,NUM)

      RETURN

      END
