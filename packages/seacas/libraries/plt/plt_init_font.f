C    Copyright (C) 2009-2017 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C
C        * Neither the name of NTESS nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
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

      SUBROUTINE PLT_SANSERIF()
      DIMENSION FNTDAT(4,4501)

      call initserif(fntdat)

      CHSP = 4.0
      POINT = 21.0
      call plt_init_font(chsp, point, fntdat)
      return
      end

      SUBROUTINE PLT_STICK()
      DIMENSION FNTDAT(4,1936)

      call initstick(fntdat)

      CHSP = 2.0
      POINT = 7.0
      call plt_init_font(chsp, point, fntdat)
      return
      end

      SUBROUTINE PLT_ROMAN()
      DIMENSION FNTDAT(4,4122)

      call initroman(fntdat);

      CHSP = 4.0
      POINT = 21.0
      call plt_init_font(chsp, point, fntdat)
      return
      end

      subroutine plt_init_font(chsp, point, fntdat)
      real fntdat(4,*)

      REAL TEXTP(40)
      COMMON /TEXT/TEXTP

      INTEGER IDEX(200,2)
      INTEGER NVECT(200,2)
      REAL XSIZE(200,2)
      REAL YSIZE(200,2)
      REAL X0(2300,2)
      REAL Y0(2300,2)
      REAL X1(2300,2)
      REAL Y1(2300,2)
      COMMON /FONT/IDEX,NVECT,XSIZE,YSIZE,X0,Y0,X1,Y1

      logical first, last
      data first /.false./

      DIMENSION RVAL(4)

      TEXTP(39) = CHSP
      TEXTP(40) = POINT

      I = 0
      DO 2060 IC = 1, 2
        LAST = .FALSE.
        NTOTVC = 0
        FIRST = .TRUE.
        NLOCVC = 0
        XMAX = 0.
        YMAX = 0.
        MAXVECT = 0

 2080   CONTINUE
        i = i + 1
        rval(1) = fntdat(1,i)
        rval(2) = fntdat(2,i)
        rval(3) = fntdat(3,i)
        rval(4) = fntdat(4,i)

        IF (RVAL(1) .EQ. -99.) THEN
          IF (FIRST) THEN
            FIRST = .FALSE.
            ICH = RVAL(4)
            IDEX(ICH,ic) = 1
            GOTO 2090
          ELSE
            IF (RVAL(4) .EQ. 999.) THEN
              LAST = .TRUE.
              GOTO 2100
            ENDIF
            IF (ICH .GT. 32) THEN
              XSIZE(ICH,ic) = (XMAX + CHSP)/POINT
            ELSE
              XSIZE(ICH,ic) = (XMAX)/POINT
            ENDIF
            YSIZE(ICH,ic) = YMAX/POINT
            NVECT(ICH,ic) = NLOCVC
            IF (MAXVECT .LT. NLOCVC) THEN
              MAXVECT = NLOCVC
              MAXCHR = ICH
            ENDIF
            ICH = RVAL(4)
            IDEX(ICH,ic) = NTOTVC + 1
            XMAX = 0.
            YMAX = 0.
            NLOCVC = 0
          ENDIF
        ELSE
          NTOTVC = NTOTVC + 1
          NLOCVC = NLOCVC + 1
          X0(NTOTVC,ic) = RVAL(1)/POINT
          Y0(NTOTVC,ic) = RVAL(2)/POINT
          X1(NTOTVC,ic) = RVAL(3)/POINT
          Y1(NTOTVC,ic) = RVAL(4)/POINT
          IF (RVAL(1) .GT. XMAX) THEN
            XMAX = RVAL(1)
          ENDIF
          IF (RVAL(3) .GT. XMAX) THEN
            XMAX = RVAL(3)
          ENDIF
          IF (RVAL(2) .GT. YMAX) THEN
            YMAX = RVAL(2)
          ENDIF
          IF (RVAL(4) .GT. YMAX) THEN
            YMAX = RVAL(4)
          ENDIF
        ENDIF
        IF (LAST) THEN
          GOTO 2100
        ENDIF
 2090   GOTO 2080

 2100   CONTINUE
 100    CONTINUE
        IF (ICH .GT. 32) THEN
          XSIZE(ICH,ic) = (XMAX + CHSP)/POINT
        ELSE
          XSIZE(ICH,ic) = (XMAX)/POINT
        ENDIF
        YSIZE(ICH,ic) = YMAX/POINT
        NVECT(ICH,ic) = NLOCVC
        DO 2110 I = 48, 57
          XSIZE(I,ic) = 1.
 2110   CONTINUE
        XSIZE(43,ic) = 1.
        XSIZE(45,ic) = 1.
        NVECT(32,ic) = 0
 2060 CONTINUE
      END

