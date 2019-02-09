C    Copyright(C) 1988-2017 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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

      SUBROUTINE CAVITY (A, CRD, IDESS, NEESS, NNESS, IPEESS, IPNESS,
     *   LTEESS, LTNESS, FACESS, DISP, NUMNP, NDIM, NUMESS,
     *   TIME, ITMSEL, TITLE, CENT, CENTER)
C
      include 'nu_io.blk'
      DIMENSION A(*), CRD(NUMNP,NDIM), IDESS(*), NEESS(*),
     *   NNESS(*), IPEESS(*), IPNESS(*), LTEESS(*), LTNESS(*),
     *   FACESS(*), TIME(*), DISP(NUMNP,NDIM), CENT(3)
      LOGICAL ITMSEL(*)
      CHARACTER*80 TITLE
      include 'nu_logs.blk'
      include 'nu_ptim.blk'
      include 'nu_cav.blk'
      LOGICAL ERROR, CENTER
C
      CALL GETCAV (ERROR, IDESS, NUMESS)
      IF (ERROR) RETURN
C
      TVOL = 0.0
      DO 10 NCAV = 1, NUMCAV
         IFLG = IFND(NCAV)
         IPTR = IPNESS(IFLG)
         IF (NDIM .EQ. 3) THEN
            CALL VOL3D( CRD, LTNESS(IPTR), NEESS(IFLG), VOLUME,
     *         NDIM, NUMESS, CENT, NUMNP, CENTER)
         ELSE
            CALL VOL2D( CRD, LTNESS(IPTR), NEESS(IFLG), VOLUME,
     *         NDIM, NUMESS, AXI, CENT, NUMNP, CENTER)
         END IF
C
         TVOL = TVOL + VOLUME
   10 CONTINUE
      DO 20 IO=IOMIN, IOMAX
         WRITE (IO,30) (ICAV(I),I=1,NUMCAV)
         IF (NDIM .EQ. 2) THEN
            WRITE (IO, 40) CENT(1),CENT(2)
         ELSE
            WRITE (IO, 50) CENT(1),CENT(2),CENT(3)
         END IF
         WRITE (IO,60) TVOL
   20 CONTINUE
   30 FORMAT (/' Cavity Flag(s): ',8I8)
   40 FORMAT ( ' Apex at X =',1PE15.8,', Y =',1PE15.8)
   50 FORMAT ( ' Apex at X =',1PE15.8,', Y =',1PE15.8,', Z =',1PE15.8)
   60 FORMAT (/' Undeformed Volume of Cavity is ',1PE15.8)

C
C ... REWIND EXODUS FILE TO BEGINNING OF TIMESTEPS
C
      IF (EXODUS .AND. ISDIS) THEN
         TIMEL = STMIN
         CALL GETDSP (CRD, DISP, NDIM, NUMNP, TIME, ITMSEL, 'R', ISTAT)
         IF (ISTAT .NE. 0) GO TO 140
         DO 70 IO=IOMIN, IOMAX
            WRITE (IO, 80)
   70    CONTINUE
   80    FORMAT (/,
     *      4X,'                 Cavity           Total',
     *     '            Timestep         Rate of',/
     *      4X,'Time             Volume           Change',
     *     '           Change           Change',/
     *      4X,'----             ------           ------',
     *     '           --------         -------')
C
         DELLAS = 0.0
   90    CONTINUE
         CALL GETDSP (CRD, DISP, NDIM, NUMNP, TIME, ITMSEL, 'S', ISTAT)
         IF (ISTAT .NE. 0) GO TO 140
         DELVOL = 0.0
         DO 100 NCAV = 1, NUMCAV
            IFLG = IFND(NCAV)
            IPTR = IPNESS(IFLG)

C     NOTE: Positive delcav = shrink in cavity volume

            IF (NDIM .EQ. 3) THEN
               CALL DVOL3D(CRD, DISP, LTNESS(IPTR),
     *            NEESS(IFLG), DELCAV, NDIM, NUMNP)
            ELSE
               CALL DVOL2D(CRD, DISP, LTNESS(IPTR),
     *            NEESS(IFLG), DELCAV, NDIM, AXI, NUMNP)
            END IF
C
            DELVOL =  DELVOL + DELCAV
  100    CONTINUE
         DELDEL = DELVOL - DELLAS
         IF (TREAD .EQ. TIMEL) THEN
            DO 110 IO=IOMIN, IOMAX
               WRITE (IO, 130) TREAD, TVOL-DELVOL, -DELVOL,
     *            -DELDEL
  110       CONTINUE
         ELSE
            DELRAT = DELDEL / (TREAD - TIMEL)
            DO 120 IO=IOMIN, IOMAX
               WRITE (IO, 130) TREAD, TVOL-DELVOL, -DELVOL,
     *            -DELDEL, -DELRAT
  120       CONTINUE
         END IF
  130    FORMAT (1X,5(1PE15.8,2X))
         DELLAS = DELVOL
         TIMEL  = TREAD
         GO TO 90
      END IF
  140 CONTINUE
      CALL INIINT (NUMCAV, 0, ICAV)
      NUMCAV = 0
      RETURN
      END
