C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Governement retains certain rights in this software.
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
C        * Neither the name of Sandia Corporation nor the names of its
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

C $Id: cavity.f,v 1.7 2007/03/21 20:12:37 gdsjaar Exp $
C $Log: cavity.f,v $
C Revision 1.7  2007/03/21 20:12:37  gdsjaar
C Several commands which can work on the deformed geometry were only
C checking whether the file was an exodus file (had timesteps) when
C requesting deformed coordinates.  Changed to also check whether the
C file had valid displacments also.
C
C Revision 1.6  1999/02/16 21:37:58  gdsjaar
C Converted to read exodusII database format.  Somewhat tested, not
C ready for production yet.
C
C Revision 1.5  1993/04/08 20:07:28  gdsjaar
C Changed cavity volume output to have correct sign
C
c Revision 1.4  1992/12/11  22:34:12  gdsjaar
c Fixed problem with incorrect determination of cavity center in 2d
c
c Revision 1.3  1992/07/20  22:38:09  gdsjaar
c Multiple cavity volume was using VOLUME (last cavity volume) instead
c of TVOL (total cavity volume).
c
c Revision 1.2  1992/07/20  22:22:23  gdsjaar
c Initialize variable DELLAS - unset before
c
c Revision 1.1.1.1  1991/02/21  15:42:22  gdsjaar
c NUMBERS: Greg Sjaardema, initial Unix release
c
c Revision 1.1  1991/02/21  15:42:21  gdsjaar
c Initial revision
c
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
   40 FORMAT ( ' Apex at X =',1PE10.3,', Y =',1PE10.3)
   50 FORMAT ( ' Apex at X =',1PE10.3,', Y =',1PE10.3,', Z =',1PE10.3)
   60 FORMAT (/' Undeformed Volume of Cavity is ',1PE10.3)

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
     *      4X,'           Cavity      Total      Timestep    Rate of',/
     *      4X,'Time       Volume      Change      Change      Change',/
     *      4X,'----       ------      ------     --------    -------')
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
  130    FORMAT (1X,5(1PE10.3,2X))
         DELLAS = DELVOL
         TIMEL  = TREAD
         GO TO 90
      END IF
  140 CONTINUE
      CALL INIINT (NUMCAV, 0, ICAV)
      NUMCAV = 0
      RETURN
      END
