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

C $Id: tablet.f,v 1.2 1998/07/14 18:20:07 gdsjaar Exp $
C $Log: tablet.f,v $
C Revision 1.2  1998/07/14 18:20:07  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:17:08  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:17:07  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]TABLET.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE TABLET (MP, ML, MS, MR, MSNAP, MCOM, ICOM, JCOM, CIN,
     &   RIN, IIN, KIN, IDUMP, N, IPOINT, COOR, IPBOUN, ILINE, LTYPE,
     &   NINT, FACTOR, LCON, iLBOUN, ISBOUN, ISIDE, NLPS, IFLINE,
     &   ILLIST, IREGN, IMAT, NSPR, IFSIDE, ISLIST, IRPB, IPBF, NPPF,
     &   IFPB, LISTPB, ILBF, NLPF, IFLB, LISTLB, ISBF, NSPF, IFSB,
     &   LISTSB, LINKP, LINKL, LINKS, LINKR, LINKM, LINKPB, LINKLB,
     &   LINKSB, REXTRM, IHOLDP, IHOLDL, IHOLDS, IHOLDR, IHOLDM, IHOLD2,
     &   IHOLD3, IWTPBF, IWTLBF, IWTSBF, IRGFLG, TITLE, NOROOM, DRWTAB,
     &   XX1, YY1, SCALE, CT, ST, X1, X2, Y1, Y2, ALPHA, DEV1, SNAP,
     &   SNAPDX, NSNAP, VAXVMS, TBZOOM, AXIST, WROTE, BATCH, VERSN,
     &   TIME1)
C***********************************************************************
C
C  SUBROUTINE TABLET = SUBROUTINE TO CONTROL DIGITIZE FUNCTIONS
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     FASTQ = A PROGRAM TO QUICKLY PREPARE QMESH INPUT
C
C***********************************************************************
C
C  SUBROUTINES CALLED:
C     DREAD  = SETS ALL PARAMETERS UP FOR READING FROM A DIGI-PAD
C     DPREAD = READS INPUT FROM A DIGI-PAD DIGITIZER
C     CLOSE  = FINDS THE CLOSEST EXISTING POINT TO THE MOUSE
C
C***********************************************************************
C
C  VARIABLES USED:
C     IANS   = LOGICAL RESPONSE FROM YES-NO QUESTION
C     TITLE  = MESH TITLE
C     DRWTAB = .TRUE. IF THE TABLET IS INITIALIZED TO A DRAWING
C     XX1    = DIGITIZATION PAD X COORDINATE OF POINT 1 (PAD INIT)
C     YY1    = DIGITIZATION PAD Y COORDINATE OF POINT 1 (PAD INIT)
C     X1     = USER X COORDINATE OF POINT 1 (PAD INIT)
C     Y1     = USER Y COORDINATE OF POINT 1 (PAD INIT)
C     X2     = USER X COORDINATE OF POINT 2 (PAD INIT)
C     Y2     = USER Y COORDINATE OF POINT 2 (PAD INIT)
C     X      = THE X LOCATION IN USER COORDINATES
C     Y      = THE Y LOCATION IN USER COORDINATES
C     SCALE  = THE SCALE FACTOR FROM DIGITIZED TO USER COORDINATES
C     CT     = THE COSINE OF THE ANGLE OF THE DRAWING ON THE PAD
C     ST     = THE SINE OF THE ANGLE OF THE DRAWING ON THE PAD
C     CHANGE = .TRUE. IF THE POINT USED IS THE CLOSEST POINT ONLY
C     SLIDE  = .TRUE. IF THE NEXT POINT IS TO HAVE THE CLOSEST POINT'S
C              COORDINATES, BUT NEW NUMBERING (SLIDE LINE USE)
C     NOROOM = .TRUE. IF THE AMOUNT OF DATA EXCEEDS DIMENSIONED LIMITS
C
C***********************************************************************
C
      DIMENSION IPOINT(MP), COOR(2, MP), IPBOUN(MP)
      DIMENSION ILINE(ML), LTYPE(ML), NINT(ML), FACTOR(ML), LCON(3, ML)
      DIMENSION ILBOUN(ML), ISBOUN(ML)
      DIMENSION ISIDE(MS), NLPS(MS), IFLINE(MS), ILLIST(MS*3)
      DIMENSION IREGN(MR), IMAT(MR), NSPR(MR), IFSIDE(MR), ISLIST(MR*4)
      DIMENSION IRPB(MR)
      DIMENSION IPBF(MP), NPPF(MP), IFPB(MP), LISTPB(2, MP)
      DIMENSION ILBF(ML), NLPF(ML), IFLB(ML), LISTLB(2, ML)
      DIMENSION ISBF(ML), NSPF(ML), IFSB(ML), LISTSB(2, ML)
      DIMENSION IWTPBF(3, MP), IWTLBF(3, ML), IWTSBF(3, ML)
      DIMENSION LINKP(2, MP), LINKL(2, ML), LINKS(2, MS), LINKR(2, MR)
      DIMENSION LINKM(2, (MS + MR)), LINKPB(2, MP), LINKLB(2, ML)
      DIMENSION LINKSB(2, ML)
      DIMENSION IHOLDP(2, MP), IHOLDL(ML*2), IHOLDR(2, MR)
      DIMENSION IHOLDM(2, (MS + MR)), IHOLD2(2, ML), IHOLD3(2, ML)
      DIMENSION IHOLDS(2, MS), IRGFLG(MR)
      DIMENSION N(29), REXTRM(4, MR), SNAPDX(2, MSNAP), NSNAP(2)
      DIMENSION KIN(MCOM), IIN(MCOM), RIN(MCOM)
C
      CHARACTER*72 TITLE, CIN(MCOM)
      CHARACTER DEV1*3, INTRNL*8, VERSN*9
C
      LOGICAL IANS, DRWTAB, ERR, NOROOM
      LOGICAL ALPHA
      LOGICAL SNAP, VAXVMS, TBZOOM, DRAWN, AXIST, WROTE, BATCH
C
      DRAWN=.FALSE.
C
C  GET THE BODY EXTREMES
C
      CALL GETEXT(MP, ML, MS, MR, N, IPOINT, COOR, ILINE, LTYPE,
     &   LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP,
     &   LINKL, LINKS, LINKR, REXTRM, XMIN1, XMAX1, YMIN1, YMAX1)
C
C  GET THE DEFAULT ZOOM AND GRID DEFINITIONS IF NOTHING HAS BEEN DEFINED
C
      IF (.NOT.TBZOOM) THEN
C
C  SET THE BODY EXTREMES AS THE ZOOM EXTREMES
C
         X1 = XMIN1
         X2 = XMAX1
         Y1 = YMIN1
         Y2 = YMAX1
         WRITE (*, 10010) X1, X2, Y1, Y2
C
C  GET THE DEFAULT TABLET INITIALIZATION
C
         CALL TABINT (X1, X2, Y1, Y2, CT, ST, SCALE, XX1, YY1, XX2, YY2,
     &      DRWTAB)
         TBZOOM = .TRUE.
      ELSE
         X1OLD = X1
         X2OLD = X2
         Y1OLD = Y1
         Y2OLD = Y2
      ENDIF
C
C  GET THE DEFAULT GRID IF NO GRID IS DEFINED
C
      IF (SNAP .AND. (NSNAP(1) .LT. 2 .OR. NSNAP(2) .LT. 2)) THEN
         CALL SNAPXY (MP, MSNAP, N(1), IPOINT, COOR, LINKP, SNAPDX,
     &      NSNAP)
      END IF
C
C  ENTER DIGITIZING OPTION
C
  100 CONTINUE
      IF (ICOM .GT. JCOM) THEN
         CALL MESAGE (' ')
         CALL FREFLD (IZ, IZ, 'ENTER TABLET OPTION: ', MCOM, IOSTAT,
     &      JCOM, KIN, CIN, IIN, RIN)
         ICOM = 1
      ENDIF
C
C  SPAWN A PROCESS
C
      IF ((CIN(ICOM)(1:2) .EQ. 'SP') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'sp')) THEN
         ICOM = ICOM + 1
         CALL SPAWN (VAXVMS)
C
C  SET THE SNAP FLAG ON
C
      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'S') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 's')) THEN
         ICOM = ICOM + 1
         IF (SNAP) THEN
            SNAP = .FALSE.
            CALL MESAGE (' ')
            CALL MESAGE ('SNAP COORDINATE DIGITIZAITON - DISABLED')
         ELSE
            SNAP = .TRUE.
            CALL MESAGE (' ')
            CALL MESAGE ('SNAP COORDINATE DIGITIZAITON - ENABLED')
            IF ((NSNAP(1) .LT. 2) .OR. (NSNAP(2) .LT. 2)) THEN
               CALL MESAGE ('PROPOSED SNAP GRID SPACING')
               SDX = (X2 - X1)/20.
               SDY = (Y2 - Y1)/15.
               WRITE (*, 10020) X1 - SDX, X2 + SDX, SDX, Y1 - SDX,
     &            Y2 + SDX, SDY
               CALL INTRUP ('IS THIS EXCEPTABLE',
     &            IANS, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN)
               IF (IANS) THEN
                  INDEX = 1
                  CALL UNISNP (MSNAP, SNAPDX, NSNAP, INDEX, X1 - SDX,
     &               X2 + SDX, SDX)
                  INDEX = 2
                  CALL UNISNP (MSNAP, SNAPDX, NSNAP, INDEX, Y1 - SDY,
     &               Y2 + SDY, SDY)
               ELSE
                  CALL MESAGE ('PLEASE DEFINE THE GRID USING UNIFORM'//
     &               ', UX, UY, X, OR Y OPTIONS')
               ENDIF
            ENDIF
         ENDIF
C
C  SHOW THE BUTTON DEFINITIONS
C
      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'A') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'a')) THEN
         ICOM = ICOM + 1
         IF(AXIST)THEN
            AXIST=.FALSE.
            CALL MESAGE('AXIS DRAWING - OFF')
         ELSE
            AXIST=.TRUE.
            CALL MESAGE('AXIS DRAWING - ON')
         ENDIF

C
C  SHOW THE BUTTON DEFINITIONS
C
      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'B') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'b')) THEN
         ICOM = ICOM + 1
         CALL HELP_FQ (2)
C
C  ADD UNIFORM Y SNAP GRID SPACINGS
C
      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'UY') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'uy')) THEN
         ICOM = ICOM + 1
         IF (ICOM .GT. JCOM) THEN
            CALL FREFLD (IZ, IZ, 'ENTER GRID YMIN, YMAX, AND Y GRID '//
     &         'SPACING:', MCOM, IOSTAT, JCOM, KIN, CIN, IIN, RIN)
            ICOM = 1
         ENDIF
         IF ((JCOM - ICOM + 1) .GE. 3) THEN
            INDEX = 2
            CALL UNISNP (MSNAP, SNAPDX, NSNAP, INDEX, RIN(ICOM),
     &         RIN(ICOM + 1), RIN(ICOM + 2))
            ICOM = ICOM + 3
            SNAP = .TRUE.
         ELSE
            CALL MESAGE ('NOT ENOUGH INFORMATION DEFINED TO ENTER'//
     &         ' UNIFORM Y GRID')
            CALL MESAGE ('NO ADDITIONAL Y GRID DEFINED')
         ENDIF
C
C  ADD UNIFORM X SNAP GRID SPACINGS
C
      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'UX') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'ux')) THEN
         ICOM = ICOM + 1
         IF (ICOM .GT. JCOM) THEN
            CALL FREFLD (IZ, IZ, 'ENTER GRID XMIN, XMAX, AND X GRID '//
     &         'SPACING:', MCOM, IOSTAT, JCOM, KIN, CIN, IIN, RIN)
            ICOM = 1
         ENDIF
         IF ((JCOM - ICOM + 1) .GE. 3) THEN
            INDEX = 1
            CALL UNISNP (MSNAP, SNAPDX, NSNAP, INDEX, RIN(ICOM),
     &         RIN(ICOM + 1), RIN(ICOM + 2))
            ICOM = ICOM + 3
            SNAP = .TRUE.
         ELSE
            CALL MESAGE ('NOT ENOUGH INFORMATION DEFINED TO ENTER'//
     &         ' UNIFORM X GRID')
            CALL MESAGE ('NO ADDITIONAL X GRID DEFINED')
         ENDIF
C
C  ADD UNIFORM SNAP GRID SPACINGS
C
      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'U') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'u')) THEN
         ICOM = ICOM + 1
         IF (ICOM .GT. JCOM) THEN
            CALL FREFLD (IZ, IZ, 'ENTER GRID XMIN, XMAX, YMIN, YMAX, '//
     &         'AND GRID SPACING:', MCOM, IOSTAT, JCOM, KIN, CIN, IIN,
     &         RIN)
            ICOM = 1
         ENDIF
         IF ((JCOM - ICOM + 1) .GE. 5) THEN
            SNAP = .TRUE.
            INDEX = 1
            CALL UNISNP (MSNAP, SNAPDX, NSNAP, INDEX, RIN(ICOM),
     &         RIN(ICOM + 1), RIN(ICOM + 4))
            INDEX = 2
            CALL UNISNP (MSNAP, SNAPDX, NSNAP, INDEX, RIN(ICOM + 2),
     &         RIN(ICOM + 3), RIN(ICOM + 4))
            ICOM = ICOM + 5
         ELSE
            CALL MESAGE ('NOT ENOUGH INFORMATION DEFINED TO ENTER'//
     &         ' UNIFORM GRID')
            CALL MESAGE ('NO GRID DEFINED')
         ENDIF
C
C  CLEAR ALL X GRID DEFINITIONS
C
      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'XC') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'xc')) THEN
         ICOM = ICOM + 1
         NSNAP(1) = 0
         CALL MESAGE ('ALL X SNAP GRID DEFINITIONS HAVE BEEN CLEARED')
         CALL MESAGE (' ')
C
C  ADD X SNAP GRID SPACINGS
C
      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'X') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'x')) THEN
         ICOM = ICOM + 1
         IF (ICOM .GT. JCOM) THEN
            CALL MESAGE ('ENTER X GRID VALUES IN ANY ORDER:')
            CALL MESAGE ('SEPERATE VALUES BY DELIMITERS OR RETURN KEY')
            CALL MESAGE ('HIT RETURN TO END INPUT')
         ENDIF
  110    CONTINUE
         IF (ICOM .GT. JCOM) THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN, IIN,
     &         RIN)
            ICOM = 1
         ENDIF
         IF (JCOM .GE. ICOM) THEN
            SNAP = .TRUE.
            INDEX = 1
            DO 120 I = ICOM, JCOM
               IF (KIN(I) .GT. 0) THEN
                  ICOM = ICOM + 1
                  CALL ADDSNP (MSNAP, SNAPDX, NSNAP, INDEX, RIN(I), ERR)
                  IF (ERR) THEN
                     WRITE (*, 10000) 'X', RIN(I - 1)
                     GO TO 130
                  ENDIF
               ELSE
                  GO TO 130
               ENDIF
  120       CONTINUE
            GO TO 110
  130       CONTINUE
         ENDIF
C
C  CLEAR ALL Y GRID DEFINITIONS
C
      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'YC') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'yc')) THEN
         ICOM = ICOM + 1
         NSNAP(2) = 0
         CALL MESAGE ('ALL Y SNAP GRID DEFINITIONS HAVE BEEN CLEARED')
         CALL MESAGE (' ')
C
C  ADD Y SNAP GRID SPACINGS
C
      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'Y') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'y')) THEN
         ICOM = ICOM + 1
         IF (ICOM .GT. JCOM) THEN
            CALL MESAGE ('ENTER Y GRID VALUES IN ANY ORDER:')
            CALL MESAGE ('SEPERATE VALUES BY DELIMITERS OR RETURN KEY')
            CALL MESAGE ('HIT RETURN TO END INPUT')
         ENDIF
  140    CONTINUE
         IF (ICOM .GT. JCOM) THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN, IIN,
     &         RIN)
            ICOM = 1
         ENDIF
         IF (JCOM .GE. ICOM) THEN
            SNAP = .TRUE.
            INDEX = 2
            DO 150 I = ICOM, JCOM
               IF (KIN(I) .GT. 0) THEN
                  ICOM = ICOM + 1
                  CALL ADDSNP (MSNAP, SNAPDX, NSNAP, INDEX, RIN(I), ERR)
                  IF (ERR) THEN
                     WRITE (*, 10000) 'Y', RIN(I - 1)
                     GO TO 160
                  ENDIF
               ELSE
                  GO TO 160
               ENDIF
  150       CONTINUE
            GO TO 140
  160       CONTINUE
         ENDIF
C
C  SET ZOOM LIMITS FOR PLOTTING
C
      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'P') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'p')) THEN
         ICOM = ICOM + 1
         CALL MESAGE ('GRID POINTS ADDED FOR ALL POINT '//
     &      'COORDINATE VALUES')
         CALL SNAPXY (MP, MSNAP, N(1), IPOINT, COOR, LINKP, SNAPDX,
     &      NSNAP)
         SNAP = .TRUE.
C
C  SET ZOOM LIMITS FOR PLOTTING
C
      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'Z') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'z')) THEN
         ICOM = ICOM + 1
         CALL ZOOMLT(MCOM, ICOM, JCOM, CIN, RIN, IIN, KIN,
     &      IDUMP, DRAWN, ALPHA, DEV1, X1, X2, Y1, Y2, X11, X22, Y11,
     &      Y22, XMIN1, XMAX1, YMIN1, YMAX1, XMIN, XMAX, YMIN, YMAX)
         DRAWN = .FALSE.
         IF (DRWTAB) THEN
            CALL MESAGE ('TABLET EXTREMES REMAIN TIED TO DRAWING'//
     &         ' INITIALIZATION')
            CALL MESAGE ('SCREEN PLOTTING ZOOM CHANGED')
            X1OLD = XMIN
            X2OLD = XMAX
            Y1OLD = YMIN
            Y2OLD = YMAX
         ELSE
            X1 = XMIN
            X2 = XMAX
            Y1 = YMIN
            Y2 = YMAX
            CALL TABINT (X1, X2, Y1, Y2, CT, ST, SCALE, XX1, YY1,
     &         XX2, YY2, DRWTAB)
            CALL MESAGE ('SCREEN PLOTTING ZOOM CHANGED')
            CALL MESAGE ('TABLET EXTREMES RESET TO ZOOM LIMITS')
         ENDIF
C
C  INITIALIZE DIGITIZING PAD
C
      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'I') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'i')) THEN
         ICOM = ICOM + 1
         CALL INITDG (MCOM, ICOM, JCOM, CIN, RIN, IIN, KIN, IDUMP, XX1,
     &      YY1, SCALE, CT, ST, X1, X2, Y1, Y2, DRWTAB, SNAP)
         IF (DRWTAB) THEN
            X1OLD = X1
            X2OLD = X2
            Y1OLD = Y1
            Y2OLD = Y2
            XMIN1 = X1
            XMAX1 = X2
            YMIN1 = Y1
            YMAX1 = Y2
         END IF
C
C  CLEAR ALL GRID DEFINITIONS
C
      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'C') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'c')) THEN
         ICOM = ICOM + 1
         NSNAP(1) = 0
         NSNAP(2) = 0
         CALL MESAGE ('ALL SNAP GRID DEFINITIONS HAVE BEEN CLEARED')
         CALL MESAGE (' ')
C
C  DIGITIZING OPTION
C
      ELSE IF ((CIN(ICOM)(1:1) .EQ. 'D') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'd')) THEN
         ICOM = ICOM + 1
C
C  GENERATE A DEFAULT SNAP GRID IF NEEDED
C
         IF ((SNAP).AND.((NSNAP(1) .LT. 2) .OR. (NSNAP(2) .LT. 2))) THEN
            NSNAP(1) = 0
            NSNAP(2) = 0
            CALL MESAGE ('PROPOSED SNAP GRID SPACING')
            SDX = (X2 - X1)/20.
            WRITE (INTRNL, '(E8.1)') SDX
            READ (INTRNL, '(E8.1)') SDX
            SDY = (Y2 - Y1)/15.
            WRITE (INTRNL, '(E8.1)') SDY
            READ (INTRNL, '(E8.1)') SDY
            SDX = MAX(SDX, SDY)
            SDY = SDX
            WRITE (*, 10020) X1 - SDX, X2 + SDX, SDX, Y1 - SDX,
     &         Y2 + SDX, SDY
            CALL INTRUP ('IS THIS EXCEPTABLE',
     &         IANS, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN)
            IF (IANS) THEN
               INDEX = 1
               CALL UNISNP (MSNAP, SNAPDX, NSNAP, INDEX, X1 - SDX,
     &            X2 + SDX, SDX)
               INDEX = 2
               CALL UNISNP (MSNAP, SNAPDX, NSNAP, INDEX, Y1 - SDY,
     &            Y2 + SDY, SDY)
            ELSE
               CALL INTRUP ('WOULD YOU CARE TO CONTINUE DIGITIZING '//
     &            'WITHOUT GRID SNAP', IANS, MCOM, ICOM, JCOM, CIN,
     &            IIN, RIN, KIN)
               IF (IANS) THEN
                  SNAP = .FALSE.
               ELSE
                  CALL MESAGE ('PLEASE DEFINE THE GRID USING '//
     &               'UNIFORM, UX, UY, X, OR Y OPTIONS')
                  GO TO 100
               ENDIF
            ENDIF
         ENDIF
C
C  NOW ENTER THE MOUSE CONTROL
C
         CALL DIGIT (MP, ML, MS, MR, MSNAP, MCOM, ICOM, JCOM, CIN, RIN,
     &      IIN, KIN, IDUMP, N, IPOINT, COOR, IPBOUN, ILINE, LTYPE,
     &      NINT, FACTOR, LCON, ILBOUN, ISBOUN, ISIDE, NLPS, IFLINE,
     &      ILLIST, IREGN, IMAT, NSPR, IFSIDE, ISLIST, IRPB, IPBF, NPPF,
     &      IFPB, LISTPB, ILBF, NLPF, IFLB, LISTLB, ISBF, NSPF, IFSB,
     &      LISTSB, LINKP, LINKL, LINKS, LINKR, LINKM, LINKPB, LINKLB,
     &      LINKSB, IHOLDP, IHOLDL, IHOLDS, IHOLDR, IHOLDM, IHOLD2,
     &      IHOLD3, IWTPBF, IWTLBF, IWTSBF, IRGFLG, TITLE, NOROOM, XX1,
     &      YY1, SCALE, CT, ST, X1, X2, Y1, Y2, X11, X22, Y11, Y22,
     &      XMIN1, XMAX1, YMIN1, YMAX1, XMIN2, XMAX2, YMIN2, YMAX2,
     &      X1OLD, X2OLD, Y1OLD, Y2OLD, ALPHA, DEV1, SNAP, SNAPDX,
     &      NSNAP, DRWTAB, AXIST)
         DRAWN = .TRUE.
         WROTE = .FALSE.
C
C  GO GET MORE ROOM IF NEEDED AND GO STRAIGHT BACK INTO DIGITIZING
C
         IF (NOROOM) THEN
            JCOM = 1
            ICOM = 1
            CIN(1) = 'DIG'
            RETURN
         ENDIF
C
C  RETURN FROM DIGITIZING
C
      ELSE IF (CIN(ICOM)(1:1) .EQ. ' ') THEN
         ICOM = ICOM + 1
         RETURN

C
C  EXIT OPTION - EXITS FASTQ
C
      ELSE IF ((CIN(ICOM)(1:2) .EQ. 'EX') .OR.
     &   (CIN(ICOM)(1:2) .EQ. 'ex')) THEN
         ICOM = ICOM + 1
         CALL STRLNG (CIN(ICOM), LEN)
         IF (((LEN .GT. 1) .AND. (CIN(ICOM)(2:2) .NE. 'X')) .OR.
     &      ((LEN .GT. 1) .AND. (CIN(ICOM)(2:2) .NE. 'x'))) THEN
            CALL HELP_FQ (14)
         ELSE
            CALL FEXIT (WROTE, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN,
     &         TIME1, BATCH, VERSN)
         ENDIF
C
C  WRITE OUT THE HELP MESSAGE
C
      ELSE
         ICOM = ICOM + 1
         CALL HELP_FQ (14)
      ENDIF
      GO TO 100
C
10000 FORMAT(' THE LAST SUCESSFULL ', A1, ' INPUT WAS: ', G14.7)
10010 FORMAT(' THE TABLET (AND PLOTTING) LIMITS ARE DEFAULTED TO:', /
     &   '         XMIN: ', G14.7, /,
     &   '         XMAX: ', G14.7, /,
     &   '         YMIN: ', G14.7, /,
     &   '         YMAX: ', G14.7)
10020 FORMAT(' FOR THE DEFAULT GRID, THE MINIMUM X IS: ', G14.7, /,
     &   '                       THE MAXIMUM X IS: ', G14.7, /,
     &   '                          THE X STEP IS: ', G14.7, /,
     &   '                       THE MINIMUM Y IS: ', G14.7, /,
     &   '                       THE MAXIMUM Y IS: ', G14.7, /,
     &   '                          THE Y STEP IS: ', G14.7)
C
      END
