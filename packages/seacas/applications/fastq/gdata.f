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

C $Id: gdata.f,v 1.3 1999/07/01 13:48:26 gdsjaar Exp $
C $Log: gdata.f,v $
C Revision 1.3  1999/07/01 13:48:26  gdsjaar
C Fixed another negative array index.
C
C Revision 1.2  1998/07/14 18:18:59  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:07:51  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:07:49  gdsjaar
c Initial revision
c
C
CC* FILE: [.MAIN]GDATA.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GDATA (MP, ML, MS, MR, MSC, MCOM, ICOM, JCOM, CIN,
     &   RIN, IIN, KIN, IDUMP, N, IPOINT, COOR, IPBOUN, ILINE, LTYPE,
     &   NINT, FACTOR, LCON, ILBOUN, ISBOUN, ISIDE, NLPS, IFLINE,
     &   ILLIST, IBARST, JMAT, JCENT, NLPB, JFLINE, JLLIST, IREGN,
     &   IMAT, NSPR, IFSIDE, ISLIST, LINKP, LINKL, LINKS, LINKB, LINKR,
     &   LINKSC, REXTRM, RSIZE, SCHEME, DEFSCH, DEFSIZ, TITLE, LABP,
     &   LABL, LABR, AXISD, LABMD, LABI, LABF, LABPB, LABLB, LABSBD,
     &   LABSC, LABSZ, FULL, IDEV, ALPHA, DEV1, VAXVMS, VERSN, WROTE,
     &   TIME1, HARDPL, BATCH)
C***********************************************************************
C
C  GDATA = SUBROUTINE TO INPUT LIGHT TABLE POINTS
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     FASTQ = A PROGRAM TO QUICKLY PREPARE QMESH INPUT
C
C***********************************************************************
C
      DIMENSION IPOINT (MP), COOR (2, MP), IPBOUN (MP)
      DIMENSION ILINE (ML), LTYPE (ML), NINT (ML)
      DIMENSION FACTOR (ML), LCON (3, ML)
      DIMENSION ILBOUN (ML), ISBOUN (ML)
      DIMENSION ISIDE (MS), NLPS (MS), IFLINE (MS), ILLIST (MS*3)
      DIMENSION IBARST (MS), JMAT (MS), JCENT (MS)
      DIMENSION NLPB (MS), JFLINE (MS), JLLIST (MS*3)
      DIMENSION IREGN (MR), IMAT (MR), NSPR (MR), IFSIDE (MR)
      DIMENSION ISLIST (MR*4)
      DIMENSION SCHEME (MSC), RSIZE (MR)
      DIMENSION LINKP (2, MP), LINKL (2, ML), LINKS (2, MS)
      DIMENSION LINKB (2, MS), LINKR (2, MR), LINKSC (2, MR)
      DIMENSION REXTRM (4, MR), N (29)
      DIMENSION KIN (MCOM), IIN (MCOM), RIN (MCOM)
      DIMENSION IDEV (2), III (1)
C
      CHARACTER*72 SCHEME, DEFSCH, TITLE, CIN (MCOM)
      CHARACTER*3 DEV1, VERSN*9
C
      LOGICAL DRAWN, FLAG, GOPLOT, ALPHA
      LOGICAL ADDLNK, VAXVMS, WROTE
      LOGICAL LABP, LABL, LABR, AXISD, LABMD
      LOGICAL LABI, LABF, LABPB, LABLB, LABSBD
      LOGICAL FULL, LABSC, LABSZ, GETMAX, TEST
      LOGICAL NUMPLT, HARDPL, BATCH, FOUND
C
      IZ = 0
      DRAWN = .FALSE.
      ADDLNK = .FALSE.
      GETMAX = .FALSE.
      TEST = .FALSE.
      NUMPLT = .FALSE.
C
C  FLAG ALL THE DATA TO BE PLOTTED
C
      FLAG = .TRUE.
      CALL FLAGD (MP, N (18), LINKP, IPOINT, FLAG)
      CALL FLAGD (ML, N (19), LINKL, ILINE, FLAG)
      CALL FLAGD (MR, N (22), LINKR, IREGN, FLAG)
C
C  GET THE REGION AND BODY EXTREMES
C
      CALL GETEXT (MP, ML, MS, MR, N, IPOINT, COOR, ILINE, LTYPE,
     &   LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP,
     &   LINKL, LINKS, LINKR, REXTRM, XMIN, XMAX, YMIN, YMAX)
C
      XMIN1 = XMIN
      XMAX1 = XMAX
      YMIN1 = YMIN
      YMAX1 = YMAX
      HXMIN = XMIN
      HYMIN = YMIN
      HXMAX = XMAX
      HYMAX = YMAX
C
C  ENTER GRAPHICS OPTION
C
  100 CONTINUE
      IF (ICOM .GT. JCOM) THEN
         CALL MESAGE (' ')
         CALL FREFLD (IZ, IZ, 'ENTER GRAPHICS OPTION: ', MCOM,
     &      IOSTAT, JCOM, KIN, CIN, IIN, RIN)
         ICOM = 1
      ENDIF
C
C  TOGGLE DRAWING OF THE AXIS
C
      IF ( (CIN (ICOM) (1:1) .EQ. 'A') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'a')) THEN
         ICOM = ICOM + 1
         IF (AXISD) THEN
            AXISD = .FALSE.
            CALL MESAGE ('AXIS DRAWING - OFF')
         ELSE
            AXISD = .TRUE.
            CALL MESAGE ('AXIS DRAWING - ON')
         ENDIF
C
C  TOGGLE THE FACTOR NUMBERS
C
      ELSEIF ( (CIN (ICOM) (1:2) .EQ. 'FA') .OR.
     &   (CIN (ICOM) (1:2) .EQ. 'fa')) THEN
         ICOM = ICOM + 1
         IF (LABF) THEN
            LABF = .FALSE.
            CALL MESAGE ('FACTOR LABELS - OFF')
         ELSE
            LABF = .TRUE.
            CALL MESAGE ('FACTOR LABELS - ON')
         ENDIF
C
C  TOGGLE THE FULL DISPLAY OF PROPERTIES
C
      ELSEIF ( (CIN (ICOM) (1:1) .EQ. 'F') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'f')) THEN
         ICOM = ICOM + 1
         IF (FULL) THEN
            FULL = .FALSE.
            CALL MESAGE ('FULL DISPLAY OF PROPERTIES - OFF')
         ELSE
            FULL = .TRUE.
            CALL MESAGE ('FULL DISPLAY OF PROPERTIES - ON')
         ENDIF
C
C  TOGGLE THE SCHEME DISPLAY
C
      ELSEIF ( (CIN (ICOM) (1:2) .EQ. 'SC') .OR.
     &   (CIN (ICOM) (1:2) .EQ. 'sc')) THEN
         ICOM = ICOM + 1
         IF (LABSC) THEN
            LABSC = .FALSE.
            CALL MESAGE ('SCHEME LABELS - OFF')
         ELSE
            LABSC = .TRUE.
            CALL MESAGE ('SCHEME LABELS - ON')
         ENDIF
C
C  TOGGLE THE SCHEME DISPLAY
C
      ELSEIF ( (CIN (ICOM) (1:2) .EQ. 'SI') .OR.
     &   (CIN (ICOM) (1:2) .EQ. 'si')) THEN
         ICOM = ICOM + 1
         IF (LABSZ) THEN
            LABSZ = .FALSE.
            CALL MESAGE ('ELEMENT SIZE LABELS - OFF')
         ELSE
            LABSZ = .TRUE.
            CALL MESAGE ('ELEMENT SIZE LABELS - ON')
         ENDIF
C
C  TOGGLE THE MATERIAL NUMBERS
C
      ELSEIF ( (CIN (ICOM) (1:1) .EQ. 'M') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'm')) THEN
         ICOM = ICOM + 1
         IF (LABMD) THEN
            LABMD = .FALSE.
            CALL MESAGE ('MATERIAL LABELS - OFF')
         ELSE
            LABMD = .TRUE.
            CALL MESAGE ('MATERIAL LABELS - ON')
         ENDIF
C
C  ENTER LINE INTERVALS
C
      ELSEIF ( (CIN (ICOM) (1:2) .EQ. 'II') .OR.
     &   (CIN (ICOM) (1:2) .EQ. 'ii')) THEN
         ICOM = ICOM + 1
         IF (ICOM .GT. JCOM) THEN
            CALL MESAGE ('ENTER LINE INTERVALS IN THE FOLLOWING '//
     &         'FORMAT:')
            CALL MESAGE ('[ LINE NO.  (OR NEG SIDE NO.),  INTERVALS ]')
            CALL MESAGE ('HIT RETURN TO END INPUT')
         ENDIF
  110    CONTINUE
         IF (ICOM .GT. JCOM) THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN,
     &         IIN, RIN)
            ICOM = 1
         ENDIF
         CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2, IFOUND)
         III (1) = I1
         IF (IFOUND .GT. 0) THEN
            CALL ININTR (ML, MS, 1, I2, III, N (19), N (20), NINT,
     &         NLPS, IFLINE, ILLIST, LINKL, LINKS, ADDLNK)
            GOTO 110
         ENDIF
C
C  ENTER LINE FACTORS
C
      ELSEIF ( (CIN (ICOM) (1:2) .EQ. 'IF') .OR.
     &   (CIN (ICOM) (1:2) .EQ. 'if')) THEN
         ICOM = ICOM + 1
         IF (ICOM .GT. JCOM) THEN
            CALL MESAGE ('ENTER LINE FACTORS IN THE FOLLOWING FORMAT:')
            CALL MESAGE ('[ LINE NO.  (OR NEG. SIDE NO., ) FACTOR ]')
            CALL MESAGE ('HIT RETURN TO END INPUT')
         ENDIF
  120    CONTINUE
         IF (ICOM .GT. JCOM) THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN,
     &         IIN, RIN)
            ICOM = 1
         ENDIF
         CALL GETI1R (MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN, I1, R1,
     &      IFOUND)
         III (1) = I1
         IF (IFOUND .GT. 0) THEN
            CALL INFACT (ML, MS, 1, R1, III, N (19), N (20), FACTOR,
     &         NLPS, IFLINE, ILLIST, LINKL, LINKS, ADDLNK)
            GOTO 120
         ENDIF
C
C  TOGGLE THE INTERVAL NUMBERS
C
      ELSEIF ( (CIN (ICOM) (1:1) .EQ. 'I') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'i')) THEN
         ICOM = ICOM + 1
         IF (LABI) THEN
            LABI = .FALSE.
            CALL MESAGE ('INTERVAL LABELS - OFF')
         ELSE
            LABI = .TRUE.
            CALL MESAGE ('INTERVAL LABELS - ON')
         ENDIF
C
C  TOGGLE THE NODE BOUNDARY NUMBERS
C
      ELSEIF ( (CIN (ICOM) (1:1) .EQ. 'N') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'n')) THEN
         ICOM = ICOM + 1
         IF (LABLB) THEN
            LABLB = .FALSE.
            CALL MESAGE ('LINE BOUNDARY LABELS - OFF')
         ELSE
            LABLB = .TRUE.
            CALL MESAGE ('LINE BOUNDARY LABELS - ON')
         ENDIF
C
C  FLAG LINES TO BE PROCESSED
C
      ELSEIF ( (CIN (ICOM) (1:2) .EQ. 'LP') .OR.
     &   (CIN (ICOM) (1:2) .EQ. 'lp')) THEN
         ICOM = ICOM + 1
         HXMIN = XMIN
         HYMIN = YMIN
         HXMAX = XMAX
         HYMAX = YMAX
         GETMAX = .TRUE.
         FOUND = .FALSE.
         FLAG = .FALSE.
         CALL FLAGD (MP, N (18), LINKP, IPOINT, FLAG)
         CALL FLAGD (ML, N (19), LINKL, ILINE, FLAG)
         CALL FLAGD (MR, N (22), LINKR, IREGN, FLAG)
         GOPLOT = .FALSE.
         CALL MESAGE ('PLOT LINES FROM <I1> TO <I2>')
         CALL MESAGE ('HIT RETURN TO END INPUT')
  130    CONTINUE
         IF (ICOM .GT. JCOM) THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN,
     &         IIN, RIN)
            ICOM = 1
         ENDIF
         CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2, IFOUND)
         IF (IFOUND .GT. 0) THEN
            IF (I1 .GT. 0) THEN
               CALL CHECK (I1, I2, N (19))
C
C  FLAG ALL DATA ASSOCIATED WITH THE LINES
C
               DO 150 I = I1, I2
                  CALL LTSORT (ML, LINKL, I, KK, ADDLNK)
                  IF (KK .GT. 0) THEN
                     GOPLOT = .TRUE.
C
                     DO 140 L = 1, 3
                        IF (LCON (L, KK) .GT. 0) THEN
                           CALL LTSORT (MP, LINKP, LCON (L, KK),
     &                        LL, ADDLNK)
                           IF (LL .GT. 0) THEN
                              IPOINT (LL) = -IABS (IPOINT (LL))
                              IF (.NOT.FOUND) THEN
                                 XMAX = COOR (1, LL)
                                 XMIN = COOR (1, LL)
                                 YMAX = COOR (2, LL)
                                 YMIN = COOR (2, LL)
                                 FOUND = .TRUE.
                              ENDIF
                           ENDIF
                        ENDIF
  140                CONTINUE
C
                     CALL DLINE (MP, ML, COOR, LINKP, ILINE (KK),
     &                  LTYPE (KK), LCON (1, KK), LCON (2, KK),
     &                  LCON (3, KK), NUMPLT, X1, Y1, TEST, GETMAX,
     &                  XMIN, XMAX, YMIN, YMAX)
                     ILINE (KK) = -IABS (ILINE (KK))
                  ENDIF
  150          CONTINUE
               GOTO 130
            ENDIF
         ENDIF
         GETMAX = .FALSE.
         IF (GOPLOT) THEN
C
C  PLOT THE LINE DATA THAT HAS BEEN FLAGGED
C
            IF (ALPHA) THEN
               CALL MESAGE ('NO PLOTTING POSSIBLE ON ALPHANUMERIC '//
     &            'TERMINAL')
            ELSE
               CALL PDATA (MP, ML, MR, MSC, IPOINT, COOR, IPBOUN,
     &            ILINE, LTYPE, NINT, LCON, FACTOR, ILBOUN, ISBOUN,
     &            IREGN, IMAT, LINKP, LINKL, LINKR, LINKSC, RSIZE,
     &            SCHEME, DEFSCH, DEFSIZ, REXTRM, N, LABP, LABL, LABR,
     &            FULL, LABMD, LABI, LABF, LABPB, LABLB, LABSBD, LABSC,
     &            LABSZ, AXISD, TITLE, XMIN, XMAX, YMIN, YMAX, XX1,
     &            YY1, XX2, YY2, DEV1, VERSN)
               DRAWN = .TRUE.
            ENDIF
         ELSE
            XMIN = HXMIN
            YMIN = HYMIN
            XMAX = HXMAX
            YMAX = HYMAX
         ENDIF
         XMIN1 = XMIN
         YMIN1 = YMIN
         XMAX1 = XMAX
         YMAX1 = YMAX
C
C  TOGGLE THE LINE NUMBERS
C
      ELSEIF ( (CIN (ICOM) (1:1) .EQ. 'L') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'l')) THEN
         ICOM = ICOM + 1
         IF (LABL) THEN
            LABL = .FALSE.
            CALL MESAGE ('LINE LABELS - OFF')
         ELSE
            LABL = .TRUE.
            CALL MESAGE ('LINE LABELS - ON')
         ENDIF
C
C  TOGGLE THE POINT BOUNDARY NUMBERS
C
      ELSEIF ( (CIN (ICOM) (1:2) .EQ. 'PB') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'pb')) THEN
         ICOM = ICOM + 1
         IF (LABPB) THEN
            LABPB = .FALSE.
            CALL MESAGE ('POINT BOUNDARY LABELS - OFF')
         ELSE
            LABPB = .TRUE.
            CALL MESAGE ('POINT BOUNDARY LABELS - ON')
         ENDIF
C
C  TOGGLE THE SIDE BOUNDARY NUMBERS
C
      ELSEIF ( (CIN (ICOM) (1:2) .EQ. 'EB') .OR.
     &   (CIN (ICOM) (1:2) .EQ. 'eb')) THEN
         ICOM = ICOM + 1
         IF (LABSBD) THEN
            LABSBD = .FALSE.
            CALL MESAGE ('ELEMENT BOUNDARY LABELS - OFF')
         ELSE
            LABSBD = .TRUE.
            CALL MESAGE ('ELEMENT BOUNDARY LABELS - ON')
         ENDIF
C
C  TOGGLE THE POINT NUMBERS
C
      ELSEIF ( (CIN (ICOM) (1:2) .EQ. 'PO') .OR.
     &   (CIN (ICOM) (1:2) .EQ. 'po')) THEN
         ICOM = ICOM + 1
         IF (LABP) THEN
            LABP = .FALSE.
            CALL MESAGE ('POINT LABELS - OFF')
         ELSE
            LABP = .TRUE.
            CALL MESAGE ('POINT LABELS - ON')
         ENDIF
C
C  TOGGLE THE REGION NUMBERS
C
      ELSEIF ( (CIN (ICOM) (1:2) .EQ. 'RE') .OR.
     &   (CIN (ICOM) (1:2) .EQ. 're')) THEN
         ICOM = ICOM + 1
         IF (LABR) THEN
            LABR = .FALSE.
            CALL MESAGE ('REGION LABELS - OFF')
         ELSE
            LABR = .TRUE.
            CALL MESAGE ('REGION LABELS - ON')
         ENDIF
C
C  SPAWN A PROCESS
C
      ELSEIF ( (CIN (ICOM) (1:3) .EQ. 'SPA') .OR.
     &   (CIN (ICOM) (1:3) .EQ. 'spa')) THEN
         ICOM = ICOM + 1
         CALL SPAWN (VAXVMS)
C
C  FLAG SIDES TO BE PROCESSED
C
      ELSEIF ( (CIN (ICOM) (1:2) .EQ. 'SP') .OR.
     &   (CIN (ICOM) (1:2) .EQ. 'sp')) THEN
         ICOM = ICOM + 1
         HXMIN = XMIN
         HYMIN = YMIN
         HXMAX = XMAX
         HYMAX = YMAX
         XMIN = 100000.
         YMIN = 100000.
         XMAX = -100000.
         YMAX = -100000.
         FOUND = .FALSE.
         GETMAX = .TRUE.
         FLAG = .FALSE.
         CALL FLAGD (MP, N (18), LINKP, IPOINT, FLAG)
         CALL FLAGD (ML, N (19), LINKL, ILINE, FLAG)
         CALL FLAGD (MR, N (22), LINKR, IREGN, FLAG)
         GOPLOT = .FALSE.
         CALL MESAGE ('PLOT SIDES FROM <I1> TO <I2>')
         CALL MESAGE ('HIT RETURN TO END INPUT')
  160    CONTINUE
         IF (ICOM .GT. JCOM) THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN,
     &         IIN, RIN)
            ICOM = 1
         ENDIF
         CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2, IFOUND)
         IF (IFOUND .GT. 0) THEN
            IF (I1 .GT. 0) THEN
               CALL CHECK (I1, I2, N (20))
C
C  FLAG ALL DATA ASSOCIATED WITH THE SIDES
C
               DO 190 I = I1, I2
                  CALL LTSORT (MS, LINKS, I, JJ, ADDLNK)
                  IF (JJ .GT. 0) THEN
                     GOPLOT = .TRUE.
                     DO 180 K = IFLINE (JJ), IFLINE (JJ)+NLPS (JJ)-1
                        CALL LTSORT (ML, LINKL, ILLIST (K), KK, ADDLNK)
                        IF (KK .GT. 0) THEN
                           CALL DLINE (MP, ML, COOR, LINKP, ILINE (KK),
     &                        LTYPE (KK), LCON (1, KK), LCON (2, KK),
     &                        LCON (3, KK), NUMPLT, X1, Y1, TEST,
     &                        GETMAX, XMIN, XMAX, YMIN, YMAX)
                           ILINE (KK) = -IABS (ILINE (KK))
                           DO 170 L = 1, 3
                              IF (LCON (L, KK) .GT. 0) THEN
                                 CALL LTSORT (MP, LINKP, LCON (L, KK),
     &                              LL, ADDLNK)
                                 IF (LL .GT. 0) THEN
                                    IPOINT (LL) = -IABS (IPOINT (LL))
                                    IF (.NOT.FOUND) THEN
                                       XMAX = COOR (1, LL)
                                       XMIN = COOR (1, LL)
                                       YMAX = COOR (2, LL)
                                       YMIN = COOR (2, LL)
                                       FOUND = .TRUE.
                                    ENDIF
                                 ENDIF
                              ENDIF
  170                      CONTINUE
                        ENDIF
  180                CONTINUE
                  ENDIF
  190          CONTINUE
               GOTO 160
            ENDIF
         ENDIF
         GETMAX = .FALSE.
         IF (GOPLOT) THEN
C
C  PLOT THE SIDE DATA THAT HAS BEEN FLAGGED
C
            IF (ALPHA) THEN
               CALL MESAGE ('NO PLOTTING POSSIBLE ON ALPHANUMERIC ' //
     &            'TERMINAL')
            ELSE
               CALL PDATA (MP, ML, MR, MSC, IPOINT, COOR, IPBOUN,
     &            ILINE, LTYPE, NINT, LCON, FACTOR, ILBOUN, ISBOUN,
     &            IREGN, IMAT, LINKP, LINKL, LINKR, LINKSC, RSIZE,
     &            SCHEME, DEFSCH, DEFSIZ, REXTRM, N, LABP, LABL, LABR,
     &            FULL, LABMD, LABI, LABF, LABPB, LABLB, LABSBD, LABSC,
     &            LABSZ, AXISD, TITLE, XMIN, XMAX, YMIN, YMAX, XX1,
     &            YY1, XX2, YY2, DEV1, VERSN)
               DRAWN = .TRUE.
            ENDIF
         ELSE
            XMIN = HXMIN
            YMIN = HYMIN
            XMAX = HXMAX
            YMAX = HYMAX
         ENDIF
         XMIN1 = XMIN
         YMIN1 = YMIN
         XMAX1 = XMAX
         YMAX1 = YMAX
C
C  SHOW STATUS OF ALL TOGGLES
C
      ELSEIF ( (CIN (ICOM) (1:1) .EQ. 'S') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 's')) THEN
         ICOM = ICOM + 1
         CALL MESAGE (' ')
         CALL MESAGE ('THE CURRENT STATUS OF ALL PLOTTING TOGGLES IS:')
         IF (AXISD) THEN
            CALL MESAGE ('   AXIS PLOTTING                - ON')
         ELSE
            CALL MESAGE ('   AXIS PLOTTING                - OFF')
         ENDIF
         IF (LABP) THEN
            CALL MESAGE ('   LABELING OF POINTS           - ON')
         ELSE
            CALL MESAGE ('   LABELING OF POINTS           - OFF')
         ENDIF
         IF (LABPB) THEN
            CALL MESAGE ('   LABELING OF POINBC FLAGS     - ON')
         ELSE
            CALL MESAGE ('   LABELING OF POINBC FLAGS     - OFF')
         ENDIF
         IF (LABL) THEN
            CALL MESAGE ('   LABELING OF LINES            - ON')
         ELSE
            CALL MESAGE ('   LABELING OF LINES            - OFF')
         ENDIF
         IF (LABI) THEN
            CALL MESAGE ('   LABELING OF LINE INTERVALS   - ON')
         ELSE
            CALL MESAGE ('   LABELING OF LINE INTERVALS   - OFF')
         ENDIF
         IF (LABF) THEN
            CALL MESAGE ('   LABELING OF LINE FACTORS     - ON')
         ELSE
            CALL MESAGE ('   LABELING OF LINE FACTORS     - OFF')
         ENDIF
         IF (LABLB) THEN
            CALL MESAGE ('   LABELING OF NODEBC FLAGS     - ON')
         ELSE
            CALL MESAGE ('   LABELING OF NODEBC FLAGS     - OFF')
         ENDIF
         IF (LABSBD) THEN
            CALL MESAGE ('   LABELING OF ELEMBC FLAGS     - ON')
         ELSE
            CALL MESAGE ('   LABELING OF ELEMBC FLAGS     - OFF')
         ENDIF
         IF (LABR) THEN
            CALL MESAGE ('   LABELING OF REGIONS          - ON')
         ELSE
            CALL MESAGE ('   LABELING OF REGIONS          - OFF')
         ENDIF
         IF (LABMD) THEN
            CALL MESAGE ('   LABELING OF BLOCK ID  (MAT)   - ON')
         ELSE
            CALL MESAGE ('   LABELING OF BLOCK ID  (MAT)   - OFF')
         ENDIF
         IF (LABSC) THEN
            CALL MESAGE ('   LABELING OF REGION SCHEMES   - ON')
         ELSE
            CALL MESAGE ('   LABELING OF REGION SCHEMES   - OFF')
         ENDIF
         IF (LABSZ) THEN
            CALL MESAGE ('   LABELING OF REGION ELEM SIZE - ON')
         ELSE
            CALL MESAGE ('   LABELING OF REGION ELEM SIZE - OFF')
         ENDIF
         IF (FULL) THEN
            CALL MESAGE ('   FULL LABELING OF PROPERTIES  - ON')
         ELSE
            CALL MESAGE ('   FULL LABELING OF PROPERTIES  - OFF')
         ENDIF
         CALL MESAGE ('*-------------------- NOTE -------------------*')
         CALL MESAGE ('    PLOTTING ORDER AT POINTS IS:               ')
         CALL MESAGE ('        POINT NO./POINBC FLAG                  ')
         CALL MESAGE ('    PLOTTING ORDER AT LINE CENTERS IS:         ')
         CALL MESAGE ('       LINE NO./INTERVALS/FACTORS/NODEBC/ELEMBC')
         CALL MESAGE ('    PLOTTING ORDER AT REGION CENTERS IS:       ')
         CALL MESAGE ('        REGION NO./BLOCK ID NO./SCHEME         ')
         CALL MESAGE ('*-------------------- NOTE -------------------*')
C
C  FLAG REGIONS TO BE PROCESSED
C
      ELSEIF ( (CIN (ICOM) (1:1) .EQ. 'R') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'r')) THEN
         ICOM = ICOM + 1
         HXMIN = XMIN
         HYMIN = YMIN
         HXMAX = XMAX
         HYMAX = YMAX
         XMIN = 100000.
         YMIN = 100000.
         XMAX = -100000.
         YMAX = -100000.
         FLAG = .FALSE.
         CALL FLAGD (MP, N (18), LINKP, IPOINT, FLAG)
         CALL FLAGD (ML, N (19), LINKL, ILINE, FLAG)
         CALL FLAGD (MR, N (21), LINKB, IBARST, FLAG)
         CALL FLAGD (MR, N (22), LINKR, IREGN, FLAG)
         GOPLOT = .FALSE.
         CALL MESAGE ('PLOT REGIONS FROM <I1> TO <I2>')
         CALL MESAGE ('HIT RETURN TO END INPUT')
  200    CONTINUE
         IF (ICOM .GT. JCOM) THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN,
     &         IIN, RIN)
            ICOM = 1
         ENDIF
         CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2, IFOUND)

         IF (IFOUND .GT. 0) THEN
            IF (I1 .GT. 0) THEN
               CALL CHECK (I1, I2, N (22))
C
C  FLAG ALL DATA ASSOCIATED WITH THE REGIONS
C
               DO 250 I = I1, I2
                  CALL LTSORT (MR, LINKR, I, II, ADDLNK)
                  IF (II .GT. 0) THEN
C
C  FIND THE MAXIMUM AND MINIMUM
C
                     XMIN = AMIN1 (XMIN, REXTRM (1, II))
                     XMAX = AMAX1 (XMAX, REXTRM (2, II))
                     YMIN = AMIN1 (YMIN, REXTRM (3, II))
                     YMAX = AMAX1 (YMAX, REXTRM (4, II))
                     GOPLOT = .TRUE.
                     IREGN (II) = -IABS (IREGN (II))
                     DO 240 J = IFSIDE (II), IFSIDE (II) + NSPR (II)-1
C
C  FLAG SIDE DATA
C
                        IF ( ISLIST (J) .GT. 0) then
                          CALL LTSORT (MS, LINKS, ISLIST(J), JJ, ADDLNK)
                          if (JJ .GT. 0) THEN
                            DO 220 K = IFLINE (JJ), IFLINE (JJ) +
     &                        NLPS (JJ)-1
                              CALL LTSORT (ML, LINKL, ILLIST (K), KK,
     &                           ADDLNK)
                              IF (KK .GT. 0) THEN
                                 ILINE (KK) = -IABS (ILINE (KK))
                                 DO 210 L = 1, 3
                                    IF (LCON (L, KK) .GT. 0) THEN
                                       CALL LTSORT (MP, LINKP,
     &                                    LCON (L, KK),  LL, ADDLNK)
                                       IF (LL .GT. 0) THEN
                                          IPOINT (LL) =
     &                                       -IABS (IPOINT (LL))
                                       ENDIF
                                    ENDIF
  210                            CONTINUE
                              ENDIF
  220                      CONTINUE
                         end if
C
C  FLAG LINE DATA
C
                        ELSE
                           JJ = IABS (ISLIST (J))
                           CALL LTSORT (ML, LINKL, JJ, KK, ADDLNK)
                           IF (KK .GT. 0) THEN
                              ILINE (KK) = -IABS (ILINE (KK))
                              DO 230 L = 1, 3
                                 IF (LCON (L, KK) .GT. 0) THEN
                                    CALL LTSORT (MP, LINKP,
     &                                 LCON (L, KK),  LL, ADDLNK)
                                    IF (LL .GT. 0) THEN
                                       IPOINT (LL) = -IABS (IPOINT (LL))
                                    ENDIF
                                 ENDIF
  230                         CONTINUE
                           ENDIF
                        ENDIF
  240                CONTINUE
                  ENDIF
  250          CONTINUE
               GOTO 200
            ENDIF
         ENDIF
         IF (GOPLOT) THEN
C
C  PLOT THE REGION DATA THAT HAS BEEN FLAGGED
C
            IF (ALPHA) THEN
               CALL MESAGE ('NO PLOTTING POSSIBLE ON ALPHA-NUMERIC'//
     &            ' TERMINAL')
            ELSE
               CALL PDATA (MP, ML, MR, MSC, IPOINT, COOR, IPBOUN,
     &            ILINE, LTYPE, NINT, LCON, FACTOR, ILBOUN, ISBOUN,
     &            IREGN, IMAT, LINKP, LINKL, LINKR, LINKSC, RSIZE,
     &            SCHEME, DEFSCH, DEFSIZ, REXTRM, N, LABP, LABL, LABR,
     &            FULL, LABMD, LABI, LABF, LABPB, LABLB, LABSBD, LABSC,
     &            LABSZ, AXISD, TITLE, XMIN, XMAX, YMIN, YMAX, XX1,
     &            YY1, XX2, YY2, DEV1, VERSN)
               DRAWN = .TRUE.
            ENDIF
         ELSE
            XMIN = HXMIN
            YMIN = HYMIN
            XMAX = HXMAX
            YMAX = HYMAX
         ENDIF
         XMIN1 = XMIN
         YMIN1 = YMIN
         XMAX1 = XMAX
         YMAX1 = YMAX
C
C  FLAG BARSETS TO BE PLOTTED
C
      ELSEIF ( (CIN (ICOM) (1:1) .EQ. 'B') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'b')) THEN
         ICOM = ICOM + 1
         HXMIN = XMIN
         HYMIN = YMIN
         HXMAX = XMAX
         HYMAX = YMAX
         XMIN = 100000.
         YMIN = 100000.
         XMAX = -100000.
         YMAX = -100000.
         FLAG = .FALSE.
         CALL FLAGD (MP, N (18), LINKP, IPOINT, FLAG)
         CALL FLAGD (ML, N (19), LINKL, ILINE, FLAG)
         CALL FLAGD (MR, N (21), LINKB, IBARST, FLAG)
         CALL FLAGD (MR, N (22), LINKR, IREGN, FLAG)
         GOPLOT = .FALSE.
         CALL MESAGE ('PLOT BARSETS FROM <I1> TO <I2>')
         CALL MESAGE ('HIT RETURN TO END INPUT')
  260    CONTINUE
         IF (ICOM .GT. JCOM) THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN,
     &         IIN, RIN)
            ICOM = 1
         ENDIF
         CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2, IFOUND)
         IF (IFOUND .GT. 0) THEN
            IF (I1 .GT. 0) THEN
               CALL CHECK (I1, I2, N (21))
C
C  FLAG ALL LINES ASSOCIATED WITH THE BARSETS
C
               DO 290 I = I1, I2
                  CALL LTSORT (MS, LINKB, I, II, ADDLNK)
                  IF (II .GT. 0) THEN
                     GOPLOT = .TRUE.
                     IBARST (II) = -IABS (IBARST (II))
                     DO 280 J = JFLINE (II), JFLINE (II)+NLPB (II)-1
                        JJ = IABS (JLLIST (J))
                        CALL LTSORT (ML, LINKL, JJ, KK, ADDLNK)
                        IF (KK .GT. 0) THEN
                           ILINE (KK) = -IABS (ILINE (KK))
                           DO 270 L = 1, 3
                              IF (LCON (L, KK) .GT. 0) THEN
                                 CALL LTSORT (MP, LINKP, LCON (L, KK),
     &                              LL, ADDLNK)
                                 IF (LL .GT. 0) THEN
                                    IPOINT (LL) = -IABS (IPOINT (LL))
                                    XMIN = AMIN1 (XMIN, COOR (1, LL))
                                    XMAX = AMAX1 (XMAX, COOR (1, LL))
                                    YMIN = AMIN1 (YMIN, COOR (2, LL))
                                    YMAX = AMAX1 (YMAX, COOR (2, LL))
                                 ENDIF
                              ENDIF
  270                      CONTINUE
                        ENDIF
  280                CONTINUE
                  ENDIF
  290          CONTINUE
               GOTO 260
            ENDIF
         ENDIF
         IF (GOPLOT) THEN
C
C  PLOT THE BARSET DATA THAT HAS BEEN FLAGGED
C
            IF (ALPHA) THEN
               CALL MESAGE ('NO PLOTTING POSSIBLE ON ALPHANUMERIC '//
     &            'TERMINAL')
            ELSE
               CALL PDATA (MP, ML, MR, MSC, IPOINT, COOR, IPBOUN,
     &            ILINE, LTYPE, NINT, LCON, FACTOR, ILBOUN, ISBOUN,
     &            IREGN, IMAT, LINKP, LINKL, LINKR, LINKSC, RSIZE,
     &            SCHEME, DEFSCH, DEFSIZ, REXTRM, N, LABP, LABL, LABR,
     &            FULL, LABMD, LABI, LABF, LABPB, LABLB, LABSBD, LABSC,
     &            LABSZ, AXISD, TITLE, XMIN, XMAX, YMIN, YMAX, XX1,
     &            YY1, XX2, YY2, DEV1, VERSN)
               DRAWN = .TRUE.
            ENDIF
         ELSE
            XMIN = HXMIN
            YMIN = HYMIN
            XMAX = HXMAX
            YMAX = HYMAX
         ENDIF
         XMIN1 = XMIN
         YMIN1 = YMIN
         XMAX1 = XMAX
         YMAX1 = YMAX
C
C  ENTER ZOOM LOCATION
C
      ELSEIF ( (CIN (ICOM) (1:1) .EQ. 'Z') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'z')) THEN
         ICOM = ICOM + 1
         CALL ZOOMLT (MCOM, ICOM, JCOM, CIN, RIN, IIN, KIN, IDUMP,
     &      DRAWN, ALPHA, DEV1, X1, X2, Y1, Y2, XX1, XX2, YY1, YY2,
     &      XMIN1, XMAX1, YMIN1, YMAX1, XMIN, XMAX, YMIN, YMAX)
         DRAWN = .FALSE.
C
C  RETURN FROM DATA PLOTTING
C
      ELSEIF (CIN (ICOM) (1:1) .EQ. ' ') THEN
         ICOM = ICOM + 1
         FLAG = .FALSE.
         CALL FLAGD (MP, N (18), LINKP, IPOINT, FLAG)
         CALL FLAGD (ML, N (19), LINKL, ILINE, FLAG)
         CALL FLAGD (MS, N (21), LINKB, IBARST, FLAG)
         CALL FLAGD (MR, N (22), LINKR, IREGN, FLAG)
         RETURN
C
C GENERATE A HARDCOPY QMS PLOT
C
      ELSEIF ( ( (CIN (ICOM) (1:1) .EQ. 'H') .OR.
     &   (CIN (ICOM) (1:2) .EQ. 'h')) .AND.
     &   (CIN (ICOM) (2:2).NE.'E') .AND.
     &   (CIN (ICOM) (2:2).NE.'e')) THEN
         ICOM = ICOM + 1
         CALL VDIQES (10002, KAVAL2)
         IF (KAVAL2 .EQ. 1) THEN
            IF (.NOT.ALPHA)CALL VDESCP (10002, 0, 0)
            CALL PDATA (MP, ML, MR, MSC, IPOINT, COOR, IPBOUN, ILINE,
     &         LTYPE, NINT, LCON, FACTOR, ILBOUN, ISBOUN, IREGN, IMAT,
     &         LINKP, LINKL, LINKR, LINKSC, RSIZE, SCHEME, DEFSCH,
     &         DEFSIZ, REXTRM, N, LABP, LABL, LABR, FULL, LABMD, LABI,
     &         LABF, LABPB, LABLB, LABSBD, LABSC, LABSZ, AXISD, TITLE,
     &         XMIN, XMAX, YMIN, YMAX, XX1, YY1, XX2, YY2, 'XXX', VERSN)
            IF (.NOT.ALPHA)CALL VDESCP (10001, 0, 0)
            CALL MESAGE ('HARDCOPY PLOT GENERATED')
            HARDPL = .TRUE.
         ELSE
            CALL MESAGE ('HARDCOPY DEVICE NOT AVAILABLE')
         ENDIF
C
C  PLOT THE CURRENT ACTIVE ITEMS
C
      ELSEIF ( (CIN (ICOM) (1:1) .EQ. 'P') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'p')) THEN
         ICOM = ICOM + 1
         IF (ALPHA) THEN
            CALL MESAGE ('NO PLOTTING POSSIBLE ON ALPHANUMERIC '//
     &         'TERMINAL')
         ELSE
            CALL PDATA (MP, ML, MR, MSC, IPOINT, COOR, IPBOUN, ILINE,
     &         LTYPE, NINT, LCON, FACTOR, ILBOUN, ISBOUN, IREGN, IMAT,
     &         LINKP, LINKL, LINKR, LINKSC, RSIZE, SCHEME, DEFSCH,
     &         DEFSIZ, REXTRM, N, LABP, LABL, LABR, FULL, LABMD, LABI,
     &         LABF, LABPB, LABLB, LABSBD, LABSC, LABSZ, AXISD, TITLE,
     &         XMIN, XMAX, YMIN, YMAX, XX1, YY1, XX2, YY2, DEV1, VERSN)
            DRAWN = .TRUE.
         ENDIF
C
C  EXIT OPTION - EXITS FASTQ
C
      ELSEIF ( (CIN (ICOM) (1:2) .EQ. 'EX') .OR.
     &   (CIN (ICOM) (1:2) .EQ. 'ex')) THEN
         ICOM = ICOM + 1
         CALL STRLNG (CIN(ICOM), LEN)
         IF (((LEN .GT. 1) .AND. (CIN(ICOM)(2:2) .NE. 'X')) .OR.
     &      ((LEN .GT. 1) .AND. (CIN(ICOM)(2:2) .NE. 'x'))) THEN
            CALL HELP_FQ(5)
         ELSE
            CALL FEXIT (WROTE, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN,
     &         TIME1, BATCH, VERSN)
         ENDIF
         GOTO 100
C
C  WRITE OUT THE HELP MESSAGE
C
      ELSE
         ICOM = ICOM + 1
         CALL HELP_FQ (5)
      ENDIF
      GOTO 100
C
      END
