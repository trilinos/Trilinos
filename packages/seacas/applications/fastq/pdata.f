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

C $Id: pdata.f,v 1.5 2007/07/24 13:10:18 gdsjaar Exp $
C $Log: pdata.f,v $
C Revision 1.5  2007/07/24 13:10:18  gdsjaar
C Fix problem with boundary condition memory overwrite.
C
C Remove old ls5 and r25 terminal tests
C
C Revision 1.4  1999/06/21 22:43:40  gdsjaar
C Fixed more uninitialized variables; one was causing core dump on g77
C compiled executable.
C
C VERSN was not consistently defined -- now 10 characters everywhere
C
C Updated so full version string output
C
C Added capability to debug memory using unit specified in EXT99
C variable. Similar to STRTUP in SUPLIB
C
C Cleaned up some other code
C
C Upped version
C
C Revision 1.3  1998/07/14 18:19:31  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.2  1992/12/08 22:13:54  gdsjaar
C Changed color of point label output from yellow to red
C
c Revision 1.1.1.1  1990/11/30  11:13:14  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:13:13  gdsjaar
c Initial revision
c
C
CC* FILE: [.MAIN]PDATA.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE PDATA (MP, ML, MR, MSC, IPOINT, COOR, IPBOUN, ILINE,
     &   LTYPE, NINT, LCON, FACTOR, ILBOUN, ISBOUN, IREGN, IMAT, LINKP,
     &   LINKL, LINKR, LINKSC, RSIZE, SCHEME, DEFSCH, DEFSIZ, REXTRM,
     &   N, LABP, LABL, LABR, FULL, LABMD, LABI, LABF, LABPB, LABLB,
     &   LABSBD, LABSC, LABSZ, AXISD, TITLE, XMIN, XMAX, YMIN, YMAX,
     &   XX1, YY1, XX2, YY2, DEV1, VERSN)
C***********************************************************************
C
C  SUBROUTINE PDATA = PLOTS FLAGGED POINTS,  LINES,  AND REGIONS
C
C***********************************************************************
C
      DIMENSION IPOINT (MP), COOR (2, MP), IPBOUN (MP)
      DIMENSION ILINE (ML), LTYPE (ML), NINT (ML), LCON (3, ML)
      DIMENSION FACTOR (ML)
      DIMENSION ILBOUN (ML), ISBOUN (ML)
      DIMENSION IREGN (MR), IMAT (MR), REXTRM (4, MR), RSIZE (MR)
      DIMENSION SCHEME (MSC)
      DIMENSION LINKP (2, MP), LINKL (2, ML), LINKR (2, MR)
      DIMENSION LINKSC (2, MR)
      DIMENSION N (29), XDUM (2), YDUM (2)
C
      CHARACTER*72 DUMMY, SCHEME, DEFSCH, TITLE, DEV1*3, CDUMMY*4
      CHARACTER*8 DATE, TIME, VERSN*10
C
      LOGICAL LABP, LABL, LABR, AXISD, LABMD, LABI, LABF
      LOGICAL LABPB, LABLB, LABSBD
      LOGICAL ADDLNK, CPUIFC, TEST, FULL, LABSC
      LOGICAL GETMAX, ADD, LABSZ
C
C  INITIALIZE THE PLOTTING SURFACE
C
      TEST = .FALSE.
      GETMAX = .FALSE.
      IF (TEST)OPEN (UNIT = 12, FILE = 'HP7580.DAT', STATUS = 'NEW')
      ADDLNK = .FALSE.
      CALL PLTBGN
      CALL PLTSTV (2, 160.)
      XDIMR = ABS (XMAX - XMIN)
      YDIMR = ABS (YMAX - YMIN)
      XDIMD = 1.
      YDIMD = .75
      CALL MPVIEW (0., XDIMD, 0., YDIMD)
      XRAT = XDIMR/XDIMD
      YRAT = YDIMR/YDIMD
      IF (XRAT.LT.YRAT) THEN
         XDIMR = XDIMD*YRAT
         XX1 =  (XMIN + XMAX - XDIMR)*.5
         XX2 =  (XMIN + XMAX + XDIMR)*.5
         XDIMR = XX2 - XX1
         YY1 = YMIN
         YY2 = YMAX
      ELSE
         YDIMR = YDIMD*XRAT
         YY1 =  (YMIN + YMAX - YDIMR)*.5
         YY2 =  (YMIN + YMAX + YDIMR)*.5
         YDIMR = YY2 - YY1
         XX1 = XMIN
         XX2 = XMAX
      ENDIF
C
C  SET UP SCALING EXTREMES FOR AXIS
C
      IF (TEST) THEN
         WRITE (12, 10000)'IN;SP6;;IP - 5710, -10060, 15710, 10060;'
         WRITE (12, 10010)
     &      'SC', IFIX (XX1*1000), ', ', IFIX (YY1*1000), ', ',
     &      IFIX (XX2*1000), ', ', IFIX (YY2*1000), ';'
      ENDIF
      IF (AXISD) THEN
         XDUM (1) = XX1 -  (XDIMR*.05)
         XDUM (2) = XX2 +  (XDIMR*.05)
         YDUM (1) = YY1 -  (YDIMR*.05)
         YDUM (2) = YY2 +  (YDIMR*.05)
         SHRINK = .2
      ELSE
         SHRINK = .1
      ENDIF
C
C  SHRINK TO FIT A BORDER ON THE PLOT
C
      XX1 = XX1 -  (XDIMR*SHRINK)
      XX2 = XX2 +  (XDIMR*SHRINK)
      YY1 = YY1 -  (YDIMR*SHRINK)
      YY2 = YY2 +  (YDIMR*SHRINK)
      CALL MPORT2 (XX1, XX2, YY1, YY2)
      CALL PLTFRM (0)
C
C  PLOT THE TITLE AND THE TRACE
C
      CALL STRLNG (TITLE, LEN)
      IF ( (LEN.GT.1) .OR. (TITLE (1:1).NE.' ')) THEN
         CALL PLTXHL (TITLE (1:LEN), XLEN)
         XBEGIN = AMAX1 (0.,  (XDIMD*.5 - XLEN*.5))
         CALL PLTXTH (XBEGIN, YDIMD*.95, TITLE (1:LEN))
      ENDIF
      DUMMY(1:10) = ' DRAWN BY '
      DUMMY(11:20) = VERSN
      DUMMY(21:22) = '  '
      CALL EXDATE (DATE)
      DUMMY(23:30) = DATE
      DUMMY(31:32) = '  '
      CALL EXTIME (TIME)
      DUMMY(33:40) = TIME
      CALL PLTXTH (0., 0., DUMMY(1:40))
C
C  DRAW THE AXIS IF REQUIRED,  AND SET CLIPPING WITHIN AXIS
C
      IF (AXISD)CALL SETAXS (XDUM, YDUM)
      IF (CPUIFC (.TRUE.))GOTO 130
C
C  PLOT THE POINTS FLAGGED
C
      IF ( (LABP) .OR. (LABPB)) THEN
         DO 100 I = 1, N (18)
            IF (CPUIFC (.TRUE.))GOTO 130
            CALL LTSORT (MP, LINKP, I, II, ADDLNK)
            IF (II.GT.0) THEN
               IF (IPOINT (II).LT.0) THEN
                  INUM =  - IPOINT (II)
                  CALL MP2PT (1, COOR (1, II), COOR (2, II),
     &               X1, Y1, MASK)
                  IF (MOD (MASK, 2).NE.0) THEN
C
C  PLOT THE POINT LABELS
C
                     IF (LABP) THEN
                        CALL PLTSTD (1, 1.)
                        CALL GETDUM (INUM, DUMMY, LEN)
                        CALL PLTXTH (X1, Y1, DUMMY (1:LEN))
                        CALL PLTXHE (X1, Y1)
                     ENDIF
C
C  PLOT THE POINBC FLAGS
C
                     IF ( ( (LABPB) .OR. ( (FULL) .AND. (LABP))) .AND.
     &                  (IPBOUN (II).GT.0)) THEN
                        CALL PLTSTD (1, 5.)
                        IF (LABP) THEN
                           CALL PLTXTH (X1, Y1, '/')
                           CALL PLTXHE (X1, Y1)
                        ENDIF
                        CALL GETDUM (IPBOUN (II), DUMMY, LEN)
                        CALL PLTXTH (X1, Y1, DUMMY (1:LEN))
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
  100    CONTINUE
      ENDIF
C
C  PLOT ALL LINES THAT HAVE BEEN FLAGGED
C
      DO 110 I = 1, N (19)
         IF (CPUIFC (.TRUE.))GOTO 130
         CALL LTSORT (ML, LINKL, I, II, ADDLNK)
         IF (II.GT.0) THEN
            IF (LABL) THEN
               ADD = .TRUE.
            ELSE
               ADD = .FALSE.
            ENDIF
            CALL PLTSTD (1, 7.)
            IF (ILINE (II).LT.0) THEN
               IF ( (LABL) .OR. (LABLB) .OR. (LABSBD) .OR. (LABF)
     &            .OR. (LABI)) THEN
                  KNUM =  - ILINE (II)
               ELSE
                  KNUM = 0
               ENDIF
               LT = LTYPE (II)
               IP1 = LCON (1, II)
               IP2 = LCON (2, II)
               IP3 = LCON (3, II)
               CALL LTSORT (MP, LINKP, IP1, IPNTR1, ADDLNK)
               CALL LTSORT (MP, LINKP, IP2, IPNTR2, ADDLNK)
               IF (IP3.NE.0) THEN
                  CALL LTSORT (MP, LINKP, IABS (IP3), IPNTR3, ADDLNK)
               ELSE
                  IPNTR3 = 0
               ENDIF
               IF ((IPNTR1.GT.0) .AND. (IPNTR2.GT.0) .AND.
     &            ((LT.EQ.1) .OR. (IPNTR3.GT.0)) ) THEN
                  CALL DLINE (MP, ML, COOR, LINKP, KNUM, LT, IP1, IP2,
     &               IP3, LABL, X1, Y1, TEST, GETMAX, DUM1, DUM2, DUM3,
     &               DUM4)
C
C  PLOT INTERVAL NUMBERS
C
                  IF ( ( (FULL) .AND. (LABL)) .OR. (LABI)) THEN
                     CALL PLTSTD (1, 5.)
                     IF (ADD) THEN
                        CALL PLTXHE (X1, Y1)
                        CALL PLTXTH (X1, Y1, '/')
                        CALL PLTXHE (X1, Y1)
                     ENDIF
                     CALL GETDUM (NINT (II), DUMMY, LEN)
                     IF (TEST) THEN
                        CALL PLTD2G (X1, Y1, XR, YR)
                        CALL PLTG2D (X1, Y1, XR, YR)
                        WRITE (12, 10020)'PU;PA', IFIX (XR*1000.),
     &                     ', ', IFIX (YR*1000.), ';LB',
     &                     DUMMY (1:LEN), CHAR (3)
                     ENDIF
                     CALL PLTXTH (X1, Y1, DUMMY (1:LEN))
                     ADD = .TRUE.
                  ENDIF
C
C  PLOT THE LINE FACTOR
C
                  IF ( ( (FULL) .AND. (LABL)) .OR. (LABF)) THEN
                     IF (ADD) THEN
                        CALL PLTSTD (1, 1.)
                        CALL PLTXHE (X1, Y1)
                        CALL PLTXTH (X1, Y1, '/')
                        CALL PLTXHE (X1, Y1)
                     ENDIF
                     CALL GTXDUM (FACTOR (II), DUMMY, LEN)
                     CALL PLTXTH (X1, Y1, DUMMY (1:LEN))
                     ADD = .TRUE.
                  ENDIF
C
C  PLOT THE LINEBC FLAGS
C
                  IF ( ( ( (FULL) .AND. (LABL)) .OR. (LABLB)) .AND.
     &               (ILBOUN (II).GT.0)) THEN
                     CALL PLTSTD (1, 2.)
                     IF (ADD) THEN
                        CALL PLTXHE (X1, Y1)
                        CALL PLTXTH (X1, Y1, '/')
                        CALL PLTXHE (X1, Y1)
                     ENDIF
                     CALL GETDUM (ILBOUN (II), DUMMY, LEN)
                     CALL PLTXTH (X1, Y1, DUMMY (1:LEN))
                     ADD = .TRUE.
                  ENDIF
C
C  PLOT THE SIDEBC FLAGS
C
                  IF ( ( ( (FULL) .AND. (LABL)) .OR. (LABSBD)) .AND.
     &               (ISBOUN (II).GT.0)) THEN
                     CALL PLTSTD (1, 3.)
                     IF (ADD) THEN
                        CALL PLTXHE (X1, Y1)
                        CALL PLTXTH (X1, Y1, '/')
                        CALL PLTXHE (X1, Y1)
                     ENDIF
                     CALL GETDUM (ISBOUN (II), DUMMY, LEN)
                     CALL PLTXTH (X1, Y1, DUMMY (1:LEN))
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
  110 CONTINUE
C
C  PLOT ALL REGIONS FLAGGED
C
      IF ( (LABR) .OR. (LABMD) .OR. (LABSC) .OR. (LABSZ)) THEN
         IF (CPUIFC (.TRUE.))GOTO 130
         DO 120 I = 1, N (22)
            CALL LTSORT (MR, LINKR, I, II, ADDLNK)
            IF (II.GT.0) THEN
               IF (IREGN (II).LT.0) THEN
                  ADD = .FALSE.
                  INUM =  - IREGN (II)
                  XMID =  (REXTRM (1, II) + REXTRM (2, II))/2.
                  YMID =  (REXTRM (3, II) + REXTRM (4, II))/2.
                  CALL MP2PT (1, XMID, YMID, X1, Y1, MASK)
                  IF ( (MOD (MASK, 2).NE.0)) THEN
C
C  PLOT THE REGION NUMBER
C
                     IF (LABR) THEN
                        CALL PLTSTD (1, 2.)
                        CALL GETDUM (INUM, DUMMY, LEN)
                        CALL PLTXTH (X1, Y1, DUMMY (1:LEN))
                        ADD = .TRUE.
                     ENDIF
C
C  PLOT OUT THE MATERIAL NUMBER
C
                     IF (((FULL) .AND. (LABR)) .OR. (LABMD)) THEN
                        CALL PLTSTD (1, 1.)
                        IF (ADD) THEN
                           CALL PLTXHE (X1, Y1)
                           CALL PLTXTH (X1, Y1, '/')
                           CALL PLTXHE (X1, Y1)
                        ENDIF
                        ADD = .TRUE.
                        CALL GETDUM (IMAT (II), DUMMY, LEN)
                        CALL PLTXTH (X1, Y1, DUMMY (1:LEN))
                     ENDIF
C
C  PLOT OUT THE SIZE NUMBER FOR THE REGION
C
                     IF (((FULL) .AND. (LABR)) .OR. (LABSZ)) THEN
                        CALL PLTSTD (1, 1.)
                        IF (ADD) THEN
                           CALL PLTXHE (X1, Y1)
                           CALL PLTXTH (X1, Y1, '/')
                           CALL PLTXHE (X1, Y1)
                        ENDIF
                        ADD = .TRUE.
                        CALL GTXDUM (RSIZE (II), DUMMY, LEN)
                        CALL PLTXTH (X1, Y1, DUMMY (1:LEN))
                     ENDIF
C
C  PLOT OUT THE SCHEME
C
                     IF (((FULL) .AND. (LABR)) .OR. (LABSC)) THEN
                        CALL PLTSTD (1, 7.)
                        IF (ADD) THEN
                           CALL PLTXHE (X1, Y1)
                           CALL PLTXTH (X1, Y1, '/')
                           CALL PLTXHE (X1, Y1)
                        ENDIF
                        CALL LTSORT (MR, LINKSC, INUM, IPNTR, ADDLNK)
                        IF ( (INUM.LE.N (24)) .AND. (IPNTR.GT.0)) THEN
                           CALL STRLNG (SCHEME (IPNTR), LEN)
                           IF (TEST) THEN
                              CALL PLTD2G (X1, Y1, XR, YR)
                              WRITE (12, 10020)'PU;PA', IFIX (XR*1000.),
     &                           ', ', IFIX (YR*1000.), ';LB',
     &                           SCHEME (IPNTR) (1:LEN), CHAR (3)
                           ENDIF
                           CALL PLTXTH (X1, Y1, SCHEME (IPNTR) (1:LEN))
                        ELSE
                           CALL STRLNG (DEFSCH, LEN)
                           IF (TEST) THEN
                              CALL PLTD2G (X1, Y1, XR, YR)
                              WRITE (12, 10020)'PU;PA', IFIX (XR*1000.),
     &                           ', ', IFIX (YR*1000.), ';LB',
     &                           DEFSCH (1:LEN), CHAR (3)
                           ENDIF
                           CALL PLTXTH (X1, Y1, DEFSCH (1:LEN))
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
  120    CONTINUE
      ENDIF
  130 CONTINUE
      CALL PLTSTD (1, 7.)
      CALL PLTBEL
      CALL PLTFLU
C
      RETURN
C
10000 FORMAT (A)
10010 FORMAT (A2, I10, A1, I10, A1, I10, A1, I10, A1)
10020 FORMAT (A5, I10, A1, I10, A3, A, A1)
C
      END
