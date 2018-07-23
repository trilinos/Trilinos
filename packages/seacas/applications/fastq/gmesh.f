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

C $Id: gmesh.f,v 1.2 1998/07/14 18:19:07 gdsjaar Exp $
C $Log: gmesh.f,v $
C Revision 1.2  1998/07/14 18:19:07  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:09:00  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:08:57  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]GMESH.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GMESH (NPNODE, NPELEM, MXNFLG, MXSFLG, NPNBC, NPSBC,
     &   MAXKXN, MR, NPREGN, MCOM, ICOM, JCOM, CIN, RIN, IIN, KIN, NNN,
     &   KKK, NUMMAT, NNXK, IPART, NODES, NNFLG, NNPTR, NSFLG, NVPTR,
     &   NVLEN, NSIDEN, MAPDXG, XN, YN, NXK, MAT, ILOOK, MAPGXD, CENTK,
     &   MATMAP, WTNODE, WTSIDE, NBCNOD, NNLIST, NBCSID, NVLIST, TITLE,
     &   IDUMP, AXIS, AREACG, LABE, LABO, LABN, LABNB, LABSB, LABM,
     &   LABW, IDEV, ALPHA, DEV1, EIGHT, NINE, VAXVMS, VERSN, WROTE,
     &   TIME1, HARDPL, BATCH)
C***********************************************************************
C
C  SUBROUTINE GMESH = SETS UP GRAPHICS FOR THE GENERATED MESH
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     MESH = ALLOWS MESH GENERATION AND DISPLAY
C
C***********************************************************************
C
C  VARIABLES USED:
C     TITLE  = MESH TITLE
C     LABE   = .TRUE. IF ELEMENT NUMBERS ARE TO BE PLOTTED
C     LABN   = .TRUE. IF NODE NUMBERS ARE TO BE PLOTTED
C     LABNB  = .TRUE. IF NODE BOUNDARY NUMBERS ARE TO BE PLOTTED
C     LABSB  = .TRUE. IF ELEMENT BOUNDARY NUMBERS ARE TO BE PLOTTED
C     LABM   = .TRUE. IF MATERIAL NUMBERS ARE TO BE PLOTTED
C     AXIS   = .TRUE. IF THE AXIS IS TO BE DRAWN
C     AREACG = .TRUE. IF THE AREA AND C.G. ARE CALCULATED AND THE C.G.
C              IS DISPLAYED
C
C***********************************************************************
C
      DIMENSION IPART (3, NPREGN), CENTK (2, NPELEM)
      DIMENSION ILOOK (NNXK * MAXKXN)
      DIMENSION XN (NPNODE), YN (NPNODE), NXK (NNXK, NPELEM)
      DIMENSION MAT (NPELEM)
      DIMENSION MAPGXD (NPNODE)
      DIMENSION NODES (NPNBC), NSIDEN (NPSBC), WTNODE (NPNBC)
      DIMENSION WTSIDE (NPSBC)
      DIMENSION NNFLG (MXNFLG), NNPTR (MXNFLG)
      DIMENSION NSFLG (MXSFLG), NVPTR (MXSFLG), NVLEN (MXSFLG)
      DIMENSION MAPDXG (NPNODE), MATMAP (3, NPREGN)
      DIMENSION KIN (MCOM), IIN (MCOM), RIN (MCOM)
      DIMENSION IDEV (2)
C
      LOGICAL LABE, LABN, LABNB, LABSB, GOPLOT, SETFLG, DRAWN
      LOGICAL AXIS, LABM, AREACG
      LOGICAL LABW, ALPHA, EIGHT, NINE, OLD, LABO, VAXVMS, REGPLT
      LOGICAL WROTE, HARDPL
      LOGICAL BATCH
C
      CHARACTER*72 TITLE, CIN (MCOM)
      CHARACTER*3 DEV1, VERSN*9
C
      IZ = 0
C
C  CALCULATE THE CENTER OF EACH ELEMENT FOR CLIPPING CONSIDERATIONS
C
      DO 100 I = 1, KKK
         IF (NXK (3, I) .EQ. 0) THEN
            CENTK (1, I) = .5 * (XN (NXK (1, I)) + XN (NXK (2, I)))
            CENTK (2, I) = .5 * (YN (NXK (1, I)) + YN (NXK (2, I)))
         ELSE IF (NXK (4, I) .EQ. 0) THEN
            CENTK (1, I) = .33333 * (XN (NXK (1, I)) +
     &         XN (NXK (2, I)) + XN (NXK (3, I)))
            CENTK (2, I) = .33333 * (YN (NXK (1, I)) +
     &         YN (NXK (2, I)) + YN (NXK (3, I)))
         ELSE IF ((EIGHT).OR. (NINE)) THEN
            CENTK (1, I) = .5 * (XN (NXK (1, I)) + XN (NXK (5, I)))
            CENTK (2, I) = .5 * (YN (NXK (1, I)) + YN (NXK (5, I)))
         ELSE
            CENTK (1, I) = .5 * (XN (NXK (1, I)) + XN (NXK (3, I)))
            CENTK (2, I) = .5 * (YN (NXK (1, I)) + YN (NXK (3, I)))
         END IF
  100 CONTINUE
      DRAWN = .FALSE.
C
C  FIND THE BODY MIN AND MAX
C
      CALL MINMAX_FQ (NPNODE, NNN, XN, YN, XMIN, XMAX, YMIN, YMAX)
      XMIN1 = XMIN
      XMAX1 = XMAX
      YMIN1 = YMIN
      YMAX1 = YMAX
C
C  ENTER GRAPHICS OPTION
C
  110 CONTINUE
      IF (ICOM .GT. JCOM) THEN
         CALL MESAGE (' ')
         CALL FREFLD (IZ, IZ, 'ENTER MESH GRAPHICS OPTION: ', MCOM,
     &      IOSTAT, JCOM, KIN, CIN, IIN, RIN)
         ICOM = 1
      END IF
C
C  TOGGLE DRAWING OF THE AXIS
C
      IF ((CIN (ICOM) (1:1) .EQ. 'A')
     &   .OR. (CIN (ICOM) (1:1) .EQ. 'a')) THEN
         ICOM = ICOM + 1
         IF (AXIS) THEN
            AXIS = .FALSE.
            CALL MESAGE ('AXIS DRAWING - OFF')
         ELSE
            AXIS = .TRUE.
            CALL MESAGE ('AXIS DRAWING - ON')
         END IF
C
C  TOGGLE CALCULATION OF AREA AND C.G.
C
      ELSE IF ((CIN (ICOM) (1:1) .EQ. 'C')
     &   .OR. (CIN (ICOM) (1:1) .EQ. 'c')) THEN
         ICOM = ICOM + 1
         IF (AREACG) THEN
            AREACG = .FALSE.
            CALL MESAGE ('AREA AND C.G. REPORT - OFF')
         ELSE
            AREACG = .TRUE.
            CALL MESAGE ('AREA AND C.G. REPORT - ON')
         END IF
C
C  TOGGLE NODAL BOUNDARY DISPLAY
C
      ELSE IF ((CIN (ICOM) (1:2) .EQ. 'NB') .OR.
     &   (CIN (ICOM) (1:2) .EQ. 'nb')) THEN
         ICOM = ICOM + 1
         IF (LABNB) THEN
            LABNB = .FALSE.
            CALL MESAGE ('NODAL BOUNDARY DISPLAY - OFF')
         ELSE
            LABNB = .TRUE.
            CALL MESAGE ('NODAL BOUNDARY DISPLAY - ON')
         END IF
C
C  TOGGLE ELEMENT SIDE BOUNDARY DISPLAY
C
      ELSE IF ((CIN (ICOM) (1:2) .EQ. 'EB') .OR.
     &   (CIN (ICOM) (1:2) .EQ. 'eb')) THEN
         ICOM = ICOM + 1
         IF (LABSB) THEN
            LABSB = .FALSE.
            CALL MESAGE ('ELEMENT SIDE BOUNDARY DISPLAY - OFF')
         ELSE
            LABSB = .TRUE.
            CALL MESAGE ('ELEMENT SIDE BOUNDARY DISPLAY - ON')
         END IF
C
C  TOGGLE WEIGHTING FACTOR DISPLAY
C
      ELSE IF ((CIN (ICOM) (1:1) .EQ. 'W') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'w')) THEN
         ICOM = ICOM + 1
         IF (LABW) THEN
            LABW = .FALSE.
            CALL MESAGE ('BOUNDARY WEIGHTING DISPLAY - OFF')
         ELSE
            LABW = .TRUE.
            CALL MESAGE ('BOUNDARY WEIGHTING DISPLAY - ON')
         END IF
C
C  TOGGLE ELEMENT NUMBERING
C
      ELSE IF ((CIN (ICOM) (1:2) .EQ. 'EN') .OR.
     &   (CIN (ICOM) (1:2) .EQ. 'en')) THEN
         ICOM = ICOM + 1
         IF (LABE) THEN
            LABE = .FALSE.
            CALL MESAGE ('ELEMENT NUMBERS - OFF')
         ELSE
            LABE = .TRUE.
            CALL MESAGE ('ELEMENT NUMBERS - ON')
            LABO = .FALSE.
         END IF
C
C  TOGGLE NODE NUMBERING
C
      ELSE IF ((CIN (ICOM) (1:1) .EQ. 'N') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'n')) THEN
         ICOM = ICOM + 1
         IF (LABN) THEN
            LABN = .FALSE.
            CALL MESAGE ('NODE NUMBERS - OFF')
         ELSE
            LABN = .TRUE.
            CALL MESAGE ('NODE NUMBERS - ON')
         END IF
C
C  TOGGLE MATERIAL NUMBER DISPLAY
C
      ELSE IF ((CIN (ICOM) (1:2) .EQ. 'MN') .OR.
     &   (CIN (ICOM) (1:2) .EQ. 'mn')) THEN
         ICOM = ICOM + 1
         IF (LABM) THEN
            LABM = .FALSE.
            CALL MESAGE ('MATERIAL NUMBERING - OFF')
         ELSE
            LABM = .TRUE.
            CALL MESAGE ('MATERIAL NUMBERING - ON')
         END IF
C
C  TOGGLE OPTIMIZED ORDER DISPLAY
C
      ELSE IF ((CIN (ICOM) (1:1) .EQ. 'O') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'o')) THEN
         ICOM = ICOM + 1
         IF (LABO) THEN
            LABO = .FALSE.
            CALL MESAGE ('OPTIMIZED ORDER NUMBERING - OFF')
         ELSE
            LABO = .TRUE.
            CALL MESAGE ('OPTIMIZER ORDER NUMBERING - ON')
            LABE = .FALSE.
         END IF
C
C  PLOT ALL ACTIVE ELEMENTS  (ZOOM STILL APPLIES)
C
      ELSE IF ((CIN (ICOM) (1:1) .EQ. 'P') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'p')) THEN
         ICOM = ICOM + 1
         IF (ALPHA) THEN
            CALL MESAGE ('NO PLOTTING POSSIBLE ON ALPHA-NUMERIC '//
     &         'TERMINAL')
         ELSE
            CALL PMESH (NPNODE, NPELEM, MXNFLG, MXSFLG, NPNBC, NPSBC,
     &         MAXKXN, NNN, KKK, NNXK, NBCNOD, NBCSID, NNLIST, NVLIST,
     &         NODES, NSIDEN, NNFLG, NNPTR, NSFLG, NVPTR, NVLEN, XN, YN,
     &         NXK, MAT, MAPDXG, MAPGXD, WTNODE, WTSIDE, AXIS, AREACG,
     &         LABE, LABO, LABN, LABNB, LABSB, LABM, LABW, CENTK, ILOOK,
     &         XMIN, XMAX, YMIN, YMAX, XX1, XX2, YY1, YY2, TITLE, DEV1,
     &         EIGHT, NINE, VERSN, VAXVMS)
            DRAWN = .TRUE.
         END IF
C
C  PLOT A LIMITED NUMBER OF ELEMENTS BY ELEMENT NUMBER
C
      ELSE IF ((CIN (ICOM) (1:2) .EQ. 'EP') .OR.
     &   (CIN (ICOM) (1:2) .EQ. 'ep')) THEN
         ICOM = ICOM + 1
         SETFLG = .FALSE.
         IF (LABO) THEN
            OLD = .FALSE.
         ELSE
            OLD = .TRUE.
         END IF
         CALL FLAGK (NPELEM, NNXK, NXK, MAPDXG, 1, KKK, SETFLG, OLD)
         SETFLG = .TRUE.
         GOPLOT = .FALSE.
         CALL MESAGE ('PLOT ELEMENTS FROM <I1> TO <I2>')
         CALL MESAGE ('HIT RETURN TO END INPUT')
  120    CONTINUE
         IF (ICOM .GT. JCOM) THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN, IIN,
     &         RIN)
            ICOM = 1
         END IF
         CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2, IFOUND)
         IF (IFOUND .GT. 0) THEN
            IF (I1 .GT. 0) THEN
               CALL CHECK (I1, I2, KKK)
               CALL FLAGK (NPELEM, NNXK, NXK, MAPDXG, I1, I2, SETFLG,
     &            OLD)
               GOPLOT = .TRUE.
               GO TO 120
            END IF
         END IF
         IF (GOPLOT) THEN
            IF (ALPHA) THEN
               CALL MESAGE ('NO PLOTTING POSSIBLE ON ALPHA-NUMERIC '//
     &            ' TERMINAL')
            ELSE
               CALL MNMXK (NPELEM, NPNODE, NNXK, NXK, XN, YN, CENTK,
     &            KKK, XMIN, XMAX, YMIN, YMAX)
               CALL PMESH (NPNODE, NPELEM, MXNFLG, MXSFLG, NPNBC, NPSBC,
     &            MAXKXN, NNN, KKK, NNXK, NBCNOD, NBCSID, NNLIST,
     &            NVLIST, NODES, NSIDEN, NNFLG, NNPTR, NSFLG, NVPTR,
     &            NVLEN, XN, YN, NXK, MAT, MAPDXG, MAPGXD, WTNODE,
     &            WTSIDE, AXIS, AREACG, LABE, LABO, LABN, LABNB, LABSB,
     &            LABM, LABW, CENTK, ILOOK, XMIN, XMAX, YMIN, YMAX, XX1,
     &            XX2, YY1, YY2, TITLE, DEV1, EIGHT, NINE, VERSN,
     &            VAXVMS)
               DRAWN = .TRUE.
            END IF
         END IF
C
C  PLOT ELEMENTS BY REGION (S) OR BARSET (S) CHOSEN
C
      ELSE IF ((CIN (ICOM) (1:1) .EQ. 'R') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'r') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'B') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'b')) THEN
         SETFLG = .FALSE.
         OLD = .TRUE.
         CALL FLAGK (NPELEM, NNXK, NXK, MAPDXG, 1, KKK, SETFLG, OLD)
         SETFLG = .TRUE.
         GOPLOT = .FALSE.
         IF ((CIN (ICOM) (1:1) .EQ. 'R') .OR.
     &      (CIN (ICOM) (1:1) .EQ. 'r')) THEN
            CALL MESAGE ('PLOT REGIONS FROM <I1> TO <I2>')
            REGPLT = .TRUE.
         ELSE
            CALL MESAGE ('PLOT BARSETS FROM <I1> TO <I2>')
            REGPLT = .FALSE.
         END IF
         ICOM = ICOM + 1
         CALL MESAGE ('HIT RETURN TO END INPUT')
  130    CONTINUE
         IF (ICOM .GT. JCOM) THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN, IIN,
     &         RIN)
            ICOM = 1
         END IF
         CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2, IFOUND)
         IF (IFOUND .GT. 0) THEN
            IF (I1 .GT. 0) THEN
               IF (I2 .LT. I1)I2 = I1
               DO 140 I = 1, NPREGN
                  IF ((IPART (1, I) .GE. I1)
     &               .AND. (IPART (1, I) .LE. I2)) THEN
                     IF ((REGPLT .AND. (NXK (4, IPART (2, I)) .GT. 0))
     &                  .OR. (NXK (4, IPART (2, I)) .EQ. 0)) THEN
                        CALL FLAGK (NPELEM, NNXK, NXK, MAPDXG,
     &                     IPART (2, I), IPART (3, I), SETFLG, OLD)
                        GOPLOT = .TRUE.
                     END IF
                  END IF
  140          CONTINUE
               GO TO 130
            END IF
         END IF
         IF (GOPLOT) THEN
            IF (ALPHA) THEN
               CALL MESAGE ('NO PLOTTING POSSIBLE ON ALPHA-NUMERIC '//
     &            'TERMINAL')
            ELSE
               CALL MNMXK (NPELEM, NPNODE, NNXK, NXK, XN, YN, CENTK,
     &            KKK, XMIN, XMAX, YMIN, YMAX)
               CALL PMESH (NPNODE, NPELEM, MXNFLG, MXSFLG, NPNBC, NPSBC,
     &            MAXKXN, NNN, KKK, NNXK, NBCNOD, NBCSID, NNLIST,
     &            NVLIST, NODES, NSIDEN, NNFLG, NNPTR, NSFLG, NVPTR,
     &            NVLEN, XN, YN, NXK, MAT, MAPDXG, MAPGXD, WTNODE,
     &            WTSIDE, AXIS, AREACG, LABE, LABO, LABN, LABNB, LABSB,
     &            LABM, LABW, CENTK, ILOOK, XMIN, XMAX, YMIN, YMAX, XX1,
     &            XX2, YY1, YY2, TITLE, DEV1, EIGHT, NINE, VERSN,
     &            VAXVMS)
               DRAWN = .TRUE.
               XMIN1 = XMIN
               XMAX1 = XMAX
               YMIN1 = YMIN
               YMAX1 = YMAX
            END IF
         END IF
C
C  PLOT ELEMENTS BY MATERIAL NUMBER
C
      ELSE IF ((CIN (ICOM) (1:1) .EQ. 'M') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'm')) THEN
         ICOM = ICOM + 1
         SETFLG = .FALSE.
         OLD = .TRUE.
         CALL FLAGK (NPELEM, NNXK, NXK, MAPDXG, 1, KKK, SETFLG, OLD)
         GOPLOT = .FALSE.
         SETFLG = .TRUE.
         CALL MESAGE ('PLOT MATERIALS FROM <I1> TO <I2>')
         CALL MESAGE ('HIT RETURN TO END INPUT')
  150    CONTINUE
         IF (ICOM .GT. JCOM) THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN, IIN,
     &         RIN)
            ICOM = 1
         END IF
         CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2, IFOUND)
         IF (IFOUND .GT. 0) THEN
            IF (I1 .GT. 0) THEN
               IF (I2 .LT. I1)I2 = I1
               DO 170 IM = I1, I2
                  DO 160 I = 1, NUMMAT
                     IF (IM .EQ. MATMAP (1, I)) THEN
                        CALL FLAGK (NPELEM, NNXK, NXK, MAPDXG,
     &                     MATMAP (2, I), MATMAP (3, I), SETFLG, OLD)
                        GOPLOT = .TRUE.
                     END IF
  160             CONTINUE
  170          CONTINUE
               GO TO 150
            END IF
         END IF
         IF (GOPLOT) THEN
            IF (ALPHA) THEN
               CALL MESAGE ('NO PLOTTING POSSIBLE ON ALPHA-NUMERIC '//
     &            'TERMINAL')
            ELSE
               CALL MNMXK (NPELEM, NPNODE, NNXK, NXK, XN, YN, CENTK,
     &            KKK, XMIN, XMAX, YMIN, YMAX)
               CALL PMESH (NPNODE, NPELEM, MXNFLG, MXSFLG, NPNBC, NPSBC,
     &            MAXKXN, NNN, KKK, NNXK, NBCNOD, NBCSID, NNLIST,
     &            NVLIST, NODES, NSIDEN, NNFLG, NNPTR, NSFLG, NVPTR,
     &            NVLEN, XN, YN, NXK, MAT, MAPDXG, MAPGXD, WTNODE,
     &            WTSIDE, AXIS, AREACG, LABE, LABO, LABN, LABNB, LABSB,
     &            LABM, LABW, CENTK, ILOOK, XMIN, XMAX, YMIN, YMAX, XX1,
     &            XX2, YY1, YY2, TITLE, DEV1, EIGHT, NINE, VERSN,
     &            VAXVMS)
               DRAWN = .TRUE.
            END IF
         END IF
C
C  SPAWN A PROCESS
C
      ELSE IF ((CIN (ICOM) (1:2) .EQ. 'SP') .OR.
     &   (CIN (ICOM) (1:2) .EQ. 'sp')) THEN
         ICOM = ICOM + 1
         CALL SPAWN (VAXVMS)
C
C  SHOW STATUS OF ALL TOGGLES
C
      ELSE IF ((CIN (ICOM) (1:1) .EQ. 'S') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 's')) THEN
         ICOM = ICOM + 1
         CALL MESAGE (' ')
         CALL MESAGE ('THE CURRENT STATUS OF ALL PLOTTING TOGGLES IS:')
         IF (AXIS) THEN
            CALL MESAGE ('   AXIS PLOTTING                - ON')
         ELSE
            CALL MESAGE ('   AXIS PLOTTING                - OFF')
         END IF
         IF (AREACG) THEN
            CALL MESAGE ('   AREA AND C.G. REPORT         - ON')
         ELSE
            CALL MESAGE ('   AREA AND C.G. REPORT         - OFF')
         END IF
         IF (LABN) THEN
            CALL MESAGE ('   LABELING OF NODES            -  ON')
         ELSE
            CALL MESAGE ('   LABELING OF NODES            -  OFF')
         END IF
         IF (LABNB) THEN
            CALL MESAGE ('   LABELING OF NODAL BOUND.     -  ON')
         ELSE
            CALL MESAGE ('   LABELING OF NODAL BOUND.     -  OFF')
         END IF
         IF (LABW) THEN
            CALL MESAGE ('   LABELING OF NODE WEIGHTING   -  ON')
         ELSE
            CALL MESAGE ('   LABELING OF NODE WEIGHTING   -  OFF')
         END IF
         IF (LABE) THEN
            CALL MESAGE ('   LABELING OF ELEMENTS         -  ON')
         ELSE
            CALL MESAGE ('   LABELING OF ELEMENTS         -  OFF')
         END IF
         IF (LABSB) THEN
            CALL MESAGE ('   LABELING OF ELEM SIDE BOUND. - ON')
         ELSE
            CALL MESAGE ('   LABELING OF ELEM SIDE BOUND. - OFF')
         END IF
         IF (LABM) THEN
            CALL MESAGE ('   LABELING OF BLOCK ID  (MAT)  - ON')
         ELSE
            CALL MESAGE ('   LABELING OF BLOCK ID  (MAT)  - OFF')
         END IF
         CALL MESAGE (' ')
         CALL MESAGE ('*----------------- NOTE -----------------*')
         CALL MESAGE ('    PLOTTING ORDER AT NODES IS:           ')
         CALL MESAGE ('        NODE NO./NODAL BOUND. FLAG/WEIGHT ')
         CALL MESAGE ('    PLOTTING ORDER AT ELEMENT CENTER IS:  ')
         CALL MESAGE ('        ELEMENT NO./BLOCK ID  (MAT) NO.    ')
         CALL MESAGE ('*----------------- NOTE -----------------*')
C
C  GET A QMS PLOT FILE OF THE CURRENT SCREEN
C
      ELSE IF (( (CIN (ICOM) (1:1) .EQ. 'H') .OR.
     &   (CIN (ICOM) (1:2) .EQ. 'h')) .AND.
     &   (CIN (ICOM) (2:2).NE.'E') .AND.
     &   (CIN (ICOM) (2:2).NE.'e')) THEN
         ICOM = ICOM + 1
         CALL VDIQES (10002, KAVAL2)
         IF (KAVAL2 .EQ. 1) THEN
            IF (.NOT.ALPHA)CALL VDESCP (10002, 0, 0)
            CALL PMESH (NPNODE, NPELEM, MXNFLG, MXSFLG, NPNBC, NPSBC,
     &         MAXKXN, NNN, KKK, NNXK, NBCNOD, NBCSID, NNLIST, NVLIST,
     &         NODES, NSIDEN, NNFLG, NNPTR, NSFLG, NVPTR, NVLEN, XN, YN,
     &         NXK, MAT, MAPDXG, MAPGXD, WTNODE, WTSIDE, AXIS, AREACG,
     &         LABE, LABO, LABN, LABNB, LABSB, LABM, LABW, CENTK, ILOOK,
     &         XMIN, XMAX, YMIN, YMAX, XX1, XX2, YY1, YY2, TITLE, 'XXX',
     &         EIGHT, NINE, VERSN, VAXVMS)
            IF (.NOT.ALPHA) THEN
               CALL PLTFLU
               CALL VDESCP (10001, 0, 0)
            END IF
            CALL MESAGE ('HARDCOPY PLOT GENERATED')
            HARDPL = .TRUE.
         ELSE
            CALL MESAGE ('HARDCOPY DEVICE NOT AVAILABLE')
         END IF
C
C  ENTER ZOOM LOCATION
C
      ELSE IF ((CIN (ICOM) (1:1) .EQ. 'Z') .OR.
     &   (CIN (ICOM) (1:1) .EQ. 'z')) THEN
         ICOM = ICOM + 1
         CALL ZOOMLT (MCOM, ICOM, JCOM, CIN, RIN, IIN, KIN, IDUMP,
     &      DRAWN, ALPHA, DEV1, X1, X2, Y1, Y2, XX1, XX2, YY1, YY2,
     &      XMIN1, XMAX1, YMIN1, YMAX1, XMIN, XMAX, YMIN, YMAX)
         DRAWN = .FALSE.
C
C  EXIT OPTION - EXITS FASTQ
C
      ELSE IF ((CIN (ICOM) (1:2) .EQ. 'EX') .OR.
     &   (CIN (ICOM) (1:2) .EQ. 'ex')) THEN
         ICOM = ICOM + 1
         CALL STRLNG (CIN(ICOM), LEN)
         IF (((LEN .GT. 1) .AND. (CIN(ICOM)(2:2) .NE. 'X')) .OR.
     &      ((LEN .GT. 1) .AND. (CIN(ICOM)(2:2) .NE. 'x'))) THEN
            CALL HELP_FQ (11)
         ELSE
            CALL FEXIT (WROTE, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN,
     &         TIME1, BATCH, VERSN)
         ENDIF
         GO TO 110
C
C  RETURN TO MESH ROUTINE
C
      ELSE IF (CIN (ICOM) (1:1) .EQ. ' ') THEN
         ICOM = ICOM + 1
         DO 180 I = 1, KKK
            NXK (1, I) = IABS (NXK (1, I))
  180    CONTINUE
         RETURN
C
C  GET HELP MESSAGE
C
      ELSE
         ICOM = ICOM + 1
         CALL HELP_FQ (11)
      END IF
      GO TO 110
C
      END
