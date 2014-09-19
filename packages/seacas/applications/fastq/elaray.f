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

C $Id: elaray.f,v 1.3 2000/11/13 15:39:04 gdsjaar Exp $
C $Log: elaray.f,v $
C Revision 1.3  2000/11/13 15:39:04  gdsjaar
C Cleaned up unused variables and labels.
C
C Removed some real to int conversion warnings.
C
C Revision 1.2  1991/04/10 19:56:49  gdsjaar
C Fixed some logical variables
C
c Revision 1.1.1.1  1990/11/30  11:06:29  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:06:28  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]ELARAY.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE ELARAY (XNOLD, YNOLD, NXKOLD, MMPOLD, LINKEG, LISTEG,
     &   MLINK, NPROLD, NPNOLD, NPEOLD, NNXK, XMIN, XMAX, YMIN, YMAX,
     &   IDIVIS)
C***********************************************************************
C
C  SUBROUTINE ELARAY = PUTS ELEMENTS INTO AN ARRAY BASED ON THEIR
C                      PHYSICAL LOCATION
C
C***********************************************************************
C
      DIMENSION XNOLD(NPNOLD), YNOLD(NPNOLD)
      DIMENSION NXKOLD(NNXK, NPEOLD), MMPOLD(3, NPROLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD)
C
      LOGICAL LCROSS, INSIDE
C
C  FIND THE EXTREMES FOR THE MESH DATA
C
      XMIN = XNOLD(1)
      XMAX = XNOLD(1)
      YMIN = YNOLD(1)
      YMAX = YNOLD(1)
      DO 100 I = 2, NPNOLD
         XMIN = AMIN1 (XMIN, XNOLD(I))
         XMAX = AMAX1 (XMAX, XNOLD(I))
         YMIN = AMIN1 (YMIN, YNOLD(I))
         YMAX = AMAX1 (YMAX, YNOLD(I))
  100 CONTINUE
C
C  SET UP THE SIZE OF THE ARRAY BASED ON THE MLINK DIMENSION
C        IF MLINK = 55 THEN THERE ARE 5 COLUMNS AND 5 ROWS
C                 = 66 THEN THERE ARE 6 COLUMNS AND 6 ROWS, ETC.
C
      IF (MLINK .EQ. 22) THEN
         IDIVIS = 2
      ELSE IF (MLINK .EQ. 33) THEN
         IDIVIS = 3
      ELSE IF (MLINK .EQ. 44) THEN
         IDIVIS = 4
      ELSE IF (MLINK .EQ. 55) THEN
         IDIVIS = 5
      ELSE IF (MLINK .EQ. 66) THEN
         IDIVIS = 6
      ELSE IF (MLINK .EQ. 77) THEN
         IDIVIS = 7
      ELSE IF (MLINK .EQ. 88) THEN
         IDIVIS = 8
      ELSE IF (MLINK .EQ. 99) THEN
         IDIVIS = 9
      ENDIF
C
C  NOW THE ELEMENTS MUST BE SORTED INTO ANY ARRAY SPACE THAT THE ELEMENT
C  CROSSES.  THE ARRAY IS LOGICALLY A SQUARE, BUT PHYSICALLY CAN BE
C  RECTANGULAR SINCE THE X AND Y EXTREMES MAY FORM ANY SIZE RECTANGLE.
C  ROWS FIRST IN THE ARRAY AND THEN COLUMNS.
C
      XDELTA = (XMAX - XMIN) / FLOAT (IDIVIS)
      YDELTA = (YMAX - YMIN) / FLOAT (IDIVIS)
      KOUNT = 0
      DO 160 J = IDIVIS, 1, -1
         IF (J .EQ. 1) THEN
            YL = YMIN
         ELSE
            YL = YMIN + (YDELTA * FLOAT(J - 1))
         ENDIF
         IF (J .EQ. IDIVIS) THEN
            YU = YMAX
         ELSE
            YU = YMIN + (YDELTA * FLOAT(J))
         ENDIF
         DO 150 I = 1, IDIVIS
            IF (I .EQ. 1) THEN
               XL = XMIN
            ELSE
               XL = XMIN + (XDELTA * FLOAT(I - 1))
            ENDIF
            IF (I .EQ. IDIVIS) THEN
               XU = XMAX
            ELSE
               XU = XMIN + (XDELTA * FLOAT(I))
            ENDIF
            INDEX = ((IDIVIS - J + 1) * 10) + I
            LINKEG (1, INDEX) = KOUNT + 1
            LINKEG (2, INDEX) = 0
C
C  ONLY CHECK ELEMENTS OF THE SAME MATERIAL ID (BLOCK ID)
C
            DO 140 KELEM = 1, NPEOLD
               DO 120 ICON = 1, 4
                  X1 = XNOLD (NXKOLD (ICON, KELEM))
                  Y1 = YNOLD (NXKOLD (ICON, KELEM))
C
C  TEST TO SEE IF THE NODE FITS IN THE GRID
C
                  IF ( ((X1 .LE. XU) .AND. (X1 .GE. XL)) .AND.
     &               ((Y1 .LE. YU) .AND. (Y1 .GE. YL))  ) THEN
                     KOUNT = KOUNT + 1
                     IF (KOUNT .GT. NPEOLD*4) THEN
                        CALL MESAGE ('** ERROR - NOT ENOUGH ROOM '//
     &                     'IN LISTEG, SUBROUTINE ELARAY **')
                        GOTO 170
                     ENDIF
                     LINKEG (2, INDEX) = LINKEG (2, INDEX) + 1
                     LISTEG (KOUNT) = KELEM
                     GOTO 130
                  ENDIF
C
C  TEST TO SEE IF THE EDGE OF THE ELEMENT CROSSES THE GRID
C
                  IF (ICON .EQ. 4) THEN
                     JCON = 1
                  ELSE
                     JCON = ICON + 1
                  ENDIF
                  X2 = XNOLD (NXKOLD (JCON, KELEM))
                  Y2 = YNOLD (NXKOLD (JCON, KELEM))
                  CALL INTSCT (X1, Y1, X2, Y2, XL, YL, XU, YL, U, W,
     &               LCROSS)
                  IF (.NOT. LCROSS) CALL INTSCT (X1, Y1, X2, Y2,
     &               XU, YL, XU, YU, U, W, LCROSS)
                  IF (.NOT. LCROSS) CALL INTSCT (X1, Y1, X2, Y2,
     &               XU, YU, XL, YU, U, W, LCROSS)
                  IF (.NOT. LCROSS) CALL INTSCT (X1, Y1, X2, Y2,
     &               XL, YU, XL, YL, U, W, LCROSS)
                  IF (LCROSS) THEN
                     KOUNT = KOUNT + 1
                     IF (KOUNT .GT. NPEOLD*4) THEN
                        CALL MESAGE ('** ERROR - NOT ENOUGH ROOM '//
     &                     'IN LISTEG, SUBROUTINE ELARAY **')
                        GOTO 170
                     ENDIF
                     LINKEG (2, INDEX) = LINKEG (2, INDEX) + 1
                     LISTEG (KOUNT) = KELEM
                     GOTO 130
                  ENDIF
C
C  OTHERWISE TEST TO SEE IF THE ELEMENT COMPLETELY ENCLOSES THE GRID
C
                  XEMIN = XNOLD (NXKOLD (1, KELEM))
                  XEMAX = XNOLD (NXKOLD (1, KELEM))
                  YEMIN = YNOLD (NXKOLD (1, KELEM))
                  YEMAX = YNOLD (NXKOLD (1, KELEM))
                  DO 110 IC = 2, 4
                     XEMIN = AMIN1 (XEMIN, XNOLD (NXKOLD (IC, KELEM)))
                     XEMAX = AMAX1 (XEMAX, XNOLD (NXKOLD (IC, KELEM)))
                     YEMIN = AMIN1 (YEMIN, YNOLD (NXKOLD (IC, KELEM)))
                     YEMAX = AMAX1 (YEMAX, YNOLD (NXKOLD (IC, KELEM)))
  110             CONTINUE
                  IF ((XL .GT. XEMIN) .OR. (XU .LT. XEMAX) .OR.
     &               (YL .GT. YEMIN) .OR. (YU .LT. YEMAX)) THEN
                     INSIDE = .FALSE.
                  ELSE
                     CALL INVMAP (X1, Y1, XL, YL, XU, YL, XU, YU, XL,
     &                  YU, XI, ETA, INSIDE)
                  ENDIF
                  IF (INSIDE) THEN
                     KOUNT = KOUNT + 1
                     IF (KOUNT .GT. NPEOLD*4) THEN
                        CALL MESAGE ('** ERROR - NOT ENOUGH ROOM '//
     &                     'IN LISTEG, SUBROUTINE ELARAY **')
                        GOTO 170
                     ENDIF
                     LINKEG (2, INDEX) = LINKEG (2, INDEX) + 1
                     LISTEG (KOUNT) = KELEM
                     GOTO 130
                  ENDIF
C
  120          CONTINUE
  130          CONTINUE
C
  140       CONTINUE
C
  150    CONTINUE

  160 CONTINUE
C
  170 CONTINUE
      RETURN
C
      END
