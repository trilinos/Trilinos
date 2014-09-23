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

C $Id: getsiz.f,v 1.1 1990/11/30 11:08:45 gdsjaar Exp $
C $Log: getsiz.f,v $
C Revision 1.1  1990/11/30 11:08:45  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]GETSIZ.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETSIZ (XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR,
     &   MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN,
     &   REYMAX, IDIVIS, SIZMIN, EMAX, EMIN, X, Y, SIZE)
C***********************************************************************
C
C  SUBROUTINE GETSIZ = GETS THE SIZE OF AN ELEMENT EDGE BASED ON THE
C                      OLD MESH SIZE AT THE GIVEN X,Y LOCATION AND THE
C                      RELATIVE MEASURE OF THE ERROR ESTIMATOR AT THAT
C                      LOCATION
C
C***********************************************************************
C
      DIMENSION XNOLD(NPNOLD), YNOLD(NPNOLD)
      DIMENSION NXKOLD(NNXK, NPEOLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD), BMESUR(NPNOLD)
C
      LOGICAL INSIDE, BAD
C
C  ASSUME A LINEAR REDUCTION FACTOR FROM R0 TO R1 WHERE R0 IS THE
C  DESIRED REDUCTION A 0. NORMALIZED ERROR MEASURE AND R1 IS A DESIRED
C  REDUCTION AT 1.0 NORMALIZED ERROR MEASURE
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/12/90
CC* MODIFICATION: CHANGED R0 TO BE 1.5 (WAS 1.2)
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/17/90
CC* MODIFICATION: CHANGED R0 TO BE 1.2 (WAS 1.5)
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/17/90
CC* MODIFICATION: CHANGED R1 TO BE 1/6 (WAS 1/7)
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: CHANGED CALl TO PASS IN R0 AND R1 AND RENAMED THEM
CC**              TO BE EMAX AND EMIN RESPECTIVELY - THEY ARE NOW SET
CC**              BY THE USER IN THE RDGENS ROUTINE
C
C  FIND THE ELEMENT THAT THIS POINT FALLS INTO
C
      DELX = (REXMAX - REXMIN) / FLOAT(IDIVIS)
      DELY = (REYMAX - REYMIN) / FLOAT(IDIVIS)
      IX = INT((X - REXMIN) / DELX) + 1
      IF (X .GE. REXMAX) IX = IDIVIS
      IY = INT((REYMAX - Y) / DELY) + 1
      IF (Y .LE. REYMIN) IY = IDIVIS
      INDEX = (IY * 10) + IX
      IBEGIN = LINKEG (1, INDEX)
      IEND = IBEGIN + LINKEG(2,INDEX) - 1
      DO 110 I = IBEGIN, IEND
         KELEM = LISTEG(I)
         XEMIN = XNOLD (NXKOLD (1, KELEM))
         XEMAX = XNOLD (NXKOLD (1, KELEM))
         YEMIN = YNOLD (NXKOLD (1, KELEM))
         YEMAX = YNOLD (NXKOLD (1, KELEM))
         DO 100 IC = 2, 4
            XEMIN = AMIN1 (XEMIN, XNOLD (NXKOLD (IC, KELEM)))
            XEMAX = AMAX1 (XEMAX, XNOLD (NXKOLD (IC, KELEM)))
            YEMIN = AMIN1 (YEMIN, YNOLD (NXKOLD (IC, KELEM)))
            YEMAX = AMAX1 (YEMAX, YNOLD (NXKOLD (IC, KELEM)))
  100    CONTINUE
         IF ((X .LT. XEMIN) .OR. (X .GT. XEMAX) .OR.
     &      (Y .LT. YEMIN) .OR. (Y .GT. YEMAX)) THEN
            INSIDE = .FALSE.
         ELSE
            CALL INVMAP (X, Y,
     &         XNOLD(NXKOLD (1, KELEM)), YNOLD(NXKOLD (1, KELEM)),
     &         XNOLD(NXKOLD (2, KELEM)), YNOLD(NXKOLD (2, KELEM)),
     &         XNOLD(NXKOLD (3, KELEM)), YNOLD(NXKOLD (3, KELEM)),
     &         XNOLD(NXKOLD (4, KELEM)), YNOLD(NXKOLD (4, KELEM)),
     &         XI, ETA, INSIDE)
         ENDIF
         IF (INSIDE) THEN
            KIN = KELEM
            GOTO 170
         ENDIF
  110 CONTINUE
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 8/1/90
CC* MODIFICATION: ADDED THE CHECK TO SEE IF A POINT BEING SIZED IS ACTUALLY
CC*               OUTSIDE THE INITIAL MESH BUT INSIDE THE GEOMETRY.
C
C  THERE IS A POSSIBILITY THAT THE POINT IS ON AN ARC WHICH IS NOT
C  INCLUDED IN THE ORIGINAL MESH - THIS MUST BE CHECKED.
C
      DTEST = 1.E15
      DO 130 I = IBEGIN, IEND
         KELEM = LISTEG(I)
         DO 120 IC = 1, 4
            JC = IC + 1
            IF (JC .EQ. 5) JC = 1
            X1 = XNOLD (NXKOLD (IC, KELEM))
            X2 = XNOLD (NXKOLD (JC, KELEM))
            Y1 = YNOLD (NXKOLD (IC, KELEM))
            Y2 = YNOLD (NXKOLD (JC, KELEM))
C
C  GET THE PARAMETERS FOR THE LINE
C
            CALL DLPARA (X1, Y1, X2, Y2, XM1, B1, BAD)
C
C  GET DISTANCE FOR VERTICAL LINE
C
            IF (BAD) THEN
               DTRY = ABS(X1 - X)
               XTRY = X1
               YTRY = Y
C
C  GET DISTANCE FOR HORIZONTAL LINE
C
            ELSE IF (ABS(XM1) .LT. .000001) THEN
               DTRY = ABS(Y1 - Y)
               XTRY = X
               YTRY = Y1
C
C  GET PERPENDICULAR DISTANCE TO ARBITRARY LINE
C
            ELSE
               XM2 = -1./XM1
               B2 = Y - (XM2*X)
               XTRY = (B2 - B1)/(XM1 - XM2)
               YTRY = (XM1*XTRY) + B1
               DTRY = SQRT((X - XTRY)**2 + (Y - YTRY)**2)
            END IF
C
C  CHECK THE INTERSECTION TO MAKE SURE THAT IT CUTS THE LINE SEGMENT
C  WE HAVE
C
            IF ((XTRY .GE. AMIN1(X1, X2)) .AND.
     &         (XTRY .LE. AMAX1(X1, X2)) .AND.
     &         (YTRY .GE. AMIN1(Y1, Y2)) .AND.
     &         (YTRY .LE. AMAX1(Y1, Y2)) ) THEN
C
C  NOW GET THE SHORTEST INTERSECTION AND GET NEEDED SIZE VALUE BASED ON
C  THE XTRY AND YTRY LOCATION
C
               IF (DTRY .LT. DTEST) THEN
                  DTEST = DTRY
                  INSIDE = .TRUE.
                  DT = SQRT( (X1 - X2)**2 + (Y1 - Y2)**2)
                  D1 = SQRT( (XTRY - X1)**2 + (YTRY - Y1)**2)
                  RATIO = D1 / DT
                  IF (IC .EQ. 1) THEN
                     XI = RATIO
                     ETA = 0.
                  ELSEIF (IC .EQ. 2) THEN
                     XI = 0.
                     ETA = RATIO
                  ELSEIF (IC .EQ. 3) THEN
                     XI = 1.0 - RATIO
                     ETA = 0.
                  ELSE
                     XI = 0.
                     ETA = 1.0 - RATIO
                  ENDIF
                  KIN = KELEM
                  ICIN = IC
                  JCIN = JC
               ENDIF
            ENDIF
  120    CONTINUE
  130 CONTINUE
C
C  NOW CHECK THE ELEMENT THAT HAS BEEN FOUND AND MAKE SURE THAT IT IS
C  A ELEMENT ALONG THE SIDE OF THE MESH AND THAT THE EDGE CLOSEST IS
C  NOT SHARED BY ANY OTHER ELEMENT.
C
      IF (INSIDE) THEN
         DO 150 I = 1, NPEOLD
            IF (I .NE. KIN) THEN
               DO 140 IC = 1, 4
                  JC = IC + 1
                  IF (JC .EQ. 5) JC = 1
                  IF ((IC .EQ. JCIN) .AND. (JC .EQ. ICIN)) THEN
                     CALL MESAGE ('** ERROR WITH ELEMENT SIDE FOUND'//
     &                  ' BEING INTERIOR TO MESH IN GETSIZ **')
                     INSIDE = .FALSE.
                     GOTO 160
                  ENDIF
  140          CONTINUE
            ENDIF
  150    CONTINUE
  160    CONTINUE
      ENDIF
C
C  THE ELMENT HAS BEEN FOUND - NOW INTERPOLATE THE STRESS VALUE FOR
C  THIS LEVEL
C
  170 CONTINUE
      IF (INSIDE) THEN
         N1 = NXKOLD (1, KIN)
         N2 = NXKOLD (2, KIN)
         N3 = NXKOLD (3, KIN)
         N4 = NXKOLD (4, KIN)
         E1 = BMESUR(N1)
         E2 = BMESUR(N2)
         E3 = BMESUR(N3)
         E4 = BMESUR(N4)
         ERROR = E1 + ((E2 - E1) * XI) + ((E4 - E1) * ETA) +
     &      ((E1 - E2 + E3 - E4) * XI * ETA)
         D1 = SQRT ( ((XNOLD(N2) - XNOLD(N1)) ** 2) +
     &      ((YNOLD(N2) - YNOLD(N1)) ** 2) )
         D2 = SQRT ( ((XNOLD(N3) - XNOLD(N2)) ** 2) +
     &      ((YNOLD(N3) - YNOLD(N2)) ** 2) )
         D3 = SQRT ( ((XNOLD(N4) - XNOLD(N3)) ** 2) +
     &      ((YNOLD(N4) - YNOLD(N3)) ** 2) )
         D4 = SQRT ( ((XNOLD(N1) - XNOLD(N4)) ** 2) +
     &      ((YNOLD(N1) - YNOLD(N4)) ** 2) )
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 8/6/90
CC* MODIFICATION: CHANGED FROM A SIZE INTERPOLATION TO THE
CC**              MINIMUM SIDE LENGTH AS THE EXISTING ELEMENT SIZE
C
         REDUC = EMAX - (ERROR * EMAX) + (ERROR * EMIN)
         SIZE = AMAX1 ((AMIN1 (D1, D2, D3, D4) * REDUC), SIZMIN)
      ELSE
C
C  ERROR HAS OCCURRED IN FINDING THE ELEMENT
C
         CALL MESAGE ('** ERROR - ENCLOSING ELEMENT NOT FOUND IN '//
     &      'GETSIZ **')
      ENDIF
C
      RETURN
C
      END
