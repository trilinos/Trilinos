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

C $Id: linlen.f,v 1.3 2000/11/13 15:39:04 gdsjaar Exp $
C $Log: linlen.f,v $
C Revision 1.3  2000/11/13 15:39:04  gdsjaar
C Cleaned up unused variables and labels.
C
C Removed some real to int conversion warnings.
C
C Revision 1.2  1991/03/21 15:44:53  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:11:13  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:11:11  gdsjaar
c Initial revision
c
C
CC* FILE: [.PAVING]LINLEN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE LINLEN (MP, COOR, LINKP, KNUM, LNUM, KT, I3, J1, J2,
     &   J3, DIST, ERR)
C***********************************************************************
C
C  SUBROUTINE LINLEN = CALCULATES THE LENGTH OF A GIVEN LINE
C
C***********************************************************************
C
C  VARIABLES USED:
C     NID   = AN ARRAY OF UNIQUE NODE IDENTIFIERS.
C     REAL  = .TRUE. FOR AN ACTUAL GENERATION
C           = .FALSE. FOR A TRIAL GENERATION
C     ERR   = .TRUE. IF AN ERROR WAS ENCOUNTERED
C     J1    = POINTER FOR THE FIRST POINT
C     J2    = POINTER FOR THE SECOND POINT
C     J3    = POINTER FOR THE THIRD POINT
C     MAXNP = MAXIMUM NUMBER OF NODES ON THE PERIMETER
C             NOTE: MAXNP MUST BE ADJUSTED FOR THE CURRENT
C                   LOCATION IN X, Y, & NID
C     KT    = THE LINE TYPE:
C           = 1 FOR STRAIGHT LINES
C           = 2 FOR CORNER LINES
C           = 3 FOR ARC WITH CENTER GIVEN
C           = 4 FOR ARC WITH THIRD POINT ON THE ARC
C           = 5 FOR PARABOLA
C           = 6 FOR ARC WITH RADIUS GIVEN
C
C***********************************************************************
C
      DIMENSION COOR (2, MP), LINKP (2, MP)
C
      LOGICAL ERR
C
      PI = ATAN2(0.0, -1.0)
C
      DIST = 0.
      ERR = .TRUE.
C
C  STRAIGHT LINE GENERATION
C
      IF (KT.EQ.1) THEN
         YDIFF = COOR (2, J2) -COOR (2, J1)
         XDIFF = COOR (1, J2) -COOR (1, J1)
         DIST = SQRT (YDIFF ** 2 + XDIFF ** 2)
         IF (DIST.EQ.0.) THEN
            WRITE (*, 10000) KNUM
            RETURN
         ENDIF
C
C  CORNER GENERATION
C
      ELSEIF (KT.EQ.2) THEN
         XDA = COOR (1, J3) -COOR (1, J1)
         YDA = COOR (2, J3) -COOR (2, J1)
         XDB = COOR (1, J2) -COOR (1, J3)
         YDB = COOR (2, J2) -COOR (2, J3)
         DA = SQRT (XDA ** 2 + YDA ** 2)
         DB = SQRT (XDB ** 2 + YDB ** 2)
         IF ((DA.EQ.0.) .OR. (DB.EQ.0.) )THEN
            WRITE (*, 10000) KNUM
            RETURN
         ENDIF
         DIST = DA+DB
C
C  CIRCULAR ARC
C
      ELSEIF ((KT.EQ.3) .OR. (KT.EQ.4) .OR. (KT.EQ.6) )THEN
         XSTART = COOR (1, J1)
         YSTART = COOR (2, J1)
         CALL ARCPAR (MP, KT, KNUM, COOR, LINKP, J1, J2, J3, I3,
     &      XCEN, YCEN, THETA1, THETA2, TANG, R1, R2, ERR, ICCW, ICW,
     &      XK, XA)
C
C  GENERATE THE CIRCLE
C
         ANG = THETA1
         DEL = TANG/30
         DO 100 I = 2, 29
            ANG = ANG+DEL
            RADIUS = XA * EXP (XK * ANG)
            XEND = XCEN+COS (ANG) * RADIUS
            YEND = YCEN+SIN (ANG) * RADIUS
            DIST = DIST+SQRT ((XEND-XSTART) ** 2 + (YEND-YSTART) ** 2)
            XSTART = XEND
            YSTART = YEND
  100    CONTINUE
         XEND = COOR (1, J2)
         YEND = COOR (2, J2)
         DIST = DIST+SQRT ((XEND-XSTART) ** 2 + (YEND-YSTART) ** 2)
C
C  ELIPSE
C
      ELSEIF (KT .EQ. 7) THEN
         XSTART = COOR (1, J1)
         YSTART = COOR (2, J1)
         CALL ELPSPR (MP, KT, KNUM, COOR, LINKP, J1, J2, J3,
     &      I3, XCEN, YCEN, THETA1, THETA2, TANG, IDUM1, IDUM2,
     &      AVALUE, BVALUE, ERR)
C
C  GENERATE THE ELIPSE
C
         ANG = THETA1
         DEL = TANG/30
         DO 110 I = 2, 29
            ANG = ANG+DEL
            RADIUS = SQRT ( (AVALUE **2 * BVALUE **2) /
     &         ( (BVALUE **2 * COS (ANG) **2) +
     &         (AVALUE **2 * SIN (ANG) **2) ) )
            XEND = XCEN+COS (ANG) * RADIUS
            YEND = YCEN+SIN (ANG) * RADIUS
            DIST = DIST+SQRT ((XEND-XSTART) ** 2 + (YEND-YSTART) ** 2)
            XSTART = XEND
            YSTART = YEND
  110    CONTINUE
         XEND = COOR (1, J2)
         YEND = COOR (2, J2)
         DIST = DIST+SQRT ((XEND-XSTART) ** 2 + (YEND-YSTART) ** 2)
C
C     PARABOLA
C
      ELSEIF (KT.EQ.5) THEN
C
C  CHECK LEGITIMACY OF DATA
C
         XMID =  (COOR (1, J1) +COOR (1, J2) ) * 0.5
         YMID =  (COOR (2, J1) +COOR (2, J2) ) * 0.5
         DOT =  (COOR (1, J2) -COOR (1, J1) ) * (COOR (1, J3) -XMID)
     &      + (COOR (2, J2) -COOR (2, J1) ) * (COOR (2, J3) -YMID)
         PERP = SQRT ((COOR (1, J2) -COOR (1, J1) ) ** 2 +
     &      (COOR (2, J2) - COOR (2, J1) ) ** 2) *
     &      SQRT ((COOR (1, J3) -XMID) ** 2 + (COOR (2, J3)
     &      -YMID) ** 2)
         IF (DOT.GE.0.05 * PERP) THEN
            WRITE (*, 10030) KNUM
            RETURN
         ENDIF
C
C  GETARC LENGTH
C
         HALFW = SQRT ((COOR (1, J2) -COOR (1, J1) ) ** 2 +
     &      (COOR (2, J2) - COOR (2, J1) ) ** 2 )  * 0.5
         IF (HALFW.EQ.0.) THEN
            WRITE (*, 10000) KNUM
            RETURN
         ENDIF
         HEIGHT = SQRT ((XMID-COOR (1, J3) ) ** 2
     &      + (YMID-COOR (2, J3) ) ** 2)
         COEF = HEIGHT/HALFW ** 2
         TCOEF = 2.0 * COEF
C
C  PARC IS A STATEMENT FUNCTION
C
         PLEFT = PARC (-TCOEF * HALFW, TCOEF)
         ARCTOT = 2.0 * PARC (TCOEF * HALFW, TCOEF)
         ARCDEL = ARCTOT/30
         ARCNXT = ARCDEL
         ARCNOW = 0.0
         THETA = ATAN2 (COOR (2, J2) -COOR (2, J1) , COOR (1, J2)
     &      - COOR (1, J1) )
C
C  CORRECT FOR ORIENTATION
C
         CROSS =  (COOR (1, J3) -XMID) *  (COOR (2, J2) -COOR (2, J1) )-
     &      (COOR (2, J3) -YMID) *  (COOR (1, J2) -COOR (1, J1) )
         IF (CROSS.LT.0.0) THETA = THETA+PI
         SINT = SIN (THETA)
         COST = COS (THETA)
C
C  FIND POINTS APPROXIMATELY BY INTEGRATION
C
         XL = -HALFW
         FL = SQRT (1.0+ (TCOEF * XL) ** 2)
         KOUNT = 1
         DELX = 2.0 * HALFW/200.0
         XSTART = COOR (1, J1)
         YSTART = COOR (2, J1)
         DO 120 I = 1, 100
            FM = SQRT (1.0+ (TCOEF * (XL+DELX) ) ** 2)
            XR = - HALFW + FLOAT (I) * 2.0 * DELX
            FR = SQRT (1.0+ (TCOEF * XR) ** 2)
            ARCOLD = ARCNOW
            ARCNOW = ARCNOW+DELX * (FL+4.0 * FM+FR) / 3.0
            IF (ARCNOW.GE.ARCNXT) THEN
C
C  COMPUTE POSITION IN LOCAL COORDINATE SYSTEM
C
               FRAC =  (ARCNXT-ARCOLD) / (ARCNOW-ARCOLD)
               XK = XL+FRAC * 2.0 * DELX
               YK = COEF * XK ** 2
C
C  CORRECT FOR ORIENTATION PROBLEM
C
               IF (CROSS.LT.0.0) XK = -XK
C
C  ROTATE IN LINE WITH GLOBAL COORDINATE SYSTEM
C
               ROTX = XK * COST - YK * SINT
               ROTY = YK * COST + XK * SINT
C
C  RESTORE XK
C
               IF (CROSS.LT.0.0) XK = -XK
C
C  TRANSLATE
C
               XEND = ROTX+COOR (1, J3)
               YEND = ROTY+COOR (2, J3)
               DIST = DIST+SQRT ((XEND-XSTART) ** 2 + (YEND-YSTART) **2)
               KOUNT = KOUNT+1
               XSTART = XEND
               YSTART = YEND
C
C  PREPARE FOR NEXT POINT
C
               IF (KOUNT.GE.29) GOTO 130
               ARCNXT = ARCNXT+ARCDEL
C
C  RESTART INTEGRATION
C
               XR = XK
               FR = SQRT (1.0+ (TCOEF * XR) ** 2)
C
C  CORRECT FOR INTEGRATION ERROR
C
               ARCNOW = PARC (TCOEF * XR, TCOEF) -PLEFT
            ENDIF
            XL = XR
            FL = FR
  120    CONTINUE
  130    CONTINUE
         XEND = COOR (1, J2)
         YEND = COOR (2, J2)
         DIST = DIST+SQRT ((XEND-XSTART) ** 2+ (YEND-YSTART) ** 2)
      ENDIF
C
C     NORMAL EXIT
C
      ERR = .FALSE.
      RETURN
C
10000 FORMAT (' ZERO LINE LENGTH ENCOUNTERED FOR LINE', I5)
10030 FORMAT (' POINTS GIVEN FOR LINE', I5, ' DO NOT DEFINE A PARABOLA')
C
      END
