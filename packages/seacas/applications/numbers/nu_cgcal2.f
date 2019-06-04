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

C $Id: cgcal2.f,v 1.3 2000/07/06 18:07:42 gdsjaar Exp $
C $Log: cgcal2.f,v $
C Revision 1.3  2000/07/06 18:07:42  gdsjaar
C Fix assumption that variables are saved between subroutine calls
C
C Revision 1.2  1999/02/16 21:37:58  gdsjaar
C Converted to read exodusII database format.  Somewhat tested, not
C ready for production yet.
C
C Revision 1.1.1.1  1991/02/21 15:42:24  gdsjaar
C NUMBERS: Greg Sjaardema, initial Unix release
C
c Revision 1.1  1991/02/21  15:42:23  gdsjaar
c Initial revision
c
      SUBROUTINE CGCAL2(CRD,IX,MAT,MASS,VOL,DENS,VOLM,CG,ZITOT,XXX,
     *    XG,XI,XINI,AJ,NNODES,NDIM,NQUAD,VOLMN,IELM,NELBLK,
     *    AXI,NUMNP)
C
      DIMENSION CRD(NUMNP,*), IX(NNODES,*), MAT(6,*), MASS(*),
     *    DENS(*), VOLM(*), CG(*), ZITOT(*),VOLMN(4,*),IELM(4,*),
     *    XXX(NDIM+1,NQUAD,*),XG(NDIM,*), XI(NDIM,*), XINI(*),
     *    AJ(2,*)
      DIMENSION ZI(6), ZMOM(3)
      DIMENSION CCC(2,4)
C
      LOGICAL AXI
      REAL MASS, MASSE
      PI = ATAN2(0.0, -1.0)
C
C ... VOLMN(1,*) = MINIMUM VOLUME (AREAS FOR 2-D)
C     VOLMN(2,*) = MAXIMUM VOLUME
C     VOLMN(3,*) = TOTAL VOLUME
C
      DO 10 I=1, NELBLK
          VOLMN(1,I) = 1.0E30
          VOLMN(2,I) = 0.0
          VOLMN(3,I) = 0.0
          IELM (3,I) = 0
          MASS(I)    = 0.0
          VOLM(I)    = 0.0
   10 CONTINUE
      DO 20 I=1,3
          ZITOT(I)   = 0.0
          ZITOT(I+3) = 0.0
          ZMOM(I)    = 0.0
   20 CONTINUE
      ZMAS = 0.0
      VOL  = 0.0
C
C ... GET QUADRATURE POINT LOCATIONS, EVALUATE SHAPE FUNCTIONS
C
      CALL QUAD(XXX, XI, XG, NDIM, NNODES, NQUAD, WT)
C
      DO 110 IBLK = 1, NELBLK
          IF (MAT(5,IBLK) .NE. 1) GOTO 110
          IELBEG = MAT(3,IBLK)
          IELEND = MAT(4,IBLK)
          MIEL   = IBLK
          DO 100 IEL = IELBEG, IELEND
C
C ... CALCULATE AREA, VOLUME, AND MOMENTS OF INERTIA OF ELEMENT
C
              DO 30 I=1,3
                  ZI(I)   = 0.0
                  ZI(I+3) = 0.0
                  CG(I)   = 0.0
   30         CONTINUE
              VOLUME = 0.0
C
              DO 40 I=1,4
                  CCC(1,I) = CRD(IX(I,IEL),1)
                  CCC(2,I) = CRD(IX(I,IEL),2)
   40         CONTINUE
              DO 80 NG=1,NQUAD
                  DET = 0.0
                  DO 60 J=1,2
                      XINI(J) = 0.0
                      DO 50 K=1,2
                          AJ(K,J) = 0.0
   50                 CONTINUE
   60             CONTINUE
C
                  DO 70 I=1,4
                      XINI(1) = XINI(1)+XXX(1,I,NG) * CCC(1,I)
                      AJ(1,1) = AJ(1,1)+XXX(2,I,NG) * CCC(1,I)
                      AJ(2,1) = AJ(2,1)+XXX(3,I,NG) * CCC(1,I)

                      XINI(2) = XINI(2)+XXX(1,I,NG) * CCC(2,I)
                      AJ(1,2) = AJ(1,2)+XXX(2,I,NG) * CCC(2,I)
                      AJ(2,2) = AJ(2,2)+XXX(3,I,NG) * CCC(2,I)
   70             CONTINUE
C
                  DET = ( AJ(1,1)*AJ(2,2) - AJ(2,1)*AJ(1,2) ) * WT
                  DETW = DET * DENS(MIEL)
C
                  IF (AXI) THEN
C
C ... CG(1) IS THE VOLUME FOR AXI 2-D, VOLUME IS ACTUALLY C/S AREA
C
                      CG(1) = CG(1) + DET * XINI(1)
                      CG(2) = CG(2) + DETW * XINI(1) * XINI(2)
                      ZI(2) = ZI(2) + DETW * XINI(1)**3
                      ZI(1) = ZI(1) + DETW * XINI(1)*XINI(2)**2
                      VOLUME  = VOLUME  + DET
                  ELSE
                      CG(1) = CG(1) + DETW * XINI(1)
                      CG(2) = CG(2) + DETW * XINI(2)
                      ZI(1) = ZI(1) + DETW * XINI(2)**2
                      ZI(2) = ZI(2) + DETW * XINI(1)**2
                      ZI(3) = ZI(3) + DETW * XINI(1)*XINI(2)
                      VOLUME = VOLUME + DET
                  END IF
C
   80         CONTINUE
C
C ... DETERMINE MIN/MAX ELEMENT VOLUMES FOR EACH MATERIAL AND
C        COUNT NUMBER OF ELEMENTS FOR EACH MATERIAL
C
              IELM(3,MIEL)      = IELM(3,MIEL) + 1
              VOLMN(3,MIEL)     = VOLMN(3,MIEL) + VOLUME
              IF (VOLUME .LT. VOLMN(1,MIEL)) THEN
                  VOLMN(1,MIEL) = VOLUME
                  IELM(1,MIEL)  = IEL
              ELSE IF (VOLUME .GT. VOLMN(2,MIEL)) THEN
                  VOLMN(2,MIEL) = VOLUME
                  IELM(2,MIEL)  = IEL
              ENDIF
C
              IF (AXI) THEN
                  VOLUME = 2. * PI * CG(1)
                  ZI(2)  = ZI(2) * 2. * PI
                  ZI(1)  = ZI(1) * 2. * PI + ZI(2) / 2.0
                  ZI(3)  = ZI(1)
              END IF
              DO 90 I=1,3
                  ZITOT(I) = ZITOT(I) + ZI(I)
   90         CONTINUE
              MASSE = VOLUME * DENS(MIEL)
              MASS(MIEL)= MASS(MIEL) + MASSE
              VOLM(MIEL)= VOLM(MIEL) + VOLUME
              VOL  = VOL  + VOLUME
              ZMAS = ZMAS + MASSE
              IF (AXI) THEN
                  CG(1)   = 0.0
                  ZMOM(2) = ZMOM(2) + CG(2) * 2. * PI
              ELSE
                  ZMOM(1) = ZMOM(1) + CG(1)
                  ZMOM(2) = ZMOM(2) + CG(2)
              END IF
  100     CONTINUE
  110 CONTINUE
      FIX = SIGN(0.5, ZMAS) + SIGN(0.5, -ZMAS)
      DO 120 I=1,3
          CG(I) = ZMOM(I) / (ZMAS + FIX)
  120 CONTINUE
      IF (AXI) THEN
          ZITOT(1) = ZITOT(1) - ZMAS * CG(2)**2
          ZITOT(3) = ZITOT(3) - ZMAS * CG(2)**2
      ELSE
          ZITOT(1) = ZITOT(1) - ZMAS * CG(2)**2
          ZITOT(2) = ZITOT(2) - ZMAS * CG(1)**2
          ZITOT(4) = ZITOT(3) - ZMAS * CG(1) * CG(2)
          ZITOT(3) = ZITOT(1) + ZITOT(2)
      END IF
C
      RETURN
      END
