C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.  
C 
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 

C=====================================================================
*DECK,MKEI
      SUBROUTINE MKEI(IST,ISTP,TIMES,IDBLKA,ICONA,NDLSTA,
     &                XA,YA,ZA,VELXA,VELYA,VELZA,
     &                EMSSA,DENSA,RNMSA,
     &                TMXA,TMYA,TMZA,TKEA,TPSQA,TJ2A,
     &                SIGXXA,SIGYYA,SIGZZA,SIGXYA,SIGYZA,SIGZXA,
     &                IDBLKB,ICONB,NDLSTB,
     &                XB,YB,ZB,VELXB,VELYB,VELZB,
     &                EMSSB,DENSB,RNMSB,
     &                TMXB,TMYB,TMZB,TKEB,TPSQB,TJ2B,ICOMPL,
     &                SIGXXB,SIGYYB,SIGZZB,SIGXYB,SIGYZB,SIGZXB)
C
C     ****************************************************************
C
C Set up arrays for computing momenta and kinetic energy
C Write results to text file
C
C
C Called by MAPVAR
C
C Calls MKE
C
C     ****************************************************************
C IST    INT   time step counter for loop in MAPVAR
C ISTP   INT   Time step
C TIMES  REAL  Array of times (TIMES(ISTP)=real time
C IDBLKA INT   Donor mesh element block I.D.
C ICONA  INT   Donor mesh connectivity array
C NDLSTA INT   Donor mesh nodes in element block
C XA     REAL  Donor mesh coordinates
C YA     REAL  Donor mesh coordinates
C ZA     REAL  Donor mesh coordinates
C VEL*A  REAL  Donor mesh velocities
C EMSSA  REAL  Donor mesh element mass
C RNMSA  REAL  Donor mesh nodal mass
C TMXA   REAL  Donor mesh sum over elt block x-mom each time step
C TMYA   REAL  Donor mesh sum over elt block y-mom each time step
C TMZA   REAL  Donor mesh sum over elt block z-mom each time step
C TKEA   REAL  Donor mesh sum over elt block ke each time step
C TPSQA  REAL  Donor mesh sum over elt block pressure squared each ts
C TKEA   REAL  Donor mesh sum over elt block J2 each time step
C SIG**A REAL  Donor mesh stress components
C IDBLKB INT   Recipient mesh element block I.D.
C ICONB  INT   Recipient mesh connectivity array
C NDLSTB INT   Recipient mesh nodes in element block
C XB     REAL  Recipient mesh coordinates
C YB     REAL  Recipient mesh coordinates
C ZB     REAL  Recipient mesh coordinates
C VEL*B  REAL  Recipient mesh velocities
C EMSSB  REAL  Recipient mesh element mass
C RNMSB  REAL  Recipient mesh nodal mass
C TMXB   REAL  Rec mesh sum over elt block x-mom each time step
C TMYB   REAL  Rec mesh sum over elt block y-mom each time step
C TMZB   REAL  Rec mesh sum over elt block z-mom each time step
C TKEB   REAL  Rec mesh sum over elt block ke each time step
C TPSQA  REAL  REC mesh sum over elt block pressure squared each ts
C TKEA   REAL  REC mesh sum over elt block J2 each time step
C ICOMPL INT   Flag to indicate completion of rec mesh element block
C SIG**B REAL  Rec mesh stress components
C
C
C     ****************************************************************
C
      include 'exodusII.inc'
C
      include 'amesh.blk'
      include 'bmesh.blk'
      include 'ebbyeb.blk'
      include 'ex2tp.blk'
      include 'tapes.blk'
      include 'varnpt.blk'
      include 'varept.blk'
C
      DIMENSION TIMES(*),ICONA(NELNDA,*),NDLSTA(*)
      DIMENSION XA(*),YA(*),ZA(*),VELXA(*),VELYA(*),VELZA(*)
      DIMENSION EMSSA(*),DENSA(*),RNMSA(*)
      DIMENSION TMXA(*),TMYA(*),TMZA(*),TKEA(*),TPSQA(*),TJ2A(*)
      DIMENSION SIGXXA(*),SIGYYA(*),SIGZZA(*),
     &          SIGXYA(*),SIGYZA(*),SIGZXA(*)
      DIMENSION ICONB(NELNDB,*),NDLSTB(*)
      DIMENSION XB(*),YB(*),ZB(*),VELXB(*),VELYB(*),VELZB(*)
      DIMENSION EMSSB(*),DENSB(*),RNMSB(*)
      DIMENSION TMXB(*),TMYB(*),TMZB(*),TKEB(*),TPSQB(*),TJ2B(*)
      DIMENSION SIGXXB(*),SIGYYB(*),SIGZZB(*),
     &          SIGXYB(*),SIGYZB(*),SIGZXB(*)
      DIMENSION XX(27),YY(27),ZZ(27)
C
C     ****************************************************************
C
C get vel's and elmass for donor mesh-A
C
      CALL EXGNV(NTP2EX,ISTP,IXVEL,NODESA,VELXA,IERR)
      CALL EXGNV(NTP2EX,ISTP,IYVEL,NODESA,VELYA,IERR)
      IF (NDIMA .EQ. 3)THEN
        CALL EXGNV(NTP2EX,ISTP,IZVEL,NODESA,VELZA,IERR)
      END IF
C
      IF(IELMS .NE. 0)THEN
        CALL EXGEV(NTP2EX,ISTP,IELMS,IDBLKA,NUMEBA,EMSSA,IERR)
      ELSE IF (IDENS .NE. 0)THEN
        CALL EXGEV(NTP2EX,ISTP,IDENS,IDBLKA,NUMEBA,DENSA,IERR)
        NNODES = NELNDA
        IF (ITYPE .EQ. 6)NNODES = 4
        DO 100 IEL = 1, NUMEBA
          DO 105 I = 1, NNODES
            XX(I) = XA(ICONA(I,IEL))
            YY(I) = YA(ICONA(I,IEL))
            IF (NDIMA .EQ. 3)THEN
              ZZ(I) = ZA(ICONA(I,IEL))
            ELSE
              ZZ(I) = 0.
            END IF
 105      CONTINUE
          CALL VOL(ITYPE,XX,YY,ZZ,VOLUME)
          EMSSA(IEL) = DENSA(IEL) * VOLUME
 100    CONTINUE
      ELSE
        CALL ERROR('MKEI','NEITHER ELMASS NOR DENSITY AVAILABLE',
     &             'CANNOT COMPUTE MOMENTA AND KINETIC ENERGY',
     &             0,' ',0,' ',0,1)
      END IF
      IF (ISXX .NE. 0)THEN
        CALL EXGEV(NTP2EX,ISTP,ISXX,IDBLKA,NUMEBA,SIGXXA,IERR)
        CALL EXGEV(NTP2EX,ISTP,ISYY,IDBLKA,NUMEBA,SIGYYA,IERR)
        CALL EXGEV(NTP2EX,ISTP,ISZZ,IDBLKA,NUMEBA,SIGZZA,IERR)
        CALL EXGEV(NTP2EX,ISTP,ISXY,IDBLKA,NUMEBA,SIGXYA,IERR)
      END IF
      IF (NDIMA .EQ. 3 .AND. ISYZ .NE.0)THEN
        CALL EXGEV(NTP2EX,ISTP,ISYZ,IDBLKA,NUMEBA,SIGYZA,IERR)
        CALL EXGEV(NTP2EX,ISTP,ISZX,IDBLKA,NUMEBA,SIGZXA,IERR)
      END IF
C
      CALL MKE(NELNDA,NUMEBA,NUMNDA,ICONA,NDLSTA,ITYPE,
     &         VELXA,VELYA,VELZA,EMSSA,RNMSA,
     &         RMXA,RMYA,RMZA,RKEA,PSQA,RJ2A,
     &         SIGXXA,SIGYYA,SIGZZA,SIGXYA,SIGYZA,SIGZXA)
C
      IF (ICOMPL .EQ. 1)THEN
        TMXA(IST)  = TMXA(IST)  + RMXA
        TMYA(IST)  = TMYA(IST)  + RMYA
        TMZA(IST)  = TMZA(IST)  + RMYA
        TKEA(IST)  = TKEA(IST)  + RKEA
        TPSQA(IST) = TPSQA(IST) + PSQA
        TJ2A(IST)  = TJ2A(IST)  + RJ2A
      END IF
C
C repeat for recipient mesh
C get vel's and elmass for mesh-B
C
      CALL EXGNV(NTP4EX,IST,IXVEL,NODESB,VELXB,IERR)
      CALL EXGNV(NTP4EX,IST,IYVEL,NODESB,VELYB,IERR)
      IF (NDIMA .EQ. 3)THEN
        CALL EXGNV(NTP4EX,IST,IZVEL,NODESB,VELZB,IERR)
      END IF
C
      IF(IELMS .NE. 0)THEN
        CALL EXGEV(NTP4EX,IST,IELMS,IDBLKB,NUMEBB,EMSSB,IERR)
      ELSE IF (IDENS .NE. 0)THEN
        CALL EXGEV(NTP4EX,IST,IDENS,IDBLKB,NUMEBB,DENSB,IERR)
        NNODES = NELNDB
        IF (ITYPE .EQ. 6)NNODES = 4
        DO 200 IEL = 1, NUMEBB
          DO 205 I = 1, NNODES
            XX(I) = XB(ICONB(I,IEL))
            YY(I) = YB(ICONB(I,IEL))
            IF (NDIMB .EQ. 3)THEN
              ZZ(I) = ZB(ICONB(I,IEL))
            ELSE
              ZZ(I) = 0.
            END IF
 205      CONTINUE
          CALL VOL(ITYPE,XX,YY,ZZ,VOLUME)
          EMSSB(IEL) = DENSB(IEL) * VOLUME
 200    CONTINUE
      ELSE
        CALL ERROR('MKEI','NEITHER ELMASS NOR DENSITY AVAILABLE',
     &             'CANNOT COMPUTE MOMENTA AND KINETIC ENERGY',
     &             0,' ',0,' ',0,1)
       END IF 
      IF (ISXX .NE. 0)THEN
        CALL EXGEV(NTP4EX,IST,ISXX,IDBLKB,NUMEBB,SIGXXB,IERR)
        CALL EXGEV(NTP4EX,IST,ISYY,IDBLKB,NUMEBB,SIGYYB,IERR)
        CALL EXGEV(NTP4EX,IST,ISZZ,IDBLKB,NUMEBB,SIGZZB,IERR)
        CALL EXGEV(NTP4EX,IST,ISXY,IDBLKB,NUMEBB,SIGXYB,IERR)
      END IF
      IF (NDIMB .EQ. 3 .AND. ISYZ .NE.0)THEN
        CALL EXGEV(NTP4EX,IST,ISYZ,IDBLKB,NUMEBB,SIGYZB,IERR)
        CALL EXGEV(NTP4EX,IST,ISZX,IDBLKB,NUMEBB,SIGZXB,IERR)
      END IF
C
      CALL MKE(NELNDB,NUMEBB,NUMNDB,ICONB,NDLSTB,ITYPE,
     &         VELXB,VELYB,VELZB,EMSSB,RNMSB,
     &         RMXB,RMYB,RMZB,RKEB,PSQB,RJ2B,
     &         SIGXXB,SIGYYB,SIGZZB,SIGXYB,SIGYZB,SIGZXB)
C
      IF (ICOMPL .EQ. 1)THEN
        TMXB(IST) = TMXB(IST) + RMXB
        TMYB(IST) = TMYB(IST) + RMYB
        TMZB(IST) = TMZB(IST) + RMYB
        TKEB(IST) = TKEB(IST) + RKEB
        TPSQB(IST) = TPSQB(IST) + PSQB
        TJ2B(IST)  = TJ2B(IST)  + RJ2B
C
C write stuff out here
C
        WRITE(NTPOUT,1010)IDBLKA,IDBLKB
 1010   FORMAT(/,5X,'DONOR MESH ID ',I7,5X,'RECIPIENT MESH ID ',I7)
        WRITE(NTPOUT,1020)TIMES(ISTP)
 1020   FORMAT(5X,'TIME ',1PE16.9)
        WRITE(NTPOUT,1030)TMXA(IST),TMXB(IST)
 1030   FORMAT(5X,'X-MOMENTUM,       DONOR ',1PE14.6,5X,
     &            'RECIPIENT ',1PE14.6)
        WRITE(NTPOUT,1040)TMYA(IST),TMYB(IST)
 1040   FORMAT(5X,'Y-MOMENTUM,       DONOR ',1PE14.6,5X,
     &            'RECIPIENT ',1PE14.6)
        WRITE(NTPOUT,1050)TMZA(IST),TMZB(IST)
 1050   FORMAT(5X,'Z-MOMENTUM,       DONOR ',1PE14.6,5X,
     &            'RECIPIENT ',1PE14.6)
        WRITE(NTPOUT,1060)TKEA(IST),TKEB(IST)
 1060   FORMAT(5X,'KINETIC ENERGY,   DONOR ',1PE14.6,5X,
     &            'RECIPIENT ',1PE14.6)
        WRITE(NTPOUT,1070)TPSQA(IST),TPSQB(IST)
 1070   FORMAT(5X,'PRESSURE SQUARED, DONOR ',1PE14.6,5X,
     &            'RECIPIENT ',1PE14.6)
        WRITE(NTPOUT,1080)TJ2A(IST),TJ2B(IST)
 1080   FORMAT(5X,'J2,               DONOR ',1PE14.6,5X,
     &            'RECIPIENT ',1PE14.6)
      END IF
      RETURN
      END
