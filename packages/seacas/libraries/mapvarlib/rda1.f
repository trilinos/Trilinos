C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
*DECK,RDA1
      SUBROUTINE RDA1(XA,YA,ZA,DISXA,DISYA,DISZA)

C     ******************************************************************

C     SUBROUTINE TO EXTRACTREAD THE CRITICAL INPUT AND SIZING PARAMETERS
C     FROM THE GENESIS FILE FOR MESH-A

C     READS MESH A, WRITES MESH C DATA AS APPROPRIATE

C     Calls function LENSTR
C     Calls subroutine ERROR

C     Called by MAPVAR

C     ******************************************************************

C  XA,etc     REAL  Coordinates of mesh-A nodes (1:nodesa)
C  DISXA,etc  REAL  Displacements of mesh-A nodes (1:nodesa)

C     ******************************************************************

      include 'aexds1.blk'
      include 'aexds2.blk'
      include 'amesh.blk'
      include 'contrl.blk'
      include 'ex2tp.blk'
      include 'rundat.blk'
      include 'steps.blk'
      include 'varnpt.blk'
      include 'varept.blk'

      DIMENSION xa(*),ya(*),za(*)
      DIMENSION DISXA(*),DISYA(*),DISZA(*)

      PARAMETER (MAXQA=240)

C     ******************************************************************

C nodal point coordinates and names

      CALL EXGCON (NTP2EX,NAMECO,IERR)

C Convert to upper case

      DO 10 I = 1, NDIMA
        CALL EXUPCS(NAMECO(I))
 10   CONTINUE

      CALL EXPCON (NTP4EX,NAMECO,IERR)
      CALL EXGCOR (NTP2EX,XA,YA,ZA,IERR)

C QA

      NQAREC = EXINQI(NTP2EX,EXQA)

      IF (NQAREC.GT.MAXQA) THEN
        CALL ERROR ('RDA1','TOO MANY QA RECORDS IN
     1MESH-A FILE','NO. RECORDS',NQAREC,'NO. RECORDS ALLOWED',240,
     2'ONLY LAST 239 QA RECORDS RETAINED',
     3'IF NOT ACCEPTABLE SEE CODE SPONSOR TO INCREASE ARRAY QALINE',0)
        NQAREC = MAXQA
      END IF

      CALL EXGQA (NTP2EX,QALINE,IERR)
      IF (NQAREC .EQ. MAXQA)THEN
        DO IQ = 1,4
           DO NQ = 2,MAXQA
              QALINE(IQ,NQ) = QALINE(IQ,NQ-1)
           end do
        end do
      ELSE
        NQAREC = NQAREC + 1
      END IF

      QALINE(1,NQAREC) = QAINFO(1)
      QALINE(2,NQAREC) = QAINFO(3)
      QALINE(3,NQAREC) = QAINFO(5)
      QALINE(4,NQAREC) = QAINFO(6)
      CALL EXPQA (NTP4EX,NQAREC,QALINE,IERR)

C VARIABLE NAMES

      CALL EXGVP (NTP2EX,"G",NVARGP,IERR)

C Do some error checking on number of variables - got me once

      NUMNAM = NVARGP
      IF (NUMNAM .GT. MXVAR)CALL ERROR('RDA1','TOO MANY VARIABLE NAMES
     1IN MESH-A DATA BASE','NUMBER OF VARIABLE NAMES ENCOUNTERED SO
     2FAR',NUMNAM,'NUMBER ALLOWED - FIXED DIMENSION',MXVAR,'SEE CODE
     3SPONSOR FOR INCREASE IN --NAMVAR--',' ',1)

      if (nvargp .gt. 0) then
         CALL EXGVAN (NTP2EX,"G",NVARGP,NAMVAR,IERR)

C     Convert to upper case

         DO I = 1, NVARGP
            CALL EXUPCS(NAMVAR(I))
         end do

         CALL EXPVP (NTP4EX,"G",NVARGP,IERR)
         CALL EXPVAN (NTP4EX,"G",NVARGP,NAMVAR,IERR)
      end if

      CALL EXGVP (NTP2EX,"E",NVAREL,IERR)

      NUMNAM = NUMNAM + NVAREL
      IF (NUMNAM .GT. MXVAR)CALL ERROR('RDA1','TOO MANY VARIABLE NAMES
     1IN MESH-A DATA BASE','NUMBER OF VARIABLE NAMES ENCOUNTERED SO
     2FAR',NUMNAM,'NUMBER ALLOWED - FIXED DIMENSION',MXVAR,'SEE CODE
     3SPONSOR FOR INCREASE IN --NAMVAR--',' ',1)

      if (nvarel .gt. 0) then
         CALL EXGVAN (NTP2EX,"E",NVAREL,NAMVAR(NVARGP+1),IERR)

C     Convert to upper case

         DO I = 1, NVAREL
            CALL EXUPCS(NAMVAR(NVARGP+I))
         end do

         CALL EXPVP (NTP4EX,"E",NVAREL,IERR)
         CALL EXPVAN (NTP4EX,"E",NVAREL,NAMVAR(NVARGP+1),IERR)
      end if

      CALL EXGVP (NTP2EX,"N",NVARNP,IERR)

      NUMNAM = NUMNAM + NVARNP
      IF (NUMNAM .GT. MXVAR)CALL ERROR('RDA1','TOO MANY VARIABLE NAMES
     1IN MESH-A DATA BASE','NUMBER OF VARIABLE NAMES ENCOUNTERED SO
     2FAR',NUMNAM,'NUMBER ALLOWED - FIXED DIMENSION',MXVAR,'SEE CODE
     3SPONSOR FOR INCREASE IN --NAMVAR--',' ',1)

      if (nvarnp .gt. 0) then
         CALL EXGVAN (NTP2EX,"N",NVARNP,NAMVAR(NVARGP+NVAREL+1),IERR)

C     Convert to upper case

         DO I = 1, NVARNP
            CALL EXUPCS(NAMVAR(NVARGP+NVAREL+I))
         end do

         CALL EXPVP (NTP4EX,"N",NVARNP,IERR)
         CALL EXPVAN (NTP4EX,"N",NVARNP,NAMVAR(NVARGP+NVAREL+1),IERR)
      end if

      LC1 = LENSTR (NAMECO(1))
      LC2 = LENSTR (NAMECO(2))
      LC3 = 2
      IF (NDIMA .EQ. 3) LC3 = LENSTR (NAMECO(3))
      ISTART = 1+NVARGP+NVAREL
      DO I = ISTART, NUMNAM
        LN = LENSTR (NAMVAR(I))
        IF (NAMVAR(I)(1:1) .EQ. 'D')THEN
          IF(NAMVAR(I)(LN:LN) .EQ. NAMECO(1)(LC1:LC1))
     &                             IXDIS=I-ISTART+1
          IF(NAMVAR(I)(LN:LN) .EQ. NAMECO(2)(LC2:LC2))
     &                             IYDIS=I-ISTART+1
          IF (NDIMA .GE. 3)THEN
            IF(NAMVAR(I)(LN:LN) .EQ. NAMECO(3)(LC3:LC3))
     &                               IZDIS=I-ISTART+1
          END IF
        END IF
      end do

      IF (IDEF .NE. 0 .AND. IXDIS .NE. 0 .AND. IYDIS .NE. 0)THEN

C  Work in deformed coordinates

        CALL EXGNV (NTP2EX,ISTEP,IXDIS,NODESA,DISXA,IERR)
        CALL EXGNV (NTP2EX,ISTEP,IYDIS,NODESA,DISYA,IERR)
        IF (NDIMA .GE. 3) THEN
          CALL EXGNV (NTP2EX,ISTEP,IZDIS,NODESA,DISZA,IERR)
        END IF
        DO I = 1, NODESA
          XA(I) = XA(I) + DISXA(I)
          YA(I) = YA(I) + DISYA(I)
          IF (NDIMA .GE. 3) ZA(I) = ZA(I) + DISZA(I)
       end do
      ELSE

C No displacements in Mesh-A data, can't do deformed processing

        IDEF = 0
      END IF
      IF (IACCU .EQ. 1)THEN

C ********************************************************************
C accuracy checK
C ********************************************************************

C find needed variables
C 1st velocities
C coordinate names - velocity will start with "v"
C                    and end with last character
C                    of coordinate name

        LC1 = LENSTR(NAMECO(1))
        LC2 = LENSTR(NAMECO(2))
        LC3 = 2
        IF (NDIMA .EQ. 3)LC3 = LENSTR(NAMECO(3))

C search nodal variables, get ptrs to vel's and elmass if available

        IXVEL = 0
        IYVEL = 0
        IZVEL = 0
        IELMS = 0
        IDENS = 0
        ISTART = 1+NVARGP+NVAREL
        DO 90 INAM = ISTART, NUMNAM
          IF (NAMVAR(INAM)(1:1) .EQ. 'V')THEN
            LN = LENSTR(NAMVAR(INAM))
            IF (NAMVAR(INAM)(LN:LN) .EQ. NAMECO(1)(LC1:LC1))THEN
              IXVEL = INAM - ISTART + 1
              GO TO 90
            END IF
            IF (NAMVAR(INAM)(LN:LN) .EQ. NAMECO(2)(LC2:LC2))THEN
              IYVEL = INAM - ISTART + 1
              GO TO 90
            END IF
            IF (NDIMA .EQ. 3)THEN
              IF (NAMVAR(INAM)(LN:LN) .EQ. NAMECO(3)(LC3:LC3))THEN
                IZVEL = INAM - ISTART + 1
              END IF
            END IF
          END IF
   90   CONTINUE
        ISTART = 1+NVARGP
        IEND = NVARGP+NVAREL
        DO 100 INAM = ISTART, IEND
          IF (NAMVAR(INAM)(1:6) .EQ. 'ELMASS')THEN
            IELMS = INAM - ISTART + 1
            GO TO 100
          ELSE IF (NAMVAR(INAM)(1:4) .EQ. 'DENS')THEN
            IDENS = INAM - ISTART + 1
            GO TO 100
          ELSE IF (NAMVAR(INAM)(1:3) .EQ. 'SIG')THEN
            LN = LENSTR(NAMVAR(INAM))
            IF (NAMVAR(INAM)(LN-1:LN) .EQ. 'XX')THEN
              ISXX = INAM - ISTART +1
              GO TO 100
            ELSE IF (NAMVAR(INAM)(LN-1:LN) .EQ. 'YY')THEN
              ISYY = INAM - ISTART +1
              GO TO 100
            ELSE IF (NAMVAR(INAM)(LN-1:LN) .EQ. 'ZZ')THEN
              ISZZ = INAM - ISTART +1
              GO TO 100
            ELSE IF (NAMVAR(INAM)(LN-1:LN) .EQ. 'XY')THEN
              ISXY = INAM - ISTART +1
              GO TO 100
            ELSE IF (NAMVAR(INAM)(LN-1:LN) .EQ. 'YZ')THEN
              ISYZ = INAM - ISTART +1
              GO TO 100
            ELSE IF (NAMVAR(INAM)(LN-1:LN) .EQ. 'ZX')THEN
              ISZX = INAM - ISTART +1
              GO TO 100
            END IF
          ELSE IF (NAMVAR(INAM)(1:4) .EQ. 'USIG')THEN
            LN = LENSTR(NAMVAR(INAM))
            IF (NAMVAR(INAM)(LN-1:LN) .EQ. 'XX')THEN
              ISXX = INAM - ISTART +1
              GO TO 100
            ELSE IF (NAMVAR(INAM)(LN-1:LN) .EQ. 'YY')THEN
              ISYY = INAM - ISTART +1
              GO TO 100
            ELSE IF (NAMVAR(INAM)(LN-1:LN) .EQ. 'ZZ')THEN
              ISZZ = INAM - ISTART +1
              GO TO 100
            ELSE IF (NAMVAR(INAM)(LN-1:LN) .EQ. 'XY')THEN
              ISXY = INAM - ISTART +1
              GO TO 100
            ELSE IF (NAMVAR(INAM)(LN-1:LN) .EQ. 'YZ')THEN
              ISYZ = INAM - ISTART +1
              GO TO 100
            ELSE IF (NAMVAR(INAM)(LN-1:LN) .EQ. 'ZX')THEN
              ISZX = INAM - ISTART +1
              GO TO 100
            END IF
          END IF
  100   CONTINUE
      END IF

      RETURN
      END
