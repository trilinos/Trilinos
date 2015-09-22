C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Government retains certain rights in this software.
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

C $Id: adjmsh.f,v 1.2 1999/06/21 22:43:40 gdsjaar Exp $
C $Log: adjmsh.f,v $
C Revision 1.2  1999/06/21 22:43:40  gdsjaar
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
C Revision 1.1.1.1  1990/11/30 11:03:22  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:03:20  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]ADJMSH.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE ADJMSH (MS, MR, NPNODE, NPELEM, MXNFLG, MXSFLG, NPREGN,
     &   NPNBC, NPSBC, MCOM, ICOM, JCOM, CIN, RIN, IIN, KIN,
     &   NNN, KKK, NNXK, NODES, NELEMS, NNFLG, NNPTR, NNLEN, NSFLG,
     &   NSPTR, NSLEN, NVPTR, NVLEN, NSIDEN, MAPDXG, XN, YN, NXK, MAT,
     &   MAPGXD, MATMAP, WTNODE, WTSIDE, NBCNOD, NNLIST, NBCSID, NSLIST,
     &   NVLIST, NUMMAT, LINKM, TITLE, ERR, EIGHT, NINE, VERSN)
C***********************************************************************
C
C  SUBROUTINE ADJMSH = ADJUSTS A GENESIS DATABASE OUTPUT
C
C***********************************************************************
C
      DIMENSION XN (NPNODE), YN (NPNODE), NXK (NNXK, NPELEM)
      DIMENSION MAT (NPELEM)
      DIMENSION NODES (NPNBC), NELEMS (NPSBC), NSIDEN (NPSBC)
      DIMENSION NNFLG (MXNFLG), NNLEN (MXNFLG)
      DIMENSION NNPTR (MXNFLG), WTNODE (NPNBC)
      DIMENSION NSFLG (MXSFLG), NSLEN (MXSFLG)
      DIMENSION NSPTR (MXSFLG), WTSIDE (NPSBC)
      DIMENSION NVLEN (MXSFLG), NVPTR (MXSFLG), LINKM (2,  (MS+MR))
      DIMENSION MAPDXG (NPNODE), MAPGXD (NPNODE), MATMAP (3, NPREGN)
      DIMENSION KIN (MCOM), IIN (MCOM), RIN (MCOM)
C
      LOGICAL FOUND, ERR
C
      CHARACTER*72 TITLE, CIN (MCOM)
      CHARACTER*10 VERSN
C
      CALL MESAGE (' ')
      CALL MESAGE
     &   ('*********************************************************')
      CALL MESAGE
     &   ('** MESH ADJUST OPTION IS CURRENTLY LIMITED TO DELETING **')
      CALL MESAGE
     &   ('**      ELEMENTS SIDE BOUNDARY FLAGS BY MATERIAL       **')
      CALL MESAGE
     &   ('*********************************************************')
      CALL MESAGE (' ')
C
C  ADJUST SIDE BOUNDARY FLAGS BY MATERIALS
C
      CALL MESAGE ('ENTER DATA IN THE FOLLOWING FORMAT:')
      CALL MESAGE ('[ MATERIAL NUMBER, FLAG ID ]')
      CALL MESAGE ('HIT RETURN TO END INPUT')
  100 CONTINUE
      IF (ICOM .GT. JCOM) THEN
         CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN, IIN,
     &      RIN)
         ICOM = 1
      END IF
      IF ((ICOM .GT. JCOM) .OR. (CIN (ICOM) (1:1) .EQ. ' ')) THEN
         ICOM = ICOM + 1
         GOTO 190
      ELSE
         I1 = IIN (ICOM)
         ICOM = ICOM + 1
         IF ((ICOM .LE. JCOM) .AND. (KIN (ICOM) .GT. 0)) THEN
            I2 = IIN (ICOM)
            ICOM = ICOM + 1
         ELSE
            ICOM = ICOM + 1
            CALL MESAGE ('** NOT ENOUGH INFORMATION IS SUPPLIED **')
            GOTO 100
         ENDIF
      ENDIF
C
C  NOW THAT THE MATERIAL (I1) AND THE FLAG ID (I2) ARE ENTERED
C  FIRST CHECK TO MAKE SURE THAT THAT MATERIAL IS PRESENT
C
      DO 110 I = 1, NUMMAT
         IF (MATMAP (1, I) .EQ. I1) THEN
            J1 = MATMAP (2, I)
            J2 = MATMAP (3, I)
            GOTO 120
         ENDIF
  110 CONTINUE
      CALL MESAGE('** THAT MATERIAL IS NOT PRESENT IN THE MESH **')
      GOTO 100
C
  120 CONTINUE
C
C  NOW FIND THE ELEMENT SIDE FLAG
C
      DO 130 I = 1, NBCSID
         IF (NSFLG (I) .EQ. I2) THEN
            II = I
            GOTO 140
         ENDIF
  130 CONTINUE
      CALL MESAGE ('** THAT ELEMENT BOUNDARY FLAG IS NOT IN THE '//
     &   'MESH **')
      GOTO 100
C
  140 CONTINUE
C
C  NOW SEARCH THE LOOP FOR ELEMENTS ATTACHED TO THAT BOUNDARY FLAG
C  OF THE SPECIFIED MATERIAL
C
      IBEGIN = NSPTR (II)
      IEND = NSPTR (II) + NSLEN (I) - 1
C
      FOUND = .FALSE.
      KOUNT = 0
C
      DO 180 I = IBEGIN, IEND
         IF ((NELEMS (I - KOUNT) .GE. J1) .AND.
     &      (NELEMS (I - KOUNT) .LE. J2)) THEN
C
C  AN ELEMENT SIDE FLAG HAS BEEN FOUND - NOW DELETE IT
C
            FOUND = .TRUE.
C
            DO 150 J = I - KOUNT, NSLIST - 1
               NELEMS (J) = NELEMS (J + 1)
  150       CONTINUE
            NSLIST = NSLIST - 1
C
            DO 160 J = (((I - KOUNT) * 2) -1), NVLIST - 2
               NSIDEN (J) = NSIDEN (J + 2)
               WTSIDE (J) = WTSIDE (J + 2)
  160       CONTINUE
            NVLIST = NVLIST - 2
C
            NSLEN (II) = NSLEN (II) - 1
            NVLEN (II) = NVLEN (II) - 2
            DO 170 J = II + 1, NBCSID
               NSPTR (J) = NSPTR (J) - 1
               NVPTR (J) = NVPTR (J) - 2
  170       CONTINUE
C
            KOUNT = KOUNT + 1
         ENDIF
  180 CONTINUE
      IF (.NOT. FOUND) THEN
         CALL MESAGE ('** NO MATCHES OF ELEMENTS WITH THAT BOUNDARY '//
     &      'FLAG AND MATERIAL **')
      ENDIF
      GOTO 100
C
  190 CONTINUE
      RETURN
C
      END
