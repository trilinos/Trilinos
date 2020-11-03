C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ADJMSH (MS, MR, NPNODE, NPELEM, MXNFLG, MXSFLG, NPREGN,
     &   NPNBC, NPSBC, MCOM, ICOM, JCOM, CIN, RIN, IIN, KIN,
     &   NNN, KKK, NNXK, NODES, NELEMS, NNFLG, NNPTR, NNLEN, NSFLG,
     &   NSPTR, NSLEN, NVPTR, NVLEN, NSIDEN, MAPDXG, XN, YN, NXK, MAT,
     &   MAPGXD, MATMAP, WTNODE, WTSIDE, NBCNOD, NNLIST, NBCSID, NSLIST,
     &   NVLIST, NUMMAT, LINKM, TITLE, ERR, EIGHT, NINE, VERSN)
C***********************************************************************

C  SUBROUTINE ADJMSH = ADJUSTS A GENESIS DATABASE OUTPUT

C***********************************************************************

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

      LOGICAL FOUND, ERR

      CHARACTER*72 TITLE, CIN (MCOM)
      CHARACTER*10 VERSN

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

C  ADJUST SIDE BOUNDARY FLAGS BY MATERIALS

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

C  NOW THAT THE MATERIAL (I1) AND THE FLAG ID (I2) ARE ENTERED
C  FIRST CHECK TO MAKE SURE THAT THAT MATERIAL IS PRESENT

      DO 110 I = 1, NUMMAT
         IF (MATMAP (1, I) .EQ. I1) THEN
            J1 = MATMAP (2, I)
            J2 = MATMAP (3, I)
            GOTO 120
         ENDIF
  110 CONTINUE
      CALL MESAGE('** THAT MATERIAL IS NOT PRESENT IN THE MESH **')
      GOTO 100

  120 CONTINUE

C  NOW FIND THE ELEMENT SIDE FLAG

      DO 130 I = 1, NBCSID
         IF (NSFLG (I) .EQ. I2) THEN
            II = I
            GOTO 140
         ENDIF
  130 CONTINUE
      CALL MESAGE ('** THAT ELEMENT BOUNDARY FLAG IS NOT IN THE '//
     &   'MESH **')
      GOTO 100

  140 CONTINUE

C  NOW SEARCH THE LOOP FOR ELEMENTS ATTACHED TO THAT BOUNDARY FLAG
C  OF THE SPECIFIED MATERIAL

      IBEGIN = NSPTR (II)
      IEND = NSPTR (II) + NSLEN (I) - 1

      FOUND = .FALSE.
      KOUNT = 0

      DO 180 I = IBEGIN, IEND
         IF ((NELEMS (I - KOUNT) .GE. J1) .AND.
     &      (NELEMS (I - KOUNT) .LE. J2)) THEN

C  AN ELEMENT SIDE FLAG HAS BEEN FOUND - NOW DELETE IT

            FOUND = .TRUE.

            DO 150 J = I - KOUNT, NSLIST - 1
               NELEMS (J) = NELEMS (J + 1)
  150       CONTINUE
            NSLIST = NSLIST - 1

            DO 160 J = (((I - KOUNT) * 2) -1), NVLIST - 2
               NSIDEN (J) = NSIDEN (J + 2)
               WTSIDE (J) = WTSIDE (J + 2)
  160       CONTINUE
            NVLIST = NVLIST - 2

            NSLEN (II) = NSLEN (II) - 1
            NVLEN (II) = NVLEN (II) - 2
            DO 170 J = II + 1, NBCSID
               NSPTR (J) = NSPTR (J) - 1
               NVPTR (J) = NVPTR (J) - 2
  170       CONTINUE

            KOUNT = KOUNT + 1
         ENDIF
  180 CONTINUE
      IF (.NOT. FOUND) THEN
         CALL MESAGE ('** NO MATCHES OF ELEMENTS WITH THAT BOUNDARY '//
     &      'FLAG AND MATERIAL **')
      ENDIF
      GOTO 100

  190 CONTINUE
      RETURN

      END
