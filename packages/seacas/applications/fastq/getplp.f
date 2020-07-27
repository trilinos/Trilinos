C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GETPLP (NPNODE, NPELEM, MAXKXN, NNXK, MXLIST, KXN,
     &   NXK, NUID, IP1, LINE, IP2, LIST, NLIST, NNN, LASTN, NOROOM,
     &   ERR)
C***********************************************************************

C  SUBROUTINE GETPLP = PRODUCES THE LIST OF NODES FROM POINT IP1
C                      THROUGH LINE TO POINT IP2

C***********************************************************************

C  NOTE:
C     THIS LIST WILL BE (LIST (I), I=1,NLIST) AND THESE WILL BE INDICES
C     INTO THE NODE TABLE

C***********************************************************************

      DIMENSION NXNLST (20)
      DIMENSION KXN (NNXK, MAXKXN), NXK (NNXK, NPELEM), NUID (NPNODE)
      DIMENSION LIST (MXLIST)

      LOGICAL ERR, ALL, NOROOM

      ERR = .FALSE.
      NOROOM = .FALSE.

C  FIND FIRST POINT

      IF (NLIST .EQ. 0) THEN
         N = INDX (NNN, NUID, IP1)
         IF (N .LE. 0) THEN
            WRITE ( * , 10000)IP1
            ERR = .TRUE.
            RETURN
         ENDIF
         NLIST = NLIST + 1
         IF (NLIST .GT. MXLIST) THEN
            NOROOM = .TRUE.
            RETURN
         ENDIF
         LIST (NLIST) = N
      ELSE
         NLIST = 0
         N = LASTN
      ENDIF

C  FOLLOW THE LINE

      IF (LINE .LE. 0) RETURN
      NPREV = 0
  100 CONTINUE
      ALL = .FALSE.
      CALL GETNXN (NPNODE, NPELEM, MAXKXN, NNXK, KXN, NXK, NUID, N,
     &   NXNLST, NUMN, ALL, ERR)
      IF (ERR) RETURN
      DO 110 I = 1, NUMN
         NEW = NXNLST (I)
         IF ( (NEW .NE. NPREV) .AND. (NUID (NEW) .GE. 1000000000)) THEN
            L =  (NUID (NEW) - 1000000000) / 100000
            IF (L .EQ. LINE) THEN
               NLIST = NLIST + 1
               IF (NLIST .GT. MXLIST) THEN
                  NOROOM = .TRUE.
                  RETURN
               ENDIF
               LIST (NLIST) = NEW
               NPREV = N
               N = NEW
               GOTO 100
            ENDIF
         ENDIF
  110 CONTINUE

C  LINE FINISHED  -  FIND IP2

      IF (IP2 .LE. 0) RETURN
      DO 120 I = 1, NUMN
         NEW = NXNLST (I)
         IF (NUID (NEW) .EQ. IP2) THEN
            NLIST = NLIST + 1
            IF (NLIST .GT. MXLIST) THEN
               NOROOM = .TRUE.
               RETURN
            ENDIF
            LIST (NLIST) = NEW
            LASTN = NEW
            RETURN
         ENDIF
  120 CONTINUE

C  LINE DID NOT MATCH UP RIGHT

      WRITE ( * , 10010)IP1, LINE, IP2
      ERR = .TRUE.
      RETURN

10000 FORMAT (' POINT', I5, ' IS NOT IN THE MESH')
10010 FORMAT (' P-L-P SEQUENCE OF', I5, ' -', I5, ' -', I5,
     &   'IS AN ILLEGAL SEQUENCE')
      END
