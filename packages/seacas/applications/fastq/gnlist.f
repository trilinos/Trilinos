C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GNLIST (MXLIST, NNUID, MSC, NPNODE, NPELEM, MAXKXN,
     &   NNXK, KXN, NXK, NUID, XN, YN, LIST, NLIST, NUMBER, KCRD, NNN,
     &   ERR, NOROOM)
C***********************************************************************

C  SUBRUOTINE GNLIST = GETS INITIAL NODE LIST TO BEGIN CUTHILL - MCKEE
C                      PROCESS

C***********************************************************************

C  NOTE:
C     AS MANY CARDS OF AS MANY TYPES AS DESIRED MAY BE USED IN
C     ANY ORDER.  IF A NODE IS REFERENCED MORE THAN ONCE A WARNING
C     WILL BE PRINTED AND ONLY THE FIRST REFERENCE WILL BE USED.
C      (IT MAY BE NECESSARY TO MULTIPLY REFERENCE A NODE IN THE
C     CASE OF MORE THAN ONE P - L - P CARD.)

C***********************************************************************

C     EXAMPLE INPUT CARDS
C     COL.1    5   ETC.
C         X-Y              3.5       4.0
C         NODE               7
C         NODE               7  100100002  100100003  100100004    8
C         P-L-P         1    1    2
C         P-L-P         1   77    3   66    5

C***********************************************************************

      DIMENSION LIST (MXLIST), XN (NPNODE), YN (NPNODE)
      DIMENSION KXN (NNXK, MAXKXN), NXK (NNXK, NPELEM), NUID (NNUID)

      CHARACTER*80 NUMBER (MSC)

      LOGICAL ERR, NOROOM

C  INITIALIZE

      ERR=.FALSE.
      NLIST=0

C  NEXT DATA CARD

      DO 150 K=1, KCRD

C  X - Y

         IF ( (NUMBER (K) (1:3) .EQ. 'X-Y') .OR.
     &      (NUMBER (K) (1:3) .EQ. 'x-y')) THEN
            READ (NUMBER (K) (11:20), ' (E10.0)')XVAL
            READ (NUMBER (K) (21:30), ' (E10.0)')YVAL
            INEAR=1
            DIST=SQRT ( (XN (1) - XVAL)**2 +  (YN (1) - YVAL)**2)
            DO 100 I=2, NNN
               D=SQRT ( (XN (I) - XVAL)**2 +  (YN (I) - YVAL)**2)
               IF (D .LT. DIST) THEN
                  DIST=D
                  INEAR=I
               ENDIF
  100       CONTINUE
            NNOW=MAX0 (NLIST, 1)
            IF (IOCCUR (NNOW, LIST, INEAR) .EQ. 1) THEN
               WRITE (*, 10000)NUID (INEAR)
            ELSE
               NLIST=NLIST + 1
               LIST (NLIST)=INEAR
            ENDIF

C  NODE ID

         ELSEIF ( (NUMBER (K) (1:3) .EQ. 'NOD') .OR.
     &      (NUMBER (K) (1:3) .EQ. 'nod')) THEN
            DO 110 I=11, 71, 10
               J=I + 9
               READ (NUMBER (K) (I:J), ' (I10)')IVAL
               IF (IVAL .GT. 0) THEN
                  NEW=INDX (NNN, NUID, IVAL)
                  IF (NEW .EQ. 0) THEN
                     WRITE (*, 10010)IVAL
                     CALL MESAGE ('THIS NODE WILL BE SKIPPED')
                  ELSEIF (IOCCUR (NLIST, LIST, NEW) .EQ. 0) THEN
                     NLIST=NLIST + 1
                     LIST (NLIST)=NEW
                  ELSE
                     WRITE (*, 10000)IVAL
                  ENDIF
               ENDIF
  110       CONTINUE

C  P-L-P

         ELSEIF ( (NUMBER (K) (1:3) .EQ. 'P-L') .OR.
     &      (NUMBER (K) (1:3) .EQ. 'p-l')) THEN
            NUMNEW=0
            DO 130 J=6, 66, 10
               READ (NUMBER (K) (J:J + 4), ' (I5)')IP1
               READ (NUMBER (K) (J + 5:J + 9), ' (I5)')LINE
               READ (NUMBER (K) (J + 10:J + 14), ' (I5)')IP2
               IF (IP1 .GT. 0) THEN
                  MXLST1=MXLIST - NLIST
                  CALL GETPLP (NPNODE, NPELEM, MAXKXN, NNXK, MXLST1,
     &               KXN, NXK, NUID, IP1, LINE, IP2, LIST (NLIST + 1),
     &               NUMNEW, NNN, LASTN, NOROOM, ERR)
                  IF (NOROOM) THEN
                     CALL MESAGE ('DIMENSIONS MUST BE INCREASED')
                     RETURN
                  ELSEIF (ERR) THEN
                     RETURN
                  ENDIF
                  NLIST1=NLIST + 1
                  NLISTN=NLIST + NUMNEW
                  IF (NUMNEW .GT. 0) THEN
                     DO 120 I=NLIST1, NLISTN
                        IF ( (NLIST .EQ. 0) .OR.
     &                     (IOCCUR (NLIST, LIST, LIST (I)) .EQ. 0)) THEN
                           NLIST=NLIST + 1
                           LIST (NLIST)=LIST (I)
                        ELSE
                           WRITE (*, 10000)LIST (I)
                        ENDIF
  120                CONTINUE
                  ENDIF
               ELSE
                  GOTO 140
               ENDIF
  130       CONTINUE
  140       CONTINUE
         ENDIF
  150 CONTINUE

      RETURN

10000 FORMAT (' NODE', I10, ' IS ALREADY IN THE LIST')
10010 FORMAT (' NODE', I10,
     &   ' IS NOT AN IDENTIFIER OF A NODE IN THIS MESH')

      END
