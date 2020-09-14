C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SAVREG (MXND, MAXNBC, MAXSBC, XN, YN, NUID, LXK, NXL,
     &   LXN, LSTNBC, LSTSBC, KNBC, KSBC, NNN, KKK, NUMREG, IUNIT, BAR,
     &   M1)
C***********************************************************************

C  SUBROUTINE SAVREG = SAVES THE NODE AND ELEMENT DESCRIPTIONS AS WELL
C                      AS THE BOUNDARY CONDITIONS

C***********************************************************************

C  NOTE:
C     THE MESH TABLES ARE EFFECTIVELY DESTROYED BY THIS ROUTINE

C***********************************************************************

      DIMENSION XN (MXND), YN (MXND), NUID (MXND)
      DIMENSION LXK (4, MXND), NXL (2, MXND*3), LXN (4, MXND)
      DIMENSION LSTNBC (MAXNBC), LSTSBC (MAXSBC), NODES (4)

      LOGICAL CCW, BAR

      CCW = .TRUE.
      IF (.NOT.BAR) THEN

C  DEFINE NUID-S FOR INTERIOR NODES.
C  SKIP DELETED NODES AND CONTINUATIONS.

         K = 0
         DO 100 I = 1, NNN
            IF ((NUID (I) .EQ. 0) .AND. (LXN (1, I) .GT. 0)) THEN
               K = K+1
               NUID (I) = NUMREG*100000+K
            ENDIF
  100    CONTINUE

C  GET COUNTER-CLOCKWISE NODE LISTS.
C   (THESE LISTS WILL OVERWRITE THE LXK ARRAY.)
C  DELETED ELEMENTS WILL BE SKIPPED.

         J = 0
         IDEL = 0
         DO 130 K = 1, KKK
            IF (LXK (1, K) .LE. 0) THEN
               DO 110 JJ = 2, KSBC, 3
                  IF (LSTSBC (JJ) .GE. (K - IDEL)) THEN
                     LSTSBC (JJ) = LSTSBC (JJ) - 1
                  ENDIF
  110          CONTINUE
               IDEL = IDEL + 1
            ELSE
               CALL GNXKA (MXND, XN, YN, K, NODES, AREA, LXK, NXL, CCW)
               J = J+1
               DO 120 I = 1, 4
                  N = NODES (I)
                  LXK (I, J) = IABS (NUID (N))
  120          CONTINUE
            ENDIF
  130    CONTINUE
         KKK = J
      ELSE
         DO 140 I = 1, KKK
            LXK (1, I) = IABS (NUID (LXK (1, I)))
            LXK (2, I) = IABS (NUID (LXK (2, I)))
  140    CONTINUE
      ENDIF

C  COLLAPSE THE NODE ARRAYS TO ELIMINATE DELETED NODES,
C  CONTINUATIONS,  AND NODES ALREADY WRITTEN OUT.

      K = 0
      DO 150 I = 1, NNN
         IF ( ( (LXN (1, I) .GT. 0) .OR. (BAR))
     &      .AND. (NUID (I) .GT. 0) ) THEN
            K = K+1
            XN (K) = XN (I)
            YN (K) = YN (I)
            NUID (K) = NUID (I)
         ENDIF
  150 CONTINUE
      NNN = K

C  WRITE HEADER,  NODE LIST,  ELEMENT LIST,  AND BOUNDARY LISTS

      WRITE (IUNIT)KKK, NNN, KNBC, KSBC, NUMREG, BAR, M1
      IF (NNN .GE. 1) WRITE (IUNIT) (NUID (I), XN (I), YN (I),
     &   I = 1, NNN)
      WRITE (IUNIT) ((LXK (I, J), I = 1, 4), J = 1, KKK)
      IF (KNBC .GT. 0)WRITE (IUNIT) (LSTNBC (I), I = 1, KNBC)
      IF (KSBC .GT. 0)WRITE (IUNIT) (LSTSBC (I), I = 1, KSBC)

      RETURN

      END
