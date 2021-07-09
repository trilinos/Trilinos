C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE NSPLIT (MXND, MLN, LNODES, ANGLE, NSTART, KANG,
     &   INODE, NNODE, NWANT, MAXSIZ)
C***********************************************************************

C  SUBROUTINE NSPLIT = SPLITS UP THE KANG CONSECUTIVE NODES STARTING
C                      AT NSTART INTO NWANT INTERVALS (OR AS CLOSE
C                      AS AN BE EXPECTED).  THE MAXIMUM NWANT SHOULD
C                      BE IS 4.

C***********************************************************************

      DIMENSION LNODES (MLN, MXND), ANGLE (MXND), INODE (4)

      LOGICAL MAXSIZ

      NNODE = 0

      IF (KANG .LE. NWANT) THEN
         NNOW = NSTART
         DO 100 I = 1, KANG
            INODE (I) = NNOW
            NNOW = LNODES (3, NNOW)
  100    CONTINUE
         NNODE = KANG

      ELSEIF (NWANT .EQ. 1) THEN
         NNODE = 1
         IF (KANG .EQ. 2) THEN
            IF (MAXSIZ) THEN
               IF (ANGLE (NSTART) .GT. ANGLE (LNODES (3, NSTART)) ) THEN
                  INODE (1) = NSTART
               ELSE
                  INODE (1) = LNODES (3, NSTART)
               ENDIF
            ELSE
               IF (ANGLE (NSTART) .GT. ANGLE (LNODES (3, NSTART)) ) THEN
                  INODE (1) = LNODES (3, NSTART)
               ELSE
                  INODE (1) = NSTART
               ENDIF
            ENDIF
         ELSE
            INODE (1) = JUMPLP (MXND, MLN, LNODES, NSTART, KANG / 2)
         ENDIF

      ELSEIF (NWANT .EQ. 2) THEN
         NNODE = 2
         NJUMP = NINT (DBLE(KANG + 1) / 4.)
         INODE (1) = JUMPLP (MXND, MLN, LNODES, NSTART,
     &      NJUMP - 1)
         INODE (2) = JUMPLP (MXND, MLN, LNODES, NSTART,
     &      KANG - NJUMP)

      ELSEIF (NWANT .EQ. 3) THEN
         NNODE = 3
         NJUMP1 = NINT (DBLE(KANG + 1) / 6.)
         NJUMP2 = NINT (DBLE(KANG + 1) / 2.)
         INODE (1) = JUMPLP (MXND, MLN, LNODES, NSTART,
     &      NJUMP1 - 1)
         INODE (2) = JUMPLP (MXND, MLN, LNODES, NSTART,
     &      NJUMP2 - 1)
         INODE (3) = JUMPLP (MXND, MLN, LNODES, NSTART,
     &      KANG - NJUMP1)

      ELSEIF (NWANT .EQ. 4) THEN
         NNODE = 4
         XKANG = KANG + 1
         NJUMP1 = NINT (XKANG / 8.) - 1
         NJUMP2 = NINT (XKANG / 2.) - NINT (XKANG / 8.) -1
         NJUMP3 = NINT (XKANG / 2.) + NINT (XKANG / 8.) -1
         INODE (1) = JUMPLP (MXND, MLN, LNODES, NSTART,
     &      NJUMP1)
         INODE (2) = JUMPLP (MXND, MLN, LNODES, NSTART,
     &      NJUMP2)
         INODE (3) = JUMPLP (MXND, MLN, LNODES, NSTART,
     &      NJUMP3)
         INODE (4) = JUMPLP (MXND, MLN, LNODES, NSTART,
     &      KANG - NJUMP1 - 1)
      ENDIF

      RETURN

      END
