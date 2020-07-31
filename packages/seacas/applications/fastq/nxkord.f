C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE NXKORD (NODES, N1)
C***********************************************************************

C  SUBROUTINE NXKORD = ROTATES THE LIST OF FOUR NODES SO N1 APPEARS
C                      FIRST IF IT IS IN THE LIST

C***********************************************************************

      DIMENSION NODES (4)

      DO 100 I = 1, 4
         IF (NODES (I) .EQ. N1) THEN
            IF (I .EQ. 1) THEN
               RETURN
            ELSEIF (I .EQ. 2) THEN
               NSAVE = NODES (1)
               NODES (1) = NODES (2)
               NODES (2) = NODES (3)
               NODES (3) = NODES (4)
               NODES (4) = NSAVE
            ELSEIF (I .EQ. 3) THEN
               NSAVE = NODES (1)
               NSAVE2 = NODES (2)
               NODES (1) = NODES (3)
               NODES (2) = NODES (4)
               NODES (3) = NSAVE
               NODES (4) = NSAVE2
            ELSE
               NSAVE = NODES (4)
               NODES (4) = NODES (3)
               NODES (3) = NODES (2)
               NODES (2) = NODES (1)
               NODES (1) = NSAVE
            ENDIF
            RETURN
         ENDIF
  100 CONTINUE

      RETURN

      END
