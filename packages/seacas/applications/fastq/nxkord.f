C $Id: nxkord.f,v 1.1 1990/11/30 11:13:00 gdsjaar Exp $
C $Log: nxkord.f,v $
C Revision 1.1  1990/11/30 11:13:00  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]NXKORD.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE NXKORD (NODES, N1)
C***********************************************************************
C
C  SUBROUTINE NXKORD = ROTATES THE LIST OF FOUR NODES SO N1 APPEARS
C                      FIRST IF IT IS IN THE LIST
C
C***********************************************************************
C
      DIMENSION NODES (4)
C
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
C
      RETURN
C
      END
