C $Id: mkused.f,v 1.1 1990/11/30 11:12:17 gdsjaar Exp $
C $Log: mkused.f,v $
C Revision 1.1  1990/11/30 11:12:17  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]MKUSED.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE MKUSED (MXNL, MP, ML, LISTL, IPOINT, NINT, LINKP,
     &   LINKL, LCON, NL)
C***********************************************************************
C
C  SUBROUTINE MKUSED = MARKS ALL LINES AND POINTS USED IN THE PERIMETER
C
C***********************************************************************
C
C  NOTE:
C     THE MESH TABLES ARE EFFECTIVELY DESTROYED BY THIS ROUTINE
C
C***********************************************************************
C
      DIMENSION LISTL (MXNL), IPOINT (MP), NINT (ML)
      DIMENSION LINKP (2, MP), LINKL (2, ML)
      DIMENSION LCON (3, ML)
C
      LOGICAL ADDLNK
C
      ADDLNK = .FALSE.
C
C  FLAG ALL LINES AND POINTS AS BEING USED
C
      DO 100 I = 1, NL
         CALL LTSORT (ML, LINKL, LISTL (I), IL, ADDLNK)
         NINT (IL) =  - IABS (NINT (IL))
         CALL LTSORT (MP, LINKP, LCON (1, IL), IP1, ADDLNK)
         CALL LTSORT (MP, LINKP, LCON (2, IL), IP2, ADDLNK)
         IPOINT (IP1) =  - IABS (IPOINT (IP1))
         IPOINT (IP2) =  - IABS (IPOINT (IP2))
  100 CONTINUE
C
      RETURN
C
      END
