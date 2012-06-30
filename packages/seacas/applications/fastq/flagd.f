C $Id: flagd.f,v 1.1 1990/11/30 11:07:33 gdsjaar Exp $
C $Log: flagd.f,v $
C Revision 1.1  1990/11/30 11:07:33  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]FLAGD.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE FLAGD (MDIM, N, LINK, INUM, FLAG)
C***********************************************************************
C
C  SUBROUTINE FLAGD = FLAGS THE DATA TO BE PLOTTED
C
C***********************************************************************
C
      DIMENSION LINK(2,MDIM), INUM(MDIM)
C
      LOGICAL FLAG, ADDLNK
C
      ADDLNK = .FALSE.
C
      DO 100 I = 1, N
         CALL LTSORT (MDIM, LINK, I, II, ADDLNK)
         IF (II .GT. 0) THEN
            IF (FLAG) THEN
               INUM(II) = -IABS (INUM (II))
            ELSE
               INUM(II) = IABS (INUM (II))
            ENDIF
         ENDIF
  100 CONTINUE
      RETURN
      END
