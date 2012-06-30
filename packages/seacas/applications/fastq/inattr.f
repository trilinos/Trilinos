C $Id: inattr.f,v 1.1 1990/11/30 11:09:14 gdsjaar Exp $
C $Log: inattr.f,v $
C Revision 1.1  1990/11/30 11:09:14  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]INATTR.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INATTR (MS, MR, MA, N17, N23, JJ, RIN, IFOUND, ATTRIB,
     &   LINKM, NOROOM)
C***********************************************************************
C
C  SUBROUTINE INATTR = INPUTS MATERIAL ATTRIBUTES INTO THE DATABASE
C
C***********************************************************************
C
      DIMENSION ATTRIB (MA, MR+MS), LINKM (2,  (MS+MR))
C
      LOGICAL NOROOM, ADDLNK
C
      NOROOM = .TRUE.
      ADDLNK = .FALSE.
C
C  UPDATE THE COUNTER IF NEEDED
C
      IF (JJ.GT.N23)N23 = JJ
C
C  ADD THE ATTRIBUTES INTO THE DATABASE
C
      N17 = N17 + 1
      J = N17
      IF (J .GT. (MS + MR))RETURN
      CALL LTSORT (MS + MR, LINKM, JJ, IPNTR, ADDLNK)
      ADDLNK = .TRUE.
      IF (IPNTR .LE. 0)THEN
         J = -J
         MINUSJ  =  -J
         CALL LTSORT (MS + MR, LINKM, JJ, MINUSJ, ADDLNK)
         J = IABS (MINUSJ)
      ELSE
         CALL LTSORT (MS + MR, LINKM, JJ, J, ADDLNK)
      ENDIF
      IF (IFOUND .GT. MA)THEN
         IEND = MA
         WRITE (*, 10000)IFOUND, MA
      ELSE
         IEND = IFOUND
      ENDIF
      DO 100 I = 1, IEND
         ATTRIB (J, I) = RIN (I)
  100 CONTINUE
C
      NOROOM = .FALSE.
C
      RETURN
C
10000 FORMAT (' FOR MATERIAL NUMBER:', I5,
     &   ' NUMBER OF ATTRIBUTES READ:',  I5, /,
     &   '                               EXCEEDS MAX ALLOWED OF:',
     &   I5)
C
      END
