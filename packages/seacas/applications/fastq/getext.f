C $Id: getext.f,v 1.2 1999/06/17 19:02:22 gdsjaar Exp $
C $Log: getext.f,v $
C Revision 1.2  1999/06/17 19:02:22  gdsjaar
C Fixed several problems related to holes.  In several places, a
C nonpositive integer was being used to index into an array.  This seems
C to fix all of those cases.  I'm not sure if I fixed the true cause of
C these errors or just the symptom though...
C
C Revision 1.1.1.1  1990/11/30 11:08:04  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:08:02  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]GETEXT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETEXT (MP, ML, MS, MR, N, IPOINT, COOR, ILINE, LTYPE,
     &   LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP,
     &   LINKL, LINKS, LINKR, REXTRM, BXMIN, BXMAX, BYMIN, BYMAX)
C***********************************************************************
C
C  SUBROUTINE GETEXT = GETS THE REGION AND BODY EXTREMES
C
C***********************************************************************
C
      DIMENSION IPOINT (MP), COOR (2, MP)
      DIMENSION ILINE (ML), LTYPE (ML), LCON (3, ML)
      DIMENSION NLPS (MS), IFLINE (MS), ILLIST (MS * 3)
      DIMENSION NSPR (MR), IFSIDE (MR), ISLIST (MR * 4)
      DIMENSION LINKP (2, MP), LINKL (2, ML), LINKS (2, MS)
      DIMENSION LINKR (2, MR)
      DIMENSION REXTRM (4, MR), N (29)
C
      LOGICAL FOUND, GETMAX, ADDLNK
      LOGICAL NUMPLT, TEST
C
C  GET THE POINTS EXTREMES
C
      ADDLNK = .FALSE.
      GETMAX = .TRUE.
      FOUND = .FALSE.
      DO 100 I = 1, N (1)
         CALL LTSORT (MP, LINKP, IABS (IPOINT (I)), II, ADDLNK)
         IF (II .GT. 0) THEN
            IF (FOUND) THEN
               BXMAX = AMAX1 (COOR (1, II), BXMAX)
               BYMAX = AMAX1 (COOR (2, II), BYMAX)
               BXMIN = AMIN1 (COOR (1, II), BXMIN)
               BYMIN = AMIN1 (COOR (2, II), BYMIN)
            ELSE
               BXMAX = COOR (1, II)
               BXMIN = COOR (1, II)
               BYMAX = COOR (2, II)
               BYMIN = COOR (2, II)
               FOUND = .TRUE.
            ENDIF
         ENDIF
  100 CONTINUE
C
C  GET ALL THE LINES EXTREMES
C
      IF (FOUND) THEN
         DO 110 I = 1, N (2)
            CALL LTSORT (ML, LINKL, IABS (ILINE (I)), II, ADDLNK)
            IF (II .GT. 0) THEN
               CALL DLINE (MP, ML, COOR, LINKP, ILINE (II),
     &            LTYPE (II), LCON (1, II), LCON (2, II), LCON (3, II),
     &            NUMPLT, X1, Y1, TEST, GETMAX, BXMIN, BXMAX, BYMIN,
     &            BYMAX)
            ENDIF
  110    CONTINUE
      ELSE
         BXMIN = 0.
         BXMAX = 20.
         BYMIN = 0.
         BYMAX = 15.
      ENDIF
C
C  CALCULATE THE EXTREMES FOR EACH REGION
C
      DO 120 I = 1, N (22)
         CALL LTSORT (MR, LINKR, I, II, ADDLNK)
         IF (II .GT. 0) THEN
            CALL REGEXT (MP, ML, MS, MR, N, II, COOR, ILINE, LTYPE,
     &         LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP,
     &         LINKL, LINKS, LINKR, XMIN, XMAX, YMIN, YMAX)
            REXTRM (1, II) = XMIN
            REXTRM (2, II) = XMAX
            REXTRM (3, II) = YMIN
            REXTRM (4, II) = YMAX
         ENDIF
  120 CONTINUE
C
      RETURN
C
      END
