C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Government retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C    
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C    
C        * Neither the name of Sandia Corporation nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
C    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    

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
