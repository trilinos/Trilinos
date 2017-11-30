C    Copyright(C) 2014-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C    
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C    
C    * Neither the name of NTESS nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
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

C $Id: movep.F,v 1.1 1998/04/10 04:30:18 gdsjaar Exp $
C $Log: movep.F,v $
C Revision 1.1  1998/04/10 04:30:18  gdsjaar
C Added the ifdef needsdoubleescape to the graphics symbol
C plotting. Required moving some .f to .F files and specifying in Imakefile
C
C Revision 1.1.1.1  1990/11/30 11:12:29  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:12:28  gdsjaar
c Initial revision
c 
C

CC* FILE: [.MAIN]MOVEP.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE MOVEP (MP, ML, MS, MR, N, COOR, LTYPE, LCON, NLPS,
     &   IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP, LINKL, LINKS,
     &   LINKR, IMOVE, ALPHA, X, Y)
C***********************************************************************
C
C  SUBROUTINE MOVEP = MOVES A POINT, AND REDRAWS ALL LINES DEPENDANT ON
C                     THAT POINT
C
C***********************************************************************
C
      DIMENSION COOR (2, MP), LTYPE (ML), LCON (3, ML)
      DIMENSION NLPS (MS), IFLINE (MS)
      DIMENSION ILLIST (MS * 3)
      DIMENSION NSPR (MR), IFSIDE (MR), ISLIST (MR * 4)
      DIMENSION LINKP (2, MP), LINKL (2, ML), LINKS (2, MS)
      DIMENSION LINKR (2, MR)
      DIMENSION N (29)
C
      LOGICAL NUMPLT, ALPHA, ADDLNK, TEST, GETMAX
C
      KNUM = 0
      ADDLNK = .FALSE.
      TEST = .FALSE.
      GETMAX = .FALSE.
      IF (.NOT.ALPHA)CALL PLTSTD (1, 0.)
C
C  ERASE ANY LINES USING THAT POINT
C
      DO 120 I = 1, N (19)
         CALL LTSORT (ML, LINKL, I, II, ADDLNK)
         IF (II .GT. 0) THEN
            DO 100 J = 1, 3
               IF (IMOVE .EQ. LCON (J, II)) THEN
                  LT = LTYPE (II)
                  IP1 = LCON (1, II)
                  IP2 = LCON (2, II)
                  IP3 = LCON (3, II)
                  CALL LTSORT (MP, LINKP, IP1, IPNTR1, ADDLNK)
                  CALL LTSORT (MP, LINKP, IP2, IPNTR2, ADDLNK)
                  IF (IP3 .NE. 0) THEN
                     CALL LTSORT (MP, LINKP, IABS (IP3), IPNTR3, ADDLNK)
                  ELSE
                     IPNTR3 = 0
                  ENDIF
                  IF ( (IPNTR1 .GT. 0) .AND. (IPNTR2 .GT. 0) .AND.
     &               ( (LT .EQ. 1) .OR. (IPNTR3 .GT. 0))) THEN
                     IF (.NOT.ALPHA) THEN
                        CALL DLINE (MP, ML, COOR, LINKP, KNUM, LT, IP1,
     &                     IP2, IP3, NUMPLT, DUM1, DUM2, TEST, GETMAX,
     &                     DUM3, DUM4, DUM5, DUM6)
                        CALL PLTFLU
                     ENDIF
                  ENDIF
                  GOTO 110
               ENDIF
  100       CONTINUE
  110       CONTINUE
         ENDIF
  120 CONTINUE
C
C  ERASE THE POINT
C
      ADDLNK = .FALSE.
      CALL LTSORT (MP, LINKP, IMOVE, II, ADDLNK)
      IF (.NOT.ALPHA) THEN
#if NeedsDoubleEscape
         CALL MPD2SY (1, COOR (1, II), COOR (2, II), '\\CX')
         CALL MPD2SY (1, COOR (1, II), COOR (2, II), '\\CSQ')
#else
         CALL MPD2SY (1, COOR (1, II), COOR (2, II), '\CX')
         CALL MPD2SY (1, COOR (1, II), COOR (2, II), '\CSQ')
#endif
      ENDIF
C
C  REDEFINE THE POINT
C
      COOR (1, II) = X
      COOR (2, II) = Y
C
      IF (.NOT.ALPHA)CALL PLTSTD (1, 7.)
C
C  REDRAW ANY LINES USING THAT POINT
C
      DO 150 I = 1, N (19)
         CALL LTSORT (ML, LINKL, I, II, ADDLNK)
         IF (II .GT. 0) THEN
            DO 130 J = 1, 3
               IF (IMOVE .EQ. LCON (J, II)) THEN
                  LT = LTYPE (II)
                  IP1 = LCON (1, II)
                  IP2 = LCON (2, II)
                  IP3 = LCON (3, II)
                  CALL LTSORT (MP, LINKP, IP1, IPNTR1, ADDLNK)
                  CALL LTSORT (MP, LINKP, IP2, IPNTR2, ADDLNK)
                  IF (IP3 .NE. 0) THEN
                     CALL LTSORT (MP, LINKP, IABS (IP3), IPNTR3, ADDLNK)
                  ELSE
                     IPNTR3 = 0
                  ENDIF
                  IF ( (IPNTR1 .GT. 0) .AND. (IPNTR2 .GT. 0) .AND.
     &               ( (LT .EQ. 1) .OR. (IPNTR3 .GT. 0))) THEN
                     IF (.NOT.ALPHA) THEN
                        CALL DLINE (MP, ML, COOR, LINKP, KNUM, LT, IP1,
     &                     IP2, IP3, NUMPLT, DUM1, DUM2, TEST, GETMAX,
     &                     DUM3, DUM4, DUM5, DUM6)
                        CALL PLTFLU
                     ENDIF
                  ENDIF
                  GOTO 140
               ENDIF
  130       CONTINUE
  140       CONTINUE
         ENDIF
  150 CONTINUE
C
C  REDRAW THE POINT
C
      IF (.NOT.ALPHA)CALL PLTSTD (1, 3.)
C
      ADDLNK = .FALSE.
      CALL LTSORT (MP, LINKP, IMOVE, II, ADDLNK)
      IF (.NOT.ALPHA) THEN
#if NeedsDoubleEscape
         CALL MPD2SY (1, COOR (1, II), COOR (2, II), '\\CX')
#else
         CALL MPD2SY (1, COOR (1, II), COOR (2, II), '\CX')
#endif
         CALL PLTSTD (1, 7.)
         CALL PLTFLU
      ENDIF
C
C
      RETURN
      END
