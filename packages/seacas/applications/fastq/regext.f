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

C $Id: regext.f,v 1.2 1999/06/17 19:02:22 gdsjaar Exp $
C $Log: regext.f,v $
C Revision 1.2  1999/06/17 19:02:22  gdsjaar
C Fixed several problems related to holes.  In several places, a
C nonpositive integer was being used to index into an array.  This seems
C to fix all of those cases.  I'm not sure if I fixed the true cause of
C these errors or just the symptom though...
C
C Revision 1.1.1.1  1990/11/30 11:14:37  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:14:36  gdsjaar
c Initial revision
c
C
CC* FILE: [.MAIN]REGEXT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE REGEXT (MP, ML, MS, MR, N, II, COOR, ILINE, LTYPE,
     &   LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP,
     &   LINKL, LINKS, LINKR, XMIN, XMAX, YMIN, YMAX)
C***********************************************************************
C
C  SUBROUTINE REGEXT = GETS THE REGION EXTREMES
C
C***********************************************************************
C
      DIMENSION COOR (2, MP)
      DIMENSION ILINE (ML), LTYPE (ML), LCON (3, ML)
      DIMENSION NLPS (MS), IFLINE (MS), ILLIST (MS * 3)
      DIMENSION NSPR (MR), IFSIDE (MR), ISLIST (MR * 4)
      DIMENSION LINKP (2, MP), LINKL (2, ML), LINKS (2, MS)
      DIMENSION LINKR (2, MR)
      DIMENSION N (29)
C
      LOGICAL FOUND, GETMAX, ADDLNK
      LOGICAL NUMPLT, TEST
C
      ADDLNK = .FALSE.
      GETMAX = .TRUE.
      FOUND = .FALSE.
C
      DO 110 J = IFSIDE (II), IFSIDE (II) + NSPR (II) - 1
C
C  GET SIDE EXTREMES
C
         IF ( ISLIST(J) .GT. 0) THEN
            CALL LTSORT (MS, LINKS, ISLIST (J), IPNTR, ADDLNK)
            IF (IPNTR .GT. 0) THEN
               DO 100 K = IFLINE (IPNTR), IFLINE (IPNTR) +
     &              NLPS (IPNTR) - 1
                 CALL LTSORT (ML, LINKL, ILLIST (K), KK, ADDLNK)
                 IF (KK .GT. 0) THEN
                    IF (.NOT.FOUND) THEN
                       CALL LTSORT (MP, LINKP, IABS (LCON (1, KK)),
     &                      IPNT, ADDLNK)
                       IF (IPNT .GT. 0) THEN
                          XMAX = COOR (1, IPNT)
                          XMIN = COOR (1, IPNT)
                          YMAX = COOR (2, IPNT)
                          YMIN = COOR (2, IPNT)
                          FOUND = .TRUE.
                       ENDIF
                    ENDIF
                    IF (FOUND) THEN
                       CALL DLINE (MP, ML, COOR, LINKP, ILINE (KK),
     &                      LTYPE (KK), LCON (1, KK), LCON (2, KK),
     &                      LCON (3, KK), NUMPLT, X1, Y1, TEST, GETMAX,
     &                      XMIN, XMAX, YMIN, YMAX)
                    ENDIF
                 ENDIF
 100           CONTINUE
            END IF
C
C  GET LINE EXTREMES
C
         ELSEIF (ISLIST (J) .LT. 0) THEN
            JJ = IABS (ISLIST (J))
            CALL LTSORT (ML, LINKL, JJ, KK, ADDLNK)
            IF (KK .GT. 0) THEN
               IF (.NOT.FOUND) THEN
                  CALL LTSORT (MP, LINKP, IABS (LCON (1, KK)), IPNT,
     &               ADDLNK)
                  IF (IPNT .GT. 0) THEN
                     XMAX = COOR (1, IPNT)
                     XMIN = COOR (1, IPNT)
                     YMAX = COOR (2, IPNT)
                     YMIN = COOR (2, IPNT)
                     FOUND = .TRUE.
                  ENDIF
               ENDIF
               IF (FOUND) THEN
                  CALL DLINE (MP, ML, COOR, LINKP, ILINE (KK),
     &               LTYPE (KK), LCON (1, KK), LCON (2, KK),
     &               LCON (3, KK), NUMPLT, X1, Y1, TEST, GETMAX, XMIN,
     &               XMAX, YMIN, YMAX)
               ENDIF
            ENDIF
         ENDIF
  110 CONTINUE
C
      RETURN
C
      END
