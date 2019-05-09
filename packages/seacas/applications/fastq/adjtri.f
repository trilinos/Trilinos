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

C $Id: adjtri.f,v 1.2 1998/07/14 18:18:18 gdsjaar Exp $
C $Log: adjtri.f,v $
C Revision 1.2  1998/07/14 18:18:18  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:03:27  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:03:26  gdsjaar
c Initial revision
c
C
CC* FILE: [.PAVING]ADJTRI.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE ADJTRI (MXND, MLN, LNODES, XN, YN, ZN, NUID, LXK, KXL,
     &   NXL, LXN, NNN, NAVAIL, IAVAIL, NODE, KELEM, ANG, TOLER1,
     &   TOLER2, N1, N2, N3, KREG, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     &   KKK, LLL, DEV1, DONE, CHECK, GRAPH, VIDEO, NOROOM, ERR, KKKADD)
C***********************************************************************
C
C  SUBROUTINE ADJTRI = ADJUSTS A TRIANGULAR SHAPED ELEMENT WHERE
C                      POSSIBLE
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY TRIDEL
C
C***********************************************************************
C
C  THERE ARE THREE POSSIBILITIES FOR CHANGE:
C     1) ANYTHING OVER TOLER1 GETS THE CORRESPONDING ELEMENT
C        DELETED
C     2) ANYTHING OVER TOLER2 AND HOOKED TO ANOTHER 3-LINE NODE GETS
C        THE CORRESPONDING ELEMENT DELETED
C     3) AN ELONGATED ELEMENT OVER 150 DEGREES GETS A 3 ELEMENT
C        REPLACEMENT FOR THE TWO ELEMENTS THERE
C
C***********************************************************************
C
      DIMENSION LXK(4, MXND), NXL(2, 3*MXND), KXL(2, 3*MXND)
      DIMENSION LXN(4, MXND), XN(MXND), YN(MXND), ZN(MXND), NUID(MXND)
      DIMENSION LNODES (MLN, MXND)
C
      CHARACTER*3 DEV1
C
      LOGICAL NOROOM, ERR, DONE, GRAPH, CHECK, VIDEO
C
C  CHECK FOR CASE 1
C
      IF (ANG .GT. TOLER1) THEN
         IF (GRAPH) THEN
            CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN,
     &         YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
            CALL LCOLOR ('PINK ')
            CALL D2NODE (MXND, XN, YN, NODE, N1)
            CALL D2NODE (MXND, XN, YN, NODE, N2)
            CALL LCOLOR ('WHITE')
            CALL SFLUSH
         ENDIF
  100    CONTINUE
         CALL DELEM (MXND, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &      NNN, NAVAIL, IAVAIL, NODE, KELEM, NODE1, NODE3, DONE,
     &      CHECK, NOROOM, ERR)
         IF ((NOROOM) .OR. (ERR)) GOTO 120
         IF (LXN (3, NODE1) .LE. 0) THEN
            NODE = NODE1
            KELEM = KXL (1, LXN (1, NODE))
            CHECK = .FALSE.
            GOTO 100
         ELSEIF (LXN (3, NODE3) .LE. 0) THEN
            NODE = NODE3
            KELEM = KXL (1, LXN (1, NODE))
            CHECK = .FALSE.
            GOTO 100
         ENDIF
         CHECK = .TRUE.
         IF ((ERR) .OR. (DONE)) GOTO 120
      ENDIF
C
C  CHECK FOR CASE 2
C
      IF ( (ANG .GT. TOLER2) .AND. (LXN (4, N3) .EQ. 0) .AND.
     &   (LXN (3, N3) .GT. 0) .AND. (LXN (2, N3) .GT. 0)) THEN
         IF (GRAPH) THEN
            CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN,
     &         YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
            CALL LCOLOR ('PINK ')
            CALL D2NODE (MXND, XN, YN, NODE, N1)
            CALL D2NODE (MXND, XN, YN, NODE, N2)
            CALL LCOLOR ('WHITE')
            CALL SFLUSH
         ENDIF
  110    CONTINUE
         CALL DELEM (MXND, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &      NNN, NAVAIL, IAVAIL, NODE, KELEM, NODE1, NODE3, DONE,
     &      CHECK, NOROOM, ERR)
         IF ((NOROOM) .OR. (ERR)) GOTO 120
         IF (LXN (3, NODE1) .LE. 0) THEN
            NODE = NODE1
            KELEM = KXL (1, LXN (1, NODE))
            CHECK = .FALSE.
            GOTO 110
         ELSEIF (LXN (3, NODE3) .LE. 0) THEN
            NODE = NODE3
            KELEM = KXL (1, LXN (1, NODE))
            CHECK = .FALSE.
            GOTO 110
         ENDIF
         CHECK = .TRUE.
         IF ((ERR) .OR. (DONE)) GOTO 120
      ENDIF
C
C  CHECK FOR CASE 3
C
      CALL LONGEL (MXND, MLN, LNODES, XN, YN, NUID, LXK, KXL, NXL,
     &   LXN, NNN, NAVAIL, IAVAIL, NODE, KELEM, ANG, TOLER2,
     &   N1, N2, KREG, XMIN, XMAX, YMIN, YMAX, KKK, LLL, DONE, GRAPH,
     &   VIDEO, NOROOM, ERR, KKKADD)
C
  120 CONTINUE
      RETURN
C
      END
