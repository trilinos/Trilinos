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

C $Id: rdmesh.f,v 1.1 1990/11/30 11:14:30 gdsjaar Exp $
C $Log: rdmesh.f,v $
C Revision 1.1  1990/11/30 11:14:30  gdsjaar
C Initial revision
C
C
CC* FILE: [.RENUM]RDMESH.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE RDMESH (NPNODE, NPELEM, NPNBC, NPSBC, NPREGN, MS, MR,
     &   NNUID, NNXK, IUNIT, NNN, KKK, IPART, LSTNBC, LSTSBC, NUID, XN,
     &   YN, NXK, MAT, MATMAP, NUMMAT, ISIDE, NLPS, IFLINE, ILLIST,
     &   LINKS, LINKR, IMAT, LINKB, JMAT, NNNBC, NNSBC, ERR)
C***********************************************************************
C
C  SUBROUTINE RDMESH = THE CURRENT MESH STORED ON DISK
C
C***********************************************************************
C
      DIMENSION IPART(3, NPREGN), ISIDE(MS), NLPS(MS), IFLINE(MS)
      DIMENSION ILLIST(MS*3), LINKS(2, MS), LINKB(2, MS), LINKR(2, MR)
      DIMENSION IMAT(MR), JMAT(MS)
C
      DIMENSION NUID(NNUID), XN(NPNODE), YN(NPNODE)
      DIMENSION NXK(NNXK, NPELEM), MAT(NPELEM)
C
      DIMENSION LSTNBC(NPNBC), LSTSBC(NPSBC)
      DIMENSION MATMAP(3, NPREGN)
C
      LOGICAL ERR, BAR, ADDLNK
C
      ERR = .FALSE.
C
C  READ THE MESH TAPE
C
      REWIND IUNIT
      NNN = 0
      KKK = 0
      NNNBC = 0
      NNSBC = 0
      ADDLNK = .FALSE.
      DO 100 I = 1, NPREGN
         IPART(1, I) = 0
         IPART(2, I) = 0
         IPART(3, I) = 0
  100 CONTINUE
C
      NUMMAT = 0
      DO 220 IR = 1, NPREGN
         READ(IUNIT, END = 230) KKKREG, NNNREG, NNBCRG, NSBCRG, KREG,
     &      BAR, M1
C
C  READ THE NODES
C
         N1 = NNN + 1
         NNN = NNN + NNNREG
         IF (NNNREG .GE. 1)READ(IUNIT, END = 230) (NUID(I), XN(I),
     &      YN(I), I = N1, NNN)
C
C  READ THE ELEMENTS
C
         K1 = KKK + 1
         IPART(1, IR) = KREG
         IPART(2, IR) = K1
         KKK = KKK + KKKREG
         IPART(3, IR) = KKK
         READ(IUNIT, END = 230) ((NXK(I, K), I = 1, 4), K = K1, KKK)
C
C  ZERO THE MIDSIDE NODE LOCATIONS IN THE NXK ARRAY
C
         DO 120 I = 5, NNXK
            DO 110 K = K1, KKK
               NXK(I, K) = 0
  110       CONTINUE
  120    CONTINUE
C
C  SET UP THE MATERIAL ARRAY AND MAXIMUM NUMBER OF MATERIALS
C
         IF (BAR) THEN
            CALL LTSORT (MS, LINKB, KREG, IPNTR, ADDLNK)
            KMAT = ABS (JMAT(IPNTR))
         ELSE
            CALL LTSORT (MR, LINKR, KREG, IPNTR, ADDLNK)
            KMAT = IMAT(IPNTR)
         END IF
C
C  SEE IF ALTERNATING MATERIALS WITHIN A REGION ARE ENABLED
C
         IF (KMAT .LT. 0) THEN
            CALL LTSORT (MS, LINKS, IABS(KMAT), JPNTR, ADDLNK)
            IF ((JPNTR .GT. 0) .AND. (NLPS(JPNTR) .GE. 2)) THEN
C
C  ADD MATERIAL NUMBER BY ROW OF ELEMENTS
C
               MATPNT = IFLINE(JPNTR)
               DO 140 K = K1, KKK - M1 + 1, M1
                  DO 130 L = K, K + M1 - 1
                     MAT(L) = ILLIST(MATPNT)
  130             CONTINUE
C
C  UPDATE THE POINTER TO THE NEXT MATERIAL
C
                  MATPNT = MATPNT + 1
                  IF (MATPNT .GT. IFLINE(JPNTR) + NLPS(JPNTR) - 1)
     &               MATPNT = IFLINE(JPNTR)
  140          CONTINUE
            ELSE
               WRITE(*, 10000) IABS(KMAT), ISIDE(JPNTR)
            END IF
C
C  PUT THE NEW MATERIALS INTO THE MATERIAL ARRAYS
C
            DO 170 MATPNT = IFLINE(JPNTR),
     &         IFLINE(JPNTR) + NLPS(JPNTR) - 1
               DO 150 K = 1, NUMMAT
                  IF (MATMAP(1, K) .EQ. ILLIST(MATPNT)) GO TO 160
  150          CONTINUE
               NUMMAT = NUMMAT + 1
               MATMAP(1, NUMMAT) = ILLIST(MATPNT)
  160          CONTINUE
  170       CONTINUE
         ELSE
C
C  JUST INPUT THE ONE MATERIAL
C
            DO 180 K = K1, KKK
               MAT(K) = KMAT
  180       CONTINUE
            DO 190 K = 1, NUMMAT
               IF (MATMAP(1, K) .EQ. KMAT) GO TO 200
  190       CONTINUE
            NUMMAT = NUMMAT + 1
            MATMAP(1, NUMMAT) = KMAT
         END IF
  200    CONTINUE
C
C  READ THE NODAL BOUNDARY CONDITIONS
C
         NNNBC1 = NNNBC + 1
         NNNBC = NNNBC + NNBCRG
         IF (NNBCRG .GE. 1) READ(IUNIT)(LSTNBC(I), I = NNNBC1, NNNBC)
C
C  READ THE SIDE BOUNDARY CONDITIONS
C
         NNSBC1 = NNSBC + 1
         NNSBC = NNSBC + NSBCRG
         IF (NSBCRG .GE. 1) THEN
            READ(IUNIT)(LSTSBC(I), I = NNSBC1, NNSBC)
            DO 210 I = NNSBC1 + 1, NNSBC, 3
               LSTSBC(I) = LSTSBC(I) + K1 - 1
  210       CONTINUE
         END IF
  220 CONTINUE
C
      RETURN
C
  230 CONTINUE
      CALL MESAGE ('PREMATURE END OF FILE ON MESH READ')
      CALL MESAGE ('CHECK MESH PROCESSING OUTPUT TO DETERMINE')
      CALL MESAGE ('INCOMPLETE MESH DESCRIPTION')
      ERR = .TRUE.
      NNN = 0
      KKK = 0
      RETURN
C
10000 FORMAT(' THE ALTERNATING MATERIAL NUMBERS FOR REGION(S) WITH', /,
     &   ' NEGATIVE MATERIAL NUMBER:', I5, ' DO NOT CORRESPOND TO A',/,
     &   ' VALID SIDE NUMBER:', I5, ' WITH AT LEAST TWO LINES.')
C
      END
