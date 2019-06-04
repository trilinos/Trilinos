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

C $Id: inbrst.f,v 1.2 1999/06/21 22:43:40 gdsjaar Exp $
C $Log: inbrst.f,v $
C Revision 1.2  1999/06/21 22:43:40  gdsjaar
C Fixed more uninitialized variables; one was causing core dump on g77
C compiled executable.
C
C VERSN was not consistently defined -- now 10 characters everywhere
C
C Updated so full version string output
C
C Added capability to debug memory using unit specified in EXT99
C variable. Similar to STRTUP in SUPLIB
C
C Cleaned up some other code
C
C Upped version
C
C Revision 1.1.1.1  1990/11/30 11:09:24  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:09:22  gdsjaar
c Initial revision
c
C
CC* FILE: [.MAIN]INBRST.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INBRST (MS, MR, N5, N6, N21, N23, JJ, IMTRL, JC, IIN,
     &   IFOUND, IBARST, JMAT, JCENT, NLPB, JFLINE, JLLIST, LINKB,
     &   LINKM, NHOLDM, IHOLDM, NHOLDB, IHOLDB, MERGE, NOROOM)
C***********************************************************************
C
C  SUBROUTINE INBRST = INPUTS A BAR SET INTO THE DATABASE
C
C***********************************************************************
C
      DIMENSION IBARST (MS), JMAT (MS), JCENT (MS), NLPB (MS)
      DIMENSION JFLINE (MS)
      DIMENSION JLLIST (3 * MS), LINKB (2, MS), LINKM (2, MS + MR)
      DIMENSION IHOLDM (2,  (MS + MR)), IHOLDB (2, MS)
      DIMENSION IIN (IFOUND)
C
      LOGICAL MERGE, NOROOM, ADDLNK
C
      IZ = 0
      NOROOM = .TRUE.
      N22 = 0
C
C  ZERO OUT THE LINK ARRAY IF NEEDED
C
      IF (JJ .GT. N21) THEN
         N21 = JJ
C
C  FIND THE CORRECT BAR SET NUMBER IF MERGING
C
C
C  SET UP POINTERS FOR MERGING DATA
C
      ELSEIF (MERGE) THEN
         JHOLD = JJ
         CALL LTSORT (MS, LINKB, JJ, IPNTR, ADDLNK)
         IF (IPNTR .GT. 0) THEN
            IF (JHOLD .GT. NHOLDB)NHOLDB = JHOLD
            CALL LTSORT (MS, IHOLDB, JHOLD, IPNTR, ADDLNK)
            IF (IPNTR .GT. 0) THEN
               JJ = IPNTR
            ELSE
               JJ = N22 + 1
               N22 = JJ
               WRITE ( * , 10000)JHOLD, JJ
               ADDLNK = .TRUE.
               CALL LTSORT (MS, IHOLDB, JHOLD, JJ, ADDLNK)
            ENDIF
         ENDIF
      ENDIF
C
C  INPUT THE BAR SET DATA INTO THE DATABASE
C
      N5 = N5 + 1
      J = N5
      IF (J .GT. MS)RETURN
      ADDLNK = .TRUE.
      CALL LTSORT (MS, LINKB, JJ, J, ADDLNK)
      IBARST (J) = JJ
      JCENT (J) = JC
      JFLINE (J) = N6 + 1
      DO 100 I = 1, IFOUND
         JJ = IIN (I)
         IF (JJ .EQ. 0)GOTO 110
         N6 = N6 + 1
         IF (N6 .GT. MS * 3)RETURN
         JLLIST (N6) = JJ
  100 CONTINUE
  110 CONTINUE
      NLPB (J) = N6 - JFLINE (J) + 1
      IF (NLPB (J) .LT. 1) THEN
         WRITE ( * , 10010)J
         CALL LTSORT (MS, LINKB, JJ, IZ, ADDLNK)
      ENDIF
      ADDLNK = .FALSE.
C
C
C  LINK UP THE MATERIAL
C
C
C  ZERO THE LINK ARRAY IF NEEDED
C
      IF (IMTRL .GT. N23) THEN
         N23 = IMTRL
C
C  SET UP POINTERS FOR MERGING DATA
C
      ELSEIF (MERGE) THEN
         JHOLD = IMTRL
         CALL LTSORT (MS + MR, LINKM, IMTRL, IPNTR, ADDLNK)
         IF (IPNTR .NE. 0) THEN
            IF (JHOLD .GT. NHOLDM)NHOLDM = JHOLD
            ADDLNK = .FALSE.
            CALL LTSORT ( (MS + MR), IHOLDM, JHOLD, IPNTR, ADDLNK)
            IF (IPNTR .GT. 0) THEN
               IMTRL = IPNTR
            ELSE
               IMTRL = N23 + 1
               N23 = IMTRL
               WRITE ( * , 10010)JHOLD, IMTRL
               ADDLNK = .TRUE.
               CALL LTSORT ( (MS + MR), IHOLDM, JHOLD, IMTRL, ADDLNK)
            ENDIF
         ENDIF
      ENDIF
C
C  ADD THE MATERIAL INTO THE DATABASE
C
      NOROOM = .FALSE.
      ADDLNK = .FALSE.
      CALL LTSORT (MS + MR, LINKM, IMTRL, IPNTR, ADDLNK)
      ADDLNK = .TRUE.
      IF (IPNTR .GT. 0) THEN
         CALL MESAGE (' ')
         WRITE ( * , 10020)IMTRL, IBARST (J)
         CALL LTSORT (MS, LINKB, IBARST (J), IZ, ADDLNK)
         RETURN
      ELSEIF (IPNTR .EQ. 0) THEN
         IMINUS =  - 1
         CALL LTSORT (MS + MR, LINKM, IMTRL, IMINUS, ADDLNK)
      ENDIF
      JMAT (J) = IMTRL
C
      RETURN
C
10000 FORMAT ('   OLD BAR SET NO:', I5, '   TO NEW BAR SET NO:', I5)
10010 FORMAT (' BAR SET:', I5, ' HAS LESS THAN ONE LINE',  / ,
     &   ' THIS BAR SET WILL NOT BE INPUT INTO DATABASE')
10020 FORMAT (' MATERIAL:', I5, ' FOR BAR SET:', I5,
     &   ' HAS BEEN DESIGNATED',
     &   / , ' AS A REGION  (4 NODE ELEMENT) MATERIAL.',  / ,
     &   ' ELEMENTS WITH 2 AND 4 NODES CANNOT SHARE MATERIAL ID''S',
     &   / , ' THIS BAR SET WILL NOT BE INPUT INTO DATABASE')
      END
