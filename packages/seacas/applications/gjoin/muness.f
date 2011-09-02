C Copyright (c) 2008 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C 
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 

C=======================================================================
      SUBROUTINE MUNESS (NUMESS, ISTAT, LESSEL, LESSDL, 
     &  IDESS, NEESS, NEDSS, IXEESS, IXEDSS,
     &  LTEESS, LTSSS, FACSS, 
     &  LTEX, LTSX, TDX, IXESS, IXDSS, NEX, NDX, ISCR)
C=======================================================================
C $Id: muness.f,v 1.2 2001/06/26 17:38:54 gdsjaar Exp $
C
C   --*** MUNESS *** (GJOIN) Compress and rearrange element side sets
C   --   Written by Amy Gilkey - revised 02/25/88
C   --
C   --MUNESS processes the element side sets according to the set status.
C   --Sets may be combined or deleted.
C   --
C   --Parameters:
C   --   NUMESS - IN/OUT - the number of element side sets
C   --   ISTAT - IN - the status of each set:
C   --      0 = same
C   --      - = delete
C   --      n = combine with set n
C   --   LESSEL - IN/OUT - the length of the element side sets element list
C   --   LESSDL - IN/OUT - the length of the element side sets dist-fact list
C   --   IDESS - IN/OUT - the element side set ID for each set
C   --   NEESS - IN/OUT - the number of elements for each set
C   --   NEDSS - IN/OUT - the number of dist-fact for each set
C   --   IXEESS - IN/OUT - the index of the first element for each set
C   --   IXEDSS - IN/OUT - the index of the first dist-fact for each set
C   --   LTEESS - IN/OUT - the elements for all sets
C   --   LTSSS - IN/OUT - the sides for all sets
C   --   FACSS - IN/OUT - the dist-fact for all sets
C   --   LTEX - SCRATCH - sized to hold the set elements
C   --   LTSX - SCRATCH - sized to hold the set sides
C   --   TDX - SCRATCH - sized to hold the set dist-fact
C   --   IXESS - SCRATCH - size = NUMESS
C   --   IXDSS - SCRATCH - size = NUMESS -- dist-fact
C   --   NEX - SCRATCH - size = NUMESS
C   --   NDX - SCRATCH - size = NUMESS -- dist-face
C   --   ISCR - SCRATCH - size = NUMESS

      INTEGER ISTAT(*)
      INTEGER IDESS(*)
      INTEGER NEESS(*), NEDSS(*)
      INTEGER IXEESS(*)
      INTEGER LTEESS(*), LTEX(*)
      INTEGER LTSSS(*), LTSX(*)
      INTEGER IXESS(*), IXEDSS(*)
      INTEGER NEX(*), NDX(*)
      INTEGER ISCR(*)
      REAL    FACSS(*), TDX(*)

      IF (NUMESS .LE. 0) RETURN

      JESS = 0
      JNE = 1
      JND = 1
      DO 110 IESS = 1, NUMESS

         IF (ISTAT(IESS) .EQ. 0) THEN
            NINSET = 1
            ISCR(NINSET) = IESS
         ELSE IF (ISTAT(IESS) .EQ. IESS) THEN
            CALL GETALL (IESS, NUMESS, ISTAT, NINSET, ISCR)
         ELSE
            NINSET = 0
         END IF

         IF (NINSET .GT. 0) THEN
            JESS = JESS + 1
            IXESS(JESS) = IESS
            NEX(JESS) = 0
            NDX(JESS) = 0
         END IF
         DO 100 ISET = 1, NINSET
            N = ISCR(ISET)
            CALL MOVINT (NEESS(N), LTEESS(IXEESS(N)), LTEX(JNE))
            CALL MOVINT (NEESS(N), LTSSS(IXEESS(N)),  LTSX(JNE))
            CALL MOVREA (NEDSS(N), FACSS(IXEDSS(N)),  TDX(JND))
            JNE = JNE + NEESS(N)
            JND = JND + NEDSS(N)
            NEX(JESS) = NEX(JESS) + NEESS(N)
            NDX(JESS) = NDX(JESS) + NEDSS(N)
  100    CONTINUE
  110 CONTINUE

      CALL ORDIX (JESS, IXESS, NUMESS, IDESS, ISCR, IDESS)
      CALL MOVINT (JESS, NEX, NEESS)
      CALL MOVINT (JESS, NDX, NEDSS)
      NUMESS = JESS
      JNE = 1
      JND = 1
      DO 120 IESS = 1, NUMESS
         IXEESS(IESS) = JNE
         IXEDSS(IESS) = JND
         JNE = JNE + NEESS(IESS)
         JND = JND + NEDSS(IESS)
  120 CONTINUE
      LESSEL = JNE - 1
      LESSDL = JND - 1
      CALL MOVINT (LESSEL, LTEX, LTEESS)
      CALL MOVINT (LESSEL, LTSX, LTSSS)
      CALL MOVREA (LESSDL, TDX,  FACSS)
      RETURN
      END
