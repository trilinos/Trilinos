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
      SUBROUTINE ZMESS (NUMESS, LESSEL, LESSDL, 
     &   IDESS, NEESS, NEDSS, IXEESS, IXEDSS, LTEESS, LTSSS, LTSNC, FAC)
C=======================================================================
C $Id: zmess.f,v 1.3 2002/01/28 19:44:47 gdsjaar Exp $

C   --*** ZMESS *** (GJOIN) Compress element side sets
C   --   Written by Amy Gilkey - revised 01/20/88
C   --
C   --ZMESS compresses the element side sets by removing deleted elements
C   --and their nodes.  Assumes that the elements and nodes are already
C   --renumbered and ordered.
C   --
C   --Parameters:
C   --   NUMESS - IN/OUT - the number of element side sets
C   --   LESSEL - IN/OUT - the length of the element side sets element list
C   --   IDESS - IN/OUT - the element side set ID for each set
C   --   NEESS - IN/OUT - the number of elements for each set
C   --   NEDSS - IN/OUT - the number of dist-fac for each set
C   --   IXEESS - IN/OUT - the index of the first element for each set
C   --   IXEDSS - IN/OUT - the index of the first dist-fac for each set
C   --   LTEESS - IN/OUT - the elements for all sets
C   --   LTSSS - IN/OUT - the sides for all sets
C   --   LTSNC - IN/OUT - the face count for each element/side in the list
C   --   FACESS - IN/OUT - the distribution factors for all sets????????????

      INTEGER IDESS(*)   ! NUMESS
      INTEGER NEESS(*)   ! NUMESS
      INTEGER NEDSS(*)   ! NUMESS
      INTEGER IXEESS(*)  ! NUMESS
      INTEGER IXEDSS(*)  ! NUMESS
      INTEGER LTEESS(*)  ! LESSEL
      INTEGER LTSSS(*)   ! LESSEL
      INTEGER LTSNC(*)   ! LESSEL
      REAL    FAC(*)     ! LESSDL

      IF (NUMESS .LE. 0) RETURN

      JESS = 0
      JNE = 0
      JDF = 0
      DO 120 IESS = 1, NUMESS
         JNELST = JNE
         JDFLST = JDF
         nd1 = ixedss(iess) ! Index of First distribution factor for this list
         IDFB = ND1
         DO 110 N = 1, NEESS(IESS)
C ... N     is the 'local' index within the current set.
C     NE    is the 'global' index within the concatenated (LTEESS, LTSSS, LTSNC) lists
C     ND1   is the FAC index of the first df for the current list.
C     ICNT  is the number of DF for element N in the current list
           NE = N + IXEESS(IESS) - 1 
           ICNT = LTSNC(NE)
C IDFB = index of first df for local element N, global element NE
C IDFE = index of last  df for local element N, global element NE
           IDFE = IDFB + ICNT - 1
           IF (LTEESS(NE) .GT. 0) THEN
               JNE = JNE + 1
               LTEESS(JNE) = LTEESS(NE)
               LTSSS(JNE)  = LTSSS(NE)
               LTSNC(JNE)  = LTSNC(NE)
               do 100 nd = idfb, idfe
                 JDF = JDF + 1
                 fac(JDF) = fac(ND)
 100           continue
            END IF
            IDFB = IDFE + 1
  110    CONTINUE
         N = JNE - JNELST
         IF (N .GT. 0) THEN
C ... There is at least 1 element remaining in the list...
            JESS = JESS + 1  ! increment sideset count
            IDESS(JESS) = IDESS(IESS) ! copy the sideset id
            NEESS(JESS) = N ! Set the elements per list count
            IXEESS(JESS) = JNELST + 1 ! set the index 
            NEDSS(JESS) = JDF - JDFLST ! set the df per list count
            IXEDSS(JESS) = JDFLST + 1
         END IF
  120 CONTINUE
      if (idfe .ne. lessdl) stop 'ZMESS: Internal error'
      NUMESS = JESS
      LESSEL = JNE
      LESSDL = JDF
      RETURN
      END
