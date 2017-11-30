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

C $Id: inline.f,v 1.2 1998/12/08 14:26:04 gdsjaar Exp $
C $Log: inline.f,v $
C Revision 1.2  1998/12/08 14:26:04  gdsjaar
C Detect whether negative line intervals entered. Output warning message
C and fix (make positive).
C
C Upped version to 2.10
C
C Revision 1.1.1.1  1990/11/30 11:09:53  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:09:51  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]INLINE.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INLINE (ML, N2, N19, JJ, LTYP, IP1, IP2, IP3, NN,
     &   FACT, NHOLDL, IHOLDL, ILINE, LTYPE, NINT, FACTOR, LCON,
     &   ILBOUN, ISBOUN, LINKL, MERGE, NOROOM)
C***********************************************************************
C
C  SUBROUTINE INLINE = INPUTS A LINE INTO THE DATABASE
C
C***********************************************************************
C
      DIMENSION ILINE (ML), LTYPE (ML), NINT (ML), FACTOR (ML)
      DIMENSION LCON (3, ML)
      DIMENSION ILBOUN (ML), ISBOUN (ML), LINKL (2, ML), IHOLDL (2, ML)
C
      LOGICAL MERGE, NOROOM, ADDLNK
C
      NOROOM = .TRUE.
      JHOLD = JJ
C
C  ADJUST THE COUNTER IF NEEDED
C
      IF (JJ .GT. N19) THEN
         N19 = JJ
C
C  GET THE CORRECT LINE NUMBER IF MERGING
C
      ELSEIF (MERGE) THEN
         ADDLNK = .FALSE.
         CALL LTSORT (ML, LINKL, JJ, IPNTR, ADDLNK)
         IF (IPNTR .GT. 0) THEN
            IF (JHOLD .GT. NHOLDL)NHOLDL = JHOLD
            CALL LTSORT (ML, IHOLDL, JHOLD, IPNTR, ADDLNK)
            IF (IPNTR .GT. 0) THEN
               JJ = IPNTR
            ELSE
               JJ = N19 + 1
               N19 = JJ
               WRITE ( * , 10000)JHOLD, JJ
               ADDLNK = .TRUE.
               CALL LTSORT (ML, IHOLDL, JHOLD, JJ, ADDLNK)
            ENDIF
         ENDIF
      ENDIF
C
C  INPUT THE LINE DATA
C
      N2 = N2 + 1
      J = N2
      IF (J .GT. ML)RETURN
      ADDLNK = .TRUE.
      CALL LTSORT (ML, LINKL, JJ, J, ADDLNK)
      ILINE (J) = JJ
      LTYPE (J) = LTYP
      LCON (1, J) = IP1
      LCON (2, J) = IP2
      LCON (3, J) = IP3
      IF (NN .LT. 0) THEN
        NINT(J) = -NN
        WRITE (*, 10010) J
      ELSE
        NINT (J) = NN
      END IF
      FACTOR (J) = FACT
      IF (FACTOR (J) .LE. 0.)FACTOR (J) = 1.
      ILBOUN (J) = 0
      ISBOUN (J) = 0
      NOROOM = .FALSE.
      RETURN
C
10000 FORMAT ('   OLD LINE NO:', I5, '   TO NEW LINE NO:', I5)
10010 FORMAT ('WARNING: Intervals on line ', I5, ' are negative.',
     &  ' Changed to positive.')
      END
