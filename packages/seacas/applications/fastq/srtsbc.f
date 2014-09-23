C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Governement retains certain rights in this software.
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

C $Id: srtsbc.f,v 1.2 2000/11/13 15:39:05 gdsjaar Exp $
C $Log: srtsbc.f,v $
C Revision 1.2  2000/11/13 15:39:05  gdsjaar
C Cleaned up unused variables and labels.
C
C Removed some real to int conversion warnings.
C
C Revision 1.1.1.1  1990/11/30 11:16:37  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:16:36  gdsjaar
c Initial revision
c 
C
CC* FILE: [.RENUM]SRTSBC.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SRTSBC (MXSFLG, NPSBC, NPELEM, NNXK, NXK, NSFLG, NSLEN,
     &   NSPTR, NVLEN, NVPTR, LSTSBC, NELEMS, NSIDEN, NSIDES, NNSBC,
     &   NSLIST, NVLIST, NBCSID)
C***********************************************************************
C
C  SUBROUTINE SRTSBC = SORTS THE LIST OF SIDE BOUNDARY CARDS
C
C***********************************************************************
C
C  VARIABLES USED:
C     NSFLG  = THE ARRAY OF FLAG VALUES
C     NSLEN  = NUMBER OF ELEMENTS IN ILIST ASSOCIATED WITH EACH FLAG
C     NSPTR  = POINTER TO THE FIRST ELEMENT IN LIST FOR EACH FLAG
C     ILIST  = THE ELEMENT LIST
C     KKK    = THE NUMBER OF ELEMENTS IN THE MESH
C     MXSFLG = THE NUMBER OF ENTRIES IN THE BOUNDARY LIST
C     FOUND  = .TRUE. IF A NEW UNIQUE FLAG HAS BEEN FOUND
C
C***********************************************************************
C
      DIMENSION NXK (NNXK, NPELEM)
      DIMENSION NSFLG (MXSFLG), NSLEN (MXSFLG), NSPTR (MXSFLG)
      DIMENSION NVLEN (MXSFLG), NVPTR (MXSFLG)
      DIMENSION LSTSBC (NPSBC), NELEMS (NPSBC)
      DIMENSION NSIDES (NPSBC), NSIDEN (NPSBC)
C
      LOGICAL FOUND
C
      IFLAG  = -1
      NSLIST = 0
      NBCSID = 0
      IBEGIN = 1
C
  100 CONTINUE
      FOUND = .FALSE.
C
      DO 110 I = IBEGIN, NNSBC, 3
         IF (LSTSBC (I) .LT. 0) THEN
            IF (FOUND) THEN
               IF (IFLAG .EQ. ABS (LSTSBC (I))) THEN
                  NSLIST = NSLIST + 1
                  NELEMS (NSLIST) = LSTSBC (I + 1)
                  NSIDES (NSLIST) = LSTSBC (I + 2)
                  NSLEN (NBCSID) = NSLEN (NBCSID) + 1
                  LSTSBC (I) = 0
               ENDIF
            ELSE
               FOUND = .TRUE.
               NBCSID = NBCSID + 1
               IFLAG =  - LSTSBC (I)
               NSFLG (NBCSID) = IFLAG
               NSLEN (NBCSID) = 1
               NSLIST = NSLIST + 1
               NSPTR (NBCSID) = NSLIST
               NELEMS (NSLIST) = LSTSBC (I + 1)
               NSIDES (NSLIST) = LSTSBC (I + 2)
               LSTSBC (I) = 0
               IBEGIN = I
            ENDIF
         ENDIF
  110 CONTINUE
C
      IF (FOUND) THEN
         GOTO 100
      ELSE
C
C  PUT ALL THE NODES ATTACHED TO THE ELEMENT BCC INTO THE
C  NSIDEN LIST
C
         NVLIST = 0
         DO 130 I = 1, NBCSID
            ISTART = NSPTR (I)
            IEND = NSPTR (I) + NSLEN (I) - 1
            NVPTR (I) = NVLIST + 1
            DO 120 J = ISTART, IEND
               J1 = NSIDES (J)
               J2 = J1 + 1
               IF (J2 .EQ. 5) J2 = 1
               NVLIST = NVLIST + 1
               NSIDEN (NVLIST) = NXK (J1, NELEMS (J))
               NVLIST = NVLIST + 1
               NSIDEN (NVLIST) = NXK (J2, NELEMS (J))
  120       CONTINUE
            NVLEN (I) = NVLIST - NVPTR (I) + 1
  130    CONTINUE
         RETURN
      ENDIF
C
      END
