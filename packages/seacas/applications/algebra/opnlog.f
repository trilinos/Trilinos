C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
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
C    * Neither the name of Sandia Corporation nor the names of its
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

C=======================================================================
      SUBROUTINE OPNLOG (LOGU)
C=======================================================================
C $Id: opnlog.f,v 1.8 2008/03/14 13:45:28 gdsjaar Exp $
C   --*** OPNLOG *** (BLOT) Open log file and write header
C   --   Written by Amy Gilkey - revised 12/21/87
C   --
C   --OPNLOG opens the log file and writes the command line as the header
C   --for the log file.
C   --
C   --Parameters:
C   --   NLOG - IN - the log file number
C   --
C   --Common Variables:
C   --   Uses QAINFO of /PROGQA/
C   --   Uses NDBIN, NDBOUT of /DBASE/

      include 'params.blk'
      include 'progqa.blk'
      include 'dbase.blk'

      CHARACTER*256 INLINE
      CHARACTER*256 STR

      NLOG = LOGU
      CALL OPNFIL (NLOG, 'U', 'L', 0, IERR)
      IF (IERR .NE. 0) THEN
         CALL PRTERR ('WARNING', 'Log file cannot be opened')
         NLOG = -1
         GOTO 100
      END IF

      INLINE = '$$$ ' // QAINFO(1)
      L = LENSTR (INLINE) + 1

      CALL EXNAME(NDBIN, STR, LFIL)
      IF (L .LT. LEN (INLINE)) INLINE(L+1:) = STR(:LFIL)
      L = LENSTR (INLINE) + 1

      CALL EXNAME(NDBOUT, STR, LFIL)
      IF (L .LT. LEN (INLINE)) INLINE(L+1:) = STR(:LFIL)
      L = LENSTR (INLINE) + 1

      WRITE (NLOG, '(A)') INLINE(:L-1)

  100 CONTINUE
      LOGU = NLOG
      RETURN
      END
