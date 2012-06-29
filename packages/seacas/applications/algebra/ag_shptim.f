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
      SUBROUTINE SHPTIM (NPTIMS, IPTIMS, TIMES)
C=======================================================================

C   --*** SHPTIM *** (TIMSEL) Display time request
C   --   Written by Amy Gilkey - revised 02/01/88
C   --
C   --SHPTIM displays the selected times.
C   --
C   --Parameters:
C   --   NPTIMS - IN - the number of selected times
C   --   IPTIMS - IN - the selected time steps
C   --   TIMES  - IN - the database times
C
C   --Routines Called:
C   --   LENSTR - (STRLIB) Find string length
C   --   NUMSTR - (STRLIB) Convert numbers to engineering notation
C   --   SQZSTR - (STRLIB) Delete extra blanks from string


      INTEGER NPTIMS
      INTEGER IPTIMS(*)
      REAL TIMES(*)

      CHARACTER*5 ISTRA
      CHARACTER*20 RSTR(3)
      REAL RNUM(3)

      WRITE (*, *)

      NPTIMW = NPTIMS

      CALL INTSTR (1, 0, NPTIMS, ISTRA, LSTRA)

      WRITE (*, 10010)
     &      ISTRA(:LSTRA), ' selected time steps'

      IF (NPTIMS .GT. 0) THEN
         MXSTEP = 0
         DO 100 I = 1, NPTIMS
            MXSTEP = MAX (MXSTEP, IPTIMS(I))
  100    CONTINUE
         CALL INTSTR (1, 0, MXSTEP, ISTRA, LSTRA)
         RNUM(2) = TIMES(IPTIMS(1))
         IF ((RNUM(2) .EQ. 0.0) .AND. (NPTIMS .GT. 1))
     &      RNUM(2) = TIMES(IPTIMS(2))
         RNUM(3) = TIMES(IPTIMS(NPTIMS))

         DO 110 I = 1, NPTIMS
            CALL INTSTR (1, LSTRA, IPTIMS(I), ISTRA, L)
            RNUM(1) = TIMES(IPTIMS(I))
            CALL NUMSTR (3, 4, RNUM, RSTR, LSTR)
            WRITE (*, 10000) I,
     &            ISTRA(:L), RSTR(1)(:LSTR)
10000       FORMAT (1X, I5, ')', 3X, '(step ', A, ')', 3X, A,
     &            :, 3X, A)

  110    CONTINUE
      END IF

      RETURN
10010  FORMAT (1X, 5A)
      END
