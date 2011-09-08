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
      SUBROUTINE NXTTIM ()
C=======================================================================

C   --*** NXTTIM *** (ALGEBRA) Move current time variables to last time
C   --   Written by Amy Gilkey - revised 12/10/87
C   --
C   --NXTTIM moves the values of the variables for the current time to the
C   --location for the last time.  The move is implemented by swaping
C   --the location pointers for the last and current time.
C   --
C   --Parameters:
C   --
C   --Common Variables:
C   --   Sets ISTVAR of /VAR../
C   --   Uses NUMINP, IXLHS of /VAR../

      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)
      
      include 'namlen.blk'
      include 'var.blk'

c      LOGICAL WHOTIM

c      LOGICAL ISHIST

C   --If previous time requested, swap pointers for current and previous time

      DO 100 IVAR = 1, NUMINP
c         ISHIST = ((TYPVAR(IVAR) .EQ. 'H') .OR. (TYPVAR(IVAR) .EQ. 'T'))
c         IF (ISHIST .OR. WHOTIM) THEN
         IF (ISTVAR(ILSTTM,IVAR) .GT. 0) THEN
            I = ISTVAR(ICURTM,IVAR)
            ISTVAR(ICURTM,IVAR) = ISTVAR(ILSTTM,IVAR)
            ISTVAR(ILSTTM,IVAR) = I
         END IF
c         END IF
  100 CONTINUE

      DO 110 IVAR = IXLHS, MAXVAR
c         ISHIST = ((TYPVAR(IVAR) .EQ. 'H') .OR. (TYPVAR(IVAR) .EQ. 'T'))
c         IF (ISHIST .OR. WHOTIM) THEN
         IF (ISTVAR(ILSTTM,IVAR) .GT. 0) THEN
            I = ISTVAR(ICURTM,IVAR)
            ISTVAR(ICURTM,IVAR) = ISTVAR(ILSTTM,IVAR)
            ISTVAR(ILSTTM,IVAR) = I
         END IF
c         END IF
  110 CONTINUE

      RETURN
      END
