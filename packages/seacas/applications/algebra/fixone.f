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
      SUBROUTINE FIXONE (MAXNE, VARVAL)
C=======================================================================

C   --*** FIXONE *** (ALGEBRA) Handle first time step variables
C   --   Written by Amy Gilkey - revised 12/10/87
C   --
C   --FIXONE is called after the first time step is read in.  It copies
C   --the current variables to the last and first time step variables
C   --(if needed).  It also corrects the current variable pointer if
C   --only the first time step for a variable is needed.
C   --
C   --Parameters:
C   --   MAXNE  - IN     - the VARVAL dimension
C   --   VARVAL - IN/OUT - the input data and copied data
C   --
C   --Common Variables:
C   --   Uses ISTVAR of /VAR../
C   --   Uses NVARNP, NVAREL of /DBNUMS/

      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)
      
      include 'namlen.blk'
      include 'var.blk'

      REAL VARVAL(MAXNE,*)

      DO 100 IVAR = 1, NUMINP
            IF (ISTVAR(IONETM,IVAR) .GT. 0) THEN
               IF (ISTVAR(ICURTM,IVAR) .LT. 0) THEN
                  ISTVAR(ICURTM,IVAR) = 0
               ELSE
                  IFROM = ISTVAR(ICURTM,IVAR)
                  ITO   = ISTVAR(IONETM,IVAR)
                  CALL CPYVAR (TYPVAR(IVAR), MAXNE,
     &                         VARVAL(1,IFROM), VARVAL(1,ITO))
               END IF
            END IF
            IF (ISTVAR(ILSTTM,IVAR) .GT. 0) THEN
               IFROM = ISTVAR(ICURTM,IVAR)
               ITO = ISTVAR(ILSTTM,IVAR)
               CALL CPYVAR (TYPVAR(IVAR), MAXNE,
     &                      VARVAL(1,IFROM), VARVAL(1,ITO))
            END IF
  100 CONTINUE

      DO 120 IVAR = IXLHS, MAXVAR
            DO 110 ITM = 1, 3
               IF (ITM .NE. ICURTM) THEN
                  IF (ISTVAR(ITM,IVAR) .GT. 0) THEN
                     IFROM = ISTVAR(ICURTM,IVAR)
                     ITO = ISTVAR(ITM,IVAR)
                     CALL CPYVAR (TYPVAR(IVAR), MAXNE,
     &                            VARVAL(1,IFROM), VARVAL(1,ITO))
                  END IF
               END IF
  110       CONTINUE
  120 CONTINUE

      RETURN
      END
