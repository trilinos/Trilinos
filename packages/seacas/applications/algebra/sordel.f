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
      SUBROUTINE SORDEL (ITOTAL, IBEGIN, IEND)
C=======================================================================

C   --*** SORDEL *** (ALGEBRA) Sort deleted variables to end of entries
C   --   Written by Amy Gilkey - revised 12/11/87
C   --
C   --SORDEL sorts the deleted variables so that they appear at the start
C   --of all non-deleted entries.
C   --
C   --Parameters:
C   --   ITOTAL - IN - the ending /VAR../ index of the variables to sort
C   --   IBEGIN - IN - the starting /VAR../ index of the deleted variables
C   --   IEND - OUT - the ending /VAR../ index of the deleted variables
C   --
C   --Common Variables:
C   --   Sets NAMVAR, TYPVAR, IDVAR, ISTVAR, IEVVAR of /VAR../

      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)
      include 'namlen.blk'
      include 'var.blk'

      CHARACTER*(maxnam) NAMTMP
      CHARACTER TYPTMP
      INTEGER ISTTMP(3)

      IEND = IBEGIN - 1

      DO 100 NVAR = IBEGIN, ITOTAL

         IF (ISTVAR(ICURTM,NVAR) .EQ. -1) THEN

            IEND = IEND + 1

            IF (IEND .NE. NVAR) THEN
               NAMTMP = NAMVAR(NVAR)
               TYPTMP = TYPVAR(NVAR)
               IDTMP = IDVAR(NVAR)
               ISTTMP(1) = ISTVAR(1,NVAR)
               ISTTMP(2) = ISTVAR(2,NVAR)
               ISTTMP(3) = ISTVAR(3,NVAR)
               IEVTMP = IEVVAR(NVAR)

               NAMVAR(NVAR) = NAMVAR(IEND)
               TYPVAR(NVAR) = TYPVAR(IEND)
               IDVAR(NVAR) = IDVAR(IEND)
               ISTVAR(1,NVAR) = ISTVAR(1,IEND)
               ISTVAR(2,NVAR) = ISTVAR(2,IEND)
               ISTVAR(3,NVAR) = ISTVAR(3,IEND)
               IEVVAR(NVAR) = IEVVAR(IEND)

               NAMVAR(IEND) = NAMTMP
               TYPVAR(IEND) = TYPTMP
               IDVAR(IEND) = IDTMP
               ISTVAR(1,IEND) = ISTTMP(1)
               ISTVAR(2,IEND) = ISTTMP(2)
               ISTVAR(3,IEND) = ISTTMP(3)
               IEVVAR(IEND) = IEVTMP
            END IF
         END IF
  100 CONTINUE

      RETURN
      END
