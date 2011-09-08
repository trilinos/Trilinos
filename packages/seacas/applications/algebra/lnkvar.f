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
      SUBROUTINE LNKVAR (*)
C=======================================================================

C   --*** LNKVAR *** (ALGEBRA) Link equation references with variables
C   --   Written by Amy Gilkey - revised 11/30/87
C   --
C   --LNKVAR links the variables in the equation with the /VAR../ arrays.
C   --Each equation variable (from /ENT../) is found in the /VAR../
C   --arrays, and the /VAR../ index is stored in the INXENT array.
C   --Note that since the /VAR../ arrays were constructed from the
C   --equation variables, each variable will be found.
C   --
C   --Parameters:
C   --   * - return statement if an error is found; message is printed
C   --
C   --Common Variables:
C   --   Sets INXENT, VALENT of /ENT../
C   --   Uses NAMENT, TYPENT of /ENT../
C   --   Uses NAMVAR, TYPVAR, IDVAR, ISTVAR, NUMINP, IXLHS of /VAR../

      include 'namlen.blk'
      include 'params.blk'
      include 'numeqn.blk'
      include 'ent.blk'
      include 'var.blk'
      
      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)

      NUMLHS = MAXVAR - IXLHS + 1

      DO 110 NEQN = 1, NUMEQN
         DO 100 NENT = 1, NUMENT(NEQN)
            IF (TYPENT(NENT,NEQN) .EQ. 'V') THEN

C            --If searching for a LHS variable (INXENT < 0),
C            --   then search the assigned variables, NAMVAR(IXLHS..MAXVAR)
C            --   else search the input variables, NAMVAR(1..NUMINP)

               IF (INXENT(NENT,NEQN) .LT. 0) THEN
                  IV = LOCSTR (NAMENT(NENT,NEQN),
     &               NUMLHS, NAMVAR(IXLHS))
                  IF (IV .GT. 0) IV = IXLHS + IV - 1
               ELSE
                  IV = LOCSTR (NAMENT(NENT,NEQN), NUMINP, NAMVAR(1))
               END IF

               IF (IV .LE. 0) THEN
                  CALL PRTERR ('PROGRAM', 'Variable "'
     &               // NAMENT(NENT,NEQN) // '" cannot be found')
                  GOTO 120
               END IF

               INXENT(NENT,NEQN) = IV

               IF ((TYPVAR(IV) .EQ. 'T')
     &            .OR. (TYPVAR(IV) .EQ. 'G')) THEN
c               IF ((TYPVAR(IV) .EQ. 'T')
c     &            .OR. (TYPVAR(IV) .EQ. 'H')
c     &            .OR. (TYPVAR(IV) .EQ. 'G')) THEN
                  IF (NAMVAR(IV) .EQ. '.GLOBAL') THEN
c                  IF ((NAMVAR(IV) .EQ. '.HISTORY')
c     &               .OR. (NAMVAR(IV) .EQ. '.GLOBAL')) THEN
                     VALENT(NENT,NEQN) =
     &                  VALENT(NENT,NEQN) + IDVAR(IV) - 1
                  ELSE
                     VALENT(NENT,NEQN) = IDVAR(IV)
                  END IF
               END IF

            END IF
  100    CONTINUE
  110 CONTINUE
      RETURN

  120 CONTINUE
      RETURN 1
      END
