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
      SUBROUTINE PXN2 (RETMIN, PARTYP, LENARY, NELBLK, IXELB, ISEVOK,
     &           NUMIX, IXNODE, IXELEM, SIG11, SIG22, SIG12, RESULT)
C=======================================================================

C   --*** PXN2 *** (ALGEBRA) Calculate 2x2 principal values
C   --   Written by Amy Gilkey - revised 05/17/88
C   --
C   --PXN2 determines minimum or maximum principal values for a 2 by 2
C   --array.
C   --
C   --Parameters:
C   --   RETMIN - IN - true iff minimum versus maximum to be calculated
C   --   PARTYP - IN - the parameter type (element must be handled differently)
C   --   LENARY - IN - the length of the parameter arrays
C   --   NELBLK - IN - the number of element blocks
C   --   IXELB  - IN - the cumulative element counts for each element block
C   --   ISEVOK - IN - the element variable truth table for each element block
C   --   NUMIX  - IN - the number of selected values; <0 if all
C   --   IXNODE - IN - the indices of the selected nodes (only if NUMIX >= 0)
C   --   IXELEM - IN - the indices of the selected elements (only if NUMIX >= 0)
C   --   SIG11, SIG22, SIG12 - IN - the tensor components
C   --   RESULT - OUT - the returned principal values

      LOGICAL RETMIN
      CHARACTER PARTYP
      INTEGER IXELB(0:NELBLK)
      LOGICAL ISEVOK(*)
      INTEGER IXNODE(*)
      INTEGER IXELEM(*)
      REAL SIG11(*), SIG22(*), SIG12(*)
      REAL RESULT(*)

      IF (NUMIX .GE. 0) THEN
         IF (PARTYP .NE. 'E') THEN
            DO 100 I = 1, NUMIX
               J = IXNODE(I)
               SROOT = SQRT
     &            ((0.5 * (SIG11(J)-SIG22(J)))**2 + SIG12(J)**2)

               IF (RETMIN) THEN
                  RESULT(J) = 0.5 * (SIG11(J)+SIG22(J)) - SROOT
               ELSE
                  RESULT(J) = 0.5 * (SIG11(J)+SIG22(J)) + SROOT
               END IF
  100       CONTINUE
         ELSE
            DO 120 IELB = 1, NELBLK
               IF (ISEVOK(IELB)) THEN
                  DO 110 I = IXELB(IELB-1)+1, IXELB(IELB)
                     J = IXELEM(I)
                     SROOT = SQRT
     &                  ((0.5 * (SIG11(J)-SIG22(J)))**2 + SIG12(J)**2)

                     IF (RETMIN) THEN
                        RESULT(J) = 0.5 * (SIG11(J)+SIG22(J)) - SROOT
                     ELSE
                        RESULT(J) = 0.5 * (SIG11(J)+SIG22(J)) + SROOT
                     END IF
  110             CONTINUE
               END IF
  120       CONTINUE
         END IF

      ELSE
C        NUMIX < 0; ALL VALUES SELECTED
C        NUMEL=NUMELO; NUMNP=NUMNPO; IXELB(J)=IXELBO(J), J=1,NELBLK
         IF (PARTYP .NE. 'E') THEN
            DO 130 J = 1, LENARY
               SROOT = SQRT
     &            ((0.5 * (SIG11(J)-SIG22(J)))**2 + SIG12(J)**2)

               IF (RETMIN) THEN
                  RESULT(J) = 0.5 * (SIG11(J)+SIG22(J)) - SROOT
               ELSE
                  RESULT(J) = 0.5 * (SIG11(J)+SIG22(J)) + SROOT
               END IF
  130       CONTINUE

         ELSE
            DO 150 IELB = 1, NELBLK
               IF (ISEVOK(IELB)) THEN
                  DO 140 J = IXELB(IELB-1)+1, IXELB(IELB)
                     SROOT = SQRT
     &                  ((0.5 * (SIG11(J)-SIG22(J)))**2 + SIG12(J)**2)

                     IF (RETMIN) THEN
                        RESULT(J) = 0.5 * (SIG11(J)+SIG22(J)) - SROOT
                     ELSE
                        RESULT(J) = 0.5 * (SIG11(J)+SIG22(J)) + SROOT
                     END IF
  140             CONTINUE
               END IF
  150       CONTINUE
         END IF
      END IF

      RETURN
      END
