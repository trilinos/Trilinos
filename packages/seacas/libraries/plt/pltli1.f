C Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
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

C $Id: pltli1.f,v 1.1 1993/07/16 16:48:37 gdsjaar Exp $ 
C $Log: pltli1.f,v $
C Revision 1.1  1993/07/16 16:48:37  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE PLTLI1(PLL,PUR,N,XV,YV,NO,XVO,YVO)
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='PLTLI1')
      DIMENSION PLL(2),PUR(2),XV(*),YV(*),XVO(*),YVO(*)
      DIMENSION P(2),S(2)
      LOGICAL INSIDE

      NMAX = NO
      NO = 0
      S(1) = XV(N)
      S(2) = YV(N)
      DO 2220 I = 1,N
         P(1) = XV(I)
         P(2) = YV(I)
         INSIDE = P(2) .GE. PLL(2)
         IF (INSIDE) THEN
            INSIDE = S(2) .GE. PLL(2)
            IF (INSIDE) THEN
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1)
               YVO(NO) = P(2)

            ELSE
               TEMP = PUR(1) - PLL(1)
               FP = (S(2)-PLL(2))*TEMP
               FQ = (P(2)-PLL(2))*TEMP
               XL = FQ/ (FQ-FP)
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1) + XL* (S(1)-P(1))
               YVO(NO) = P(2) + XL* (S(2)-P(2))
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1)
               YVO(NO) = P(2)
            END IF

         ELSE
            INSIDE = S(2) .GE. PLL(2)
            IF (INSIDE) THEN
               TEMP = PUR(1) - PLL(1)
               FP = (S(2)-PLL(2))*TEMP
               FQ = (P(2)-PLL(2))*TEMP
               XL = FQ/ (FQ-FP)
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1) + XL* (S(1)-P(1))
               YVO(NO) = P(2) + XL* (S(2)-P(2))
            END IF

         END IF

         S(1) = P(1)
         S(2) = P(2)
 2220 CONTINUE
      RETURN

      END
