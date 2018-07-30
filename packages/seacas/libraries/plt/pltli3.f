C Copyright (C) 2009-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C     * Neither the name of NTESS nor the names of its
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

C $Id: pltli3.f,v 1.2 2000/10/25 18:55:02 gdsjaar Exp $ 
C $Log: pltli3.f,v $
C Revision 1.2  2000/10/25 18:55:02  gdsjaar
C In the pltli? functions, check for N==0 before doing any array
C accesses.
C
C Also changed all references to 'mask' to be arrays where they were
C scalars since downstream code seems to treat them as arrays.
C
C Revision 1.1  1993/07/16 16:48:39  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE PLTLI3(PLL,PUR,N,XV,YV,NO,XVO,YVO)
      DIMENSION PLL(2),PUR(2),XV(*),YV(*),XVO(*),YVO(*)
      DIMENSION P(2),S(2)
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='PLTLI3')
      LOGICAL INSIDE

      NMAX = NO
      NO = 0
      if (n .gt. 0) then
        S(1) = XV(N)
        S(2) = YV(N)
        DO 2260 I = 1,N
          P(1) = XV(I)
          P(2) = YV(I)
          INSIDE = P(2) .LE. PUR(2)
          IF (INSIDE) THEN
            INSIDE = S(2) .LE. PUR(2)
            IF (INSIDE) THEN
              NO = NO + 1
              IF (NO.GT.NMAX) THEN
                CALL PLTFLU
                CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *            ,2)
                RETURN

              END IF

              XVO(NO) = P(1)
              YVO(NO) = P(2)

            ELSE
              TEMP = PLL(1) - PUR(1)
              FP = (S(2)-PUR(2))*TEMP
              FQ = (P(2)-PUR(2))*TEMP
              XL = FQ/ (FQ-FP)
              NO = NO + 1
              IF (NO.GT.NMAX) THEN
                CALL PLTFLU
                CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *            ,2)
                RETURN

              END IF

              XVO(NO) = P(1) + XL* (S(1)-P(1))
              YVO(NO) = P(2) + XL* (S(2)-P(2))
              NO = NO + 1
              IF (NO.GT.NMAX) THEN
                CALL PLTFLU
                CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *            ,2)
                RETURN

              END IF

              XVO(NO) = P(1)
              YVO(NO) = P(2)
            END IF

          ELSE
            INSIDE = S(2) .LE. PUR(2)
            IF (INSIDE) THEN
              TEMP = PLL(1) - PUR(1)
              FP = (S(2)-PUR(2))*TEMP
              FQ = (P(2)-PUR(2))*TEMP
              XL = FQ/ (FQ-FP)
              NO = NO + 1
              IF (NO.GT.NMAX) THEN
                CALL PLTFLU
                CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *            ,2)
                RETURN

              END IF

              XVO(NO) = P(1) + XL* (S(1)-P(1))
              YVO(NO) = P(2) + XL* (S(2)-P(2))
            END IF

          END IF

          S(1) = P(1)
          S(2) = P(2)
 2260   CONTINUE
      end if
      RETURN

      END
