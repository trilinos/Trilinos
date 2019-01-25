C    Copyright(C) 2014-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
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
C    * Neither the name of NTESS nor the names of its
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

C $Id: sortst.f,v 1.1 1990/11/30 11:16:14 gdsjaar Exp $
C $Log: sortst.f,v $
C Revision 1.1  1990/11/30 11:16:14  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]SORTST.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SORTST(N, ARRIN, ITOPP, IBOTP, IRANGE, INDX)
C***********************************************************************
C
C SUBROUTINE SORTST = THIS SUBROUTINE SORTS A REAL ARRAY IN ASCENDING
C                     ORDER  AND DOES IT JUST CHANGING THE INDEX ARRAY.
C
C***********************************************************************
C
C VARIABLES   IN : N ......NUMBER OF ELEMENTS IN THE ARRAY
C                  ARRIN ..REAL ARRAY WITH DATA TO BE SORTED
C                  ITOPP...TOP POINTER IN THE X ARRAY
C                  IBOTP...BOTTON POINTER IN THE X ARRAY
C             OUT: INDX ...INDEX ARRAY WITH ITS SORTED VALUES IN
C                            ASCENDING ORDER.
C                  IRANGE..RANGE BETWEEN THE ARRAY 'INDX' WAS SORTED
C
C WRITTEN BY:  HORACIO RECALDE                 DATE: JAN 20, 1988
C MODIFIED BY: MB STEPHENSON                   DATA: MAR 08, 1989
C   REPLACED EXCHANGE SORT WITH HEAPSORT FOR EFFICIENCY
C***********************************************************************
C
      REAL ARRIN(N)
      INTEGER INDX(N)
C
C...  CHECK POINTERS
C
      ITOP = ITOPP
      IBOT = IBOTP
      IF (ITOP .GT. N) ITOP = ITOP - N
      IF (IBOT .GT. N) IBOT = IBOT - N
C
C---  CALCULATE THE RANGE AND INITIALIZE INDEX ARRAY
C
      IF (ITOP .EQ. IBOT) THEN
         IRANGE = 1
         INDX(1) = ITOP
         RETURN
      ELSE IF (ITOP .LT. IBOT) THEN
         IRANGE = IBOT - ITOP + 1
      ELSE
         IRANGE = IBOT - ITOP + N + 1
      ENDIF
C
      DO 100 J = 1, IRANGE
         INDX(J) = ITOP
         ITOP = ITOP + 1
         IF (ITOP .GT. N) ITOP = 1
  100 CONTINUE
C
C---  PERFORM A HEAPSORT ON THE ELEMENTS
C           (SEE NUMERICAL RECEIPTS, PG. 233)
C           NOTE:  THERE MUST BE AT LEAST 2 ELEMENTS IN THE ARRAY
C
      L = IRANGE/2 + 1
      IR = IRANGE
  110 CONTINUE
      IF (L .GT. 1) THEN
         L = L - 1
         INDXT = INDX(L)
         Q = ARRIN(INDXT)
      ELSE
         INDXT = INDX(IR)
         Q = ARRIN(INDXT)
         INDX(IR) = INDX(1)
         IR = IR - 1
         IF (IR .EQ. 1) THEN
            INDX(1) = INDXT
            RETURN
         END IF
      END IF
C
      I = L
      J = L + L
  120 CONTINUE
      IF (J .LE. IR) THEN
         IF (J .LT. IR) THEN
            IF (ARRIN(INDX(J)) .LT. ARRIN(INDX(J + 1))) J = J + 1
         END IF
         IF (Q .LT. ARRIN(INDX(J))) THEN
            INDX(I) = INDX(J)
            I = J
            J = J + J
         ELSE
            J = IR + 1
         END IF
         GO TO 120
      END IF
      INDX(I) = INDXT
      GO TO 110
C
      END
