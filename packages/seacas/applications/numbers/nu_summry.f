C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Government retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C    
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C    
C        * Neither the name of Sandia Corporation nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
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

C $Id: summry.f,v 1.2 2000/07/06 16:49:57 gdsjaar Exp $
C $Log: summry.f,v $
C Revision 1.2  2000/07/06 16:49:57  gdsjaar
C Changed real*4 to real
C
C Revision 1.1.1.1  1991/02/21 15:45:55  gdsjaar
C NUMBERS: Greg Sjaardema, initial Unix release
C
c Revision 1.1  1991/02/21  15:45:54  gdsjaar
c Initial revision
c
      SUBROUTINE SUMMRY (TYPE, NUM, SELECT, VALUE, SUMR, ISUMR, IOFF)
      CHARACTER*(*) TYPE
      LOGICAL SELECT(*)
      REAL VALUE(*), SUMR(*)
      INTEGER ISUMR(*)
C
C ... TYPE = 'A' - calculate stats based on absolute values
C              NOTE: VALUE will be modified if this option used
C          = ' ' - calculate stats based on true value
C
C ... SUMR(1) = MINIMUM   ISUMR(1) = ELEMENT NUMBER
C ... SUMR(2) = MAXIMUM   ISUMR(2) = ELEMENT NUMBER
C ... SUMR(3) = AVERAGE
C ... SUMR(4) = STD. DEV.
C
      SUMR(1) =  1.0E30
      SUMR(2) = -1.0E30
      NUMSEL = 0

      RMEAN  = 0.0
      STDDEV = 0.0

      IF (TYPE(:1) .EQ. 'A') THEN
         DO 5 I=1, NUM
            VALUE(I) = ABS(VALUE(I))
    5    CONTINUE
      END IF

      DO 10 I = 1, NUM
         IF (SELECT(I)) THEN
            NUMSEL = NUMSEL + 1
            TMEAN  = RMEAN + (VALUE(I) - RMEAN) / NUMSEL
            STDDEV = STDDEV + (VALUE(I) - RMEAN) * (VALUE(I) - TMEAN)
            RMEAN  = TMEAN

            IF (VALUE(I) .LT. SUMR(1)) THEN
               SUMR(1) = VALUE(I)
               ISUMR(1) = I
            END IF

            IF (VALUE(I) .GT. SUMR(2)) THEN
               SUMR(2) = VALUE(I)
               ISUMR(2) = I
            END IF

         END IF
   10 CONTINUE

      SUMR(3) = RMEAN
      SUMR(4) = SQRT(STDDEV / MAX(1.0, FLOAT(NUMSEL-1)))

      ISUMR(1) = ISUMR(1) + IOFF
      ISUMR(2) = ISUMR(2) + IOFF

      RETURN
      END
