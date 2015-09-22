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

C $Id: pltncf.f,v 1.1 1993/07/16 16:48:59 gdsjaar Exp $ 
C $Log: pltncf.f,v $
C Revision 1.1  1993/07/16 16:48:59  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE PLTNCF(X,TYPE,FN,NE)
      CHARACTER*(*) TYPE
      CHARACTER*1 TTYPE
      REAL FNICE(17)
      DATA FNICE/-10.,-8.,-6.,-5.,-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.,5.,6.,
     *     8.,10./

      TTYPE = TYPE
      CALL CHRUP(TTYPE,TTYPE)
      IF (X.EQ.0.) THEN
         FN = 0.
         RETURN

      END IF

      F1 = X/10.**NE
      IF (TTYPE.EQ.'O') THEN
         DO 2870 I = 1,17
            IF (F1.LE.FNICE(I)) THEN
               GO TO 2880

            END IF

 2870    CONTINUE
 2880    CONTINUE
      END IF

      IF (TTYPE.EQ.'U') THEN
         DO 2890 I = 17,1,-1
            IF (F1.GE.FNICE(I)) THEN
               GO TO 2900

            END IF

 2890    CONTINUE
 2900    CONTINUE
      END IF

      IF (I.GT.17) THEN
         I = 17
      END IF

      IF (I.LT.1) THEN
         I = 1
      END IF

      FN = FNICE(I)
      RETURN

      END
