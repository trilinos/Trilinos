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

C $Id: gapout.f,v 1.1 1991/02/21 15:43:16 gdsjaar Exp $
C $Log: gapout.f,v $
C Revision 1.1  1991/02/21 15:43:16  gdsjaar
C Initial revision
C
      SUBROUTINE GAPOUT (DIRCOS, MASSLV, NUMNIQ, NDIM, IFLGM, IFLGS,
     *   GMTHD)
      include 'nu_io.blk'
      DIMENSION DIRCOS(NDIM + 2,*), MASSLV(2,*)
      CHARACTER*8 GMTHD, COSLAB(3), DISLAB(2)
      LOGICAL ISABRT
      DATA COSLAB /'Cosine X','Cosine Y','Cosine Z'/
      DATA DISLAB /' Normal ','Tangent '/
      LENPAG = 50
C
      DO 30 IO=IOMIN,IOMAX
         DO 20 IPAG = 1, NUMNIQ, LENPAG
            IF (ISABRT()) RETURN
            WRITE (IO, 40) IFLGM, IFLGS, GMTHD, (COSLAB(I),I=1,NDIM),
     *         DISLAB(1), DISLAB(2)
C
            DO 10 I=IPAG, MIN(IPAG+LENPAG-1,NUMNIQ)
               WRITE (IO, 50) I, (MASSLV(J,I),J=1,2),
     *            (DIRCOS(J,I),J=1,NDIM+2)
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
C
   40 FORMAT ('1',
     *   ' Master Flag = ',I6,', Slave Flag = ',I6,', Method: ',A8,//,
     *   '     #   Master   Slave   ',5(A8,3X))
   50 FORMAT (1X,I5,':',2I8,5(1PE11.3))
      RETURN
      END
