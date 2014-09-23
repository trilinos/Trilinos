C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Governement retains certain rights in this software.
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

C $Id: shwblk.f,v 1.1 1991/02/21 15:45:39 gdsjaar Exp $
C $Log: shwblk.f,v $
C Revision 1.1  1991/02/21 15:45:39  gdsjaar
C Initial revision
C
      SUBROUTINE SHWBLK (NELBLK, MAT, NSELND, NSELEL)
      DIMENSION MAT(6, NELBLK)
      CHARACTER*16 TYPE
      CHARACTER*80 STRTMP
      include 'nu_io.blk'
C
      DO 10 IO=IOMIN, IOMAX
         WRITE (IO, 50)
   10 CONTINUE
      DO 30 ITMP=1, NELBLK
         I = MAT(6, ITMP)
         IF (MAT(5,I) .EQ. 1) THEN
            TYPE = 'Selected'
         ELSE
            TYPE = 'Not Selected'
         END IF
         DO 20 IO=IOMIN, IOMAX
            WRITE (IO, 60) MAT(1,I), MAT(3,I), MAT(4,I), MAT(2,I), I,
     *         TYPE
   20    CONTINUE
   30 CONTINUE
      WRITE (STRTMP, 70) NSELND, NSELEL
      CALL SQZSTR(STRTMP, LSTR)
      DO 40 IO=IOMIN, IOMAX
         WRITE (IO, 80) STRTMP(:LSTR)
   40 CONTINUE

      RETURN
   50 FORMAT (//
     *   '     Material  First     Last   Number of     Block'/,
     *   '      Number  Element   Element  Elements ')
   60 FORMAT (5I10,5X,A16)
   70 FORMAT (I8,' nodes and ',I8,' elements selected')
   80 FORMAT (/5X,A)
      END
