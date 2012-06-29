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

C $Id: pltav2.f,v 1.4 2000/10/25 18:55:01 gdsjaar Exp $ 
C $Log: pltav2.f,v $
C Revision 1.4  2000/10/25 18:55:01  gdsjaar
C In the pltli? functions, check for N==0 before doing any array
C accesses.
C
C Also changed all references to 'mask' to be arrays where they were
C scalars since downstream code seems to treat them as arrays.
C
C Revision 1.3  1993/07/19 17:06:31  gdsjaar
C Changed hex constants back to preceding X, --needed on cray. Works
C either way on other systems.
C
c Revision 1.2  1993/07/16  17:33:07  gdsjaar
c Integer constant too big on sun, replaced it with hexadecimal notation
c
c Revision 1.1  1993/07/16  16:47:40  gdsjaar
c Changed plt to library rather than single source file.
c 
C=======================================================================
      SUBROUTINE PLTAV2(UMAP,N,X1,Y1,X2,Y2,TH,XL)
      REAL UMAP(*)
      INTEGER N
      REAL X1(*),Y1(*)
      REAL X2(*),Y2(*)
      REAL TH
      REAL XL
      REAL PX(32),PY(32),QX(32),QY(32)
      INTEGER MASK(1)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, X'80000000'/

      MASK(1) = -1
      CALL PLTMV2(UMAP,N,MASK,X1,Y1,X2,Y2,PX,PY,QX,QY)
      DO 2000 I = 1,N
         IF (IAND(MASK(1),IZBIT(I)).NE.0) THEN
            CALL PLTARR(PX(I),PY(I),QX(I),QY(I),TH,XL)
         END IF

 2000 CONTINUE
      RETURN

      END
