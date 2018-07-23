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

C $Id: mpmul4.f,v 1.4 1993/07/19 17:06:30 gdsjaar Exp $ 
C $Log: mpmul4.f,v $
C Revision 1.4  1993/07/19 17:06:30  gdsjaar
C Changed hex constants back to preceding X, --needed on cray. Works
C either way on other systems.
C
c Revision 1.3  1993/07/16  22:11:18  gdsjaar
c Unrolled do loops to speed up execution.
c
c Revision 1.2  1993/07/16  17:33:06  gdsjaar
c Integer constant too big on sun, replaced it with hexadecimal notation
c
c Revision 1.1  1993/07/16  16:47:20  gdsjaar
c Changed plt to library rather than single source file.
c 
C=======================================================================
      SUBROUTINE MPMUL4(N,MASK,ARR1,ARR2,ARR3,ARR4,MAT,RES1,RES2,RES3,
     *                  RES4)
      DIMENSION ARR1(*),ARR2(*),ARR3(*),ARR4(*),RES1(*),
     *          RES2(*),RES3(*),RES4(*),MAT(4,4)
      REAL MAT
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, X'80000000'/

      IF (MASK.EQ.0) THEN
         RETURN

      END IF

      DO 3160 I = 1,N
         IF (IAND(MASK,IZBIT(I)).NE.0) THEN
           RES1(I) = MAT(1,1)*ARR1(I) + MAT(2,1)*ARR2(I) +
     *       MAT(3,1)*ARR3(I) + MAT(4,1)*ARR4(I)
           RES2(I) = MAT(1,2)*ARR1(I) + MAT(2,2)*ARR2(I) +
     *       MAT(3,2)*ARR3(I) + MAT(4,2)*ARR4(I)
           RES3(I) = MAT(1,3)*ARR1(I) + MAT(2,3)*ARR2(I) +
     *       MAT(3,3)*ARR3(I) + MAT(4,3)*ARR4(I)
           RES4(I) = MAT(1,4)*ARR1(I) + MAT(2,4)*ARR2(I) +
     *       MAT(3,4)*ARR3(I) + MAT(4,4)*ARR4(I)
         END IF

 3160 CONTINUE
      RETURN

      END
