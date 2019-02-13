C    Copyright(C) 1988-2017 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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

C $Id: locout.f,v 1.1 1991/02/21 15:43:58 gdsjaar Exp $
C $Log: locout.f,v $
C Revision 1.1  1991/02/21 15:43:58  gdsjaar
C Initial revision
C
      SUBROUTINE LOCOUT (TYPE, NDIM, NODEL, TOLER, SORT, P1, P2, BOUND)
      DIMENSION P1(NDIM), P2(NDIM), TOLER(2)
      CHARACTER*(*) NODEL, BOUND, SORT, TYPE
      CHARACTER*16  BNAME
      include 'nu_io.blk'
C
      IF (NODEL(:1) .EQ. 'E') THEN
         NODEL = 'Elements'
      ELSE
         NODEL = 'Nodes'
      ENDIF
C
      IF (BOUND(:3) .EQ. 'BOU') THEN
         BNAME = 'Bounded Search'
      ELSE
         BNAME = 'Unbounded Search'
      END IF
C
      DO 10 IO=IOMIN, IOMAX
         WRITE (IO, 20) NODEL(:LENSTR(NODEL)), TOLER(1), TOLER(2),
     *      TYPE(:LENSTR(TYPE))
   20 FORMAT (//' Locating all ',A,' at a distance ',1PE15.8,
     *   ' plus/minus ',1PE15.8,/' from the ',A)
C
      IF (TYPE .EQ. 'LINE') THEN
         IF (NDIM .EQ. 2) THEN
            WRITE (IO, 30) (P1(I),I=1,NDIM), (P2(I),I=1,NDIM), BNAME
         ELSE
            WRITE (IO, 40) (P1(I),I=1,NDIM), (P2(I),I=1,NDIM), BNAME
         END IF
      ELSE IF (TYPE .EQ. 'POINT') THEN
         IF (NDIM .EQ. 2) THEN
            WRITE (IO, 50) (P1(I), I=1, NDIM)
         ELSE
            WRITE (IO, 60) (P1(I), I=1, NDIM)
         END IF
      ELSE IF (TYPE .EQ. 'PLANE') THEN
      A = P2(1)
      B = P2(2)
      C = P2(3)
      D = A * P1(1) + B * P1(2) + C * P1(3)
C
            WRITE (IO, 70) A, B, C, D
      END IF
      WRITE (IO, 80) SORT(:LENSTR(SORT))
   10 CONTINUE

   30 FORMAT (' from Point (',2(1PE11.3),')',/
     *        ' to   Point (',2(1PE11.3),')',2X,A)
   40 FORMAT (' from Point (',3(1PE11.3),')',/
     *        ' to   Point (',3(1PE11.3),')',2X,A)
   50 FORMAT (' (',2(1PE11.3),')')
   60 FORMAT (' (',3(1PE11.3),')')
   70 FORMAT (' ',1PE15.8,' X + ',1PE15.8,' Y + ',1PE15.8,
     *    ' Z = ',1PE15.8)
   80 FORMAT (' Sorted on field ',A,/)
C
      RETURN
      END
