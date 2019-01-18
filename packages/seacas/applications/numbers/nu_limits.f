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

C $Id: limits.f,v 1.4 1999/02/16 21:38:00 gdsjaar Exp $
C $Log: limits.f,v $
C Revision 1.4  1999/02/16 21:38:00  gdsjaar
C Converted to read exodusII database format.  Somewhat tested, not
C ready for production yet.
C
C Revision 1.3  1993/07/21 22:34:56  gdsjaar
C Removed unused variable--error
C
c Revision 1.2  1991/02/21  16:37:58  gdsjaar
c Moved ENGNOT function out of write statements
c
c Revision 1.1.1.1  1991/02/21  15:43:50  gdsjaar
c NUMBERS: Greg Sjaardema, initial Unix release
c
c Revision 1.1  1991/02/21  15:43:48  gdsjaar
c Initial revision
c
      SUBROUTINE LIMITS (XYZMIN, XYZMAX, CRD, IX, MAT, NDIM, NEBLK,
     *   NNODES, EXODUS, TIME, ITMSEL, CORDSP, NUMNP)
      DIMENSION CRD(NUMNP, *), IX(NNODES,*), MAT(6,*),
     *   XYZMIN(NDIM,NEBLK), XYZMAX(NDIM,NEBLK),
     *   TIME(*), CORDSP(NUMNP, *)
      LOGICAL ITMSEL(*), ISABRT

      CHARACTER*16 ENGNOT, ENG1
      DIMENSION OVMIN (3), OVMAX (3)
      LOGICAL EXODUS
      include 'nu_io.blk'
      include 'nu_ptim.blk'
C
C ... IF NOT EXODUS, THEN CORDSP CONTAINS COORDINATES
C
      IF (EXODUS) THEN
         CALL GETDSP (CRD, CORDSP, NDIM, NUMNP, TIME, ITMSEL, 
     *      'R', ISTAT)
         IF (ISTAT .NE. 0) GO TO 130
      END IF

   10 CONTINUE
      IF (EXODUS) THEN
         CALL GETDSP (CRD, CORDSP, NDIM, NUMNP, TIME, ITMSEL, 
     *      'A', ISTAT)
         IF (ISTAT .NE. 0) GO TO 130
      END IF
      IF (ISABRT()) RETURN

      IF (NDIM .EQ. 2) THEN
         DO 40 IBLK = 1, NEBLK
            IF (MAT(5,IBLK) .NE. 1) GOTO 40
            IELBEG = MAT(3,IBLK)
            IELEND = MAT(4,IBLK)
C
            XYZMIN(1,IBLK) = CORDSP(IX(1,IELBEG),1)
            XYZMIN(2,IBLK) = CORDSP(IX(1,IELBEG),2)
            XYZMAX(1,IBLK) = CORDSP(IX(1,IELBEG),1)
            XYZMAX(2,IBLK) = CORDSP(IX(1,IELBEG),2)
C
            DO 30 IEL = IELBEG, IELEND
               DO 20 I = 1, NNODES
                  XYZMIN(1,IBLK) = MIN(XYZMIN(1,IBLK),
     *               CORDSP(IX(I,IEL),1))
                  XYZMIN(2,IBLK) = MIN(XYZMIN(2,IBLK),
     *               CORDSP(IX(I,IEL),2))
                  XYZMAX(1,IBLK) = MAX(XYZMAX(1,IBLK),
     *               CORDSP(IX(I,IEL),1))
                  XYZMAX(2,IBLK) = MAX(XYZMAX(2,IBLK),
     *               CORDSP(IX(I,IEL),2))
   20          CONTINUE
   30       CONTINUE
   40    CONTINUE
C
      ELSE
         DO 70 IBLK = 1, NEBLK
            IF (MAT(5,IBLK) .NE. 1) GOTO 70
            IELBEG = MAT(3,IBLK)
            IELEND = MAT(4,IBLK)
C
            XYZMIN(1,IBLK) = CORDSP(IX(1,IELBEG),1)
            XYZMIN(2,IBLK) = CORDSP(IX(1,IELBEG),2)
            XYZMIN(3,IBLK) = CORDSP(IX(1,IELBEG),3)
            XYZMAX(1,IBLK) = CORDSP(IX(1,IELBEG),1)
            XYZMAX(2,IBLK) = CORDSP(IX(1,IELBEG),2)
            XYZMAX(3,IBLK) = CORDSP(IX(1,IELBEG),3)
C
            DO 60 IEL = IELBEG, IELEND
               DO 50 I = 1, NNODES
                  XYZMIN(1,IBLK) = MIN(XYZMIN(1,IBLK),
     *               CORDSP(IX(I,IEL),1))
                  XYZMIN(2,IBLK) = MIN(XYZMIN(2,IBLK),
     *               CORDSP(IX(I,IEL),2))
                  XYZMIN(3,IBLK) = MIN(XYZMIN(3,IBLK),
     *               CORDSP(IX(I,IEL),3))
                  XYZMAX(1,IBLK) = MAX(XYZMAX(1,IBLK),
     *               CORDSP(IX(I,IEL),1))
                  XYZMAX(2,IBLK) = MAX(XYZMAX(2,IBLK),
     *               CORDSP(IX(I,IEL),2))
                  XYZMAX(3,IBLK) = MAX(XYZMAX(3,IBLK),
     *               CORDSP(IX(I,IEL),3))
   50          CONTINUE
   60       CONTINUE
   70    CONTINUE
      END IF
      DO 90 I=1,NDIM
         OVMIN(I) = XYZMIN(I,1)
         OVMAX(I) = XYZMAX(I,1)
         DO 80 J=2,NEBLK
            OVMIN(I) = MIN(OVMIN(I), XYZMIN(I,J))
            OVMAX(I) = MAX(OVMAX(I), XYZMAX(I,J))
   80    CONTINUE
   90 CONTINUE
C
      ENG1 = ENGNOT(TREAD,2)
      DO 120 IO=IOMIN, IOMAX
         IF (EXODUS) WRITE (IO, 140) ENG1
         IF (NDIM .EQ. 2) THEN
            WRITE (IO, 150)
            DO 100 ITMP=1,NEBLK
               I = MAT(6, ITMP)
               IF (MAT(5,I) .NE. 1) GOTO 100
               WRITE (IO, 170) MAT(1,I),(XYZMIN(J,I),J=1,NDIM),
     *            (XYZMAX(J,I),J=1,NDIM),
     *            ((XYZMAX(J,I)-XYZMIN(J,I)),J=1,NDIM)
  100       CONTINUE
         ELSE
            WRITE (IO, 160)
            DO 110 ITMP=1,NEBLK
               I = MAT(6, ITMP)
               IF (MAT(5,I) .NE. 1) GOTO 110
               WRITE (IO, 180) MAT(1,I),(XYZMIN(J,I),J=1,NDIM),
     *            (XYZMAX(J,I),J=1,NDIM),
     *            ((XYZMAX(J,I)-XYZMIN(J,I)),J=1,NDIM)
  110       CONTINUE
         END IF
         IF (NEBLK .GT. 1) THEN
            IF (NDIM .EQ. 2) THEN
               WRITE (IO, 190) (OVMIN(I),I=1,NDIM),(OVMAX(I),I=1,NDIM),
     *              ((OVMAX(I)-OVMIN(I)),I=1,NDIM)
            ELSE
               WRITE (IO, 200) (OVMIN(I),I=1,NDIM),(OVMAX(I),I=1,NDIM),
     *              ((OVMAX(I)-OVMIN(I)),I=1,NDIM)
            END IF
         END IF
  120 CONTINUE
      IF (EXODUS) GO TO 10
  130 CONTINUE
      RETURN

  140 FORMAT(//5X,'Time = ',A16)
  150 FORMAT(/5X,'Material',T22,'X',T36,'Y')
  160 FORMAT(/5X,'Material',T22,'X',T36,'Y',T50,'Z')
  170 FORMAT (I10,3X,2(2X,1PE12.5),'  Minimum',/
     *   ,13X,2(2X,1PE12.5),'  Maximum',/
     *   ,13X,2(2X,1PE12.5),'  Range',/)
  180 FORMAT (I10,3X,3(2X,1PE12.5),'  Minimum',/
     *   ,13X,3(2X,1PE12.5),'  Maximum',/
     *   ,13X,3(2X,1PE12.5),'  Range',/)
  190 FORMAT (5X,'Limits for Total Body:',/
     *   ,13X,2(2X,1PE12.5),'  Minimum',/
     *   ,13X,2(2X,1PE12.5),'  Maximum',/
     *   ,13X,2(2X,1PE12.5),'  Range',/)
  200 FORMAT (5X,'Limits for Total Body:',/
     *   ,13X,3(2X,1PE12.5),'  Minimum',/
     *   ,13X,3(2X,1PE12.5),'  Maximum',/
     *   ,13X,3(2X,1PE12.5),'  Range',/)
      END

