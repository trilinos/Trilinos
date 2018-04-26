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

C $Id: lxrnl.f,v 1.3 1993/07/16 18:28:52 gdsjaar Exp $ 
C $Log: lxrnl.f,v $
C Revision 1.3  1993/07/16 18:28:52  gdsjaar
C Changed real*8 to double precision
C
c Revision 1.2  1993/07/16  18:07:49  gdsjaar
c Added external pltblk statements so that linkers would pull in block
c data subprogram to initialize constants.
c
c Revision 1.1  1993/07/16  16:46:47  gdsjaar
c Changed plt to library rather than single source file.
c 
C=======================================================================
      LOGICAL FUNCTION LXRNL(VAL,N,CH)
      DOUBLE PRECISION VAL(*)
      CHARACTER CH
      LOGICAL LDUM,LXGTCH,LXGTWH,LXREAL
      DOUBLE PRECISION XX,XY

      I = 1
 2320 IF (.NOT. (I.LE.N)) GO TO 2340
      VAL(I) = 0.
      I = I + 1
      GO TO 2320

 2340 CONTINUE
      N = 0
      LXRNL = .TRUE.
 2350 CONTINUE
      LDUM = LXGTWH(CH)
      IF (LXREAL(VAL(N+1),CH)) THEN
         N = N + 1
         LDUM = LXGTWH(CH)
         IF (LXGTCH('#',CH)) THEN
            RETURN

         END IF

         IF (LXGTCH(',',CH)) THEN
            GO TO 2360

         END IF

         IF (CH.EQ.CHAR(0)) THEN
            RETURN

         END IF

         IF (LXGTCH('*',CH)) THEN
            LDUM = LXGTWH(CH)
            XX = VAL(N)
            IF (.NOT.LXREAL(XY,CH)) THEN
               LXRNL = .FALSE.
               RETURN

            END IF

            M = INT(XX + .1)
            N0 = N
 2380       IF (.NOT. (N.LT.M+N0)) GO TO 2400
            VAL(N) = XY
            N = N + 1
            GO TO 2380

 2400       CONTINUE
            LDUM = LXGTWH(CH)
            N = N0 + MAX(M-1,0)
            IF (LXGTCH(',',CH)) THEN
               GO TO 2360

            END IF

            IF (LXGTCH('#',CH)) THEN
               RETURN

            END IF

            IF (CH.EQ.CHAR(0)) THEN
               RETURN

            END IF

         END IF

      ELSE IF (LXGTCH(',',CH)) THEN
         VAL(N+1) = 0.
         N = N + 1

      ELSE IF (LXGTCH('#',CH)) THEN
         RETURN

      ELSE IF (CH.EQ.CHAR(0)) THEN
         IF (N.EQ.0) THEN
            RETURN

         END IF

         N = N + 1
         VAL(N) = 0.
         RETURN

      ELSE
         LXRNL = .FALSE.
         RETURN

      END IF

 2360 GO TO 2350

      END
