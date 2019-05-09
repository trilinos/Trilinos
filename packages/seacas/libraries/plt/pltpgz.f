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

C $Id: pltpgz.f,v 1.2 1997/11/18 20:03:15 gdsjaar Exp $
C $Log: pltpgz.f,v $
C Revision 1.2  1997/11/18 20:03:15  gdsjaar
C Fixed problem accessing array outside of valid bounds.
C Fixes coredump problem on DEC
C
C Revision 1.1  1993/07/16 16:49:06  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      FUNCTION PLTPGZ(N,X,Y,Z,XQ,YQ)
      DIMENSION X(3),Y(3),Z(3),D(2)

      X21 = X(2) - X(1)
      X31 = X(3) - X(1)
      Y21 = Y(2) - Y(1)
      Y31 = Y(3) - Y(1)
      Z21 = Z(2) - Z(1)
      Z31 = Z(3) - Z(1)
      A = Y31*Z21 - Z31*Y21
      B = Z31*X21 - X31*Z21
      C = X31*Y21 - Y31*X21
      IF (C.EQ.0.) THEN
         DO 2200 I = 1,N - 2
            D(1) = X(I+1) - X(I)
            D(2) = Y(I+1) - Y(I)
            DDD = D(1)*D(1) + D(2)*D(2)
            IF (DDD.NE.0.) THEN
               GO TO 2210

            END IF

 2200    CONTINUE
 2210    CONTINUE
         IF (DDD.EQ.0.) THEN
            PLTPGZ = 0.0
            RETURN

         END IF

         DDP = D(1)*XQ + D(2)*YQ
         DDP1 = D(1)*X(1) + D(2)*Y(1)
         ALPHA = (DDP-DDP1)/DDD
         PLTPGZ = Z(1) + ALPHA* (Z(2)-Z(1))
         RETURN

      END IF

      PLTPGZ = ((A* (X(1)-XQ)+B* (Y(1)-YQ))/C) + Z(1)
      RETURN

      END
