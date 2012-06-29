C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
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

C $Log: gaussf.f,v $
C Revision 1.2  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:01:23  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:50:56  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE GAUSSF (IS3DIM, NLNKF, LINKF1, XN, YN, ZN,
     &   XGAUSS, YGAUSS, ZGAUSS)
C=======================================================================

C   --*** GAUSSF *** (DETOUR) Calculate gauss point for a face
C   --   Written by Amy Gilkey - revised 10/20/87
C   --
C   --GAUSSF calculates the 4 element gauss points given the coordinates
C   --of the 4 nodes.
C   --
C   --Parameters:
C   --   IS3DIM - IN - true iff 3D versus 2D
C   --   NLNKF - IN - the number of nodes per face; must be 4
C   --   LINKF1 - IN - the connectivity for the face
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   XGAUSS, YGAUSS, ZGAUSS - OUT - the gauss point coordinates

      PARAMETER (S = .577350269189626)

      LOGICAL IS3DIM
      INTEGER LINKF1(NLNKF)
      REAL XN(*), YN(*), ZN(*)
      REAL XGAUSS(4), YGAUSS(4), ZGAUSS(4)

      REAL WT(4,4)
      LOGICAL FIRST
      SAVE FIRST, WT
C      --FIRST - true iff first time in routine
C      --WT - the gauss point constants

      DATA FIRST /.TRUE./

      IF (FIRST) THEN
         FIRST = .FALSE.

C      --Determine gauss point constants

         DO 100 N = 1, 4
            IF (N .EQ. 1) THEN
               SI = -S
               TI = -S
            ELSE IF (N .EQ. 2) THEN
               SI = +S
               TI = -S
            ELSE IF (N .EQ. 3) THEN
               SI = +S
               TI = +S
            ELSE IF (N .EQ. 4) THEN
               SI = -S
               TI = +S
            END IF
            WT(1,N) = .25 * (1-SI) * (1-TI)
            WT(2,N) = .25 * (1+SI) * (1-TI)
            WT(3,N) = .25 * (1+SI) * (1+TI)
            WT(4,N) = .25 * (1-SI) * (1+TI)
  100    CONTINUE
      END IF

C   --Determine gauss point coordinates

      DO 110 N = 1, 4
         XGAUSS(N) = WT(1,N)*XN(LINKF1(1)) + WT(2,N)*XN(LINKF1(2))
     &      + WT(3,N)*XN(LINKF1(3)) + WT(4,N)*XN(LINKF1(4))
         YGAUSS(N) = WT(1,N)*YN(LINKF1(1)) + WT(2,N)*YN(LINKF1(2))
     &      + WT(3,N)*YN(LINKF1(3)) + WT(4,N)*YN(LINKF1(4))
         IF (IS3DIM) THEN
            ZGAUSS(N) = WT(1,N)*ZN(LINKF1(1)) + WT(2,N)*ZN(LINKF1(2))
     &         + WT(3,N)*ZN(LINKF1(3)) + WT(4,N)*ZN(LINKF1(4))
         END IF
  110 CONTINUE

      RETURN
      END
