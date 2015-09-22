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

C $Log: elecor.f,v $
C Revision 1.2  2009/03/25 12:36:43  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:00:12  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:49:44  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE ELECOR (NDIM, NELBLK, LEN, NLNK, LINK,
     &   XN, YN, ZN, XE, YE, ZE)
C=======================================================================

C   --*** ELECOR *** (BLOT) Calculate element centers
C   --   Written by Amy Gilkey, revised 10/23/87
C   --
C   --ELECOR calculates the element centers from the nodal coordinates.
C   --
C   --Parameters:
C   --   NDIM - IN - the number of dimensions
C   --   NELBLK - IN - the number of element blocks
C   --   LEN - IN - the cumulative element counts by element block
C   --   NLNK - IN - the number of nodes per element
C   --   LINK - IN - the connectivity array; connectivity all zero if
C   --      element is undefined
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   XE, YE, ZE - OUT - the element center coordinates

      INTEGER LEN(0:*), LINK(*)
      INTEGER NLNK(*)
      REAL XN(*), YN(*), ZN(*)
      REAL XE(*), YE(*), ZE(*)

      DO 170 IELB = 1, NELBLK
         IF (NLNK(IELB) .LE. 0) GOTO 170

         DIVLNK = 1.0 / NLNK(IELB)

         IXL0 = IDBLNK (IELB, 0, LEN, NLNK) - 1
         IF (NDIM .EQ. 2) THEN
            IF (NLNK(IELB) .NE. 8) THEN
               DO 110 IEL = LEN(IELB-1)+1, LEN(IELB)
                  X0 = 0.0
                  Y0 = 0.0
                  IF (LINK(IXL0+1) .NE. 0) THEN
                     DO 100 I = 1, NLNK(IELB)
                        X0 = X0 + XN(LINK(IXL0+I))
                        Y0 = Y0 + YN(LINK(IXL0+I))
  100                CONTINUE
                  END IF
                  XE(IEL) = X0 * DIVLNK
                  YE(IEL) = Y0 * DIVLNK
                  IXL0 = IXL0 + NLNK(IELB)
  110          CONTINUE

            ELSE IF (NLNK(IELB) .EQ. 8) THEN
C            --Note that 8-node elements are numbered with corners at 1-3-5-7
               DO 140 IEL = LEN(IELB-1)+1, LEN(IELB)
                  X0 = 0.0
                  Y0 = 0.0
                  IF (LINK(IXL0+1) .NE. 0) THEN
                     DO 120 I = 1, NLNK(IELB), 2
                        X0 = X0 - 0.25 * XN(LINK(IXL0+I))
                        Y0 = Y0 - 0.25 * YN(LINK(IXL0+I))
  120                CONTINUE
                     DO 130 I = 2, NLNK(IELB), 2
                        X0 = X0 + 0.5 * XN(LINK(IXL0+I))
                        Y0 = Y0 + 0.5 * YN(LINK(IXL0+I))
  130                CONTINUE
                  END IF
                  XE(IEL) = X0
                  YE(IEL) = Y0
                  IXL0 = IXL0 + NLNK(IELB)
  140          CONTINUE
            END IF

         ELSE IF (NDIM .EQ. 3) THEN
            DO 160 IEL = LEN(IELB-1)+1, LEN(IELB)
               X0 = 0.0
               Y0 = 0.0
               Z0 = 0.0
               IF (LINK(IXL0+1) .NE. 0) THEN
                  DO 150 I = 1, NLNK(IELB)
                     X0 = X0 + XN(LINK(IXL0+I))
                     Y0 = Y0 + YN(LINK(IXL0+I))
                     Z0 = Z0 + ZN(LINK(IXL0+I))
  150             CONTINUE
               END IF
               XE(IEL) = X0 * DIVLNK
               YE(IEL) = Y0 * DIVLNK
               ZE(IEL) = Z0 * DIVLNK
               IXL0 = IXL0 + NLNK(IELB)
  160       CONTINUE
         END IF
  170 CONTINUE

      RETURN
      END
